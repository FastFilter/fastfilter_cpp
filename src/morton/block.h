/*
Copyright (c) 2019 Advanced Micro Devices, Inc.
 
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
 
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Author: Alex D. Breslow 
        Advanced Micro Devices, Inc.
        AMD Research
*/
#ifndef _BLOCK_H
#define _BLOCK_H

#include <stdint.h>

#include <cstring>
#include <sstream>
#include <iostream>

#include "util.h"
#include "vector_types.h"

#ifndef INLINE
#define INLINE __attribute__((always_inline)) inline
#endif

namespace Morton{

const uint64_t bytes_to_bits = 8lu;

template<uint64_t block_size_bits, 
         uint64_t entry_size_bits, 
         typename atom_t> // What composes the block, e.g., 8 atoms
struct Block{ // Assuming block size is a multiple of atom_t's size in bytes
  #include "vector_types.h" // Needs to get included here to properly 
                            // define an atom_t vector
  private:
    const static uint64_t atom_size_bits = bytes_to_bits * sizeof(atom_t);
    // Limits entry size to 63 bits
    const static uint64_t read_mask = (static_cast<atom_t>(1) << entry_size_bits) - 1;
    const static uint64_t atoms_per_block = block_size_bits / (atom_size_bits);
    // You need this for cases where you don't use the final several bits of 
    // an atom because the entries don't evenly divide into it.
    const static uint64_t entries_per_atom = atom_size_bits / entry_size_bits;
    atom_t _block_storage[atoms_per_block]{}; // FIXME: <= 64 bits
    const static bool explicit_print = true;
  public:
    inline atom_t operator[](uint32_t atom_index) const{
      return _block_storage[atom_index];
    }

    // Items don't cross atom boundaries.
    inline void add(uint64_t index, atom_t item){
      // Which atom we need to examine
      uint64_t atom_index = index / entries_per_atom;
      // Where to start reading in the atom
      uint64_t atom_bit_offset = (index % entries_per_atom) * entry_size_bits;
      _block_storage[atom_index] &= ~(read_mask << atom_bit_offset); // clear
      _block_storage[atom_index] |= (item << atom_bit_offset); // set
    }

    INLINE bool read_bit(uint64_t raw_index_in_bits) const{
      uint64_t atom_index = raw_index_in_bits / (sizeof(atom_t) * 8);
      uint64_t atom_offset = raw_index_in_bits % (sizeof(atom_t) * 8);
      return (_block_storage[atom_index] >> atom_offset) & 0x1;
    }

    INLINE atom_t read_byte(uint64_t byte_index) const{
      uint64_t atom_index = byte_index / (sizeof(atom_t));
      uint64_t local_byte_index = byte_index % (sizeof(atom_t));
      return (_block_storage[atom_index] >> 
        (sizeof(atom_t) * local_byte_index)) & 0xff;
    }

    INLINE void set_bit(uint64_t raw_index_in_bits, atom_t value){
      uint64_t atom_index = raw_index_in_bits / (sizeof(atom_t) * 8);
      uint64_t atom_offset = raw_index_in_bits % (sizeof(atom_t) * 8);
      atom_t mask = (~static_cast<atom_t>(0)) ^ (value << atom_offset);
      _block_storage[atom_index] &= mask;
    }

    // Bit's initial state 0 and value of 0 ---> no change (0)
    // Bit's initial state 0 and value of 1 ---> set to 1  (0 ---> 1)
    // Bit's initial state 1 and any value  ---> no change (1)  
    INLINE void sticky_set_bit(uint64_t raw_index_in_bits, uint64_t value){
      uint64_t atom_index = raw_index_in_bits / (sizeof(atom_t) * 8);
      uint64_t atom_offset = raw_index_in_bits % (sizeof(atom_t) * 8);
      _block_storage[atom_index] |= value << atom_offset;
    }

    // Templated version of the function below it
    // Note that like the version below it, it assumes that you only need to 
    // read a single word of type T from _block_storage.
    template<class T> // Set to type of desired read size
    inline void add_t(uint64_t raw_offset_bits, uint64_t field_width_bits, 
      uint64_t index, T item){
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      uint8_t word_index = global_index / (sizeof(T) * 8ULL);
      T read_mask_specific = (static_cast<T>(1) << field_width_bits) - 
        static_cast<T>(1);
      uint8_t word_bit_offset = global_index % (sizeof(T) * 8ULL);
      T* block_storage = reinterpret_cast<T*>(_block_storage);
      block_storage[word_index] &= ~(read_mask_specific << word_bit_offset);
      block_storage[word_index] |= item << word_bit_offset;
    }

    // This version of add assumes that the read doesn't need to cross 
    // atom boundaries.  There is but a single atom that we need to access.
    // FIXME: As written, it currently only supports accessing atom 0 if the 
    // field width doesn't evenly divide a block.
    inline void add(uint64_t raw_offset_bits, uint64_t field_width_bits,
      uint64_t index, atom_t item){
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      // Get atom index of where the item is located.  It could be one or two 
      // atoms if the entry size is less than or equal to atom_size_bits.
      // If it's greater than an atom in size, then it could span more than 
      // 2 atoms.  I'm not going to support that here.
      uint8_t atom = global_index / atom_size_bits;

      // May be different than the default read mask if field_width_bits 
      // differs from what is stored in the struct
      atom_t read_mask_specific = (static_cast<atom_t>(1) << field_width_bits) - static_cast<atom_t>(1);
      // The entry's LSB and its index in atom0
      uint8_t atom_bit_offset = global_index % atom_size_bits;

      // Clear the bits for the field
      _block_storage[atom] &= ~(read_mask_specific << atom_bit_offset);
      
      // Store the item
      _block_storage[atom] |= item << atom_bit_offset;
    }

    // A more general implementation of adding elements that allows for 
    // the crossing of word boundaries, having elements begin at a raw 
    // bit offset from the start, and specifying an arbitrary field width in 
    // bits.  This flexibility allows using the block store for storing 
    // several arrays within a block, with potentially different sized elements 
    // in each.  An example would be an array of 3 bit counters followed 
    // immediately or offset by several bits by an array of 5 bit entries. 
    INLINE void add_cross(uint64_t raw_offset_bits, 
      uint64_t field_width_bits, 
      uint64_t index, atom_t item){
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      // Get atom index of where the item is located.  It could be one or two 
      // atoms if the entry size is less than or equal to atom_size_bits.
      // If it's greater than an atom in size, then it could span more than 
      // 2 atoms.  I'm not going to support that here.
      atom_t atom0 = global_index / atom_size_bits;
      atom_t atom1 = (global_index + field_width_bits) / atom_size_bits;
      // How many bits spilled over into atom1
      uint16_t spillover = (global_index + field_width_bits) % atom_size_bits;
      spillover = (atom1 - atom0) ? spillover : 0; // Compare and swap ins?

      // May be different than the default read mask if field_width_bits 
      // differs from what is stored in the struct
      atom_t read_mask_specific = (static_cast<atom_t>(1) << field_width_bits) - 1;

      // The entry's LSB and its index in atom0
      uint16_t atom0_bit_offset = global_index % atom_size_bits;

      _block_storage[atom0] &=  ~(read_mask_specific << atom0_bit_offset); // clear

      atom_t shifted_read_mask = read_mask_specific >> (field_width_bits - spillover);
      _block_storage[atom1] &= ~(shifted_read_mask); // clear if spillover

      _block_storage[atom0] |= (item << atom0_bit_offset); // Store
      _block_storage[atom1] |= (item >> (field_width_bits - spillover));
    }

    // Entries can cross atom boundaries.
    INLINE void add_cross(uint64_t index, atom_t item){
      add_cross(0, entry_size_bits, index, item);
    }

    template<class T>
    inline T read_zeroth_word(uint64_t raw_offset_bits,
      uint64_t field_width_bits, uint64_t index) const{
      constexpr T one = static_cast<T>(1);
      uint_fast8_t lsb = raw_offset_bits + field_width_bits * index;
      const T read_mask_specific = (one << field_width_bits) - one;
      const T* block_storage = reinterpret_cast<const T*>(_block_storage);
      return (block_storage[0] >> lsb) & read_mask_specific;
    }

    //__attribute__((always_inline)) 
    inline atom_t read_atom0(uint64_t raw_offset_bits, 
      uint64_t field_width_bits, uint64_t index) const{
      constexpr atom_t one = static_cast<atom_t>(1);
      // Index of the least significant bit that we want to read
      uint_fast8_t lsb = raw_offset_bits + field_width_bits * index;
      // Mask out content to the left of where we want to read
      const atom_t read_mask_specific = (one << field_width_bits) - one;
      return (_block_storage[0] >> lsb) & read_mask_specific;
    }
 
    // This version of read assumes that the read doesn't need to cross 
    // atom boundaries.  There is but a single atom that we need to access.
    // FIXME: As written, it currently only supports accessing atom 0 if the 
    // field width doesn't evenly divide a block.
    inline atom_t read(uint64_t raw_offset_bits, uint64_t field_width_bits,
      uint64_t index) const{
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      // Get atom index of where the item is located.  It could be one or two 
      // atoms if the entry size is less than or equal to atom_size_bits.
      // If it's greater than an atom in size, then it could span more than 
      // 2 atoms.  I'm not going to support that here.
      uint8_t atom = global_index / atom_size_bits;

      // May be different than the default read mask if field_width_bits 
      // differs from what is stored in the struct
      atom_t read_mask_specific = (static_cast<atom_t>(1) << field_width_bits) - static_cast<atom_t>(1);
      // The entry's LSB and its index in atom0
      uint8_t atom_bit_offset = global_index % atom_size_bits;
      // Return the item
      return (_block_storage[atom] >> atom_bit_offset) & 
        read_mask_specific;
    }

    // TODO: Should return an object that is big enough to hold an entry 
    inline atom_t read(uint64_t index) const{
      // Which atom we need to examine
      uint64_t atom_index = index / entries_per_atom;
      // Where to start reading in the atom
      uint64_t atom_bit_offset = (index % entries_per_atom) * entry_size_bits;
      return (_block_storage[atom_index] >> atom_bit_offset) & read_mask;
    }

		// FIXME: Not rigorously tested
    // Note that this code assumes that you're reading a block of items 
    // that is no bigger than 1 plus an atom's size in bits. This code is 
    // useful for reading a block of items that may span multiple atoms 
    // due to their initial offset but which could fit into a single word.
    // It's much more efficient than calling read_cross multiple times on 
    // adjacent items.  The index tells you where to start reading at.
    inline atom_t read_cross_many(uint64_t raw_offset_bits,
      uint64_t field_width_bits, uint64_t index, uint64_t items_read) const{
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      // Get atom index of where the item is located.  It could be one or two 
      // atoms if the entry size is less than or equal to atom_size_bits.
      // If it's greater than an atom in size, then it could span more than 
      // 2 atoms.  I'm not going to support that here.
      uint16_t atom0 = global_index / atom_size_bits;
      uint16_t atom1 = (global_index + field_width_bits * items_read) / atom_size_bits;
      atom_t multi_item_read_mask = (static_cast<atom_t>(1) << (field_width_bits * items_read)) - 1;
      // How many bits spilled over into atom1
      uint16_t spillover = (global_index + field_width_bits * items_read) % atom_size_bits;
      spillover = (atom1 - atom0) ? spillover : 0; // Compare and swap ins?

      // The entry's LSB and its index in atom0
      atom_t atom0_bit_offset = global_index % atom_size_bits;
      // This is either the full entry or the bits of it present in atom 0.
      atom_t lower = (_block_storage[atom0] >> atom0_bit_offset) & multi_item_read_mask;

      atom_t shifted_read_mask = multi_item_read_mask >> (field_width_bits * items_read - spillover);
      // This is nothing or the remaining bits of the entry in atom1.
      atom_t upper = (_block_storage[atom1] & (shifted_read_mask));
      return (upper << (field_width_bits * items_read - spillover)) | lower; 
    }

    // Increments a field by 1
    // !!!NOTE: Doesn't do overflow testing FIXME?  I'm going to rely on the 
    // application to do the check. If you're getting a counter overflow, 
    // something is seriously wrong.
    inline void inc(uint64_t raw_offset_bits, 
      uint64_t field_width_bits, uint64_t index){
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      atom_t atom = global_index / atom_size_bits;
      uint64_t counter_lsb = global_index % atom_size_bits;
      _block_storage[atom] += static_cast<atom_t>(1) << counter_lsb;
    }

    // Key is to propagate the carry bit across words
    inline void inc_cross(uint64_t raw_offset_bits, uint64_t field_width_bits, 
      uint64_t index){
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      uint16_t atom = global_index / atom_size_bits;
      // Bit position of the counter's LSB
      uint16_t counter_lsb = global_index % atom_size_bits;
      
      // Need to get the initial value of the MSB of the atom to know if we 
      // need to perform a carry operation
      uint8_t msb = _block_storage[atom] >> (atom_size_bits - 1);
      _block_storage[atom] += (static_cast<atom_t>(1) << counter_lsb);
      uint8_t msb_after_inc = _block_storage[atom] >> (atom_size_bits - 1);

      // Only perform a carry if the MSB bit flips from 1 to 0. 
      // We don't have to worry about 
      // overflowing into an unmapped page because if the counter spans two 
      // atoms, then the next atom is also valid (if you're using this 
      // function correctly :) ).
      uint8_t carry = (msb == 1) & (msb_after_inc == 0); 
      _block_storage[atom + 1] += carry; // IMPORTANT: This could segfault or 
                                         // corrupt memory if the page after 
                                         // the last atom isn't mapped.
    }

    // A variant of read_cross that is more parametrizable
    INLINE atom_t read_cross(uint64_t raw_offset_bits, 
      uint64_t field_width_bits, uint64_t index) const{
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      // Get atom index of where the item is located.  It could be one or two 
      // atoms if the entry size is less than or equal to atom_size_bits.
      // If it's greater than an atom in size, then it could span more than 
      // 2 atoms.  I'm not going to support that here.
      uint8_t atom0 = global_index / atom_size_bits;
      uint8_t atom1 = (global_index + field_width_bits) / atom_size_bits;
      // How many bits spilled over into atom1
      uint8_t spillover = (global_index + field_width_bits) % atom_size_bits;
      spillover = (atom1 - atom0) ? spillover : 0; // Compare and swap ins?

      // May be different than the default read mask if field_width_bits 
      // differs from what is stored in the struct
      atom_t read_mask_specific = (static_cast<atom_t>(1) << field_width_bits) - static_cast<atom_t>(1);

      // The entry's LSB and its index in atom0
      uint8_t atom0_bit_offset = global_index % atom_size_bits;
      // This is either the full entry or the bits of it present in atom 0.
      atom_t lower = (_block_storage[atom0] >> atom0_bit_offset) & 
        read_mask_specific;

      atom_t shifted_read_mask = read_mask_specific >> (field_width_bits - spillover);
      // This is nothing or the remaining bits of the entry in atom1.
      atom_t upper = (_block_storage[atom1] & (shifted_read_mask));
      return (upper << (field_width_bits - spillover)) | lower; 
    }

    // Simple case where the read is not occurring on an array that doesn't 
    // have an initial bit offset and the field width is just entry_size_bits
    INLINE atom_t read_cross(uint64_t index) const{
      return read_cross(0, entry_size_bits, index);
    }

    // A simpler version to understand
    // It oddly gives the same performance.
    INLINE void add_cross_left_displace_v2(uint64_t raw_offset_bits, 
      uint64_t field_width_bits, uint64_t index, atom_t item){
      constexpr atom_t one = 1;
      uint64_t global_index = index * field_width_bits + raw_offset_bits;

      // Which atom to begin shifting at
      uint8_t atom0 = global_index / atom_size_bits;
      // The first atom is special.  Part shifts and part likely doesn't.
      // Grab the part that's going to shift 
      uint8_t atom0_bit_offset = global_index % atom_size_bits;

      Block<block_size_bits, entry_size_bits, atom_t> shifted_block = *this;
      
      for(uint64_t atom = 0; atom < atoms_per_block; atom++){
        shifted_block._block_storage[atom] <<= field_width_bits;
      }
      for(uint64_t atom = 1; atom < atoms_per_block; atom++){
        shifted_block._block_storage[atom] |= _block_storage[atom - 1] >> (atom_size_bits - field_width_bits);
      } 

      // Only keep the part that doesn't move (everything to the right of the 
      // the offset position
      atom_t atom0_read_mask = (one << atom0_bit_offset) - one;

      // 1) Every atom before where the change took place stays in position
      // No change (a no-op)

      // 2) Blend atom0 where the change occurred 
      // Low order bits don't change
      // High order bits after proceeding the shift do change
      _block_storage[atom0] &= atom0_read_mask;
      _block_storage[atom0] |= ((~atom0_read_mask) & shifted_block._block_storage[atom0]);

      // 3) Every atom after where the change took place is just shifted
      for(uint64_t atom = atom0 + 1; atom < atoms_per_block; atom++){
        _block_storage[atom] = shifted_block[atom];
      }
      
      // 4) Write the new element into position (clears slot prior to writing)
      add_cross(raw_offset_bits, field_width_bits, index, item);
    }

    INLINE void add_left_displace(uint64_t raw_offset_bytes,
      uint64_t field_width_bytes, uint64_t index, atom_t item, 
      uint64_t region_size_bytes){
      
      uint8_t* source_address = reinterpret_cast<uint8_t*>(&_block_storage[0]) 
        + raw_offset_bytes + field_width_bytes * index;
      uint8_t* destination_address = source_address + field_width_bytes;
      uint64_t bytes_moved = region_size_bytes - field_width_bytes - index * field_width_bytes; 
      memmove(destination_address, source_address, bytes_moved);
      memcpy(source_address, &item, field_width_bytes);
    }

    INLINE void add_cross_left_displace(uint64_t raw_offset_bits, 
      uint64_t field_width_bits, uint64_t index, atom_t item){
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      // Which atom to begin shifting at
      uint8_t atom0 = global_index / atom_size_bits;
      // The first atom is special.  Part shifts and part likely doesn't.
      // Grab the part that's going to shift 
      uint8_t atom0_bit_offset = global_index % atom_size_bits;

      Block<block_size_bits, entry_size_bits, atom_t> shifted_block = *this;
      
      for(uint64_t atom = 0; atom < atoms_per_block; atom++){
        shifted_block._block_storage[atom] <<= field_width_bits;
      }
      for(uint64_t atom = 1; atom < atoms_per_block; atom++){
        shifted_block._block_storage[atom] |= _block_storage[atom - 1] >> (atom_size_bits - field_width_bits);
      } 

      // Clear the entire atom
      constexpr atom_t full_clear_mask = static_cast<atom_t>(0);
      
      // Preserve the entire atom
      constexpr atom_t full_preserve_mask = ~full_clear_mask;
 
      // Only keep the part that doesn't move (everything to the right of the 
      // the offset position
      atom_t atom0_read_mask = (static_cast<atom_t>(1) << atom0_bit_offset) - 
        static_cast<atom_t>(1);
      for(uint64_t atom = 0; atom < atoms_per_block; atom++){
        atom_t alternate_clear_mask = (atom == atom0) ? atom0_read_mask :
          full_preserve_mask;
        atom_t clear_mask = (atom > atom0) ? full_clear_mask : alternate_clear_mask;
        atom_t atom_read_mask = ~clear_mask;
        _block_storage[atom] &= clear_mask;
        _block_storage[atom] |= ((atom_read_mask) & shifted_block._block_storage[atom]);
      }
      // Write the new element into position
      add_cross(raw_offset_bits, field_width_bits, index, item);
    }

    // Specialized version for when the field and the things that you are 
    // deleting and shifting are all multiples of a byte and so is the alignment
    INLINE void del_right_displace(uint64_t raw_offset_bytes, uint64_t 
      field_width_bytes, uint64_t index, uint64_t region_size_bytes){
      // Base address for the region that we're shifting
      uint8_t* base_address = reinterpret_cast<uint8_t*>(&_block_storage[0]) 
        + raw_offset_bytes;
      uint8_t* source_address = base_address + field_width_bytes * (index + 1);
      uint8_t* destination_address = source_address - field_width_bytes;
      uint64_t bytes_moved = region_size_bytes - field_width_bytes * (index + 1);
      memmove(destination_address, source_address, bytes_moved); 
      memset(base_address + region_size_bytes - field_width_bytes, 0x0, 
        field_width_bytes);  
    }


    // General version for when can fingerprints cross word boundaries
    INLINE void del_cross_right_displace(uint64_t raw_offset_bits,
      uint64_t field_width_bits, uint64_t index){
      uint64_t global_index = index * field_width_bits + raw_offset_bits;
      // Which atom to begin shifting at
      uint8_t atom0 = global_index / atom_size_bits;
      // The first atom is special.  Part shifts and part likely doesn't.
      // Grab the part that's going to shift 
      uint8_t atom0_bit_offset = global_index % atom_size_bits;

      Block<block_size_bits, entry_size_bits, atom_t> shifted_block = *this;
     
      // Shift right 
      for(uint64_t atom = 0; atom < atoms_per_block; atom++){
        shifted_block._block_storage[atom] >>= field_width_bits;
      }
      // this -> _block_storage[atom]'s low order bits needed to OR'ed into
      // shifted_block._block_starge[atom - 1]'s high order bits
      for(uint64_t atom = 1; atom < atoms_per_block; atom++){
        shifted_block._block_storage[atom - 1] |= _block_storage[atom] << (atom_size_bits - field_width_bits);
      } 

      // What follows is a merge operation that is exactly the same as 
      // what appears in add_cross_left_displace
      // ***** START OF MERGE ****
      // Clear the entire atom
      constexpr atom_t full_clear_mask = static_cast<atom_t>(0);
      
      // Preserve the entire atom
      constexpr atom_t full_preserve_mask = ~full_clear_mask;
 
      // Only keep the part that doesn't move (everything to the right of the 
      // the offset position
      atom_t atom0_read_mask = (static_cast<atom_t>(1) << atom0_bit_offset) - 
        static_cast<atom_t>(1);
      for(uint64_t atom = 0; atom < atoms_per_block; atom++){
        atom_t alternate_clear_mask = (atom == atom0) ? atom0_read_mask :
          full_preserve_mask;
        atom_t clear_mask = (atom > atom0) ? full_clear_mask : alternate_clear_mask;
        atom_t atom_read_mask = ~clear_mask;
        _block_storage[atom] &= clear_mask;
        _block_storage[atom] |= ((atom_read_mask) & shifted_block._block_storage[atom]);
      } // **** END OF MERGE ****
    }

    // Zeroes bits from index raw_offset_bits to raw_offset_bits + len_in_bits - 1
    inline void clear_swath(uint64_t raw_offset_bits, uint64_t len_in_bits){
      constexpr uint64_t one = 1; 
      uint8_t* zeroth_byte_pointer = reinterpret_cast<uint8_t*>(&_block_storage[0]);
      int64_t first_full_byte = (raw_offset_bits / 8) + ((raw_offset_bits % 8) != 0);
      int64_t leading_partial_byte = raw_offset_bits / 8;
      int64_t last_full_byte = ((raw_offset_bits + len_in_bits) / 8) - ((raw_offset_bits % 8) != 0);
      int64_t trailing_partial_byte = (raw_offset_bits + len_in_bits) / 8;
      uint8_t* first_full_byte_ptr = zeroth_byte_pointer + first_full_byte; 
      if(last_full_byte >= first_full_byte) memset(first_full_byte_ptr, 0x0, last_full_byte - first_full_byte + one);
      uint8_t* leading_partial_byte_pointer = zeroth_byte_pointer + leading_partial_byte;
      *leading_partial_byte_pointer &= (one << (raw_offset_bits % 8)) - one;
      uint8_t* trailing_partial_byte_pointer = zeroth_byte_pointer + trailing_partial_byte;
      *trailing_partial_byte_pointer &= ~((one << ((raw_offset_bits + len_in_bits) % 8)) - one);
    } 

    inline void clear(uint64_t index){
      // Which atom we need to examine
      uint64_t atom_index = index / entries_per_atom;
      // Where to start reading in the atom
      uint64_t atom_bit_offset = (index % entries_per_atom) * entry_size_bits;
      _block_storage[atom_index] &= ~(read_mask << atom_bit_offset); // clear
    }

    // Entries can cross atom boundaries.
    inline void clear_cross(uint64_t index){
      uint64_t global_index = index * entry_size_bits;
      // Get atom index of where the item is located.  It could be one or two 
      // atoms if the entry size is less than or equal to atom_size_bits.
      // If it's greater than an atom in size, then it could span more than 
      // 2 atoms.  I'm not going to support that here.
      atom_t atom0 = global_index / atom_size_bits;
      atom_t atom1 = (global_index + entry_size_bits) / atom_size_bits;
      // How many bits spilled over into atom1
      uint64_t spillover = (global_index + entry_size_bits) % atom_size_bits;
      spillover = (atom1 - atom0) ? spillover : 0; // Compare and swap ins?

      // The entry's LSB and its index in atom0
      uint64_t atom0_bit_offset = global_index % atom_size_bits;

      _block_storage[atom0] &=  ~(read_mask << atom0_bit_offset); // clear

      uint64_t shifted_read_mask = read_mask >> (entry_size_bits - spillover);
      _block_storage[atom1] &= ~(shifted_read_mask); // clear if spillover
    }

    // Prints out the block storage bit by bit with a space every "spacing"
    // digits
    std::string block_storage_as_bit_string(uint32_t spacing){
      std::stringstream ss;
      for(int32_t atom_id = atoms_per_block - 1; atom_id > -1; 
      atom_id--){
        for(int32_t i = sizeof(atom_t) * 8 - 1; i > -1; i--){
          ss << ((_block_storage[atom_id] >> i) & 1 ? '1' : '0');
          if((i + atom_size_bits * atom_id) % spacing == 0){
            ss << ' ';
          }
        }
      }
      return ss.str();
    }

    friend std::ostream& operator<<(std::ostream &os, 
      const Block<block_size_bits, entry_size_bits, atom_t>& b){
     
      // Both of these print the most significant bits to the left.  So the 
      // zeroth entry in the zeroth atom is the furthest to the right. 
      if(b.explicit_print){
        for(int64_t i = atoms_per_block - 1; i >= 0; i--){
          for(int64_t j = atom_size_bits - 1; j >= 0; j--){
            char digit = ((b._block_storage[i] >> (j)) & 1) ? '1' : '0';
            os << digit;
          }
          os << " ";
        }
      }
      else{ // Print the actual items rather than the storage.
        for(int64_t i = atoms_per_block - 1; i >= 0; i--){
          for(int64_t j = entries_per_atom - 1; j >= 0; j--){
            os << b.read(i * entries_per_atom + j) << " "; 
          }
          os << " ";
        }
      }
      return os;
    }
};

}// Morton namespace

#endif // _BLOCK_H
