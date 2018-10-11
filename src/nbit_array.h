#ifndef NBIT_ARRAY_H_
#define NBIT_ARRAY_H_

//namespace nbit_array {

template <typename ItemType>
class UIntArray {
    size_t byteCount;
    ItemType* data;
public:
    UIntArray(size_t size) {
        byteCount = sizeof(ItemType[size]);
        data = new ItemType[size]();
    }
    ~UIntArray() {
        delete[] data;
    }
    inline ItemType get(size_t index) {
        return data[index];
    }
    inline void set(size_t index, ItemType value) {
        data[index] = value;
    }
    void bulkSet(uint16_t* source, size_t length) {
        for(int i = 0; i < length; i++) {
            data[i] = source[i];
        }
    }
    inline ItemType mask(ItemType fingerprint) {
        return fingerprint;
    }
    size_t getByteCount() {
        return byteCount;
    }
};

class UInt12Array {
    size_t byteCount;
    uint8_t* data;
public:
    UInt12Array(size_t size) {
        byteCount = size * 3 / 2 + 32;
        data = new uint8_t[byteCount]();
    }
    ~UInt12Array() {
        delete[] data;
    }
    // the returned value may contain other high-order bits;
    // call mask() to clear them
    inline uint32_t get(size_t index) {
        size_t firstBytePos = (index >> 1) + index;
        uint32_t word;
        memcpy(&word, data + firstBytePos, sizeof(uint32_t));
        return word >> ((index & 1) << 2);
    }
    void bulkSet(uint16_t* source, size_t length) {
        for(size_t i = 0, j = 0; i < length;) {
            uint32_t a = source[i++];
            uint32_t b = source[i++];
            data[j++] = (uint8_t) (a);
            data[j++] = (uint8_t) ((a >> 8) | (b << 4));
            data[j++] = (uint8_t) (b >> 4);
        }
    }
    inline void set(size_t index, uint32_t value) {
        size_t wordpos = (index >> 1) * 3;
        unsigned int wp = (index & 1) * 12;
        size_t offval = (value & 0xfff) << wp;
        size_t offmask = 0xfff << wp;
        uint32_t word;
        memcpy(&word, data + wordpos, sizeof(uint32_t));
        // no need for word & ~(offmask) if it's always 0 at the beginning
        // word = (word) | offval;
        word = (word & ~(offmask)) | offval;
        memcpy(data + wordpos, &word, sizeof(uint32_t));
    }
    inline uint32_t mask(uint32_t fingerprint) {
        return fingerprint & 0xfff;
    }
    size_t getByteCount() {
        return byteCount;
    }
};

template <typename ItemType, size_t bitsPerEntry, uint32_t bitMask = (1 << bitsPerEntry) - 1>
class NBitArray {
    size_t byteCount;
    uint8_t* data;
public:
    NBitArray(size_t size) {
        byteCount = (size * bitsPerEntry + 63 + 128) / 64 * 64 / 8;
        data = new uint8_t[byteCount]();
    }
    ~NBitArray() {
        delete[] data;
    }
    inline ItemType get(size_t index) {
        size_t bitPos = index * bitsPerEntry;
        size_t firstBytePos = (size_t) (bitPos >> 3);
        uint32_t word = __builtin_bswap32(*((uint32_t*) (data + firstBytePos))) >> 8;
        return (ItemType) ((word >> (24 - bitsPerEntry - (bitPos & 7))) & bitMask);
    }
    void bulkSet(uint16_t* source, size_t length) {
        for(size_t i = 0; i < length; i++) {
            set(i, source[i]);
        }
    }
    inline void set(size_t index, ItemType value) {
        size_t bitPos = index * bitsPerEntry;
        size_t firstBytePos = (size_t) (bitPos >> 3);
        uint32_t word = __builtin_bswap32(*((uint32_t*) (data + firstBytePos))) >> 8;
        word &= ~(bitMask << (24 - bitsPerEntry - (bitPos & 7)));
        word |= ((value & bitMask) << (24 - bitsPerEntry - (bitPos & 7)));
        data[firstBytePos] = (uint8_t) (word >> 16);
        data[firstBytePos + 1] = (uint8_t) (word >> 8);
        data[firstBytePos + 2] = (uint8_t) word;
    }
    inline ItemType mask(ItemType fingerprint) {
        return fingerprint & bitMask;
    }
    size_t getByteCount() {
        return byteCount;
    }
};

// }  // namespace n_bit_array

#endif  // NBIT_ARRAY_H_
