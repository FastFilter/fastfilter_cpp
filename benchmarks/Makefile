# Uncomment one of the following to switch between debug and opt mode
OPT = -O3 -DNDEBUG
#OPT = -g -ggdb

CXXFLAGS += -fno-strict-aliasing -Wall -std=c++11 -I. -I../src/ $(OPT) 

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),aarch64)
        CXXFLAGS += 
else
	CXXFLAGS += -march=native
endif

LDFLAGS = -Wall 

HEADERS = $(wildcard ../src/*.h) *.h



.PHONY: all

BINS = bulk-insert-and-query.exe

all: $(BINS)

clean:
	/bin/rm -f $(BINS)

%.exe: %.cc ${HEADERS}  Makefile
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
