# Uncomment one of the following to switch between debug and opt mode
OPT = -O3 -DNDEBUG
#OPT = -g -ggdb

CXXFLAGS += -fno-strict-aliasing -Wall -std=c++11 -I. -I../src/ \
    -I../src/bloom/ -I../src/cuckoo/ -I../src/gcs \
    -I../src/gqf/ -I../src/morton/ -I../src/xorfilter \
    $(OPT)

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),aarch64)
        CXXFLAGS +=
else
ifeq ($(UNAME_P),unknown)
        CXXFLAGS +=
else
	CXXFLAGS += -march=native
endif
endif

LDFLAGS = -Wall

HEADERS = $(wildcard ../src/*.h \
    ../src/bloom/*.h ../src/cuckoo/*.h ../src/gcs/*.h \
    ../src/gqf/*.h ../src/morton/*.h ../src/xorfilter/*.h \
    ) *.h

.PHONY: all

BINS = bulk-insert-and-query.exe

all: $(BINS)

clean:
	/bin/rm -f $(BINS)

%.exe: %.cc ${HEADERS}  Makefile
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
