

# Uncomment one of the following to switch between debug and opt mode
OPT = -O3 -DNDEBUG -std=c++17
#OPT = -g -ggdb -fsanitize=address -fno-omit-frame-pointer -Wextra -fsanitize=undefined

CXXFLAGS += -fno-strict-aliasing -Wall -std=c++11 -I. -I../src/ \
    -I../src/bloom/ -I../src/cuckoo/ -I../src/gcs \
    -I../src/gqf/ -I../src/morton/ -I../src/xorfilter -I../src/ribbon \
    $(OPT)

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
        CXXFLAGS += -march=native
else
        CXXFLAGS +=
endif
LDFLAGS = -Wall -Wextra

HEADERS = $(wildcard ../src/*.h \
    ../src/bloom/*.h ../src/cuckoo/*.h ../src/gcs/*.h \
    ../src/gqf/*.h ../src/morton/*.h ../src/xorfilter/*.h ../src/ribbon/*.h \
    ) *.h

.PHONY: all

BINS = bulk-insert-and-query.exe stream.exe construction-failure.exe

all: $(BINS)

clean:
	/bin/rm -f $(BINS)

%.exe: %.cc ${HEADERS}  Makefile
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)
