#Get OS
UNAME=$(shell uname)
# Compiler definition
CC=g++

# Flags for the compiles
ifeq ($(UNAME), Linux)
# do something Linux-y
CFLAGS= -O3 -Wno-unused-result -Wno-write-strings
else
CFLAGS= -O3 -Wno-write-strings
endif

UTILS = util/

all: designBarcode findIndexes

designBarcode: hybrid-ss-min.o CtEnergy.o
	$(CC) $(CFLAGS) -fopenmp Main_GenerateandFindIndependentSet_v1.cpp hybrid-ss-min.cpp HybridMin.cpp CtEnergy.cpp energy.cpp util.cpp -o designBarcode -lpthread

findIndexes: findIndexes.o
	$(CC) $(CFLAGS) findIndexes.o -o findIndexes -lpthread -lz

ifeq ($(UNAME), Linux)
findIndexes.o: findIndexes.cpp
	$(CC) -I$(UTILS) -c -std=c++0x findIndexes.cpp
else
findIndexes.o: findIndexes.cpp
	$(CC) -I$(UTILS) -c findIndexes.cpp
endif

clean:
	rm -rf *o designBarcode
	rm -rf *o findIndexes
