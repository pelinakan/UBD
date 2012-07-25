# Compiler definition
CC=g++
# Flags for the compiles
CFLAGS=-Wall -O3
UTILS = util/

all: designBarcode findIndexes

designBarcode: hybrid-ss-min.o CtEnergy.o
	$(CC) $(CFLAGS) -fopenmp Main_GenerateandFindIndependentSet_v1.cpp hybrid-ss-min.cpp HybridMin.cpp CtEnergy.cpp energy.cpp util.cpp -o designBarcode -lpthread

findIndexes: findIndexes.o
	$(CC) findIndexes.o -o findIndexes -lpthread -lz

findIndexes.o: findIndexes.cpp
	$(CC) -I$(UTILS) -c -std=c++0x findIndexes.cpp


clean:
	rm -rf *o designBarcode
	rm -rf *o findIndexes