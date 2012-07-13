# Compiler definition
CC=g++
# Flags for the compiles
CFLAGS=-Wall -O3

all: designBarcode

designBarcode: hybrid-ss-min.o CtEnergy.o
	$(CC) $(CFLAGS) -fopenmp Main_GenerateandFindIndependentSet_v1.cpp hybrid-ss-min.cpp HybridMin.cpp CtEnergy.cpp energy.cpp util.cpp -o designBarcode -lpthread

clean:
	rm -rf *o designBarcode