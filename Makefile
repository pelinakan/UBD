# Compiler definition
CC=g++
# Flags for the compiles
CFLAGS=-Wall -O2 -ggdb

all: designBarcode

designBarcode: hybrid-ss-min.o CtEnergy.o
	$(CC) -openmp Main_GenerateandFindIndependentSet_v1.cpp hybrid-ss-min.cpp CtEnergy.cpp energy.cpp util.cpp -o designBarcode -lpthread

clean:
	rm -rf *o designBarcode