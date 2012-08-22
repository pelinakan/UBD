# Compiler definition
CC=g++
# Flags for the compiles
CFLAGS=-Wall -O3
UTILS = util/

all: findIndexes

findIndexes: findIndexes.o
	$(CC) findIndexes.o -o findIndexes -lpthread -lz

findIndexes.o: findIndexes.cpp
	$(CC) -I$(UTILS) -c -std=c++0x findIndexes.cpp


clean:
	rm -rf *o findIndexes
