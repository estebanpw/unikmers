CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: unikmers



unikmers: unikmers.c
	$(CXX) $(CFLAGS) alignmentFunctions.c -lm commonFunctions.c -lm unikmers.c -lm -o $(BIN)/unikmers

clean:
	rm -rf $(BIN)/unikmers
