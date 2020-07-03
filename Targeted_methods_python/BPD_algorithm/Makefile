CC=g++
CFLAGS=-std=c++11 -O2
NAME=$(shell uname -s)


ifeq ($(NAME),Darwin)
all: BPD_algorithm_mac
else
all: BPD_algorithm
endif

clean:
	rm -f BP_algorithm

BPD_algorithm: BFalgorithm.cpp
	$(CC) -O3 BFalgorithm.cpp  -o BP_algorithm

BPD_algorithm_mac: BFalgorithm.cpp 
	$(CC) $(CFLAGS) BFalgorithm.cpp  -o BP_algorithm

