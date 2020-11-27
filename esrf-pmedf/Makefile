# Compile pmedf_readC.cpp into pmedf_readC.oct.

#CC=gcc
#LDFLAGS=-lm

all: pmedf_readC.oct

pmedf_readC.oct: pmedf_readC.cpp
	mkoctfile $< -s -lz
	rm -f pmedf_readC.o

