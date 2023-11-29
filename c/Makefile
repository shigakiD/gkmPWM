CC=gcc
CXX=g++
CXXFLAGS=-fopenmp -O3 -fPIC -fno-loop-optimize -fno-aggressive-loop-optimizations
SRC:=$(wildcard *.c)
OBJ:=$(patsubst %.c,%.o,$(SRC))
LIB:=$(wildcard *.a)

.PHONY: all 
all: gkmPWMlasso4

gkmPWMlasso4 : $(OBJ) $(LIB) main.o
	$(CXX) $(CXXFLAGS) $(OBJ) main.o $(LIB) -o gkmPWMlasso4
    
%.o : %.c
	$(CC) -c $(CXXFLAGS) $< -o $@
main.o : main.cpp
	$(CC) -c $(CXXFLAGS) $< -o $@
    
.PHONY: clean
clean:
	$(RM) *.o *~ gkmPWMlasso4
