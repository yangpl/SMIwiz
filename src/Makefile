CC = mpicc
CFLAGS = -g -Wall

fftw3_lib = $(fftw3)/lib
fftw3_inc = $(fftw3)/include

BIN = ../bin
LIB = -L$(fftw3_lib) -lfftw3 -lm -fopenmp -lpthread  -lmpi
INC = -I../include -I$(fftw3_inc)
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)



all: clean SMIwiz

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

SMIwiz:	$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN)/SMIwiz $(OBJ) $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) $(BIN)/SMIwiz


