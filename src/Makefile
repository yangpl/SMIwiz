CC=mpicc # Intel compiler
CFLAGS= -g -Wall 

#assume your FFTW has been installed and the path is setup as variable $(fftw)
fftw3_lib = $(fftw)/lib
fftw3_inc = $(fftw)/include

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


