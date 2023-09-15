CC=mpicc # Intel compiler
CFLAGS= -g -Wall 


BIN = ../bin
LIB =  -lm -fopenmp -lpthread -lfftw3 -lmpi
INC = -I../include
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


