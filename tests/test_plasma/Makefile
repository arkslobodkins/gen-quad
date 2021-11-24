CC = gcc
VERSION = -std=c99

MATH = -lm
CFLAGS=-g  -Wall -fopenmp
LIBSLA = -llapack -llapacke
LIBSBLAS= -lblas

all : plasma.x
plasma.x: Driver.o get_time.o libplasma.a libcoreblas.a
	$(CC) $(VERSION) $(CFLAGS) -o $@ $^ $(LIBSBLAS) $(LIBSLA) $(MATH)




