CC = gcc
VERSION = -std=c99

MATH = -lm
CFLAGS=-g  -Wall -fopenmp -O3 -march=native
#LIBSLA = -llapack -llapacke
LIBSBLAS= -lopenblas

all : plasma.x
plasma.x: Driver.o get_time.o libplasma.a libcoreblas.a
	$(CC) $(VERSION) $(CFLAGS) -o $@ $^ $(LIBSBLAS) $(LIBSLA) $(MATH)

clean :
	\rm -f *.o plasma.x





