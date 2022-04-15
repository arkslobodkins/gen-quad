CC = gcc
VERSION = -std=c99

LPATH = -L /opt/intel/oneapi/mkl/2022.0.2/lib/intel64/
IPATH = -I /opt/intel/oneapi/mkl/2022.0.2/include/

CFLAGS=-g  -Wall -fopenmp -O3 -march=native $(IPATH) $(LPATH)
LFLAGS= -lmkl_rt -lm

all : plasma.x
plasma.x: Driver.o get_time.o libplasma.a libcoreblas.a
	$(CC) $(VERSION) $(CFLAGS) $(LPATH) $(IPATH) -o $@ $^ $(LFLAGS)

clean :
	\rm -f *.o plasma.x





