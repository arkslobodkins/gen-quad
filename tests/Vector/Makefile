LPATH = -L /opt/intel/oneapi/mkl/2022.0.2/lib/intel64/
IPATH = -I /opt/intel/oneapi/mkl/2022.0.2/include/
VPATH = ../../src/

CC = gcc
CFLAGS = $(IPATH) -fsanitize=address

all : TestVector

TestVector : TestVector.o Vector.o
	$(CC) $(CFLAGS) -o $@ $^ $(LPATH) -lmkl_rt -lm

clean:
	rm -f *.o TestVector
