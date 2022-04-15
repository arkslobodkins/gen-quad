LPATH = -L /opt/intel/oneapi/mkl/2022.0.2/lib/intel64/
IPATH = -I /opt/intel/oneapi/mkl/2022.0.2/include/
VPATH = ../../src/

CC = gcc
CFLAGS = $(IPATH) -fsanitize=address -fopenmp -O3

all : TestFunction

TestFunction : TestFunction.o GetFunction.o Vector.o Matrix.o     \
	get_time.o Basis.o Quadrature.o Constraints.o Print.o Tensor.o \
	GENERAL_QUADRATURE.o AddDimension.o GeneralGaussTensor.o libgaussq.a
	$(CC) $(CFLAGS) -o $@ $^ $(LPATH) -lmkl_rt -lm

clean:
	rm -f *.o TestFunction
