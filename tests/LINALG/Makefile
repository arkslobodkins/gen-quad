LPATH = -L /opt/intel/oneapi/mkl/2022.0.2/lib/intel64/
IPATH = -I /opt/intel/oneapi/mkl/2022.0.2/include/ -I ../../plasma-17.1/include/
VPATH = ../../src/

CC = gcc
CFLAGS = $(IPATH) -fopenmp -O3 -fsanitize=address

all : TestQR

TestQR : TestQR.o LINALG.o Matrix.o                                     \
	Tensor.o Vector.o Print.o Quadrature.o Basis.o AddDimension.o        \
	GENERAL_QUADRATURE.o GeneralGaussTensor.o Constraints.o get_time.o   \
	libplasma.a libcoreblas.a libgaussq.a
	$(CC) $(CFLAGS) -o $@ $^ $(LPATH) -lmkl_rt -lm

clean:
	rm -f *.o TestQR
