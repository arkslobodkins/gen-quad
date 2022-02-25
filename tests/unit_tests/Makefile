CVERSION    = g++ -std=c++11
CFLAGS      = -g -Wall
GTEST_FLAGS = -lgtest -lpthread
MATH        = -lm -lopenblas

all : tests.x
tests.x: driver.o Vector.o Matrix.o Tensor.o GENERAL_QUADRATURE.o
	$(CVERSION)  $(CFLAGS) -o $@ $^ $(MATH) $(GTEST_FLAGS)

clean :
	\rm -f *.o tests.x



