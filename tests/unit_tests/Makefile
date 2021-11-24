GVERSION    = g++ -std=c++11
COMPFLAGS   = -g -Wall
GTEST_FLAGS = -lgtest -lpthread
MATH        = -lm
VPATH= ../../src

all : tests.x
tests.x: driver.o Vector.o Matrix.o
	$(GVERSION)  $(COMPFLAGS) -o $@ $^ $(MATH) $(GTEST_FLAGS)

clean :
	\rm -f *.o

realclean : clean
	\rm -f tests.x *~ *.a


