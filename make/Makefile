assertions = 1
OPENMP     = 1
debug      = 1

CXX = g++
CC  = gcc
CXXFLAGS = -O3 -mfma -march=native -isystem../plasma-17.1/include  -fsanitize=address #-fopt-info -Og
CFLAGS   = -O3 -mfma -march=native -isystem../plasma-17.1/include  -fsanitize=address #-fopt-info -Og

CVERSION   = -std=c11
CPPVERSION = -std=c++11

ifeq ($(OPENMP), 1)
CFLAGS   += -fopenmp
CXXFLAGS += -fopenmp
endif

ifeq ($(assertions), 0)
CFLAGS   += -DNDEBUG
CXXFLAGS += -DNDEBUG
endif

ifeq ($(debug), 0)
CFLAGS   += -DQUAD_DEBUG_OFF
CXXFLAGS += -DQUAD_DEBUG_OFF
else
CFLAGS   += -g -Wall -Wextra
CXXFLAGS += -g -Wall -Wextra

endif

LIBSBLAS= -lblas -llapack -llapacke
#alternatively:
#LIBSBLAS = -lopenblas
MATH = -lm
RANLIB = ranlib
VPATH = ../src:../src/Gauss_Lib

OBJECTS = Main.o LeastSquaresNewton.o  GeneralGaussTensor.o GENERAL_QUADRATURE.o GetJacobian.o  \
   GetFunction.o NodeElimination.o Output.o AddDimension.o ComputeDomain.o QuadDriver.o         \
   Constraints.o ConstrainedOptimization.o Vector.o Matrix.o Print.o Quadrature.o               \
	LINALG.o Basis.o Tensor.o get_time.o libgaussq.a libplasma.a libcoreblas.a

all : quadc quadcpp

quadc : quad_c results data
quad_c : $(OBJECTS)
	$(CC) $(CVERSION) $(CFLAGS) -o $@ $^ $(LIBSBLAS) $(MATH)

quadcpp : quad_cpp results data
quad_cpp : $(OBJECTS)
	$(CXX) $(CPPVERSION) $(CXXFLAGS) -o $@ $^ $(LIBSBLAS) $(MATH)


OBJS = GaussCore.o Jacobi.o
libgaussq.a : $(OBJS)
	$(AR) rv libgaussq.a $(OBJS)
	$(RANLIB) libgaussq.a

results:
	mkdir -p ../results
data:
	mkdir -p ../data

clean:
	\rm -f *.o

distclean : clean txtclean
	\rm -f quad_c quad_cpp *~

SUBDIR_ROOTS := ../results ../data
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
PATTERNS := *.txt *.png
TXT_FILES := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(PATTERNS)))
txtclean:
	rm -rf $(TXT_FILES)

