assertions = 1
OPENMP     = 1
debug      = 1

CXX = g++
CC  = gcc
CXXFLAGS = -g -Wall -O3 -march=native -I../plasma-17.1/include #-fsanitize=address -fopt-info
CFLAGS   = -g -Wall -O3 -march=native -I../plasma-17.1/include #-fsanitize=address #-fopt-info

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
endif

#alternatively: LIBSBLAS= -lblas -llapack -llapacke
LIBSBLAS = -lopenblas
MATH = -lm
RANLIB = ranlib
VPATH = ../src:../src/Gauss_Lib

OBJECTS = Main.o LeastSquaresNewton.o  GeneralGaussTensor.o GENERAL_QUADRATURE.o GetJacobian.o  \
	GetFunction.o NodeElimination.o Output.o AddDimension.o ComputeDomain.o QuadDriver.o         \
	Constraints.o ConstrainedOptimization.o Vector.o Matrix.o Print.o Quadrature.o               \
	LINALG.o Basis.o get_time.o libgaussq.a libplasma.a libcoreblas.a

all : quad_c quad_cpp

quad_c : quad_c.x results data
quad_c.x : $(OBJECTS)
	$(CC) $(CVERSION) $(CFLAGS) -o $@ $^ $(LIBSBLAS) $(MATH)

quad_cpp : quad_cpp.x results data
quad_cpp.x : $(OBJECTS)
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
	\rm -f quad_c.x quad_cpp.x*~

SUBDIR_ROOTS := ../results ../data
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
PATTERNS := *.txt *.png
TXT_FILES := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(PATTERNS)))
txtclean:
	rm -rf $(TXT_FILES)

