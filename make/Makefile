CC = gcc
VERSION = -std=c99

#uncomment "-DQUAD_DEBUG_OFF"   to disable debug mode
#uncomment "-NDEBUG"            to disable assert statements
CFLAGS=-g -Wall -fopenmp        \
      #-DQUAD_DEBUG_OFF         \
      #-DNDEBUG


LIBSLA = -llapack
LIBSBLAS= -lblas
MATH= -lm
RANLIB = ranlib
VPATH= ../src:../src/Gauss_Lib

all : quad results
quad: Main.o InsertionSort.o BasisIntegrals.o LeastSquaresNewton.o GeneralGaussTensor.o                \
      GetJacobian.o GetFunction.o NodeElimination.o BasisFunctions.o BasisIndices.o Output.o           \
      AddDimension.o ComputeDomain.o QuadDriver.o Constraints.o Vector.o Matrix.o glist.o Print.o      \
      Quadrature.o LINALG.o Polytope.o libgaussq.a
	$(CC) $(VERSION) $(CFLAGS) -o $@ $^ $(LIBSBLAS) $(LIBSLA)  $(MATH)

OBJS =  GaussCore.o Jacobi.o
libgaussq.a: $(OBJS)
	$(AR) rv libgaussq.a $(OBJS)
	$(RANLIB) libgaussq.a

results:
	mkdir -p ../results

clean :
	\rm -f *.o

distclean : clean txtclean
	\rm -f quad *~ *.a

SUBDIR_ROOTS := ../results
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
PATTERNS := *.txt
TXT_FILES := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(PATTERNS)))
txtclean:
	rm -rf $(TXT_FILES)

