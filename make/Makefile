CC = gcc -std=c99 -Wall

#comment   "-DCONSTR_OPT"       to disable constrained optimization\
#uncomment "-DQUAD_DEBUG_OFF"   to disable debug mode
#uncomment "-NDEBUG"            to disable assert statements
CFLAGS=-g                 \
		-DCONSTR_OPT        \
		#-DQUAD_DEBUG_OFF   \
		#-DNDEBUG


LIBSLA = -llapack
LIBSBLAS= -lblas
MATH= -lm
RANLIB = ranlib
VPATH= ../src:../src/Gauss_Lib:../results
OBJS =  GaussCore.o Jacobi.o

all : quad
quad: Main.o InsertionSort.o BasisIntegrals.o LeastSquaresNewton.o GeneralGaussTensor.o                \
	   GetJacobian.o GetFunction.o NodeElimination.o BasisFunctions.o BasisIndices.o Output.o           \
		AddDimension.o ComputeDomain.o QuadDriver.o Constraints.o Vector.o Matrix.o Print.o Quadrature.o \
	   LINALG.o libgaussq.a

	$(CC)  $(CFLAGS) -o $@ $^ $(LIBSBLAS) $(LIBSLA)  $(MATH)

libgaussq.a: $(OBJS)
	$(AR) rv libgaussq.a $(OBJS)
	$(RANLIB) libgaussq.a

clean :
	\rm -f *.o

realclean : clean
	\rm -f quad *~ *.a

SUBDIR_ROOTS := ../results
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
PATTERNS := *.txt
TXT_FILES := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(PATTERNS)))
txtclean:
	rm -rf $(TXT_FILES)

