CC = gcc
VERSION = -std=c11

#uncomment "-DQUAD_DEBUG_OFF"   to disable debug mode
#uncomment "-NDEBUG"            to disable assert statements
#comment out -fopenmp           to disable OpenMP parallelism
CFLAGS= -g -Wall -fopenmp -O3 -march=x86-64
      #-fsanitize=address              \
      #-DQUAD_DEBUG_OFF                \
      #-DNDEBUG                        \
      #-fopt-info                      \

# alternatively, LIBSBLAS= -lopenblas -llapacke if openblas is available
LIBSBLAS= -lblas -llapack -llapacke
MATH= -lm
RANLIB = ranlib
VPATH= ../src:../src/Gauss_Lib

all : quad results data
quad: Main.o LeastSquaresNewton.o  GeneralGaussTensor.o GENERAL_QUADRATURE.o GetJacobian.o  \
	GetFunction.o NodeElimination.o Output.o AddDimension.o ComputeDomain.o QuadDriver.o     \
	Constraints.o ConstrainedOptimization.o Vector.o Matrix.o Print.o Quadrature.o           \
	LINALG.o Basis.o get_time.o libgaussq.a libplasma.a libcoreblas.a
	$(CC) $(VERSION) $(CFLAGS) -o $@ $^ $(LIBSBLAS) $(MATH)

OBJS =  GaussCore.o Jacobi.o
libgaussq.a: $(OBJS)
	$(AR) rv libgaussq.a $(OBJS)
	$(RANLIB) libgaussq.a

results:
	mkdir -p ../results
data:
	mkdir -p ../data

clean :
	\rm -f *.o

distclean : clean txtclean
	\rm -f quad *~

SUBDIR_ROOTS := ../results ../data
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
PATTERNS := *.txt *.png
TXT_FILES := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(PATTERNS)))
txtclean:
	rm -rf $(TXT_FILES)

