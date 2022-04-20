OPENMP     = 1
assertions = 1
debug      = 1

LPATH = -L /opt/intel/oneapi/mkl/2022.0.2/lib/intel64/
IPATH = -I /opt/intel/oneapi/mkl/2022.0.2/include/ -I ../plasma-17.1/include/

CC  = gcc
CVERSION = -std=c11
LFLAGS = -lm -lmkl_rt
CFLAGS = -O3 -march=native $(IPATH) #-fsanitize=address #-fopt-info -Og

ifeq ($(OPENMP), 1)
CFLAGS += -fopenmp
else
CFLAGS += -Wno-unknown-pragmas
endif

ifeq ($(assertions), 0)
CFLAGS += -DNDEBUG
endif

ifeq ($(debug), 0)
CFLAGS   += -DQUAD_DEBUG_OFF
else
CFLAGS   += -g -Wall -Wextra
endif

RANLIB = ranlib
VPATH = ../src:../src/Gauss_Lib

OBJECTS = Main.o NonlinearSolve.o  GeneralGaussTensor.o GENERAL_QUADRATURE.o GetJacobian.o      \
   GetFunction.o NodeElimination.o Output.o AddDimension.o ComputeDomain.o QuadDriver.o         \
   Constraints.o ConstrainedOptimization.o Vector.o Matrix.o Print.o Quadrature.o               \
	LINALG.o Basis.o Tensor.o get_time.o libgaussq.a libplasma.a libcoreblas.a

all : quad_c

quad_c : quad results
quad : $(OBJECTS)
	$(CC) $(CVERSION) $(CFLAGS) -o $@ $^ $(LFLAGS) $(LPATH)

OBJS = GaussCore.o Jacobi.o
libgaussq.a : $(OBJS)
	$(AR) rv libgaussq.a $(OBJS)
	$(RANLIB) libgaussq.a

results:
	mkdir -p ../results

clean:
	\rm -f *.o quad quad_cpp

distclean : clean txtclean
	\rm -f quad quad_cpp *~

SUBDIR_ROOTS := ../results
DIRS := . $(shell find $(SUBDIR_ROOTS) -type d)
PATTERNS := *.txt *.png
TXT_FILES := $(foreach DIR,$(DIRS),$(addprefix $(DIR)/,$(PATTERNS)))
txtclean:
	rm -rf $(TXT_FILES)

