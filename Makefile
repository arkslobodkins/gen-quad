####################################################################
debug = 0
MKL = 1
compiler = gcc

ifeq ($(compiler), gcc)
CXX = g++
CXXFLAGS =-fopenmp -march=native
else ifeq ($(compiler), clang)
CXX = clang++
CXXFLAGS =-fopenmp=libomp -march=native
else ifeq ($(compiler), icpx)
CXX = icpx
CXXFLAGS =-qopenmp -xHost -fp-model=precise
endif


CXXFLAGS +=-std=c++14 -O3
LFLAGS = -lm

IPATH = -I ./include/
VPATH = ./src/:./alglib/src/:./include/
ODIR = ./src/objdir/

ifeq ($(debug), 1)
CXXFLAGS += -g -DGEN_QUAD_DEBUG_ON
ifeq ($(compiler), gcc)
CXXFLAGS += -D_GLIBCXX_DEBUG
endif
endif

ifeq ($(debug), 2)
CXXFLAGS += -g -DGEN_QUAD_DEBUG_ON -fsanitize=undefined -fsanitize=address
endif

####################################################################
ifeq ($(MKL), 1)

ifdef MKLROOT
CXXFLAGS += -DGEN_QUAD_USE_MKL
LFLAGS += -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
#LFLAGS += -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -ldl # flags for gnu threading if intel threading is not available
endif

endif

####################################################################
OBJ = gen_quad.o nonlinear_solve.o util.o domain.o node_elimination.o function.o basis.o quadrature.o compute_quad.o   \
quad_tensor.o quad_driver.o integration.o ap.o alglibmisc.o alglibinternal.o linalg.o specialfunctions.o mesh.o main.o
OBJS = $(addprefix $(ODIR), $(OBJ))

DEPS = gen_quad.hpp nonlinear_solve.hpp util.hpp domain.hpp node_elimination.hpp function.hpp basis.hpp quadrature.hpp \
compute_quad.hpp quad_tensor.hpp quad_driver.hpp math_func.hpp exceptions.hpp headers.hpp mesh.hpp

all: objdir results results_dir quad

$(ODIR)%.o: %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(IPATH) -c $< -o $@

quad: $(OBJS)
	$(CXX) $(CXXFLAGS) -o quad.x $^ $(LFLAGS)
objdir:
	mkdir -vp 'src/objdir'

results:
	mkdir -vp 'results'

results_dir:
	mkdir -vp 'results/times' 'results/history' 'results/quad_rules' 'results/efficiencies'

####################################################################
.PHONY: clean

clean:
	rm -f src/objdir/*.o quad.x

distclean:
	rm -f src/objdir/*.o results/times/* results/history/* results/quad_rules/* results/efficiencies/* quad.x
