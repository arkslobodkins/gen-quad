debug = 0
MKL = 1

CXX = g++
CXXFLAGS =-std=c++14 -O3 -fopenmp
LFLAGS = -lm

IPATH = -I ./include/
VPATH = ./src/:./alglib/src/:./include/
ODIR = ./src/objdir/

ifeq ($(debug), 1)
CXXFLAGS += -g -DGEN_QUAD_DEBUG_ON
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

OBJS = $(ODIR)gen_quad.o $(ODIR)nonlinear_solve.o $(ODIR)util.o $(ODIR)domain.o  \
$(ODIR)node_elimination.o $(ODIR)function.o $(ODIR)basis.o $(ODIR)quadrature.o   \
$(ODIR)compute_quad.o $(ODIR)quad_tensor.o $(ODIR)quad_driver.o                  \
$(ODIR)integration.o $(ODIR)ap.o $(ODIR)alglibmisc.o $(ODIR)alglibinternal.o     \
$(ODIR)linalg.o $(ODIR)specialfunctions.o $(ODIR)mesh.o $(ODIR)main.o

DEPS = gen_quad.hpp nonlinear_solve.hpp util.hpp domain.hpp                      \
node_elimination.hpp function.hpp basis.hpp quadrature.hpp                       \
compute_quad.hpp quad_tensor.hpp quad_driver.hpp                                 \
math_func.hpp exceptions.hpp headers.hpp mesh.hpp

all: objdir results quad

$(ODIR)gen_quad.o : gen_quad.cpp gen_quad.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)nonlinear_solve.o : nonlinear_solve.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)util.o : util.cpp util.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)domain.o : domain.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)node_elimination.o : node_elimination.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)function.o : function.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)basis.o : basis.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)quadrature.o : quadrature.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)compute_quad.o : compute_quad.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)quad_tensor.o : quad_tensor.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)quad_driver.o : quad_driver.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)integration.o : integration.cpp integration.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)ap.o : ap.cpp ap.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)alglibinternal.o : alglibinternal.cpp alglibinternal.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)alglibmisc.o : alglibmisc.cpp alglibmisc.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)linalg.o : linalg.cpp linalg.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)specialfunctions.o : specialfunctions.cpp specialfunctions.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)mesh.o : mesh.cpp mesh.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ODIR)main.o: main.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(IPATH) -c -o $@ $<

quad: $(OBJS)
	$(CXX) $(CXXFLAGS) -o quad.x $^ $(LFLAGS)
objdir:
	mkdir -vp 'src/objdir'

results:
	mkdir -vp 'results' 'results/times' 'results/history' 'results/quad_rules' 'results/efficiencies'

.PHONY: clean

clean:
	\rm -f src/objdir/*.o quad.x

distclean:
	\rm -f src/objdir/*.o results/times/* results/history/* results/quad_rules/* results/efficiencies/* quad.x
