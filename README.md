# gen-quad 
**Note**    
*Rules up to dimension 4 for multiple domains are already provided in ComputedRules directory.     
If you want to obtain rules for higher degrees and/or dimensions, follow use instructions below.*  
 
**gen-quad is a package that computes efficient quadrature rules of moderate and high  
accuracy for various polytopes.** Dimensions of interest lie in the range 2-6,  
but one can go to higher dimensions. Theoretically, the algorithm can compute quadrature  
rules of arbitrarily high degree and dimension. In practice, it is limited to available  
memory resources and adequate runtime. Mathematical details are described in our  
paper https://arxiv.org/abs/2207.10737.    
  
**Currently supported domains are:**   
INTERVAL(simply returns Gauss-Legendre rules)  
N-dimensional CUBE  
N-dimensional SIMPLEX(tetrahedron)  
N-dimensional CUBESIMPLEX(tensor product of a cube and a simplex)  
N-dimensional SIMPLEXSIMPLEX(tensor product of two simplexes)  
3-dimensional PYRAMID  
**arbitrary 2-dimensional polytope**
  
All domains (except INTERVAL) are computed using a node elimination algorithm in combination  
with recursive procedures. The node elimination is an iterative technique, which accepts initial  
quadrature rule as a guess and eliminates as many nodes as possible, thereby obtaining more  
efficient quadrature rule. Due to high computational complexity of the algorithm, runtime increases  
rapidly with the increase in degree of accuracy and dimension. To address that, **it has been parallelized using OpenMP.**  
  
Use instructions:   
**gen-quad has no dependencies besides Eigen and ALGLIB, which are included in the directory. Therefore, nothing needs to be installed.
However, if MKL is available on your system, Eigen is linked to MKL library, which speeds up linear algebra routines.
MKL support can be disabled by setting MKL=0 at the beginning of Makefile or passing it from command line: make MKL=0.
Requires C++14 compiler support.**  

To obtain quadrature rule for a particular domain, degree, and dimension,   
appropriate driver routines must be called. Examples are provided in main.cpp.  
**Final quadrature rules are stored in .txt files in results/quad_rules directory.**   

To make an executable, type "make" in the main project directory  
(executable will be produced in the same directory).   
By default "make" uses gcc compiler, but other compilers can be used also.  
For example, use "make compiler=clang" or "make compiler=icpx" for clang and intel compilers, respectively.  
Similarly, "make -j8 compiler=clang" will compile in parallel.  
To compile debug mode: "make -j8 compiler=gcc debug=1". It enables many assertions within gen-quad, as well as Eigen. 
One can also additionally enable address sanitizer by setting "debug=2", but should be used only for testing purposes, 
since it will greatly reduce performance.  
Has been tested with GCC and intel compilers on multiple Linux platforms.    
 

Additional notes:
At each iteration of the node elimination algorithm one needs to solve an underdetermined
system of equations. As a result, there are many possible solutions and small floating
point perturbations may lead to different results over time. Therefore, depending
on the hardware and number of threads used in the program, results for the same problem
may vary slightly from one execution to another. 

   
 
