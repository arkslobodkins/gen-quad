gen-quad is a package that computes efficient quadrature rules of moderate and high  
accuracy for various polytopes. Dimensions of interest lie in the range 2-6,  
but one can go to higher dimensions. Theoretically, the algorithm can compute quadrature  
rules of arbitrarily high degree and dimension. In practice, it is limited to available  
memory resources and adequate runtime.  
  
Currently supported domains are:  
N-dimensional CUBE  
N-dimensional SIMPLEX(tetrahedron)  
N-dimensional CUBESIMPLEX(tensor product of a cube and a simplex)  
N-dimensional SIMPLEXSIMPLEX(tensor product of two simplexes)  
INTERVAL(simply returns Gauss-Legendre rules)  
  
All domains (except INTERVAL) are computed using a node elimination algorithm in combination  
with recursive procedures. The node elimination is an iterative technique, which accepts initial  
quadrature rule as a guess and eliminates as many nodes as possible, thereby obtaining more   
efficient quadrature rule. Due to high computational complexity of the algorithm, runtime increases   
rapidly with the increase in degree of accuracy and dimension. To address that, it has been parallelized 
using OpenMP. 
  
  
Use instructions:  
Requires at least C99 or C++98 or more recent, as well as BLAS and LAPACK libraries. It also
requires availability of LAPACKE interface to support initialization of PLASMA routines. PLASMA 
library has been compiled and statically linked, therefore it does not need to be installed.  
Has been tested with gcc and g++ on multiple Linux platforms.  
Makefile contains instructions for enabling/disabling debugging mode.  
  
When running executable, the following command line arguments should be supplied:  
1.Shape of the domain, which is case-insensitive. For example, either cube, CUBE or cubE are acceptable.  
2.Degree of precision, which should be 1 or greater.  
3.Dimension/dimensions of the domain, depending on domain type.  
  
More specifically, inputs should be passed in the following order:  
cube degree dim                 s.t. dim >= 2.               For example, cube 8 3  
simplex degree dim              s.t  dim>= 2.                For example, simplex 3 4  
cubesimplex degree dim1 dim2    s.t. dim1 >= 1, dim2 >= 2.   For example, cubesimplex 4 1 3  
simplexsimplex degree dim1 dim2 s.t, dim1 >= 2, dim2 >= 2.   For example, simplexsimplex 6 3 2  
  
Final quadrature rules and results are stored in .txt files in results directory.  
