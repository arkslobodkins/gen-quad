/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "BasisFunctions.h"
#include "BasisIntegrals.h"
#include "BasisIndices.h"
#include "AddDimension.h"
#include "GeneralGaussTensor.h"
#include "Gauss_Lib/Jacobi.h"
#include "Quadrature.h"
#include "Print.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

static void BasisSimplexPolyhedralOne(int dim, int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);
static void BasisSimplexPolyhedralTwo(int dim, int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi);

static void BasisPrimeSimplexPolyhedralOne(int dim, int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);
static void BasisPrimeSimplexPolyhedralTwo(int dim, int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime);

typedef void(*BasisFunc)(int *dims, int deg, const INT_8 *basis_id, const double *x,  double *phi);

ATTR_UNUSED static int finite_difference_test(int dim, int *dims, int deg, const INT_8 *basis_id,
                                              const double *x, const double *phi_prime, BasisFunc func);


void LegendrePoly(int order, double x, double *p, double *dp)
{
   int k;
   double fac1 = 0.0, fac2 = 0.0;

   p[0] = 1.0;
   dp[0] = 0.0;
   p[1] = x;
   dp[1] = 1.0;

   for( k = 1; k < order-1; ++k)
   {
      fac1 = (2.0*k+1.0)/(k+1.0);
      fac2 = k/(k+1.0);
      p[k+1] = fac1*x*p[k] - fac2*p[k-1];
      dp[k+1] = fac1*(p[k] + x*dp[k]) - fac2*dp[k-1];
   }
}

#define ADD(x, y) (x)+(y)
#define SUB(x, y) (x)-(y)
void JacobiPoly(int order, double x, double alpha, double beta, double *p)
{
   double fac1, fac2, fac3;
   fac1 = fac2 = fac3 = 0.0;

   // compute some constants to be used inside a for loop to avoid redundant computations,
   // whereas some are computed repeatedly inide the loop(tuned)
   double a_plus_b               = ADD(alpha, beta);
   double a_plus_b_minus_1       = SUB(a_plus_b, 1.0);
   double a_plus_b_minus_2       = SUB(a_plus_b, 2.0);

   p[0] = 1.0;
   if(order > 1)
      p[1] = 0.5*(x-1.0) * (alpha+beta+2.0) + alpha+1.0;

   for (int k = 1; k < order-1; ++k)
   {
      double two_n_plus_a_plus_b = ADD(2.0*(k+1.0), a_plus_b);
      fac1 = (2.0*(k+1)+a_plus_b_minus_1) * (
          (2.0*(k+1.0)+alpha+beta) * (2.0*(k+1.0)+alpha+beta-2.0) * x + (alpha*alpha-beta*beta)
          );
      fac2 = -2.0 * (k+1.0+alpha-1.0) * (k+1.0+beta-1.0) * two_n_plus_a_plus_b;
      fac3 = 1.0 / (2*(k+1.0)*(k+1.0+a_plus_b) * (2*(k+1.0)+a_plus_b_minus_2));
      p[k+1] = fac3 * (fac1*p[k] + fac2*p[k-1]);
   }
}
#undef ADD
#undef SUB

void JacobiPolyPrime(int order, double x, double alpha, double beta, double *dp)
{
   double p[order];
   double a_plus_b_plus1 = 1.0+alpha+beta;

   JacobiPoly(order, x, alpha+1.0, beta+1.0, p);
   dp[0] = 0.0;
   dp[1] = 0.5 * (alpha+beta+2.0);

   for(int k = 2; k < order; ++k)
      dp[k] = 0.5*(a_plus_b_plus1+k) * p[k-1];

}


// Generates orthogonal polynomial basis for the unit cube
void CubeBasisFuncs(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi)
{
   int dim = dims[0];
   assert(dim >= 1);
   int degPlus1 = deg+1;
   int num_funcs = BasisSize(deg, dim);

   double legendre[degPlus1*dim];
   double dxlegendre[degPlus1*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < num_funcs; ++k) phi[k] = 1.0;

   for(int k = 0; k < num_funcs; ++k)
      for(int d = 0; d < dim; ++d)
         phi[k] *= legendre[basis_id[k*dim+d] + degPlus1*d];

}// BasisCube


// Generates derivatives of orthogonal polynomial basis for the unit cube
void CubeBasisFuncsDer(int *dims, int deg, const INT_8 *basis_id, const double* x, double *phiPrime)
{
   assert(dims[0] >= 1);
   int dim      = dims[0];
   int degPlus1    = deg+1;
   int num_funcs = BasisSize(deg, dim);

   double legendre[degPlus1*dim];
   double dxlegendre[degPlus1*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < num_funcs*dim; ++k) phiPrime[k] = 1.0;

   for(int d = 0; d < dim; ++d)
   {

      // dimension < d
      for(int j = 0; j < d; ++j)
      {
         double *legendre_ptr = &legendre[degPlus1*j];
         for(int k = 0; k < num_funcs; ++k)
            phiPrime[d*num_funcs+k] *= legendre_ptr[basis_id[dim*k+j]];
      }
      // dimension = d
      double *dxlegendre_ptr = &dxlegendre[degPlus1*d];
      for(int k = 0; k < num_funcs; ++k)
         phiPrime[d*num_funcs+k] *= 2.0 * dxlegendre_ptr[basis_id[dim*k+d]];

      // dimension > d
      for(int j = d+1; j < dim; ++j)
      {
         double *legendre_ptr = &legendre[degPlus1*j];
         for(int k = 0; k < num_funcs; ++k)
            phiPrime[d*num_funcs+k] *= legendre_ptr[basis_id[dim*k+j]];
      }

   }

#ifdef QUAD_DEBUG_ON
   BasisFunc phi_func = &CubeBasisFuncs;
   int flag = finite_difference_test(dim, dims, deg, basis_id, x, phiPrime, phi_func);
   if(flag == 1)
   {
      PRINT_ERR("failed finite_difference_test", __LINE__, __FILE__);
      for(int d = 0; d < dim; ++d)
         printf("x[%i] = %lf  ", d, x[d]);
      printf("\n");
   }
#endif

}// BasisPrimeCube


// Generates orthogonal polynomial basis for the unit simplex
void BasisSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi)
{
   assert(dims[0] >= 2);
   int dim      = dims[0];
   int degPlus1    = deg+1;
   int num_funcs = BasisSize(deg, dim);

   if(dim == 2)
   {
      double legendre[degPlus1];
      double dxlegendre[degPlus1];
      double jacobi[degPlus1*degPlus1];

      LegendrePoly(degPlus1, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
      for(int j = 0; j < degPlus1; ++j) JacobiPoly(degPlus1, 1.0-2*x[0], 2.0*j+1, 0.0, &jacobi[j*degPlus1]);
      for(int k = 0; k < num_funcs; ++k) phi[k] = 1.0;

      int r, n;
      for(int k = 0; k < num_funcs; ++k)
      {
         r = basis_id[dim*k];
         n = basis_id[dim*k+1];
         double xFactor = 1.0; for(int i = 0; i < r; ++i) xFactor *= x[0];
         phi[k] *= xFactor * jacobi[r*degPlus1+n] * legendre[r];
      }

      return;
   }

   // dim >=  3
   int i, j, k, d;
   double legendre[degPlus1];
   double dxlegendre[degPlus1];

   double *jacobi = (double *) malloc( SQUARE(degPlus1)*(dim-1)*size_double );
   LegendrePoly(degPlus1, (2.0*x[dim-1]-x[dim-2])/x[dim-2], legendre, dxlegendre);
   for(k = 0; k < num_funcs; ++k) phi[k] = 1.0;

   double polyX[dim-1];
   for(d = 0; d < dim-1; ++d) polyX[d] = 1.0;
   for(d = 0; d < dim-2; ++d) polyX[d] = x[dim-d-2]/x[dim-d-3];
   polyX[dim-2] = x[0];

   for(d = 1; d < dim; ++d)
   {
      int dimCur = d-1;
      double jacobiX_dimCur = 1.0;
      if(d < dim-1) {
         jacobiX_dimCur /= x[dim-d-2];
         jacobiX_dimCur = 1.0-2.0*x[dim-d-1] * jacobiX_dimCur;
      }
      else if(d == dim-1) {
         jacobiX_dimCur *= x[0];
         jacobiX_dimCur = 1.0-2.0*jacobiX_dimCur;
      }

      double alpha = 0.0;
      for(j = 0; j < degPlus1; ++j) {
         alpha = 2*j+d;
         JacobiPoly(degPlus1, jacobiX_dimCur, alpha, 0.0, &jacobi[degPlus1*dimCur*deg+degPlus1*j]);
      }
   }


   for(d = 1; d < dim; ++d)
   {
      double xTemp = polyX[d-1];

      for(k = 0; k < num_funcs; ++k)
      {
         const INT_8 *index = &basis_id[k*dim];
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= xTemp;
         double phiTemp = 1.0;
         phiTemp *= jacobi[*(index+d) + degPlus1*(d-1)*deg + degPlus1*xPower] * xFactor;
         phi[k] *= phiTemp;
      }

   }
   for(k = 0; k < num_funcs; ++k) {
      const INT_8 *index = &basis_id[k*dim];
      phi[k] *= legendre[*index];
   }

   free(jacobi);

   return;
}// BasisSimplex


// Generates derivatives of orthogonal polynomial basis for the unit simplex
void BasisPrimeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime)
{
   int dim = dims[0];
   int num_funcs = BasisSize(deg, dim);
   double h = POW_DOUBLE(10, -5)*5.0;

   double *phi_backw1 = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_backw2 = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_forw1  = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_forw2  = (double *)malloc( SIZE_DOUBLE(num_funcs) );

   double x_backw1[dim];
   double x_backw2[dim];
   double x_forw1[dim];
   double x_forw2[dim];

   for(int d = 0; d < dim; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t)
      {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t] = x[d_t];
         x_forw2[d_t] = x[d_t];
      }
      x_backw1[d] = x_backw1[d] - h;
      x_backw2[d] = x_backw2[d] - 2.0*h;
      x_forw1[d] = x_forw1[d] + h;
      x_forw2[d] = x_forw2[d] + 2.0*h;

      BasisSimplex(dims, deg, basis_id, x_backw1, phi_backw1);
      BasisSimplex(dims, deg, basis_id, x_backw2, phi_backw2);
      BasisSimplex(dims, deg, basis_id, x_forw1, phi_forw1);
      BasisSimplex(dims, deg, basis_id, x_forw2, phi_forw2);

      for(int k = 0; k < num_funcs; ++k)
         phiPrime[d*num_funcs+k] = ( 1.0/12.0*phi_backw2[k] - 2.0/3.0*phi_backw1[k] + 2.0/3.0*phi_forw1[k]-1.0/12.0*phi_forw2[k] ) / h;
   }

   free(phi_backw1);
   free(phi_backw2);
   free(phi_forw1);
   free(phi_forw2);
}// BasisPrimeSimplex


// Generates orthogonal polynomial basis for the unit CUBESIMPLEX
void BasisCubeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi)
{
   assert(dims[0] >= 1 && dims[1] >= 2);
   int dim1 = dims[0];
   int dim2 = dims[1];
   int dim  = dim1+dim2;

   int degPlus1 = deg+1;
   int num_funcs = BasisSize(deg, dim);

   double legendre[degPlus1*dim1];
   double dxlegendre[degPlus1*dim1];

   for(int d = 0; d < dim1; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < num_funcs; ++k) phi[k] = 1.0;
   for(int k = 0; k < num_funcs; ++k)
      for(int d = 0; d < dim1; ++d)
         phi[k] *= legendre[basis_id[k*dim+d]+degPlus1*d];

   double phiSimplex[num_funcs];
   BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x, phiSimplex);
   for(int k = 0; k < num_funcs; ++k) phi[k] *= phiSimplex[k];

}// BasisCubeSimplex


// Generates derivatives of orthogonal polynomial basis for the unit CUBESIMPLEX
void BasisPrimeCubeSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime)
{
   assert(dims[0] >= 1 && dims[1] >= 2);
   int dim1     = dims[0];
   int dim2     = dims[1];
   int dim      = dim1+dim2;
   int num_funcs = BasisSize(deg, dim);
   int degPlus1    = deg+1;

   double legendre[degPlus1*dim1];
   double dxlegendre[degPlus1*dim1];

   for(int d = 0; d < dim1; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < num_funcs*dim; ++k) phiPrime[k] = 1.0;

   for(int d = 0; d < dim; ++d)
   {
      for(int j = 0; j < dim1; ++j)
      {
         if(j != d)
            for(int k = 0; k < num_funcs; ++k)
               phiPrime[k+d*num_funcs] *= legendre[basis_id[k*dim+j]+degPlus1*j];

         if(j == d)
            for(int k = 0; k < num_funcs; ++k)
               phiPrime[k+d*num_funcs] *= 2.0 * dxlegendre[basis_id[k*dim+j]+degPlus1*j];
      }
   }

   double h = POW_DOUBLE(10, -5)*5.0;

   double *phi_backw1 = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_backw2 = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_forw1  = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_forw2  = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_polyhedral  = (double *)malloc( SIZE_DOUBLE(num_funcs) );

   double x_backw1[dim];
   double x_backw2[dim];
   double x_forw1[dim];
   double x_forw2[dim];

   BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x, phi_polyhedral);
   for(int d = 0; d < dim1; ++d) {
      for(int k = 0; k < num_funcs; ++k)
         phiPrime[d*num_funcs+k] *= phi_polyhedral[k];
   }

   for(int d = 0; d < dim2; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t) {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t] = x[d_t];
         x_forw2[d_t] = x[d_t];
      }
      x_backw2[d+dim1] = x_backw2[d+dim1] - 2.0*h;
      x_backw1[d+dim1] = x_backw1[d+dim1] - h;
      x_forw1[d+dim1] = x_forw1[d+dim1] + h;
      x_forw2[d+dim1] = x_forw2[d+dim1] + 2.0*h;

      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x_backw1, phi_backw1);
      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x_backw2, phi_backw2);
      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x_forw1, phi_forw1);
      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x_forw2, phi_forw2);

      for(int k = 0; k < num_funcs; ++k)
         phiPrime[(d+dim1)*num_funcs+k] *= ( 1.0/12.0*phi_backw2[k] - 2.0/3.0*phi_backw1[k] + 2.0/3.0*phi_forw1[k]-1.0/12.0*phi_forw2[k] ) / h;
   }

   free(phi_backw1);
   free(phi_backw2);
   free(phi_forw1);
   free(phi_forw2);
   free(phi_polyhedral);

}// BasisPrimeCubeSimplex


// Generates orthogonal polynomial basis for the unit SIMPLEXSIMPLEX
void BasisSimplexSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi)
{
   assert(dims[0] >= 2 && dims[1] >= 2);
   int dim = dims[0]+dims[1];
   int num_funcs = BasisSize(deg, dim);
   double *phiSimplex = (double *)malloc( SIZE_DOUBLE(num_funcs) );

   for(int k = 0; k < num_funcs; ++k) phi[k] = 1.0;

   BasisSimplexPolyhedralOne(dim, dims, deg, basis_id, x, phiSimplex);
   for(int k = 0; k < num_funcs; ++k)
      phi[k] *= phiSimplex[k];

   BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x, phiSimplex);
   for(int k = 0; k < num_funcs; ++k)
      phi[k] *= phiSimplex[k];

   free(phiSimplex);

} // BasisSimplexSimplex


// Generates derivatives of orthogonal polynomial basis for the unit SIMPLEXSIMPLEX
void BasisPrimeSimplexSimplex(int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime)
{
   assert(dims[0] >= 2 && dims[1] >= 2);

   int dim1     = dims[0];
   int dim2     = dims[1];
   int dim      = dim1+dim2;
   int num_funcs = BasisSize(deg, dim);

   for(int k = 0; k < num_funcs*dim; ++k) phiPrime[k] = 1.0;

   double h = POW_DOUBLE(10, -5)*5.0;
   double *phi_backw1      = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_backw2      = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_forw1       = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_forw2       = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_polyhedral  = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double x_backw1[dim];
   double x_backw2[dim];
   double x_forw1[dim];
   double x_forw2[dim];

   BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x, phi_polyhedral);
   for(int d = 0; d < dim1; ++d) {
      for(int k = 0; k < num_funcs; ++k)
         phiPrime[d*num_funcs+k] *= phi_polyhedral[k];
   }
   BasisSimplexPolyhedralOne(dim, dims, deg, basis_id, x, phi_polyhedral);
   for(int d = 0; d < dim2; ++d) {
      for(int k = 0; k < num_funcs; ++k)
         phiPrime[(d+dim1)*num_funcs+k] *= phi_polyhedral[k];
   }

   for(int d = 0; d < dim2; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t) {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t] = x[d_t];
         x_forw2[d_t] = x[d_t];
      }
      x_backw2[d+dim1] = x_backw2[d+dim1] - 2.0*h;
      x_backw1[d+dim1] = x_backw1[d+dim1] - h;
      x_forw1[d+dim1] = x_forw1[d+dim1] + h;
      x_forw2[d+dim1] = x_forw2[d+dim1] + 2.0*h;

      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x_backw1, phi_backw1);
      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x_backw2, phi_backw2);
      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x_forw1, phi_forw1);
      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x_forw2, phi_forw2);

      for(int k = 0; k < num_funcs; ++k)
         phiPrime[(d+dim1)*num_funcs+k] *= ( 1.0/12.0*phi_backw2[k] - 2.0/3.0*phi_backw1[k] + 2.0/3.0*phi_forw1[k]-1.0/12.0*phi_forw2[k] ) / h;
   }

   for(int d = 0; d < dim1; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t) {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t] = x[d_t];
         x_forw2[d_t] = x[d_t];
      }
      x_backw2[d] = x_backw2[d] - 2.0*h;
      x_backw1[d] = x_backw1[d] - h;
      x_forw1[d] = x_forw1[d] + h;
      x_forw2[d] = x_forw2[d] + 2.0*h;

      BasisSimplexPolyhedralOne(dim, dims, deg, basis_id, x_backw1, phi_backw1);
      BasisSimplexPolyhedralOne(dim, dims, deg, basis_id, x_backw2, phi_backw2);
      BasisSimplexPolyhedralOne(dim, dims, deg, basis_id, x_forw1, phi_forw1);
      BasisSimplexPolyhedralOne(dim, dims, deg, basis_id, x_forw2, phi_forw2);

      for(int k = 0; k < num_funcs; ++k)
         phiPrime[d*num_funcs+k] *= ( 1.0/12.0*phi_backw2[k] - 2.0/3.0*phi_backw1[k] + 2.0/3.0*phi_forw1[k]-1.0/12.0*phi_forw2[k] ) / h;
   }

   free(phi_backw1);
   free(phi_backw2);
   free(phi_forw1);
   free(phi_forw2);
   free(phi_polyhedral);
}


// Generates orthogonal polynomial basis for the first [0, dim1-1] coordinates
// over simplex in dim-dimensional space.
static void BasisSimplexPolyhedralOne(int dim, int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi)
{
   assert(dims[0] >= 2);

   int dim1 = dims[0];
   int num_funcs = BasisSize(deg, dim);
   int degPlus1 = deg+1;

   if(dim1 == 2)
   {
      double legendre[degPlus1];
      double dxlegendre[degPlus1];
      double jacobi[SQUARE(degPlus1)];

      LegendrePoly(degPlus1, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
      for(int j = 0; j < degPlus1; ++j) JacobiPoly(degPlus1, 1.0-2*x[0], 2.0*j+1, 0.0, &jacobi[j*degPlus1]);
      for(int k = 0; k < num_funcs; ++k) phi[k] = 1.0;

      int r, n;
      for(int k = 0; k < num_funcs; ++k)
      {
         r = basis_id[dim*k];
         n = basis_id[dim*k+1];
         double xFactor = 1.0; for(int i = 0; i < r; ++i) xFactor *= x[0];
         phi[k] *= xFactor * jacobi[r*degPlus1+n] * legendre[r];
      }
      return;
   }


   // dim1 > 3
   int i, j, k, d;
   double legendre[degPlus1];
   double dxlegendre[degPlus1];

   double *jacobi = (double *) malloc( SQUARE(degPlus1)*(dim1-1)*size_double );
   LegendrePoly(degPlus1, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], legendre, dxlegendre);
   for(k = 0; k < num_funcs; ++k) phi[k] = 1.0;

   double polyX[dim1-1];
   for(d = 0; d < dim1-1; ++d) polyX[d] = 1.0;
   for(d = 0; d < dim1-2; ++d) polyX[d] = x[dim1-d-2]/x[dim1-d-3];
   polyX[dim1-2] = x[0];

   double jacobiX[dim1-1]; for(d = 0; d < dim1-1; ++d) jacobiX[d] = 1.0;

   for(d = 1; d < dim1; ++d)
   {
      int dimCur = d-1;
      if(d < dim1-1)
      {
         jacobiX[dimCur] /= x[dim1-d-2];
         jacobiX[dimCur] = 1.0-2.0*x[dim1-d-1] * jacobiX[dimCur];
      }
      if(d == dim1-1)
      {
         jacobiX[dimCur] *= x[0];
         jacobiX[dimCur] = 1.0-2.0*jacobiX[dimCur];
      }

      double alpha = 0.0;
      for(j = 0; j < deg+1; ++j)
      {
         alpha = 2*j+d;
         JacobiPoly(degPlus1, jacobiX[dimCur], alpha, 0.0, &jacobi[degPlus1*(d-1)*deg+degPlus1*j]);
      }
   }

   for(k = 0; k < num_funcs; ++k)
   {
      const INT_8 *index = &basis_id[k*dim];
      double phiTemp = 1.0;
      for(d = 1; d < dim1; ++d)
      {
         double xTemp = polyX[d-1];
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= xTemp;
         phiTemp *= jacobi[*(index+d) + degPlus1*(d-1)*deg + degPlus1*xPower] * xFactor;
      }
      phi[k] *= phiTemp*legendre[*index];
   }

   free(jacobi);
   return;
}// BasisSimplexPolyhedralOne


// Generates orthogonal polynomial basis for coordinates between [dim1, dim2-1]
// over a simplex in dim-dimensional space.
static void BasisSimplexPolyhedralTwo(int dim, int *dims, int deg, const INT_8 *basis_id, const double *x, double *phi)
{
   assert(dims[1] >= 2);
   int dim1 = dims[0];
   int dim2 = dims[1];
   int dimTwo = dim1+dim2;
   int degPlus1 = deg+1;
   int num_funcs = BasisSize(deg, dim);

   if(dim2 == 2)
   {
      double legendre[degPlus1];
      double dxlegendre[degPlus1];
      double jacobi[SQUARE(degPlus1)];

      LegendrePoly(degPlus1, (2.0*x[dim1+1]-x[dim1])/x[dim1], legendre, dxlegendre);
      for(int j = 0; j < degPlus1; ++j) JacobiPoly(degPlus1, 1.0-2*x[dim1], 2.0*j+1, 0.0, &jacobi[j*degPlus1]);
      for(int k = 0; k < num_funcs; ++k) phi[k] = 1.0;

      int r, n;
      for(int k = 0; k < num_funcs; ++k)
      {
         r = basis_id[dim*k+dim1];
         n = basis_id[dim*k+dim1+1];
         double xFactor = 1.0; for(int i = 0; i < r; ++i) xFactor *= x[dim1];
         phi[k] *= xFactor*jacobi[r*degPlus1+n] * legendre[r];
      }
      return;
   }


   // dim1 > 3
   int i, j, k, d;
   double legendre[degPlus1];
   double dxlegendre[degPlus1];

   double *jacobi = (double *) malloc( SQUARE(degPlus1)*(dim2-1)*size_double );
   LegendrePoly(degPlus1, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], legendre, dxlegendre);
   for(k = 0; k < num_funcs; ++k) phi[k] = 1.0;

   double polyX[dim2-1];
   for(d = 0; d < dim2-1; ++d) polyX[d] = 1.0;
   for(d = 0; d < dim2-2; ++d) polyX[d] = x[dimTwo-d-2]/x[dimTwo-d-3];
   polyX[dim2-2] = x[dim1];

   double jacobiX[dim2-1]; for(d = 0; d < dim2-1; ++d) jacobiX[d] = 1.0;

   for(d = 1; d < dim2; ++d)
   {
      int dimCur = d-1;
      if(d < dim2-1)
      {
         jacobiX[dimCur] /= x[dimTwo-d-2];
         jacobiX[dimCur] = 1.0-2.0*x[dimTwo-d-1] * jacobiX[dimCur];
      }
      if(d == dim2-1)
      {
         jacobiX[dimCur] *= x[dim1];
         jacobiX[dimCur] = 1.0-2.0*jacobiX[dimCur];
      }

      double alpha = 0.0;
      for(j = 0; j < deg+1; ++j)
      {
         alpha = 2*j+d;
         JacobiPoly(degPlus1, jacobiX[dimCur], alpha, 0.0, &jacobi[degPlus1*(d-1)*deg+degPlus1*j]);
      }
   }

   for(k = 0; k < num_funcs; ++k)
   {
      const INT_8 *index = &basis_id[k*dim+dim1];
      double phiTemp = 1.0;
      for(d = 1; d < dim2; ++d)
      {
         double xTemp = polyX[d-1];
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= xTemp;
         phiTemp *= jacobi[*(index+d) + degPlus1*(d-1)*deg + degPlus1*xPower] * xFactor;
      }
      phi[k] *= phiTemp*legendre[*index];
   }
   free(jacobi);

   return;
}// BasisSimplexPolyhedralTwo


// Generates derivatives of orthogonal polynomial basis for the first [0, dim1-1] coordinates
// over simplex in dim-dimensional space.
static void BasisPrimeSimplexPolyhedralOne(int dim, int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime)
{
   assert(dims[0] >= 2);

   int dim1 = dims[0];
   int degPlus1 = deg+1;
   int num_funcs = BasisSize(deg, dim);

   if(dim1 == 2)
   {
      int i, j, k, d;
      double factor1, factor2;
      double legendre[degPlus1];
      double dxlegendre[degPlus1];
      double jacobi[SQUARE(degPlus1)];
      double dxjacobi[SQUARE(degPlus1)];

      LegendrePoly(degPlus1, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
      for(i = 0; i < degPlus1; ++i) JacobiPoly(degPlus1, 1.0-2*x[0], 2.0*i+1, 0.0, &jacobi[i*degPlus1]);
      for(i = 0; i < degPlus1; ++i) JacobiPolyPrime(degPlus1, 1.0-2*x[0], 2.0*i+1, 0.0, &dxjacobi[i*degPlus1]);
      for(k = 0; k < num_funcs*dim; ++k) phiPrime[k] = 1.0;

      int power1, power2, index;
      for(int k = 0; k < num_funcs; ++k)
      {
         power1 = basis_id[k*dim];
         power2 = basis_id[k*dim+1];
         index = power2+power1*degPlus1;
         factor1 = 1.0; for(j = 0; j < power1; ++j) factor1 *= x[0];
         factor2 = 1.0; for(j = 0; j < power1-1; ++j) factor2 *= x[0];

         phiPrime[k] = -2.0*x[1]/SQUARE(x[0]) * dxlegendre[power1] * factor1 * jacobi[index] +
                        legendre[power1] * ( power1 * factor2 * jacobi[index] + factor1 * dxjacobi[index] * (-2) );
         phiPrime[k+num_funcs] = factor1 * jacobi[index] * dxlegendre[power1] * 2.0 / x[0];

         for(d = 2; d < dim; ++d)
            phiPrime[k+d*num_funcs] *= factor1 * jacobi[power2+power1*degPlus1] * legendre[power1];
      }
   }


   if(dim1 >= 3)
   {
      int i, j, d;
      double jacobiX[dim1];
      double polyX[dim1];
      double primeFactor[dim1];

      double *jacobi = (double *)malloc( SQUARE(degPlus1)*(dim1)*size_double );
      double *dxjacobi = (double *)malloc( SQUARE(degPlus1)*(dim1)*size_double );
      LegendrePoly(degPlus1, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], &jacobi[0], &dxjacobi[0]);

      double *phi = (double *)malloc(num_funcs*sizeof(double));
      BasisSimplexPolyhedralOne(dim, dims, deg, basis_id, x, phi);
      for(int k = 0; k < num_funcs*dim; ++k) phiPrime[k] = 1.0;

      for(d = 0; d < dim1; ++d)   polyX[d] = 1.0;
      for(d = 1; d < dim1-1; ++d) polyX[d] = x[dim1-d-1]/x[dim1-d-2];
      polyX[0] = 0; polyX[dim1-1] = x[0];

      for(d = 0; d < dim1; ++d) jacobiX[d] = 1.0;
      for(d = 1; d < dim1; ++d)
      {
         int dimCur = d;
         jacobiX[0] = (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2];
         if(d < dim1-1)
         {
            jacobiX[dimCur] /= x[dim1-d-2];
            jacobiX[dimCur] = 1.0-2.0*x[dim1-d-1] * jacobiX[dimCur];
         }
         if(d == dim1-1)
         {
            jacobiX[dimCur] *= x[0];
            jacobiX[dimCur] = 1.0-2.0*jacobiX[dimCur];
         }

         int alpha = 0;
         for(j = 0; j < deg+1; ++j)
         {
            alpha = 2*j+d;
            JacobiPoly(degPlus1, jacobiX[dimCur], alpha, 0.0, &jacobi[SQUARE(degPlus1)*d+degPlus1*j]);
            JacobiPolyPrime(degPlus1, jacobiX[dimCur], alpha, 0.0, &dxjacobi[SQUARE(degPlus1)*d+degPlus1*j]);
         }
      }

      const INT_8 *index;
      double polyTemp[4];
      double polyFactor[5];
      for(int k = 0; k < num_funcs; ++k)
      {
         int dimPrev = 0;
         int xPowerPrev = 0;
         int jacobiIndexPrev = 0;
         index = &basis_id[k*dim];

         for(j = 0; j < dim1; ++j)
         {
            int xPower = 0; for(i = 0; i < j; ++i) xPower += *(index+i);

            polyTemp[0] = polyX[j];
            double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= polyTemp[0];
            int jacobiIndex = SQUARE(degPlus1)*j+degPlus1*xPower;
            for(d = 0; d < dim1; ++d) primeFactor[d] = jacobi[*(index+j)+jacobiIndex] * xFactor;

            if(j == 0)
            {
               primeFactor[dim1-1] = dxjacobi[*(index)+jacobiIndex] * 2.0/x[dim1-2];
               primeFactor[dim1-2] = 1.0;
            }
            else if(j == 1)
            {
               polyTemp[0] = x[dim1-3];
               polyTemp[1] = x[dim1-2];
               polyFactor[0] = x[dim1-3];
               polyFactor[1] = 1.0;
               for(i = 0; i < xPower-1; ++i)
               {
                  polyFactor[0] *= polyTemp[0];
                  polyFactor[1] *= polyTemp[1];
               }
               primeFactor[dim1-2] = jacobi[*(index+1) + jacobiIndex] * xFactor * dxjacobi[*(index) + jacobiIndexPrev] * (-2.0) * x[dim1-1] / SQUARE(x[dim1-2]) +
                                     jacobi[*(index) + jacobiIndexPrev] * ( dxjacobi[*(index+1) + jacobiIndex] * (-2.0)/x[dim1-3] * xFactor +
                                     jacobi[*(index+1) + jacobiIndex] * xPower * polyFactor[1] / polyFactor[0] );
            }
            else if(j == dim1-1)
            {
               int d1 = basis_id[k*dim + dim1-1];
               int d2 = basis_id[k*dim + dim1-2];
               polyTemp[0] = x[1]/x[0];
               polyTemp[1] = x[0];
               polyTemp[2] = x[1];

               for(int s = 0; s < 5; ++s) polyFactor[s] = 1.0;
               for(i = 0; i < xPowerPrev; ++i)   polyFactor[0] *= polyTemp[0];
               for(i = 0; i < xPower; ++i)       polyFactor[1] *= polyTemp[1];
               for(i = 0; i < xPowerPrev; ++i)   polyFactor[2] *= polyTemp[2];
               for(i = 1; i < xPower; ++i)       polyFactor[3] *= polyTemp[1];
               for(i = 0; i < xPowerPrev+1; ++i) polyFactor[4] *= polyTemp[1];

               double v1 = polyFactor[0] * jacobi[d2+jacobiIndexPrev] *
                           ( dxjacobi[d1+jacobiIndex] * (-2) * polyFactor[1] + jacobi[d1+jacobiIndex] * xPower * polyFactor[3] );

               double v2 = jacobi[d1+jacobiIndex] * polyFactor[1] *
                           ( polyFactor[0] * dxjacobi[d2+jacobiIndexPrev] * 2.0 * x[1]/SQUARE(x[0]) +
                            jacobi[d2+jacobiIndexPrev] * polyFactor[2] * (-xPowerPrev) / polyFactor[4] );

               primeFactor[0] = v1 + v2;
            }

            if(j <= dim1-2) primeFactor[dim1-j-2] = 1.0;
            if(j == dim1-2) primeFactor[0] = 1.0;

            if( (dim1 > 3) && (j >= 2) && (j != dim1-1) )
            {
               int d1 = basis_id[k*dim+j];
               int d2 = basis_id[k*dim+j-1];
               int prev = dim1-j-2;
               int cur = dim1-j-1;
               int next = dim1-j;

               polyTemp[0] = x[cur]/x[prev];
               polyTemp[1] = x[next]/x[cur];
               polyTemp[2] = x[next];
               polyTemp[3] = x[cur];

               for(int s = 0; s < 5; ++s) polyFactor[s] = 1.0;
               for(i = 0; i < xPower; ++i)       polyFactor[0] *= polyTemp[0];
               for(i = 0; i < xPowerPrev; ++i)   polyFactor[1] *= polyTemp[1];
               for(i = 0; i < xPowerPrev; ++i)   polyFactor[2] *= polyTemp[2];
               for(i = 0; i < xPowerPrev+1; ++i) polyFactor[3] *= polyTemp[3];
               for(i = 0; i < xPower-1; ++i)     polyFactor[4] *= polyTemp[0];

               double v1 = polyFactor[1] * jacobi[d2+jacobiIndexPrev] *
                           ( dxjacobi[d1+jacobiIndex] * (-2.0)/x[prev] * polyFactor[0] + jacobi[d1+jacobiIndex] * (xPower)*polyFactor[4]*1.0/x[prev] );

               double v2 = jacobi[d1+jacobiIndex] * polyFactor[0] *
                           ( polyFactor[1] * dxjacobi[d2+jacobiIndexPrev] * 2.0 * x[next] / SQUARE(x[cur]) +
                             jacobi[d2+jacobiIndexPrev] * polyFactor[2] * (-xPowerPrev) / polyFactor[3] );

               primeFactor[dim1-1-j] = v1 + v2;
            }

            for(d = 0; d < dim1; ++d)
               phiPrime[d*num_funcs+k] *= primeFactor[d];

            dimPrev = j;
            xPowerPrev = xPower;
            jacobiIndexPrev = SQUARE(degPlus1)*dimPrev + degPlus1*xPowerPrev;
         }
      }

      for(d = dim1; d < dim; ++d)
         for(int k = 0; k < num_funcs; ++k)
            phiPrime[k+d*num_funcs] *= phi[k];

      free(jacobi);
      free(dxjacobi);
      free(phi);

   }// end if dim >= 3


}// BasisPrimeSimplexPolyhedralOne


// Generates derivatives of orthogonal polynomial basis for coordinates
// between [dim1, dim2-1] over simplex in dimal-dimensional space.
static void BasisPrimeSimplexPolyhedralTwo(int dim, int *dims, int deg, const INT_8 *basis_id, const double *x, double *phiPrime)
{
   assert(dims[1] >= 2);

   int dim1 = dims[0];
   int dim2 = dims[1];
   int dimTwo = dim1+dim2;
   int degPlus1 = deg+1;
   int num_funcs = BasisSize(deg, dim);


   if(dim2 == 2)
   {
      int i, j, k, d;
      double factor1, factor2;
      double legendre[degPlus1];
      double dxlegendre[degPlus1];
      double jacobi[SQUARE(degPlus1)];
      double dxjacobi[SQUARE(degPlus1)];

      LegendrePoly(degPlus1, (2.0*x[dim1+1]-x[dim1])/x[dim1], legendre, dxlegendre);
      for(i = 0; i < degPlus1; ++i) JacobiPoly(degPlus1, 1.0-2*x[dim1], 2.0*i+1, 0.0, &jacobi[i*degPlus1]);
      for(i = 0; i < degPlus1; ++i) JacobiPolyPrime(degPlus1, 1.0-2*x[dim1], 2.0*i+1, 0.0, &dxjacobi[i*degPlus1]);
      for(k = 0; k < num_funcs*dim; ++k) phiPrime[k] = 1.0;

      int power1, power2, index;
      for(int k = 0; k < num_funcs; ++k)
      {
         power1 = basis_id[k*dim+dim1];
         power2 = basis_id[k*dim+dim1+1];
         index = power2+power1*degPlus1;
         factor1 = 1.0; for(j = 0; j < power1; ++j) factor1 *= x[dim1];
         factor2 = 1.0; for(j = 0; j < power1-1; ++j) factor2 *= x[dim1];

         phiPrime[k+dim1*num_funcs] = -2.0*x[dim1+1]/SQUARE(x[dim1]) * dxlegendre[power1] * factor1 * jacobi[index] +
                                      legendre[power1] * ( power1 * factor2 * jacobi[index] + factor1 * dxjacobi[index] * (-2) );
         phiPrime[k+(dim1+1)*num_funcs] = factor1 * jacobi[index] * dxlegendre[power1] * 2.0/x[dim1];

         for(d = 0; d < dim1; ++d)
            phiPrime[k+d*num_funcs] *= factor1 * jacobi[power2+power1*degPlus1] * legendre[power1];
         for(d = dimTwo; d < dim; ++d)
            phiPrime[k+d*num_funcs] *= factor1 * jacobi[power2+power1*degPlus1] * legendre[power1];
      }

      return;
   }


   if(dim2 >= 3)
   {
      int i, j, k, d;
      double jacobiX[dim2];
      double polyX[dim2];
      double primeFactor[dim2];

      double *jacobi = (double *)malloc( SQUARE(degPlus1)*(dim2)*size_double );
      double *dxjacobi = (double *)malloc( SQUARE(degPlus1)*(dim2)*size_double );
      LegendrePoly(degPlus1, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], &jacobi[0], &dxjacobi[0]);

      for(k = 0; k < num_funcs*dim; ++k) phiPrime[k] = 1.0;
      double *phi = (double *)malloc( SIZE_DOUBLE(num_funcs) );
      BasisSimplexPolyhedralTwo(dim, dims, deg, basis_id, x, phi);

      for(d = 0; d < dim2; ++d)   polyX[d] = 1.0;
      for(d = 1; d < dim2-1; ++d) polyX[d] = x[dimTwo-d-1]/x[dimTwo-d-2];
      polyX[0] = 0; polyX[dim2-1] = x[dim1];

      for(d = 0; d < dim2; ++d) jacobiX[d] = 1.0;
      for(d = 1; d < dim2; ++d)
      {
         int dimCur = d;
         jacobiX[0] = (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2];
         if(d < dim2-1)
         {
            jacobiX[dimCur] /= x[dimTwo-d-2];
            jacobiX[dimCur] = 1.0-2.0*x[dimTwo-d-1] * jacobiX[dimCur];
         }
         if(d == dim2-1)
         {
            jacobiX[dimCur] *= x[dim1];
            jacobiX[dimCur] = 1.0-2.0*jacobiX[dimCur];
         }

         int alpha = 0;
         for(j = 0; j < deg+1; ++j)
         {
            alpha = 2*j+d;
            JacobiPoly(degPlus1, jacobiX[dimCur], alpha, 0.0, &jacobi[SQUARE(degPlus1)*d+degPlus1*j]);
            JacobiPolyPrime(degPlus1, jacobiX[dimCur], alpha, 0.0, &dxjacobi[SQUARE(degPlus1)*d+degPlus1*j]);
         }
      }

      const INT_8 *index;
      double polyTemp[4];
      double polyFactor[5];
      for(int k = 0; k < num_funcs; ++k)
      {
         int dimPrev = 0;
         int xPowerPrev = 0;
         int jacobiIndexPrev = 0;
         index = &basis_id[k*dim+dim1];

         for(j = 0; j < dim2; ++j)
         {
            int xPower = 0; for(i = 0; i < j; ++i) xPower += *(index+i);

            polyTemp[0] = polyX[j];
            double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= polyTemp[0];
            int jacobiIndex = SQUARE(degPlus1)*j+degPlus1*xPower;
            for(d = 0; d < dim2; ++d) primeFactor[d] = jacobi[*(index+j)+jacobiIndex] * xFactor;

            if(j == 0)
            {
               primeFactor[dim2-1] = dxjacobi[*(index)+jacobiIndex] * 2.0/x[dimTwo-2];
               primeFactor[dim2-2] = 1.0;
            }
            else if(j == 1)
            {
               polyTemp[0] = x[dimTwo-3];
               polyTemp[1] = x[dimTwo-2];
               polyFactor[0] = x[dimTwo-3];
               polyFactor[1] = 1.0;
               for(i = 0; i < xPower-1; ++i)
               {
                  polyFactor[0] *= polyTemp[0];
                  polyFactor[1] *= polyTemp[1];
               }
               primeFactor[dim2-2] = jacobi[*(index+1) + jacobiIndex] * xFactor * dxjacobi[*(index) + jacobiIndexPrev] *
                                     (-2.0) * x[dimTwo-1] / SQUARE(x[dimTwo-2]) + jacobi[*(index) + jacobiIndexPrev] *
                                     ( dxjacobi[*(index+1) + jacobiIndex] * (-2.0)/x[dimTwo-3] * xFactor +
                                      jacobi[*(index+1) + jacobiIndex] * xPower * polyFactor[1] / polyFactor[0] );
            }
            else if(j == dim2-1)
            {
               int d1 = basis_id[k*dim + dimTwo-1];
               int d2 = basis_id[k*dim + dimTwo-2];
               polyTemp[0] = x[dim1+1]/x[dim1];
               polyTemp[1] = x[dim1];
               polyTemp[2] = x[dim1+1];

               for(int s = 0; s < 5; ++s) polyFactor[s] = 1.0;
               for(i = 0; i < xPowerPrev; ++i)   polyFactor[0] *= polyTemp[0];
               for(i = 0; i < xPower; ++i)       polyFactor[1] *= polyTemp[1];
               for(i = 0; i < xPowerPrev; ++i)   polyFactor[2] *= polyTemp[2];
               for(i = 1; i < xPower; ++i)       polyFactor[3] *= polyTemp[1];
               for(i = 0; i < xPowerPrev+1; ++i) polyFactor[4] *= polyTemp[1];

               double v1 = polyFactor[0] * jacobi[d2+jacobiIndexPrev] *
                           ( dxjacobi[d1+jacobiIndex] * (-2.0) * polyFactor[1] + jacobi[d1+jacobiIndex] * xPower * polyFactor[3] );

               double v2 = jacobi[d1+jacobiIndex] * polyFactor[1] *
                           ( polyFactor[0] * dxjacobi[d2+jacobiIndexPrev] * 2.0 * x[dim1+1] / SQUARE(x[dim1]) +
                             jacobi[d2+jacobiIndexPrev] * polyFactor[2] * (-xPowerPrev) / polyFactor[4] );

               primeFactor[0] = v1 + v2;
            }

            if(j <= dim2-2) primeFactor[dim2-j-2] = 1.0;
            if(j == dim2-2) primeFactor[0] = 1.0;

            if( (dim2 > 3) && (j >= 2) && (j != dim2-1) )
            {
               int d1 = basis_id[k*dim+dim1+j];
               int d2 = basis_id[k*dim+dim1+j-1];
               int prev = dimTwo-j-2;
               int cur = dimTwo-j-1;
               int next = dimTwo-j;

               polyTemp[0] = x[cur]/x[prev];
               polyTemp[1] = x[next]/x[cur];
               polyTemp[2] = x[next];
               polyTemp[3] = x[cur];

               for(int s = 0; s < 5; ++s) polyFactor[s] = 1.0;
               for(i = 0; i < xPower; ++i)       polyFactor[0] *= polyTemp[0];
               for(i = 0; i < xPowerPrev; ++i)   polyFactor[1] *= polyTemp[1];
               for(i = 0; i < xPowerPrev; ++i)   polyFactor[2] *= polyTemp[2];
               for(i = 0; i < xPowerPrev+1; ++i) polyFactor[3] *= polyTemp[3];
               for(i = 0; i < xPower-1; ++i)     polyFactor[4] *= polyTemp[0];

               double v1 = polyFactor[1] * jacobi[d2+jacobiIndexPrev] *
                           ( dxjacobi[d1+jacobiIndex] * (-2)/x[prev] * polyFactor[0] +
                             jacobi[d1+jacobiIndex] * xPower * polyFactor[4] / x[prev] );

               double v2 = jacobi[d1+jacobiIndex] * polyFactor[0] *
                           ( polyFactor[1] * dxjacobi[d2+jacobiIndexPrev] * 2.0 * x[next] / SQUARE(x[cur]) +
                            jacobi[d2+jacobiIndexPrev] * polyFactor[2] * (-xPowerPrev) / polyFactor[3] );

               primeFactor[dim2-1-j] = v1 + v2;
            }

            for(d = 0; d < dim2; ++d)
               phiPrime[(d+dim1)*num_funcs+k] *= primeFactor[d];

            dimPrev = j;
            xPowerPrev = xPower;
            jacobiIndexPrev = SQUARE(degPlus1)*dimPrev + degPlus1*xPowerPrev;
         }
      }

      for(d = 0; d < dim1; ++d)
         for(int k = 0; k < num_funcs; ++k)
            phiPrime[k+d*num_funcs] *= phi[k];
      for(d = dimTwo; d < dim; ++d)
         for(int k = 0; k < num_funcs; ++k)
            phiPrime[k+d*num_funcs] *= phi[k];

      free(jacobi);
      free(dxjacobi);
      free(phi);

      return;
   }// end if dim >= 3

}// BasisPrimeSimplexPolyhedralTwo


void BasisMonomial(int dim, int deg, const INT_8 *basis_id, const double *x, double *phi)
{
   int num_funs = BasisSize(dim, deg);
   for(int k = 0; k < num_funs; ++k)
      phi[k] = 1.0;

   for(int k = 0; k < num_funs; ++k)
   {
      for(int d = 0; d < dim; ++d)
      {
         int basis_power = basis_id[k*dim+d];
         double factor = x[d];
         double product = 1.0;

         for(int r = 0; r < basis_power; ++r)
            product *= factor;

         phi[k] *= product;
      }
   }

}


static int finite_difference_test(int dim, int *dims, int deg, const INT_8 *basis_id,
                                  const double *x, const double *phi_prime, BasisFunc phi_func)
{
   int num_funcs = BasisSize(deg, dim);

   double h = POW_DOUBLE(10, -5)*5.0;
   double tol = POW_DOUBLE(10, -8);

   double *phi_approx = (double *)malloc( SIZE_DOUBLE(num_funcs*dim) );
   double *phi_backw1 = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_backw2 = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_forw1  = (double *)malloc( SIZE_DOUBLE(num_funcs) );
   double *phi_forw2  = (double *)malloc( SIZE_DOUBLE(num_funcs) );

   double x_backw1[dim];
   double x_backw2[dim];
   double x_forw1[dim];
   double x_forw2[dim];

   int flag = 0;

   for(int d = 0; d < dim; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t)
      {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t] = x[d_t];
         x_forw2[d_t] = x[d_t];
      }
      x_backw1[d] = x_backw1[d] - h;
      x_backw2[d] = x_backw2[d] - 2.0*h;
      x_forw1[d] = x_forw1[d] + h;
      x_forw2[d] = x_forw2[d] + 2.0*h;

      phi_func(dims, deg, basis_id, x_backw1, phi_backw1);
      phi_func(dims, deg, basis_id, x_backw2, phi_backw2);
      phi_func(dims, deg, basis_id, x_forw1, phi_forw1);
      phi_func(dims, deg, basis_id, x_forw2, phi_forw2);

      for(int k = 0; k < num_funcs; ++k)
         phi_approx[d*num_funcs+k] = ( 1.0/12.0*phi_backw2[k] - 2.0/3.0*phi_backw1[k] + 2.0/3.0*phi_forw1[k]-1.0/12.0*phi_forw2[k] ) / h;

      for(int k = 0; k < num_funcs; ++k)
         if( fabs( phi_prime[d*num_funcs+k] - phi_approx[d*num_funcs+k] ) > tol )
         {
            flag = 1;
//            goto FREERETURN;
         }
   }

//   if(flag == 1 || dim == 3)
//   {
//
//   double maxError = fabs(phi_approx[0] - phi_prime[0]);
//   for(int i = 1; i < dim*num_funcs; ++i)
//      maxError = MAX(maxError, fabs(phi_approx[i] - phi_prime[i]));
//   printf("maximum error = %.16e\n", maxError);
//   }

FREERETURN:
   free(phi_backw1);
   free(phi_backw2);
   free(phi_forw1);
   free(phi_forw2);
   free(phi_approx);

   return flag;
}


void orthogonal_simplex_basis_test(int deg, int dim)
{
   int num_funcs = BasisSize(deg, dim);

   int n = floor(deg/2)+1 + floor(dim/2) + floor(deg/6) + 3; // make n large enough to get exact quadrature
   int dims_1D[1] = {1};
   quadrature *quad_1D = quadrature_init_basic(n, 1, dims_1D, deg, INTERVAL);
   Jacobi(quad_1D->num_nodes, 0.0, 0.0, quad_1D->x, quad_1D->w);

   int N = POW_INT(n, dim);
   quadrature *quad_S = quadrature_init_basic(N, dim, &dim, deg, SIMPLEX);
   GeneralizedNodesTensor(quad_1D, quad_S);
   GeneralizedWeightsTensor(quad_1D, quad_S);
   GeneralDuffy(quad_S);

   double *basis = (double *)calloc(num_funcs, size_double);
   double *quad_integrals = (double *)calloc(num_funcs, size_double);

   double *exact_integrals = (double *)calloc(num_funcs, size_double);
   BasisIntegralsSimplex(&dim, deg, exact_integrals);

   INT_8 *basis_indices = (int8_t *)malloc((num_funcs*dim)*sizeof(INT_8));
   BasisIndices(deg, dim, basis_indices);

   for(int i = 0; i < N; ++i)
   {
      BasisSimplex(&dim, deg, basis_indices, &quad_S->x[dim*i], basis);
      for(int j = 0; j < num_funcs; ++j)
         quad_integrals[j] += basis[j]*quad_S->w[i];
   }

   double max_res = fabs(quad_integrals[0] - exact_integrals[0]);
   for(int j = 1; j < num_funcs; ++j)
   {
      double res = fabs(quad_integrals[j]-exact_integrals[j]);
      max_res = MAX(max_res, res);
   }
   printf("\nTesting orthogonality of basis functions. Maximum error of basis functions = %.16e\n\n", max_res);

   free(basis_indices);
   free(basis);
   free(quad_integrals);
   free(exact_integrals);
   quadrature_free(quad_1D);
   quadrature_free(quad_S);

}
