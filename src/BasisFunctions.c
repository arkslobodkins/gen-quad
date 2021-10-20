/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "BasisFunctions.h"

#include "BasisIndices.h"
#include "GeneralGaussTensor.h"
#include "Gauss_Lib/Jacobi.h"
#include "Quadrature.h"
#include "AddDimension.h"
#include "GENERAL_QUADRATURE.h"
#include "Print.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>


static void PhiSimplexPolyhedralOne(const int_fast8_t *basis_id, const double *x,
                                    const quadParams *params, double *phi);

static void PhiSimplexPolyhedralTwo(const int_fast8_t *basis_id, const double *x,
                                    const quadParams *params, double *phi);

static void PhiPrimeSimplexPolyhedralOne(const int_fast8_t *basis_id, const double *x,
                                         const quadParams *params, double *phiPrime);

static void PhiPrimeSimplexPolyhedralTwo(const int_fast8_t *basis_id, const double *x,
                                         const quadParams *params, double *phiPrime);

typedef void(*PhiFunc)(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi);

static int finite_difference_test(const int_fast8_t *basis_id, const double *x,
                                  const quadParams *params, const double *phi_prime, PhiFunc func);

static void LegendrePoly(int order, double x, double *p, double *dp);
static void LegendrePoly(int order, double x, double *p, double *dp)
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
      fac2 =k/(k+1.0);
      p[k+1] = fac1*x*p[k] - fac2*p[k-1];
      dp[k+1] = fac1*(p[k] + x*dp[k]) - fac2*dp[k-1];
   }

}


#define SQUARE(x) ((x)*(x))
#define ADD(x, y) (x)+(y)
#define SUB(x, y) (x)-(y)

static void JacobiPoly(int order, double x, double alpha, double beta, double *p);
static void JacobiPolyPrime(int order, double x, double alpha, double beta, double *dp);
static void JacobiPoly(int order, double x, double alpha, double beta, double *p)
{
   // compute constants to be used inside a for loop for improved performance
   double fac1 = 0.0, fac2 = 0.0, fac3 = 0.0;
   double alpha_sq = SQUARE(alpha);
   double beta_sq = SQUARE(beta);
   double a_plus_b = ADD(alpha, beta);
   double a_plus_b_minus_1 = SUB(a_plus_b, 1.0);
   double a_plus_b_minus_2 = SUB(a_plus_b, 2.0);
   double alpha_sq_minus_beta_sq = SUB(alpha_sq, beta_sq);
   double a_minus_1 = SUB(alpha, 1.0);
   double b_minus_1 = SUB(beta, 1.0);

   p[0] = 1.0;
   if(order > 1)
      p[1] = 0.5*(x-1.0) * (alpha+beta+2.0) + alpha+1.0;

   for (int k = 1; k < order-1; ++k)
   {
      int k_plus_1 = k+1;
      int two_n = 2*k_plus_1;
      int two_n_plus_a_plus_b = ADD(two_n, a_plus_b);
      fac1 = (2.0*k_plus_1+a_plus_b_minus_1) * (two_n_plus_a_plus_b * (two_n_plus_a_plus_b-2.0) * x + alpha_sq_minus_beta_sq);
      fac2 = -2.0 * (k_plus_1+a_minus_1) * (k_plus_1+b_minus_1) * two_n_plus_a_plus_b;
      fac3 = 1.0 / (two_n*(k_plus_1+a_plus_b) * (two_n+a_plus_b_minus_2));
      p[k_plus_1] = fac3 * (fac1*p[k] + fac2*p[k-1]);
   }

}


static void JacobiPolyPrime(int order, double x, double alpha, double beta, double *dp)
{
   double p[order];
   double a_plus_b_plus1 = 1.0+alpha+beta;

   JacobiPoly(order, x, alpha+1.0, beta+1.0, p);
   dp[0] = 0.0;
   dp[1] = 0.5 * (alpha+beta+2.0);

   for(int k = 2; k < order; ++k)
      dp[k] = 0.5*(a_plus_b_plus1+k) * p[k-1];

}


/* PhiCube
 * Generates orthogonal polynomial basis for the unit cube
 */
void PhiCube(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi)
{
   assert(params->dim >= 1);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim = params->dim;
   int deg = params->deg;
   int order = params->deg+1;
   int num_funs = BasisSize(dim, deg);

   double legendre[order*dim];
   double dxlegendre[order*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePoly(order, 2*x[d]-1, &legendre[d*order], &dxlegendre[d*order]);

   for(int k = 0; k < num_funs; ++k) phi[k] = 1.0;

   for(int k = 0; k < num_funs; ++k)
      for(int d = 0; d < dim; ++d)
         phi[k] *= legendre[basis_id[k*dim+d] + order*d];

}// end PhiCube


/* PhiPrimeCube
 * Generates derivatives of orthogonal polynomial basis for the unit cube
 */
void PhiPrimeCube(const int_fast8_t *basis_id, const double* x, const quadParams *params, double *phiPrime)
{
   assert(params->dim >= 1);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim = params->dim;
   int order = params->deg+1;
   int num_funs = BasisSize(dim, params->deg);

   double legendre[order*dim];
   double dxlegendre[order*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePoly(order, 2*x[d]-1, &legendre[d*order], &dxlegendre[d*order]);

   for(int k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

   for(int d = 0; d < dim; ++d)
   {

      // dimension < d
      for(int j = 0; j < d; ++j)
      {
         double *legendre_ptr = &legendre[order*j];
         for(int k = 0; k < num_funs; ++k)
            phiPrime[d*num_funs+k] *= legendre_ptr[basis_id[dim*k+j]];
      }
      // dimension = d
      double *dxlegendre_ptr = &dxlegendre[order*d];
      for(int k = 0; k < num_funs; ++k)
         phiPrime[d*num_funs+k] *= 2.0 * dxlegendre_ptr[basis_id[dim*k+d]];

      // dimension > d
      for(int j = d+1; j < dim; ++j)
      {
         double *legendre_ptr = &legendre[order*j];
         for(int k = 0; k < num_funs; ++k)
            phiPrime[d*num_funs+k] *= legendre_ptr[basis_id[dim*k+j]];
      }

   }

#ifdef QUAD_DEBUG_ON
   PhiFunc phi_func = &PhiCube;
   int flag = finite_difference_test(basis_id, x, params, phiPrime, phi_func);
//   if(flag == 1)
//      PRINT_ERR("failed finite_difference_test", __LINE__, __FILE__);
#endif

}// end PhiPrimeCube


/* PhiSimplex
 * Generates orthogonal polynomial basis for the unit simplex
 */
void PhiSimplex(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi)
{
   assert(params->dim >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim = params->dim;
   int deg = params->deg;
   int order = deg+1;
   int num_funs = BasisSize(dim, deg);

   if(dim == 2)
   {
      double legendre[order];
      double dxlegendre[order];
      double jacobi[order*order];

      LegendrePoly(order, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
      for(int j = 0; j < order; ++j) JacobiPoly(order, 1.0-2*x[0], 2.0*j+1, 0.0, &jacobi[j*order]);
      for(int k = 0; k < num_funs; ++k) phi[k] = 1.0;

      int r, n;
      for(int k = 0; k < num_funs; ++k)
      {
         r = basis_id[dim*k];
         n = basis_id[dim*k+1];
         double xFactor = 1.0; for(int i = 0; i < r; ++i) xFactor *= x[0];
         phi[k] *= xFactor * jacobi[r*order+n] * legendre[r];
      }

      return;
   }

   // dim > 3
   int i, j, k, d;
   double legendre[order];
   double dxlegendre[order];

   double *jacobi = (double *) malloc( SQUARE(order)*(dim-1)*size_double );
   LegendrePoly(order, (2.0*x[dim-1]-x[dim-2])/x[dim-2], legendre, dxlegendre);
   for(k = 0; k < num_funs; ++k) phi[k] = 1.0;

   double polyX[dim-1];
   for(d = 0; d < dim-1; ++d) polyX[d] = 1.0;
   for(d = 0; d < dim-2; ++d) polyX[d] = x[dim-d-2]/x[dim-d-3];
   polyX[dim-2] = x[0];

   double jacobiX[dim-1];
   for(d = 1; d < dim; ++d)
   {
      int dimCur = d-1;
      double jacobiX_dimCur = 1.0;
      if(d < dim-1)
      {
         jacobiX_dimCur /= x[dim-d-2];
         jacobiX_dimCur = 1.0-2.0*x[dim-d-1] * jacobiX_dimCur;
      }
      else if(d == dim-1)
      {
         jacobiX_dimCur *= x[0];
         jacobiX_dimCur = 1.0-2.0*jacobiX_dimCur;
      }

      double alpha = 0.0;
      for(j = 0; j < order; ++j)
      {
         alpha = 2*j+d;
         JacobiPoly(order, jacobiX_dimCur, alpha, 0.0, &jacobi[order*(d-1)*deg+order*j]);
      }
      jacobiX[dimCur] = jacobiX_dimCur;
   }

   for(k = 0; k < num_funs; ++k)
   {
      const int_fast8_t *index = &basis_id[k*dim];
      double phiTemp = 1.0;
      for(d = 1; d < dim; ++d)
      {
         double xTemp = polyX[d-1];
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= xTemp;
         phiTemp *= jacobi[*(index+d) + order*(d-1)*deg + order*xPower] * xFactor;
      }

      phi[k] *= phiTemp * legendre[*index];
   }
   free(jacobi);

   return;
}// end PhiSimplex


/* PhiPrimeSimplex
 * Generates derivatives of orthogonal polynomial basis for the unit simplex
 */
void PhiPrimeSimplex(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phiPrime)
{
   assert(params->dim >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim = params->dim;
   int deg = params->deg;
   int order = deg+1;
   int num_funs = BasisSize(dim, deg);


   if(dim == 2)
   {
      int i, j, k;
      double factor1, factor2;
      double legendre[order];
      double dxlegendre[order];
      double jacobi[SQUARE(order)];
      double dxjacobi[SQUARE(order)];
      for(int k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

      LegendrePoly(order, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
      for(i = 0; i < order; ++i) JacobiPoly(order, 1.0-2*x[0], 2.0*i+1, 0.0, &jacobi[i*order]);
      for(i = 0; i < order; ++i) JacobiPolyPrime(order, 1.0-2*x[0], 2.0*i+1, 0.0, &dxjacobi[i*order]);

      int power1, power2, index;
      for(k = 0; k < num_funs; ++k)
      {
         power1 = basis_id[k*dim];
         power2 = basis_id[k*dim+1];
         index = power2+power1*order;
         factor1 = 1.0; for(j = 0; j < power1; ++j) factor1 *= x[0];
         factor2 = 1.0; for(j = 0; j < power1-1; ++j) factor2 *= x[0];

         phiPrime[k] = -2.0 * x[1] / SQUARE(x[0]) * dxlegendre[power1] * factor1 * jacobi[index] +
                        legendre[power1] * ( power1 * factor2 * jacobi[index] + factor1 * dxjacobi[index] * (-2.0) );
         phiPrime[k+num_funs] = factor1 * jacobi[index] * dxlegendre[power1] * 2.0/x[0];
      }

   }


   if(dim >= 3)
   {
      int i, j, k, d;
      double *phi = (double *)malloc(num_funs*sizeof(double));
      double jacobiX[dim];
      double polyX[dim];
      double primeFactor[dim];

      double *jacobi = (double *)malloc( SQUARE(order)*dim*size_double );
      double *dxjacobi = (double *)malloc( SQUARE(order)*dim*size_double );
      LegendrePoly(order, (2.0*x[dim-1]-x[dim-2])/x[dim-2], &jacobi[0], &dxjacobi[0]);
      for(int k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

      for(d = 0; d < dim; ++d)   polyX[d] = 1.0;
      for(d = 1; d < dim-1; ++d) polyX[d] = x[dim-d-1]/x[dim-d-2];
      polyX[0] = 0; polyX[dim-1] = x[0];

      for(d = 0; d < dim; ++d) jacobiX[d] = 1.0;
      for(d = 1; d < dim; ++d)
      {
         int dimCur = d;
         jacobiX[0] = (2.0*x[dim-1]-x[dim-2])/x[dim-2];
         if(d < dim-1)
         {
            jacobiX[dimCur] /= x[dim-d-2];
            jacobiX[dimCur] = 1.0-2.0*x[dim-d-1] * jacobiX[dimCur];
         }
         if(d == dim-1)
         {
            jacobiX[dimCur] *= x[0];
            jacobiX[dimCur] = 1.0-2.0*jacobiX[dimCur];
         }

         int alpha = 0;
         for(j = 0; j < deg+1; ++j)
         {
            alpha = 2*j+d;
            JacobiPoly(order, jacobiX[dimCur], alpha, 0.0, &jacobi[SQUARE(order)*d+order*j]);
            JacobiPolyPrime(order, jacobiX[dimCur], alpha, 0.0, &dxjacobi[SQUARE(order)*d+order*j]);
         }

      }

      const int_fast8_t *index;
      double polyTemp[4];
      double polyFactor[5];
      for(k = 0; k < num_funs; ++k)
      {
         int dimPrev = 0;
         int xPowerPrev = 0;
         int jacobiIndexPrev = 0;
         double xFactorPrev = 0.0;
         index = &basis_id[k*dim];

         for(j = 0; j < dim; ++j)
         {
            int xPower = 0; for(i = 0; i < j; ++i) xPower += *(index+i);

            polyTemp[0] = polyX[j];
            double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= polyTemp[0];
            int jacobiIndex = SQUARE(order)*j+order*xPower;
            for(d = 0; d < dim; ++d) primeFactor[d] = jacobi[*(index+j)+jacobiIndex] * xFactor;

            if(j == 0)
            {
               primeFactor[dim-1] = dxjacobi[*(index)+jacobiIndex] * 2.0/x[dim-2];
               primeFactor[dim-2] = 1.0;
            }
            else if(j == 1)
            {
               polyTemp[0] = x[dim-3];
               polyTemp[1] = x[dim-2];
               polyFactor[0] = x[dim-3];
               polyFactor[1] = 1.0;
               for(i = 0; i < xPower-1; ++i)
               {
                  polyFactor[0] *= polyTemp[0];
                  polyFactor[1] *= polyTemp[1];
               }
               primeFactor[dim-2] = jacobi[*(index+1) + jacobiIndex] * xFactor * dxjacobi[*(index) + jacobiIndexPrev] * (-2.0) * x[dim-1] / SQUARE(x[dim-2]) +
                                    jacobi[*(index) + jacobiIndexPrev] *
                                    ( dxjacobi[*(index+1) + jacobiIndex] * (-2.0) / x[dim-3] * xFactor +
                                      jacobi[*(index+1) + jacobiIndex] * xPower * polyFactor[1] / polyFactor[0] );
            }
            else if(j == dim-1)
            {
               int d1 = basis_id[k*dim + dim-1];
               int d2 = basis_id[k*dim + dim-2];
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
                           ( dxjacobi[d1+jacobiIndex] * (-2.0) * polyFactor[1] + jacobi[d1+jacobiIndex] * xPower * polyFactor[3] );

               double v2 = jacobi[d1+jacobiIndex] * polyFactor[1] *
                           ( polyFactor[0] * dxjacobi[d2+jacobiIndexPrev] * 2.0 * x[1] / SQUARE(x[0]) +
                            jacobi[d2+jacobiIndexPrev] * polyFactor[2] * (-xPowerPrev) / polyFactor[4] );

               primeFactor[0] = v1 + v2;
            }

            if(j <= dim-2) primeFactor[dim-j-2] = 1.0;
            if(j == dim-2) primeFactor[0] = 1.0;

            if( (dim > 3) && (j >= 2) && (j != dim-1) )
            {
               int d1 = basis_id[k*dim+j];
               int d2 = basis_id[k*dim+j-1];
               int prev = dim-j-2;
               int cur = dim-j-1;
               int next = dim-j;

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
                           ( dxjacobi[d1+jacobiIndex] * (-2.0)/x[prev] * polyFactor[0] + jacobi[d1+jacobiIndex] * xPower * polyFactor[4] / x[prev] );

               double v2 = jacobi[d1+jacobiIndex] * polyFactor[0] *
                           ( polyFactor[1] * dxjacobi[d2+jacobiIndexPrev] * 2.0 * x[next]/SQUARE(x[cur]) +
                           jacobi[d2+jacobiIndexPrev] * polyFactor[2] * (-xPowerPrev) / polyFactor[3] );

               primeFactor[dim-1-j] = v1 + v2;
            }

            for(d = 0; d < dim; ++d) phiPrime[d*num_funs+k] *= primeFactor[d];

            dimPrev = j;
            xPowerPrev = xPower;
            jacobiIndexPrev = SQUARE(order)*dimPrev + order*xPowerPrev;
            xFactorPrev = xFactor;
         }
      }

      free(jacobi);
      free(dxjacobi);
      free(phi);
   }// end if dim >= 3

#ifdef QUAD_DEBUG_ON
   double *phi = (double *)malloc(num_funs*sizeof(double));
   PhiFunc phi_func = &PhiSimplex;
   int flag = finite_difference_test(basis_id, x, params, phiPrime, phi_func);
   if(flag == 1)
   {
//      PRINT_ERR("failed finite_difference_test", __LINE__, __FILE__);
//      printf("%lf, %lf\n", x[0], x[1]);
   }
   free(phi);
#endif

}// end PhiPrimeSimplex


/* PhiCubeSimplex
 * Generates orthogonal polynomial basis for the unit CUBESIMPLEX
 */
void PhiCubeSimplex(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi)
{
   assert(params->dims[0] >= 1 && params->dims[1] >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim1 = params->dims[0];
   int dim2 = params->dims[1];
   int dim = dim1+dim2;
   int deg = params->deg;
   int order = deg+1;
   int num_funs = BasisSize(dim, deg);

   double legendre[order*dim1];
   double dxlegendre[order*dim1];

   for(int d = 0; d < dim1; ++d)
      LegendrePoly(order, 2*x[d]-1, &legendre[d*order], &dxlegendre[d*order]);

   for(int k = 0; k < num_funs; ++k) phi[k] = 1.0;
   for(int k = 0; k < num_funs; ++k)
      for(int d = 0; d < dim1; ++d)
         phi[k] *= legendre[basis_id[k*dim+d]+order*d];

   double phiSimplex[num_funs];
   PhiSimplexPolyhedralTwo(basis_id, x, params, phiSimplex);
   for(int k = 0; k < num_funs; ++k) phi[k] *= phiSimplex[k];

}// end PhiCubeSimplex


/* PhiPrimeCubeSimplex
 * Generates derivatives of orthogonal polynomial basis for the unit CUBESIMPLEX
 */
void PhiPrimeCubeSimplex(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phiPrime)
{
   assert(params->dims[0] >= 1 && params->dims[1] >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int num_funs = params->num_funs;
   int dim1 = params->dims[0];
   int dim2 = params->dims[1];
   int dim = dim1+dim2;
   int order = params->deg+1;

   double legendre[order*dim1];
   double dxlegendre[order*dim1];
   double phiPrimeSimplex[num_funs*dim];

   for(int d = 0; d < dim1; ++d)
      LegendrePoly(order, 2*x[d]-1, &legendre[d*order], &dxlegendre[d*order]);

   for(int k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

   for(int d = 0; d < dim; ++d)
   {
      for(int j = 0; j < dim1; ++j)
      {
         if(j != d)
            for(int k = 0; k < num_funs; ++k)
               phiPrime[k+d*num_funs] *= legendre[basis_id[k*dim+j]+order*j];

         if(j == d)
            for(int k = 0; k < num_funs; ++k)
               phiPrime[k+d*num_funs] *= 2.0 * dxlegendre[basis_id[k*dim+j]+order*j];

      }
   }

   PhiPrimeSimplexPolyhedralTwo(basis_id, x, params, phiPrimeSimplex);
   for(int k = 0; k < num_funs*dim; ++k)
      phiPrime[k] *= phiPrimeSimplex[k];


#ifdef QUAD_DEBUG_ON
   PhiFunc phi_func = &PhiCubeSimplex;
   int flag = finite_difference_test(basis_id, x, params, phiPrime, phi_func);
   if(flag == 1)
      PRINT_ERR("failed finite_difference_test", __LINE__, __FILE__);
#endif

}// end PhiPrimeCubeSimplex


/* PhiSimplexSimplex
 * Generates orthogonal polynomial basis for the unit SIMPLEXSIMPLEX
 */
void PhiSimplexSimplex(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi)
{
   assert(params->dims[0] >= 2 && params->dims[1] >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int num_funs = params->num_funs;
   double *phiSimplex = (double *)malloc( SIZE_DOUBLE(num_funs) );

   for(int k = 0; k < num_funs; ++k) phi[k] = 1.0;

   PhiSimplexPolyhedralOne(basis_id, x, params, phiSimplex);
   for(int k = 0; k < num_funs; ++k)
      phi[k] *= phiSimplex[k];

   PhiSimplexPolyhedralTwo(basis_id, x, params, phiSimplex);
   for(int k = 0; k < num_funs; ++k)
      phi[k] *= phiSimplex[k];

   free(phiSimplex);

}


/* PhiPrimeSimplexSimplex
 * Generates derivatives of orthogonal polynomial basis for the unit SIMPLEXSIMPLEX
 */
void PhiPrimeSimplexSimplex(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phiPrime)
{
   assert(params->dims[0] >= 2 && params->dims[1] >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim = params->dim;
   int num_funs = params->num_funs;
   double *phiPrimeSimplex = (double *)malloc(SIZE_DOUBLE(num_funs*dim));

   for(int k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

   PhiPrimeSimplexPolyhedralOne(basis_id, x, params, phiPrimeSimplex);
   for(int k = 0; k < num_funs*dim; ++k)
      phiPrime[k] *= phiPrimeSimplex[k];

   PhiPrimeSimplexPolyhedralTwo(basis_id, x, params, phiPrimeSimplex);
   for(int k = 0; k < num_funs*dim; ++k)
      phiPrime[k] *= phiPrimeSimplex[k];

#ifdef QUAD_DEBUG_ON
   PhiFunc phi_func = &PhiSimplexSimplex;
   int flag = finite_difference_test(basis_id, x, params, phiPrime, phi_func);
   if(flag == 1)
      PRINT_ERR("failed finite_difference_test", __LINE__, __FILE__);
#endif

   free(phiPrimeSimplex);
}


/* PhiCubeSimplexSimplex
 * Generates orthogonal polynomial basis for the unit CUBESIMPLEXSIMPLEX
 */
void PhiCubeSimplexSimplex(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi)
{
//   assert(params->num_funs == BasisSize(params->dim, params->deg));
//   int k, d;
//   int p = params->deg;
//   int order = p+1;
//   int m = params->num_funs;
//   int dimSimplex1 = params->dims[0];
//   int dimSimplex2 = params->dims[1];
//   int dimCube = params->dims[2];
//   int twoDims = dimSimplex1+dimSimplex2;
//   int dim = params->dim;
//   quadParams *newParams = (quadParams *)malloc(sizeof(quadParams));
//   newParams->dims = (int *)malloc(3*sizeof(int));
//   SetParams(dim, 3, newParams->dims, p, newParams);
//   int size = BasisSize(dim, p); //number of basis functions in dim-dimensions
//   double phiSimplex[m];
//
//   for(k = 0; k < m; ++k)
//      phi[k] = 1.0;
//
//   PhiSimplexPolyhedralOne(basis_id, x, newParams, phiSimplex);
//   for(k = 0; k < m; ++k)
//      phi[k] *= phiSimplex[k];
//
//   PhiSimplexPolyhedralTwo(basis_id, x, newParams, phiSimplex);
//   for(k = 0; k < m; ++k)
//      phi[k] *= phiSimplex[k];
//
//   double* legendre = (double *)malloc(order*dimCube*sizeof(double));
//   double* dxlegendre = (double *)malloc(order*dimCube*sizeof(double));
//   for(d = 0; d < dimCube; ++d)
//      LegendrePoly(order, 2*x[twoDims+d]-1, &legendre[(d)*order], &dxlegendre[d*order]);
//
//   for(k = 0; k < m; ++k)
//   {
//      for(d = 0; d < dimCube; ++d)
//         phi[k] = phi[k] * legendre[basis_id[k*dim+twoDims+d]+order*d];
//   }
//   free(legendre);
//   free(dxlegendre);
//   free(newParams->dims);
//   free(newParams);
}//end PhiCubeSimplexSimplex


/* PhiPrimeCubeSimplexSimplex
 * Generates derivatives of orthogonal polynomial basis for the unit CUBESIMPLEXSIMPLEX
 */
void PhiPrimeCubeSimplexSimplex(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phiPrime)
{
   assert(params->num_funs == BasisSize(params->dim, params->deg));
   int j, k, d;
   int p = params->deg;
   int order = p+1;
   int m = params->num_funs;
   int dimSimplex1 = params->dims[0];
   int dimSimplex2 = params->dims[1];
   int dimCube = params->dims[2];
   int twoDims = dimSimplex1+dimSimplex2;
   int dim = params->dim;
   int size = BasisSize(dim, p);
   double phiPrimeSimplex[m*dim];

   for(k = 0; k < m*dim; ++k)
      phiPrime[k] = 1.0;

   PhiPrimeSimplexPolyhedralOne(basis_id, x, params, phiPrimeSimplex);
   for(k = 0; k < m*dim; ++k)
      phiPrime[k] = phiPrime[k] * phiPrimeSimplex[k];

   PhiPrimeSimplexPolyhedralTwo(basis_id, x, params, phiPrimeSimplex);
   for(k = 0; k < m*dim; ++k)
      phiPrime[k] = phiPrime[k] * phiPrimeSimplex[k];

   double* legendre = (double *)malloc(order*dimCube*sizeof(double));
   double* dxlegendre = (double *)malloc(order*dimCube*sizeof(double));
   for(d = 0; d < dimCube; ++d)
      LegendrePoly(order, 2*x[twoDims+d]-1, &legendre[d*order], &dxlegendre[d*order]);

   for(d = 0; d < dim; ++d)
   {
      for(k = 0; k < m; ++k)
      {
         for(j = twoDims; j < dim; ++j)
         {
            if(j != d)
               phiPrime[k+d*m] = phiPrime[k+d*m] * legendre[basis_id[k*dim+j] + order*(j-twoDims)];

            if(j == d)
               phiPrime[k+d*m] = 2*phiPrime[k+d*m] * dxlegendre[basis_id[k*dim+j] + order*(j-twoDims)];
         }
      }
   }

   free(legendre);
   free(dxlegendre);

#ifdef QUAD_DEBUG_ON
   void (*phi_func)(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi) = &PhiCubeSimplexSimplex;
   int flag = finite_difference_test(basis_id, x, params, phiPrime, phi_func);
#endif
}// end PhiPrimeCubeSimplexSimplex


/* PhiSimplexPolyhedralOne
 * Generates orthogonal polynomial basis for the first [0, dim1-1] coordinates
 * over simplex in dim-dimensional space.
 */
static void PhiSimplexPolyhedralOne(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi)
{
   assert(params->dims[0] >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim1 = params->dims[0];
   int dim = params->dim;
   int deg = params->deg;
   int num_funs = params->num_funs;
   int order = deg+1;

   if(dim1 == 2)
   {
      double legendre[order];
      double dxlegendre[order];
      double jacobi[SQUARE(order)];

      LegendrePoly(order, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
      for(int j = 0; j < order; ++j) JacobiPoly(order, 1.0-2*x[0], 2.0*j+1, 0.0, &jacobi[j*order]);
      for(int k = 0; k < num_funs; ++k) phi[k] = 1.0;

      int r, n;
      for(int k = 0; k < num_funs; ++k)
      {
         r = basis_id[dim*k];
         n = basis_id[dim*k+1];
         double xFactor = 1.0; for(int i = 0; i < r; ++i) xFactor *= x[0];
         phi[k] *= xFactor * jacobi[r*order+n] * legendre[r];
      }
      return;
   }


   // dim1 > 3
   int i, j, k, d;
   double legendre[order];
   double dxlegendre[order];

   double *jacobi = (double *) malloc( SQUARE(order)*(dim1-1)*size_double );
   LegendrePoly(order, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], legendre, dxlegendre);
   for(k = 0; k < num_funs; ++k) phi[k] = 1.0;

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
         JacobiPoly(order, jacobiX[dimCur], alpha, 0.0, &jacobi[order*(d-1)*deg+order*j]);
      }
   }

   for(k = 0; k < num_funs; ++k)
   {
      const int_fast8_t *index = &basis_id[k*dim];
      double phiTemp = 1.0;
      for(d = 1; d < dim1; ++d)
      {
         double xTemp = polyX[d-1];
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= xTemp;
         phiTemp *= jacobi[*(index+d) + order*(d-1)*deg + order*xPower] * xFactor;
      }
      phi[k] *= phiTemp*legendre[*index];
   }

   free(jacobi);
   return;
}// end PhiSimplexPolyhedralOne


/* PhiSimplexPolyhedralTwo
 * Generates orthogonal polynomial basis for coordinates between [dim1, dim2-1]
 * over a simplex in dim-dimensional space.
 */
static void PhiSimplexPolyhedralTwo(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi)
{
   assert(params->dims[1] >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim1 = params->dims[0];
   int dim2 = params->dims[1];
   int dimTwo = dim1+dim2;
   int dim = params->dim;
   int deg = params->deg;
   int num_funs = params->num_funs;
   int order = deg+1;

   if(dim2 == 2)
   {
      double legendre[order];
      double dxlegendre[order];
      double jacobi[SQUARE(order)];

      LegendrePoly(order, (2.0*x[dim1+1]-x[dim1])/x[dim1], legendre, dxlegendre);
      for(int j = 0; j < order; ++j) JacobiPoly(order, 1.0-2*x[dim1], 2.0*j+1, 0.0, &jacobi[j*order]);
      for(int k = 0; k < num_funs; ++k) phi[k] = 1.0;

      int r, n;
      for(int k = 0; k < num_funs; ++k)
      {
         r = basis_id[dim*k+dim1];
         n = basis_id[dim*k+dim1+1];
         double xFactor = 1.0; for(int i = 0; i < r; ++i) xFactor *= x[dim1];
         phi[k] *= xFactor*jacobi[r*order+n] * legendre[r];
      }
      return;
   }


   // dim1 > 3
   int i, j, k, d;
   double legendre[order];
   double dxlegendre[order];

   double *jacobi = (double *) malloc( SQUARE(order)*(dim2-1)*size_double );
   LegendrePoly(order, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], legendre, dxlegendre);
   for(k = 0; k < num_funs; ++k) phi[k] = 1.0;

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
         JacobiPoly(order, jacobiX[dimCur], alpha, 0.0, &jacobi[order*(d-1)*deg+order*j]);
      }
   }

   for(k = 0; k < num_funs; ++k)
   {
      const int_fast8_t *index = &basis_id[k*dim+dim1];
      double phiTemp = 1.0;
      for(d = 1; d < dim2; ++d)
      {
         double xTemp = polyX[d-1];
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= xTemp;
         phiTemp *= jacobi[*(index+d) + order*(d-1)*deg + order*xPower] * xFactor;
      }
      phi[k] *= phiTemp*legendre[*index];
   }
   free(jacobi);

   return;
}// end PhiSimplexPolyhedralTwo


/* PhiPrimeSimplexPolyhedralOne
 * Generates derivatives of orthogonal polynomial basis for the first [0, dim1-1] coordinates
 * over simplex in dimal-dimensional space.
 */
static void PhiPrimeSimplexPolyhedralOne(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phiPrime)
{
   assert(params->dims[0] >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim = params->dim;
   int dim1 = params->dims[0];
   int deg = params->deg;
   int order = deg+1;
   int num_funs = BasisSize(dim, deg);

   if(dim1 == 2)
   {
      int i, j, k, d;
      double factor1, factor2;
      double legendre[order];
      double dxlegendre[order];
      double jacobi[SQUARE(order)];
      double dxjacobi[SQUARE(order)];

      LegendrePoly(order, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
      for(i = 0; i < order; ++i) JacobiPoly(order, 1.0-2*x[0], 2.0*i+1, 0.0, &jacobi[i*order]);
      for(i = 0; i < order; ++i) JacobiPolyPrime(order, 1.0-2*x[0], 2.0*i+1, 0.0, &dxjacobi[i*order]);
      for(k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

      int power1, power2, index;
      for(int k = 0; k < num_funs; ++k)
      {
         power1 = basis_id[k*dim];
         power2 = basis_id[k*dim+1];
         index = power2+power1*order;
         factor1 = 1.0; for(j = 0; j < power1; ++j) factor1 *= x[0];
         factor2 = 1.0; for(j = 0; j < power1-1; ++j) factor2 *= x[0];

         phiPrime[k] = -2.0*x[1]/SQUARE(x[0]) * dxlegendre[power1] * factor1 * jacobi[index] +
                        legendre[power1] * ( power1 * factor2 * jacobi[index] + factor1 * dxjacobi[index] * (-2) );
         phiPrime[k+num_funs] = factor1 * jacobi[index] * dxlegendre[power1] * 2.0 / x[0];

         for(d = 2; d < dim; ++d)
            phiPrime[k+d*num_funs] *= factor1 * jacobi[power2+power1*order] * legendre[power1];
      }
   }


   if(dim1 >= 3)
   {
      int i, j, d;
      double jacobiX[dim1];
      double polyX[dim1];
      double primeFactor[dim1];

      double *jacobi = (double *)malloc( SQUARE(order)*(dim1)*size_double );
      double *dxjacobi = (double *)malloc( SQUARE(order)*(dim1)*size_double );
      LegendrePoly(order, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], &jacobi[0], &dxjacobi[0]);

      double *phi = (double *)malloc(num_funs*sizeof(double));
      PhiSimplexPolyhedralOne(basis_id, x, params, phi);
      for(int k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

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
            JacobiPoly(order, jacobiX[dimCur], alpha, 0.0, &jacobi[SQUARE(order)*d+order*j]);
            JacobiPolyPrime(order, jacobiX[dimCur], alpha, 0.0, &dxjacobi[SQUARE(order)*d+order*j]);
         }
      }

      const int_fast8_t *index;
      double polyTemp[4];
      double polyFactor[5];
      for(int k = 0; k < num_funs; ++k)
      {
         int dimPrev = 0;
         int xPowerPrev = 0;
         int jacobiIndexPrev = 0;
         double xFactorPrev = 0.0;
         index = &basis_id[k*dim];

         for(j = 0; j < dim1; ++j)
         {
            int xPower = 0; for(i = 0; i < j; ++i) xPower += *(index+i);

            polyTemp[0] = polyX[j];
            double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= polyTemp[0];
            int jacobiIndex = SQUARE(order)*j+order*xPower;
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
               phiPrime[d*num_funs+k] *= primeFactor[d];

            dimPrev = j;
            xPowerPrev = xPower;
            jacobiIndexPrev = SQUARE(order)*dimPrev + order*xPowerPrev;
            xFactorPrev = xFactor;
         }
      }

      for(d = dim1; d < dim; ++d)
         for(int k = 0; k < num_funs; ++k)
            phiPrime[k+d*num_funs] *= phi[k];

      free(jacobi);
      free(dxjacobi);
      free(phi);

   }// end if dim >= 3


}// end PhiPrimeSimplexPolyhedralOne


/* PhiPrimeSimplexPolyhedralTwo
 * Generates derivatives of orthogonal polynomial basis for coordinates between [dim1, dim2-1]
 * over simplex in dimal-dimensional space.
 */
static void PhiPrimeSimplexPolyhedralTwo(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phiPrime)
{
   assert(params->dims[1] >= 2);
   assert(params->num_funs == BasisSize(params->dim, params->deg));

   int dim = params->dim;
   int dim1 = params->dims[0];
   int dim2 = params->dims[1];
   int dimTwo = dim1+dim2;
   int deg = params->deg;
   int order = deg+1;
   int num_funs = BasisSize(dim, deg);


   if(dim2 == 2)
   {
      int i, j, k, d;
      double factor1, factor2;
      double legendre[order];
      double dxlegendre[order];
      double jacobi[SQUARE(order)];
      double dxjacobi[SQUARE(order)];

      LegendrePoly(order, (2.0*x[dim1+1]-x[dim1])/x[dim1], legendre, dxlegendre);
      for(i = 0; i < order; ++i) JacobiPoly(order, 1.0-2*x[dim1], 2.0*i+1, 0.0, &jacobi[i*order]);
      for(i = 0; i < order; ++i) JacobiPolyPrime(order, 1.0-2*x[dim1], 2.0*i+1, 0.0, &dxjacobi[i*order]);
      for(k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

      int power1, power2, index;
      for(int k = 0; k < num_funs; ++k)
      {
         power1 = basis_id[k*dim+dim1];
         power2 = basis_id[k*dim+dim1+1];
         index = power2+power1*order;
         factor1 = 1.0; for(j = 0; j < power1; ++j) factor1 *= x[dim1];
         factor2 = 1.0; for(j = 0; j < power1-1; ++j) factor2 *= x[dim1];

         phiPrime[k+dim1*num_funs] = -2.0*x[dim1+1]/SQUARE(x[dim1]) * dxlegendre[power1] * factor1 * jacobi[index] +
                                      legendre[power1] * ( power1 * factor2 * jacobi[index] + factor1 * dxjacobi[index] * (-2) );
         phiPrime[k+(dim1+1)*num_funs] = factor1 * jacobi[index] * dxlegendre[power1] * 2.0/x[dim1];

         for(d = 0; d < dim1; ++d)
            phiPrime[k+d*num_funs] *= factor1 * jacobi[power2+power1*order] * legendre[power1];
         for(d = dimTwo; d < dim; ++d)
            phiPrime[k+d*num_funs] *= factor1 * jacobi[power2+power1*order] * legendre[power1];
      }

      return;
   }


   if(dim2 >= 3)
   {
      int i, j, k, d;
      double jacobiX[dim2];
      double polyX[dim2];
      double primeFactor[dim2];

      double *jacobi = (double *)malloc( SQUARE(order)*(dim2)*size_double );
      double *dxjacobi = (double *)malloc( SQUARE(order)*(dim2)*size_double );
      LegendrePoly(order, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], &jacobi[0], &dxjacobi[0]);

      for(k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;
      double *phi = (double *)malloc( SIZE_DOUBLE(num_funs) );
      PhiSimplexPolyhedralTwo(basis_id, x, params, phi);

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
            JacobiPoly(order, jacobiX[dimCur], alpha, 0.0, &jacobi[SQUARE(order)*d+order*j]);
            JacobiPolyPrime(order, jacobiX[dimCur], alpha, 0.0, &dxjacobi[SQUARE(order)*d+order*j]);
         }
      }

      const int_fast8_t *index;
      double polyTemp[4];
      double polyFactor[5];
      for(int k = 0; k < num_funs; ++k)
      {
         int dimPrev = 0;
         int xPowerPrev = 0;
         int jacobiIndexPrev = 0;
         double xFactorPrev = 0.0;
         index = &basis_id[k*dim+dim1];

         for(j = 0; j < dim2; ++j)
         {
            int xPower = 0; for(i = 0; i < j; ++i) xPower += *(index+i);

            polyTemp[0] = polyX[j];
            double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= polyTemp[0];
            int jacobiIndex = SQUARE(order)*j+order*xPower;
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
               phiPrime[(d+dim1)*num_funs+k] *= primeFactor[d];

            dimPrev = j;
            xPowerPrev = xPower;
            jacobiIndexPrev = SQUARE(order)*dimPrev + order*xPowerPrev;
            xFactorPrev = xFactor;
         }
      }

      for(d = 0; d < dim1; ++d)
         for(int k = 0; k < num_funs; ++k)
            phiPrime[k+d*num_funs] *= phi[k];
      for(d = dimTwo; d < dim; ++d)
         for(int k = 0; k < num_funs; ++k)
            phiPrime[k+d*num_funs] *= phi[k];

      free(jacobi);
      free(dxjacobi);
      free(phi);

      return;
   }// end if dim >= 3

}// end PhiPrimeSimplexPolyhedralTwo


/* PhiMonomial
 * Generates monomial basis
 */
void PhiMonomial(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi)
{
   int deg = params->deg;
   int dim = params->dim;
   int num_funs = BasisSize(dim, deg);

   for(int k = 0; k < num_funs; ++k) phi[k] = 1.0;

   for(int k = 0; k < num_funs; ++k)
   {
      for(int d = 0; d < dim; ++d)
      {
         double product = 1.0;
         double factor = x[d];
         for(int r = 0; r < basis_id[k*dim+d]; ++r) product *= factor;

         phi[k] *= product;
      }
   }

}// end PhiCube


/* PhiPrimeCube
 * Generates derivatives of monomial basis for the unit cube
 */
void PhiPrimeMonomial(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phiPrime)
{
   int deg = params->deg;
   int dim = params->dim;
   int num_funs = BasisSize(dim, deg);

   for(int k = 0; k < num_funs*dim; ++k) phiPrime[k] = 1.0;

   for(int d = 0; d < dim; ++d)
   {
      for(int k = 0; k < num_funs; ++k)
      {
         double factor = 0.0, product = 0.0;

         // dimension < d
         for(int j = 0; j < d; ++j)
         {
            product = 1.0;
            factor = x[j];
            for(int r = 0; r < basis_id[k*dim+j]; ++r) product *= factor;
            phiPrime[k+d*num_funs] *= product;
         }

         // dimension > d
         for(int j = d+1; j < dim; ++j)
         {
            product = 1.0;
            factor = x[j];
            for(int r = 0; r < basis_id[k*dim+j]; ++r) product *= factor;
            phiPrime[k+d*num_funs] *= product;
         }

         // dimension == d
         if(basis_id[k*dim+d] > 0)
         {
            product = 1.0;
            factor = x[d];
            for(int r = 0; r < basis_id[k*dim+d]-1; ++r) product *= factor;
            phiPrime[k+d*num_funs] *= basis_id[k*dim+d] * product;
         }
         else if(basis_id[k*dim+d] == 0) phiPrime[k+d*num_funs] = 0.0;

      }
   }

}// end PhiPrimeMonomial


static int finite_difference_test(const int_fast8_t *basis_id, const double *x, const quadParams *params, const double *phi_prime,
                                  void(*phi_func)(const int_fast8_t *basis_id, const double *x, const quadParams *params, double *phi))
{
   int num_funs = params->num_funs;
   int dim = params->dim;

   double h = POW(10, -6);
   double tol = POW(10, -5);

   double *phi_approx = (double *)malloc( SIZE_DOUBLE(num_funs*dim) );
   double *phi_backw1 = (double *)malloc( SIZE_DOUBLE(num_funs) );
   double *phi_backw2 = (double *)malloc( SIZE_DOUBLE(num_funs) );
   double *phi_forw1  = (double *)malloc( SIZE_DOUBLE(num_funs) );
   double *phi_forw2  = (double *)malloc( SIZE_DOUBLE(num_funs) );

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

      phi_func(basis_id, x_backw1, params, phi_backw1);
      phi_func(basis_id, x_backw2, params, phi_backw2);
      phi_func(basis_id, x_forw1, params, phi_forw1);
      phi_func(basis_id, x_forw2, params, phi_forw2);

      for(int k = 0; k < num_funs; ++k)
         phi_approx[d*num_funs+k] = ( 1.0/12.0*phi_backw2[k] - 2.0/3.0*phi_backw1[k] + 2.0/3.0*phi_forw1[k]-1.0/12.0*phi_forw2[k] ) / h;

      for(int k = 0; k < num_funs; ++k)
         if( fabs( phi_prime[d*num_funs+k] - phi_approx[d*num_funs+k] ) > tol )
         {
            flag = 1;
            goto FREERETURN;
         }
   }

FREERETURN:
   free(phi_backw1);
   free(phi_backw2);
   free(phi_forw1);
   free(phi_forw2);
   free(phi_approx);

   return flag;
}


double  orthogonal_simplex_basis_test(const int_fast8_t *basis_id, const quadParams *params)
{
   int deg = params->deg;
   int dim = params->dim;
   int num_funs = params->num_funs;

   int n = floor(deg/2)+1 + floor(dim/2) + floor(deg/6) + 3; // make n large enough to get exact quadrature
   int dims_1D[1] = {1};
   quadrature *quad_1D = quadrature_init(n, 1, dims_1D, deg, INTERVAL);
   Jacobi(quad_1D->k, 0.0, 0.0, quad_1D->x, quad_1D->w);

   int N = POW(n, dim);
   quadrature *quad_S = quadrature_init(N, dim, params->dims, deg, SIMPLEX);
   GeneralizedNodesTensor( (const_quadrature *)quad_1D, quad_S );
   GeneralizedWeightsTensor( (const_quadrature *)quad_1D, quad_S );
   GeneralDuffy(quad_S);

   double *phi = (double *)malloc(num_funs*size_double);
   double *integrals = (double *)calloc(num_funs, size_double);

   for(int i = 0; i < N; ++i)
   {
      PhiSimplex(basis_id, &quad_S->x[dim*i], params, phi);
      for(int j = 0; j < num_funs; ++j)
         integrals[j] += phi[j]*quad_S->w[i];
   }

   double ie1 = 1.0; for(int i = 1; i < dim; ++i) ie1 = ie1/(i+1);

   double max_res = fabs(integrals[0] - ie1);
   double res = max_res;
   for(int j = 1; j < num_funs; ++j)
   {
      res = fabs(integrals[j]);
      max_res = MAX(max_res, res);
   }

   quadrature_free(quad_1D);
   quadrature_free(quad_S);
   free(phi);
   free(integrals);

   return max_res;
}
