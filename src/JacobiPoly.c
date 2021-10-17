/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "JacobiPoly.h"

#define SQUARE(x) ((x)*(x))
#define ADD(x, y) (x)+(y)
#define SUB(x, y) (x)-(y)

void JacobiPoly(int order, double x, double alpha, double beta, double *p)
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
