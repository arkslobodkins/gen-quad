/* Evaluate the Legendre polynomials and their first derivatives
 * by the recursion formula
 */

#include "LegendrePoly.h"

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
      fac2 =k/(k+1.0);
      p[k+1] = fac1*x*p[k] - fac2*p[k-1];
      dp[k+1] = fac1*(p[k] + x*dp[k]) - fac2*dp[k-1];
   }

}

