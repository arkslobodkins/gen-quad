#include "GENERAL_QUADRATURE.h"
#include <strings.h>
#include <assert.h>
#include <omp.h>

char *get_domain_string(DOMAIN_TYPE D)
{
   switch(D)
   {
      case INTERVAL:
         return (char *)"LINE";
      case CUBE:
         return (char *)"CUBE";
      case SIMPLEX:
         return (char *)"SIMPLEX";
      case CUBESIMPLEX:
         return (char *)"CUBESIMPLEX";
      case SIMPLEXSIMPLEX:
         return  (char *)"SIMPLEXSIMPLEX";
      default:
         return (char *)"0";
   }
}

bool string_to_domain(const char *shape, DOMAIN_TYPE *D)
{
   if( strcasecmp("INTERVAL", shape) == 0 )                *D = INTERVAL;
   else if( strcasecmp("CUBE", shape) == 0 )               *D = CUBE;
   else if( strcasecmp("SIMPLEX", shape) == 0 )            *D = SIMPLEX;
   else if( strcasecmp("CUBESIMPLEX", shape) == 0 )        *D = CUBESIMPLEX;
   else if (strcasecmp("SIMPLEXSIMPLEX", shape) == 0 )     *D = SIMPLEXSIMPLEX;
   else                                                     return false;

   return true;
}

int IntPower(int x, int power)
{
   int result = 1;
   for(int i = 0; i < power; ++i)
      result *= x;

   return result;
}

double DoubleIntPower(double x, int power)
{
   double result = 1.0;
   for(int i = 0; i < power; ++i)
      result *= x;

   return result;
}

int factorial(int n)
{
    if (n == 0) return 1;
    return n * factorial(n - 1);
}

int binomial(int k, int n)
{
   return factorial(n)/(factorial(k)*factorial(n-k));
}

double expNDim(int dim, double x[])
{
   double expVal = 1.0;
   for(int i = 0; i < dim; ++i)
      expVal *= expl(x[i]);

   return expVal;
}

long double expIntegral1D(long double c)
{
   return 1.0L/c * (expl(c)-1.0L);
}

double expIntegralNDimCube(int dim)
{
   double expI = 1.0;
   for(int i = 0; i < dim; ++i)
      expI *= expIntegral1D(1.0);

   return expI;
}

double expIntegralNDimSimplex(int dim)
{
   assert(dim >= 2);

   double expI = 0.0;
   double sign = 1.0;
   for(int i = dim; i >= 0; --i) {
         expI += (double)binomial(i, dim)*expl(i)*sign;
      sign *= -1.0;
   }
   expI /= factorial(dim);

   return expI;
}

#ifdef _OPENMP
GQ_BOOL OMP_CONDITION(int deg, int dim)
{
   if(dim == 3  && deg >= 12)      return GQ_TRUE;
   else if(dim == 4  && deg  >= 7) return GQ_TRUE;
   else if(dim == 5  && deg  >= 6) return GQ_TRUE;
   else if(dim == 6  && deg  >= 5) return GQ_TRUE;
   else if(dim > 6)                return GQ_TRUE;
   else                            return GQ_FALSE;
}

GQ_BOOL PLASMA_CONDITION()
{
   if(omp_get_max_threads() > 4) return GQ_TRUE;
   else                          return GQ_FALSE;
}
#endif
