#include "GENERAL_QUADRATURE.h"
#include <strings.h>

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

#ifdef _OPENMP
GQ_BOOL OMP_CONDITION(int deg, int dim)
{
   if(dim == 3  && deg >= 12)      return GQ_TRUE;
   else if(dim == 4  && deg  >= 6) return GQ_TRUE;
   else if(dim == 5  && deg  >= 5) return GQ_TRUE;
   else if(dim > 5)                return GQ_TRUE;
   else                            return GQ_FALSE;
}
#endif
