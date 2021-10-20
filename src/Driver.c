/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GENERAL_QUADRATURE.h"
#include "ComputeDomain.h"
#include "ImplementDomain.h"
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

int is_int(const char *number);
int is_pos_int(int, char *);

// main
// Receives dimension sizes and degree of precision as command line arguments.
// Performs tests to test appropriateness of the input arguments. If inputs are
// valid, it proceeds to recursive initial guess procedure and Node Elimination algorithm.
int main(int argc, char *argv[])
{
   printf("\n");

   const char *SHAPE;
   int deg, dims[3];

   if(argv[1] == NULL)
   {
      fprintf(stderr, "No input parameters specified, exiting program.\n\n");
      exit(EXIT_FAILURE);
   }
   else
      SHAPE = argv[1];

   DOMAIN_TYPE D;
   if( strcasecmp("INTERVAL", SHAPE) == 0 )  D = INTERVAL;
   else if( strcasecmp("CUBE", SHAPE) == 0 )  D = CUBE;
   else if( strcasecmp("SIMPLEX", SHAPE) == 0 )  D = SIMPLEX;
   else if( strcasecmp("CUBESIMPLEX", SHAPE) == 0 )  D = CUBESIMPLEX;
   else if (strcasecmp("SIMPLEXSIMPLEX", SHAPE) == 0 )  D = SIMPLEXSIMPLEX;
   else if (strcasecmp("CUBESIMPLEXSIMPLEX", SHAPE) == 0 )  D = CUBESIMPLEXSIMPLEX;
   else
   {
      fprintf(stderr, "Domain is not specified properly. Acceptable parameters:\n"
            " 'INTERVAL', 'CUBE', 'SIMPLEX', 'CUBESIMPLEX', 'SIMPLEXSIMPLEX', 'CUBESIMPLEXSIMPLEX.");

      fprintf(stderr, "Exiting program.\n\n");
      exit(EXIT_FAILURE);
   }

   int test;
   switch(D)
   {

      case(INTERVAL):
         if( argv[2] == NULL )
         {
            fprintf(stderr, "Not enough input parameters for %s, exiting program.\n\n", SHAPE);
            exit(EXIT_FAILURE);
         }
         else
         {
            test = 1;
            deg = atoi(argv[2]);
            if( !is_pos_int(2, argv[2]) )
            {
               fprintf(stderr, "Degree of precision should be integer >= 1.\n");
               test = 0;
            }
            if(test == 0)
            {
               fprintf(stderr, "Exiting program.\n\n");
               exit(EXIT_FAILURE);
            }
         }

         if(argv[3] != NULL)
         {
            fprintf(stderr, "Received redundant parameters for %s, exiting program.\n\n", SHAPE);
            exit(EXIT_FAILURE);
         }
         break;



      case(CUBE):
         if( (argv[2] == NULL) || (argv[3] == NULL) )
         {
            fprintf(stderr, "Not enough input parameters for %s, exiting program.\n\n", SHAPE);
            exit(EXIT_FAILURE);
         }
         else
         {
            test = 1;
            deg = atoi(argv[2]);
            if( !is_pos_int(2, argv[2]) )
            {
               fprintf(stderr, "Degree of precision should be integer >= 1.\n");
               test = 0;
            }
            if( !is_pos_int(3, argv[3]) || (atoi(argv[3]) < 2 ) )
            {
               fprintf(stderr, "Dimension for %s should be integer >= 2.\n", SHAPE);
               test = 0;
            }
            if(test == 0)
            {
               fprintf(stderr, "Exiting program.\n\n");
               exit(EXIT_FAILURE);
            }
         }

         if(argv[4] != NULL)
         {
            fprintf(stderr, "Received redundant parameters for %s, exiting program.\n\n", SHAPE);
            exit(EXIT_FAILURE);
         }
         break;



      case(SIMPLEX):
         if( (argv[2] == NULL) || (argv[3] == NULL) )
         {
            fprintf(stderr, "Not enough input parameters for %s, exiting program.\n\n", SHAPE);
            exit(EXIT_FAILURE);
         }
         else
         {
            deg = atoi(argv[2]);
            test = 1;
            if ( !is_pos_int(2, argv[2]) )
            {
               fprintf(stderr, "Degree of precision should be integer >= 1.\n");
               test = 0;
            }
            if ( !is_pos_int(3, argv[3]) || (atoi(argv[3]) < 2) )
            {
               fprintf(stderr, "Dimension for %s should be integer >= 2.\n", SHAPE);
               test = 0;
            }
            if(test == 0)
            {
               fprintf(stderr, "Exiting program.\n\n");
               exit(EXIT_FAILURE);
            }
         }

         if(argv[4] != NULL)
         {
            fprintf(stderr, "Received redundant parameters for %s, exiting program.\n\n", SHAPE);
               exit(EXIT_FAILURE);
         }
         break;



      case(CUBESIMPLEX):
         if( (argv[2] == NULL) || (argv[3] == NULL) || (argv[4] == NULL) )
         {
            fprintf(stderr, "Not enough input parameters for %s, exiting program.\n\n", SHAPE);
            exit(EXIT_FAILURE);
         }
         else
         {
            deg = atoi(argv[2]);
            test = 1;
            if ( !is_pos_int(2, argv[2]) )
            {
               fprintf(stderr, "Degree of precision should be integer >= 1.\n");
               test = 0;
            }
            if ( !is_pos_int(3, argv[3]) )
            {
               fprintf(stderr, "First dimension for %s should be integer >= 1.\n", SHAPE);
               test = 0;
            }
            if ( !is_pos_int(4, argv[4]) || (atoi(argv[4]) < 2) )
            {
               fprintf(stderr, "Second dimension for %s should be integer >= 2.\n", SHAPE);
               test = 0;
            }
            if(test == 0)
            {
               fprintf(stderr, "Exiting program.\n\n");
               exit(EXIT_FAILURE);
            }
         }
         if(argv[5] != NULL)
         {
            fprintf(stderr, "Received redundant parameters for %s, exiting program.\n\n", SHAPE);
               exit(EXIT_FAILURE);
         }
         break;


      case(SIMPLEXSIMPLEX):
         if( (argv[2] == NULL) || (argv[3] == NULL) || (argv[4] == NULL) )
         {
            fprintf(stderr, "Not enough input parameters for %s, exiting program.\n\n", SHAPE);
            exit(EXIT_FAILURE);
         }
         else
         {
            deg = atoi(argv[2]);
            test = 1;
            if (!is_pos_int(2, argv[2]))
            {
               fprintf(stderr, "Degree of precision should be integer >= 1.\n");
               test = 0;
            }
            if ( !is_pos_int(3, argv[3]) || (atoi(argv[3]) < 2) )
            {
               fprintf(stderr, "First dimension for %s should be integer >= 2.\n", SHAPE);
               test = 0;
            }
            if ( !is_pos_int(4, argv[4]) || (atoi(argv[4]) < 2) )
            {
               fprintf(stderr, "Second dimension for %s should be integer >= 2.\n", SHAPE);
               test = 0;
            }
            if(test == 0)
            {
               fprintf(stderr, "Exiting program.\n\n");
               exit(EXIT_FAILURE);
            }
         }
         if(argv[5] != NULL)
         {
            fprintf(stderr, "Received redundant parameters for %s, exiting program.\n\n", SHAPE);
            exit(EXIT_FAILURE);
         }
         break;


      case(CUBESIMPLEXSIMPLEX):
         if( (argv[2] == NULL) || (argv[3] == NULL) || (argv[4] == NULL) || (argv[5] == NULL) )
         {
            fprintf(stderr, "Not enough input parameters for %s, exiting program.\n", SHAPE);
            fprintf(stderr, "\n");
            exit(EXIT_FAILURE);
         }
         else
         {
            deg = atoi(argv[2]);
            test = 1;
            if ( !is_pos_int(2, argv[2]) )
            {
               fprintf(stderr, "Degree of precision should be integer >= 1.\n");
               test = 0;
            }
            if ( !is_pos_int(3, argv[3]) )
            {
               fprintf(stderr, "First dimension for %s should be integer >= 1.\n", SHAPE);
               test = 0;
            }
            if ( !is_pos_int(4, argv[4]) || (atoi(argv[4]) < 2) )
            {
               fprintf(stderr, "Second dimension for %s should be integer >= 2.\n", SHAPE);
               test = 0;
            }
            if ( !is_pos_int(5, argv[5]) || (atoi(argv[5]) < 2) )
            {
               fprintf(stderr, "Third dimension for %s should be integer >= 2.\n", SHAPE);
               test = 0;
            }
            if (test == 0)
            {
               fprintf(stderr, "Exiting program.\n\n");
               exit(EXIT_FAILURE);
            }
         }
         if(argv[6] != NULL)
         {
            fprintf(stderr, "Received redundant parameters for %s, exiting program.\n\n", SHAPE);
               exit(EXIT_FAILURE);
         }
         break;
   }

   int counter = 0;
   int i = 3;
   //convert inputs to dimension arguments
   while( (argv[i] != NULL) && (counter < 3) )
   {
      dims[counter] = atoi(argv[i]);
      ++i ;
      ++counter;
   }

   time_t start = clock();
   ImplementDomain(deg, dims, D);
   time_t end = clock();
   double total = (double)(end - start)/CLOCKS_PER_SEC;
   extern double LSQ_TIME;
   printf("total time for LAPACK routine in LeastSquaresNewton = %le\n", LSQ_TIME);
   printf("runtime = %le\n", total);

   return EXIT_SUCCESS;
} //end main


int is_int(const char *number)
{
   int i;
   int neg = 0;

   //checking for negative numbers
   if (number[0] == '-')
   {
      i = 1;
      neg = 1;

      for (i = 1; number[i] != 0; ++i)
      {
         if ( (number[i] > '9') || (number[i] < '0') )
            if ( !isdigit(number[i]) )
               return 0;
      }
   }
   else
   {
      for (i = 0; number[i] != 0; ++i)
      {
         if ( (number[i] > '9') || (number[i] < '0') )
            if ( !isdigit(number[i]) )
               return 0;
      }
   }

   if(neg == 1)
      return -1;

   return 1;
}


int is_pos_int(int index, char *number)
{
   int pos_int = 1;
   if( (is_int(number) != 1) && (is_int(number) != -1) )
      pos_int = 0;
   else if( (is_int(number) == -1) || (atoi(number) == 0) )
      pos_int = 0;
   return pos_int;
}
