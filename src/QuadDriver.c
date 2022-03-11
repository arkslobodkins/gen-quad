/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "QuadDriver.h"
#include "ComputeDomain.h"
#include "GENERAL_QUADRATURE.h"
#include "get_time.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mkl.h>

static void TimesToScreen(double total_time);
static void TimesToFile(DOMAIN_TYPE D, int deg, int dim, double total_time);

static bool checkIntervalParams(int argc, char **argv);
static bool checkCubeParams(int argc, char **argv);
static bool checkSimplexParams(int argc, char **argv);
static bool checkCubeSimplexParams(int argc, char **argv);
static bool checkSimplexSimplexParams(int argc, char **argv);
static int is_int(const char *number);
static int is_pos_int(char *);


void QuadDriver(int argc, char **argv)
{
   printf("\n");

   if( argv[1] == NULL )
   {
      fprintf(stderr, "No input parameters, exiting program.\n\n");
      exit(EXIT_FAILURE);
   }

   const char *SHAPE = argv[1];
   DOMAIN_TYPE D;
   if( !string_to_domain(SHAPE, &D) )
   {
      fprintf(stderr, "Domain is not specified properly. Acceptable parameters:\n"
            " 'INTERVAL', 'CUBE', 'SIMPLEX', 'CUBESIMPLEX', 'SIMPLEXSIMPLEX'.");
      fprintf(stderr, "Invalid input, exiting program.\n\n");
      exit(EXIT_FAILURE);
   }

   switch(D)
   {
   case(INTERVAL):
      if( !checkIntervalParams(argc, argv) ) {
         fprintf(stderr, "Invalid input, exiting program.\n\n");
         exit(EXIT_FAILURE);
      }
      break;

   case(CUBE):
      if( !checkCubeParams(argc, argv) ) {
         fprintf(stderr, "Invalid input, exiting program.\n\n");
         exit(EXIT_FAILURE);
      }
      break;
   case(SIMPLEX):
      if( !checkSimplexParams(argc, argv) ) {
         fprintf(stderr, "Invalid input, exiting program.\n\n");
         exit(EXIT_FAILURE);
      }
      break;
   case(CUBESIMPLEX):
      if( !checkCubeSimplexParams(argc, argv) ) {
         fprintf(stderr, "Invalid input, exiting program.\n\n");
         exit(EXIT_FAILURE);
      }
      break;
   case(SIMPLEXSIMPLEX):
      if( !checkSimplexSimplexParams(argc, argv) ) {
         printf("Invalid input, exiting program.\n\n");
         exit(EXIT_FAILURE);
      }
      break;
   }

   int deg = atoi(argv[2]);
   int counter = 0;
   int dims[3]; dims[0] = 0; dims[1] = 0; dims[2] = 0;
   int i = 3;
   //convert inputs to dimension arguments
   while( (argv[i] != NULL) && (counter < 3) )
   {
      dims[counter] = atoi(argv[i]);
      ++i ;
      ++counter;
   }


#ifdef QUAD_DEBUG_ON
   printf("DEBUG MODE ON\n\n");
#else
   printf("DEBUG MODE OFF\n\n");
#endif
#ifdef _OPENMP
   mkl_set_num_threads(1);
   printf("OPENMP enanbled with %i threads\n\n", omp_get_max_threads());
#endif

   double time_start_wall = get_cur_time();
   switch(D)
   {
      case INTERVAL:
         ComputeInterval(deg);
         break;
      case CUBE:
         ComputeCube(deg, dims[0]);
         break;
      case SIMPLEX:
         ComputeSimplex(deg, dims[0]);
         break;
      case CUBESIMPLEX:
         ComputeCubeSimplex(deg, dims[0], dims[1]);
         break;
      case SIMPLEXSIMPLEX:
         ComputeSimplexSimplex(deg, dims[0], dims[1]);
         break;
   }
   double time_end_wall = get_cur_time();
   double total_time = time_end_wall - time_start_wall;
   int dim = dims[0]+dims[1]+dims[2];
   TimesToScreen(total_time);
   TimesToFile(D, deg, dim, total_time);

}



static void TimesToScreen(double total_time)
{
   extern double LSQ_TIME;
   extern double JACOBIAN_TIME;
   extern double FUNCTION_TIME;
   extern double PREDICTOR_TIME;
   double MAIN_ROUTINES_TIME = LSQ_TIME + JACOBIAN_TIME + FUNCTION_TIME + PREDICTOR_TIME;
   printf("wall clock time for least squares routine in LeastSquaresNewton = %le\n", LSQ_TIME);
   printf("wall clock time for GetJacobian routine = %le\n", JACOBIAN_TIME);
   printf("wall clock time for GetFunction routine = %le\n", FUNCTION_TIME);
   printf("wall clock time for predictor routine   = %le\n", PREDICTOR_TIME);
   printf("\n");
   printf("main routines took %f percent of the total computation\n", 100.0*MAIN_ROUTINES_TIME/total_time);
   printf("total wall clock time = %le\n", total_time);
}

static void TimesToFile(DOMAIN_TYPE D, int deg, int dim, double total_time)
{
   FILE *file;
   char str[80];
   const char *shape = get_domain_string(D);
   sprintf(str, "../results/times_%s_dim%i_deg%i.txt", shape, dim, deg);
   file = fopen(str, "w");

   extern double LSQ_TIME;
   extern double JACOBIAN_TIME;
   extern double FUNCTION_TIME;
   extern double PREDICTOR_TIME;
   double MAIN_ROUTINES_TIME = LSQ_TIME + JACOBIAN_TIME + FUNCTION_TIME + PREDICTOR_TIME;
   fprintf(file, "wall clock time for least squares routine in LeastSquaresNewton = %le\n", LSQ_TIME);
   fprintf(file, "wall clock time for GetJacobian routine = %le\n", JACOBIAN_TIME);
   fprintf(file, "wall clock time for GetFunction routine = %le\n", FUNCTION_TIME);
   fprintf(file, "wall clock time for predictor routines  = %le\n", PREDICTOR_TIME);
   fprintf(file,"\n");
   fprintf(file,"main routines took %f percent of the total computation\n", 100.0*MAIN_ROUTINES_TIME/total_time);
   fprintf(file, "total wall clock time = %le\n", total_time);

   fclose(file);
}



static bool checkIntervalParams(int argc, char **argv)
{
   if( argc < 3 )
   {
      fprintf(stderr, "Not enough input parameters for INTERVAL.\n");
      return false;
   }
   else
   {
      if( !is_pos_int(argv[2]) )
      {
         fprintf(stderr, "Degree of precision should be integer >= 1.\n");
         return false;
      }
   }

   if( argc > 3 )
   {
      fprintf(stderr, "Received redundant parameters for INTERVAL.\n");
      return false;
   }

   return true;
}

static bool checkCubeParams(int argc, char **argv)
{
   if( argc < 4 )
   {
      fprintf(stderr, "Not enough input parameters for CUBE.\n");
      return false;
   }
   else
   {
      int test = 1;
      if( !is_pos_int(argv[2]) )
      {
         fprintf(stderr, "Degree of precision should be integer >= 1.\n");
         test = 0;
      }
      if( !is_pos_int(argv[3]) || (atoi(argv[3]) < 2 ) )
      {
         fprintf(stderr, "Dimension for CUBE should be integer >= 2.\n");
         test = 0;
      }
      if( test == 0 ) return false;
   }

   if( argc > 4 )
   {
      fprintf(stderr, "Received redundant parameters for CUBE.\n");
      return false;
   }

   return true;
}

static bool checkSimplexParams(int argc, char **argv)
{
   if( argc < 4 )
   {
      fprintf(stderr, "Not enough input parameters for SIMPLEX.\n");
      return false;
   }
   else
   {
      int test = 1;
      if( !is_pos_int(argv[2]) )
      {
         fprintf(stderr, "Degree of precision should be integer >= 1.\n");
         test = 0;
      }
      if( !is_pos_int(argv[3]) || (atoi(argv[3]) < 2) )
      {
         fprintf(stderr, "Dimension for SIMPLEX should be integer >= 2.\n");
         test = 0;
      }
      if( test == 0 ) return false;
   }

   if( argc > 4 )
   {
      fprintf(stderr, "Received redundant parameters for SIMPLEX.\n");
      return false;
   }

   return true;
}

static bool checkCubeSimplexParams(int argc, char **argv)
{
   if( argc < 5 )
   {
      fprintf(stderr, "Not enough input parameters for CUBESIMPLEX.\n");
      return false;
   }
   else
   {
      int test = 1;
      if( !is_pos_int(argv[2]) )
      {
         fprintf(stderr, "Degree of precision should be integer >= 1.\n");
         test = 0;
      }
      if( !is_pos_int(argv[3]) )
      {
         fprintf(stderr, "First dimension for CUBESIMPLEX should be integer >= 1.\n");
         test = 0;
      }
      if( !is_pos_int(argv[4]) || (atoi(argv[4]) < 2) )
      {
         fprintf(stderr, "Second dimension for CUBESIMPLEX should be integer >= 2.\n");
         test = 0;
      }
      if( test == 0 ) return false;
   }

   if( argc > 5 )
   {
      fprintf(stderr, "Received redundant parameters for CUBESIMPLEX.\n");
      return false;
   }

   return true;
}

static bool checkSimplexSimplexParams(int argc, char **argv)
{
   if( argc < 5 )
   {
      fprintf(stderr, "Not enough input parameters for SIMPLEXSIMPLEX, exiting program.\n");
      return false;
   }
   else
   {
      int test = 1;
      if( !is_pos_int(argv[2]) )
      {
         fprintf(stderr, "Degree of precision should be integer >= 1.\n");
         test = 0;
      }
      if( !is_pos_int(argv[3]) || (atoi(argv[3]) < 2) )
      {
         fprintf(stderr, "First dimension for SIMPLEXSIMPLEX should be integer >= 2.\n");
         test = 0;
      }
      if( !is_pos_int(argv[4]) || (atoi(argv[4]) < 2) )
      {
         fprintf(stderr, "Second dimension for SIMPLEXSIMPLEX should be integer >= 2.\n");
         test = 0;
      }
      if( test == 0 ) return false;
   }

   if( argc > 5 )
   {
      fprintf(stderr, "Received redundant parameters for SIMPLEXSIMPLEX.\n");
      return false;
   }

   return true;
}

static int is_int(const char *number)
{
   bool neg = false;

   //checking for negative numbers
   if( number[0] == '-' )
   {
      neg = true;
      for(int i = 1; number[i] != 0; ++i)
         if( !isdigit(number[i]) )
            return 0;
   }
   else
   {
      for(int i = 0; number[i] != 0; ++i)
         if( !isdigit(number[i]) )
            return 0;
   }

   if(neg) return -1;
   else    return 1;
}

static int is_pos_int(char *number)
{
   if( (is_int(number) != 1) || (atoi(number) == 0) )
      return 0;
   return 1;
}
