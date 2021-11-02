#include "Print.h"
#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#include <stdlib.h>
#include <stdio.h>

void PrintNodes(const_quadrature *q, const char *name)
{
   int dim = q->dim;

   printf("nodes for %s\n", name);
   for(int i = 0; i < q->k; ++i)
   {
      for(int j = 0; j < dim; ++j)
      {
         printf("x[%i][%i] = %.16lf" "  ", i, j, q->x[dim*i+j]);
      }
      printf("\n");
   }
   printf("\n");
}

void PrintNodeAndWeight(int id, const_quadrature *q, const char *name)
{
   int dim = q->dim;
   printf("%i. node and weight for %s\n", id, name);

   printf("w[%i] = %.16f ", id, q->w[id]);
   for(int j = 0; j < dim; ++j)
      printf("x[%i][%i] = %.16lf" "  ", id, j, q->x[dim*id+j]);

   printf("\n\n");
}

void PrintNodesAndWeights(const_quadrature *q, const char *name)
{
   int dim = q->dim;

   printf("nodes and weights for %s\n", name);
   for(int i = 0; i < q->k; ++i)
   {
      printf("w[%i] = %.6f ", i, q->w[i]);
      for(int j = 0; j < dim; ++j)
      {
         printf("x[%i][%i] = %.12lf" "  ", i, j, q->x[dim*i+j]);
      }
      printf("\n");
   }
   printf("\n");
}

void PrintNodeInfo(int iters, double error_norm, const_quadrature *q, const char *name)
{
   printf("quadrature info for %s during iteration of LSQ\n", name);
   printf("iterations = %i error_norm = %.16le InDomain = %i\n\n", iters, error_norm, QuadInDomain(q));
}

void PrintElimInfo(int dim, int num_nodes, int opt, double opt_factor)
{
   printf("dimension = %i, current number of nodes = %i, optimal = %i, efficiency = %f\n", dim, num_nodes, opt, opt_factor);

}

void PrintInt(int x, const char *name)
{
   printf("%s = %i\n", name, x);
}

void PrintBool(bool x, const char *name)
{
   printf("%s = %d\n", name, x);
}

void PrintFloat(float x, const char *name)
{
   printf("%s = %f\n", name, x);
}

void PrintDouble(double x, const char *name)
{
   printf("%s = %.16e\n", name, x);
}

void Print(const char *x)
{
   printf("%s\n", x);
}

void PRINT_ERR(const char *x, int line_num, const char *file_name)
{
      fprintf(stderr, "\nerror on line %d in %s, %s\n\n", line_num, file_name, x);
}


