#include "Print.h"
#include "GENERAL_QUADRATURE.h"

#include <stdlib.h>
#include <stdio.h>

void PrintNodes(const_quadrature *q, const char *name)
{
   int dim = q->params->dim;

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
   int dim = q->params->dim;
   printf("%i. node and weight for %s\n", id, name);

   printf("w[%i] = %.16f ", id, q->w[id]);
   for(int j = 0; j < dim; ++j)
      printf("x[%i][%i] = %.16lf" "  ", id, j, q->x[dim*id+j]);

   printf("\n\n");
}


void PrintNodesAndWeights(const_quadrature *q, const char *name)
{
   int dim = q->params->dim;

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


void PrintNodeInfo(int iters, double error_norm, const _DomainFuncs dom_funcs, const_quadrature *q, const char *name)
{
   printf("quadrature info for %s during iteration of LSQ\n", name);
   printf("iterations = %i error_norm = %.16le InDomain = %i\n\n", iters, error_norm, dom_funcs.inDomain(q));
}


void PrintElimInfo(int dim, int num_nodes, int opt)
{
   printf("dimension = %i, current number of nodes = %i, optimal = %i\n", dim, num_nodes, opt);

}


void print_int(int x, const char *name)
{
   printf("%s = %i\n", name, x);
}

void print_bool(bool x, const char *name)
{
   printf("%s = %d\n", name, x);
}

void print_float(float x, const char *name)
{
   printf("%s = %f\n", name, x);
}

void print_double(double x, const char *name)
{
   printf("%s = %.16e\n", name, x);
}

void print_string(const char *x, int format)
{
   if (format == 0)
      printf("%s\n", x);
   if (format == 1)
      printf("\n%s\n\n", x);
   if (format == 2)
      printf("\n\n%s\n\n\n", x);
   else if(format == 3)
      printf("\n%s\n", x);
   else if(format == 4)
      printf("%s\n\n", x);
}


void print_int_err(int x, const char *name, int line_num, const char *file_name)
{
   fprintf(stderr, "error on line %d in %s, %s = %i\n", line_num, file_name, name, x);
}

void print_bool_err(bool x, const char *name, int line_num, const char *file_name)
{
   fprintf(stderr, "error on line %d in %s, %s = %d\n", line_num, file_name, name, x);
}

void print_float_err(float x, const char *name, int line_num, const char *file_name)
{
   fprintf(stderr, "error on line %d in %s, %s = %f\n", line_num, file_name, name, x);
}

void print_double_err(double x, const char *name, int line_num, const char *file_name)
{
   fprintf(stderr, "error on line %d in %s, %s = %.16e\n", line_num, file_name, name, x);
}

void print_string_err(const char *x, int format, int line_num, const char *file_name)
{
   if (format == 0)
      fprintf(stderr, "error on line %d in %s, %s\n", line_num, file_name, x);
   if (format == 1)
      fprintf(stderr, "\nerror on line %d in %s, %s\n\n", line_num, file_name, x);
   if (format == 2)
      fprintf(stderr, "\n\nerror on line %d in %s, %s\n\n\n", line_num, file_name, x);
   else if(format == 3)
      fprintf(stderr, "\nerror on line %d in %s, %s\n", line_num, file_name, x);
   else if(format == 4)
      fprintf(stderr, "error on line %d in %s, %s \n\n", line_num, file_name, x);
}


