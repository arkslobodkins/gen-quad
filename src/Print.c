#include "Print.h"
#include "Quadrature.h"

#include <stdlib.h>
#include <stdio.h>

void PrintNodes(const quadrature *q, const char *name)
{
   int dim = q->dim;

   printf("nodes for %s\n", name);
   for(int i = 0; i < q->num_nodes; ++i)
   {
      for(int j = 0; j < dim; ++j)
      {
         printf("x[%i][%i] = %.16lf" "  ", i, j, q->x[dim*i+j]);
      }
      printf("\n");
   }
   printf("\n");
}

void PrintNodeAndWeight(int id, const quadrature *q, const char *name)
{
   int dim = q->dim;
   printf("%i. node and weight for %s\n", id, name);

   printf("w[%i] = %.16f ", id, q->w[id]);
   for(int j = 0; j < dim; ++j)
      printf("x[%i][%i] = %.16lf" "  ", id, j, q->x[dim*id+j]);

   printf("\n\n");
}

void PrintNodesAndWeights(const quadrature *q, const char *name)
{
   int dim = q->dim;
   printf("nodes and weights for %s\n", name);
   for(int i = 0; i < q->num_nodes; ++i)
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

void PrintNodesAndWeightsToFile(const quadrature *q, const char *name, FILE *file)
{
   int dim = q->dim;
   fprintf(file, "nodes and weights for %s\n", name);
   for(int i = 0; i < q->num_nodes; ++i)
   {
      fprintf(file, "w[%i] = %.6f ", i, q->w[i]);
      for(int j = 0; j < dim; ++j)
      {
         fprintf(file, "x[%i][%i] = %.12lf" "  ", i, j, q->x[dim*i+j]);
      }
      fprintf(file, "\n");
   }
   fprintf(file, "\n");
}

void PrintQuadLSQInfo(int iters, double error_norm, const quadrature *q, const char *name)
{
   printf("quadrature info for %s during iteration of LSQ\n", name);
   printf("iterations = %i error_norm = %.16le InDomain = %i\n\n",
         iters, error_norm, QuadInDomain(q));
}

void PrintElimInfo(int dim, int num_nodes, int opt, double opt_factor)
{
   printf("dimension = %i, current number of nodes = %i, optimal = %i, efficiency = %f\n",
         dim, num_nodes, opt, opt_factor);
}

void PrintInt(int x, const char *name)
{
   printf("%s = %i\n", name, x);
}

void PrintInts(int *x, int n, const char *name)
{
   for(int i = 0; i < n; ++i)
      printf("%s[%i] = %i\n", name, i, x[i]);
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

void PrintDoubles(double *x, int n, const char *name)
{
   for(int i = 0; i < n; ++i)
      printf("%s[%i] = %.16e\n", name, i, x[i]);
}

void Print(const char *x)
{
   printf("%s\n", x);
}

void PRINT_ERR(const char *x, int line_num, const char *file_name)
{
      fprintf(stderr, "\nerror on line %d in %s, %s\n\n", line_num, file_name, x);
}

void PRINT_WARN(const char *x, int line_num, const char *file_name)
{
      fprintf(stderr, "\nwarning on line %d in %s, %s\n\n", line_num, file_name, x);
}


void PrintHistElem(void *data)
{
   hist_data *h_data = (hist_data *)data;
   printf("nodes_tot = %d\n", h_data->nodes_tot);
   printf("success_node = %d\n", h_data->success_node);
   printf("success_its id = %d\n", h_data->success_its);
}

const char *ERR_STRING(int return_val)
{
   switch(return_val)
   {
   case GQ_SUCCESS:
      return STR_GQ_SUCCESS;
      break;
   case NULL_VAL:
      return STR_NULL_VAL;
      break;
   case ALLOC_FAIL:
      return STR_ALLOC_FAIL;
      break;
   case INV_INPUT:
      return STR_INV_INPUT;
      break;
   case INF_VAL:
      return STR_INF_VAL;
      break;
   case NAN_VAL:
      return STR_NAN_VAL;
      break;
   case QUAD_HUGE_ERR:
      return STR_QUAD_HUGE_ERR;
      break;
   case LAPACK_ERR:
      return STR_LAPACK_ERR;
      break;
   case PLASMA_ERR:
      return STR_PLASMA_ERR;
      break;
   case DIV_BY_ZERO:
      return STR_DIV_BY_ZERO;
      break;
   case NOT_CONVERGE:
      return STR_NOT_CONVERGE;
      break;
   case DIVERGE_ERR:
      return STR_DIVERGE_ERR;
      break;
   case LARGE_RES:
      return STR_LARGE_RES;
      break;
   case CONSTR_ERROR:
      return STR_CONSTR_ERROR;
   default:
      return "OTHER";
   }

}

