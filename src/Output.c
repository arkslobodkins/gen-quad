/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */
// to do: ensure results directory exists

#include "Output.h"
#include "GENERAL_QUADRATURE.h"
#include <stdio.h>

static char *get_domain_string(DOMAIN_TYPE D);


// Prints history and quadrature parameters to a file.
void Output(int num_nodes_initial, double res, const_quadrature *q, const elim_history hist)
{
   char str[50];
   int i;
   int eliminations = hist.tot_elims;
   int deg = q->deg;
   int dim = q->dim;
   int m = q->num_funcs;
   char *shape = get_domain_string(q->D);
   sprintf(str, "../results/history_%s_dim%i_deg%i.txt", shape, dim, deg);
   FILE* FID;
   FID = fopen(str, "w");
   fprintf(FID, "degree of precision = %i\n", deg);
   fprintf(FID, "number of basis functions = %i\n", m);
   fprintf(FID, "dimension = %i\n", dim);
   fprintf(FID, "initial number of nodes = %i\n", num_nodes_initial);
   fprintf(FID, "final number of nodes = %i\n", q->k);
   fprintf(FID, "final residual = %0.16f\n\n", res);

   for(i = 0; i < eliminations; ++i)
   {
      fprintf(FID, "total number of nodes[%i] = %i\n", i, hist.nodes_tot[i]);
      fprintf(FID, "success_node[%i] = %i\n", i, hist.success_node[i]);
      fprintf(FID, "Converged in %i iterations\n\n", hist.success_its[i]);
   }

   fclose(FID);
}


// Prints final quadrature to a file
void DumpCubatureRule(const_quadrature *quad)
{
   int i = -1, j = -1;
   char str[50] = "0";
   int dim = quad->dim;
   int deg = quad->deg;
   const char *shape = get_domain_string(quad->D);

   sprintf(str, "../results/quadrature_%s_dim%i_deg%i.txt", shape, dim, deg);
   FILE *FID = fopen(str, "w");
   fprintf(FID, "%i %i\n", dim, quad->k);

   for(i = 0; i < quad->k; ++i)
   {
      for(j = 0; j < dim; ++j)
         fprintf(FID, "%.16e  ", quad->x[dim*i+j]);

      fprintf(FID, "%.16e\n", quad->w[i]);
   }

   fclose(FID);
}


static char *get_domain_string(DOMAIN_TYPE D)
{
   switch(D)
   {
      case INTERVAL:
         return (char *)"line";
      case CUBE:
         return (char *)"cube";
      case SIMPLEX:
         return (char *)"simplex";
      case CUBESIMPLEX:
         return (char *)"cubesimplex";
      case SIMPLEXSIMPLEX:
         return  (char *)"simplexsimplex";
      case CUBESIMPLEXSIMPLEX:
         return (char *)"cubesimplexsimplex";
      default:
         return (char *)"0";
   }
}
