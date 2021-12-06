/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */
// to do: ensure results directory exists

#include "Output.h"
#include <stdio.h>


// Prints history and quadrature parameters to a file.
void HistoryToFile(const quadrature *q, int arr_size, history **hist)
{
   char str[80];
   int i,j;
   int deg = q->deg;
   int dim = q->dim;
   char *shape = get_domain_string(q->D);
   sprintf(str, "../results/history_%s_dim%i_deg%i.txt", shape, dim, deg);
   FILE* FID;
   FID = fopen(str, "w");

   for(i = 0; i < arr_size; ++i)
   {
      fprintf(FID, "*********************************************************\n\n");
      fprintf(FID, "dimension = %i\n", hist[i]->dim);
      fprintf(FID, "degree of precision = %i\n", hist[i]->degree);
      fprintf(FID, "DOMAIN_TYPE = %s\n", get_domain_string(hist[i]->D));
      fprintf(FID, "number of basis functions = %i\n", hist[i]->num_funcs);
      fprintf(FID, "initial number of nodes = %i\n", hist[i]->nodes_initial);
      fprintf(FID, "final number of nodes = %i\n", hist[i]->nodes_final);
      fprintf(FID, "total eliminations = %i\n", hist[i]->total_elims);
      fprintf(FID, "final residual = %.16e\n\n", hist[i]->res);
      for(j = 0; j < hist[i]->total_elims; ++j)
      {
         fprintf(FID, "total number of nodes[%i] = %i\n", j, hist[i]->hist_array[j].nodes_tot);
         fprintf(FID, "success_node[%i] = %i\n", j, hist[i]->hist_array[j].success_node);
         fprintf(FID, "Converged in %i iterations\n\n", hist[i]->hist_array[j].success_its);
      }
      fprintf(FID, "\n");
   }

   fclose(FID);
}

// Prints final quadrature to a file
void QuadratureToFile(const quadrature *quad)
{
   int i = -1, j = -1;
   char str[50] = "0";
   int dim = quad->dim;
   int deg = quad->deg;
   const char *shape = get_domain_string(quad->D);

   sprintf(str, "../results/quadrature_%s_dim%i_deg%i.txt", shape, dim, deg);
   FILE *FID = fopen(str, "w");
   fprintf(FID, "%i %i\n", dim, quad->num_nodes);

   for(i = 0; i < quad->num_nodes; ++i)
   {
      for(j = 0; j < dim; ++j)
         fprintf(FID, "%.16e  ", quad->x[dim*i+j]);

      fprintf(FID, "%.16e\n", quad->w[i]);
   }

   fclose(FID);
}

