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
void Output(const quadrature *q, int arr_size, history **hist_arr)
{
   char str[50];
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
      fprintf(FID, "dimension = %i\n", hist_arr[i]->dim);
      fprintf(FID, "degree of precision = %i\n", hist_arr[i]->degree);
      fprintf(FID, "DOMAIN_TYPE = %s\n", get_domain_string(hist_arr[i]->D));
      fprintf(FID, "number of basis functions = %i\n", hist_arr[i]->num_funcs);
      fprintf(FID, "initial number of nodes = %i\n", hist_arr[i]->nodes_initial);
      fprintf(FID, "final number of nodes = %i\n", hist_arr[i]->nodes_final);
      fprintf(FID, "total eliminations = %i\n", hist_arr[i]->list->size);
      fprintf(FID, "final residual = %.16e\n\n", hist_arr[i]->res);
      node *curr = hist_arr[i]->list->first;
      for(j = 0; j < hist_arr[i]->list->size; ++j)
      {
         hist_data *cur_d = (hist_data *)(curr->data);
         fprintf(FID, "total number of nodes[%i] = %i\n", j, cur_d->nodes_tot);
         fprintf(FID, "success_node[%i] = %i\n", j, cur_d->success_node);
         fprintf(FID, "Converged in %i iterations\n\n", cur_d->success_its);
         curr = curr->next;
      }
      fprintf(FID, "\n");
   }

   fclose(FID);
}


// Prints final quadrature to a file
void DumpCubatureRule(const quadrature *quad)
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
         return (char *)"LINE";
      case CUBE:
         return (char *)"CUBE";
      case SIMPLEX:
         return (char *)"SIMPLEX";
      case CUBESIMPLEX:
         return (char *)"CUBESIMPLEX";
      case SIMPLEXSIMPLEX:
         return  (char *)"SIMPLEXSIMPLEX";
      case CUBESIMPLEXSIMPLEX:
         return (char *)"CUBESIMPLEXSIMPLEX";
      default:
         return (char *)"0";
   }
}
