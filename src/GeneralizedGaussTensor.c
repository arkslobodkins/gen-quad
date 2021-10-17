/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GeneralizedGaussTensor.h"

#include "GaussTensor.h"
#include "Quadrature.h"
#include "TestIntegral.h"
#include "GENERAL_QUADRATURE.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>


void GeneralizedNodesTensor(const_quadrature *quad_1D, quadrature *quad_gen)
{
   assert(quad_1D->params->dim == 1);
   assert(quad_gen->params->dim > 1);
   int dim = quad_gen->params->dim; assert(quad_gen->k = POW(quad_1D->k, dim));


   int id1 = -1, id2 = -1;
   int p = quad_1D->params->deg;
   int n = quad_1D->k;
   int n_sq = SQUARE(n);

   int dims_temp1[1] = {2}; quadrature *quad_temp1 = quadrature_init(n_sq, 2, dims_temp1, p, CUBE);

   int count = pow(n, dim-2);
   NodesTensor( (const_quadrature *)quad_1D, (const_quadrature *)quad_1D, quad_temp1 );
   for(int i = 0; i < count; ++i)
   {
      id1 = dim-2 + i*dim*n_sq;
      id2 = dim-1 + i*dim*n_sq;
      for(int k = 0; k < n_sq; ++k)
      {
         quad_gen->x[id1+k*dim] = quad_temp1->x[2*k];
         quad_gen->x[id2+k*dim] = quad_temp1->x[2*k+1];
      }
   }
   quadrature_free(quad_temp1);


   for(int j = dim-3; j >= 0; --j)
   {
      count = POW(n, j);
      for(int i = 0; i < count; ++i)
      {
         int dims_temp2[1] = {dim-j-1}; quadrature *quad_temp2 = quadrature_init(POW(n, dim-j-1), dim-j-1, dims_temp2, p, CUBE);
         int dims_temp3[1] = {dim-j};   quadrature *quad_temp3 = quadrature_init(POW(n, dim-j), dim-j,  dims_temp3, p, CUBE);

         for(int s = 0; s < pow(n, dim-j-1); ++s)
            quad_temp2->x[s] = quad_gen->x[s*dim+j+1];

         NodesTensor((const_quadrature *)quad_1D, (const_quadrature *)quad_temp2, quad_temp3);
         id1 = i*(int)pow(n, dim-j)*dim + j;
         for(int k = 0; k < pow(n, dim-j); ++k)
            quad_gen->x[k*dim+id1] = quad_temp3->x[2*k];

         quadrature_free(quad_temp2);
         quadrature_free(quad_temp3);
      }
   }

}


void GeneralizedWeightsTensor(const_quadrature *quad_1D, quadrature *quad_gen)
{
   assert(quad_1D->params->dim == 1);
   assert(quad_gen->params->dim > 1);
   int dim = quad_gen->params->dim;
   assert(quad_gen->k = pow(quad_1D->k, dim));

   int n = quad_1D->k;
   int n_sq = SQUARE(n);
   int deg = quad_gen->params->deg;

   int dims_temp1[3] = {2};
   quadrature *quad_temp1 = quadrature_init(n_sq, 2, dims_temp1, deg, CUBE);
   int counter = pow(n, dim-2);
   for(int i = 0; i < counter; ++i)
   {
      WeightsTensor((const_quadrature *)quad_1D, (const_quadrature *)quad_1D, quad_temp1);
      for(int k = 0; k < n_sq; ++k)
         quad_gen->w[i*n_sq+k] = quad_temp1->w[k];
   }
   quadrature_free(quad_temp1);


   int dims_temp2[3] = {3}; quadrature *quad_temp2 = quadrature_init(pow(n, 3), 3, dims_temp2, deg, CUBE);
   int dims_temp3[3] = {3}; quadrature *quad_temp3 = quadrature_init(pow(n, 3), 3, dims_temp3, deg, CUBE);

   for(int j = dim-3; j >= 0; --j)
   {
      int dims_temp2[1] = {dim-j-1}; quadrature_realloc(pow(n, dim-j-1), dim-j-1, dims_temp2, deg, quad_temp2);
      int dims_temp3[1] = {dim-j};   quadrature_realloc(pow(n, dim-j), dim-j, dims_temp3, deg, quad_temp3);

      for(int i = 0; i < pow(n, dim-j-1); ++i)
         quad_temp2->w[i] = quad_gen->w[i];

      counter = pow(n, j);
      for(int i = 0; i < counter; ++i)
      {
         WeightsTensor((const_quadrature *)quad_1D, (const_quadrature *)quad_temp2, quad_temp3);
         for(int k = 0; k < pow(n, dim-j-1) * n; ++k)
            quad_gen->w[i*(int)pow(n, dim-j)+k] = quad_temp3->w[k];

      }
   }
   quadrature_free(quad_temp2);
   quadrature_free(quad_temp3);

}
