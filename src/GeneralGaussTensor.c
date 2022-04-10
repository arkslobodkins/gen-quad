/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "GeneralGaussTensor.h"
#include "Quadrature.h"
#include "Print.h"

#include <assert.h>
#include <stdio.h>

void NodesTensor2D(const quadrature *quad1, const quadrature *quad2, quadrature *quad_new)
{
   assert(quad1->num_nodes * quad2->num_nodes == quad_new->num_nodes);

   int n1 = quad1->num_nodes;
   int n2 = quad2->num_nodes;
   double *x1 = quad1->x;
   double *x2 = quad2->x;
   double *x_new = quad_new->x;

   for(int i = 0; i < n1; ++i) {
      for(int j = 0; j < n2; ++j) {
         x_new[2*(n2*i+j)]   = x1[i];
         x_new[2*(n2*i+j)+1] = x2[j];
      }
   }
}

void WeightsTensor2D(const quadrature *quad1, const quadrature *quad2, quadrature *quad_new)
{
   assert(quad1->num_nodes * quad2->num_nodes == quad_new->num_nodes);

   int n1 = quad1->num_nodes;
   int n2 = quad2->num_nodes;
   double *w1 = quad1->w;
   double *w2 = quad2->w;
   double *w_new = quad_new->w;

   for(int i = 0; i < n1; ++i)
      for(int j = 0; j < n2; ++j)
         w_new[i*n2+j] = w1[i] * w2[j];
}

void GeneralizedNodesTensor(const quadrature *quad_1D, quadrature *quad_gen)
{
   assert(quad_1D->dim == 1);
   assert(quad_gen->dim > 1);
   assert(quad_gen->num_nodes == POW_INT(quad_1D->num_nodes, quad_gen->dim));

   int dim = quad_gen->dim;
   int p = quad_1D->deg;
   int n = quad_1D->num_nodes;
   int n_sq = SQUARE(n);
   int id1, id2;
   double *x_gen = quad_gen->x;

   int dims_temp1[1] = {2};
   quadrature *quad_temp1 = quadrature_init_basic(n_sq, 2, dims_temp1, p, CUBE);
   double *x_t1   = quad_temp1->x;

   NodesTensor2D(quad_1D, quad_1D, quad_temp1);
   int count = POW_INT(n, dim-2);
   for(int i = 0; i < count; ++i) {
      id1 = dim-2 + i*dim*n_sq;
      id2 = dim-1 + i*dim*n_sq;
      for(int k = 0; k < n_sq; ++k) {
         x_gen[id1+k*dim] = x_t1[2*k];
         x_gen[id2+k*dim] = x_t1[2*k+1];
      }
   }
   quadrature_free(quad_temp1);


   for(int j = dim-3; j >= 0; --j)
   {
      count = POW_INT(n, j);
      for(int i = 0; i < count; ++i)
      {
         int dims_temp2[1] = {dim-j-1};
         quadrature *quad_temp2 = quadrature_init_basic(POW_INT(n, dim-j-1), dim-j-1, dims_temp2, p, CUBE);
         double *x_t2 = quad_temp2->x;

         int dims_temp3[1] = {dim-j};
         quadrature *quad_temp3 = quadrature_init_basic(POW_INT(n, dim-j), dim-j,  dims_temp3, p, CUBE);
         double *x_t3 = quad_temp3->x;

         for(int s = 0; s < POW_INT(n, dim-j-1); ++s)
            x_t2[s] = x_gen[s*dim+j+1];

         id1 = POW_INT(n, dim-j) * i * dim + j;
         NodesTensor2D(quad_1D, quad_temp2, quad_temp3);
         for(int k = 0; k < POW_INT(n, dim-j); ++k)
            x_gen[k*dim+id1] = x_t3[2*k];

         quadrature_free(quad_temp2);
         quadrature_free(quad_temp3);
      }
   }

}

void GeneralizedWeightsTensor(const quadrature *quad_1D, quadrature *quad_gen)
{
   assert(quad_1D->dim == 1);
   assert(quad_gen->dim > 1);
   assert(quad_gen->num_nodes == POW_INT(quad_1D->num_nodes, quad_gen->dim));

   int dim = quad_gen->dim;
   int n = quad_1D->num_nodes;
   int n_sq = SQUARE(n);
   int deg = quad_gen->deg;
   double *w_gen = quad_gen->w;

   int dims_temp1[3] = {2};
   quadrature *quad_temp1 = quadrature_init_basic(n_sq, 2, dims_temp1, deg, CUBE);
   double *w1 = quad_temp1->w;

   int counter = POW_INT(n, dim-2);
   for(int i = 0; i < counter; ++i) {
      WeightsTensor2D(quad_1D, quad_1D, quad_temp1);
      for(int k = 0; k < n_sq; ++k)
         w_gen[i*n_sq+k] = w1[k];
   }
   quadrature_free(quad_temp1);

   if(dim == 2) return;

   int dims_temp2[1] = {3};
   quadrature *quad_temp2 = quadrature_init_basic(POW_INT(n, 3), 3, dims_temp2, deg, CUBE);
   double *w2 = quad_temp2->w;

   int dims_temp3[1] = {3};
   quadrature *quad_temp3 = quadrature_init_basic(POW_INT(n, 3), 3, dims_temp3, deg, CUBE);
   double *w3 = quad_temp3->w;

   for(int j = dim-3; j >= 0; --j)
   {
      dims_temp2[0] = dim-j-1;
      quadrature_reinit_basic(POW_INT(n, dim-j-1), dim-j-1, dims_temp2, deg, quad_temp2);
      w2 = quad_temp2->w;

      dims_temp3[0] = dim-j;
      quadrature_reinit_basic(POW_INT(n, dim-j), dim-j, dims_temp3, deg, quad_temp3);
      w3 = quad_temp3->w;

      for(int i = 0; i < POW_INT(n, dim-j-1); ++i)
         w2[i] = w_gen[i];

      counter = POW_INT(n, j);
      for(int i = 0; i < counter; ++i) {
         WeightsTensor2D(quad_1D, quad_temp2, quad_temp3);
         for(int k = 0; k < POW_INT(n, dim-j-1) * n; ++k)
            w_gen[POW_INT(n, dim-j)*i+k] = w3[k];
      }
   }

   quadrature_free(quad_temp2);
   quadrature_free(quad_temp3);
}

void TestGeneralizedTensor()
{
   int deg = 1;
   int d1 = 1;
   int *dims1 = &d1;
   int nodes1D = 3;
   quadrature *quad1D = quadrature_init_basic(nodes1D, 1, dims1, deg, INTERVAL);
   quad1D->x[0] = 1.0; quad1D->x[1] = 2.0; quad1D->x[2] = 3.0;
   quad1D->w[0] = -1.0; quad1D->w[1] = -2.0; quad1D->w[2] = -3.0;

   int dimTensor = 3;
   int *dimsT = &dimTensor;
   int nodesGen = POW_INT(nodes1D, dimTensor);
   quadrature *quadTensor = quadrature_init_basic(nodesGen, dimTensor, dimsT, deg, CUBE);

   GeneralizedNodesTensor(quad1D, quadTensor);
   GeneralizedWeightsTensor(quad1D, quadTensor);
   FILE* file;
   file = fopen("../results/TensorTest.txt", "w");

   fprintf(file, "Testing Generalized Tensor product:");
   PrintNodesAndWeightsToFile(quad1D, "quad1D", file);
   PrintNodesAndWeightsToFile(quadTensor, "quadTensor", file);
   fclose(file);

   quadrature_free(quad1D);
   quadrature_free(quadTensor);
}
