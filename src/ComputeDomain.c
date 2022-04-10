/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "ComputeDomain.h"

#include "AddDimension.h"
#include "GeneralGaussTensor.h"
#include "NodeElimination.h"
#include "Quadrature.h"
#include "Output.h"
#include "Print.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

int MAX_DIM;

static history *hist_create(int n);
static void hist_free(history *h_list);
static void hist_save_start(quadrature *q, history *hist);
static void hist_save_end(quadrature *q, history *hist);

static quadrature* ComputeCubeNext(int deg, int dim_cur, int dim_next,
                                   const quadrature *q_cur, int hist_start, history **hist);
static quadrature* ComputeSimplexNext(int deg, int dim_cur, int dim_next,
                                      const quadrature *q_cur, int hist_start, history **hist);


void ComputeInterval(int deg)
{
   assert(deg > 0);

   quadrature *q_gauss = quadrature_gauss_legendre(deg);
   QuadratureToFile(q_gauss);
   quadrature_free(q_gauss);
}


void ComputeCubeFull(int deg, int dim)
{
   assert(deg > 0);
   assert(dim > 1);
   MAX_DIM = dim;

   quadrature *q_init = quadrature_full_cube_tensor(deg, dim);
   quadrature *q_new  = quadrature_make_full_copy(q_init);

   history *hist_cube = (history *)malloc(sizeof(history));
   hist_cube = hist_create(q_new->num_nodes);

   NodeElimination(q_init, q_new, hist_cube);

   quadrature_free(q_init);
   quadrature_free(q_new);
   hist_free(hist_cube);
}


void ComputeSimplexFull(int deg, int dim)
{
   assert(deg > 0);
   assert(dim > 1);
   MAX_DIM = dim;

   quadrature *q_init = quadrature_full_simplex_tensor(deg, dim);
   quadrature *q_new  = quadrature_make_full_copy(q_init);

   history *hist_cube = (history *)malloc(sizeof(history));
   hist_cube = hist_create(q_new->num_nodes);

   NodeElimination(q_init, q_new, hist_cube);

   quadrature_free(q_init);
   quadrature_free(q_new);
   hist_free(hist_cube);
}


void ComputeCube(int deg, int dim)
{
   assert(deg > 0);
   assert(dim > 1);
   MAX_DIM = dim;

#ifdef QUAD_DEBUG_ON
   double basisRes = orthogonal_cube_basis_test(deg, dim);
   if(basisRes > POW_DOUBLE(10.0, -13))
   {
      PRINT_ERR("Residual of orthogonal cube basis functions is too large", __LINE__, __FILE__);
      return;
   }
#endif


   int hist_dims = dim-1;
   history **hist_cube = (history **)malloc(hist_dims*sizeof(history *));

   quadrature *q_gauss = quadrature_gauss_legendre(deg);
   quadrature *q_new = ComputeCubeNext(deg, 1, dim, q_gauss, 0, hist_cube);

   double res = QuadTestIntegral(q_new, monomial);
   printf("Monomial basis residual = %.16e\n", res);
   res = QuadTestIntegralExp(q_new);
   printf("Exponential residual = %.16e\n\n", res);

   HistoryToFile(q_new, dim-1, hist_cube);
   QuadratureToFile(q_new);
//   BoundaryCubeStats(q_new);

   quadrature_free(q_gauss);
   quadrature_free(q_new);

   for(int i = 0; i < hist_dims; ++i)
      hist_free(hist_cube[i]);
   free(hist_cube);
}// end ComputeCube


void ComputeSimplex(int deg, int dim)
{
   assert(deg > 0);
   assert(dim > 1);
   MAX_DIM = dim;

#ifdef QUAD_DEBUG_ON
   double basisRes = orthogonal_simplex_basis_test(deg, dim);
   if(basisRes > POW_DOUBLE(10.0, -13))
   {
      PRINT_ERR("Residual of orthogonal simplex basis functions is too large", __LINE__, __FILE__);
      return;
   }
#endif

   int hist_dims = dim-1;
   history **hist_simplex = (history **)malloc((hist_dims)*sizeof(history *));

   quadrature *q_gauss = quadrature_gauss_legendre(deg);
   quadrature *q_new = ComputeSimplexNext(deg, 1, dim, q_gauss, 0, hist_simplex);

   double res = QuadTestIntegral(q_new, monomial);
   printf("Monomial basis residual = %.16e\n", res);
   res = QuadTestIntegralExp(q_new);
   printf("Exponential residual = %.16e\n\n", res);

   HistoryToFile(q_new, dim-1, hist_simplex);
   QuadratureToFile(q_new);

   quadrature_free(q_gauss);
   quadrature_free(q_new);
   for(int i = 0; i < hist_dims; ++i)
      hist_free(hist_simplex[i]);
   free(hist_simplex);
}// end ComputeSimplex


void ComputeCubeSimplex(int deg, int dim1, int dim2)
{
   assert(deg > 0);
   assert((dim1 > 0) && (dim2 > 1));

   int dim = dim1 + dim2;
   MAX_DIM = dim;
   int dim_simplex_max = dim2;

   quadrature *q_new_s = NULL, *q_init_cs = NULL, *q_new_cs = NULL;
   history **hist_cs = (history **)malloc((dim-1)*sizeof(history *));

   quadrature *q_gauss = quadrature_gauss_legendre(deg);

   // implement recursive scheme for computing the initial guess and run Node Elimination algorithm
   q_new_s = ComputeSimplexNext(deg, 1, dim_simplex_max, q_gauss, 0, hist_cs);
   for(int d = 2; d <= dim; ++d)
   {
      int dim_simplex = MIN(d, dim_simplex_max); // increase SIMPLEX dimension, unless maximum dimension is reached
      int dim_cube = d - dim_simplex;            // when dim_simplex is at maximum, increase dimension for CUBE

      if(dim_cube > 0) // set up initial guess for CUBESIMPLEX
      {

         int d_init[2] = {dim_cube, dim_simplex};
         int d_new[2] = {dim_cube, dim_simplex};
         if(dim_cube == 1) {
            int n_nodes_next = q_new_s->num_nodes * q_gauss->num_nodes;
            q_init_cs = quadrature_init_full(n_nodes_next, d, d_init, deg, CUBESIMPLEX);
            AddLineFirst(q_gauss, q_new_s, q_init_cs);

            q_new_cs = quadrature_init_full(n_nodes_next, d, d_new, deg, CUBESIMPLEX);
         }
         else if(dim_cube > 1) {
            int n_nodes_next = q_new_cs->num_nodes * q_gauss->num_nodes;
            quadrature_free(q_init_cs);
            q_init_cs = quadrature_init_full(n_nodes_next, d, d_init, deg, CUBESIMPLEX);
            AddLineFirst(q_gauss, q_new_cs, q_init_cs);

            quadrature_free(q_new_cs);
            q_new_cs = quadrature_init_full(n_nodes_next, d, d_new, deg, CUBESIMPLEX);
         }

         hist_cs[d-2] = hist_create(q_new_cs->num_nodes);
         hist_save_start(q_init_cs, hist_cs[d-2]);
         NodeElimination(q_init_cs, q_new_cs, hist_cs[d-2]);
         hist_save_end(q_new_cs, hist_cs[d-2]);
      }
   }

   double res = QuadTestIntegral(q_new_cs, monomial);
   printf("Monomial basis residual = %.16e\n", res);
   res = QuadTestIntegralExp(q_new_cs);
   printf("Exponential residual = %.16e\n\n", res);
   HistoryToFile(q_new_cs, dim-1, hist_cs);
   QuadratureToFile(q_new_cs);

   quadrature_free(q_gauss);
   quadrature_free(q_new_s);
   quadrature_free(q_init_cs);
   quadrature_free(q_new_cs);
   for(int i = 0; i < dim-1; ++i)
      hist_free(hist_cs[i]);
   free(hist_cs);
}// end ComputeCubeSimplex


void ComputeCubeSimplexTensor(int deg, int dim1, int dim2)
{
   assert(deg > 0);
   assert((dim1 > 0) && (dim2 > 1));

   int dim = dim1 + dim2;
   MAX_DIM = dim;
   int dim_cube = dim1;
   int dim_simplex = dim2;

   quadrature *q_gauss = quadrature_gauss_legendre(deg);

   int hist_dims = dim-1;
   history **hist_cs = (history **)malloc((hist_dims)*sizeof(history *));

   quadrature *q_new_cube = NULL;
   if(dim_cube > 1) q_new_cube = ComputeCubeNext(deg, 1, dim_cube, q_gauss, 0, hist_cs);
   else             q_new_cube = quadrature_gauss_legendre(deg);
   quadrature *q_new_simplex = ComputeSimplexNext(deg, 1, dim_simplex, q_gauss, dim_cube-1, hist_cs);

   // set up parameters for CUBESIMPLEX
   int n_nodes_next_cs = q_new_cube->num_nodes * q_new_simplex->num_nodes;
   int d_cs[2] = {dim_cube, dim_simplex};
   quadrature *q_init_cs = quadrature_init_full(n_nodes_next_cs, dim, d_cs, deg, CUBESIMPLEX);
   quadrature *q_new_cs  = quadrature_init_full(n_nodes_next_cs, dim, d_cs, deg, CUBESIMPLEX);

   // compute initial guess for CUBESIMPLEX
   MixedTensor(q_new_cube, q_new_simplex, q_init_cs);
   quadrature_assign(q_init_cs, q_new_cs);

   // perform node elimination for CUBESIMPLEX
   hist_cs[hist_dims-1] = hist_create(q_init_cs->num_nodes);
   hist_save_start(q_init_cs, hist_cs[hist_dims-1]);
   NodeElimination(q_init_cs, q_new_cs, hist_cs[hist_dims-1]);
   hist_save_end(q_new_cs, hist_cs[hist_dims-1]);

   double res = QuadTestIntegral(q_new_cs, monomial);
   printf("Monomial basis residual = %.16e\n", res);
   res = QuadTestIntegralExp(q_new_cs);
   printf("Exponential residual = %.16e\n\n", res);
   HistoryToFile(q_new_cs, hist_dims, hist_cs);
   QuadratureToFile(q_new_cs);

   quadrature_free(q_gauss);
   quadrature_free(q_init_cs);
   quadrature_free(q_new_cube);
   quadrature_free(q_new_simplex);
   quadrature_free(q_new_cs);
   for(int i = 0; i < hist_dims; ++i)
      hist_free(hist_cs[i]);
   free(hist_cs);
}// end ComputeCubeSimplexTensor


void ComputeSimplexSimplex(int deg, int dim1, int dim2)
{
   assert(deg > 0);
   assert((dim1 > 1) && (dim2 > 1));

   int dim = dim1 + dim2;
   MAX_DIM = dim;
   // flip dimensions in a favourable order
   int dim_s1_max = dim1;
   int dim_s2_max = dim2;
   if(dim_s1_max < dim_s2_max) {
      int temp   = dim_s2_max;
      dim_s2_max = dim_s1_max;
      dim_s1_max = temp;
   }

   quadrature *q_gauss = quadrature_gauss_legendre(deg);

   int hist_dims = dim_s1_max;
   history **hist_cs = (history **)malloc((hist_dims)*sizeof(history *));

   quadrature *q_new_s1 = NULL;
   quadrature *q_new_s2 = ComputeSimplexNext(deg, 1, dim_s2_max, q_gauss, 0, hist_cs);
   if(dim_s2_max != dim_s1_max)
      q_new_s1 = ComputeSimplexNext(deg, dim_s2_max, dim_s1_max, q_new_s2, dim_s2_max-1, hist_cs);
   else
      q_new_s1 = quadrature_make_full_copy(q_new_s2);

   // set up parameters for SIMPLEXSIMPLEX
   int n_nodes_next_ss = q_new_s1->num_nodes * q_new_s2->num_nodes;
   int d_ss[2] = {dim_s1_max, dim_s2_max};
   quadrature *q_init_ss = quadrature_init_full(n_nodes_next_ss, dim, d_ss, deg, SIMPLEXSIMPLEX);
   quadrature *q_new_ss  = quadrature_init_full(n_nodes_next_ss, dim, d_ss, deg, SIMPLEXSIMPLEX);

   // compute initial guess for SIMPLEXSIMPLEX
   MixedTensor(q_new_s1, q_new_s2, q_init_ss);
   quadrature_assign(q_init_ss, q_new_ss);

   // perform node elimination for SIMPLEXSIMPLEX
   hist_cs[hist_dims-1]= hist_create(q_init_ss->num_nodes);
   hist_save_start(q_init_ss, hist_cs[hist_dims-1]);
   NodeElimination(q_init_ss, q_new_ss, hist_cs[hist_dims-1]);
   hist_save_end(q_new_ss, hist_cs[hist_dims-1]);

   double res = QuadTestIntegral(q_new_ss, monomial);
   printf("Monomial basis residual = %.16e\n", res);
   res = QuadTestIntegralExp(q_new_ss);
   printf("Exponential residual = %.16e\n\n", res);
   HistoryToFile(q_new_ss, hist_dims, hist_cs);
   QuadratureToFile(q_new_ss);

   quadrature_free(q_gauss);
   quadrature_free(q_init_ss);
   quadrature_free(q_new_s1);
   quadrature_free(q_new_s2);
   quadrature_free(q_new_ss);
   for(int i = 0; i < hist_dims; ++i)
      hist_free(hist_cs[i]);
   free(hist_cs);
}// end ComputeSimplexSimplexTensor


static quadrature *ComputeCubeNext(int deg, int dim_cur, int dim_next, const quadrature *q_cur, int hist_start, history **hist)
{
   assert(deg > 0);
   assert(dim_cur > 0);
   assert(dim_next > 1);
   assert(dim_cur != dim_next);

   quadrature *q_init_next = NULL;
   quadrature *q_new_next = NULL;

   int count = 0;
   for(int d = dim_cur+1; d <= dim_next; ++d)
   {
      quadrature *q_gauss = quadrature_gauss_legendre(deg);

      if(d == dim_cur+1) {
         int n_nodes_next = q_gauss->num_nodes * q_cur->num_nodes;
         int d_init[1] = {d};
         q_init_next = quadrature_init_full(n_nodes_next, d, d_init, deg, CUBE);
         AddLineFirst(q_gauss, q_cur, q_init_next);

         int d_new[1] = {d};
         q_new_next = quadrature_init_full(n_nodes_next, d, d_new, deg, CUBE);
      }
      else if(d > dim_cur+1) {
         int n_nodes_next = q_new_next->num_nodes * q_gauss->num_nodes;
         int d_init[1] = {d};
         quadrature_free(q_init_next);
         q_init_next = quadrature_init_full(n_nodes_next, d, d_init, deg, CUBE);
         AddLineFirst(q_gauss, q_new_next, q_init_next);

         int d_new[1] = {d};
         quadrature_free(q_new_next);
         q_new_next = quadrature_init_full(n_nodes_next, d, d_new, deg, CUBE);
         quadrature_assign(q_init_next, q_new_next);
      }
      quadrature_free(q_gauss);

      hist[hist_start+count] = hist_create(q_init_next->num_nodes);
      hist_save_start(q_init_next, hist[hist_start+count]);
      NodeElimination(q_init_next, q_new_next, hist[hist_start+count]);
      hist_save_end(q_new_next, hist[hist_start+count]);
      ++count;
   }

   quadrature_free(q_init_next);

   return q_new_next;
}


static quadrature *ComputeSimplexNext(int deg, int dim_cur, int dim_next, const quadrature *q_cur, int hist_start, history **hist)
{
   assert(deg > 0);
   assert(dim_cur > 0);
   assert(dim_next > 1);
   assert(dim_cur != dim_next);

   quadrature *q_init_next = NULL;
   quadrature *q_new_next  = NULL;

   int count = 0;
   for(int d = dim_cur+1; d <= dim_next; ++d)
   {
      quadrature *q_gauss = quadrature_gauss_jacobi(deg, d-1.0, 0.0);

      if(d == dim_cur+1) {
         int n_nodes_next = q_gauss->num_nodes * q_cur->num_nodes;
         int d_init[1] = {d};
         q_init_next = quadrature_init_full(n_nodes_next, d, d_init, deg, SIMPLEX);
         AddLineSimplex(q_gauss, q_cur, q_init_next);

         int d_new[1] = {d};
         q_new_next = quadrature_init_full(n_nodes_next, d, d_new, deg, SIMPLEX);
      }
      else if(d > dim_cur+1) {
         int n_nodes_next = q_new_next->num_nodes * q_gauss->num_nodes;
         int d_init[1] = {d};
         quadrature_free(q_init_next);
         q_init_next = quadrature_init_full(n_nodes_next, d, d_init, deg, SIMPLEX);
         AddLineSimplex(q_gauss, q_new_next, q_init_next);

         int d_new[1] = {d};
         quadrature_free(q_new_next);
         q_new_next = quadrature_init_full(n_nodes_next, d, d_new, deg, SIMPLEX);
         quadrature_assign(q_init_next, q_new_next);
      }
      quadrature_free(q_gauss);

      hist[hist_start+count] = hist_create(q_init_next->num_nodes);
      hist_save_start(q_init_next, hist[hist_start+count]);
      NodeElimination(q_init_next, q_new_next, hist[hist_start+count]);
      hist_save_end(q_new_next, hist[hist_start+count]);
      ++count;
   }

   quadrature_free(q_init_next);

   return q_new_next;
}


static history *hist_create(int n)
{
   history *h = (history *)malloc(sizeof(history));
   h->hist_array = (hist_data *)malloc(n*sizeof(hist_data));
   h->total_elims = 0;

   for(int i = 0; i < n; ++i)
      h->hist_array[i].num_solutions = 0;

   return h;
}

static void hist_free(history *h)
{
   free(h->hist_array);
   free(h);
}

static void hist_save_start(quadrature *q, history *hist)
{
   hist->dim           = q->dim;
   hist->degree        = q->deg;
   hist->nodes_initial = q->num_nodes;
   hist->nodes_optimal = ceil( 1.0 * q->basis->numFuncs / (q->dim+1) );
   hist->num_funcs     = q->basis->numFuncs;
   hist->D             = q->D;
   hist->efficiency    = 0.0;
}

static void hist_save_end(quadrature *q, history *hist)
{
   hist->nodes_final = q->num_nodes;
   hist->res = QuadTestIntegral(q, orthogonal);
}



