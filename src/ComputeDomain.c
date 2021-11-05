/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "ComputeDomain.h"

#include "AddDimension.h"
#include "BasisFunctions.h"
#include "BasisIndices.h"
#include "GeneralGaussTensor.h"
#include "Gauss_Lib/Jacobi.h"
#include "NodeElimination.h"
#include "Quadrature.h"
#include "Output.h"
#include "GENERAL_QUADRATURE.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

int MAX_DIM;

static history *hist_create();
static void hist_free(history *h_list);
static void hist_save_start(quadrature *q, history *hist);
static void hist_save_end(quadrature *q, history *hist);

static quadrature* ComputeSimplexNext(int deg, int dim_cur, int dim_next, const quadrature *q_cur, int hist_start, history **hist);


void ComputeInterval(int deg)
{
   assert(deg > 0);

   // generate Gaussian nodes and weights in 1-d on [0, 1]
   int n = ceil( (deg+1)/2.0 );
   int d_gauss[1] = {1};
   quadrature *q_gauss = quadrature_init_basic(n, 1, d_gauss, deg, INTERVAL);
   double p1 = 0.0, p2 = 0.0;
   Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);
   DumpCubatureRule(q_gauss);
   quadrature_free(q_gauss);
}


// ComputeCube
// Sets up the initial guess for the Node Elimination algorithm.
// Initial guess is computed via recursive scheme that reuses
// quadrature rules obtained in lower dimensions. More specifically,
// it computes the tensor product of (d-1)-dimensional quadrature
// rule over the unit cube and interval to obtain the initial guess for
// quadrature over the d-dimensional unit cube. When Node Elimination
// for the dimension assigned by a user is finalized, the final quadrature
// and history are printed to results directory.
void ComputeCube(int deg, int dim)
{
   assert(deg > 0);
   assert(dim > 1);

   MAX_DIM = dim;
   int hist_dims = dim-1;
   history **hist_cube = (history **)malloc(hist_dims*sizeof(history *));
   for(int i = 0; i < hist_dims; ++i)
      hist_cube[i] = hist_create();

   quadrature *q_init = NULL;
   quadrature *q_new  = NULL;

   // generate Gaussian nodes and weights in 1-d on [0, 1]
   int n = ceil( (deg+1)/2.0 );
   int d_gauss[1] = {1};
   quadrature *q_gauss = quadrature_init_basic(n, 1, d_gauss, deg, INTERVAL);
   double p1 = 0.0, p2 = 0.0;
   Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);

   // implement recursive scheme for computing the initial guess and run Node Elimination algorithm
   for(int d = 2; d <= dim; ++d)
   {
      if(d == 2) // set up initial guess for 2-dim CUBE
      {
         int n_nodes_next = q_gauss->k * q_gauss->k;
         int d_init[1] = {d};
         q_init = quadrature_init_full(n_nodes_next, d, d_init, deg, CUBE);
         AddLineFirst(q_gauss, q_gauss, q_init);

         int d_new[1] = {d};
         q_new = quadrature_init_full(n_nodes_next, d, d_new, deg, CUBE);
      }
      else if(d > 2) // set up initial guess for d-dim CUBE
      {
         int n_nodes_next = q_new->k * q_gauss->k;
         int d_init[1] = {d};
         quadrature_realloc(n_nodes_next, d, d_init, deg, q_init);
         AddLineFirst(q_gauss, q_new, q_init);

         int d_new[1] = {d};
         quadrature_realloc(n_nodes_next, d, d_new, deg, q_new);
      }

      hist_save_start(q_init, hist_cube[d-2]);
      NodeElimination(q_init, q_new, hist_cube[d-2]->list);
      hist_save_end(q_new, hist_cube[d-2]);
   }

   Output(q_new, dim-1, hist_cube);
   DumpCubatureRule(q_new);

   quadrature_free(q_gauss);
   quadrature_free(q_init);
   quadrature_free(q_new);

   for(int i = 0; i < hist_dims; ++i)
      hist_free(hist_cube[i]);
   free(hist_cube);
}// end ComputeCube



// ComputeSimplex
// Sets up the initial guess for the Node Elimination algorithm.
// Initial guess is computed via recursive scheme that reuses
// quadrature rules obtained in lower dimensions. More specifically,
// it computes the tensor product of (d-1)-dimensional quadrature
// rule over the unit simplex and interval and uses variant of Duffy
// Transformation to obtain the initial guess for quadrature over the
// d-dimensional unit simplex. When Node Elimination for the dimension
// assigned by a user is finalized, the final quadrature and history
// are printed to results directory.
void ComputeSimplex(int deg, int dim)
{
   assert(deg > 0);
   assert(dim > 1);

   {
      // test orthogonality of basis functions
      double x[dim];
      for(int d = 0; d < dim; ++d)
         x[d] = (double)rand() / (double) RAND_MAX;

      int num_funs = BasisSize(deg, dim);
      int8_t *basis_id = (int8_t *)malloc((num_funs*dim)*sizeof(int8_t));
      double *phi = (double *)malloc(num_funs*sizeof(double));
      BasisIndices(deg, dim, basis_id);

      int s_dims = dim;
      BasisSimplex(&s_dims, deg, basis_id, x, phi);
      double inf_norm = orthogonal_simplex_basis_test(&s_dims, deg, basis_id);
      printf("\nTesting orthogonality of basis functions. Maximum error of basis functions = %.16e\n\n", inf_norm);

      free(basis_id);
      free(phi);
   }

   MAX_DIM = dim;
   history **hist_simplex = (history **)malloc((dim-1)*sizeof(history *));
   for(int i = 0; i < dim-1; ++i)
   {
      hist_simplex[i] = hist_create();
      hist_simplex[i]->dim = i+2;
   }

   double p1 = 0.0, p2 = 0.0;
   int d_gauss[1] = {1};
   int n2 = ceil( (deg+1)/2.0 );
   quadrature *q_gauss = quadrature_init_basic(n2, 1,  d_gauss, 2*n2-1, INTERVAL);
   Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);
   quadrature *q_new = ComputeSimplexNext(deg, 1, dim, q_gauss, 0, hist_simplex);

   Output(q_new, dim-1, hist_simplex);
   DumpCubatureRule(q_new);

   quadrature_free(q_gauss);
   quadrature_free(q_new);
   for(int i = 0; i < dim-1; ++i)
      hist_free(hist_simplex[i]);
   free(hist_simplex);
}// end ComputeSimplex


// ComputeCubeSimplex
// Sets up the initial guess for the Node Elimination algorithm.
// Initial guess is computed via recursive scheme that reuses
// quadrature rules obtained in lower dimensions. At the first stage,
// it computes the tensor product of (d-1)-dimensional quadrature
// rule over the unit simplex and interval and uses variant of Duffy
// Transformation to obtain the initial guess for quadrature over the
// d-dimensional unit simplex. Later on, it computes tensor product of
// (d-1)-dimensional CUBESIMPLEX and interval to obtain initial guess for quadrature
// over d-dimensional CUBESIMPLEX. When Node Elimination for the dimensions
// assigned by a user is finalized, the final quadrature and history are printed to results directory.
void ComputeCubeSimplex(int deg, int dim1, int dim2)
{
   assert(deg > 0);
   assert((dim1 > 0) && (dim2 > 1));

   int dim = dim1 + dim2;
   MAX_DIM = dim;
   int dim_simplex_max = dim2;

   quadrature *q_init_s = NULL, *q_new_s = NULL, *q_init_cs = NULL, *q_new_cs = NULL;
   history **hist_cs = (history **)malloc((dim-1)*sizeof(history *));
   for(int i = 0; i < dim-1; ++i)
   {
      hist_cs[i]= hist_create();
      hist_cs[i]->dim = i+2;
   }

   // implement recursive scheme for computing the initial guess and run Node Elimination algorithm
   for(int d = 2; d <= dim; ++d)
   {
      int dim_simplex = MIN(d, dim_simplex_max); // increase SIMPLEX dimension, unless maximum dimension is reached
      int dim_cube = d - dim_simplex;            // when dim_simplex is at maximum, increase dimension for CUBE

      if(dim_cube == 0)
      {
         if(d == 2) // set up initial guess for 2-d SIMPLEX
         {
            int n1 = ceil( (deg+2)/2.0 );
            int d_gauss_1[1] = {1};
            quadrature *q_gauss_1 = quadrature_init_basic(n1, 1, d_gauss_1, 2*n1-1, INTERVAL);

            int n2 = ceil( (deg+1)/2.0 );
            int d_gauss_2[1] = {1};
            quadrature *q_gauss_2 = quadrature_init_basic(n2, 1, d_gauss_2, deg, INTERVAL);

            // generate Gaussian nodes and weights in 1-d on [0, 1]
            double p1 = 0.0, p2 = 0.0;
            Jacobi(q_gauss_1->k, p1, p2, q_gauss_1->x, q_gauss_1->w);
            Jacobi(q_gauss_2->k, p1, p2, q_gauss_2->x, q_gauss_2->w);

            int d_init[1] = {d};
            int n_nodes_next = q_gauss_1->k * q_gauss_2->k;
            q_init_s = quadrature_init_full(n_nodes_next, d, d_init, deg, SIMPLEX);
            AddLineSimplex(q_gauss_1, q_gauss_2, q_init_s);

            int d_new[1] = {d};
            q_new_s = quadrature_init_full(n_nodes_next, d, d_new, deg, SIMPLEX);

            quadrature_free(q_gauss_1);
            quadrature_free(q_gauss_2);
         }
         else if(d > 2) // set up initial guess for d-dim SIMPLEX
         {
            int n = ceil( (deg+d)/2.0 );
            int d_gauss[3] = {1};
            quadrature *q_gauss = quadrature_init_basic(n, 1, d_gauss, 2*n-1, INTERVAL);
            // generate Gaussian nodes and weights in 1-d on [0, 1]
            double p1 = 0.0, p2 = 0.0;
            Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);

            int n_nodes_next = q_new_s->k * q_gauss->k;
            int d_init[1] = {d};
            quadrature_realloc(n_nodes_next, d, d_init, deg, q_init_s);
            AddLineSimplex(q_gauss, q_new_s, q_init_s);

            int d_new[1] = {d};
            quadrature_realloc(n_nodes_next, d, d_new, deg, q_new_s);

            quadrature_free(q_gauss);
         }

         hist_save_start(q_init_s, hist_cs[d-2]);
         NodeElimination(q_init_s, q_new_s, hist_cs[d-2]->list);
         hist_save_end(q_new_s, hist_cs[d-2]);
      }
      else if(dim_cube > 0) // set up initial guess for CUBESIMPLEX
      {
         int n = ceil( (deg+1)/2.0 );
         int d_gauss[1] = {1};
         quadrature *q_gauss = quadrature_init_basic(n, 1, d_gauss, deg, INTERVAL);
         double p1 = 0.0, p2 = 0.0;
         Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);

         int d_init[2] = {dim_cube, dim_simplex};
         int d_new[2] = {dim_cube, dim_simplex};
         if(dim_cube == 1)
         {
            int n_nodes_next = q_new_s->k * q_gauss->k;
            q_init_cs = quadrature_init_full(n_nodes_next, d, d_init, deg, CUBESIMPLEX);
            AddLineFirst(q_gauss, q_new_s, q_init_cs);

            q_new_cs = quadrature_init_full(n_nodes_next, d, d_new, deg, CUBESIMPLEX);
         }
         else if(dim_cube > 1)
         {
            int n_nodes_next = q_new_cs->k * q_gauss->k;
            quadrature_realloc(n_nodes_next, d, d_init, deg, q_init_cs);
            AddLineFirst(q_gauss, q_new_cs, q_init_cs);

            quadrature_realloc(n_nodes_next, d, d_new, deg, q_new_cs);
         }
         quadrature_free(q_gauss);

         hist_save_start(q_init_cs, hist_cs[d-2]);
         NodeElimination(q_init_cs, q_new_cs, hist_cs[d-2]->list);
         hist_save_end(q_new_cs, hist_cs[d-2]);
      }

   }

   Output(q_new_cs, dim-1, hist_cs);
   DumpCubatureRule(q_new_cs);

   quadrature_free(q_init_s);
   quadrature_free(q_new_s);
   quadrature_free(q_init_cs);
   quadrature_free(q_new_cs);
   for(int i = 0; i < dim-1; ++i)
      hist_free(hist_cs[i]);
   free(hist_cs);
}// end ComputeCubeSimplex


// ComputeSimplexSimplex
// Sets up the initial guess for the Node Elimination algorithm.
// Initial guess is computed via recursive scheme that reuses
// quadrature rules obtained in lower dimensions. It computes
// the tensor product of (d-1)-dimensional quadrature rule over
// the unit simplex and interval and uses variant of Duffy Transformation
// to obtain the initial guess for quadrature over the d-dimensional unit simplex.
// When quadrature for dim1-dimensional simplex and dim2-dimensional simplices
// are computed by Node Elimination, the tensor product of dim1 and dim2 simplices
// is used as the final initial guess for SIMPLEXSIMPLEX. When Node Elimination for the dimensions assigned by a user
// is finalized, the final quadrature and history are printed to results directory.
void ComputeSimplexSimplex(int deg, int dim1, int dim2)
{
   assert(deg > 0);
   assert((dim1 > 1) && (dim2 > 1));

   int dim = dim1 + dim2;
   MAX_DIM = dim;
   // flip dimensions in a favourable order
   int dim_s1_max = dim1;
   int dim_s2_max = dim2;
   if(dim_s1_max < dim_s2_max)
   {
      int temp   = dim_s2_max;
      dim_s2_max = dim_s1_max;
      dim_s1_max = temp;
   }

   double p1 = 0.0, p2 = 0.0;
   int d_gauss[1] = {1};
   int n = ceil( (deg+1)/2.0 );
   quadrature *q_gauss = quadrature_init_basic(n, 1,  d_gauss, 2*n-1, INTERVAL);
   Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);

   int hist_dims = dim_s1_max;
   history **hist_ss = (history **)malloc((hist_dims)*sizeof(history *));
   for(int i = 0; i < hist_dims; ++i)
      hist_ss[i]= hist_create();

   quadrature *q_new_s1 = NULL;
   quadrature *q_new_s2 = ComputeSimplexNext(deg, 1, dim_s2_max, q_gauss, 0, hist_ss);
   if(dim_s2_max != dim_s1_max)
      q_new_s1 = ComputeSimplexNext(deg, dim_s2_max, dim_s1_max, q_new_s2, dim_s2_max-1, hist_ss);
   else
      q_new_s1 = quadrature_make_full_copy(q_new_s2);

   // set up parameters for SIMPLEXSIMPLEX
   int n_nodes_next_ss = q_new_s1->k * q_new_s2->k;
   int d_ss[2] = {dim_s1_max, dim_s2_max};
   quadrature *q_init_ss = quadrature_init_full(n_nodes_next_ss, dim, d_ss, deg, SIMPLEXSIMPLEX);
   quadrature *q_new_ss  = quadrature_init_full(n_nodes_next_ss, dim, d_ss, deg, SIMPLEXSIMPLEX);

   // compute initial guess for SIMPLEXSIMPLEX
   MixedTensor(q_new_s1, q_new_s2, q_init_ss);
   quadrature_assign(q_init_ss, q_new_ss);

   // perform node elimination for SIMPLEXSIMPLEX
   hist_save_start(q_init_ss, hist_ss[hist_dims-1]);
   NodeElimination(q_init_ss, q_new_ss, hist_ss[hist_dims-1]->list);
   hist_save_end(q_new_ss, hist_ss[hist_dims-1]);

   Output(q_new_ss, hist_dims, hist_ss);
   DumpCubatureRule(q_new_ss);

   quadrature_free(q_gauss);
   quadrature_free(q_init_ss);
   quadrature_free(q_new_s1);
   quadrature_free(q_new_s2);
   quadrature_free(q_new_ss);
   for(int i = 0; i < hist_dims; ++i)
      hist_free(hist_ss[i]);
   free(hist_ss);
}// end ComputeSimplexSimplex


void ComputeCubeSimplexSimplex(int degree, int dim1, int dim2, int dim3)
{

}


static quadrature *ComputeSimplexNext(int deg, int dim_cur, int dim_next, const quadrature *q_cur, int hist_start, history **hist)
{
   assert(deg > 0);
   assert(dim_cur > 0);
   assert(dim_next > 1);
   assert(dim_cur != dim_next);

   quadrature *q_init_next    = NULL;
   quadrature *q_new_next     = NULL;

   int count = 0;
   for(int d = dim_cur+1; d <= dim_next; ++d)
   {
      double p1 = 0.0, p2 = 0.0;
      int d_gauss[1] = {1};
      int n = ceil( (deg+d)/2.0 );
      quadrature *q_gauss = quadrature_init_basic(n, 1,  d_gauss, 2*n-1, INTERVAL);
      Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);

      if(d == dim_cur+1)
      {
         int n_nodes_next = q_gauss->k * q_cur->k;
         int d_init[1] = {d};
         q_init_next = quadrature_init_full(n_nodes_next, d, d_init, deg, SIMPLEX);
         AddLineSimplex(q_gauss, q_cur, q_init_next);

         int d_new[1] = {d};
         q_new_next = quadrature_init_full(n_nodes_next, d, d_new, deg, SIMPLEX);
      }
      else if(d > dim_cur+1)
      {
         int n_nodes_next = q_new_next->k * q_gauss->k;
         int d_init[1] = {d};
         quadrature_realloc(n_nodes_next, d, d_init, deg, q_init_next);
         AddLineSimplex(q_gauss, q_new_next, q_init_next );

         int d_new[1] = {d};
         quadrature_realloc(n_nodes_next, d, d_new, deg, q_new_next);
         quadrature_assign(q_init_next, q_new_next);
      }
      quadrature_free(q_gauss);

      hist_save_start(q_init_next, hist[hist_start+count]);
      NodeElimination(q_init_next, q_new_next, hist[hist_start+count]->list);
      hist_save_end(q_new_next, hist[hist_start+count]);
      ++count;
   }

   quadrature_free(q_init_next);

   return q_new_next;
}


static history *hist_create()
{
   history *h_list = (history *)malloc(sizeof(history));
   h_list->list = glist_create();
   return h_list;
}

static void hist_free(history *h_list)
{
   glist_free(h_list->list);
   free(h_list);
}

static void hist_save_start(quadrature *q, history *hist)
{
   hist->dim = q->dim;
   hist->degree = q->deg;
   hist->nodes_initial = q->k;
   hist->num_funcs = q->num_funcs;
   hist->D = q->D;
}

static void hist_save_end(quadrature *q, history *hist)
{
   hist->nodes_final = q->k;
   hist->res = QuadTestIntegral(q);
}


