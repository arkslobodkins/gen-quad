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
static void history_init(int k, elim_history *hist);
static void history_realloc(int k, elim_history *hist);
static void history_free(elim_history hist);


void ComputeInterval(int deg)
{
   assert(deg > 0);

   int n = ceil( (deg+1)/2.0 );
   int d_gauss[1] = {1};
   quadrature *q_gauss = quadrature_init_basic(n, 1, d_gauss, deg, INTERVAL);

   // generate Gaussian nodes and weights in 1-d on [0, 1]
   double p1 = 0.0, p2 = 0.0;
   Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);
   DumpCubatureRule((const_quadrature *)q_gauss);

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
// and analysis is printed to quadRule.txt and results.txt.
void ComputeCube(int deg, int dim)
{
   assert(deg > 0);
   assert(dim > 1);

   MAX_DIM = dim;
   int n_nodes_init = -1, n_nodes_next = -1;

   // generate Gaussian nodes and weights in 1-d on [0, 1]
   int n = ceil( (deg+1)/2.0 );
   int d_gauss[1] = {1};
   quadrature *q_gauss = quadrature_init_basic(n, 1, d_gauss, deg, INTERVAL);
   double p1 = 0.0, p2 = 0.0;
   Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);

   quadrature *q_init = NULL, *q_new = NULL;
   elim_history hist_cube = {0};
   // implement recursive scheme for computing the initial guess and run Node Elimination algorithm
   for(int d = 2; d <= dim; ++d)
   {
      if(d == 2) // set up initial guess for 2-dim CUBE
      {
         n_nodes_next = q_gauss->k * q_gauss->k;
         int d_init[1] = {d};
         q_init = quadrature_init_full(n_nodes_next, d, d_init, deg, CUBE);
         AddLineFirst( (const_quadrature *)q_gauss, (const_quadrature *)q_gauss, q_init );

         int d_new[1] = {d};
         q_new = quadrature_init_full(n_nodes_next, d, d_new, deg, CUBE);
         history_init(n_nodes_next, &hist_cube);
      }
      else if(d > 2) // set up initial guess for d-dim CUBE
      {
         n_nodes_next = q_new->k * q_gauss->k;
         int d_init[1] = {d};
         quadrature_realloc(n_nodes_next, d, d_init, deg, q_init);
         AddLineFirst( (const_quadrature *)q_gauss, (const_quadrature *)q_new, q_init );

         int d_new[1] = {d};
         quadrature_realloc(n_nodes_next, d, d_new, deg, q_new);
         history_realloc(n_nodes_next, &hist_cube);
      }

      if(d == dim) n_nodes_init = q_init->k;
      NodeElimination((const_quadrature *)q_init, q_new, &hist_cube);
   }

   // print results to results.txt and quadRule.txt
   double res = QuadTestIntegral( (const_quadrature *)q_new );
   Output(n_nodes_init, res, (const_quadrature *)q_new, hist_cube);
   DumpCubatureRule((const_quadrature *)q_new);

   quadrature_free(q_gauss);
   quadrature_free(q_init);
   quadrature_free(q_new);
   history_free(hist_cube);
}// ComputeCube




// ComputeSimplex
// Sets up the initial guess for the Node Elimination algorithm.
// Initial guess is computed via recursive scheme that reuses
// quadrature rules obtained in lower dimensions. More specifically,
// it computes the tensor product of (d-1)-dimensional quadrature
// rule over the unit simplex and interval and uses variant of Duffy
// Transformation to obtain the initial guess for quadrature over the
// d-dimensional unit simplex. When Node Elimination for the dimension
// assigned by a user is finalized, the final quadrature and analysis
// is printed to quadRule.txt and results.txt.
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
   int n_nodes_init = -1, n_nodes_next = -1;
   quadrature *q_init = NULL, *q_new = NULL;
   elim_history hist_simplex = {0};
   // implement recursive scheme for computing the initial guess and run Node Elimination algorithm
   for(int d = 2; d <= dim; ++d)
   {
      if(d == 2) // set up initial guess for 2-d SIMPLEX
      {
         double p1 = 0.0, p2 = 0.0;
         int d_gauss[1] = {1};

         int n1 = ceil( (deg+2)/2.0 );
         quadrature *q_gauss_1 = quadrature_init_basic(n1, 1,  d_gauss, 2*n1-1, INTERVAL);
         Jacobi(q_gauss_1->k, p1, p2, q_gauss_1->x, q_gauss_1->w);

         int n2 = ceil( (deg+1)/2.0 );
         quadrature *q_gauss_2 = quadrature_init_basic(n2, 1,  d_gauss, 2*n2-1, INTERVAL);
         Jacobi(q_gauss_2->k, p1, p2, q_gauss_2->x, q_gauss_2->w);

         n_nodes_next = q_gauss_1->k * q_gauss_2->k;
         int d_init[1] = {d};
         q_init = quadrature_init_full(n_nodes_next, d,  d_init, deg, SIMPLEX);
         AddLineSimplex( (const_quadrature *)q_gauss_1, (const_quadrature *)q_gauss_2, q_init);

         int d_new[1] = {d};
         q_new = quadrature_init_full(n_nodes_next, d,  d_new, deg, SIMPLEX);
         history_init(n_nodes_next, &hist_simplex);

         quadrature_free(q_gauss_1);
         quadrature_free(q_gauss_2);
      }
      else if(d > 2) // set up initial guess for d-dim SIMPLEX
      {
         double p1 = 0.0, p2 = 0.0;
         int n = ceil( (deg+d)/2.0 );
         int d_gauss[1] = {1};
         quadrature *q_gauss = quadrature_init_basic(n, 1,  d_gauss, 2*n-1, INTERVAL);
         Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);

         n_nodes_next = q_new->k * q_gauss->k;
         int d_init[1] = {d};
         quadrature_realloc(n_nodes_next, d, d_init, deg, q_init);
         AddLineSimplex( (const_quadrature *)q_gauss, (const_quadrature *)q_new, q_init );

         int d_new[1] = {d};
         quadrature_realloc(n_nodes_next, d, d_new, deg, q_new);
         quadrature_assign((const_quadrature *)q_init, q_new);
         history_realloc(n_nodes_next, &hist_simplex);

         quadrature_free(q_gauss);
      }

      if(d == dim) n_nodes_init = q_init->k;
      NodeElimination((const_quadrature *)q_init, q_new, &hist_simplex);
   }

   // print results to results.txt and quadRule.txt
   double res = QuadTestIntegral( (const_quadrature *)q_new );
   Output(n_nodes_init, res, (const_quadrature *)q_new, hist_simplex);
   DumpCubatureRule((const_quadrature *)q_new);

   quadrature_free(q_init);
   quadrature_free(q_new);
   history_free(hist_simplex);
}// ComputeSimplex


// ComputeCubeSimplex
// ets up the initial guess for the Node Elimination algorithm.
// Initial guess is computed via recursive scheme that reuses
// quadrature rules obtained in lower dimensions. At the first stage,
// it computes the tensor product of (d-1)-dimensional quadrature
// rule over the unit simplex and interval and uses variant of Duffy
// Transformation to obtain the initial guess for quadrature over the
// d-dimensional unit simplex. Later on, it computes tensor product of
// (d-1)-dimensional CUBESIMPLEX and interval to obtain initial guess for quadrature
// over d-dimensional CUBESIMPLEX. When Node Elimination for the dimensions
// assigned by a user is finalized, the final quadrature and analysis
// is printed to quadRule.txt and results.txt.
void ComputeCubeSimplex(int deg, int dim1, int dim2)
{
   assert(deg > 0);
   assert( (dim1 > 0) && (dim2 > 1) );

   int dim = dim1 + dim2;
   MAX_DIM = dim;
   int dim_cube = -1, dim_simplex = -1;
   int dim_simplex_max = dim2;
   int n_nodes_init = -1, n_nodes_next = -1;
   quadrature *q_init_s = NULL, *q_new_s = NULL, *q_init_cs = NULL, *q_new_cs = NULL;
   elim_history hist = {0};

   // implement recursive scheme for computing the initial guess and run Node Elimination algorithm
   for(int d = 2; d <= dim; ++d)
   {
      dim_simplex = MIN(d, dim_simplex_max); // increase SIMPLEX dimension, unless maximum dimension is reached
      dim_cube = d - dim_simplex;            // when dim_simplex is at maximum, increase dimension for CUBE


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
            n_nodes_next = q_gauss_1->k * q_gauss_2->k;
            q_init_s = quadrature_init_full(n_nodes_next, d, d_init, deg, SIMPLEX);
            AddLineSimplex( (const_quadrature *)q_gauss_1, (const_quadrature *)q_gauss_2, q_init_s);

            int d_new[1] = {d};
            q_new_s = quadrature_init_full(n_nodes_next, d, d_new, deg, SIMPLEX);
            history_init(n_nodes_next, &hist);

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

            n_nodes_next = q_new_s->k * q_gauss->k;
            int d_init[1] = {d};
            quadrature_realloc(n_nodes_next, d, d_init, deg, q_init_s);
            AddLineSimplex( (const_quadrature *)q_gauss, (const_quadrature *)q_new_s, q_init_s );

            int d_new[1] = {d};
            quadrature_realloc(n_nodes_next, d, d_new, deg, q_new_s);
            history_realloc(n_nodes_next, &hist);

            quadrature_free(q_gauss);
         }
         n_nodes_init = q_init_s->k;
         NodeElimination((const_quadrature *)q_init_s, q_new_s, &hist);
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
            n_nodes_next = q_new_s->k * q_gauss->k;
            q_init_cs = quadrature_init_full(n_nodes_next, d, d_init, deg, CUBESIMPLEX);
            AddLineFirst( (const_quadrature *)q_gauss, (const_quadrature *)q_new_s, q_init_cs );

            q_new_cs = quadrature_init_full(n_nodes_next, d, d_new, deg, CUBESIMPLEX);
            history_realloc(n_nodes_next, &hist);
         }
         else if(dim_cube > 1)
         {
            n_nodes_next = q_new_cs->k * q_gauss->k;
            quadrature_realloc(n_nodes_next, d, d_init, deg, q_init_cs);
            AddLineFirst( (const_quadrature *)q_gauss, (const_quadrature *)q_new_cs, q_init_cs );

            quadrature_realloc(n_nodes_next, d, d_new, deg, q_new_cs);
            history_realloc(n_nodes_next, &hist);
         }
         quadrature_free(q_gauss);

         if(d == dim) n_nodes_init = q_init_cs->k;
         NodeElimination((const_quadrature *)q_init_cs, q_new_cs, &hist);
      }


   }

   // print results to results.txt and quadRule.txt
   double res = QuadTestIntegral( (const_quadrature *)q_new_cs );
   Output(n_nodes_init, res, (const_quadrature *)q_new_cs, hist);
   DumpCubatureRule((const_quadrature *)q_new_cs);

   quadrature_free(q_init_s);
   quadrature_free(q_new_s);
   quadrature_free(q_init_cs);
   quadrature_free(q_new_cs);
   history_free(hist);
}// ComputeCubeSimplex


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
// is finalized, the final quadrature and analysis is printed to quadRule.txt and results.txt.
void ComputeSimplexSimplex(int deg, int dim1, int dim2)
{
   assert(deg > 0);
   assert( (dim1 > 1) && (dim2 > 1) );

   int dim = dim1 + dim2;
   MAX_DIM = dim;
   int dim_s1 = -1;
   // flip dimensions in a favourable order
   int dim_s1_max = dim1;
   int dim_s2_max = dim2;
   if(dim_s1_max < dim_s2_max)
   {
      int temp = dim_s2_max;
      dim_s2_max = dim_s1_max;
      dim_s1_max = temp;
   }

   int  n_nodes_init_ss = -1, n_nodes_next_s1 = -1, n_nodes_next_s2 = -1, n_nodes_next_ss = -1;
   quadrature *q_init_s1 = NULL, *q_init_ss = NULL, *q_new_s1 = NULL, *q_new_s2 = NULL, *q_new_ss = NULL;
   elim_history hist = {0};
   // implement recursive scheme for computing the initial guess and run Node Elimination algorithm
   for(int d = 2; d <= dim_s1_max; ++d)
   {
      dim_s1 = d;


      if(dim_s1 == 2) // set up initial guess for 2-d SIMPLEX
      {
         double p1 = 0.0, p2 = 0.0;
         int d_gauss[1] = {1};

         int n1 = ceil( (deg+2)/2.0 );
         quadrature *q_gauss_1 = quadrature_init_basic(n1, 1, d_gauss, 2*n1-1, INTERVAL);
         Jacobi(q_gauss_1->k, p1, p2, q_gauss_1->x, q_gauss_1->w);
         int n2 = ceil( (deg+1)/2.0 );
         quadrature *q_gauss_2 = quadrature_init_basic(n2, 1, d_gauss, deg, INTERVAL);
         Jacobi(q_gauss_2->k, p1, p2, q_gauss_2->x, q_gauss_2->w);

         n_nodes_next_s1 = q_gauss_1->k * q_gauss_2->k;
         int d_init[1] = {dim_s1};
         q_init_s1 = quadrature_init_full(n_nodes_next_s1, dim_s1, d_init, deg, SIMPLEX);
         AddLineSimplex( (const_quadrature *)q_gauss_1, (const_quadrature *)q_gauss_2, q_init_s1 );

         int d_new[1] = {dim_s1};
         q_new_s1 = quadrature_init_full(n_nodes_next_s1, dim_s1, d_new, deg, SIMPLEX);
         history_init(n_nodes_next_s1, &hist);

         quadrature_free(q_gauss_1);
         quadrature_free(q_gauss_2);
      }


      else if(dim_s1 > 2) // set up initial guess for d-dim SIMPLEX
      {
         int d_gauss[1] = {1};
         int n = ceil( (deg+dim_s1)/2.0 );
         quadrature *q_gauss = quadrature_init_basic(n, 1, d_gauss, deg, INTERVAL);
         double p1 = 0.0, p2 = 0.0;
         Jacobi(q_gauss->k, p1, p2, q_gauss->x, q_gauss->w);

         n_nodes_next_s1 = q_new_s1->k * q_gauss->k;
         int d_init[1] = {dim_s1};
         quadrature_realloc(n_nodes_next_s1, dim_s1, d_init, deg, q_init_s1);
         AddLineSimplex( (const_quadrature *)q_gauss, (const_quadrature *)q_new_s1, q_init_s1 );

         int d_s1[1] = {dim_s1};
         quadrature_realloc(n_nodes_next_s1, dim_s1, d_s1, deg, q_new_s1);
         history_realloc(n_nodes_next_s1, &hist);

         quadrature_free(q_gauss);
      }


      NodeElimination((const_quadrature *)q_init_s1, q_new_s1, &hist);
      n_nodes_next_s1 = q_new_s1->k;

      if(dim_s1 == dim_s2_max)
      {
         n_nodes_next_s2 = q_new_s1->k;
         int d_s2[1] = {d};
         q_new_s2 = quadrature_init_full(n_nodes_next_s2, d, d_s2, deg, SIMPLEX);
         quadrature_assign((const_quadrature *)q_new_s1, q_new_s2);
      }
   }

   // set up parameters for SIMPLEXSIMPLEX
   n_nodes_next_ss = n_nodes_next_s1 * n_nodes_next_s2;
   int d_ss[2] = {dim_s1_max, dim_s2_max};
   q_init_ss = quadrature_init_full(n_nodes_next_ss, dim, d_ss, deg, SIMPLEXSIMPLEX);
   q_new_ss = quadrature_init_full(n_nodes_next_ss, dim, d_ss, deg, SIMPLEXSIMPLEX);

   // compute initial guess for SIMPLEXSIMPLEX
   for(int i = 0; i < n_nodes_next_s1; ++i)
   {
      for(int j = 0; j < n_nodes_next_s2; ++j)
      {
         int id = i*n_nodes_next_s2*dim + j*dim;
         for(int d1 = 0; d1 < dim_s1_max; ++d1)
            q_init_ss->x[id+d1] = q_new_s1->x[i*dim_s1_max+d1];

         int count = 0;
         for(int d2 = dim_s1_max; d2 < dim; ++d2)
         {
            q_init_ss->x[id+d2] = q_new_s2->x[j*dim_s2_max+count];
            ++count;
         }
      }
   }
   WeightsTensor2D( (const_quadrature *)q_new_s1, (const_quadrature *)q_new_s2, q_init_ss );
   quadrature_assign((const_quadrature *)q_init_ss, q_new_ss);
   history_realloc(q_init_ss->k, &hist);

   // perform node elimination for SIMPLEXSIMPLEX
   n_nodes_init_ss = q_init_ss->k;
   NodeElimination((const_quadrature *)q_init_ss, q_new_ss, &hist);

   // print results to results.txt and quadRule.txt
   double res = QuadTestIntegral( (const_quadrature *)q_new_ss );
   Output(n_nodes_init_ss, res, (const_quadrature *)q_new_ss, hist);
   DumpCubatureRule((const_quadrature *)q_new_ss);

   quadrature_free(q_init_s1);
   quadrature_free(q_init_ss);
   quadrature_free(q_new_s1);
   quadrature_free(q_new_s2);
   quadrature_free(q_new_ss);
   history_free(hist);
}// ComputeSimplexSimplex


void ComputeCubeSimplexSimplex(int degree, int dim1, int dim2, int dim3)
{

}


static void history_init(int k, elim_history *hist)
{
   hist->tot_elims    = 0;
   hist->nodes_tot    = (int *)malloc(k*sizeof(int));
   hist->success_node = (int *)malloc(k*sizeof(int));
   hist->success_its  = (int *)malloc(k*sizeof(int));
}


static void history_realloc(int k, elim_history *hist)
{
   hist->tot_elims    = 0;
   hist->nodes_tot    = (int *)realloc(hist->nodes_tot, k*sizeof(int));
   hist->success_node = (int *)realloc(hist->success_node, k*sizeof(int));
   hist->success_its  = (int *)realloc(hist->success_its, k*sizeof(int));
}


static void history_free(elim_history hist)
{
   free(hist.nodes_tot);
   free(hist.success_its);
   free(hist.success_node);
}
