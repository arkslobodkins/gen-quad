/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "InDomain.h"

#include "Constraints.h"
#include "Quadrature.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>


bool InCubeElem(const_quadrature *q, int elem)
{
   int dim = q->params->dim;
   constraints *cons = constraints_cube_init(q->params->dims);
   get_constraints_cube(cons);
   Matrix M = cons->M;
   Vector b = cons->b;

   double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
   const double *X_ixdim = &q->x[elem*dim];

   for(int r = 0; r < M.rows; ++r)
      for(int c = 0; c < M.cols; ++c)
         lhs[r] += M.id[r][c] * X_ixdim[c];

   for(int r = 0; r < M.rows; ++r)
   {
      if(lhs[r] > b.id[r])
      {
         constraints_cube_free(cons);
         return false;
      }
   }

   constraints_cube_free(cons);
   return true;
}

//bool InCubeElem_d(double *x, int dim)
//{
//   int dim = q->params->dim;
//   constraints cons = constraints_cube_init(q->params->dims);
//   get_constraints_cube(cons);
//   Matrix M = cons.M;
//   Vector b = cons.b;
//
//   double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
//   const double *X_ixdim = &q->x[elem*dim];
//
//   for(int r = 0; r < M.rows; ++r)
//      for(int c = 0; c < M.cols; ++c)
//         lhs[r] += M.id[r][c] * X_ixdim[c];
//
//   for(int r = 0; r < M.rows; ++r)
//   {
//      if(lhs[r] > b.id[r])
//      {
//         constraints_cube_free(cons);
//         return false;
//      }
//   }
//
//   constraints_cube_free(cons);
//   return true;
//}


bool InSimplexElem(const_quadrature *q, int elem)
{

   return true;
}
bool InCubeSimplexElem(const_quadrature *q, int elem)
{

   return true;
}
bool InSimplexSimplexElem(const_quadrature *q, int elem)
{

   return true;
}
bool InCubeSimplexSimplexElem(const_quadrature *q, int elem)
{

   return true;
}

/* InCube
 * Tests whether nodes are inside unit cube.
 */
bool InCube(const_quadrature *q)
{
   int dim = q->params->dim;
   constraints *cons = constraints_cube_init(q->params->dims);
   get_constraints_cube(cons);
   Matrix M = cons->M;
   Vector b = cons->b;

   for(int i = 0; i < q->k; ++i)
   {
      double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
      const double *X_ixdim = &q->x[i*dim];

      for(int r = 0; r < M.rows; ++r)
         for(int c = 0; c < M.cols; ++c)
            lhs[r] += M.id[r][c] * X_ixdim[c];

      for(int r = 0; r < M.rows; ++r)
      {
         if(lhs[r] > b.id[r])
         {
            constraints_cube_free(cons);
            return false;
         }
      }
   }

   constraints_cube_free(cons);
   return true;
}


/* InSimplex
 * Tests whether nodes are inside unit simplex.
 */
bool InSimplex(const_quadrature *q)
{
   int dim = q->params->dim;
   constraints *cons = constraints_simplex_init(q->params->dims);
   get_constraints_simplex(cons);
   Matrix M = cons->M;
   Vector b = cons->b;

   for(int i = 0; i < q->k; ++i)
   {
      double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
      const double *X_ixdim = &q->x[i*dim];

      for(int r = 0; r < M.rows; ++r)
         for(int c = 0; c < M.cols; ++c)
            lhs[r] += M.id[r][c] * X_ixdim[c];

      for(int r = 0; r < M.rows; ++r)
      {
         if(lhs[r] > b.id[r])
         {
            constraints_simplex_free(cons);
            return false;
         }
      }
   }

   constraints_simplex_free(cons);
   return true;
}


/* InCubeSimplex
 * Tests whether nodes are inside (unit cube of dimension[0]) X (unit simplex of dimension[1]).
 */
bool InCubeSimplex(const_quadrature *q)
{
   int dim = q->params->dim;
   constraints *cons = constraints_cubesimplex_init(q->params->dims);
   get_constraints_cubesimplex(cons);
   Matrix M = cons->M;
   Vector b = cons->b;

   for(int i = 0; i < q->k; ++i)
   {
      double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
      const double *X_ixdim = &q->x[i*dim];

      for(int r = 0; r < M.rows; ++r)
         for(int c = 0; c < M.cols; ++c)
            lhs[r] += M.id[r][c] * X_ixdim[c];

      for(int r = 0; r < M.rows; ++r)
      {
         if(lhs[r] > b.id[r])
         {
            constraints_cubesimplex_free(cons);
            return false;
         }
      }
   }

   constraints_cubesimplex_free(cons);
   return true;
}


/* InSimplexSimplex
 * Tests whether nodes are inside (unit simplex of dimension[0]) X (unit simplex of dimension[1]).
 */
bool InSimplexSimplex(const_quadrature *q)
{
   int dim = q->params->dim;
   constraints *cons = constraints_simplexsimplex_init(q->params->dims);
   get_constraints_simplexsimplex(cons);
   Matrix M = cons->M;
   Vector b = cons->b;

   for(int i = 0; i < q->k; ++i)
   {
      double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
      const double *X_ixdim = &q->x[i*dim];

      for(int r = 0; r < M.rows; ++r)
         for(int c = 0; c < M.cols; ++c)
            lhs[r] += M.id[r][c] * X_ixdim[c];

      for(int r = 0; r < M.rows; ++r)
      {
         if(lhs[r] > b.id[r])
         {
            constraints_simplexsimplex_free(cons);
            return false;
         }
      }
   }

   constraints_simplexsimplex_free(cons);
   return true;
}


/* InCubeSimplexSimplex
 * Tests whether nodes are inside (unit simplex of dimension[0]) X (unit simplex of dimension[1]) X (unit cube of dimension[2]).
 */
bool InCubeSimplexSimplex(const_quadrature *q)
{
   int dim = q->params->dim;
   constraints *cons = constraints_cubesimplexsimplex_init(q->params->dims);
   get_constraints_cubesimplexsimplex(cons);
   Matrix M = cons->M;
   Vector b = cons->b;

   for(int i = 0; i < q->k; ++i)
   {
      double lhs[M.rows]; memset(lhs, 0, M.rows*size_double);
      const double *X_ixdim = &q->x[i*dim];

      for(int r = 0; r < M.rows; ++r)
         for(int c = 0; c < M.cols; ++c)
            lhs[r] += M.id[r][c] * X_ixdim[c];

      for(int r = 0; r < M.rows; ++r)
      {
         if(lhs[r] > b.id[r])
         {
            constraints_cubesimplexsimplex_free(cons);
            return false;
         }
      }
   }

   constraints_cubesimplexsimplex_free(cons);
   return true;
}


