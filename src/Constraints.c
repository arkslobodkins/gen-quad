/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "Constraints.h"

#include "GENERAL_QUADRATURE.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>




/**************************************
\* Routines for interval constraints \*
***************************************/

constraints *constraints_interval_init(int dims[1])
{
   assert(dims[0] == 1);

   int n_rows = TWO;
   int n_cols = ONE;

   constraints *cons = (constraints *)malloc(sizeof(constraints));
   cons->dim = ONE;
   cons->dims = dims;
   cons->M = Matrix_init(n_rows, n_cols);
   cons->b = Vector_init(n_rows);

   return cons;
}


void constraints_interval_realloc(constraints *constr, int dims[1])
{
   return;
}


void get_constraints_interval(constraints *cons)
{
   int n_rows = cons->M.rows;
   int n_cols = cons->M.cols;

   for(int i = 0; i < n_rows; ++i) memset(cons->M.id[i], 0, SIZE_DOUBLE(n_cols));
   memset(cons->b.id, 0, SIZE_DOUBLE(n_rows));

   cons->M.id[0][0] = -1.0;
   cons->M.id[1][0] = 1.0;

   cons->b.id[0] = 0.0;
   cons->b.id[1] = 1.0;
}


void constraints_interval_free(constraints *constr)
{
	Matrix_free(constr->M);
	Vector_free(constr->b);
   free(constr);
}




/**********************************
\* Routines for cube constraints \*
**********************************/

constraints *constraints_cube_init(int dims[1])
{
   assert(dims[0] > 0);

   int dim = dims[0];
   int n_rows = TWO * dim;
   int n_cols = dim;

   constraints *cons = (constraints *)malloc(sizeof(constraints));
   cons->dim = dim;
   cons->dims = dims;
   cons->M = Matrix_init(n_rows, n_cols);
   cons->b = Vector_init(n_rows);

   return cons;
}


void constraints_cube_realloc(constraints *cons, int dims[1])
{
   assert(dims[0] > 0);
   int dim = dims[0];
   int n_rows = TWO * dim;
   int n_cols = dim;

   Matrix_realloc(n_rows, n_cols, &cons->M);
   Vector_realloc(n_rows, &cons->b);
}


void get_constraints_cube(constraints *cons)
{
   int n_rows = cons->M.rows;
   int n_cols = cons->M.cols;

   for(int i = 0; i < n_rows; ++i) memset(cons->M.id[i], 0, SIZE_DOUBLE(n_cols));
   memset(cons->b.id, 0, SIZE_DOUBLE(n_rows));

   for(int j = 0; j < n_cols; ++j)
   {
      int two_j = 2*j;
      cons->M.id[two_j][j] = -1.0;
      cons->M.id[two_j+1][j] = 1.0;
   }

   for(int i = 0; i < n_rows; i += 2)
   {
      cons->b.id[i] = 0.0;
      cons->b.id[i+1] = 1.0;
   }
}


void constraints_cube_free(constraints *constr)
{
	Matrix_free(constr->M);
	Vector_free(constr->b);
   free(constr);
}




/*************************************
\* Routines for simplex constraints \*
*************************************/

constraints *constraints_simplex_init(int dims[1])

{
   assert(dims[0] > 1);

   int dim = dims[0];
   int n_rows = dim + ONE;
   int n_cols = dim;

   constraints *cons = (constraints *)malloc(sizeof(constraints));
   cons->dim = dim;
   cons->dims = dims;
   cons->M = Matrix_init(n_rows, n_cols);
   cons->b = Vector_init(n_rows);

   return cons;
}


void constraints_simplex_realloc(constraints *cons, int dims[1])
{
   assert(dims[0] > 1);
   int dim = dims[0];
   int n_rows = dim + ONE;
   int n_cols = dim;

   Matrix_realloc(n_rows, n_cols, &cons->M);
   Vector_realloc(n_rows, &cons->b);
}


void get_constraints_simplex(constraints *cons)
{

   int n_rows = cons->M.rows;
   int n_cols = cons->M.cols;

   for(int i = 0; i < n_rows; ++i) memset(cons->M.id[i], 0, SIZE_DOUBLE(n_cols));
   memset(cons->b.id, 0, SIZE_DOUBLE(n_rows));

   cons->b.id[0] = 1.0;

   cons->M.id[0][0] = 1.0;
   for(int j = 1; j < n_rows-1; ++j)
   {
      cons->M.id[j][j-1] = -1.0;
      cons->M.id[j][j] = +1.0;
   }
   cons->M.id[n_rows-1][n_cols-1] = -1.0;
}


void constraints_simplex_free(constraints *constr)
{
   Matrix_free(constr->M);
   Vector_free(constr->b);
   free(constr);
}




/*****************************************
\* Routines for cubesimplex constraints \*
*****************************************/

constraints *constraints_cubesimplex_init(int dims[2])
{
   assert( (dims[0] > 0) && (dims[1] > 1) );

   int dim = dims[0]+dims[1];

   constraints *cons = (constraints *)malloc(sizeof(constraints));
   cons->dim = dim;
   cons->dims = dims;

   int n_rows = (TWO * dims[0]) + (dims[1] + ONE);
   int n_cols = dim;
   cons->M = Matrix_init(n_rows, n_cols);
   cons->b = Vector_init(n_rows);

   return cons;
}


void constraints_cubesimplex_realloc(constraints *cons, int dims[2])
{
   assert( (dims[0] > 0) && (dims[1] > 1) );
   int dim = dims[0]+dims[1];
   int n_rows = (TWO * dims[0]) + (dims[1] + ONE);
   int n_cols = dim;

   Matrix_realloc(n_rows, n_cols, &cons->M);
   Vector_realloc(n_rows, &cons->b);
}


void get_constraints_cubesimplex(constraints *cons)
{

   int n_rows = cons->M.rows;
   int n_cols = cons->M.cols;

   for(int i = 0; i < n_rows; ++i) memset(cons->M.id[i], 0, SIZE_DOUBLE(n_cols));
   memset(cons->b.id, 0, SIZE_DOUBLE(n_rows));

   int dim_cube = cons->dims[0];
   constraints *cons_cube = constraints_cube_init(&dim_cube);
   get_constraints_cube(cons_cube);

   for(int i = 0; i < cons_cube->M.rows; ++i)
      for(int j = 0; j < cons_cube->M.cols; ++j)
         cons->M.id[i][j] = cons_cube->M.id[i][j];

   for(int i = 0; i < cons_cube->b.len; ++i)
      cons->b.id[i] = cons_cube->b.id[i];


   int dim_simplex = cons->dims[1];
   constraints *cons_simplex = constraints_simplex_init(&dim_simplex);
   get_constraints_simplex(cons_simplex);

   for(int i = 0; i < cons_simplex->M.rows; ++i)
      for(int j = 0; j < cons_simplex->M.cols; ++j)
         cons->M.id[cons_cube->M.rows +i][cons_cube->M.cols +j] = cons_simplex->M.id[i][j];

   for(int i = 0; i < cons_simplex->b.len; ++i)
      cons->b.id[cons_cube->b.len+i] = cons_simplex->b.id[i];


   constraints_cube_free(cons_cube);
   constraints_simplex_free(cons_simplex);
}


void constraints_cubesimplex_free(constraints *constr)
{
   Matrix_free(constr->M);
   Vector_free(constr->b);
   free(constr);
}




/********************************************
\* Routines for simplexsimplex constraints \*
********************************************/

constraints *constraints_simplexsimplex_init(int dims[2])
{
   assert( (dims[0] > 1) && (dims[1] > 1) );

   int dim = dims[0] + dims[1];
   constraints *cons = (constraints *)malloc(sizeof(constraints));
   cons->dim = dim;
   cons->dims = dims;

   int n_rows = (dims[0] + ONE) + (dims[1] + ONE);
   int n_cols = dim;
   cons->M = Matrix_init(n_rows, n_cols);
   cons->b = Vector_init(n_rows);

   return cons;
}


void constraints_simplexsimplex_realloc(constraints *cons, int dims[2])
{
   assert( (dims[0] > 1) && (dims[1] > 1) );
   int dim = dims[0] + dims[1];
   int n_rows = (dims[0] + ONE) + (dims[1] + ONE);
   int n_cols = dim;

   Matrix_realloc(n_rows, n_cols, &cons->M);
   Vector_realloc(n_rows, &cons->b);
}


void get_constraints_simplexsimplex(constraints *cons)
{

   int n_rows = cons->M.rows;
   int n_cols = cons->M.cols;

   for(int i = 0; i < n_rows; ++i) memset(cons->M.id[i], 0, SIZE_DOUBLE(n_cols));
   memset(cons->b.id, 0, SIZE_DOUBLE(n_rows));

   int dim_s1 = cons->dims[0];
   constraints *cons_s1 = constraints_simplex_init(&dim_s1);
   get_constraints_simplex(cons_s1);

   for(int i = 0; i < cons_s1->M.rows; ++i)
      for(int j = 0; j < cons_s1->M.cols; ++j)
         cons->M.id[i][j] = cons_s1->M.id[i][j];

   for(int i = 0; i < cons_s1->b.len; ++i)
      cons->b.id[i] = cons_s1->b.id[i];


   int dim_s2 = cons->dims[1];
   constraints *cons_s2 = constraints_simplex_init(&dim_s2);
   get_constraints_simplex(cons_s2);

   for(int i = 0; i < cons_s2->M.rows; ++i)
      for(int j = 0; j < cons_s2->M.cols; ++j)
         cons->M.id[cons_s1->M.rows +i][cons_s1->M.cols +j] = cons_s2->M.id[i][j];

   for(int i = 0; i < cons_s2->b.len; ++i)
      cons->b.id[cons_s1->b.len+i] = cons_s2->b.id[i];


   constraints_simplex_free(cons_s1);
   constraints_simplex_free(cons_s2);

}


void constraints_simplexsimplex_free(constraints *constr)
{
   Matrix_free(constr->M);
   Vector_free(constr->b);
   free(constr);
}




/************************************************
\* Routines for cubesimplexsimplex constraints \*
************************************************/

constraints *constraints_cubesimplexsimplex_init(int dims[3])
{
   assert( (dims[0] > 1) && (dims[1] > 1)  && (dims[2] > 0) );

	int dim = dims[0] + dims[1] + dims[2];
   constraints *cons = (constraints *)malloc(sizeof(constraints));
   cons->dim = dim;
   cons->dims = dims;

   int n_rows =  (dims[0] + ONE) + (dims[1]+ ONE) + (TWO * dims[2]); // double check the order
   int n_cols = dim;
	cons->M = Matrix_init(n_rows, n_cols);
   cons->b = Vector_init(n_rows);

   return cons;
}


void constraints_cubesimplexsimplex_realloc(constraints *cons, int dims[3])
{
   assert( (dims[0] > 1) && (dims[1] > 1)  && (dims[2] > 0) );
	int dim = dims[0] + dims[1] + dims[2];
   int n_rows =  (dims[0] + ONE) + (dims[1]+ ONE) + (TWO * dims[2]); // double check the order
   int n_cols = dim;

   Matrix_realloc(n_rows, n_cols, &cons->M);
   Vector_realloc(n_rows, &cons->b);
}


void get_constraints_cube_simplex_simplex(constraints *cons)
{

}


void get_constraints_cubesimplexsimplex(constraints *cons)
{

}


void constraints_cubesimplexsimplex_free(constraints *constr)
{
   Matrix_free(constr->M);
   Vector_free(constr->b);
   free(constr);
}

