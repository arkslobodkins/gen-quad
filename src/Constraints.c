/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "GENERAL_QUADRATURE.h"
#include "Constraints.h"
#include "Vector.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

static void get_full_constr(constraints *constr);


/**************************************
\* Routines for interval constraints \*
***************************************/

constraints *constraints_interval_init(int dims[1])
{
   assert(dims[0] == 1);

   constraints *constr = (constraints *)malloc(sizeof(constraints));
   if(constr == NULL) return NULL;

   int n_rows = TWO;
   int n_cols = ONE;
   constr->dim = ONE;
   constr->dims = dims;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void constraints_interval_realloc(constraints *constr, int dims[1])
{
   return; // nothing to be reallocated, only defined to stay consistent with other domains
}

void get_constraints_interval(constraints *constr)
{
   memset(constr->M.id, 0, SIZE_DOUBLE(constr->M.len));
   memset(constr->b.id, 0, SIZE_DOUBLE(constr->b.len));
   constr->M.rid[0][0] = -1.0;
   constr->M.rid[1][0] = 1.0;
   constr->b.id[0] = 0.0;
   constr->b.id[1] = 1.0;

   get_full_constr(constr);
}

void constraints_interval_free(constraints *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   free(constr); constr = NULL;
}



/**********************************
\* Routines for cube constraints \*
**********************************/

constraints *constraints_cube_init(int dims[1])
{
   assert(dims[0] > 1);

   constraints *constr = (constraints *)malloc(sizeof(constraints));
   if(constr == NULL) return NULL;

   int dim = dims[0];
   int n_rows = TWO * dim;
   int n_cols = dim;
   constr->dim = dim;
   constr->dims = dims;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void constraints_cube_realloc(constraints *constr, int dims[1])
{
   assert(dims[0] > 1);

   int dim = dims[0];
   int n_rows = TWO * dim;
   int n_cols = dim;
   RMatrix_realloc(n_rows, n_cols, &constr->M);
   Vector_realloc(n_rows, &constr->b);
   RMatrix_realloc(n_rows+1, n_cols+1, &constr->M_FULL);
   Vector_realloc(n_rows+1, &constr->b_FULL);
}

void get_constraints_cube(constraints *constr)
{
   int n_rows = constr->M.rows;
   int n_cols = constr->M.cols;
   memset(constr->M.id, 0, SIZE_DOUBLE(n_rows*n_cols));
   memset(constr->b.id, 0, SIZE_DOUBLE(n_rows));

   for(int j = 0; j < n_cols; ++j) {
      int two_j = 2*j;
      constr->M.rid[two_j][j] = -1.0;
      constr->M.rid[two_j+1][j] = 1.0;
   }
   for(int i = 0; i < n_rows; i += 2) {
      constr->b.id[i] = 0.0;
      constr->b.id[i+1] = 1.0;
   }

   get_full_constr(constr);
}

void constraints_cube_free(constraints *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   free(constr); constr = NULL;
}



/*************************************
\* Routines for simplex constraints \*
*************************************/

constraints *constraints_simplex_init(int dims[1])
{
   assert(dims[0] > 1);

   constraints *constr = (constraints *)malloc(sizeof(constraints));
   if(constr == NULL) return NULL;

   int dim = dims[0];
   int n_rows = dim + ONE;
   int n_cols = dim;
   constr->dim = dim;
   constr->dims = dims;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void constraints_simplex_realloc(constraints *constr, int dims[1])
{
   assert(dims[0] > 1);

   int dim = dims[0];
   int n_rows = dim + ONE;
   int n_cols = dim;
   RMatrix_realloc(n_rows, n_cols, &constr->M);
   Vector_realloc(n_rows, &constr->b);
   RMatrix_realloc(n_rows+1, n_cols+1, &constr->M_FULL);
   Vector_realloc(n_rows+1, &constr->b_FULL);
}

void get_constraints_simplex(constraints *constr)
{
   int n_rows = constr->M.rows;
   int n_cols = constr->M.cols;

   memset(constr->M.id, 0, SIZE_DOUBLE(n_rows*n_cols));
   memset(constr->b.id, 0, SIZE_DOUBLE(n_rows));
   constr->b.id[0] = 1.0;

   constr->M.rid[0][0] = 1.0;
   for(int j = 1; j < n_rows-1; ++j) {
      constr->M.rid[j][j-1] = -1.0;
      constr->M.rid[j][j] = +1.0;
   }
   constr->M.rid[n_rows-1][n_cols-1] = -1.0;

   get_full_constr(constr);
}

void constraints_simplex_free(constraints *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   free(constr); constr = NULL;
}



/*****************************************
\* Routines for cubesimplex constraints \*
*****************************************/

constraints *constraints_cubesimplex_init(int dims[2])
{
   assert( (dims[0] > 0) && (dims[1] > 1) );

   constraints *constr = (constraints *)malloc(sizeof(constraints));
   if(constr == NULL) return NULL;

   int dim = dims[0]+dims[1];
   int n_rows = (TWO * dims[0]) + (dims[1] + ONE);
   int n_cols = dim;
   constr->dim = dim;
   constr->dims = dims;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void constraints_cubesimplex_realloc(constraints *constr, int dims[2])
{
   assert( (dims[0] > 0) && (dims[1] > 1) );

   int dim = dims[0]+dims[1];
   int n_rows = (TWO * dims[0]) + (dims[1] + ONE);
   int n_cols = dim;
   RMatrix_realloc(n_rows, n_cols, &constr->M);
   Vector_realloc(n_rows, &constr->b);
   RMatrix_realloc(n_rows+1, n_cols+1, &constr->M_FULL);
   Vector_realloc(n_rows+1, &constr->b_FULL);
}

void get_constraints_cubesimplex(constraints *constr)
{
   int n_rows = constr->M.rows;
   int n_cols = constr->M.cols;
   memset(constr->M.id, 0, SIZE_DOUBLE(n_rows*n_cols));
   memset(constr->b.id, 0, SIZE_DOUBLE(n_rows));

   int dim_cube = constr->dims[0];
   constraints *cons_cube = constraints_cube_init(&dim_cube);
   get_constraints_cube(cons_cube);

   for(int i = 0; i < cons_cube->M.rows; ++i)
      for(int j = 0; j < cons_cube->M.cols; ++j)
         constr->M.rid[i][j] = cons_cube->M.rid[i][j];

   for(int i = 0; i < cons_cube->b.len; ++i)
      constr->b.id[i] = cons_cube->b.id[i];



   int dim_simplex = constr->dims[1];
   constraints *cons_simplex = constraints_simplex_init(&dim_simplex);
   get_constraints_simplex(cons_simplex);

   for(int i = 0; i < cons_simplex->M.rows; ++i)
      for(int j = 0; j < cons_simplex->M.cols; ++j)
         constr->M.rid[cons_cube->M.rows +i][cons_cube->M.cols +j] = cons_simplex->M.rid[i][j];

   for(int i = 0; i < cons_simplex->b.len; ++i)
      constr->b.id[cons_cube->b.len+i] = cons_simplex->b.id[i];


   get_full_constr(constr);

   constraints_cube_free(cons_cube);
   constraints_simplex_free(cons_simplex);
}

void constraints_cubesimplex_free(constraints *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   free(constr); constr = NULL;
}



/********************************************
\* Routines for simplexsimplex constraints \*
********************************************/

constraints *constraints_simplexsimplex_init(int dims[2])
{
   assert( (dims[0] > 1) && (dims[1] > 1) );

   constraints *constr = (constraints *)malloc(sizeof(constraints));
   if(constr == NULL) return NULL;

   int dim = dims[0] + dims[1];
   int n_rows = (dims[0] + ONE) + (dims[1] + ONE);
   int n_cols = dim;
   constr->dim = dim;
   constr->dims = dims;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void constraints_simplexsimplex_realloc(constraints *constr, int dims[2])
{
   assert( (dims[0] > 1) && (dims[1] > 1) );

   int dim = dims[0] + dims[1];
   int n_rows = (dims[0] + ONE) + (dims[1] + ONE);
   int n_cols = dim;
   RMatrix_realloc(n_rows, n_cols, &constr->M);
   Vector_realloc(n_rows, &constr->b);
   RMatrix_realloc(n_rows+1, n_cols+1, &constr->M_FULL);
   Vector_realloc(n_rows+1, &constr->b_FULL);
}

void get_constraints_simplexsimplex(constraints *constr)
{
   memset(constr->M.id, 0, SIZE_DOUBLE(constr->M.len));
   memset(constr->b.id, 0, SIZE_DOUBLE(constr->b.len));

   int dim_s1 = constr->dims[0];
   constraints *cons_s1 = constraints_simplex_init(&dim_s1);
   get_constraints_simplex(cons_s1);

   for(int i = 0; i < cons_s1->M.rows; ++i)
      for(int j = 0; j < cons_s1->M.cols; ++j)
         constr->M.rid[i][j] = cons_s1->M.rid[i][j];

   for(int i = 0; i < cons_s1->b.len; ++i)
      constr->b.id[i] = cons_s1->b.id[i];



   int dim_s2 = constr->dims[1];
   constraints *cons_s2 = constraints_simplex_init(&dim_s2);
   get_constraints_simplex(cons_s2);

   for(int i = 0; i < cons_s2->M.rows; ++i)
      for(int j = 0; j < cons_s2->M.cols; ++j)
         constr->M.rid[cons_s1->M.rows +i][cons_s1->M.cols +j] = cons_s2->M.rid[i][j];

   for(int i = 0; i < cons_s2->b.len; ++i)
      constr->b.id[cons_s1->b.len+i] = cons_s2->b.id[i];


   get_full_constr(constr);

   constraints_simplex_free(cons_s1);
   constraints_simplex_free(cons_s2);
}

void constraints_simplexsimplex_free(constraints *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   free(constr); constr = NULL;
}



/************************************************
\* Routines for cubesimplexsimplex constraints \*
************************************************/

constraints *constraints_cubesimplexsimplex_init(int dims[3])
{
   assert( (dims[0] > 1) && (dims[1] > 1)  && (dims[2] > 0) );

   constraints *constr = (constraints *)malloc(sizeof(constraints));
   if(constr == NULL) return NULL;

   int dim = dims[0] + dims[1] + dims[2];
   int n_rows =  (dims[0] + ONE) + (dims[1]+ ONE) + (TWO * dims[2]); // double check the order
   int n_cols = dim;
   constr->dim = dim;
   constr->dims = dims;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void constraints_cubesimplexsimplex_realloc(constraints *constr, int dims[3])
{
   assert( (dims[0] > 1) && (dims[1] > 1)  && (dims[2] > 0) );

   int dim = dims[0] + dims[1] + dims[2];
   int n_rows =  (dims[0] + ONE) + (dims[1]+ ONE) + (TWO * dims[2]); // double check the order
   int n_cols = dim;
   RMatrix_realloc(n_rows, n_cols, &constr->M);
   Vector_realloc(n_rows, &constr->b);
   RMatrix_realloc(n_rows+1, n_cols+1, &constr->M_FULL);
   Vector_realloc(n_rows+1, &constr->b_FULL);
}

void get_constraints_cube_simplex_simplex(constraints *constr)
{
}

void get_constraints_cubesimplexsimplex(constraints *constr)
{
}

void constraints_cubesimplexsimplex_free(constraints *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   free(constr); constr = NULL;
}

static void get_full_constr(constraints *constr)
{
   RMatrix A = constr->M;
   Vector b = constr->b;
   RMatrix A_FULL = constr->M_FULL;
   Vector b_FULL = constr->b_FULL;

   for(int i = 0; i < A.rows; ++i)
      for(int j = 0; j < A.cols; ++j)
         A_FULL.rid[i+1][j+1] = A.rid[i][j];
   A_FULL.rid[0][0] = -1.0;

   for(int i = 0; i < b.len; ++i)
      b_FULL.id[i+1] = b.id[i];
   b_FULL.id[0] = 0.0;

   constr->M_FULL = A_FULL;
   constr->b_FULL = b_FULL;
}


void TestConstraints()
{
   printf("Testing CUBE constraints:\n");
   int dimCube = 3;
   constraints *constrCube = constraints_cube_init(&dimCube);
   get_constraints_cube(constrCube);
   PrintConstraints(constrCube);
   constraints_cube_free(constrCube);

   printf("Testing SIMPLEX constraints:\n");
   int dimSimplex = 3;
   constraints *constrSimplex = constraints_simplex_init(&dimSimplex);
   get_constraints_simplex(constrSimplex);
   PrintConstraints(constrSimplex);
   constraints_simplex_free(constrSimplex);

   printf("Testing CUBESIMPLEX constraints:\n");
   int dimsCS[2]; dimsCS[0] = 2; dimsCS[1] = 2;
   constraints *constrCS = constraints_cubesimplex_init(dimsCS);
   get_constraints_cubesimplex(constrCS);
   PrintConstraints(constrCS);
   constraints_cubesimplex_free(constrCS);

   printf("Testing SIMPLEXSIMPLEX constraints:\n");
   int dimsSS[2]; dimsSS[0] = 3; dimsSS[1] = 2;
   constraints *constrSS = constraints_simplexsimplex_init(dimsSS);
   get_constraints_simplexsimplex(constrSS);
   PrintConstraints(constrSS);
   constraints_simplex_free(constrSS);
}

void PrintConstraints(constraints *constr)
{
   printf("printing constraint Matrix:\n");
   PrintRMatrix(constr->M);
   printf("printing full constraint Matrix:\n");
   PrintRMatrix(constr->M_FULL);
   printf("printing constraint vector:\n");
   VPrint(constr->b);
   printf("printing full constraint vector:\n");
   VPrint(constr->b_FULL);
}
