/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "Constraints.h"
#include "Vector.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


constraints* constraints_init(void *params, constrInterface *interface)
{
   constraints *constr  = interface->constr_init(params);
   constr->interface    = (constrInterface *)malloc(sizeof(constrInterface));
   *(constr->interface) = *interface;
   return constr;
}

void constraints_get(constraints *constr)
{
   constr->interface->constr_get(constr);
}

void constraints_free(constraints *constr)
{
   constrInterface *interface = constr->interface;
   constr->interface->constr_free(constr);
   if(interface) {
      free(interface);
      interface = NULL;
   }
}

static void get_full_constr(constraints *constr)
{
   RMatrix A = constr->M;
   Vector b  = constr->b;
   RMatrix A_FULL = constr->M_FULL;
   Vector b_FULL  = constr->b_FULL;

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

constrInterface SetIntervalConstrInterface()
{
   constrInterface interface;
   interface.constr_init = (constr_init_ptr)&constraints_interval_init;
   interface.constr_get  = (constr_get_ptr)&get_constraints_interval;
   interface.constr_free = (constr_free_ptr)&constraints_interval_free;
   return interface;
}

constrInterface SetCubeConstrInterface()
{
   constrInterface interface;
   interface.constr_init = (constr_init_ptr)&constraints_cube_init;
   interface.constr_get  = (constr_get_ptr)&get_constraints_cube;
   interface.constr_free = (constr_free_ptr)&constraints_cube_free;
   return interface;
}

constrInterface SetSimplexConstrInterface()
{
   constrInterface interface;
   interface.constr_init = (constr_init_ptr)&constraints_simplex_init;
   interface.constr_get  = (constr_get_ptr)&get_constraints_simplex;
   interface.constr_free = (constr_free_ptr)&constraints_simplex_free;
   return interface;
}

constrInterface SetCubeSimplexConstrInterface()
{
   constrInterface interface;
   interface.constr_init = (constr_init_ptr)&constraints_cubesimplex_init;
   interface.constr_get  = (constr_get_ptr)&get_constraints_cubesimplex;
   interface.constr_free = (constr_free_ptr)&constraints_cubesimplex_free;
   return interface;
}

constrInterface SetSimplexSimplexConstrInterface()
{
   constrInterface interface;
   interface.constr_init = (constr_init_ptr)&constraints_simplexsimplex_init;
   interface.constr_get  = (constr_get_ptr)&get_constraints_simplexsimplex;
   interface.constr_free = (constr_free_ptr)&constraints_simplexsimplex_free;
   return interface;
}


/**************************************
\* Routines for interval constraints \*
**************************************/

constraintsInterval* constraints_interval_init(__attribute__unused dimParamsInterval *params)
{
   constraintsInterval *constr = (constraintsInterval *)malloc(sizeof(constraintsInterval));

   int n_rows   = TWO;
   int n_cols   = ONE;
   constr->dim  = ONE;
   constr->M      = RMatrix_init(n_rows, n_cols);
   constr->b      = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void get_constraints_interval(constraintsInterval *constr)
{
   memset(constr->M.id, 0, SIZE_DOUBLE(constr->M.len));
   memset(constr->b.id, 0, SIZE_DOUBLE(constr->b.len));
   constr->M.rid[0][0] = -1.0;
   constr->M.rid[1][0] = 1.0;
   constr->b.id[0] = 0.0;
   constr->b.id[1] = 1.0;

   get_full_constr((constraints *)constr);
}

void constraints_interval_free(constraintsInterval *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   free(constr);
}



/**********************************
\* Routines for cube constraints \*
**********************************/

constraintsCube* constraints_cube_init(dimParamsCube *params)
{
   assert(params->dim > 0);

   constraintsCube *constr = (constraintsCube *)malloc(sizeof(constraintsCube));
   constr->dimParams = (dimParamsCube *)malloc(sizeof(dimParamsCube));

   int dim = params->dim;
   int n_rows = TWO * dim;
   int n_cols = dim;
   constr->dim = dim;
   constr->dimParams->dim = dim;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void get_constraints_cube(constraintsCube *constr)
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

   get_full_constr((constraints *)constr);
}

void constraints_cube_free(constraintsCube *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   if(constr->dimParams) {
      free(constr->dimParams);
      constr->dimParams = NULL;
   }
   free(constr);
}



/*************************************
\* Routines for simplex constraints \*
*************************************/

constraintsSimplex* constraints_simplex_init(dimParamsSimplex *params)
{
   assert(params->dim > 1);

   constraintsSimplex *constr = (constraintsSimplex *)malloc(sizeof(constraintsSimplex));
   constr->dimParams = (dimParamsSimplex *)malloc(sizeof(dimParamsSimplex));

   int dim = params->dim;
   int n_rows = dim + ONE;
   int n_cols = dim;
   constr->dim = dim;
   constr->dimParams->dim = dim;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void get_constraints_simplex(constraintsSimplex *constr)
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

   get_full_constr((constraints *)constr);
}

void constraints_simplex_free(constraintsSimplex *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   if(constr->dimParams) {
      free(constr->dimParams);
      constr->dimParams = NULL;
   }
   free(constr);
}



/*****************************************
\* Routines for cubesimplex constraints \*
*****************************************/

constraintsCubeSimplex* constraints_cubesimplex_init(dimParamsCubeSimplex *params)
{
   assert((params->dims[0] > 0));
   assert((params->dims[1] > 1));

   constraintsCubeSimplex *constr = (constraintsCubeSimplex *)malloc(sizeof(constraintsCubeSimplex));
   constr->dimParams = (dimParamsCubeSimplex *)malloc(sizeof(dimParamsCubeSimplex));

   int dim = params->dims[0]+params->dims[1];
   int n_rows = (TWO * params->dims[0]) + (params->dims[1] + ONE);
   int n_cols = dim;
   constr->dimParams->dims[0] = params->dims[0];
   constr->dimParams->dims[1] = params->dims[1];
   constr->dim = dim;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void get_constraints_cubesimplex(constraintsCubeSimplex *constr)
{
   int n_rows = constr->M.rows;
   int n_cols = constr->M.cols;
   memset(constr->M.id, 0, SIZE_DOUBLE(n_rows*n_cols));
   memset(constr->b.id, 0, SIZE_DOUBLE(n_rows));

   dimParamsCube paramsCube = { constr->dimParams->dims[0] };
   constraintsCube *constr_cube = constraints_cube_init(&paramsCube);
   get_constraints_cube(constr_cube);

   dimParamsSimplex paramsSimplex = { constr->dimParams->dims[1] };
   constraintsSimplex *constr_simplex = constraints_simplex_init(&paramsSimplex);
   get_constraints_simplex(constr_simplex);

   for(int i = 0; i < constr_cube->M.rows; ++i)
      for(int j = 0; j < constr_cube->M.cols; ++j)
         constr->M.rid[i][j] = constr_cube->M.rid[i][j];
   for(int i = 0; i < constr_simplex->M.rows; ++i)
      for(int j = 0; j < constr_simplex->M.cols; ++j)
         constr->M.rid[constr_cube->M.rows +i][constr_cube->M.cols +j] = constr_simplex->M.rid[i][j];

   for(int i = 0; i < constr_cube->b.len; ++i)
      constr->b.id[i] = constr_cube->b.id[i];
   for(int i = 0; i < constr_simplex->b.len; ++i)
      constr->b.id[constr_cube->b.len+i] = constr_simplex->b.id[i];

   get_full_constr((constraints *)constr);

   constraints_cube_free(constr_cube);
   constraints_simplex_free(constr_simplex);
}

void constraints_cubesimplex_free(constraintsCubeSimplex *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   if(constr->dimParams) {
      free(constr->dimParams);
      constr->dimParams = NULL;
   }
   free(constr);
}



/********************************************
\* Routines for simplexsimplex constraints \*
********************************************/

constraintsSimplexSimplex* constraints_simplexsimplex_init(dimParamsSimplexSimplex *params)
{
   assert((params->dims[0] > 1));
   assert((params->dims[1] > 1));

   constraintsSimplexSimplex *constr= (constraintsSimplexSimplex *)malloc(sizeof(constraintsSimplexSimplex));
   constr->dimParams = (dimParamsSimplexSimplex*)malloc(sizeof(dimParamsSimplexSimplex));

   int dim = params->dims[0] + params->dims[1];
   int n_rows = (params->dims[0] + ONE) + (params->dims[1] + ONE);
   int n_cols = dim;
   constr->dimParams->dims[0] = params->dims[0];
   constr->dimParams->dims[1] = params->dims[1];
   constr->dim = dim;
   constr->M = RMatrix_init(n_rows, n_cols);
   constr->b = Vector_init(n_rows);
   constr->M_FULL = RMatrix_init(n_rows+1, n_cols+1);
   constr->b_FULL = Vector_init(n_rows+1);

   return constr;
}

void get_constraints_simplexsimplex(constraintsSimplexSimplex *constr)
{
   memset(constr->M.id, 0, SIZE_DOUBLE(constr->M.len));
   memset(constr->b.id, 0, SIZE_DOUBLE(constr->b.len));

   dimParamsSimplex paramsS1 = { constr->dimParams->dims[0] };
   constraintsSimplex *constrS1 = constraints_simplex_init(&paramsS1);
   get_constraints_simplex(constrS1);

   dimParamsSimplex paramsS2 = { constr->dimParams->dims[1] };
   constraintsSimplex *constrS2 = constraints_simplex_init(&paramsS2);
   get_constraints_simplex(constrS2);

   for(int i = 0; i < constrS1->M.rows; ++i)
      for(int j = 0; j < constrS1->M.cols; ++j)
         constr->M.rid[i][j] = constrS1->M.rid[i][j];
   for(int i = 0; i < constrS2->M.rows; ++i)
      for(int j = 0; j < constrS2->M.cols; ++j)
         constr->M.rid[constrS1->M.rows +i][constrS1->M.cols +j] = constrS2->M.rid[i][j];

   for(int i = 0; i < constrS1->b.len; ++i)
      constr->b.id[i] = constrS1->b.id[i];
   for(int i = 0; i < constrS2->b.len; ++i)
      constr->b.id[constrS1->b.len+i] = constrS2->b.id[i];

   get_full_constr((constraints *)constr);

   constraints_simplex_free(constrS1);
   constraints_simplex_free(constrS2);
}

void constraints_simplexsimplex_free(constraintsSimplexSimplex *constr)
{
   if(constr == NULL) return;

   RMatrix_free(constr->M);
   Vector_free(constr->b);
   RMatrix_free(constr->M_FULL);
   Vector_free(constr->b_FULL);
   if(constr->dimParams) {
      free(constr->dimParams);
      constr->dimParams = NULL;
   }
   free(constr);
}


void TestConstraints()
{
   printf("Testing CUBE constraints:\n");
   dimParamsCube dimCube = {3};
   constraintsCube *constrCube = constraints_cube_init(&dimCube);
   get_constraints_cube(constrCube);
   PrintConstraints((constraints *)constrCube);
   constraints_cube_free(constrCube);

   printf("Testing SIMPLEX constraints:\n");
   dimParamsSimplex dimSimplex = {3};
   constraintsSimplex *constrSimplex = constraints_simplex_init(&dimSimplex);
   get_constraints_simplex(constrSimplex);
   PrintConstraints((constraints *)constrSimplex);
   constraints_simplex_free(constrSimplex);

   printf("Testing CUBESIMPLEX constraints:\n");
   dimParamsCubeSimplex dimsCS; dimsCS.dims[0] = 2; dimsCS.dims[1] = 2;
   constraintsCubeSimplex *constrCS = constraints_cubesimplex_init(&dimsCS);
   get_constraints_cubesimplex(constrCS);
   PrintConstraints((constraints*)constrCS);
   constraints_cubesimplex_free(constrCS);

   printf("Testing SIMPLEXSIMPLEX constraints:\n");
   dimParamsSimplexSimplex dimsSS; dimsSS.dims[0] = 3; dimsSS.dims[1] = 2;
   constraintsSimplexSimplex *constrSS = constraints_simplexsimplex_init(&dimsSS);
   get_constraints_simplexsimplex(constrSS);
   PrintConstraints((constraints *)constrSS);
   constraints_simplexsimplex_free(constrSS);
}

void PrintConstraints(constraints *constr)
{
   printf("printing constraint Matrix:\n");
   RMatrixPrint(constr->M);
   printf("printing full constraint Matrix:\n");
   RMatrixPrint(constr->M_FULL);
   printf("printing constraint vector:\n");
   VPrint(constr->b);
   printf("printing full constraint vector:\n");
   VPrint(constr->b_FULL);
}

