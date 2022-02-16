#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct constraints constraints;

typedef constraints*(*constr_init_ptr)(void *params);
typedef void(*constr_get_ptr)(constraints *constr);
typedef void(*constr_free_ptr)(constraints *constr);

typedef struct
{
   constr_init_ptr constr_init;
   constr_get_ptr constr_get;
   constr_free_ptr constr_free;
} constrInterface;

typedef struct
{
   void *empty;
} dimParamsInterval;

typedef struct
{
   int dim;
} dimParamsCube;

typedef struct
{
   int dim;
} dimParamsSimplex;

typedef struct
{
   int dims[2];
} dimParamsCubeSimplex;

typedef struct
{
   int dims[2];
} dimParamsSimplexSimplex;

struct constraints
{
   constrInterface *interface;
   void *dimParams;
   int dim;
   RMatrix M;
   Vector b;
   RMatrix M_FULL;
   Vector b_FULL;
};

typedef struct
{
   constrInterface *interface;
   dimParamsInterval *dimParams;
   int dim;
   RMatrix M;
   Vector b;
   RMatrix M_FULL;
   Vector b_FULL;
} constraintsInterval;

typedef struct
{
   constrInterface *interface;
   dimParamsCube *dimParams;
   int dim;
   RMatrix M;
   Vector b;
   RMatrix M_FULL;
   Vector b_FULL;
} constraintsCube;

typedef struct
{
   constrInterface *interface;
   dimParamsSimplex *dimParams;
   int dim;
   RMatrix M;
   Vector b;
   RMatrix M_FULL;
   Vector b_FULL;
} constraintsSimplex;

typedef struct
{
   constrInterface *interface;
   dimParamsCubeSimplex *dimParams;
   int dim;
   RMatrix M;
   Vector b;
   RMatrix M_FULL;
   Vector b_FULL;
} constraintsCubeSimplex;

typedef struct
{
   constrInterface *interface;
   dimParamsSimplexSimplex *dimParams;
   int dim;
   RMatrix M;
   Vector b;
   RMatrix M_FULL;
   Vector b_FULL;
} constraintsSimplexSimplex;

constrInterface SetIntervalConstrInterface();
constrInterface SetCubeConstrInterface();
constrInterface SetSimplexConstrInterface();
constrInterface SetCubeSimplexConstrInterface();
constrInterface SetSimplexSimplexConstrInterface();

constraints* constraints_init(void *params, constrInterface *interface);
void constraints_get(constraints *constr);
void constraints_free(constraints *constr);

constraintsInterval* constraints_interval_init(dimParamsInterval *params);
void get_constraints_interval(constraintsInterval *constr);
void constraints_interval_free(constraintsInterval *constr);

constraintsCube* constraints_cube_init(dimParamsCube *params);
void get_constraints_cube(constraintsCube *constr);
void constraints_cube_free(constraintsCube *constr);

constraintsSimplex* constraints_simplex_init(dimParamsSimplex *params);
void get_constraints_simplex(constraintsSimplex *constr);
void constraints_simplex_free(constraintsSimplex *constr);

constraintsCubeSimplex* constraints_cubesimplex_init(dimParamsCubeSimplex *params);
void get_constraints_cubesimplex(constraintsCubeSimplex *constr);
void constraints_cubesimplex_free(constraintsCubeSimplex *constr);

constraintsSimplexSimplex* constraints_simplexsimplex_init(dimParamsSimplexSimplex *params);
void get_constraints_simplexsimplex(constraintsSimplexSimplex *constr);
void constraints_simplexsimplex_free(constraintsSimplexSimplex *constr);

void TestConstraints();
void PrintConstraints(constraints *constr);

#ifdef __cplusplus
}
#endif

#endif
