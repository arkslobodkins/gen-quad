#ifndef GENERAL_QUADRATURE_H
#define GENERAL_QUADRATURE_H

#include "Matrix.h"
#include "Vector.h"
#include "glist.h"

#include <stdbool.h>
#include <stdint.h>
#include <math.h>


#ifdef __cplusplus
extern "C" {
#endif

// Set debug mode to ON by default. Compile with -DQUAD_DEBUG_OFF to disable debugging mode.
#ifndef QUAD_DEBUG_OFF
#define QUAD_DEBUG_ON
#endif

#ifdef __GNUC__
#define ATTR_UNUSED __attribute__ ((unused))
#else
#define ATTR_UNUSED
#endif


#define ONE 1
#define TWO 2
#define THREE 3
#define QUAD_HUGE 10.0
#define CLOSE_TO_ZERO POW_DOUBLE(10, -18)
#define QUAD_TOL  POW_DOUBLE(10.0, -15)
#define BOUND_TOL POW_DOUBLE(10.0, -12)
#define POW_DOUBLE(x,y) pow(x,y)
#define POW_INT(x,y)    IntPower(x,y)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define is_greater_than_zero(x) ((x) > 0 ? true : false)
#define SQUARE(x) ((x)*(x))
#define SQRT(x) sqrt((double)(x))
#define ij2(i, j, num_cols) ((i)*(num_cols)+j) // maps indices from 2-d layout to memory
#define size_int sizeof(int)
#define size_double sizeof(double)
#define SIZE_INT(x) ((x)*sizeof(int))
#define SIZE_DOUBLE(x) ((x)*sizeof(double))
#define size_quadrature sizeof(quadrature)


typedef enum { ON, OFF } bool_enum;
typedef enum { INTERVAL, CUBE, SIMPLEX, CUBESIMPLEX, SIMPLEXSIMPLEX, CUBESIMPLEXSIMPLEX } DOMAIN_TYPE;
typedef enum { NODE, WEIGHT, NONE } NODE_OR_WEIGHT;


typedef struct
{
   int tot_elims;
   int *nodes_tot;
   int *success_node;
   int *success_its;
} elim_history;

typedef struct
{
   int nodes_tot;
   int success_node;
   int success_its;
} hist_data;

typedef struct
{
   glist *list;
   int dim;
   int degree;
   int nodes_initial;
   int nodes_final;
   int num_funcs;
   double res;
   DOMAIN_TYPE D;
} history;

typedef struct
{
   int            nodeId;
   int            eqnId;
   double         tMin;
   bool_enum      ACTIVE;
   NODE_OR_WEIGHT N_OR_W;
} ConstrNodeData;


typedef struct
{
   int            boundaryNodeId;
   int            eqnId;
   double         tMin;
   bool_enum      ACTIVE;
   NODE_OR_WEIGHT N_OR_W;
} ConstrVectData;


typedef struct quadrature quadrature;
typedef struct const_quadrature const_quadrature;

typedef void   (*EvalBasis)       (int *dims, int deg, const int_fast8_t *basis, const double *x, double *phi);
typedef void   (*EvalBasisDer)    (int *dims, int deg, const int_fast8_t *basis, const double *x, double *phiPrime);
typedef void(*BasisIntegrals)(int *dims, int deg, double *integrals);
typedef bool(*InDomain)(const_quadrature *quad);
typedef bool(*InDomainElem)(const_quadrature *quad, int elem);
typedef bool(*InConstraint)(const_quadrature *q);
typedef bool(*InConstraintElem)(const_quadrature *quad, int elem);
typedef bool(*PosWeights)(const_quadrature *q);
typedef bool(*PosWeightsElem)(const_quadrature *q, int elem);
typedef bool(*OnTheBoundary)(const_quadrature *q, int elem);
typedef bool(*EqnOnTheBoundary)(const_quadrature *q, int elem, int eqn);
typedef double(*TestIntegral)(const_quadrature *q);
typedef void(*SetFuncs)(quadrature *q);
typedef void(*SetParams)(int dim, int num_dims, int *dims, int deg, quadrature *q);
typedef void(*FreePtr)(quadrature *quad);

typedef struct constraints constraints;
typedef constraints*(*constraints_init)(int *dims);
typedef void (*constraints_realloc)(constraints *cons, int *dims);
typedef void (*get_constraints)    (constraints *cons);
typedef void (*constraints_free)   (constraints *cons);



struct constraints
{
   int dim;
   int *dims;
   RMatrix M;
   Vector b;
   RMatrix M_FULL;
   Vector b_FULL;
};


struct quadrature
{
   DOMAIN_TYPE D;
   int dim;
   int num_dims;
   int *dims;
   int deg;
   int num_funcs;

   int k;
   double *w;
   double *x;
   Vector z;
   constraints *constr;

   int setFuncsConstrFlag;
   SetFuncs setFuncs; // Function that sets up all domain functions to point to appropriate domain

   EvalBasis           evalBasis;
   EvalBasisDer        evalBasisDer;
   BasisIntegrals      basisIntegrals;
   constraints_init    constr_init;
   constraints_realloc constr_realloc;
   get_constraints     get_constr;
   constraints_free    constr_free;

   FreePtr free_ptr;
};


#define Q_SUCCESS       0
#define NULL_VAL       -1
#define ALLOC_FAIL     -2
#define ILL_INPUT      -3
#define INF_VAL        -4
#define NAN_VAL        -5
#define QUAD_HUGE_ERR  -6
#define LAPACK_ERR     -7
#define DIV_BY_ZERO    -8


#ifdef __cplusplus
}
#endif

#endif

