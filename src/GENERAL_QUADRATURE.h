#ifndef GENERAL_QUADRATURE_H
#define GENERAL_QUADRATURE_H

#include "Matrix.h"
#include "Vector.h"

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
#define POW(x,y) pow(x,y)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define SQUARE(x) ((x)*(x))
#define SQRT(x) sqrt((double)(x))
#define QUAD_TOL POW(10, -15)
#define ij2(i, j, num_cols) i*num_cols+j // maps indices from 2-d layout to memory
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
   NODE_OR_WEIGHT N_OR_W;
   int nodeId;
   int eqnId;
   double tMin;
   int flag;
} ConstrNodeData;


typedef struct
{
   bool_enum ACTIVE_CONSTRAINTS;
   NODE_OR_WEIGHT N_OR_W;
   int boundaryNodeId;
   int  eqnId;
} ConstrVectData;


typedef struct quadrature quadrature;
typedef struct const_quadrature const_quadrature;

typedef void(*EvalBasis)(int *dims, int deg, const int_fast8_t *basis, const double *x, double *phi);
typedef void(*EvalBasisDer)(int *dims, int deg, const int_fast8_t *basis, const double *x, double *phiPrime);
typedef void(*BasisIntegrals)(int *dims, int deg, double *integrals);
typedef bool(*InDomain)(const_quadrature *quad);
typedef bool(*InDomainElem)(const_quadrature *quad, int elem);
typedef bool (*InConstraint)(const_quadrature *q);
typedef bool (*InConstraintElem)(const_quadrature *quad, int elem);
typedef bool (*PosWeights)(const_quadrature *q);
typedef bool (*PosWeightsElem)(const_quadrature *q, int elem);
typedef bool (*OnTheBoundary)(const_quadrature *q, int elem);
typedef bool (*EqnOnTheBoundary)(const_quadrature *q, int elem, int eqn);
typedef double(*TestIntegral)(const_quadrature *q);

typedef struct constraints constraints;
typedef constraints*(*constraints_init)(int *dims);
typedef void (*constraints_realloc)(constraints *cons, int *dims);
typedef void (*get_constraints)(constraints *cons);
typedef void (*constraints_free)(constraints *cons);

typedef void(*SetFuncs)(quadrature *q);
typedef void (*SetParams)(int dim, int num_dims, int *dims, int deg, quadrature *q);


struct constraints
{
   int dim;
   int *dims;
   Matrix M;
   Vector b;
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
   double *z;
   constraints *cons;
   Matrix FULL_A;
   Vector FULL_b;

   SetFuncs setFuncs; // Function that sets up all domain functions to point to appropriate domain
   int setFuncsFlag;

   EvalBasis evalBasis;
   EvalBasisDer evalBasisDer;
   BasisIntegrals basisIntegrals;
   InDomain inDomain;
   InDomainElem inDomainElem;
   InConstraint inConstraint;
   InConstraintElem inConstraintElem;
   PosWeights posWeights;
   PosWeightsElem posWeightsElem;
   OnTheBoundary onTheBoundary;
   EqnOnTheBoundary eqnOnTheBoundary;
   TestIntegral testIntegral;
   constraints_init constr_init;
   constraints_realloc constr_realloc;
   get_constraints get_constr;
   constraints_free constr_free;
};

struct const_quadrature
{
   const DOMAIN_TYPE D;
   const int dim;
   const int num_dims;
   const int *dims;
   const int deg;
   const int num_funcs;

   const int k;
   const double *w;
   const double *x;
   const double *z;
   const constraints *cons;
   const Matrix FULL_A;
   const Vector FULL_b;

   const SetFuncs setFuncs; // Function that sets up all domain functions to point to appropriate domain
   int setFuncsFlag;

   const EvalBasis evalBasis;
   const EvalBasisDer evalBasisDer;
   const BasisIntegrals basisIntegrals;
   const InDomain inDomain;
   const InDomainElem inDomainElem;
   const InConstraint inConstraint;
   const InConstraintElem inConstraintElem;
   const PosWeights posWeights;
   const PosWeightsElem posWeightsElem;
   const OnTheBoundary onTheBoundary;
   const EqnOnTheBoundary eqnOnTheBoundary;
   const TestIntegral testIntegral;
   const constraints_init constr_init;
   const constraints_realloc constr_realloc;
   const get_constraints get_constr;
   const constraints_free constr_free;
};


#ifdef __cplusplus
}
#endif

#endif

