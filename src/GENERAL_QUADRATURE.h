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


#define ij2(i, j, numCols) i*numCols+j // maps indices from 2-d layout to memory
#define POW(x,y) pow(x,y)
#define MAX(a, b) a > b ? a : b
#define MIN(a, b) a < b ? a : b
#define SQUARE(x) ((x)*(x))
#define SQRT(x) sqrt(x)
#define QUAD_TOL POW(10, -15)
#define ONE 1
#define TWO 2
#define THREE 3
#define size_int sizeof(int)
#define size_int sizeof(int)
#define size_double sizeof(double)
#define SIZE_INT(x) (x)*sizeof(int)
#define SIZE_DOUBLE(x) (x)*sizeof(double)
#define size_quadrature sizeof(quadrature)
#define size_quadParams sizeof(quadParams)


typedef enum { ON, OFF } BOOLEAN;
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
   int node_id;
   int eqn_id;
   double t_min;
   int flag;
} constr_node_data;


typedef struct
{
   BOOLEAN ACTIVE_CONSTRAINTS;
   NODE_OR_WEIGHT N_OR_W;
   int boundary_node_id;
   int  eqn_id;
   Vector additional_constr_eqn;
} constr_vect_data;


typedef struct quadrature quadrature;
typedef struct const_quadrature const_quadrature;

typedef struct quadParams quadParams;
typedef void(*EvalBasis)(const int_fast8_t *basis, const double *x, const quadParams *params, double *phi);
typedef void(*EvalBasisDer)(const int_fast8_t *basis, const double *x, const quadParams *params, double *phiPrime);
typedef void(*BasisIntegrals)(const quadParams *params, double *integrals);
typedef bool(*InDomain)(const_quadrature *quad);
typedef bool(*InDomainElem)(const_quadrature *quad, int elem);
typedef bool(*InDomainPoint)(double *x, int dim);

typedef struct constraints constraints;
typedef constraints*(*constraints_init)(int *dims);
typedef void (*constraints_realloc)(constraints *cons, int *dims);
typedef void (*get_constraints)(constraints *cons);
typedef void (*constraints_free)(constraints *cons);
typedef void (*_SetParams)(int dim, int num_dims, int *dims, int deg, quadParams *params);
//typedef bool (*InConstraint)(const InDomain inDomain, const_quadrature *quad);
//typedef bool (*InConstraintElem)(const InDomainElem InDomainElem, const_quadrature *quad, int elem);

struct quadParams
{
   int dim;
   int num_dims;
   int *dims;
   int deg;
   int num_funs;
};

struct constraints
{
   int dim;
   int *dims;
   Matrix M;
   Vector b;
};

typedef struct
{
   EvalBasis evalBasis;
   EvalBasisDer evalBasisDer;
   BasisIntegrals basisIntegrals;
   InDomain inDomain;
   InDomainElem inDomainElem;
//   InDomainPoint inDomainPoint;

   constraints_init constr_init;
   constraints_realloc constr_realloc;
   get_constraints get_constr;
   constraints_free constr_free;
} _DomainFuncs;

typedef _DomainFuncs *DomainFuncs;
typedef void(*SetDomain)(quadrature *q);


typedef struct
{
   constraints_init constr_init;
   get_constraints get_constr;
   constraints_free constr_free;
} ConstraintFuncs;
typedef ConstraintFuncs *constraintFuncs;
typedef void(*SetConstraint)(ConstraintFuncs *constrFuncs);


struct quadrature
{
   DOMAIN_TYPE D;
   int k;
   double *w;
   double *x;
   double *z;
   quadParams *params;
   constraints *cons;

   _SetParams setParams; // Function that sets up params structure for the appropriate domain
   SetDomain setDomain;  // Function that sets up all domain functions to point to appropriate domain

   DomainFuncs domFuncs; // structure that points to all the functions below
   EvalBasis evalBasis;
   EvalBasisDer evalBasisDer;
   BasisIntegrals basisIntegrals;
   InDomain inDomain;
   //InDomainElem inDomainElem;
   //InDomainPoint inDomainPoint;
   constraints_init constr_init;
   constraints_realloc constr_realloc;
   get_constraints get_constr;
   constraints_free constr_free;
};


struct const_quadrature
{
   const DOMAIN_TYPE D;
   const int k;
   const double *w;
   const double *x;
   const double *z;
   const quadParams *params;
   const constraints *cons;

   const _SetParams setParams; // Function that sets up params structure for the appropriate domain
   const SetDomain setDomain;  // Function that sets up all domain functions to point to appropriate domain

   const DomainFuncs domFuncs; // structure that points to all the functions below
   const EvalBasis evalBasis;
   const EvalBasisDer evalBasisDer;
   const BasisIntegrals basisIntegrals;
   const InDomain inDomain;
   //const InDomainElem inDomainElem;
   //const InDomainPoint inDomainPoint;
   const constraints_init constr_init;
   const constraints_realloc constr_realloc;
   const get_constraints get_constr;
   const constraints_free constr_free;
};


#ifdef __cplusplus
}
#endif

#endif

