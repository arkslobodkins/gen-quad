#ifndef GENERAL_QUADRATURE_H
#define GENERAL_QUADRATURE_H

#include "Matrix.h"
#include "Vector.h"
#include "Tensor.h"
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
#define __attribute__unused __attribute__ ((unused))
#elif defined __INTEL_COMPILER
#define __attribute__unused #pragma unused
#else
#define __attribute__unused
#endif

#define GQ_TRUE 1
#define GQ_FALSE 0
typedef int GQ_BOOL;

#define ONE 1
#define TWO 2
#define THREE 3
#define QUAD_HUGE 10.0
#define INT_8 int_fast8_t
#define CLOSE_TO_ZERO POW_DOUBLE(10, -18)
#define QUAD_TOL  POW_DOUBLE(10.0, -14)
#define BOUND_TOL POW_DOUBLE(10.0, -12)
#define BOUND_CORRECTION POW_DOUBLE(10.0, -14)
#define POW_DOUBLE(x,y) pow(x,y)
#define POW_INT(x,y)    IntPower(x,y)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define is_greater_than_zero(x) ((x) > 0 ? true : false)
#define is_even(x) (x)%2 == 0 ? true : false
#define is_odd(x) (x)%2 == 1 ? true : false
#define SQUARE(x) ((x)*(x))
#define SQRT(x) sqrt((double)(x))
#define SIZE_INT(x) ((x)*sizeof(int))
#define SIZE_DOUBLE(x) ((x)*sizeof(double))

typedef enum { SVD, ELIM } ElimType;
typedef enum { ON, OFF } bool_enum;
typedef enum { orthogonal, monomial } BASIS_TYPE;
typedef enum { INTERVAL, CUBE, SIMPLEX, CUBESIMPLEX, SIMPLEXSIMPLEX } DOMAIN_TYPE;

#define StaticVectorInit(vecName, __len)  \
   Vector vecName;                        \
   double __ ##vecName[__len];            \
   vecName.len = __len;                   \
   vecName.id = __##vecName;

char *get_domain_string(DOMAIN_TYPE D);
char *ElimToString(ElimType elimType);
bool string_to_domain(const char *shape, DOMAIN_TYPE *D);
int IntPower(int x, int power);
int factorial(int n);
int binomial(int k, int n);
long double expIntegral1D(long double c);
double expNDim(int dim, double x[]);
double expIntegralNDimCube(int dim);
double expIntegralNDimSimplex(int dim);
#ifdef _OPENMP
GQ_BOOL OMP_CONDITION(int deg, int dim);
GQ_BOOL PLASMA_CONDITION();
#endif

typedef struct
{
   int nodes_tot;
   int success_node;
   int success_its;
   ElimType elim_type;
} hist_data;

typedef struct
{
   hist_data *hist_array;
   int dim;
   int degree;
   int nodes_initial;
   int nodes_final;
   int num_funcs;
   double res;
   DOMAIN_TYPE D;
   int total_elims;
} history;


#define GQ_SUCCESS      0
#define NULL_VAL       -1
#define ALLOC_FAIL     -2
#define INV_INPUT      -3
#define INF_VAL        -4
#define NAN_VAL        -5
#define QUAD_HUGE_ERR  -6
#define LAPACK_ERR     -7
#define PLASMA_ERR     -8
#define DIV_BY_ZERO    -9
#define NOT_CONVERGE   -10
#define DIVERGE_ERR    -11
#define LARGE_RES      -12
#define CONSTR_ERROR   -13

#define STR_GQ_SUCCESS "GQ_SUCCESS"
#define STR_NULL_VAL "NULL_VAL"
#define STR_ALLOC_FAIL "ALLOC_FAIL"
#define STR_INV_INPUT "INV_INPUT"
#define STR_INF_VAL "INF_VAL"
#define STR_NAN_VAL "NAN_VAL"
#define STR_QUAD_HUGE_ERR "QUAD_HUGE_ERR"
#define STR_LAPACK_ERR "LAPACK_ERR"
#define STR_PLASMA_ERR "PLASMA_ERR"
#define STR_DIV_BY_ZERO "DIV_BY_ZERO"
#define STR_NOT_CONVERGE "NOT_CONVERGE"
#define STR_DIVERGE_ERR "DIVERGE_ERR"
#define STR_LARGE_RES "LARGE_RES"
#define STR_CONSTR_ERROR "CONSTR_ERROR"

#define STR_QUAD_NOT_FULL_INIT "supplied quadrature object was not fully initialized"


#ifdef __cplusplus
}
#endif

#endif

