#ifndef POLYTOPE_H
#define POLYTOPE_H


#include "Constraints.h"


#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
   int_fast8_t *basis_indices;
   double *basis_functions;
   double *basis_integrals;
} Basis;

typedef struct Polytope Polytope;

typedef Polytope*(*_PolytopeInit)(void *);
typedef void (*_PolytopeFree)(void *);
typedef int (*_dim)(void *);
typedef int (*_deg)(void *);
typedef DOMAIN_TYPE(*_type)(void *);
typedef const char*(*_string)();
typedef double(*_get_area)(Polytope *);


typedef struct
{
  _PolytopeInit polytope_init;
  _PolytopeFree polytope_free;
  _dim dim;
  _deg deg;
  _type type;
  _string string;
  _get_area get_area;

} PolytopeInterface;

struct Polytope
{
   Basis *basis;
   void *constr;
   void *domain_params;

   int basisSize;

   const PolytopeInterface *interface;
};

Polytope * PolytopeInit(void *init_params, PolytopeInterface *interface);
void PolytopeFree(Polytope* polytope);

int PolytopeDim(Polytope *polytope);
int PolytopeDeg(Polytope *polytope);
DOMAIN_TYPE PolytopeType(Polytope *polytope);
const char *PolytopeString(Polytope *Polytope);
double PolytopeGetArea(Polytope *polytope);
void ComputeBasis(Polytope * polytope);

typedef struct
{
  int dim;
  int deg;
  DOMAIN_TYPE D;
} CubeParams;

typedef struct
{
  int dim;
  int deg;
} InitialCubeParams;

typedef struct
{
   int dim;
   RMatrix M;
   Vector b;
   RMatrix M_FULL;
   Vector b_FULL;
} CubeConstr;

typedef struct
{
   Basis *basis;
   constraints *constr;
   CubeParams *params;

   int basisSize;
   const PolytopeInterface *interface;
} Cube;


PolytopeInterface SetCubeInterface();
Cube * CubeInit(InitialCubeParams *params, PolytopeInterface *interface);
constraints * CubeConstrInit(Cube *cube);

int CubeDim(Cube *cube);
int CubeDeg(Cube *cube);
DOMAIN_TYPE CubeType(Cube *cube);
const char *GetCubeString();
double CubeGetArea(Cube *cube);

void TestPolytope();

#ifdef __cplusplus
}
#endif


#endif
