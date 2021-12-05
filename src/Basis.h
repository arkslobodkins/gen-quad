#ifndef BASIS_H
#define BASIS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct BasisInterface BasisInterface;
typedef struct Basis Basis;
typedef struct CubeParams CubeParams;
typedef struct CubeBasis CubeBasis;
typedef struct SimplexParams SimplexParams;
typedef struct SimplexBasis SimplexBasis;
typedef struct CubeSimplexParams CubeSimplexParams;
typedef struct CubeSimplexBasis CubeSimplexBasis;
typedef struct SimplexSimplexParams SimplexSimplexParams;
typedef struct SimplexSimplexBasis SimplexSimplexBasis;
typedef struct MixedPolytopeBasis MixedPolytopeBasis;
typedef struct MixedParams MixedParams;


struct Basis
{
   BasisInterface *interface;
   void *params;
   void *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
};

typedef Basis*(*BasisInitPtr)(void *);
typedef void(*IndicesPtr)(int deg, int dim);
typedef void(*BasisFuncsPtr)(Basis *basis, const double *x, Vector v);
typedef void(*BasisDerPtr)(Basis *basis, const double *x);
typedef void(*BasisIntegralsPtr)(Basis *basis);
typedef void(*BasisFreePtr)(Basis *basis);

struct BasisInterface
{
   BasisInitPtr basisInit;
   BasisFuncsPtr computeFuncs;
   BasisDerPtr computeDer;
   BasisIntegralsPtr computeIntegrals;
   BasisFreePtr basisFree;
};

struct CubeParams
{
   int deg;
   int dim;
};

struct SimplexParams
{
   int deg;
   int dim;
};

struct CubeSimplexParams
{
   int deg;
   int dims[2];
};

struct SimplexSimplexParams
{
   int deg;
   int dims[2];
};

struct MixedParams
{
   int deg;
   int dims[2];
};

BasisInterface SetCubeBasisInterface();
BasisInterface SetSimplexBasisInterface();
BasisInterface SetCubeSimplexBasisInterface();
BasisInterface SetSimplexSimplexBasisInterface();

Basis* BasisInit(void *params, BasisInterface *interface);
void BasisFuncs(Basis *basis, const double *x, Vector v);
void BasisDer(Basis *basis, const double *x);
void Integrals(Basis *basis);
void BasisFree(Basis *basis);

struct CubeBasis
{
   const BasisInterface *interface;
   CubeParams *params;
   double *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
};

typedef struct
{
   Vector phi_backw1;
   Vector phi_backw2;
   Vector phi_forw1;
   Vector phi_forw2;
} AddDataSimplex;

struct SimplexBasis
{
   const BasisInterface *interface;
   SimplexParams *params;
   AddDataSimplex *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
};

typedef struct
{
   Vector basis_polytopic;
   Vector phi_backw1;
   Vector phi_backw2;
   Vector phi_forw1;
   Vector phi_forw2;
} AddDataCubeSimplex;

struct CubeSimplexBasis
{
   const BasisInterface *interface;
   CubeSimplexParams *params;
   AddDataCubeSimplex *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
};

typedef struct
{
   Vector basis_polytopic;
   Vector phi_backw1;
   Vector phi_backw2;
   Vector phi_forw1;
   Vector phi_forw2;
} AddDataSimplexSimplex;

struct SimplexSimplexBasis
{
   const BasisInterface *interface;
   SimplexSimplexParams *params;
   AddDataSimplexSimplex *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
};

struct MixedPolytopeBasis
{
   const BasisInterface *interface;
   MixedParams *params;
   void *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
};

void ComputeBasisIndices(Basis *basis);

CubeBasis* CubeBasisInit(CubeParams *cubeParams);
void ComputeCubeBasisFuncs(CubeBasis *cubeBasis, const double *x, Vector v);
void ComputeCubeBasisDer(CubeBasis *cubeBasis, const double *x);
void CubeBasisIntegrals(CubeBasis *cubeBasis);
void CubeBasisFree(CubeBasis *cubeBasis);

SimplexBasis* SimplexBasisInit(SimplexParams *simplexParams);
void SimplexBasisFuncs(SimplexBasis *simplexBasis, const double *x, Vector v);
void SimplexBasisDer(SimplexBasis *simplexBasis, const double *x);
void SimplexBasisIntegrals(SimplexBasis *simplexBasis);
void SimplexBasisFree(SimplexBasis *simplexBasis);

CubeSimplexBasis* CubeSimplexBasisInit(CubeSimplexParams *csParams);
void CubeSimplexBasisFuncs(CubeSimplexBasis *csBasis, const double *x, Vector v);
void CubeSimplexBasisDer(CubeSimplexBasis *csBasis, const double *x);
void CubeSimplexBasisIntegrals(CubeSimplexBasis *csBasis);
void CubeSimplexBasisFree(CubeSimplexBasis *csBasis);

SimplexSimplexBasis* SimplexSimplexBasisInit(SimplexSimplexParams *ssParams);
void SimplexSimplexBasisFuncs(SimplexSimplexBasis *ssBasis, const double *x, Vector v);
void SimplexSimplexBasisDer(SimplexSimplexBasis *ssBasis, const double *x);
void SimplexSimplexBasisIntegrals(SimplexSimplexBasis *ssBasis);
void SimplexSimplexBasisFree(SimplexSimplexBasis *ssBasis);

void SimplexFuncsPolytopicOne(MixedPolytopeBasis *mixedBasis, const double *x, Vector v);
void SimplexFuncsPolytopicTwo(MixedPolytopeBasis *mixedBasis, const double *x, Vector v);

void TestBasis();

#ifdef __cplusplus
}
#endif

#endif
