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

typedef struct
{
   int size[3];
   double ***id;
} Table3d;

struct Basis
{
   BasisInterface *interface;
   void *params;
   void *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *indices;
   Vector functions;
   Vector derivatives;
   Vector integrals;
   Table3d table;
};

typedef Basis*(*BasisInitPtr)(void *);
typedef Basis*(*MakeBasisCopy)(void *);
typedef void(*IndicesPtr)(int deg, int dim);
typedef void(*BasisFuncsPtr)(Basis *basis, const double *x, Vector v);
typedef void(*BasisDerPtr)(Basis *basis, const double *x, Vector v);
typedef void(*BasisIntegralsPtr)(Basis *basis, Vector v);
typedef void(*BasisIntegralsMonomialPtr)(Basis *basis, Vector v);
typedef void(*BasisFreePtr)(Basis *basis);

struct BasisInterface
{
   BasisInitPtr basisInit;
   MakeBasisCopy makeBasisCopy;
   BasisFuncsPtr computeFuncs;
   BasisDerPtr computeDer;
   BasisIntegralsPtr computeIntegrals;
   BasisIntegralsMonomialPtr computeIntegralsMonomial;
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
Basis* BasisCopy(Basis *basis);
void BasisFuncs(Basis *basis, const double *x, Vector v);
void BasisDer(Basis *basis, const double *x, Vector v);
void BasisIntegrals(Basis *basis, Vector v);
void BasisIntegralsMonomial(Basis *basis, Vector v);
void BasisMonomial(Basis *basis, const double *x, Vector phi);
void BasisFree(Basis *basis);

typedef struct
{
   INT_8 *idMap;
} AddDataCube;

struct CubeBasis
{
   BasisInterface *interface;
   CubeParams *params;
   AddDataCube *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *indices;
   Vector functions;
   Vector derivatives;
   Vector integrals;
   Table3d table;
};

typedef struct
{
   Vector phi_backw1;
   Vector phi_forw1;
   int* xPower;
   RMatrix xFactor;
   INT_8 *idMap;
} AddDataSimplex;

struct SimplexBasis
{
   BasisInterface *interface;
   SimplexParams *params;
   AddDataSimplex *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *indices;
   Vector functions;
   Vector derivatives;
   Vector integrals;
   Table3d table;
};

typedef struct
{
   Vector basis_polytopic;
   Vector phi_backw1;
   Vector phi_forw1;
   int* xPower[2];
   RMatrix xFactor[2];
   INT_8 *idMap;
} AddDataCubeSimplex;

struct CubeSimplexBasis
{
   BasisInterface *interface;
   CubeSimplexParams *params;
   AddDataCubeSimplex *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *indices;
   Vector functions;
   Vector derivatives;
   Vector integrals;
   Table3d table;
};

typedef struct
{
   Vector basis_polytopic;
   Vector phi_backw1;
   Vector phi_forw1;
   int* xPower[2];
   RMatrix xFactor[2];
   INT_8 *idMap;
} AddDataSimplexSimplex;

struct SimplexSimplexBasis
{
   BasisInterface *interface;
   SimplexSimplexParams *params;
   AddDataSimplexSimplex *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *indices;
   Vector functions;
   Vector derivatives;
   Vector integrals;
   Table3d table;
};

typedef struct
{
   Vector basis_polytopic;
   Vector phi_backw1;
   Vector phi_forw1;
   int* xPower[2];
   RMatrix xFactor[2];
   INT_8 *idMap;
} MixedPolyTopeData;

struct MixedPolytopeBasis
{
   BasisInterface *interface;
   MixedParams *params;
   MixedPolyTopeData *addData;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *indices;
   Vector functions;
   Vector derivatives;
   Vector integrals;
   Table3d table;
};

CubeBasis* CubeBasisInit(CubeParams *params);
CubeBasis* MakeCubeCopy(CubeBasis *basis);
void ComputeCubeBasisFuncs(CubeBasis *basis, const double *x, Vector v);
void ComputeCubeBasisDer(CubeBasis *basis, const double *x, Vector v);
void CubeBasisIntegrals(CubeBasis *basis, Vector v);
void CubeBasisIntegralsMonomial(CubeBasis *basis, Vector v);
void CubeBasisFree(CubeBasis *basis);

SimplexBasis* SimplexBasisInit(SimplexParams *params);
SimplexBasis* MakeSimplexCopy(SimplexBasis *basis);
void SimplexBasisFuncs(SimplexBasis *basis, const double *x, Vector v);
void SimplexBasisDer(SimplexBasis *basis, const double *x, Vector v);
void SimplexBasisIntegrals(SimplexBasis *basis, Vector v);
void SimplexBasisIntegralsMonomial(SimplexBasis *basis, Vector v);
void SimplexBasisFree(SimplexBasis *basis);

CubeSimplexBasis* CubeSimplexBasisInit(CubeSimplexParams *params);
CubeSimplexBasis* MakeCubeSimplexCopy(CubeSimplexBasis *basis);
void CubeSimplexBasisFuncs(CubeSimplexBasis *basis, const double *x, Vector v);
void CubeSimplexBasisDer(CubeSimplexBasis *basis, const double *x, Vector v);
void CubeSimplexBasisIntegrals(CubeSimplexBasis *basis, Vector v);
void CubeSimplexBasisIntegralsMonomial(CubeSimplexBasis *basis, Vector v);
void CubeSimplexBasisFree(CubeSimplexBasis *basis);

SimplexSimplexBasis* SimplexSimplexBasisInit(SimplexSimplexParams *params);
SimplexSimplexBasis* MakeSimplexSimplexCopy(SimplexSimplexBasis *basis);
void SimplexSimplexBasisFuncs(SimplexSimplexBasis *basis, const double *x, Vector v);
void SimplexSimplexBasisDer(SimplexSimplexBasis *basis, const double *x, Vector v);
void SimplexSimplexBasisIntegrals(SimplexSimplexBasis *basis, Vector v);
void SimplexSimplexBasisIntegralsMonomial(SimplexSimplexBasis *basis, Vector v);
void SimplexSimplexBasisFree(SimplexSimplexBasis *basis);

void SimplexFuncsPolytopicOne(MixedPolytopeBasis *basis, const double *x, Vector v);
void SimplexFuncsPolytopicTwo(MixedPolytopeBasis *basis, const double *x, Vector v);

double orthogonal_simplex_basis_test(int deg, int dim);
double orthogonal_cube_basis_test(int deg, int dim);
void PrintBasisIndices(Basis *basis);

#ifdef __cplusplus
}
#endif

#endif
