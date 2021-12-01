#include "GENERAL_QUADRATURE.h"

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

typedef Basis*(*BasisInitPtr)(void *);
typedef void(*IndicesPtr)(int deg, int dim);
typedef void(*BasisFuncsPtr)(Basis *basis, double *x);
typedef void(*BasisDerPtr)(Basis *basis, double *x);
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

struct Basis
{
   void *params;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
   const BasisInterface *interface;
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

BasisInterface SetCubeBasisInterface();
BasisInterface SetSimplexBasisInterface();
BasisInterface SetCubeSimplexBasisInterface();
BasisInterface SetSimplexSimplexBasisInterface();

Basis* BasisInit(void *params, BasisInterface *interface);
void BasisFuncs(Basis *basis, double *x);
void BasisDer(Basis *basis, double *x);
void Integrals(Basis *basis);
void BasisFree(Basis *basis);

struct CubeBasis
{
   CubeParams *params;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
   const BasisInterface *interface;
};

struct SimplexBasis
{
   SimplexParams *params;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
   const BasisInterface *interface;
};

struct CubeSimplexBasis
{
   CubeSimplexParams *params;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
   const BasisInterface *interface;
};

struct SimplexSimplexBasis
{
   SimplexSimplexParams *params;
   int deg;
   int dim;
   int numFuncs;
   INT_8 *basis_indices;
   Vector basis_funcs;
   Vector basis_der;
   Vector basis_integrals;
   const BasisInterface *interface;
};

void ComputeBasisIndices(Basis *basis);

CubeBasis* CubeBasisInit(CubeParams *cubeParams);
void ComputeCubeBasisFuncs(CubeBasis *cubeBasis, double *x);
void ComputeCubeBasisDer(CubeBasis *cubeBasis, double *x);
void CubeBasisIntegrals(CubeBasis *cubeBasis);
void CubeBasisFree(CubeBasis *cubeBasis);

SimplexBasis* SimplexBasisInit(SimplexParams *simplexParams);
void SimplexBasisFuncs(SimplexBasis *simplexBasis, double *x);
void SimplexBasisDer(SimplexBasis *simplexBasis, double *x);
void SimplexBasisIntegrals(SimplexBasis *simplexBasis);
void SimplexBasisFree(SimplexBasis *simplexBasis);

CubeSimplexBasis* CubeSimplexBasisInit(CubeSimplexParams *csParams);
void CubeSimplexBasisFuncs(CubeSimplexBasis *csBasis, double *x);
void CubeSimplexBasisDer(CubeSimplexBasis *csBasis, double *x);
void CubeSimplexBasisIntegrals(CubeSimplexBasis *csBasis);
void CubeSimplexBasisFree(CubeSimplexBasis *csBasis);

SimplexSimplexBasis* SimplexSimplexBasisInit(SimplexSimplexParams *ssParams);
void SimplexSimplexBasisFuncs(SimplexSimplexBasis *ssBasis, double *x);
void SimplexSimplexBasisDer(SimplexSimplexBasis *ssBasis, double *x);
void SimplexSimplexBasisIntegrals(SimplexSimplexBasis *ssBasis);
void SimplexSimplexBasisFree(SimplexSimplexBasis *ssBasis);

void TestBasis();
