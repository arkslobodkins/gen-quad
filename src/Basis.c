/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "Quadrature.h"
#include "Gauss_Lib/Jacobi.h"
#include "GeneralGaussTensor.h"
#include "AddDimension.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#pragma omp declare simd
static inline int id2(int i, int j, int N) { return i*N+j; }

static int BasisSize(int deg, int dim);
static void BasisIndices(int deg, int dim, INT_8 *F);

static Table3d Table3dCreate(int deg, int dim);
static void Table3dFree(Table3d table);
static void ComputeTable(Table3d table);
static void CopyTable(Table3d t1, Table3d t2);

static void LegendrePoly(int order, double x, double *p);
static void LegendrePolyAccurate(int order, double x, double *p);
static void LegendrePolyAndPrime(int order, double x, double *p, double *dp);
static void LegendrePolyAndPrimeAccurate(int order, double x, double *p, double *dp);

static void JacobiPoly(int order, double x, int alpha, double *p);
static void JacobiPolyWithTable(int order, double x, int alpha, double *p, Table3d table);
static void JacobiPolyWithTableAccurate(int order, double x, int alpha, double *p, Table3d table);

static void IntegralsCubePolyhedralMonomial(MixedPolytopeBasis *basis, Vector v);
static void IntegralsSimplexPolyhedralMonomialOne(MixedPolytopeBasis *basis, Vector v);
static void IntegralsSimplexPolyhedralMonomialTwo(MixedPolytopeBasis *basis, Vector v);


Basis* BasisInit(void *params, BasisInterface *interface)
{
   Basis *basis = interface->basisInit(params);
   basis->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(basis->interface) = *interface;

   int dim = basis->dim;
   int numFuncs = basis->numFuncs;

   basis->functions   = Vector_init(numFuncs);
   basis->derivatives = Vector_init(numFuncs*dim);
   basis->integrals   = Vector_init(numFuncs);
   return basis;
}

Basis* BasisCopy(Basis *basis)
{
   return basis->interface->makeBasisCopy(basis);
}

void BasisFuncs(Basis *basis, const double *x, Vector F)
{
   basis->interface->computeFuncs(basis, x, F);
}

void BasisDer(Basis *basis, const double *x, Vector dF)
{
   basis->interface->computeDer(basis, x, dF);
}

void BasisIntegrals(Basis *basis, Vector I)
{
   basis->interface->computeIntegrals(basis, I);
}

void BasisIntegralsMonomial(Basis *basis, Vector I)
{
   basis->interface->computeIntegralsMonomial(basis, I);
}

void BasisFree(Basis *basis)
{
   if(basis == NULL) return;
   BasisInterface *interface = basis->interface;

   Vector_free(basis->integrals);
   Vector_free(basis->derivatives);
   Vector_free(basis->functions);

   basis->interface->basisFree(basis);
   if(interface)
   {
      free(interface);
      interface = NULL;
   }
}

static int BasisSize(int deg, int dim)
{
   assert(dim >= 0 && deg >= 0);

   unsigned int first = deg + 1;
   unsigned int last = deg + dim;
   unsigned int product = 1;

   unsigned int i = first, counter = 1;
   for(; i <= last; ++i, ++counter)
      product = product * i/counter;

   return product;
}

static void BasisIndices(int deg, int dim, INT_8 *F)
{
   if(dim == 1) {
      for(int i = 0; i <= deg; ++i)
         F[i] = i;
      return;
   }

   int counter;
   // compute basis indices using nested recursion if dimension >= 2
   for(int dout = 2; counter = 0, dout <= dim; ++dout)
   {
      for(int i = 0; i <= deg; ++i)
      {
         int rsize = BasisSize(deg-i, dout-1);
         INT_8 recursiveF[rsize*(dout-1)];
         BasisIndices(deg-i, dout-1, recursiveF);
         if(dout == dim)
         {
            for(int k = 0; k < rsize; ++k)
            {
               for(int d = 0; d < dim-1; ++d)
                  F[counter+k*dim+d] = recursiveF[k*(dim-1)+d];
               F[counter+k*dim+dout-1] = i;
            }
            counter += dim*rsize;
         }
      }
   }
}

#pragma omp declare simd
static inline double DoubleIntPower(double x, int power)
{
    double result = 1.0;
    for (;;) {
        if (power & 1)
            result *= x;
        power >>= 1;
        if (!power)
            break;
        x *= x;
    }
    return result;
}

void BasisMonomial(Basis *basis, const double *x, Vector F)
{
   int dim      = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *ind   = basis->indices;

   VSetToOne(F);
   for(int k = 0; k < numFuncs; ++k) {
      for(int d = 0; d < dim; ++d) {
         INT_8 basis_power = ind[id2(k, d, dim)];
         F.id[k] *= DoubleIntPower(x[d], basis_power);
      }
   }
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetCubeBasisInterface()
{
   BasisInterface cubeInterface;
   cubeInterface.basisInit        = (BasisInitPtr)&CubeBasisInit;
   cubeInterface.makeBasisCopy    = (MakeBasisCopy)&MakeCubeCopy;
   cubeInterface.computeFuncs     = (BasisFuncsPtr)&ComputeCubeBasisFuncs;
   cubeInterface.computeDer       = (BasisDerPtr)&ComputeCubeBasisDer;
   cubeInterface.computeIntegrals = (BasisIntegralsPtr)&CubeBasisIntegrals;
   cubeInterface.computeIntegralsMonomial = (BasisIntegralsMonomialPtr)&CubeBasisIntegralsMonomial;
   cubeInterface.basisFree        = (BasisFreePtr)&CubeBasisFree;
   return cubeInterface;
}


CubeBasis* CubeBasisInit(CubeParams *cubeParams)
{
   int deg = cubeParams->deg;
   int dim = cubeParams->dim;
   assert(dim >= 2 && deg >= 1);

   CubeBasis *basis   = (CubeBasis *)malloc(sizeof(CubeBasis));
   basis->params      = (CubeParams *)malloc(sizeof(CubeParams));
   basis->params->deg = deg;
   basis->params->dim = dim;
   basis->deg         = deg;
   basis->dim         = dim;
   basis->numFuncs    = BasisSize(deg, dim);
   basis->addData     = NULL;
   memset(&basis->table, 0, sizeof(Table3d));

   int numFuncs = basis->numFuncs;
   basis->indices = (INT_8 *)malloc(numFuncs*dim*sizeof(INT_8));
   BasisIndices(basis->deg, basis->dim, basis->indices);

   basis->addData        = (AddDataCube *)malloc(sizeof(AddDataCube));
   basis->addData->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   INT_8 *idMap          = basis->addData->idMap;
   INT_8 *idBasis        = basis->indices;
   for(int i = 0; i < numFuncs; ++i)
      #pragma omp simd
      for(int j = 0; j < dim; ++j)
         idMap[(id2(j, i, numFuncs))] = idBasis[id2(i, j, dim)];
   return basis;
}


CubeBasis* MakeCubeCopy(CubeBasis *basis)
{
   CubeBasis *basisCopy   = (CubeBasis *)malloc(sizeof(CubeBasis));
   basisCopy->params      = (CubeParams *)malloc(sizeof(CubeParams));
   basisCopy->params->deg = basis->deg;
   basisCopy->params->dim = basis->dim;
   basisCopy->deg         = basis->deg;
   basisCopy->dim         = basis->dim;
   basisCopy->numFuncs    = basis->numFuncs;
   basisCopy->addData     = NULL;

   BasisInterface interface = SetCubeBasisInterface();
   basisCopy->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(basisCopy->interface) = interface;

   int dim = basisCopy->dim;
   int numFuncs = basisCopy->numFuncs;

   basisCopy->indices = (INT_8 *)malloc(numFuncs*dim*sizeof(INT_8));
   memcpy(basisCopy->indices, basis->indices, numFuncs*dim*sizeof(INT_8));

   basisCopy->addData = (AddDataCube *)malloc(sizeof(AddDataCube));
   basisCopy->addData->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   memcpy(basisCopy->addData->idMap, basis->addData->idMap, numFuncs*dim*sizeof(INT_8));

   basisCopy->functions   = Vector_init(numFuncs);
   basisCopy->derivatives = Vector_init(numFuncs*dim);
   basisCopy->integrals   = Vector_init(numFuncs);

   return basisCopy;
}


void ComputeCubeBasisFuncs(CubeBasis *basis, const double *x, Vector F)
{
   int deg      = basis->deg;
   int dim      = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *idMap = basis->addData->idMap;

   VSetToOne(F);
   for(int d = 0; d < dim; ++d)
   {
      double legendre[deg+1];
      LegendrePolyAccurate(deg+1, 2.*x[d]-1., legendre);
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         F.id[k] *= legendre[idMap[d*numFuncs+k]];
   }
}


void ComputeCubeBasisDer(CubeBasis *basis, const double *x, Vector dF)
{
   int deg = basis->deg;
   int dim = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *idMap = basis->addData->idMap;
   Tensor2D dF2D = VectorToTensor2D(dim, numFuncs, dF);

   double legendre[dim][deg+1];
   double dxlegendre[dim][deg+1];
   for(int d = 0; d < dim; ++d)
      LegendrePolyAndPrimeAccurate(deg+1, 2.*x[d]-1., legendre[d], dxlegendre[d]);

   VSetToOne(dF);
   for(int d = 0; d < dim; ++d)
   {
      // dimension < d
      for(int j = 0; j < d; ++j)
         #pragma omp simd
         for(int k = 0; k < numFuncs; k++)
            TID2(dF2D, d, k) *= legendre[j][idMap[j*numFuncs+k]];

      // dimension = d
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         TID2(dF2D, d, k) *= 2.0 * dxlegendre[d][idMap[d*numFuncs+k]];

      // dimension > d
      for(int j = d+1; j < dim; ++j)
         #pragma omp simd
         for(int k = 0; k < numFuncs; k++)
            TID2(dF2D, d, k) *= legendre[j][idMap[j*numFuncs+k]];
   }
}


void CubeBasisIntegrals(CubeBasis *basis, Vector I)
{
   assert(basis->numFuncs == I.len);
   VSetToZero(I);
   I.id[0] = 1.0;
}


void CubeBasisIntegralsMonomial(CubeBasis *basis, Vector I)
{
   int dim = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *ind = basis->indices;

   VSetToOne(I);
   for(int i = 0; i < numFuncs; ++i)
      for(int d = 0; d < dim; ++d)
         I.id[i] /= (ind[id2(i, d, dim)] + 1);
}


void CubeBasisFree(CubeBasis *basis)
{
   if(basis == NULL) return;

   if(basis->params) {
      free(basis->params);
      basis->params = NULL;
   }
   if(basis->addData) {
      free(basis->addData->idMap);
      free(basis->addData);
      basis->addData = NULL;
   }
   if(basis->indices) {
      free(basis->indices);
      basis->indices = NULL;
   }
   free(basis);
}


void ComputeCube2dTest(Basis *basis, const double *x, Vector F)
{
   int deg = basis->deg;

   double legendre1[deg+1];
   double legendre2[deg+1];
   LegendrePolyAccurate(deg+1, 2.*x[0]-1., legendre1);
   LegendrePolyAccurate(deg+1, 2.*x[1]-1., legendre2);

   for(int i = 0, count = 0; i < deg; ++i)
      for(int j = 0; j < deg-i; ++j)
         F.id[count++] = legendre1[j] * legendre2[i];
}

void ComputeCube3dTest(Basis *basis, const double *x, Vector F)
{
   int deg = basis->deg;
   double legendre[3][deg+1];
   for(int i = 0; i < 3; ++i)
      LegendrePolyAccurate(deg+1, 2.*x[i]-1., legendre[i]);

   for(int i = 0, count = 0; i < deg; ++i)
      for(int j = 0; j < deg-i; ++j)
         for(int k = 0; k < deg-i-j; ++k)
            F.id[count++] = legendre[0][k] * legendre[1][j] * legendre[2][i];
}


void ComputeCube4dTest(Basis *basis, const double *x, Vector F)
{
   int deg = basis->deg;
   double legendre[4][deg+1];
   for(int i = 0; i < 4; ++i)
      LegendrePolyAccurate(deg+1, 2.*x[i]-1., legendre[i]);

   for(int i = 0, count = 0; i < deg; ++i)
      for(int j = 0; j < deg-i; ++j)
         for(int k = 0; k < deg-i-j; ++k)
            for(int l = 0; l < deg-i-j-k; ++l)
               F.id[count++] = legendre[0][l] * legendre[1][k] * legendre[2][j] * legendre[3][i];
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetSimplexBasisInterface()
{
   BasisInterface simplexInterface;
   simplexInterface.basisInit        = (BasisInitPtr)&SimplexBasisInit;
   simplexInterface.makeBasisCopy    = (MakeBasisCopy)&MakeSimplexCopy;
   simplexInterface.computeFuncs     = (BasisFuncsPtr)&SimplexBasisFuncs;
   simplexInterface.computeDer       = (BasisDerPtr)&SimplexBasisDer;
   simplexInterface.computeIntegrals = (BasisIntegralsPtr)&SimplexBasisIntegrals;
   simplexInterface.basisFree        = (BasisFreePtr)&SimplexBasisFree;
   simplexInterface.computeIntegralsMonomial = (BasisIntegralsMonomialPtr)&SimplexBasisIntegralsMonomial;
   return simplexInterface;
}


SimplexBasis* SimplexBasisInit(SimplexParams *simplexParams)
{
   int deg = simplexParams->deg;
   int dim = simplexParams->dim;
   assert(dim >= 2 && deg >= 1);

   SimplexBasis *basis = (SimplexBasis *)malloc(sizeof(SimplexBasis));
   basis->params       = (SimplexParams *)malloc(sizeof(SimplexParams));
   basis->params->deg  = deg;
   basis->params->dim  = dim;
   basis->deg          = deg;
   basis->dim          = dim;
   basis->numFuncs     = BasisSize(deg, dim);

   basis->indices = (INT_8 *) malloc(basis->numFuncs*dim*sizeof(INT_8));
   BasisIndices(basis->deg, basis->dim, basis->indices);

   basis->addData = (AddDataSimplex *)malloc(sizeof(AddDataSimplex));
   int numFuncs = basis->numFuncs;
   AddDataSimplex *addData = basis->addData;

   addData->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   INT_8 *idMap   = addData->idMap;
   INT_8 *idBasis = basis->indices;
   for(int i = 0; i < numFuncs; ++i)
      #pragma omp simd
      for(int j = 0; j < dim; ++j)
         idMap[id2(j, i, numFuncs)] = idBasis[id2(i, j, dim)];

   addData->phi_backw1 = Vector_init(numFuncs);
   addData->phi_forw1  = Vector_init(numFuncs);

   addData->xPower  = (int *)malloc(dim*numFuncs*sizeof(int));
   addData->xFactor = RMatrix_init(dim, numFuncs);

   int *xPower = addData->xPower;
   memset(xPower, 0, dim*numFuncs*sizeof(int));
   for(int d = 1; d < dim; ++d)
      for(int k = 0; k < numFuncs; ++k)
         for(int i = 0; i < d; ++i)
            xPower[id2(d, k, numFuncs)] += basis->indices[k*dim+i];

   basis->table = Table3dCreate(deg, dim);
   ComputeTable(basis->table);
   return basis;
}


SimplexBasis* MakeSimplexCopy(SimplexBasis *b)
{
   int deg = b->deg;
   int dim = b->dim;
   int numFuncs = b->numFuncs;

   SimplexBasis *bcopy = (SimplexBasis *)malloc(sizeof(SimplexBasis));
   bcopy->params       = (SimplexParams *)malloc(sizeof(SimplexParams));
   bcopy->params->deg  = deg;
   bcopy->params->dim  = dim;
   bcopy->deg          = deg;
   bcopy->dim          = dim;
   bcopy->numFuncs     = numFuncs;

   BasisInterface interface = SetSimplexBasisInterface();
   bcopy->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(bcopy->interface) = interface;

   bcopy->indices = (INT_8 *)malloc(numFuncs*dim*sizeof(INT_8));
   memcpy(bcopy->indices, b->indices, numFuncs*dim*sizeof(INT_8));

   bcopy->addData = (AddDataSimplex *)malloc(sizeof(AddDataSimplex));
   AddDataSimplex *addDataCp = bcopy->addData;
   AddDataSimplex *addData = b->addData;

   addDataCp->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   memcpy(addDataCp->idMap, addData->idMap, numFuncs*dim*sizeof(INT_8));

   addDataCp->phi_backw1 = Vector_init(numFuncs);
   addDataCp->phi_forw1  = Vector_init(numFuncs);
   addDataCp->xFactor = RMatrix_init(dim, numFuncs);

   addDataCp->xPower = (int *)malloc(dim*numFuncs*sizeof(int));
   memcpy(addDataCp->xPower, addData->xPower, dim*numFuncs*sizeof(int));

   bcopy->table = Table3dCreate(deg, dim);
   CopyTable(b->table, bcopy->table);

   bcopy->functions   = Vector_init(numFuncs);
   bcopy->derivatives = Vector_init(numFuncs*dim);
   bcopy->integrals   = Vector_init(numFuncs);

   return bcopy;
}


void SimplexBasisFuncs(SimplexBasis *basis, const double *x, Vector F)
{
   int deg      = basis->deg;
   int dim      = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *idMap = basis->addData->idMap;

   double legendre[(deg+1)];
   double jacobi[dim-1][SQUARE(deg+1)];

   LegendrePolyAccurate(deg+1, (2.0*x[dim-1]-x[dim-2])/x[dim-2], legendre);
   double xCoord[dim-1];
   for(int d = 0; d < dim-2; ++d)
      xCoord[d] = x[dim-d-2]/x[dim-d-3];
   xCoord[dim-2] = x[0];

   double jCoord[dim];
   for(int d = 1; d < dim-1; ++d)
      jCoord[d] = 1.0-2.0*x[dim-d-1]/x[dim-d-2];
   jCoord[dim-1] = 1.0-2.0*x[0];

   for(int d = 1; d < dim; ++d) {
      for(int j = 0; j < deg+1; ++j) {
         int nextAlpha = (deg+1)*j;
         JacobiPolyWithTableAccurate(deg+1, jCoord[d], 2*j+d, &jacobi[d-1][nextAlpha], basis->table);
      }
   }

   RMatrix xFactor = basis->addData->xFactor;
   int *xPower = basis->addData->xPower;
   for(int d = 1; d < dim; ++d)
      for(int k = 0; k < numFuncs; ++k)
            xFactor.rid[d][k] = DoubleIntPower(xCoord[d-1], xPower[id2(d, k, numFuncs)]);

   for(int k = 0; k < numFuncs; ++k)
      F.id[k] = legendre[idMap[k]];
   for(int d = 1; d < dim; ++d)
      for(int k = 0; k < numFuncs; ++k)
         F.id[k] *= jacobi[d-1][idMap[d*numFuncs+k] + (deg+1)*xPower[id2(d, k, numFuncs)]] * xFactor.rid[d][k];
}


void SimplexBasisDer(SimplexBasis *basis, const double *x, Vector dF)
{
   int dim = basis->params->dim;
   int numFuncs = basis->numFuncs;
   double h = POW_DOUBLE(10, -6);

   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;

   for(int d = 0; d < dim; ++d)
   {
      double x_backw1[dim];
      double x_forw1[dim];
      for(int d_t = 0; d_t < dim; ++d_t)
      {
         x_backw1[d_t] = x[d_t];
         x_forw1[d_t] = x[d_t];
      }
      x_backw1[d] = x_backw1[d] - h;
      x_forw1[d] = x_forw1[d] + h;

      SimplexBasisFuncs(basis, x_backw1, phi_backw1);
      SimplexBasisFuncs(basis, x_forw1, phi_forw1);
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         dF.id[d*numFuncs+k] = phi_forw1.id[k] - phi_backw1.id[k];
   }
   VScale(1.0/(2.0*h), dF);
}


void SimplexBasisIntegrals(SimplexBasis *basis, Vector I)
{
   assert(basis->numFuncs == I.len);
   VSetToZero(I);
   I.id[0] = 1.0 / (double)factorial(basis->dim);
}


void SimplexBasisIntegralsMonomial(SimplexBasis *basis, Vector I)
{
   int dim      = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *ind   = basis->indices;

   VSetToOne(I);
   for(int i = 0; i < numFuncs; ++i)
   {
      for(int d = 0; d < dim; ++d)
      {
         double power = 0;
         for(int r = 0; r < dim-d; ++r)
            power += (double)ind[(i+1)*dim-(r+1)];
         I.id[i] /= (power+dim-d);
      }
   }
}


void SimplexBasisFree(SimplexBasis *basis)
{
   if(basis == NULL) return;

   Table3dFree(basis->table);
   if(basis->params) {
      free(basis->params);
      basis->params = NULL;
   }

   if(basis->addData) {
      AddDataSimplex *addData = basis->addData;
      Vector_free(addData->phi_backw1);
      Vector_free(addData->phi_forw1);
      free(addData->idMap);
      RMatrix_free(addData->xFactor);
      free(addData->xPower);
      free(addData);
      basis->addData= NULL;
   }
   if(basis->indices) {
      free(basis->indices);
      basis->indices = NULL;
   }
   free(basis);
}


void ComputeSimplex2dTest(Basis *basis, const double *x, Vector F)
{
   int deg = basis->deg;
   double legendre[deg+1];
   double jacobi[deg+1];

   LegendrePoly(deg+1, (2.*x[1]-x[0]) / x[0], legendre);
   for(int i = 0, count = 0; i < deg; ++i)
   {
      JacobiPoly(deg+1, 1.0-2.0*x[0], 2*i+1, jacobi);
      for(int j = 0; j < deg-i; ++j)
         F.id[count++] = jacobi[j] * legendre[i] * DoubleIntPower(x[0], i);
   }
}


void ComputeSimplex3dTest(Basis *basis, const double *x, Vector F)
{
   int deg = basis->deg;
   double legendre[deg+1];
   double jacobi1[deg+1];
   double jacobi2[deg+1];

   LegendrePoly(deg+1, (2.*x[2]-x[1]) / x[1], legendre);
   for(int i = 0, count = 0; i < deg; ++i)
   {
      JacobiPoly(deg+1, 1.0-2.0*x[1]/x[0], 2*i+1, jacobi1);
      for(int j = 0; j < deg-i; ++j)
      {
         JacobiPoly(deg+1, 1.0-2.0*x[0], 2*(i+j)+2, jacobi2);
         for(int k = 0; k < deg-i-j; ++k)
            F.id[count++] = jacobi2[k] * DoubleIntPower(x[0], i+j) *
                            jacobi1[j] * DoubleIntPower(x[1]/x[0], i) *
                            legendre[i];
      }
   }
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetCubeSimplexBasisInterface()
{
   BasisInterface csInterface;
   csInterface.basisInit        = (BasisInitPtr)&CubeSimplexBasisInit;
   csInterface.makeBasisCopy    = (MakeBasisCopy)&MakeCubeSimplexCopy;
   csInterface.computeFuncs     = (BasisFuncsPtr)&CubeSimplexBasisFuncs;
   csInterface.computeDer       = (BasisDerPtr)&CubeSimplexBasisDer;
   csInterface.computeIntegrals = (BasisIntegralsPtr)&CubeSimplexBasisIntegrals;
   csInterface.basisFree        = (BasisFreePtr)&CubeSimplexBasisFree;
   csInterface.computeIntegralsMonomial = (BasisIntegralsMonomialPtr)&CubeSimplexBasisIntegralsMonomial;
   return csInterface;
}


CubeSimplexBasis* CubeSimplexBasisInit(CubeSimplexParams *csParams)
{
   int deg = csParams->deg;
   int dim1 = csParams->dims[0];
   int dim2 = csParams->dims[1];
   int dim = dim1+dim2;
   assert(dim1 >= 1 && dim2 >= 2 && deg >= 1);

   CubeSimplexBasis *basis = (CubeSimplexBasis *)malloc(sizeof(CubeSimplexBasis));
   basis->params = (CubeSimplexParams *)malloc(sizeof(CubeSimplexParams));

   basis->params->deg     = deg;
   basis->params->dims[0] = dim1;
   basis->params->dims[1] = dim2;
   basis->deg             = deg;
   basis->dim             = dim;
   basis->numFuncs        = BasisSize(deg, dim);

   basis->indices = (INT_8 *) malloc(basis->numFuncs*dim*sizeof(INT_8));
   BasisIndices(basis->deg, basis->dim, basis->indices);

   basis->addData              = (AddDataCubeSimplex *)malloc(sizeof(AddDataCubeSimplex));
   AddDataCubeSimplex *addData = basis->addData;
   int numFuncs                = basis->numFuncs;
   addData->phi_backw1         = Vector_init(numFuncs);
   addData->phi_forw1          = Vector_init(numFuncs);
   addData->basis_polytopic    = Vector_init(numFuncs);

   addData->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   INT_8 *idMap   = addData->idMap;
   INT_8 *idBasis = basis->indices;
   for(int i = 0; i < numFuncs; ++i)
      #pragma omp simd
      for(int j = 0; j < dim; ++j)
         idMap[id2(j, i, numFuncs)] = idBasis[id2(i, j, dim)];

   addData->xPower[0] = NULL;
   memset(&addData->xFactor[0], 0, sizeof(RMatrix));

   addData->xPower[1]  = (int *)malloc(dim2*numFuncs*sizeof(int));
   addData->xFactor[1] = RMatrix_init(dim2, numFuncs);

   int *xPower = addData->xPower[1];
   memset(xPower, 0, dim2*numFuncs*sizeof(int));
   for(int d = 1; d < dim2; ++d)
      for(int k = 0; k < numFuncs; ++k)
         for(int i = 0; i < d; ++i)
            xPower[id2(d, k, numFuncs)] += basis->indices[k*dim+dim1+i];

   basis->table = Table3dCreate(deg, dim);
   ComputeTable(basis->table);
   return basis;
}


CubeSimplexBasis* MakeCubeSimplexCopy(CubeSimplexBasis *b)
{
   int deg = b->deg;
   int dim = b->dim;
   int dim1 = b->params->dims[0];
   int dim2 = b->params->dims[1];
   int numFuncs = b->numFuncs;

   CubeSimplexBasis *bcopy = (CubeSimplexBasis *)malloc(sizeof(CubeSimplexBasis));
   bcopy->params = (CubeSimplexParams *)malloc(sizeof(CubeSimplexParams));
   bcopy->params->deg = b->deg;
   bcopy->params->dims[0] = dim1;
   bcopy->params->dims[1] = dim2;
   bcopy->deg      = deg;
   bcopy->dim      = dim;
   bcopy->numFuncs = numFuncs;

   BasisInterface interface = SetCubeSimplexBasisInterface();
   bcopy->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(bcopy->interface) = interface;

   bcopy->indices = (INT_8 *)malloc(numFuncs*dim*sizeof(INT_8));
   memcpy(bcopy->indices, b->indices, numFuncs*dim*sizeof(INT_8));

   bcopy->addData = (AddDataCubeSimplex *)malloc(sizeof(AddDataCubeSimplex));
   AddDataCubeSimplex *addDataCp = bcopy->addData;
   AddDataCubeSimplex *addData = b->addData;

   addDataCp->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   memcpy(addDataCp->idMap, addData->idMap, numFuncs*dim*sizeof(INT_8));

   addDataCp->phi_backw1 = Vector_init(numFuncs);
   addDataCp->phi_forw1  = Vector_init(numFuncs);
   addDataCp->basis_polytopic = Vector_init(numFuncs);

   addDataCp->xPower[0] = NULL;
   memset(&addDataCp->xFactor[0], 0, sizeof(RMatrix));

   addDataCp->xFactor[1] = RMatrix_init(dim2, numFuncs);
   addDataCp->xPower[1]  = (int *)malloc(dim2*numFuncs*sizeof(int));
   memcpy(addDataCp->xPower[1], addData->xPower[1], dim2*numFuncs*sizeof(int));

   bcopy->table = Table3dCreate(deg, dim);
   CopyTable(b->table, bcopy->table);

   bcopy->functions   = Vector_init(numFuncs);
   bcopy->derivatives = Vector_init(numFuncs*dim);
   bcopy->integrals   = Vector_init(numFuncs);

   return bcopy;
}


void CubeSimplexBasisFuncs(CubeSimplexBasis *basis, const double *x, Vector F)
{
   int deg      = basis->deg;
   int dim1     = basis->params->dims[0];
   int numFuncs = basis->numFuncs;
   INT_8 *idMap = basis->addData->idMap;

   VSetToOne(F);
   for(int d = 0; d < dim1; ++d)
   {
      double legendre[deg+1];
      LegendrePoly(deg+1, 2*x[d]-1, legendre);
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         F.id[k] *= legendre[idMap[d*numFuncs+k]];
   }

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis->addData->basis_polytopic);
   Vector basis_polytopic = basis->addData->basis_polytopic;
   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      F.id[k] *= basis_polytopic.id[k];
}


void CubeSimplexBasisDer(CubeSimplexBasis *basis, const double *x, Vector dF)
{
   int deg      = basis->deg;
   int dim      = basis->dim;
   int dim1     = basis->params->dims[0];
   int dim2     = basis->params->dims[1];
   int numFuncs = basis->numFuncs;
   INT_8 *idMap = basis->addData->idMap;

   double legendre[dim1][deg+1];
   double dxlegendre[dim1][deg+1];
   for(int d = 0; d < dim1; ++d)
      LegendrePolyAndPrime(deg+1, 2*x[d]-1, legendre[d], dxlegendre[d]);

   VSetToOne(dF);
   Tensor2D dF2D = VectorToTensor2D(dim, numFuncs, dF);
   for(int d = 0; d < dim; ++d)
   {
      for(int j = 0; j < dim1; ++j)
      {
         if(j != d)
            #pragma omp simd
            for(int k = 0; k < numFuncs; ++k)
               TID2(dF2D, d, k) *= legendre[j][idMap[j*numFuncs+k]];
         if(j == d)
            #pragma omp simd
            for(int k = 0; k < numFuncs; ++k)
               TID2(dF2D, d, k) *= 2.0 * dxlegendre[j][idMap[j*numFuncs+k]];
      }
   }

   double h = POW_DOUBLE(10, -6);
   double diffFact = 1.0/(2.0*h);
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;
   Vector basis_polytopic = basis->addData->basis_polytopic;

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim1; ++d)
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         TID2(dF2D, d, k) *= basis_polytopic.id[k];

   for(int d = 0; d < dim2; ++d)
   {
      double x_backw1[dim];
      double x_forw1[dim];
      for(int d_t = 0; d_t < dim; ++d_t)
      {
         x_backw1[d_t] = x[d_t];
         x_forw1[d_t]  = x[d_t];
      }
      x_backw1[d+dim1] = x_backw1[d+dim1] - h;
      x_forw1[d+dim1]  = x_forw1[d+dim1] + h;

      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_backw1, phi_backw1);
      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_forw1, phi_forw1);
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         TID2(dF2D, d+dim1, k) *= (phi_forw1.id[k] - phi_backw1.id[k]) * diffFact;
   }
}


void CubeSimplexBasisIntegrals(CubeSimplexBasis *basis, Vector I)
{
   VSetToZero(I);
   I.id[0] = 1.0 / (double)factorial(basis->params->dims[1]);
}


void CubeSimplexBasisIntegralsMonomial(CubeSimplexBasis *basis, Vector I)
{
  int numFuncs = basis->numFuncs;
  Vector integralsCube = Vector_init(numFuncs);
  Vector integralsSimplex = Vector_init(numFuncs);
  Vector integrals = I;

  IntegralsCubePolyhedralMonomial((MixedPolytopeBasis *)basis, integralsCube);
  IntegralsSimplexPolyhedralMonomialTwo((MixedPolytopeBasis *)basis, integralsSimplex);

  #pragma omp simd
  for(int i = 0; i < numFuncs; ++i)
     integrals.id[i] = integralsCube.id[i]*integralsSimplex.id[i];

  Vector_free(integralsCube);
  Vector_free(integralsSimplex);
}


void CubeSimplexBasisFree(CubeSimplexBasis *basis)
{
   if(basis == NULL) return;

   Table3dFree(basis->table);
   if(basis->params) {
      free(basis->params);
      basis->params = NULL;
   }
   if(basis->addData) {
      AddDataCubeSimplex *addData = basis->addData;
      Vector_free(addData->basis_polytopic);
      Vector_free(addData->phi_backw1);
      Vector_free(addData->phi_forw1);
      free(addData->idMap);
      RMatrix_free(addData->xFactor[1]);
      free(addData->xPower[1]);
      free(addData);
      basis->addData = NULL;
   }
   if(basis->indices) {
      free(basis->indices);
      basis->indices = NULL;
   }
   free(basis);
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetSimplexSimplexBasisInterface()
{
   BasisInterface ssInterface;
   ssInterface.basisInit        = (BasisInitPtr)&SimplexSimplexBasisInit;
   ssInterface.makeBasisCopy    = (MakeBasisCopy)&MakeSimplexSimplexCopy;
   ssInterface.computeFuncs     = (BasisFuncsPtr)&SimplexSimplexBasisFuncs;
   ssInterface.computeDer       = (BasisDerPtr)&SimplexSimplexBasisDer;
   ssInterface.computeIntegrals = (BasisIntegralsPtr)&SimplexSimplexBasisIntegrals;
   ssInterface.computeIntegralsMonomial = (BasisIntegralsMonomialPtr)&SimplexSimplexBasisIntegralsMonomial;
   ssInterface.basisFree        = (BasisFreePtr)&SimplexSimplexBasisFree;
   return ssInterface;
}


SimplexSimplexBasis* SimplexSimplexBasisInit(SimplexSimplexParams *ssParams)
{
   int deg = ssParams->deg;
   int dim1 = ssParams->dims[0];
   int dim2 = ssParams->dims[1];
   int dim = dim1+dim2;
   assert(dim1 >= 2 && dim2 >= 2 && deg >= 1);

   SimplexSimplexBasis *basis = (SimplexSimplexBasis *)malloc(sizeof(SimplexSimplexBasis));
   basis->params = (SimplexSimplexParams *)malloc(sizeof(SimplexSimplexParams));

   basis->params->deg         = deg;
   basis->params->dims[0]     = dim1;
   basis->params->dims[1]     = dim2;
   basis->deg                 = deg;
   basis->dim                 = dim;
   basis->numFuncs            = BasisSize(deg, dim);

   basis->indices = (INT_8 *) malloc(basis->numFuncs*dim*sizeof(INT_8));
   BasisIndices(basis->deg, basis->dim, basis->indices);

   basis->addData = (AddDataSimplexSimplex *)malloc(sizeof(AddDataSimplexSimplex));
   int numFuncs = basis->numFuncs;
   AddDataSimplexSimplex *addData = basis->addData;
   addData->phi_backw1 = Vector_init(numFuncs);
   addData->phi_forw1  = Vector_init(numFuncs);
   addData->basis_polytopic = Vector_init(numFuncs);

   addData->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   INT_8 *idMap   = addData->idMap;
   INT_8 *idBasis = basis->indices;
   for(int i = 0; i < numFuncs; ++i)
      #pragma omp simd
      for(int j = 0; j < dim; ++j)
         idMap[id2(j, i, numFuncs)] = idBasis[id2(i, j, dim)];

   addData->xPower[0]  = (int *)malloc(dim1*numFuncs*sizeof(int));
   addData->xFactor[0] = RMatrix_init(dim1, numFuncs);
   addData->xPower[1]  = (int *)malloc(dim2*numFuncs*sizeof(int));
   addData->xFactor[1] = RMatrix_init(dim2, numFuncs);

   int *xPowerOne = addData->xPower[0];
   int *xPowerTwo = addData->xPower[1];
   memset(xPowerOne, 0, dim1*numFuncs*sizeof(int));
   for(int d = 1; d < dim1; ++d)
      for(int k = 0; k < numFuncs; ++k)
         for(int i = 0; i < d; ++i)
            xPowerOne[id2(d, k, numFuncs)] += basis->indices[k*dim+i];

   memset(xPowerTwo, 0, dim2*numFuncs*sizeof(int));
   for(int d = 1; d < dim2; ++d)
      for(int k = 0; k < numFuncs; ++k)
         for(int i = 0; i < d; ++i)
            xPowerTwo[id2(d, k, numFuncs)] += basis->indices[k*dim+dim1+i];

   basis->table = Table3dCreate(deg, dim);
   ComputeTable(basis->table);

   return basis;
}


SimplexSimplexBasis* MakeSimplexSimplexCopy(SimplexSimplexBasis *b)
{
   int deg = b->deg;
   int dim = b->dim;
   int dim1 = b->params->dims[0];
   int dim2 = b->params->dims[1];
   int numFuncs = b->numFuncs;

   SimplexSimplexBasis *bcopy = (SimplexSimplexBasis *)malloc(sizeof(SimplexSimplexBasis));
   bcopy->params = (SimplexSimplexParams *)malloc(sizeof(SimplexSimplexParams));
   bcopy->params->deg = deg;
   bcopy->params->dims[0] = dim1;
   bcopy->params->dims[1] = dim2;
   bcopy->deg      = deg;
   bcopy->dim      = dim;
   bcopy->numFuncs = numFuncs;

   BasisInterface interface = SetSimplexSimplexBasisInterface();
   bcopy->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(bcopy->interface) = interface;

   bcopy->indices = (INT_8 *)malloc(numFuncs*dim*sizeof(INT_8));
   memcpy(bcopy->indices, b->indices, numFuncs*dim*sizeof(INT_8));

   bcopy->addData = (AddDataSimplexSimplex *)malloc(sizeof(AddDataSimplexSimplex));
   AddDataSimplexSimplex *addDataCp = bcopy->addData;
   AddDataSimplexSimplex *addData = b->addData;

   addDataCp->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   memcpy(addDataCp->idMap, addData->idMap, numFuncs*dim*sizeof(INT_8));

   addDataCp->phi_backw1 = Vector_init(numFuncs);
   addDataCp->phi_forw1  = Vector_init(numFuncs);
   addDataCp->basis_polytopic = Vector_init(numFuncs);

   addDataCp->xFactor[0] = RMatrix_init(dim1, numFuncs);
   addDataCp->xPower[0]  = (int *)malloc(dim1*numFuncs*sizeof(int));
   memcpy(addDataCp->xPower[0], addData->xPower[0], dim1*numFuncs*sizeof(int));

   addDataCp->xFactor[1] = RMatrix_init(dim2, numFuncs);
   addDataCp->xPower[1]  = (int *)malloc(dim2*numFuncs*sizeof(int));
   memcpy(addDataCp->xPower[1], addData->xPower[1], dim2*numFuncs*sizeof(int));

   bcopy->table = Table3dCreate(deg, dim);
   CopyTable(b->table, bcopy->table);

   bcopy->functions   = Vector_init(numFuncs);
   bcopy->derivatives = Vector_init(numFuncs*dim);
   bcopy->integrals   = Vector_init(numFuncs);

   return bcopy;
}


void SimplexSimplexBasisFuncs(SimplexSimplexBasis *basis, const double *x, Vector F)
{
   int numFuncs     = basis->numFuncs;
   Vector polytopic = basis->addData->basis_polytopic;

   VSetToOne(F);
   SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x, polytopic);
   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      F.id[k] *= polytopic.id[k];

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, polytopic);
   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      F.id[k] *= polytopic.id[k];
}


void SimplexSimplexBasisDer(SimplexSimplexBasis *basis, const double *x, Vector dF)
{
   int dim1     = basis->params->dims[0];
   int dim2     = basis->params->dims[1];
   int dim      = basis->dim;
   int numFuncs = basis->numFuncs;

   double h = POW_DOUBLE(10, -6) * 5.0;
   double diffFact = 1.0/(2.0*h);
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;

   double x_backw1[dim];
   double x_forw1[dim];

   VSetToOne(dF);
   Vector basis_polytopic      = basis->addData->basis_polytopic;
   Tensor2D  dF2D              = VectorToTensor2D(dim, numFuncs, dF);
   Tensor1D  basis_polytopic1D = VectorToTensor1D(basis_polytopic);

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim1; ++d)
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         TID2(dF2D, d, k) *= TID1(basis_polytopic1D, k);

   SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim2; ++d)
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         TID2(dF2D, d+dim1, k) *= TID1(basis_polytopic1D, k);

   for(int d = 0; d < dim2; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t)
      {
         x_backw1[d_t] = x[d_t];
         x_forw1[d_t]  = x[d_t];
      }
      x_backw1[d+dim1] = x_backw1[d+dim1] - h;
      x_forw1[d+dim1]  = x_forw1[d+dim1] + h;

      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_backw1, phi_backw1);
      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_forw1, phi_forw1);

      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         TID2(dF2D, d+dim1, k) *= (phi_forw1.id[k] - phi_backw1.id[k]) * diffFact;
   }

   for(int d = 0; d < dim1; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t)
      {
         x_backw1[d_t] = x[d_t];
         x_forw1[d_t] = x[d_t];
      }
      x_backw1[d] = x_backw1[d] - h;
      x_forw1[d] = x_forw1[d] + h;

      SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x_backw1, phi_backw1);
      SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x_forw1, phi_forw1);

      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         TID2(dF2D, d, k) *= (phi_forw1.id[k] - phi_backw1.id[k]) * diffFact;
   }
}


void SimplexSimplexBasisIntegrals(SimplexSimplexBasis *basis, Vector I)
{
   VSetToZero(I);
   I.id[0] = 1.0 / (double)factorial(basis->params->dims[0]);
   I.id[0] /= (double)factorial(basis->params->dims[1]);
}


void SimplexSimplexBasisIntegralsMonomial(SimplexSimplexBasis *basis, Vector I)
{
  int numFuncs = basis->numFuncs;
  Vector integralsS1 = Vector_init(numFuncs);
  Vector integralsS2 = Vector_init(numFuncs);
  Vector integrals = I;

  IntegralsSimplexPolyhedralMonomialOne((MixedPolytopeBasis *)basis, integralsS1);
  IntegralsSimplexPolyhedralMonomialTwo((MixedPolytopeBasis *)basis, integralsS2);
  VMult(integralsS1, integralsS2, integrals);

  Vector_free(integralsS1);
  Vector_free(integralsS2);
}


void SimplexSimplexBasisFree(SimplexSimplexBasis *basis)
{
   if(basis == NULL) return;

   Table3dFree(basis->table);
   if(basis->params) {
      free(basis->params);
      basis->params = NULL;
   }
   if(basis->addData) {
      AddDataSimplexSimplex *addData = basis->addData;
      Vector_free(addData->basis_polytopic);
      Vector_free(addData->phi_backw1);
      Vector_free(addData->phi_forw1);
      free(addData->idMap);
      RMatrix_free(addData->xFactor[0]);
      RMatrix_free(addData->xFactor[1]);
      free(addData->xPower[0]);
      free(addData->xPower[1]);
      free(addData);
      basis->addData = NULL;
   }
   if(basis->indices) {
      free(basis->indices);
      basis->indices = NULL;
   }
   free(basis);
}


static void LegendrePoly(int order, double x, double *p)
{
   p[0] = 1.0;
   p[1] = x;

   for(int  k = 1; k < order-1; ++k)
   {
      double fac1 = (2.0*k+1.0)/(k+1.0);
      double fac2 = k/(k+1.0);
      p[k+1] = fac1*x*p[k] - fac2*p[k-1];
   }
}

static void LegendrePolyAccurate(int order, double x, double *p)
{
   double128 pquad[order];
   pquad[0] = 1.q;
   pquad[1] = x;

   for(int  k = 1; k < order-1; ++k)
   {
      double128 fac1 = (2.q*k+1.q)/(k+1.q);
      double128 fac2 = k/(k+1.q);
      pquad[k+1] = fac1*x*pquad[k] - fac2*pquad[k-1];
   }

   for(int  k = 0; k < order; ++k)
      p[k] = pquad[k];
}


static void LegendrePolyAndPrime(int order, double x, double *p, double *dp)
{
   p[0] = 1.0; p[1] = x;
   dp[0] = 0.0; dp[1] = 1.0;

   for(int  k = 1; k < order-1; ++k)
   {
      double fac1 = (2.0*k+1.0)/(k+1.0);
      double fac2 = k/(k+1.0);
      p[k+1] = fac1*x*p[k] - fac2*p[k-1];
      dp[k+1] = fac1*(p[k] + x*dp[k]) - fac2*dp[k-1];
   }
}

static void LegendrePolyAndPrimeAccurate(int order, double x, double *p, double *dp)
{
   double128 pquad[order];
   double128 dpquad[order];
   pquad[0] = 1.q; pquad[1] = x;
   dpquad[0] = 0.q; dpquad[1] = 1.q;

   for(int  k = 1; k < order-1; ++k)
   {
      double fac1 = (2.q*k+1.q)/(k+1.q);
      double fac2 = k/(k+1.q);
      pquad[k+1] = fac1*x*pquad[k] - fac2*pquad[k-1];
      dpquad[k+1] = fac1*(pquad[k] + x*dpquad[k]) - fac2*dpquad[k-1];
   }

   for(int  k = 0; k < order; ++k) p[k] = pquad[k];
   for(int  k = 0; k < order; ++k) dp[k] = dpquad[k];
}


static Table3d Table3dCreate(int deg, int dim)
{
   Table3d table;
   if(dim >= 3)      table.size[0] = 2*deg+2+dim-1;
   else if(dim == 2) table.size[0] = 2*deg+2;
   table.size[1] = deg+2;
   table.size[2] = 3;

   table.id = (double ***)malloc(table.size[0]*sizeof(double **));
   for(int i = 0; i < table.size[0]; ++i)
      table.id[i] = (double **)malloc(table.size[1]*sizeof(double *));
   for(int i = 0; i < table.size[0]; ++i)
      for(int j = 0; j < table.size[1]; ++j)
         table.id[i][j] = (double *)malloc(table.size[2]*sizeof(double));
   return table;
}


static void Table3dFree(Table3d table)
{
   for(int i = 0; i < table.size[0]; ++i)
      for(int j = 0; j < table.size[1]; ++j)
         free(table.id[i][j]);
   for(int i = 0; i < table.size[0]; ++i)
      free(table.id[i]);
   free(table.id);
}


static void ComputeTable(Table3d table)
{
   for(int alpha = 1; alpha < table.size[0]; ++alpha)
   {
      for (int k = 1; k < table.size[1]; ++k)
      {
         double128 fac1Part0 = (alpha+2.q*k+1.q) * (alpha+2.q*k+2.q) * (alpha+2.q*k);
         double128 fac1Part1 = (alpha+2.q*k+1.q) * (alpha*alpha);                      //excludes multiplication by x
         double128 fac2 = -2.q * (alpha+k) * k * (alpha+2.q*k+2.q);
         double128 fac3 = 1.q / (2.q*(k+1.q)*(alpha+k+1.q) * (alpha+2.q*k));
         table.id[alpha][k][0] = fac3 * fac1Part0;
         table.id[alpha][k][1] = fac3 * fac1Part1;
         table.id[alpha][k][2] = fac3 * fac2;
      }
   }
}


static void CopyTable(Table3d t1, Table3d t2)
{
   assert(t1.size[0] == t2.size[0]);
   assert(t1.size[1] == t2.size[1]);
   for(int alpha = 1; alpha < t1.size[0]; ++alpha)
   {
      for (int k = 1; k < t1.size[1]; ++k)
      {
         t2.id[alpha][k][0] = t1.id[alpha][k][0];
         t2.id[alpha][k][1] = t1.id[alpha][k][1];
         t2.id[alpha][k][2] = t1.id[alpha][k][2];
      }
   }
}


static void JacobiPoly(int order, double x, int alpha, double *p)
{
   p[0] = 1.0;
   p[1] = 0.5*(x-1.0) * (alpha+2.0) + alpha+1.0;
   for(int k = 1; k < order-1; ++k)
   {
      double fac1Part0 = (alpha+2.0*k+1.0) * (alpha+2.0*k+2.0) * (alpha+2.0*k);
      double fac1Part1 = (alpha+2.0*k+1.0) * (alpha*alpha);
      double fac2 = -2.0 * (alpha+k) * k * (alpha+2.0*k+2.0);
      double fac3 = 1.0 / (2.0*(k+1.0)*(alpha+k+1.0) * (alpha+2.0*k));
      double c1 = fac3 * fac1Part0;
      double c2 = fac3 * fac1Part1;
      double c3 = fac3 * fac2;
      p[k+1] = (c1*x + c2)*p[k] + c3*p[k-1];
   }
}


// asumes order is at least 2 or greater, i.e. polynomial degree is at least 1
static void JacobiPolyWithTable(int order, double x, int alpha, double *p, Table3d table)
{
#ifdef QUAD_DEBUG_ON
   assert(order >= 2);
#endif
   p[0] = 1.0;
   p[1] = 0.5*(x-1.0) * (alpha+2.0) + alpha+1.0;
   for(int k = 1; k < order-1; ++k)
      p[k+1] = (table.id[alpha][k][0]*x + table.id[alpha][k][1])*p[k] + table.id[alpha][k][2]*p[k-1];
}

// asumes order is at least 2 or greater, i.e. polynomial degree is at least 1
static void JacobiPolyWithTableAccurate(int order, double x, int alpha, double *p, Table3d table)
{
#ifdef QUAD_DEBUG_ON
   assert(order >= 2);
#endif
   double128 pquad[order];
   pquad[0] = 1.q;
   pquad[1] = 0.5q*(x-1.q) * (alpha+2.q) + alpha+1.q;
   for(int k = 1; k < order-1; ++k)
      pquad[k+1] = (table.id[alpha][k][0]*x + table.id[alpha][k][1])*pquad[k] + table.id[alpha][k][2]*pquad[k-1];

   for(int k = 0; k < order; ++k)
      p[k] = pquad[k];
}


static void IntegralsCubePolyhedralMonomial(MixedPolytopeBasis *basis, Vector I)
{
   int dim          = basis->dim;
   int dimCube      = basis->params->dims[0];
   int numFuncs     = basis->numFuncs;
   Vector integrals = I;
   INT_8 *indices   = basis->indices;

   for(int i = 0; i < numFuncs; ++i)
   {
      double val = 1.0;
      for(int d = 0; d < dimCube; ++d)
         val = val/(indices[i*dim+d] + 1);
      integrals.id[i] = val;
   }
}


static void IntegralsSimplexPolyhedralMonomialOne(MixedPolytopeBasis *basis, Vector I)
{
   int dim          = basis->dim;
   int numFuncs     = basis->numFuncs;
   int dimSimplex   = basis->params->dims[0];
   INT_8 *indices   = basis->indices;

   for(int i = 0; i < numFuncs; ++i)
   {
      double val = 1.0;
      for(int d = 0; d < dimSimplex; ++d)
      {
         double power = 0;
         for(int r = 0; r < dimSimplex-d; ++r)
            power += (double)indices[dim*i+dimSimplex-r-1];
         val /= (power+dimSimplex-d);
      }
      I.id[i] = val;
   }
}


static void IntegralsSimplexPolyhedralMonomialTwo(MixedPolytopeBasis *basis, Vector I)
{
   int dim          = basis->dim;
   int dimSimplex   = basis->params->dims[1];
   int numFuncs     = basis->numFuncs;
   INT_8 *indices   = basis->indices;

   for(int i = 0; i < numFuncs; ++i)
   {
      double val = 1.0;
      for(int d = 0; d < dimSimplex; ++d)
      {
         double power = 0;
         for(int r = 0; r < dimSimplex-d; ++r)
            power += (double)indices[dim*(i+1)-r-1];
         val /= (power+dimSimplex-d);
      }
      I.id[i] = val;
   }
}


void SimplexFuncsPolytopicOne(MixedPolytopeBasis *basis, const double *x, Vector F)
{
   int deg        = basis->deg;
   int dim1       = basis->params->dims[0];
   int numFuncs   = basis->numFuncs;
   INT_8 *idMap   = basis->addData->idMap;

   double jacobi[dim1-1][SQUARE(deg+1)];
   double legendre[(deg+1)];
   LegendrePoly(deg+1, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], legendre);

   double xCoord[dim1-1];
   for(int d = 0; d < dim1-2; ++d)
      xCoord[d] = x[dim1-d-2]/x[dim1-d-3];
   xCoord[dim1-2] = x[0];

   double jCoord[dim1-1];
   for(int d = 1; d < dim1-1; ++d)
      jCoord[d-1] = 1.0-2.0*x[dim1-d-1] / x[dim1-d-2];
   jCoord[dim1-2] = 1.0-2.0*x[0];

   for(int d = 1; d < dim1; ++d) {
      for(int j = 0; j < deg+1; ++j) {
         double alpha = 2*j+d;
         JacobiPolyWithTable(deg+1, jCoord[d-1], alpha, &jacobi[d-1][(deg+1)*j], basis->table);
      }
   }

   int *xPow = basis->addData->xPower[0];
   RMatrix xFactor = basis->addData->xFactor[0];
   for(int d = 1; d < dim1; ++d)
      for(int k = 0; k < numFuncs; ++k)
            xFactor.rid[d][k] = DoubleIntPower(xCoord[d-1], xPow[id2(d, k, numFuncs)]);

   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      F.id[k] = legendre[idMap[k]];

   for(int d = 1; d < dim1; ++d)
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         F.id[k] *= jacobi[d-1][idMap[d*numFuncs+k] + (deg+1)*xPow[d*numFuncs+k]] * xFactor.rid[d][k];
}


void SimplexFuncsPolytopicTwo(MixedPolytopeBasis *basis, const double *x, Vector F)
{
   int deg        = basis->deg;
   int dim        = basis->dim;
   int dim1       = basis->params->dims[0];
   int dim2       = basis->params->dims[1];
   int dimTwo     = dim1+dim2;
   int numFuncs   = basis->numFuncs;
   INT_8 *idMap   = basis->addData->idMap;

   double jacobi[dim2-1][SQUARE(deg+1)];
   double legendre[deg+1];
   LegendrePoly(deg+1, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], legendre);

   double xCoord[dim2-1];
   for(int d = 0; d < dim2-2; ++d)
      xCoord[d] = x[dimTwo-d-2]/x[dimTwo-d-3];
   xCoord[dim2-2] = x[dim1];

   double jCoord[dim];
   for(int d = 1; d < dim2; ++d)
      jCoord[d-1] = 1.0-2.0*x[dimTwo-d-1]/x[dimTwo-d-2];
   jCoord[dim2-2] = 1.0-2.0*x[dim1];

   for(int d = 1; d < dim2; ++d) {
      for(int j = 0; j < deg+1; ++j) {
         int nextAlpha = (deg+1)*j;
         JacobiPolyWithTable(deg+1, jCoord[d-1], 2*j+d, &jacobi[d-1][nextAlpha], basis->table);
      }
   }

   RMatrix xFactor = basis->addData->xFactor[1];
   int *xPower = basis->addData->xPower[1];
   for(int d = 1; d < dim2; ++d)
      for(int k = 0; k < numFuncs; ++k)
            xFactor.rid[d][k] = DoubleIntPower(xCoord[d-1], xPower[id2(d, k, numFuncs)]);

   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      F.id[k] = legendre[idMap[(dim1)*numFuncs+k]];

   for(int d = 1; d < dim2; ++d)
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         F.id[k] *= jacobi[d-1][idMap[(d+dim1)*numFuncs+k] + (deg+1)*xPower[d*numFuncs+k]] * xFactor.rid[d][k];
}


double orthogonal_cube_basis_test(int deg, int dim)
{
   quadrature *quad_C = quadrature_full_cube_tensor(deg, dim);
   double res = QuadTestIntegral(quad_C, orthogonal);
   printf("\nTesting orthogonality of cube basis functions."
          "Maximum error of basis functions = %.16e\n\n", res);
   quadrature_free(quad_C);
   return res;
}


double orthogonal_simplex_basis_test(int deg, int dim)
{
   quadrature *quad_S = quadrature_full_simplex_tensor(deg, dim);
   double res = QuadTestIntegral(quad_S, orthogonal);
   printf("\nTesting orthogonality of simplex basis functions."
          "Maximum error of basis functions = %.16e\n\n", res);
   quadrature_free(quad_S);
   return res;
}


void PrintBasisIndices(Basis *basis)
{
   INT_8 *indices = basis->indices;
   for(int i = 0; i < basis->numFuncs; ++i) {
      for(int j = 0; j < basis->dim; ++j)
         printf("%i    ", indices[i*basis->dim+j]);
      printf("\n");
   }
}


