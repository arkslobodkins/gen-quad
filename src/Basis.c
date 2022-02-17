#include "Basis.h"
#include "Quadrature.h"
#include "Gauss_Lib/Jacobi.h"
#include "GeneralGaussTensor.h"
#include "AddDimension.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define id2(i, j, N) (i)*(N)+(j)

static int BasisSize(int deg, int dim);
static void BasisIndices(int deg, int dim, INT_8 *F);

static Table3d Table3dCreate(int deg, int dim);
static void Table3dFree(Table3d table);
static void ComputeTable(Table3d table);

static void LegendrePoly(int order, double x, double *p);
static void LegendrePolyAndPrime(int order, double x, double *p, double *dp);
static void JacobiPolyWithTable(int order, double x, int alpha, double *p, Table3d table);

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

void BasisFuncs(Basis *basis, const double *x, Vector v)
{
   basis->interface->computeFuncs(basis, x, v);
}

void BasisDer(Basis *basis, const double *x, Vector v)
{
   basis->interface->computeDer(basis, x, v);
}

void BasisIntegrals(Basis *basis, Vector v)
{
   basis->interface->computeIntegrals(basis, v);
}

void BasisIntegralsMonomial(Basis *basis, Vector v)
{
   basis->interface->computeIntegralsMonomial(basis, v);
}

void BasisFree(Basis *basis)
{
   if(basis == NULL) return;
   BasisInterface *interface = basis->interface;

   Vector_free(basis->integrals);
   Vector_free(basis->derivatives);
   Vector_free(basis->functions);

   basis->interface->basisFree(basis);
   if(interface) {
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
      for(int i = 0; i <= deg; ++i) F[i] = i;
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

void BasisMonomial(Basis *basis, const double *x, Vector phi)
{
   int numFuncs    = basis->numFuncs;
   int dim         = basis->dim;
   INT_8 *basis_id = basis->indices;

   VSetToOne(phi);
   for(int k = 0; k < numFuncs; ++k) {
      for(int d = 0; d < dim; ++d) {
         INT_8 basis_power = basis_id[k*dim+d];
         phi.id[k] *= DoubleIntPower(x[d], basis_power);
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


void ComputeCubeBasisFuncs(CubeBasis *basis, const double *x, Vector v)
{
   int deg        = basis->deg;
   int dim        = basis->dim;
   int numFuncs   = basis->numFuncs;
   INT_8 *idMap   = basis->addData->idMap;
   double *phi    = v.id;

   double legendre[deg+1];
   LegendrePoly(deg+1, 2*x[0]-1, legendre);
   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      phi[k] = legendre[idMap[k]];

   for(int d = 1; d < dim; ++d)
   {
      LegendrePoly(deg+1, 2*x[d]-1, legendre);
      int nextDim = d*numFuncs;
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         phi[k] *= legendre[idMap[nextDim+k]];
   }
}


void ComputeCubeBasisDer(CubeBasis *basis, const double *x, Vector v)
{
   int deg = basis->deg;
   int dim = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *idMap = basis->addData->idMap;

   double legendre[(deg+1)*dim];
   double dxlegendre[(deg+1)*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePolyAndPrime(deg+1, 2*x[d]-1, &legendre[d*(deg+1)], &dxlegendre[d*(deg+1)]);

   VSetToOne(v);
   double *phiPrime = v.id;
   for(int d = 0; d < dim; ++d)
   {
      // dimension < d
      for(int j = 0; j < d; ++j) {
         double *legendre_ptr = &legendre[(deg+1)*j];
         for(int k = 0; k < numFuncs; k++)
            phiPrime[d*numFuncs+k] *= legendre_ptr[idMap[j*numFuncs+k]];
      }
      // dimension = d
      double *dxlegendre_ptr = &dxlegendre[(deg+1)*d];
      for(int k = 0; k < numFuncs; ++k)
         phiPrime[d*numFuncs+k] *= 2.0 * dxlegendre_ptr[idMap[d*numFuncs+k]];

      // dimension > d
      for(int j = d+1; j < dim; ++j) {
         double *legendre_ptr = &legendre[(deg+1)*j];
         for(int k = 0; k < numFuncs; k++)
            phiPrime[d*numFuncs+k] *= legendre_ptr[idMap[j*numFuncs+k]];
      }
   }
}


void CubeBasisIntegrals(CubeBasis *basis, Vector v)
{
   int len = basis->integrals.len;
   double *integrals = v.id;
   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1.0;
}


void CubeBasisIntegralsMonomial(CubeBasis *basis, Vector v)
{
   int dim = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *basisIndices = basis->indices;
   Vector integrals = v;

   for(int i = 0; i < numFuncs; ++i) {
      double val = 1.0;
      for(int d = 0; d < dim; ++d)
         val = val/(basisIndices[i*dim+d] + 1);
      integrals.id[i] = val;
   }
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
   free(basis); basis = NULL;
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


SimplexBasis* MakeSimplexCopy(SimplexBasis *basis)
{
   SimplexBasis *basisCopy = (SimplexBasis *)malloc(sizeof(SimplexBasis));
   basisCopy->params       = (SimplexParams *)malloc(sizeof(SimplexParams));
   basisCopy->params->deg  = basis->deg;
   basisCopy->params->dim  = basis->dim;
   basisCopy->deg          = basis->deg;
   basisCopy->dim          = basis->dim;
   basisCopy->numFuncs     = basis->numFuncs;

   BasisInterface interface = SetSimplexBasisInterface();
   basisCopy->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(basisCopy->interface) = interface;

   int deg = basisCopy->deg;
   int dim = basisCopy->dim;
   int numFuncs = basisCopy->numFuncs;

   basisCopy->indices = (INT_8 *)malloc(numFuncs*dim*sizeof(INT_8));
   memcpy(basisCopy->indices, basis->indices, numFuncs*dim*sizeof(INT_8));

   basisCopy->addData = (AddDataSimplex *)malloc(sizeof(AddDataSimplex));
   basisCopy->addData->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   memcpy(basisCopy->addData->idMap, basis->addData->idMap, numFuncs*dim*sizeof(INT_8));

   basisCopy->addData->phi_backw1 = Vector_init(numFuncs);
   basisCopy->addData->phi_forw1  = Vector_init(numFuncs);

   basisCopy->addData->xFactor = RMatrix_init(dim, numFuncs);

   basisCopy->addData->xPower  = (int *)malloc(dim*numFuncs*sizeof(int));
   memcpy(basisCopy->addData->xPower, basis->addData->xPower, dim*numFuncs*sizeof(int));

   basisCopy->table = Table3dCreate(deg, dim);
   for(int alpha = 1; alpha < basis->table.size[0]; ++alpha) {
      for (int k = 1; k < basis->table.size[1]; ++k) {
         basisCopy->table.id[alpha][k][0] = basis->table.id[alpha][k][0];
         basisCopy->table.id[alpha][k][1] = basis->table.id[alpha][k][1];
         basisCopy->table.id[alpha][k][2] = basis->table.id[alpha][k][2];
      }
   }

   basisCopy->functions   = Vector_init(numFuncs);
   basisCopy->derivatives = Vector_init(numFuncs*dim);
   basisCopy->integrals   = Vector_init(numFuncs);

   return basisCopy;
}


void SimplexBasisFuncs(SimplexBasis *basis, const double *x, Vector v)
{
   int deg        = basis->deg;
   int dim        = basis->dim;
   int numFuncs   = basis->numFuncs;
   double legendre[(deg+1)];
   double jacobi[SQUARE(deg+1)*(dim-1)];

   INT_8 *idMap   = basis->addData->idMap;
   INT_8 *basisId = basis->indices;

   LegendrePoly(deg+1, (2.0*x[dim-1]-x[dim-2])/x[dim-2], legendre);
   double xCoord[dim-1];
   for(int d = 0; d < dim-2; ++d) xCoord[d] = x[dim-d-2]/x[dim-d-3];
   xCoord[dim-2] = x[0];

   double jCoord[dim];
   for(int d = 1; d < dim-1; ++d)
      jCoord[d] = 1.0-2.0*x[dim-d-1]/x[dim-d-2];
   jCoord[dim-1] = 1.0-2.0*x[0];

   for(int d = 1; d < dim; ++d) {
      int nextDim = (deg+1)*(d-1)*deg;
      for(int j = 0; j < deg+1; ++j) {
         int nextAlpha = (deg+1)*j;
         JacobiPolyWithTable(deg+1, jCoord[d], 2*j+d, &jacobi[nextDim+nextAlpha], basis->table);
      }
   }

   double *phi = v.id;
   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      phi[k] = legendre[idMap[k]];

   RMatrix xFactor = basis->addData->xFactor;
   int *xPower = basis->addData->xPower;
   for(int d = 1; d < dim; ++d)
      for(int k = 0; k < numFuncs; ++k)
            xFactor.rid[d][k] = DoubleIntPower(xCoord[d-1], xPower[id2(d, k, numFuncs)]);

   for(int d = 1; d < dim; ++d) {
      int nextDim = (deg+1)*deg*(d-1);
      for(int k = 0; k < numFuncs; ++k)
         phi[k] *= jacobi[basisId[k*dim+d] + nextDim + (deg+1)*xPower[id2(d, k, numFuncs)]] * xFactor.rid[d][k];
   }
}


void SimplexBasisDer(SimplexBasis *basis, const double *x, Vector v)
{
   int dim = basis->params->dim;
   int numFuncs = basis->numFuncs;
   double h = POW_DOUBLE(10, -6);

   double x_backw1[dim];
   double x_forw1[dim];

   double *basisDer = v.id;
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;

   for(int d = 0; d < dim; ++d)
   {
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
         basisDer[d*numFuncs+k] = phi_forw1.id[k] - phi_backw1.id[k];
   }
   VectorScale(1.0/(2.0*h), v);
}


void SimplexBasisIntegrals(SimplexBasis *basis, Vector v)
{
   int dim = basis->dim;
   int len = basis->integrals.len;
   double *integrals = v.id;

   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1;
   for(int i = 1; i <= dim; ++i)
      integrals[0] /= i;

}


void SimplexBasisIntegralsMonomial(SimplexBasis *basis, Vector v)
{
   int dim             = basis->dim;
   int numFuncs        = basis->numFuncs;
   INT_8 *basisIndices = basis->indices;
   Vector integrals    = v;
   double val, power;

   for(int i = 0; i < numFuncs; ++i)
   {
      val = 1.0;
      for(int d = 0; d < dim; ++d)
      {
         power = 0;
         for(int r = 0; r < dim-d; ++r)
            power += (double)basisIndices[i*dim+dim-r-1];
         val /= (power+dim-d);
      }
      integrals.id[i] = val;
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
      Vector_free(basis->addData->phi_backw1);
      Vector_free(basis->addData->phi_forw1);
      free(basis->addData->idMap);
      RMatrix_free(basis->addData->xFactor);
      free(basis->addData->xPower);
      free(basis->addData);
      basis->addData = NULL;
   }
   if(basis->indices) {
      free(basis->indices);
      basis->indices = NULL;
   }
   free(basis); basis = NULL;
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
   assert(dim1 >= 1  && dim2 >= 2 && deg >= 1);

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


CubeSimplexBasis* MakeCubeSimplexCopy(CubeSimplexBasis *basis)
{
   CubeSimplexBasis *basisCopy = (CubeSimplexBasis *)malloc(sizeof(CubeSimplexBasis));
   basisCopy->params = (CubeSimplexParams *)malloc(sizeof(CubeSimplexParams));
   basisCopy->params->deg = basis->deg;
   basisCopy->params->dims[0] = basis->params->dims[0];
   basisCopy->params->dims[1] = basis->params->dims[1];
   basisCopy->deg         = basis->deg;
   basisCopy->dim         = basis->dim;
   basisCopy->numFuncs    = basis->numFuncs;

   BasisInterface interface = SetCubeSimplexBasisInterface();
   basisCopy->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(basisCopy->interface) = interface;

   int deg = basisCopy->deg;
   int dim = basisCopy->dim;
   int dim2 = basisCopy->params->dims[1];
   int numFuncs = basisCopy->numFuncs;

   basisCopy->indices = (INT_8 *)malloc(numFuncs*dim*sizeof(INT_8));
   memcpy(basisCopy->indices, basis->indices, numFuncs*dim*sizeof(INT_8));

   basisCopy->addData = (AddDataCubeSimplex *)malloc(sizeof(AddDataCubeSimplex));
   basisCopy->addData->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   memcpy(basisCopy->addData->idMap, basis->addData->idMap, numFuncs*dim*sizeof(INT_8));

   basisCopy->addData->phi_backw1 = Vector_init(numFuncs);
   basisCopy->addData->phi_forw1  = Vector_init(numFuncs);
   basisCopy->addData->basis_polytopic = Vector_init(numFuncs);

   basisCopy->addData->xPower[0] = NULL;
   memset(&basisCopy->addData->xFactor[0], 0, sizeof(RMatrix));

   basisCopy->addData->xFactor[1] = RMatrix_init(dim2, numFuncs);
   basisCopy->addData->xPower[1]  = (int *)malloc(dim2*numFuncs*sizeof(int));
   memcpy(basisCopy->addData->xPower[1], basis->addData->xPower[1], dim2*numFuncs*sizeof(int));

   basisCopy->table = Table3dCreate(deg, dim);
   for(int alpha = 1; alpha < basis->table.size[0]; ++alpha) {
      for (int k = 1; k < basis->table.size[1]; ++k) {
         basisCopy->table.id[alpha][k][0] = basis->table.id[alpha][k][0];
         basisCopy->table.id[alpha][k][1] = basis->table.id[alpha][k][1];
         basisCopy->table.id[alpha][k][2] = basis->table.id[alpha][k][2];
      }
   }

   basisCopy->functions   = Vector_init(numFuncs);
   basisCopy->derivatives = Vector_init(numFuncs*dim);
   basisCopy->integrals   = Vector_init(numFuncs);

   return basisCopy;
}


void CubeSimplexBasisFuncs(CubeSimplexBasis *basis, const double *x, Vector v)
{
   int deg        = basis->deg;
   int dim1       = basis->params->dims[0];
   int numFuncs   = basis->numFuncs;
   INT_8 *idMap   = basis->addData->idMap;

   VSetToOne(v);
   double *phi = v.id;
   for(int d = 0; d < dim1; ++d)
   {
      double legendre[deg+1];
      LegendrePoly(deg+1, 2*x[d]-1, legendre);
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         phi[k] *= legendre[idMap[d*numFuncs+k]];
   }

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis->addData->basis_polytopic);
   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      phi[k] *= basis->addData->basis_polytopic.id[k];
}


void CubeSimplexBasisDer(CubeSimplexBasis *basis, const double *x, Vector v)
{
   int deg          = basis->deg;
   int dim          = basis->dim;
   int dim1         = basis->params->dims[0];
   int dim2         = basis->params->dims[1];
   int numFuncs     = basis->numFuncs;
   INT_8 *idMap     = basis->addData->idMap;

   double legendre[(deg+1)*dim1];
   double dxlegendre[(deg+1)*dim1];

   for(int d = 0; d < dim1; ++d)
      LegendrePolyAndPrime(deg+1, 2*x[d]-1, &legendre[d*(deg+1)], &dxlegendre[d*(deg+1)]);

   VSetToOne(v);
   double *phiPrime = v.id;

   for(int d = 0; d < dim; ++d)
   {
      for(int j = 0; j < dim1; ++j)
      {
         if(j != d)
            #pragma omp simd
            for(int k = 0; k < numFuncs; ++k)
               phiPrime[k+d*numFuncs] *= legendre[idMap[j*numFuncs+k]+(deg+1)*j];
         if(j == d)
            #pragma omp simd
            for(int k = 0; k < numFuncs; ++k)
               phiPrime[k+d*numFuncs] *= 2.0 * dxlegendre[idMap[j*numFuncs+k]+(deg+1)*j];
      }
   }

   double h = POW_DOUBLE(10, -6);
   double diffFact = 1.0/(2.0*h);
   double x_backw1[dim];
   double x_forw1[dim];
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;
   Vector basis_polytopic = basis->addData->basis_polytopic;

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim1; ++d)
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         phiPrime[d*numFuncs+k] *= basis_polytopic.id[k];

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
         phiPrime[(d+dim1)*numFuncs+k] *= (phi_forw1.id[k] - phi_backw1.id[k])*diffFact;
   }
}


void CubeSimplexBasisIntegrals(CubeSimplexBasis *basis, Vector v)
{
   int dim2 = basis->params->dims[1];
   int len = basis->integrals.len;
   double *integrals = v.id;

   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1;
   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;
}


void CubeSimplexBasisIntegralsMonomial(CubeSimplexBasis *basis, Vector v)
{
  int numFuncs = basis->numFuncs;
  Vector integralsCube = Vector_init(numFuncs);
  Vector integralsSimplex = Vector_init(numFuncs);
  Vector integrals = v;

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
      Vector_free(basis->addData->basis_polytopic);
      Vector_free(basis->addData->phi_backw1);
      Vector_free(basis->addData->phi_forw1);
      free(basis->addData->idMap);
      RMatrix_free(basis->addData->xFactor[1]);
      free(basis->addData->xPower[1]);
      free(basis->addData);
      basis->addData = NULL;
   }
   if(basis->indices) {
      free(basis->indices);
      basis->indices = NULL;
   }
   free(basis); basis = NULL;
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
   assert(dim1 >= 2  && dim2 >= 2 && deg >= 1);

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


SimplexSimplexBasis* MakeSimplexSimplexCopy(SimplexSimplexBasis *basis)
{
   SimplexSimplexBasis *basisCopy = (SimplexSimplexBasis *)malloc(sizeof(SimplexSimplexBasis));
   basisCopy->params = (SimplexSimplexParams *)malloc(sizeof(SimplexSimplexParams));
   basisCopy->params->deg = basis->deg;
   basisCopy->params->dims[0] = basis->params->dims[0];
   basisCopy->params->dims[1] = basis->params->dims[1];
   basisCopy->deg         = basis->deg;
   basisCopy->dim         = basis->dim;
   basisCopy->numFuncs    = basis->numFuncs;

   BasisInterface interface = SetSimplexSimplexBasisInterface();
   basisCopy->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(basisCopy->interface) = interface;

   int deg = basisCopy->deg;
   int dim = basisCopy->dim;
   int dim1 = basisCopy->params->dims[0];
   int dim2 = basisCopy->params->dims[1];
   int numFuncs = basisCopy->numFuncs;

   basisCopy->indices = (INT_8 *)malloc(numFuncs*dim*sizeof(INT_8));
   memcpy(basisCopy->indices, basis->indices, numFuncs*dim*sizeof(INT_8));

   basisCopy->addData = (AddDataSimplexSimplex *)malloc(sizeof(AddDataSimplexSimplex));
   basisCopy->addData->idMap = (INT_8 *)malloc(dim*numFuncs*sizeof(INT_8));
   memcpy(basisCopy->addData->idMap, basis->addData->idMap, numFuncs*dim*sizeof(INT_8));

   basisCopy->addData->phi_backw1 = Vector_init(numFuncs);
   basisCopy->addData->phi_forw1  = Vector_init(numFuncs);
   basisCopy->addData->basis_polytopic = Vector_init(numFuncs);

   basisCopy->addData->xFactor[0] = RMatrix_init(dim1, numFuncs);
   basisCopy->addData->xPower[0]  = (int *)malloc(dim1*numFuncs*sizeof(int));
   memcpy(basisCopy->addData->xPower[0], basis->addData->xPower[0], dim1*numFuncs*sizeof(int));

   basisCopy->addData->xFactor[1] = RMatrix_init(dim2, numFuncs);
   basisCopy->addData->xPower[1]  = (int *)malloc(dim2*numFuncs*sizeof(int));
   memcpy(basisCopy->addData->xPower[1], basis->addData->xPower[1], dim2*numFuncs*sizeof(int));

   basisCopy->table = Table3dCreate(deg, dim);
   for(int alpha = 1; alpha < basis->table.size[0]; ++alpha) {
      for (int k = 1; k < basis->table.size[1]; ++k) {
         basisCopy->table.id[alpha][k][0] = basis->table.id[alpha][k][0];
         basisCopy->table.id[alpha][k][1] = basis->table.id[alpha][k][1];
         basisCopy->table.id[alpha][k][2] = basis->table.id[alpha][k][2];
      }
   }

   basisCopy->functions   = Vector_init(numFuncs);
   basisCopy->derivatives = Vector_init(numFuncs*dim);
   basisCopy->integrals   = Vector_init(numFuncs);

   return basisCopy;
}


void SimplexSimplexBasisFuncs(SimplexSimplexBasis *basis, const double *x, Vector v)
{
   int numFuncs     = basis->numFuncs;
   Vector polytopic = basis->addData->basis_polytopic;

   VSetToOne(basis->functions);
   double *phi = v.id;

   SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x, polytopic);
   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      phi[k] *= polytopic.id[k];

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, polytopic);
   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      phi[k] *= polytopic.id[k];
}


void SimplexSimplexBasisDer(SimplexSimplexBasis *basis, const double *x, Vector v)
{
   int dim1          = basis->params->dims[0];
   int dim2          = basis->params->dims[1];
   int dim           = basis->dim;
   int numFuncs      = basis->numFuncs;

   double h = POW_DOUBLE(10, -6)*5.0;
   double diffFact = 1.0/(2.0*h);
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;
   Vector basis_polytopic = basis->addData->basis_polytopic;

   double x_backw1[dim];
   double x_forw1[dim];

   VSetToOne(v);
   Vector phiPrime = v;
   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim1; ++d) {
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         phiPrime.id[d*numFuncs+k] *= basis_polytopic.id[k];
   }
   SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim2; ++d) {
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         phiPrime.id[(d+dim1)*numFuncs+k] *= basis_polytopic.id[k];
   }

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
         phiPrime.id[(d+dim1)*numFuncs+k] *= (phi_forw1.id[k] - phi_backw1.id[k]) * diffFact;
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
         phiPrime.id[(d*numFuncs+k)] *= (phi_forw1.id[k] - phi_backw1.id[k]) * diffFact;
   }
}


void SimplexSimplexBasisIntegrals(SimplexSimplexBasis *basis, Vector v)
{
   int dim1 = basis->params->dims[0];
   int dim2 = basis->params->dims[1];
   int len  = basis->integrals.len;
   double *integrals = v.id;

   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1;
   for(int i = 1; i <= dim1; ++i)
      integrals[0] /= i;

   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;
}


void SimplexSimplexBasisIntegralsMonomial(SimplexSimplexBasis *basis, Vector v)
{
  int numFuncs = basis->numFuncs;
  Vector integralsS1 = Vector_init(numFuncs);
  Vector integralsS2 = Vector_init(numFuncs);
  Vector integrals = v;

  IntegralsSimplexPolyhedralMonomialOne((MixedPolytopeBasis *)basis, integralsS1);
  IntegralsSimplexPolyhedralMonomialTwo((MixedPolytopeBasis *)basis, integralsS2);

  #pragma omp simd
  for(int i = 0; i < numFuncs; ++i)
     integrals.id[i] = integralsS1.id[i]*integralsS2.id[i];

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
      Vector_free(basis->addData->basis_polytopic);
      Vector_free(basis->addData->phi_backw1);
      Vector_free(basis->addData->phi_forw1);
      free(basis->addData->idMap);
      RMatrix_free(basis->addData->xFactor[0]);
      RMatrix_free(basis->addData->xFactor[1]);
      free(basis->addData->xPower[0]);
      free(basis->addData->xPower[1]);
      free(basis->addData);
      basis->addData = NULL;
   }
   if(basis->indices) {
      free(basis->indices);
      basis->indices = NULL;
   }
   free(basis); basis = NULL;
}


static void LegendrePoly(int order, double x, double *p)
{
   p[0] = 1.0;
   p[1] = x;

   for(int  k = 1; k < order-1; ++k) {
      double fac1 = (2.0*k+1.0)/(k+1.0);
      double fac2 = k/(k+1.0);
      p[k+1] = fac1*x*p[k] - fac2*p[k-1];
   }
}


static void LegendrePolyAndPrime(int order, double x, double *p, double *dp)
{
   p[0] = 1.0;
   p[1] = x;
   dp[0] = 0.0;
   dp[1] = 1.0;

   for(int  k = 1; k < order-1; ++k) {
      double fac1 = (2.0*k+1.0)/(k+1.0);
      double fac2 = k/(k+1.0);
      p[k+1] = fac1*x*p[k] - fac2*p[k-1];
      dp[k+1] = fac1*(p[k] + x*dp[k]) - fac2*dp[k-1];
   }
}


static Table3d Table3dCreate(int deg, int dim)
{
   Table3d table;
   if(dim >= 3)
      table.size[0] = 2*deg+2+dim-1;
   else if(dim == 2)
      table.size[0] = 2*deg+2;
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
   for(int alpha = 1; alpha < table.size[0]; ++alpha) {
      for (int k = 1; k < table.size[1]; ++k) {
         double fac1Part0 = (alpha+2.0*k+1.0) * (alpha+2.0*k+2.0) * (alpha+2.0*k);
         double fac1Part1 = (alpha+2.0*k+1.0) * (alpha*alpha);                      //excludes multiplication by x
         double fac2 = -2.0 * (alpha+k) * k * (alpha+2.0*k+2.0);
         double fac3 = 1.0 / (2.0*(k+1.0)*(alpha+k+1.0) * (alpha+2.0*k));
         table.id[alpha][k][0] = fac3 * fac1Part0;
         table.id[alpha][k][1] = fac3 * fac1Part1;
         table.id[alpha][k][2] = fac3 * fac2;
      }
   }
}


// asumes order is at least 2 or greater, i.e. polynomial degree is at least 1
static void JacobiPolyWithTable(int order, double x, int alpha, double *p, Table3d table)
{
   p[0] = 1.0;
   p[1] = 0.5*(x-1.0) * (alpha+2.0) + alpha+1.0;
   for(int k = 1; k < order-1; ++k)
      p[k+1] = (table.id[alpha][k][0]*x + table.id[alpha][k][1])*p[k] + table.id[alpha][k][2]*p[k-1];
}


static void IntegralsCubePolyhedralMonomial(MixedPolytopeBasis *basis, Vector v)
{
   int dim          = basis->dim;
   int dimCube      = basis->params->dims[0];
   int numFuncs     = basis->numFuncs;
   Vector integrals = v;
   INT_8 *indices   = basis->indices;

   for(int i = 0; i < numFuncs; ++i)
   { double val = 1.0;
      for(int d = 0; d < dimCube; ++d)
         val = val/(indices[i*dim+d] + 1);
      integrals.id[i] = val;
   }
}


static void IntegralsSimplexPolyhedralMonomialOne(MixedPolytopeBasis *basis, Vector v)
{
   int dim          = basis->dim;
   int numFuncs     = basis->numFuncs;
   int dimSimplex   = basis->params->dims[0];
   INT_8 *indices   = basis->indices;
   Vector integrals = v;

   for(int i = 0; i < numFuncs; ++i)
   {
      double val = 1.0;
      for(int d = 0; d < dimSimplex; ++d) {
         double power = 0;
         for(int r = 0; r < dimSimplex-d; ++r)
            power += (double)indices[dim*i+dimSimplex-r-1];
         val /= (power+dimSimplex-d);
      }
      integrals.id[i] = val;
   }
}


static void IntegralsSimplexPolyhedralMonomialTwo(MixedPolytopeBasis *basis, Vector v)
{
   int dim          = basis->dim;
   int dimSimplex   = basis->params->dims[1];
   int numFuncs     = basis->numFuncs;
   INT_8 *indices   = basis->indices;
   Vector integrals = v;

   for(int i = 0; i < numFuncs; ++i)
   {
      double val = 1.0;
      for(int d = 0; d < dimSimplex; ++d) {
         double power = 0;
         for(int r = 0; r < dimSimplex-d; ++r)
            power += (double)indices[dim*(i+1)-r-1];
         val /= (power+dimSimplex-d);
      }
      integrals.id[i] = val;
   }
}


void SimplexFuncsPolytopicOne(MixedPolytopeBasis *basis, const double *x, Vector v)
{
   int deg        = basis->deg;
   int dim1       = basis->params->dims[0];
   int numFuncs   = basis->numFuncs;
   INT_8 *idMap   = basis->addData->idMap;
   double *phi    = v.id;

   double jacobi[SQUARE(deg+1)*(dim1-1)];
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
      int nextDim = (deg+1)*(d-1)*deg;
      for(int j = 0; j < deg+1; ++j) {
         double alpha = 2*j+d;
         JacobiPolyWithTable(deg+1, jCoord[d-1], alpha, &jacobi[nextDim+(deg+1)*j], basis->table);
      }
   }

   int *xPow = basis->addData->xPower[0];
   RMatrix xFactor = basis->addData->xFactor[0];
   for(int d = 1; d < dim1; ++d)
      for(int k = 0; k < numFuncs; ++k)
            xFactor.rid[d][k] = DoubleIntPower(xCoord[d-1], xPow[id2(d, k, numFuncs)]);

   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      phi[k] = legendre[idMap[k]];

   for(int d = 1; d < dim1; ++d) {
      int nextDim = (deg+1)*(d-1)*deg;
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         phi[k] *= jacobi[idMap[d*numFuncs+k] + nextDim + (deg+1)*xPow[d*numFuncs+k]] * xFactor.rid[d][k];
   }
}


void SimplexFuncsPolytopicTwo(MixedPolytopeBasis *basis, const double *x, Vector v)
{
   int deg        = basis->deg;
   int dim        = basis->dim;
   int dim1       = basis->params->dims[0];
   int dim2       = basis->params->dims[1];
   int dimTwo     = dim1+dim2;
   int numFuncs   = basis->numFuncs;
   double *phi    = v.id;
   INT_8 *idMap   = basis->addData->idMap;

   double jacobi[SQUARE(deg+1)*(dim2-1)];
   double legendre[deg+1];
   LegendrePoly(deg+1, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], legendre);

   double xCoord[dim2-1];
   for(int d = 0; d < dim2-2; ++d) xCoord[d] = x[dimTwo-d-2]/x[dimTwo-d-3];
   xCoord[dim2-2] = x[dim1];

   double jCoord[dim];
   for(int d = 1; d < dim2; ++d)
      jCoord[d-1] = 1.0-2.0*x[dimTwo-d-1]/x[dimTwo-d-2];
   jCoord[dim2-2] = 1.0-2.0*x[dim1];

   for(int d = 1; d < dim2; ++d) {
      int nextDim = (deg+1)*(d-1)*deg;
      for(int j = 0; j < deg+1; ++j) {
         int nextAlpha = (deg+1)*j;
         JacobiPolyWithTable(deg+1, jCoord[d-1], 2*j+d, &jacobi[nextDim+nextAlpha], basis->table);
      }
   }

   RMatrix xFactor = basis->addData->xFactor[1];
   int *xPower = basis->addData->xPower[1];
   for(int d = 1; d < dim2; ++d)
      for(int k = 0; k < numFuncs; ++k)
            xFactor.rid[d][k] = DoubleIntPower(xCoord[d-1], xPower[id2(d, k, numFuncs)]);

   #pragma omp simd
   for(int k = 0; k < numFuncs; ++k)
      phi[k] = legendre[idMap[(dim1)*numFuncs+k]];

   for(int d = 1; d < dim2; ++d) {
      int nextDim = (deg+1)*(d-1)*deg;
      #pragma omp simd
      for(int k = 0; k < numFuncs; ++k)
         phi[k] *= jacobi[idMap[(d+dim1)*numFuncs+k] + nextDim + (deg+1)*xPower[d*numFuncs+k]] * xFactor.rid[d][k];
   }
}




double orthogonal_cube_basis_test(int deg, int dim)
{
   int numFuncs = BasisSize(deg, dim);
   int n = ceil( (deg+1)/2.0 )+1; // make n large enough to get exact quadrature
   int dims_1D[1] = {1};
   quadrature *quad_1D = quadrature_init_basic(n, 1, dims_1D, deg, INTERVAL);
   Jacobi(quad_1D->num_nodes, 0.0, 0.0, quad_1D->x, quad_1D->w);

   int N = POW_INT(n, dim);
   quadrature *quad_C = quadrature_init_full(N, dim, &dim, deg, CUBE);
   GeneralizedNodesTensor(quad_1D, quad_C);
   GeneralizedWeightsTensor(quad_1D, quad_C);
   Vector functions = quad_C->basis->functions;

   double *quad_integrals = (double *)calloc(numFuncs, sizeof(double));
   Basis *basis = quad_C->basis;
   BasisIntegrals(basis, basis->integrals);

   for(int i = 0; i < N; ++i) {
      BasisFuncs(basis, &quad_C->x[dim*i], functions);
      for(int j = 0; j < numFuncs; ++j)
         quad_integrals[j] += functions.id[j]*quad_C->w[i];
   }

   double max_res = fabs(quad_integrals[0] - basis->integrals.id[0]);
   for(int j = 1; j < numFuncs; ++j) {
      double res = fabs(quad_integrals[j]-basis->integrals.id[j]);
      max_res = MAX(max_res, res);
   }
   printf("\nTesting orthogonality of cube basis functions. Maximum error of basis functions = %.16e\n\n", max_res);

   free(quad_integrals); quad_integrals = NULL;
   quadrature_free(quad_1D);
   quadrature_free(quad_C);

   return max_res;
}


double orthogonal_simplex_basis_test(int deg, int dim)
{
   int numFuncs = BasisSize(deg, dim);
   int n = ceil( (deg+dim)/2.0 )+2; // make n large enough to get exact quadrature
   int dims_1D[1] = {1};
   quadrature *quad_1D = quadrature_init_basic(n, 1, dims_1D, deg, INTERVAL);
   Jacobi(quad_1D->num_nodes, 0.0, 0.0, quad_1D->x, quad_1D->w);

   int N = POW_INT(n, dim);
   quadrature *quad_S = quadrature_init_full(N, dim, &dim, deg, SIMPLEX);
   GeneralizedNodesTensor(quad_1D, quad_S);
   GeneralizedWeightsTensor(quad_1D, quad_S);
   GeneralDuffy(quad_S);
   Vector functions = quad_S->basis->functions;

   double *quad_integrals = (double *)calloc(numFuncs, sizeof(double));
   Basis *basis = quad_S->basis;
   BasisIntegrals(basis, basis->integrals);

   for(int i = 0; i < N; ++i) {
      BasisFuncs(basis, &quad_S->x[dim*i], functions);
      for(int j = 0; j < numFuncs; ++j)
         quad_integrals[j] += functions.id[j]*quad_S->w[i];
   }

   double max_res = fabs(quad_integrals[0] - basis->integrals.id[0]);
   for(int j = 1; j < numFuncs; ++j) {
      double res = fabs(quad_integrals[j]-basis->integrals.id[j]);
      max_res = MAX(max_res, res);
   }
   printf("\nTesting orthogonality of simplex basis functions. Maximum error of basis functions = %.16e\n\n", max_res);

   free(quad_integrals); quad_integrals = NULL;
   quadrature_free(quad_1D);
   quadrature_free(quad_S);
   return max_res;
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

