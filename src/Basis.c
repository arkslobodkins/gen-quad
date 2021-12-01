#include "Basis.h"
#include "BasisFunctions.h"
#include "BasisIndices.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

BasisInterface SetCubeBasisInterface()
{
   BasisInterface cubeInterface;
   cubeInterface.basisInit        = (BasisInitPtr)&CubeBasisInit;
   cubeInterface.computeFuncs     = (BasisFuncsPtr)&ComputeCubeBasisFuncs;
   cubeInterface.computeDer       = (BasisDerPtr)&ComputeCubeBasisDer;
   cubeInterface.computeIntegrals = (BasisIntegralsPtr)&CubeBasisIntegrals;
   cubeInterface.basisFree        = (BasisFreePtr)&CubeBasisFree;
   return cubeInterface;
}

BasisInterface SetCubeSimplexBasisInterface()
{
   BasisInterface csInterface;
   csInterface.basisInit        = (BasisInitPtr)&CubeSimplexBasisInit;
   csInterface.computeFuncs     = (BasisFuncsPtr)&CubeSimplexBasisFuncs;
   csInterface.computeDer       = (BasisDerPtr)&CubeSimplexBasisDer;
   csInterface.computeIntegrals = (BasisIntegralsPtr)&CubeSimplexBasisIntegrals;
   csInterface.basisFree        = (BasisFreePtr)&CubeSimplexBasisFree;
   return csInterface;
}

Basis* BasisInit(void *params, BasisInterface *interface)
{
   Basis *basis = interface->basisInit(params);
   basis->interface = interface;

   int numFuncs = basis->numFuncs;
   int dim = basis->dim;
   basis->basis_indices   = (INT_8 *) malloc(numFuncs*dim*sizeof(INT_8));
   basis->basis_funcs     = Vector_init(numFuncs);
   basis->basis_der       = Vector_init(numFuncs*dim);
   basis->basis_integrals = Vector_init(numFuncs);
   return basis;
}

void BasisFree(Basis *basis)
{
   if(basis->basis_indices != NULL) {
      free(basis->basis_indices);
      basis->basis_indices = NULL;
   }
   Vector_free(basis->basis_integrals);
   Vector_free(basis->basis_der);
   Vector_free(basis->basis_funcs);
   free(basis->basis_indices);

   basis->interface->basisFree(basis);
}

void BasisFuncs(Basis *basis, double *x)
{
   basis->interface->computeFuncs(basis, x);
}

void BasisDer(Basis *basis, double *x)
{
   basis->interface->computeDer(basis, x);
}

void Integrals(Basis *basis)
{
   basis->interface->computeIntegrals(basis);
}

CubeBasis* CubeBasisInit(CubeParams *cubeParams)
{
   int deg = cubeParams->deg;
   int dim = cubeParams->dim;
   assert(dim >= 0 && deg >= 0);

   CubeBasis *cubeBasis = (CubeBasis *)malloc(sizeof(CubeBasis));
   cubeBasis->params = (CubeParams *)malloc(sizeof(CubeParams));
   cubeBasis->params->deg = deg;
   cubeBasis->params->dim = dim;
   cubeBasis->deg = deg;
   cubeBasis->dim = dim;
   cubeBasis->numFuncs = BasisSize(deg, dim);
   return cubeBasis;
}

void ComputeCubeBasisFuncs(CubeBasis *cubeBasis, double *x)
{
   int deg = cubeBasis->deg;
   int dim = cubeBasis->dim;
   int numFuncs = cubeBasis->numFuncs;
   int degPlus1 = deg+1;
   INT_8 *basisId = cubeBasis->basis_indices;
   double *phi = cubeBasis->basis_funcs.id;

   double legendre[degPlus1*dim];
   double dxlegendre[degPlus1*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;
   for(int k = 0; k < numFuncs; ++k)
      for(int d = 0; d < dim; ++d)
         phi[k] *= legendre[basisId[k*dim+d] + degPlus1*d];
}

void ComputeCubeBasisDer(CubeBasis *cubeBasis, double *x)
{
   int deg = cubeBasis->deg;
   int dim = cubeBasis->dim;
   int numFuncs = cubeBasis->numFuncs;
   int degPlus1 = deg+1;
   INT_8 *basisId = cubeBasis->basis_indices;
   double *phiPrime = cubeBasis->basis_der.id;

   double legendre[degPlus1*dim];
   double dxlegendre[degPlus1*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < numFuncs*dim; ++k) phiPrime[k] = 1.0;

   for(int d = 0; d < dim; ++d)
   {
      // dimension < d
      for(int j = 0; j < d; ++j)
      {
         double *legendre_ptr = &legendre[degPlus1*j];
         for(int k = 0; k < numFuncs; ++k)
            phiPrime[d*numFuncs+k] *= legendre_ptr[basisId[dim*k+j]];
      }
      // dimension = d
      double *dxlegendre_ptr = &dxlegendre[degPlus1*d];
      for(int k = 0; k < numFuncs; ++k)
         phiPrime[d*numFuncs+k] *= 2.0 * dxlegendre_ptr[basisId[dim*k+d]];

      // dimension > d
      for(int j = d+1; j < dim; ++j)
      {
         double *legendre_ptr = &legendre[degPlus1*j];
         for(int k = 0; k < numFuncs; ++k)
            phiPrime[d*numFuncs+k] *= legendre_ptr[basisId[dim*k+j]];
      }
   }
}

void CubeBasisIntegrals(CubeBasis *cubeBasis)
{
   int len = cubeBasis->basis_integrals.len;
   double *integrals = cubeBasis->basis_integrals.id;
   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1.0;
}

void CubeBasisFree(CubeBasis *cubeBasis)
{
   free(cubeBasis->params);
   free(cubeBasis);
}

CubeSimplexBasis* CubeSimplexBasisInit(CubeSimplexParams *csParams)
{
   int deg = csParams->deg;
   int dim = csParams->dims[0]+csParams->dims[1];

   CubeSimplexBasis *csBasis = (CubeSimplexBasis *)malloc(sizeof(CubeSimplexBasis));
   csBasis->params = (CubeSimplexParams *)malloc(sizeof(CubeSimplexParams));

   csBasis->params->deg      = deg;
   csBasis->params->dims[0]  = csParams->dims[0];
   csBasis->params->dims[1]  = csParams->dims[1];
   csBasis->deg              = deg;
   csBasis->dim              = dim;
   csBasis->numFuncs = BasisSize(deg, dim);
   return csBasis;
}

void CubeSimplexBasisFree(CubeSimplexBasis *csBasis)
{
   free(csBasis->params);
   free(csBasis);
}

void CubeSimplexBasisFuncs(CubeSimplexBasis *csBasis, double *x)
{

}

void CubeSimplexBasisDer(CubeSimplexBasis *csBasis, double *x)
{

}

void CubeSimplexBasisIntegrals(CubeSimplexBasis *csBasis)
{
   int dim2 = csBasis->params->dims[1];
   int len = csBasis->basis_integrals.len;
   double *integrals = csBasis->basis_integrals.id;

   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1;
   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;
}

void ComputeBasisIndices(Basis *basis)
{
   int dim = basis->dim;
   int deg = basis->deg;
   INT_8* f = basis->basis_indices;

   if(dim == 1) {
      for(int i = 0; i <= deg; ++i) f[i] = i;
      return;
   }

   int counter;
   // compute basis indices using nested recursion if dimension >= 2
   for(int j = 2; counter = 0, j <= dim; ++j)
   {
      for(int i = 0; i <= deg; ++i)
      {
         int size = BasisSize(deg-i, j-1);
         int dimxsize = dim*size;
         INT_8 recursiveF[size*(j-1)];
         BasisIndices(deg-i, j-1, recursiveF);
         if(j == dim)
         {
            for(int k = 0; k < size; ++k)
            {
               int kxdim = k*dim;
               int kxdim_minus1 = k*(dim-1);
               for(int d = 0; d < dim-1; ++d)
                  f[counter+kxdim+d] = recursiveF[kxdim_minus1+d];
               f[counter+kxdim+j-1] = i;
            }
            counter += dimxsize;
         }
      }
   }
}

void TestBasis()
{
   CubeParams cubeParams = {2,3};
   double x[3];
   CubeParams *params = &cubeParams;
   BasisInterface interface = SetCubeBasisInterface();
   Basis *cube = BasisInit((void *)params, &interface);
   ComputeBasisIndices(cube);
   BasisFuncs(cube, x);
   BasisDer(cube, x);
   printf("basisSize = %i\n", cube->numFuncs);
   BasisFree(cube);

   CubeSimplexParams csParams;
   csParams.deg = 2; csParams.dims[0] = 1; csParams.dims[1] = 2;
   CubeSimplexParams *paramsPtr = &csParams;
   double y[3];
   BasisInterface csinterface = SetCubeSimplexBasisInterface();
   Basis *cs = BasisInit((void *)paramsPtr, &csinterface);

   BasisFree(cs);
}

