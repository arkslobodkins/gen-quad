#include "Basis.h"
#include "BasisFunctions.h"
#include "BasisIndices.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


Basis* BasisInit(void *params, BasisInterface *interface)
{
   Basis *basis = interface->basisInit(params);
   basis->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(basis->interface) = *interface;

   int numFuncs = basis->numFuncs;
   int dim = basis->dim;
   basis->basis_indices   = (INT_8 *) malloc(numFuncs*dim*sizeof(INT_8));
   ComputeBasisIndices(basis);
   basis->basis_funcs     = Vector_init(numFuncs);
   basis->basis_der       = Vector_init(numFuncs*dim);
   basis->basis_integrals = Vector_init(numFuncs);
   return basis;
}

void BasisFuncs(Basis *basis, const double *x, Vector v)
{
   basis->interface->computeFuncs(basis, x, v);
}

void BasisDer(Basis *basis, const double *x)
{
   basis->interface->computeDer(basis, x);
}

void Integrals(Basis *basis)
{
   basis->interface->computeIntegrals(basis);
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
   free(basis->interface);
   free(basis);
}


/////////////////////////////////////////////////////////////////////////////
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

CubeBasis* CubeBasisInit(CubeParams *cubeParams)
{
   int deg = cubeParams->deg;
   int dim = cubeParams->dim;
   assert(dim >= 2 && deg >= 1);

   CubeBasis *basis = (CubeBasis *)malloc(sizeof(CubeBasis));
   basis->params = (CubeParams *)malloc(sizeof(CubeParams));
   basis->params->deg = deg;
   basis->params->dim = dim;
   basis->deg         = deg;
   basis->dim         = dim;
   basis->numFuncs    = BasisSize(deg, dim);
   basis->addData     = NULL;
   return basis;
}

void ComputeCubeBasisFuncs(CubeBasis *basis, const double *x, Vector v)
{
   int deg        = basis->deg;
   int dim        = basis->dim;
   int numFuncs   = basis->numFuncs;
   INT_8 *basisId = basis->basis_indices;
   double *phi    = v.id;

   double legendre[(deg+1)*dim];
   double dxlegendre[(deg+1)*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePoly(deg+1, 2*x[d]-1, &legendre[d*(deg+1)], &dxlegendre[d*(deg+1)]);

   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;
   for(int k = 0; k < numFuncs; ++k)
      for(int d = 0; d < dim; ++d)
         phi[k] *= legendre[basisId[k*dim+d] + d*(deg+1)];
}

void ComputeCubeBasisDer(CubeBasis *basis, const double *x)
{
   int deg          = basis->deg;
   int dim          = basis->dim;
   int numFuncs     = basis->numFuncs;
   int degPlus1     = deg+1;
   INT_8 *basisId   = basis->basis_indices;
   double *phiPrime = basis->basis_der.id;

   double legendre[degPlus1*dim];
   double dxlegendre[degPlus1*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < numFuncs*dim; ++k) phiPrime[k] = 1.0;

   for(int d = 0; d < dim; ++d)
   {
      // dimension < d
      for(int j = 0; j < d; ++j) {
         double *legendre_ptr = &legendre[degPlus1*j];
         for(int k = 0; k < numFuncs; ++k)
            phiPrime[d*numFuncs+k] *= legendre_ptr[basisId[dim*k+j]];
      }
      // dimension = d
      double *dxlegendre_ptr = &dxlegendre[degPlus1*d];
      for(int k = 0; k < numFuncs; ++k)
         phiPrime[d*numFuncs+k] *= 2.0 * dxlegendre_ptr[basisId[dim*k+d]];

      // dimension > d
      for(int j = d+1; j < dim; ++j) {
         double *legendre_ptr = &legendre[degPlus1*j];
         for(int k = 0; k < numFuncs; ++k)
            phiPrime[d*numFuncs+k] *= legendre_ptr[basisId[dim*k+j]];
      }
   }
}

void CubeBasisIntegrals(CubeBasis *basis)
{
   int len = basis->basis_integrals.len;
   double *integrals = basis->basis_integrals.id;
   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1.0;
}

void CubeBasisFree(CubeBasis *basis)
{
   free(basis->params);
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetSimplexBasisInterface()
{
   BasisInterface simplexInterface;
   simplexInterface.basisInit        = (BasisInitPtr)&SimplexBasisInit;
   simplexInterface.computeFuncs     = (BasisFuncsPtr)&SimplexBasisFuncs;
   simplexInterface.computeDer       = (BasisDerPtr)&SimplexBasisDer;
   simplexInterface.computeIntegrals = (BasisIntegralsPtr)&SimplexBasisIntegrals;
   simplexInterface.basisFree        = (BasisFreePtr)&SimplexBasisFree;
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

   basis->addData = (AddDataSimplex *)malloc(sizeof(AddDataSimplex));
   int numFuncs = basis->numFuncs;
   AddDataSimplex *addData = basis->addData;
   addData->phi_backw2 = Vector_init(numFuncs);
   addData->phi_backw1 = Vector_init(numFuncs);
   addData->phi_forw1  = Vector_init(numFuncs);
   addData->phi_forw2  = Vector_init(numFuncs);
   return basis;
}

void SimplexBasisFuncs(SimplexBasis *basis, const double *x, Vector v)
{
   int i, j, k, d;
   int deg        = basis->deg;
   int dim        = basis->dim;
   int numFuncs   = basis->numFuncs;
   int degPlus1   = deg+1;
   double legendre[degPlus1];
   double dxlegendre[degPlus1];
   double jacobi[SQUARE(degPlus1)*(dim-1)];

   INT_8 *basisId = basis->basis_indices;
   double *phi    = v.id;
   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;

   LegendrePoly(degPlus1, (2.0*x[dim-1]-x[dim-2])/x[dim-2], legendre, dxlegendre);
   double xCoord[dim-1];
   for(d = 0; d < dim-1; ++d) xCoord[d] = 1.0;
   for(d = 0; d < dim-2; ++d) xCoord[d] = x[dim-d-2]/x[dim-d-3];
   xCoord[dim-2] = x[0];

   for(d = 1; d < dim; ++d)
   {
      int dimCur = d-1;
      double jCoord = 1.0;
      if(d < dim-1) {
         jCoord /= x[dim-d-2];
         jCoord = 1.0-2.0*x[dim-d-1] * jCoord;
      }
      else if(d == dim-1) {
         jCoord *= x[0];
         jCoord = 1.0-2.0*jCoord;
      }
      for(j = 0; j < degPlus1; ++j) {
         double alpha = 2*j+d;
         JacobiPoly(degPlus1, jCoord, alpha, 0.0, &jacobi[degPlus1*dimCur*deg+degPlus1*j]);
      }
   }

   for(d = 1; d < dim; ++d)
   {
      for(k = 0; k < numFuncs; ++k) {
         const INT_8 *index = &basisId[k*dim];
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = DoubleIntPower(xCoord[d-1], xPower);
         double phiTemp = jacobi[*(index+d) + degPlus1*(d-1)*deg + degPlus1*xPower] * xFactor;
         phi[k] *= phiTemp;
      }
   }

   for(k = 0; k < numFuncs; ++k) {
      const INT_8 *index = &basisId[k*dim];
      phi[k] *= legendre[*index];
   }
}

void SimplexBasisDer(SimplexBasis *basis, const double *x)
{
   int dim = basis->params->dim;
   int numFuncs = basis->numFuncs;
   double h = POW_DOUBLE(10, -5)*5.0;

   double x_backw1[dim];
   double x_backw2[dim];
   double x_forw1[dim];
   double x_forw2[dim];

   double *basisDer = basis->basis_der.id;
   Vector phi_backw2 = basis->addData->phi_backw2;
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1 = basis->addData->phi_forw1;
   Vector phi_forw2 = basis->addData->phi_forw2;

   for(int d = 0; d < dim; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t)
      {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t]  = x[d_t];
         x_forw2[d_t]  = x[d_t];
      }
      x_backw2[d] = x_backw2[d] - 2.0*h;
      x_backw1[d] = x_backw1[d] - h;
      x_forw1[d]  = x_forw1[d] + h;
      x_forw2[d]  = x_forw2[d] + 2.0*h;

      SimplexBasisFuncs(basis, x_backw2, phi_backw2);
      SimplexBasisFuncs(basis, x_backw1, phi_backw1);
      SimplexBasisFuncs(basis, x_forw1, phi_forw1);
      SimplexBasisFuncs(basis, x_forw2, phi_forw2);
      for(int k = 0; k < numFuncs; ++k)
         basisDer[d*numFuncs+k] = ( 1.0/12.0*phi_backw2.id[k] - 2.0/3.0*phi_backw1.id[k] +
                                    2.0/3.0*phi_forw1.id[k]-1.0/12.0*phi_forw2.id[k] ) / h;
   }
}

void SimplexBasisIntegrals(SimplexBasis *basis)
{
   int dim = basis->dim;
   int len = basis->basis_integrals.len;
   double *integrals = basis->basis_integrals.id;

   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1;
   for(int i = 1; i <= dim; ++i)
      integrals[0] /= i;

}

void SimplexBasisFree(SimplexBasis *basis)
{
   free(basis->params);
   Vector_free(basis->addData->phi_backw1);
   Vector_free(basis->addData->phi_backw2);
   Vector_free(basis->addData->phi_forw1);
   Vector_free(basis->addData->phi_forw2);
   free(basis->addData);
}


/////////////////////////////////////////////////////////////////////////////
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

CubeSimplexBasis* CubeSimplexBasisInit(CubeSimplexParams *csParams)
{
   int deg = csParams->deg;
   int dim1 = csParams->dims[0];
   int dim2 = csParams->dims[1];
   int dim = dim1+dim2;
   assert(dim1 >= 1  && dim2 >= 2 && deg >= 1);

   CubeSimplexBasis *basis = (CubeSimplexBasis *)malloc(sizeof(CubeSimplexBasis));
   basis->params = (CubeSimplexParams *)malloc(sizeof(CubeSimplexParams));

   basis->params->deg         = deg;
   basis->params->dims[0]     = dim1;
   basis->params->dims[1]     = dim2;
   basis->deg                 = deg;
   basis->dim                 = dim;
   basis->numFuncs            = BasisSize(deg, dim);

   basis->addData = (AddDataCubeSimplex *)malloc(sizeof(AddDataCubeSimplex));
   int numFuncs = basis->numFuncs;
   AddDataCubeSimplex *addData = basis->addData;
   addData->phi_backw2 = Vector_init(numFuncs);
   addData->phi_backw1 = Vector_init(numFuncs);
   addData->phi_forw1  = Vector_init(numFuncs);
   addData->phi_forw2  = Vector_init(numFuncs);
   addData->basis_polytopic = Vector_init(numFuncs);

   return basis;
}

void CubeSimplexBasisFuncs(CubeSimplexBasis *basis, const double *x, Vector v)
{
   int deg        = basis->deg;
   int dim        = basis->dim;
   int dim1       = basis->params->dims[0];
   int numFuncs   = basis->numFuncs;
   int degPlus1   = deg+1;
   double *phi    = v.id;
   INT_8 *basisId = basis->basis_indices;

   double legendre[degPlus1*dim1];
   double dxlegendre[degPlus1*dim1];

   for(int d = 0; d < dim1; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;
   for(int k = 0; k < numFuncs; ++k)
      for(int d = 0; d < dim1; ++d)
         phi[k] *= legendre[basisId[k*dim+d]+degPlus1*d];

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis->addData->basis_polytopic);
   for(int k = 0; k < numFuncs; ++k)
      phi[k] *= basis->addData->basis_polytopic.id[k];
}

void CubeSimplexBasisDer(CubeSimplexBasis *basis, const double *x)
{
   int deg          = basis->deg;
   int dim          = basis->dim;
   int dim1         = basis->params->dims[0];
   int dim2         = basis->params->dims[1];
   int numFuncs     = basis->numFuncs;
   int degPlus1     = deg+1;
   double *phiPrime = basis->basis_der.id;
   INT_8 *basisId = basis->basis_indices;

   double legendre[degPlus1*dim1];
   double dxlegendre[degPlus1*dim1];

   for(int d = 0; d < dim1; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

   for(int k = 0; k < numFuncs*dim; ++k) phiPrime[k] = 1.0;

   for(int d = 0; d < dim; ++d)
   {
      for(int j = 0; j < dim1; ++j) {
         if(j != d)
            for(int k = 0; k < numFuncs; ++k)
               phiPrime[k+d*numFuncs] *= legendre[basisId[k*dim+j]+degPlus1*j];
         if(j == d)
            for(int k = 0; k < numFuncs; ++k)
               phiPrime[k+d*numFuncs] *= 2.0 * dxlegendre[basisId[k*dim+j]+degPlus1*j];
      }
   }

   double h = POW_DOUBLE(10, -5)*5.0;
   double x_backw1[dim];
   double x_backw2[dim];
   double x_forw1[dim];
   double x_forw2[dim];
   Vector phi_backw2 = basis->addData->phi_backw2;
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;
   Vector phi_forw2  = basis->addData->phi_forw2;
   Vector basis_polytopic = basis->addData->basis_polytopic;

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim1; ++d) {
      for(int k = 0; k < numFuncs; ++k)
         phiPrime[d*numFuncs+k] *= basis_polytopic.id[k];
   }

   for(int d = 0; d < dim2; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t)
      {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t]  = x[d_t];
         x_forw2[d_t]  = x[d_t];
      }
      x_backw2[d+dim1] = x_backw2[d+dim1] - 2.0*h;
      x_backw1[d+dim1] = x_backw1[d+dim1] - h;
      x_forw1[d+dim1]  = x_forw1[d+dim1] + h;
      x_forw2[d+dim1]  = x_forw2[d+dim1] + 2.0*h;

      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_backw1, phi_backw1);
      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_backw2, phi_backw2);
      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_forw1, phi_forw1);
      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_forw2, phi_forw2);

      for(int k = 0; k < numFuncs; ++k)
         phiPrime[(d+dim1)*numFuncs+k] *= ( 1.0/12.0*phi_backw2.id[k] - 2.0/3.0*phi_backw1.id[k] +
                                            2.0/3.0*phi_forw1.id[k] - 1.0/12.0*phi_forw2.id[k] ) / h;
   }
}

void CubeSimplexBasisIntegrals(CubeSimplexBasis *basis)
{
   int dim2 = basis->params->dims[1];
   int len = basis->basis_integrals.len;
   double *integrals = basis->basis_integrals.id;

   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1;
   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;
}

void CubeSimplexBasisFree(CubeSimplexBasis *basis)
{
   free(basis->params);
   Vector_free(basis->addData->basis_polytopic);
   Vector_free(basis->addData->phi_backw1);
   Vector_free(basis->addData->phi_backw2);
   Vector_free(basis->addData->phi_forw1);
   Vector_free(basis->addData->phi_forw2);
   free(basis->addData);
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetSimplexSimplexBasisInterface()
{
   BasisInterface ssInterface;
   ssInterface.basisInit        = (BasisInitPtr)&SimplexSimplexBasisInit;
   ssInterface.computeFuncs     = (BasisFuncsPtr)&SimplexSimplexBasisFuncs;
   ssInterface.computeDer       = (BasisDerPtr)&SimplexSimplexBasisDer;
   ssInterface.computeIntegrals = (BasisIntegralsPtr)&SimplexSimplexBasisIntegrals;
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

   basis->addData = (AddDataSimplexSimplex *)malloc(sizeof(AddDataSimplexSimplex));
   int numFuncs = basis->numFuncs;
   AddDataSimplexSimplex *addData = basis->addData;
   addData->phi_backw2 = Vector_init(numFuncs);
   addData->phi_backw1 = Vector_init(numFuncs);
   addData->phi_forw1  = Vector_init(numFuncs);
   addData->phi_forw2  = Vector_init(numFuncs);
   addData->basis_polytopic = Vector_init(numFuncs);

   return basis;
}

void SimplexSimplexBasisFuncs(SimplexSimplexBasis *basis, const double *x, Vector v)
{
   int numFuncs     = basis->numFuncs;
   double *phi      = basis->basis_funcs.id;
   Vector polytopic = basis->addData->basis_polytopic;

   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;

   SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x, polytopic);
   for(int k = 0; k < numFuncs; ++k)
      phi[k] *= polytopic.id[k];

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, polytopic);
   for(int k = 0; k < numFuncs; ++k)
      phi[k] *= polytopic.id[k];
}

void SimplexSimplexBasisDer(SimplexSimplexBasis *basis, const double *x)
{
   int dim1          = basis->params->dims[0];
   int dim2          = basis->params->dims[1];
   int dim           = basis->dim;
   int numFuncs      = basis->numFuncs;
   Vector phiPrime   = basis->basis_der;

   double h = POW_DOUBLE(10, -5)*5.0;
   Vector phi_backw2 = basis->addData->phi_backw2;
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;
   Vector phi_forw2  = basis->addData->phi_forw2;
   Vector basis_polytopic = basis->addData->basis_polytopic;
   double x_backw1[dim];
   double x_backw2[dim];
   double x_forw1[dim];
   double x_forw2[dim];

   for(int k = 0; k < phiPrime.len; ++k) phiPrime.id[k] = 1.0;
   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim1; ++d) {
      for(int k = 0; k < numFuncs; ++k)
         phiPrime.id[d*numFuncs+k] *= basis_polytopic.id[k];
   }
   SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x, basis_polytopic);
   for(int d = 0; d < dim2; ++d) {
      for(int k = 0; k < numFuncs; ++k)
         phiPrime.id[(d+dim1)*numFuncs+k] *= basis_polytopic.id[k];
   }

   for(int d = 0; d < dim2; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t) {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t]  = x[d_t];
         x_forw2[d_t]  = x[d_t];
      }
      x_backw2[d+dim1] = x_backw2[d+dim1] - 2.0*h;
      x_backw1[d+dim1] = x_backw1[d+dim1] - h;
      x_forw1[d+dim1]  = x_forw1[d+dim1] + h;
      x_forw2[d+dim1]  = x_forw2[d+dim1] + 2.0*h;

      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_backw1, phi_backw1);
      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_backw2, phi_backw2);
      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_forw1, phi_forw1);
      SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x_forw2, phi_forw2);

      for(int k = 0; k < numFuncs; ++k)
         phiPrime.id[(d+dim1)*numFuncs+k] *= ( 1.0/12.0*phi_backw2.id[k] - 2.0/3.0*phi_backw1.id[k] +
                                               2.0/3.0*phi_forw1.id[k] - 1.0/12.0*phi_forw2.id[k] ) / h;
   }

   for(int d = 0; d < dim1; ++d)
   {
      for(int d_t = 0; d_t < dim; ++d_t) {
         x_backw1[d_t] = x[d_t];
         x_backw2[d_t] = x[d_t];
         x_forw1[d_t] = x[d_t];
         x_forw2[d_t] = x[d_t];
      }
      x_backw2[d] = x_backw2[d] - 2.0*h;
      x_backw1[d] = x_backw1[d] - h;
      x_forw1[d] = x_forw1[d] + h;
      x_forw2[d] = x_forw2[d] + 2.0*h;

      SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x_backw1, phi_backw1);
      SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x_backw2, phi_backw2);
      SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x_forw1, phi_forw1);
      SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x_forw2, phi_forw2);

      for(int k = 0; k < numFuncs; ++k)
         phiPrime.id[d*numFuncs+k] *= ( 1.0/12.0*phi_backw2.id[k] - 2.0/3.0*phi_backw1.id[k] +
                                        2.0/3.0*phi_forw1.id[k] -1.0/12.0*phi_forw2.id[k] ) / h;
   }
}

void SimplexSimplexBasisIntegrals(SimplexSimplexBasis *basis)
{
   int dim1 = basis->params->dims[0];
   int dim2 = basis->params->dims[1];
   int len  = basis->basis_integrals.len;
   double *integrals = basis->basis_integrals.id;

   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1;
   for(int i = 1; i <= dim1; ++i)
      integrals[0] /= i;

   for(int i = 1; i <= dim2; ++i)
      integrals[0] /= i;
}

void SimplexSimplexBasisFree(SimplexSimplexBasis *basis)
{
   free(basis->params);
   Vector_free(basis->addData->basis_polytopic);
   Vector_free(basis->addData->phi_backw1);
   Vector_free(basis->addData->phi_backw2);
   Vector_free(basis->addData->phi_forw1);
   Vector_free(basis->addData->phi_forw2);
   free(basis->addData);
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
//   CubeParams cubeParams = {3,3};
//   double x[3];
//   x[0] = 0.2; x[1] = 0.3; x[2] = 0.4;
//   CubeParams *params = &cubeParams;
//   BasisInterface interface = SetCubeBasisInterface();
//   Basis *cube = BasisInit((void *)params, &interface);
//   ComputeBasisIndices(cube);
//   BasisFuncs(cube, x);
//   BasisDer(cube, x);
//
//   Basis *cube2 = BasisInit((void *)params, &interface);
//   ComputeBasisIndices(cube2);
//   BasisFuncs(cube2, x);
//   CubeBasisFuncs(&cubeParams.dim, cubeParams.deg, cube2->basis_indices, x, cube2->basis_funcs.id);
//   CubeBasisFuncsDer(&cubeParams.dim, cubeParams.deg, cube2->basis_indices, x, cube2->basis_der.id);
//   Vector vc = Vector_init(cube2->basis_funcs.len);
//   VectorAddScale(1.0, cube->basis_funcs, -1.0, cube2->basis_funcs, vc);
////   VPrint(vc);
//
////   VPrint(cube2->basis_der);
////   VPrint(cube->basis_der);
//   BasisFree(cube);
//
//   CubeSimplexParams csParams;
//   csParams.deg = 3; csParams.dims[0] = 1; csParams.dims[1] = 2;
//   CubeSimplexParams *paramsPtr = &csParams;
//   double y[3];
//   y[0] = 0.2; y[1] = 0.05; y[2] = 0.5;
//   BasisInterface csinterface = SetCubeSimplexBasisInterface();
//   Basis *cs = BasisInit((void *)paramsPtr, &csinterface);
//   ComputeBasisIndices(cs);
//   BasisDer(cs, y);
//
//   Basis *cs2 = BasisInit((void *)paramsPtr, &csinterface);
//   ComputeBasisIndices(cs2);
//   BasisPrimeCubeSimplex(csParams.dims, csParams.deg, cs2->basis_indices, y, cs2->basis_der.id);
//   Vector v3 = Vector_init(cs2->basis_funcs.len*(cs2->dim));
//   VectorAddScale(1.0, cs->basis_der, -1.0, cs2->basis_der, v3);
//   VPrint(v3);

////
//   SimplexParams simplexParams;
//   simplexParams.deg = 3; simplexParams.dim = 3;
//   SimplexParams *simplexParamsPtr = &simplexParams;
//   double z[3];
//   z[0] = 0.2; z[1] = 0.05; z[2] = 0.5;
//   BasisInterface simplexInterface = SetSimplexBasisInterface();
//   Basis *simplex = BasisInit((void *)simplexParamsPtr, &simplexInterface);
//   ComputeBasisIndices(simplex);
//   BasisFuncs(simplex, z, simplex->basis_funcs);
//
//   Basis *simplex2 = BasisInit((void *)simplexParamsPtr, &simplexInterface);
//   ComputeBasisIndices(simplex2);
//   BasisSimplex(&simplexParams.dim, simplexParams.deg, simplex2->basis_indices, z, simplex2->basis_funcs.id);
//   Vector v4 = Vector_init(simplex2->basis_funcs.len);
//   VectorAddScale(1.0, simplex->basis_funcs, -1.0, simplex2->basis_funcs, v4);
//   VPrint(v4);

//   BasisFree(cs);
//   BasisFree(cs2);
}

void SimplexFuncsPolytopicOne(MixedPolytopeBasis *basis, const double *x, Vector v)
{
   int i, j, k, d;
   int deg        = basis->deg;
   int dim        = basis->dim;
   int dim1       = basis->params->dims[0];
   int numFuncs   = basis->numFuncs;
   int degPlus1   = deg+1;
   INT_8 *basisId = basis->basis_indices;
   double *phi    = v.id;

   double legendre[degPlus1];
   double dxlegendre[degPlus1];
   double jacobi[SQUARE(degPlus1)*(dim1-1)];
   for(k = 0; k < numFuncs; ++k) phi[k] = 1.0;

   LegendrePoly(degPlus1, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], legendre, dxlegendre);
   double xCoord[dim1-1];
   xCoord[dim1-2] = x[0];
   for(d = 0; d < dim1-2; ++d) xCoord[d] = x[dim1-d-2]/x[dim1-d-3];
   double jCoord[dim1-1]; for(d = 0; d < dim1-1; ++d) jCoord[d] = 1.0;

   for(d = 1; d < dim1; ++d)
   {
      int dimCur = d-1;
      if(d < dim1-1) {
         jCoord[dimCur] /= x[dim1-d-2];
         jCoord[dimCur] = 1.0-2.0*x[dim1-d-1] * jCoord[dimCur];
      }
      if(d == dim1-1) {
         jCoord[dimCur] *= x[0];
         jCoord[dimCur] = 1.0-2.0*jCoord[dimCur];
      }

      for(j = 0; j < deg+1; ++j) {
         double alpha = 2*j+d;
         JacobiPoly(degPlus1, jCoord[dimCur], alpha, 0.0, &jacobi[degPlus1*(d-1)*deg+degPlus1*j]);
      }
   }

   for(k = 0; k < numFuncs; ++k)
   {
      const INT_8 *index = &basisId[k*dim];
      double phiTemp = 1.0;
      for(d = 1; d < dim1; ++d) {
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= xCoord[d-1];
         phiTemp *= jacobi[*(index+d) + degPlus1*(d-1)*deg + degPlus1*xPower] * xFactor;
      }
      phi[k] *= phiTemp*legendre[*index];
   }
}

void SimplexFuncsPolytopicTwo(MixedPolytopeBasis *basis, const double *x, Vector v)
{
   int i, j, k, d;
   int deg        = basis->deg;
   int dim        = basis->dim;
   int dim1       = basis->params->dims[0];
   int dim2       = basis->params->dims[1];
   int dimTwo     = dim1+dim2;
   int numFuncs   = basis->numFuncs;
   int degPlus1   = deg+1;
   double *phi    = v.id;
   INT_8 *basisId = basis->basis_indices;

   double legendre[degPlus1];
   double dxlegendre[degPlus1];
   double jacobi[SQUARE(degPlus1)*(dim2-1)];
   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;

   LegendrePoly(degPlus1, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], legendre, dxlegendre);
   double xCoord[dim2-1];
   xCoord[dim2-2] = x[dim1];
   for(d = 0; d < dim2-2; ++d) xCoord[d] = x[dimTwo-d-2]/x[dimTwo-d-3];
   double jCoord[dim2-1]; for(d = 0; d < dim2-1; ++d) jCoord[d] = 1.0;

   for(d = 1; d < dim2; ++d)
   {
      int dimCur = d-1;
      if(d < dim2-1) {
         jCoord[dimCur] /= x[dimTwo-d-2];
         jCoord[dimCur] = 1.0-2.0*x[dimTwo-d-1] * jCoord[dimCur];
      }
      if(d == dim2-1) {
         jCoord[dimCur] *= x[dim1];
         jCoord[dimCur] = 1.0-2.0*jCoord[dimCur];
      }
      for(j = 0; j < deg+1; ++j) {
         double alpha = 2*j+d;
         JacobiPoly(degPlus1, jCoord[dimCur], alpha, 0.0, &jacobi[degPlus1*(d-1)*deg+degPlus1*j]);
      }
   }

   for(k = 0; k < numFuncs; ++k)
   {
      const INT_8 *index = &basisId[k*dim+dim1];
      double phiTemp = 1.0;
      for(d = 1; d < dim2; ++d) {
         int xPower = 0; for(i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = 1.0; for(i = 0; i < xPower; ++i) xFactor *= xCoord[d-1];
         phiTemp *= jacobi[*(index+d) + degPlus1*(d-1)*deg + degPlus1*xPower] * xFactor;
      }
      phi[k] *= phiTemp*legendre[*index];
   }
}


