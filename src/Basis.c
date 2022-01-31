#include "Basis.h"
#include "Quadrature.h"
#include "Gauss_Lib/Jacobi.h"
#include "GeneralGaussTensor.h"
#include "AddDimension.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void BasisIndices(int deg, int dim, INT_8 *f);
static int BasisSize(int deg, int dim);

static Table3d Table3dCreate(int deg, int dim);
static void Table3dFree(Table3d table);
static void ComputeTable(Table3d table);

static void LegendrePoly(int order, double x, double *p);
static void LegendrePolyAndPrime(int order, double x, double *p, double *dp);
static void JacobiPolyBetaZero(int order, double x, int alpha, double *p);
static void JacobiPolyWithTable(int order, double x, int alpha, double *p, Table3d table);

static void IntegralsCubePolyhedralMonomial(MixedPolytopeBasis *basis, Vector v);
static void IntegralsSimplexPolyhedralMonomialOne(MixedPolytopeBasis *basis, Vector v);
static void IntegralsSimplexPolyhedralMonomialTwo(MixedPolytopeBasis *basis, Vector v);

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

Basis* BasisInit(void *params, BasisInterface *interface)
{
   Basis *basis = interface->basisInit(params);
   basis->interface = (BasisInterface *)malloc(sizeof(BasisInterface));
   *(basis->interface) = *interface;

   int numFuncs = basis->numFuncs;
   int dim = basis->dim;
   basis->indices = (INT_8 *) malloc(numFuncs*dim*sizeof(INT_8));
   BasisIndices(basis->deg, basis->dim, basis->indices);
   basis->functions   = Vector_init(numFuncs);
   basis->derivatives = Vector_init(numFuncs*dim);
   basis->integrals   = Vector_init(numFuncs);
   return basis;
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
   if(basis->indices) {
      free(basis->indices);
      basis->indices = NULL;
   }

   basis->interface->basisFree(basis);
   if(interface) {
      free(interface);
      interface = NULL;
   }
}

void BasisMonomial(Basis *basis, const double *x, Vector phi)
{
   int numFuncs    = basis->numFuncs;
   int dim         = basis->dim;
   INT_8 *basis_id = basis->indices;

   for(int k = 0; k < numFuncs; ++k) phi.id[k] = 1.0;
   for(int k = 0; k < numFuncs; ++k) {
      for(int d = 0; d < dim; ++d) {
         INT_8 basis_power = basis_id[k*dim+d];
         double product = DoubleIntPower(x[d], basis_power);
         phi.id[k] *= product;
      }
   }
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetCubeBasisInterface()
{
   BasisInterface cubeInterface;
   cubeInterface.basisInit        = (BasisInitPtr)&CubeBasisInit;
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
   INT_8 *basisId = basis->indices;
   double *phi    = v.id;

   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;

   for(int d = 0; d < dim; ++d) {
      double legendre[deg+1];
      LegendrePoly(deg+1, 2*x[d]-1, legendre);
      for(int k = 0; k < numFuncs; ++k)
         phi[k] *= legendre[basisId[k*dim+d]];
   }
}


void ComputeCubeBasisDer(CubeBasis *basis, const double *x, Vector v)
{
   int deg          = basis->deg;
   int dim          = basis->dim;
   int numFuncs     = basis->numFuncs;
   int degPlus1     = deg+1;
   INT_8 *basisId   = basis->indices;
   double *phiPrime = v.id;

   double legendre[degPlus1*dim];
   double dxlegendre[degPlus1*dim];

   for(int d = 0; d < dim; ++d)
      LegendrePolyAndPrime(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

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


void CubeBasisIntegrals(CubeBasis *basis, Vector v)
{
   int len = basis->integrals.len;
   double *integrals = v.id;
   memset(integrals, 0, SIZE_DOUBLE(len));
   integrals[0] = 1.0;
}


void CubeBasisIntegralsMonomial(CubeBasis *basis, Vector v)
{
   int i, d;
   int dim = basis->dim;
   int numFuncs = basis->numFuncs;
   INT_8 *basisIndices = basis->indices;
   Vector integrals = v;

   for(i = 0; i < numFuncs; ++i) {
      double val = 1.0;
      for(d = 0; d < dim; ++d)
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
      free(basis->addData);
      basis->addData = NULL;
   }
   free(basis); basis = NULL;
}



/////////////////////////////////////////////////////////////////////////////
BasisInterface SetSimplexBasisInterface()
{
   BasisInterface simplexInterface;
   simplexInterface.basisInit        = (BasisInitPtr)&SimplexBasisInit;
   simplexInterface.computeFuncs     = (BasisFuncsPtr)&SimplexBasisFuncs;
   simplexInterface.computeDer       = (BasisDerPtr)&SimplexBasisDer;
   simplexInterface.computeIntegrals = (BasisIntegralsPtr)&SimplexBasisIntegrals;
   simplexInterface.computeIntegralsMonomial = (BasisIntegralsMonomialPtr)&SimplexBasisIntegralsMonomial;
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

   basis->table = Table3dCreate(deg, dim);
   ComputeTable(basis->table);
   return basis;
}


void SimplexBasisFuncs(SimplexBasis *basis, const double *x, Vector v)
{
   int deg        = basis->deg;
   int dim        = basis->dim;
   int numFuncs   = basis->numFuncs;
   int degPlus1   = deg+1;
   double legendre[degPlus1];
   double jacobi[SQUARE(degPlus1)*(dim-1)];

   INT_8 *basisId = basis->indices;

   LegendrePoly(degPlus1, (2.0*x[dim-1]-x[dim-2])/x[dim-2], legendre);
   double xCoord[dim-1];
   for(int d = 0; d < dim-2; ++d) xCoord[d] = x[dim-d-2]/x[dim-d-3];
   xCoord[dim-2] = x[0];

   double jCoord[dim];
   for(int d = 1; d < dim-1; ++d)
      jCoord[d] = 1.0-2.0*x[dim-d-1]/x[dim-d-2];
   jCoord[dim-1] = 1.0-2.0*x[0];

   for(int d = 1; d < dim; ++d) {
      int nextDim = (deg+1)*(d-1)*deg;
      for(int j = 0; j < degPlus1; ++j) {
         int nextAlpha = (deg+1)*j;
         JacobiPolyWithTable(deg+1, jCoord[d], 2*j+d, &jacobi[nextDim+nextAlpha], basis->table);
      }
   }

   double *phi = v.id;
   for(int k = 0; k < numFuncs; ++k) phi[k] = legendre[basisId[k*dim]];
   for(int d = 1; d < dim; ++d)
   {
      double coordLoc = xCoord[d-1];
      for(int k = 0; k < numFuncs; ++k)
      {
         const INT_8 *index = &basisId[k*dim];
         int xPower = 0; for(int i = 0; i < d; ++i) xPower += *(index+i);
         double xFactor = DoubleIntPower(coordLoc, xPower);
         double phiTemp = jacobi[*(index+d) + degPlus1*(d-1)*deg + degPlus1*xPower] * xFactor;
         phi[k] *= phiTemp;
      }
   }
}


void SimplexBasisDer(SimplexBasis *basis, const double *x, Vector v)
{
   int dim = basis->params->dim;
   int numFuncs = basis->numFuncs;
   double h = POW_DOUBLE(10, -5)*5.0;

   double x_backw1[dim];
   double x_backw2[dim];
   double x_forw1[dim];
   double x_forw2[dim];

   double *basisDer = v.id;
   Vector phi_backw2 = basis->addData->phi_backw2;
   Vector phi_backw1 = basis->addData->phi_backw1;
   Vector phi_forw1  = basis->addData->phi_forw1;
   Vector phi_forw2  = basis->addData->phi_forw2;

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
   int numFuncs = basis->numFuncs;
   int dim = basis->dim;
   INT_8 *basisIndices = basis->indices;
   Vector integrals = basis->integrals;
   double val, power;

   for(int i = 0; i < numFuncs; ++i)
   {
      val = 1.0;
      for(int d = 0; d < dim; ++d) {
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

   Vector_free(basis->addData->phi_backw1);
   Vector_free(basis->addData->phi_backw2);
   Vector_free(basis->addData->phi_forw1);
   Vector_free(basis->addData->phi_forw2);
   Table3dFree(basis->table);
   if(basis->params) {
      free(basis->params);
      basis->params = NULL;
   }
   if(basis->addData) {
      free(basis->addData);
      basis->addData = NULL;
   }
   free(basis); basis = NULL;
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetCubeSimplexBasisInterface()
{
   BasisInterface csInterface;
   csInterface.basisInit        = (BasisInitPtr)&CubeSimplexBasisInit;
   csInterface.computeFuncs     = (BasisFuncsPtr)&CubeSimplexBasisFuncs;
   csInterface.computeDer       = (BasisDerPtr)&CubeSimplexBasisDer;
   csInterface.computeIntegrals = (BasisIntegralsPtr)&CubeSimplexBasisIntegrals;
   csInterface.computeIntegralsMonomial = (BasisIntegralsMonomialPtr)&CubeSimplexBasisIntegralsMonomial;
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
   INT_8 *basisId = basis->indices;

   double legendre[degPlus1*dim1];

   for(int d = 0; d < dim1; ++d)
      LegendrePoly(degPlus1, 2*x[d]-1, &legendre[d*degPlus1]);

   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;
   for(int k = 0; k < numFuncs; ++k)
      for(int d = 0; d < dim1; ++d)
         phi[k] *= legendre[basisId[k*dim+d]+degPlus1*d];

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, basis->addData->basis_polytopic);
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
   int degPlus1     = deg+1;
   double *phiPrime = v.id;
   INT_8 *basisId = basis->indices;

   double legendre[degPlus1*dim1];
   double dxlegendre[degPlus1*dim1];

   for(int d = 0; d < dim1; ++d)
      LegendrePolyAndPrime(degPlus1, 2*x[d]-1, &legendre[d*degPlus1], &dxlegendre[d*degPlus1]);

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

  for(int i = 0; i < numFuncs; ++i)
     integrals.id[i] = integralsCube.id[i]*integralsSimplex.id[i];

  Vector_free(integralsCube);
  Vector_free(integralsSimplex);
}


void CubeSimplexBasisFree(CubeSimplexBasis *basis)
{
   if(basis == NULL) return;

   Vector_free(basis->addData->basis_polytopic);
   Vector_free(basis->addData->phi_backw1);
   Vector_free(basis->addData->phi_backw2);
   Vector_free(basis->addData->phi_forw1);
   Vector_free(basis->addData->phi_forw2);
   if(basis->params) {
      free(basis->params);
      basis->params = NULL;
   }
   if(basis->addData) {
      free(basis->addData);
      basis->addData = NULL;
   }
   free(basis); basis = NULL;
}


/////////////////////////////////////////////////////////////////////////////
BasisInterface SetSimplexSimplexBasisInterface()
{
   BasisInterface ssInterface;
   ssInterface.basisInit        = (BasisInitPtr)&SimplexSimplexBasisInit;
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
   double *phi      = basis->functions.id;
   Vector polytopic = basis->addData->basis_polytopic;

   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;

   SimplexFuncsPolytopicOne((MixedPolytopeBasis *)basis, x, polytopic);
   for(int k = 0; k < numFuncs; ++k)
      phi[k] *= polytopic.id[k];

   SimplexFuncsPolytopicTwo((MixedPolytopeBasis *)basis, x, polytopic);
   for(int k = 0; k < numFuncs; ++k)
      phi[k] *= polytopic.id[k];
}


void SimplexSimplexBasisDer(SimplexSimplexBasis *basis, const double *x, Vector v)
{
   int dim1          = basis->params->dims[0];
   int dim2          = basis->params->dims[1];
   int dim           = basis->dim;
   int numFuncs      = basis->numFuncs;
   Vector phiPrime   = v;

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

  for(int i = 0; i < numFuncs; ++i)
     integrals.id[i] = integralsS1.id[i]*integralsS2.id[i];

  Vector_free(integralsS1);
  Vector_free(integralsS2);
}


void SimplexSimplexBasisFree(SimplexSimplexBasis *basis)
{
   if(basis == NULL) return;

   Vector_free(basis->addData->basis_polytopic);
   Vector_free(basis->addData->phi_backw1);
   Vector_free(basis->addData->phi_backw2);
   Vector_free(basis->addData->phi_forw1);
   Vector_free(basis->addData->phi_forw2);
   if(basis->params) {
      free(basis->params);
      basis->params = NULL;
   }
   if(basis->addData) {
      free(basis->addData);
      basis->addData = NULL;
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
static void JacobiPolyBetaZero(int order, double x, int alpha, double *p)
{
   p[0] = 1.0;
   p[1] = 0.5*(x-1.0) * (alpha+2.0) + alpha+1.0;
   for (int k = 1; k < order-1; ++k) {
      double fac1 = (alpha+2.0*k+1.0) * ( (alpha+2.0*k+2.0) * (alpha+2.0*k) * x + (alpha*alpha) );
      double fac2 = -2.0 * (alpha+k) * k * (alpha+2.0*k+2.0);
      double fac3 = 1.0 / (2.0*(k+1.0)*(alpha+k+1.0) * (alpha+2.0*k));
      p[k+1] = fac3 * (fac1*p[k] + fac2*p[k-1]);
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
   int numFuncs = basis->numFuncs;
   int dimCube = basis->params->dims[0];
   int dim = basis->dim;
   Vector integrals = v;
   INT_8* indices = basis->indices;

   for(int i = 0; i < numFuncs; ++i)
   { double val = 1.0;
      for(int d = 0; d < dimCube; ++d)
         val = val/(indices[i*dim+d] + 1);
      integrals.id[i] = val;
   }
}


static void IntegralsSimplexPolyhedralMonomialOne(MixedPolytopeBasis *basis, Vector v)
{
   int numFuncs = basis->numFuncs;
   int dimSimplex = basis->params->dims[0];
   int dim = basis->dim;
   double val, power;
   INT_8* indices = basis->indices;
   Vector integrals = v;

   for(int i = 0; i < numFuncs; ++i)
   {
      val = 1.0;
      for(int d = 0; d < dimSimplex; ++d) {
         power = 0;
         for(int r = 0; r < dimSimplex-d; ++r)
            power += (double)indices[dim*i+dimSimplex-r-1];
         val /= (power+dimSimplex-d);
      }
      integrals.id[i] = val;
   }
}


static void IntegralsSimplexPolyhedralMonomialTwo(MixedPolytopeBasis *basis, Vector v)
{
   int numFuncs = basis->numFuncs;
   int dimSimplex = basis->params->dims[1];
   int dim = basis->dim;
   double val, power;
   INT_8* indices = basis->indices;
   Vector integrals = v;

   for(int i = 0; i < numFuncs; ++i)
   {
      val = 1.0;
      for(int d = 0; d < dimSimplex; ++d) {
         power = 0;
         for(int r = 0; r < dimSimplex-d; ++r)
            power += (double)indices[dim*(i+1)-r-1];
         val /= (power+dimSimplex-d);
      }
      integrals.id[i] = val;
   }
}


void SimplexFuncsPolytopicOne(MixedPolytopeBasis *basis, const double *x, Vector v)
{
   int i, j, k, d;
   int deg        = basis->deg;
   int dim        = basis->dim;
   int dim1       = basis->params->dims[0];
   int numFuncs   = basis->numFuncs;
   int degPlus1   = deg+1;
   INT_8 *basisId = basis->indices;
   double *phi    = v.id;

   double legendre[degPlus1];
   double jacobi[SQUARE(degPlus1)*(dim1-1)];
   for(k = 0; k < numFuncs; ++k) phi[k] = 1.0;

   LegendrePoly(degPlus1, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], legendre);
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
         JacobiPolyBetaZero(degPlus1, jCoord[dimCur], alpha, &jacobi[degPlus1*(d-1)*deg+degPlus1*j]);
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
   INT_8 *basisId = basis->indices;

   double legendre[degPlus1];
   double jacobi[SQUARE(degPlus1)*(dim2-1)];
   for(int k = 0; k < numFuncs; ++k) phi[k] = 1.0;

   LegendrePoly(degPlus1, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], legendre);
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
         JacobiPolyBetaZero(degPlus1, jCoord[dimCur], alpha, &jacobi[degPlus1*(d-1)*deg+degPlus1*j]);
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


static void BasisIndices(int deg, int dim, INT_8 *f)
{
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


void orthogonal_cube_basis_test(int deg, int dim)
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
}


void orthogonal_simplex_basis_test(int deg, int dim)
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
}


