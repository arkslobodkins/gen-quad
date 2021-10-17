/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "SetDomain.h"

#include "SetParams.h"
#include "Phi.h"
#include "InDomain.h"
#include "BasisIntegrals.h"
#include "GENERAL_QUADRATURE.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void SetInterval(_DomainFuncs *domFuncs)
{
   domFuncs->evalBasis = &PhiCube;
   domFuncs->evalBasisDer = &PhiPrimeCube;
   domFuncs->basisIntegrals = &IntegralsCube;
   domFuncs->inDomain = &InCube;
   domFuncs->inDomainElem = &InCubeElem;
}


void SetCube(_DomainFuncs *domFuncs)
{

   domFuncs->evalBasis = &PhiCube;
   domFuncs->evalBasisDer = &PhiPrimeCube;
   domFuncs->basisIntegrals = &IntegralsCube;
   domFuncs->inDomain = &InCube;
   domFuncs->inDomainElem = &InCubeElem;
}


void SetSimplex(_DomainFuncs *domFuncs)
{
   domFuncs->evalBasis = &PhiSimplex;
   domFuncs->evalBasisDer = &PhiPrimeSimplex;
   domFuncs->basisIntegrals = &IntegralsSimplex;
   domFuncs->inDomain = &InSimplex;
   domFuncs->inDomainElem = &InSimplexElem;
}


void SetCubeSimplex(_DomainFuncs *domFuncs)
{
   domFuncs->evalBasis = &PhiCubeSimplex;
   domFuncs->evalBasisDer = &PhiPrimeCubeSimplex;
   domFuncs->basisIntegrals = &IntegralsCubeSimplex;
   domFuncs->inDomain = &InCubeSimplex;
   domFuncs->inDomainElem = &InCubeSimplexElem;
}


void SetSimplexSimplex(_DomainFuncs *domFuncs)
{
   domFuncs->evalBasis = &PhiSimplexSimplex;
   domFuncs->evalBasisDer = &PhiPrimeSimplexSimplex;
   domFuncs->basisIntegrals = &IntegralsSimplexSimplex;
   domFuncs->inDomain = &InSimplexSimplex;
   domFuncs->inDomainElem = &InSimplexSimplexElem;
}


void SetCubeSimplexSimplex(_DomainFuncs *domFuncs)
{
   domFuncs->evalBasis = &PhiCubeSimplexSimplex;
   domFuncs->evalBasisDer = &PhiPrimeCubeSimplexSimplex;
   domFuncs->basisIntegrals = &IntegralsCubeSimplexSimplex;
   domFuncs->inDomain = &InCubeSimplexSimplex;
   domFuncs->inDomainElem = &InCubeSimplexSimplexElem;
}

