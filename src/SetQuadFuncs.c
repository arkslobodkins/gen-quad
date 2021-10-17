
/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "SetQuadFuncs.h"

#include "Phi.h"
#include "InDomain.h"
#include "BasisIntegrals.h"
#include "GENERAL_QUADRATURE.h"
#include "Constraints.h"

void SetIntervalFuncs(quadrature *q)
{
   q->evalBasis = &PhiCube;
   q->evalBasisDer = &PhiPrimeCube;
   q->basisIntegrals = &IntegralsCube;
   q->inDomain = &InCube;

   q->constr_init = &constraints_interval_init;
   q->constr_realloc = &constraints_interval_realloc;
   q->get_constr = &get_constraints_interval;
   q->constr_free = &constraints_interval_free;
}


void SetCubeFuncs(quadrature *q)
{

   q->evalBasis = &PhiCube;
   q->evalBasisDer = &PhiPrimeCube;
   q->basisIntegrals = &IntegralsCube;
   q->inDomain = &InCube;
//   q->inDomainElem = &InCubeElem;

   q->constr_init = &constraints_cube_init;
   q->constr_realloc = &constraints_cube_realloc;
   q->get_constr = &get_constraints_cube;
   q->constr_free = &constraints_cube_free;
}


void SetSimplexFuncs(quadrature *q)
{
   q->evalBasis = &PhiSimplex;
   q->evalBasisDer = &PhiPrimeSimplex;
   q->basisIntegrals = &IntegralsSimplex;
   q->inDomain = &InSimplex;
//   q->inDomainElem = &InSimplexElem;

   q->constr_init = &constraints_simplex_init;
   q->constr_realloc = &constraints_simplex_realloc;
   q->get_constr = &get_constraints_simplex;
   q->constr_free = &constraints_simplex_free;
}


void SetCubeSimplexFuncs(quadrature *q)
{
   q->evalBasis = &PhiCubeSimplex;
   q->evalBasisDer = &PhiPrimeCubeSimplex;
   q->basisIntegrals = &IntegralsCubeSimplex;
   q->inDomain = &InCubeSimplex;
//   q->inDomainElem = &InCubeSimplexElem;

   q->constr_init = &constraints_cubesimplex_init;
   q->constr_realloc = &constraints_cubesimplex_realloc;
   q->get_constr = &get_constraints_cubesimplex;
   q->constr_free = &constraints_cubesimplex_free;
}


void SetSimplexSimplexFuncs(quadrature *q)
{
   q->evalBasis = &PhiSimplexSimplex;
   q->evalBasisDer = &PhiPrimeSimplexSimplex;
   q->basisIntegrals = &IntegralsSimplexSimplex;
   q->inDomain = &InSimplexSimplex;
//   q->inDomainElem = &InSimplexSimplexElem;
//
   q->constr_init = &constraints_simplexsimplex_init;
   q->constr_realloc = &constraints_simplexsimplex_realloc;
   q->get_constr = &get_constraints_simplexsimplex;
   q->constr_free = &constraints_simplexsimplex_free;
}


void SetCubeSimplexSimplexFuncs(quadrature *q)
{
   q->evalBasis = &PhiCubeSimplexSimplex;
   q->evalBasisDer = &PhiPrimeCubeSimplexSimplex;
   q->basisIntegrals = &IntegralsCubeSimplexSimplex;
   q->inDomain = &InCubeSimplexSimplex;
//   q->inDomainElem = &InCubeSimplexSimplexElem;
//
   q->constr_init = &constraints_cubesimplexsimplex_init;
   q->constr_realloc = &constraints_cubesimplexsimplex_realloc;
   q->get_constr = &get_constraints_cubesimplexsimplex;
   q->constr_free = &constraints_cubesimplexsimplex_free;
}

