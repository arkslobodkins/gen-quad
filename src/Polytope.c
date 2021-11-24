#include "BasisIndices.h"
#include "BasisIntegrals.h"
#include "BasisFunctions.h"
#include "Polytope.h"
#include "Constraints.h"

#include <stdio.h>
#include <stdlib.h>

Polytope * PolytopeInit(void *init_params, PolytopeInterface *interface)
{
    Polytope *polytope  = interface->polytope_init(init_params);
    polytope->interface = interface;

    int dim = polytope->interface->dim(polytope);
    int deg = polytope->interface->deg(polytope);
    polytope->basisSize = BasisSize(deg, dim);
    polytope->basis_indices_flag = OFF;
    polytope->basis_indices   = NULL;

    return polytope;
}

void PolytopeFree(Polytope* polytope)
{
   if(polytope->basis_indices != NULL) {
      free(polytope->basis_indices);
      polytope->basis_indices = NULL;
   }

   polytope->interface->polytope_free(polytope);
}

const char* PolytopeString(Polytope *polytope)
{
   return polytope->interface->string();
}

DOMAIN_TYPE PolytopeType(Polytope *polytope)
{
   return polytope->interface->type(polytope);
}

int PolytopeDim(Polytope *polytope)
{
    return (polytope->interface->dim)(polytope);
}

int PolytopeDeg(Polytope *polytope)
{
    return (polytope->interface->deg)(polytope);
}

double PolytopeGetArea(Polytope *polytope)
{
    return (polytope->interface->get_area)(polytope);
}

void ComputeBasisIndices(Polytope * polytope)
{
   if(polytope->basis_indices_flag == ON) {
      fprintf(stderr, "basis indices have been computed already\n");
      return;
   }
    int dim = polytope->interface->dim(polytope);
    int deg = polytope->interface->deg(polytope);
    polytope->basis_indices = (INT_8 *)
       malloc(polytope->basisSize*dim*sizeof(INT_8));
    BasisIndices(deg, dim, polytope->basis_indices);
    polytope->basis_indices_flag = ON;
}

void ComputeBasisIntegrals(Polytope* polytope)
{
   polytope->interface->ComputeBasisIntegrals(polytope);
}

void ComputeBasisFunctions(Polytope* polytope, double *point)
{
   polytope->interface->ComputeBasisFunctions(polytope, point);
}



Cube * CubeInit(InitialCubeParams *params, PolytopeInterface *interface)
{
    Cube *cube    = (Cube *)malloc(sizeof(Cube));
    cube->constr  = CubeConstrInit(cube);

    cube->params  = (CubeParams *)malloc(sizeof(CubeParams));
    cube->params->dim = params->dim;
    cube->params->deg = params->deg;
    cube->params->D   = CUBE;

    cube->basis_functions = NULL;
    cube->basis_integrals = NULL;
    return cube;
}

double CubeGetArea(Cube *cube)
{
   return 1.0;
}

int CubeDim(Cube *cube)
{
   return cube->params->dim;
}

int CubeDeg(Cube *cube)
{
   return cube->params->deg;
}

DOMAIN_TYPE CubeType(Cube *cube)
{
   return cube->params->D;
}

const char *GetCubeString()
{
   return "CUBE";
}

void ComputeBasisIntegralsCube(Cube* cube)
{
   int dim = cube->params->dim;
   int deg = cube->params->deg;
   cube->basis_integrals = (double *)malloc(cube->basisSize*size_double);
   BasisIntegralsCube(&dim, deg, cube->basis_integrals);
   // BasisIntegrals will be fixed from dims to dim
}

void ComputeBasisFunctionsCube(Cube* cube, double *point)
{
   if(cube->basis_indices_flag != ON)
      fprintf(stderr, "basis_indices have not been computed\n");

   int dim = cube->params->dim;
   int deg = cube->params->deg;
   cube->basis_functions = (double *)malloc(cube->basisSize*size_double);
   BasisCube(&dim, deg, cube->basis_indices, point, cube->basis_functions);
}

constraints * CubeConstrInit(Cube *cube)
{
   constraints *c = (constraints *) malloc(sizeof(constraints));
   return c;
}

void CubeConstrFree(Cube *cube)
{
   if(cube->constr != NULL) { free(cube->constr); cube->constr = NULL; }
}

void CubeFree(Cube *cube)
{
   CubeConstrFree(cube);
   if(cube->basis_functions != NULL) {
      free(cube->basis_functions);
      cube->basis_functions = NULL;
   }
   if(cube->basis_integrals != NULL) {
      free(cube->basis_integrals);
      cube->basis_integrals = NULL;
   }
   if(cube->params != NULL) { free(cube->params); cube->params = NULL; }
   if(cube != NULL) { free(cube); cube = NULL; }
}

PolytopeInterface SetCubeInterface()
{
   PolytopeInterface CubeInterface;
   CubeInterface.polytope_init = (_PolytopeInit)&CubeInit;
   CubeInterface.polytope_free = (_PolytopeFree)&CubeFree;
   CubeInterface.dim           = (_dim)&CubeDim;
   CubeInterface.deg           = (_deg)&CubeDeg;
   CubeInterface.type          = (_type)&CubeType;
   CubeInterface.string        = (_string)&GetCubeString;
   CubeInterface.get_area      = (_get_area)&CubeGetArea;
   CubeInterface.ComputeBasisIntegrals =
      (_ComputeBasisIntegrals)&ComputeBasisIntegralsCube;
   CubeInterface.ComputeBasisFunctions =
      (_ComputeBasisFunctions)&ComputeBasisFunctionsCube;

   return CubeInterface;
}

void TestPolytope()
{
   PolytopeInterface CubeInterface = SetCubeInterface();
   int dim = 2; int deg = 2;
   InitialCubeParams params = {dim, deg};
   Polytope *polytope = PolytopeInit((void *)&params, &CubeInterface);

   printf("polytope dim = %i\n", PolytopeDim(polytope));
   printf("polytope deg = %i\n", PolytopeDeg(polytope));
   printf("polytope enum id = %i\n", PolytopeType(polytope));
   printf("polytope name = %s\n", PolytopeString(polytope));
   printf("polytope area = %lf\n", PolytopeGetArea(polytope));
   printf("polytope basis size = %i\n", polytope->basisSize);
   ComputeBasisIndices(polytope);
   ComputeBasisIntegrals(polytope);
   double x[2]; x[0] = 0.5; x[1] = 0.5;
   ComputeBasisFunctions(polytope, x);


   PolytopeFree(polytope);
}
