#include "BasisIndices.h"
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
    polytope->basis = (Basis *)malloc(sizeof(Basis));
    polytope->basis->basis_indices   = NULL;
    polytope->basis->basis_functions = NULL;
    polytope->basis->basis_integrals = NULL;

    return polytope;
}

void PolytopeFree(Polytope* polytope)
{
   if(polytope->basis->basis_indices != NULL) { free(polytope->basis->basis_indices);
      polytope->basis->basis_indices = NULL; }
   if(polytope->basis != NULL) { free(polytope->basis); polytope->basis = NULL; }

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

Cube * CubeInit(InitialCubeParams *params, PolytopeInterface *interface)
{
    Cube *cube    = (Cube *)malloc(sizeof(Cube));
    cube->constr  = CubeConstrInit(cube);

    cube->params  = (CubeParams *)malloc(sizeof(CubeParams));
    cube->params->dim = params->dim;
    cube->params->deg = params->deg;
    cube->params->D   = CUBE;


    return cube;
}

void ComputeBasis(Polytope * polytope)
{
    int dim = polytope->interface->dim(polytope);
    int deg = polytope->interface->deg(polytope);
    polytope->basis->basis_indices = (int_fast8_t *)
       malloc(polytope->basisSize*dim*sizeof(int_fast8_t));
    BasisIndices(deg, dim, polytope->basis->basis_indices);
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

   return CubeInterface;
}

void TestPolytope()
{
   PolytopeInterface CubeInterface = SetCubeInterface();
   int dim = 2; int deg = 5;
   InitialCubeParams params = {dim, deg};
   Polytope *polytope = PolytopeInit((void *)&params, &CubeInterface);

   printf("polytope dim = %i\n", PolytopeDim(polytope));
   printf("polytope deg = %i\n", PolytopeDeg(polytope));
   printf("polytope enum id = %i\n", PolytopeType(polytope));
   printf("polytope name = %s\n", PolytopeString(polytope));
   printf("polytope area = %lf\n", PolytopeGetArea(polytope));
   printf("polytope basis size = %i\n", polytope->basisSize);
   ComputeBasis(polytope);

   PolytopeFree(polytope);
}
