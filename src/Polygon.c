#include <stdlib.h>
#include "Constraints.h"

typedef struct
{
   double *basis_functions;
   double *basis_integrals;
   int *basis_indices;

} Basis;


typedef struct
{
  constraints* (*constraints_init)(void *polygon);
  void         (*get_constraints) (constraints *cons);
} PolygonInterface;


typedef struct
{
    void *shape_params;
    const PolygonInterface *interface;
    constraints *constr;
    Basis *basis;
} Polygon;


typedef struct
{
  int dim;
  int deg;
  DOMAIN_TYPE D;

} CubeParams;


typedef struct
{
   Basis *basis;
   constraints *constr;
   CubeParams params;
} Cube;


typedef struct
{
  constraints* (*constraints_init)(Cube *cube);
  void         (*get_constraints) (Cube *cube);
} CubeInterface;




Polygon *
PolygonCreate(void *shape_params, PolygonInterface *interface)
{
    Polygon *polygon   = (Polygon *)malloc(sizeof(Polygon));
    polygon->shape_params  = shape_params;
    polygon->interface = interface;
    return polygon;
}


constraints * PolygonConstraintsInit(Polygon *polygon)
{
    return (polygon->interface->constraints_init)(polygon->constr);
}


void PolygonGetConstraints(Polygon *polygon)
{
    return (polygon->interface->get_constraints)(polygon->constr);
}
