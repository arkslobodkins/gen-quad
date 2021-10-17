/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */

#include "SetConstraint.h"

#include "Constraints.h"
#include "GENERAL_QUADRATURE.h"

void SetConstraintCube(ConstraintFuncs *constr)
{
   constr->constr_init = &constraints_cube_init;
   constr->get_constr = &get_constraints_cube;
   constr->constr_free = &constraints_cube_free;
}


void SetConstraintSimplex(ConstraintFuncs *constr)
{
   constr->constr_init = &constraints_simplex_init;
   constr->get_constr = &get_constraints_simplex;
   constr->constr_free = &constraints_simplex_free;
}


void SetConstraintCubeSimplex(ConstraintFuncs *constr)
{
   constr->constr_init = &constraints_cubesimplex_init;
   constr->get_constr = &get_constraints_cubesimplex;
   constr->constr_free = &constraints_cubesimplex_free;
}


void SetConstraintSimplexSimplex(ConstraintFuncs *constr)
{
   constr->constr_init = &constraints_simplexsimplex_init;
   constr->get_constr = &get_constraints_simplexsimplex;
   constr->constr_free = &constraints_simplexsimplex_free;
}


void SetConstraintCubeSimplexSimplex(ConstraintFuncs *constr)
{
   constr->constr_init = &constraints_cubesimplexsimplex_init;
   constr->get_constr = &get_constraints_cubesimplexsimplex;
   constr->constr_free = &constraints_cubesimplexsimplex_free;
}
