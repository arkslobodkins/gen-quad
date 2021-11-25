#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

constraints *constraints_interval_init(int dims[1]);
void constraints_interval_realloc(constraints *constr, int dims[1]);
void get_constraints_interval(constraints *constr);
void constraints_interval_free(constraints *constr);

constraints *constraints_cube_init(int dims[1]);
void constraints_cube_realloc(constraints *constr, int dims[1]);
void get_constraints_cube(constraints *constr);
void constraints_cube_free(constraints *constr);

constraints *constraints_simplex_init(int dims[1]);
void constraints_simplex_realloc(constraints *constr, int dims[1]);
void get_constraints_simplex(constraints *constr);
void constraints_simplex_free(constraints *constr);

constraints *constraints_cubesimplex_init(int dims[1]);
void get_constraints_cubesimplex(constraints *constr);
void constraints_cubesimplex_free(constraints *constr);
void constraints_cubesimplex_realloc(constraints *constr, int dims[2]);

constraints *constraints_simplexsimplex_init(int dims[2]);
void get_constraints_simplexsimplex(constraints *constr);
void constraints_simplexsimplex_free(constraints *constr);
void constraints_simplexsimplex_realloc(constraints *constr, int dims[2]);

constraints *constraints_cubesimplexsimplex_init(int dims[3]);
void get_constraints_cubesimplexsimplex(constraints *constr);
void constraints_cubesimplexsimplex_free(constraints *constr);
void constraints_cubesimplexsimplex_realloc(constraints *constr, int dims[3]);

void TestConstraints();
void PrintConstraints(constraints *constr);

#ifdef __cplusplus
}
#endif

#endif
