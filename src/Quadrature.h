#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif


quadrature * quadrature_init(int n, int dim, int *dims, int deg, DOMAIN_TYPE D);
void quad_set_funcs_and_constr(quadrature *q);
void quadrature_realloc(int n, int dim, int *dims, int deg, quadrature *q);
void quadrature_reinit(int n, quadrature *q);
void quadrature_assign(const quadrature *quad1, quadrature *quad2);
void quadrature_remove_element(int index, quadrature *quad);
void quadrature_to_vector(const quadrature quad, Vector v);
void quadrature_get_elem(const_quadrature *q, int i, Vector v);
void vector_to_quadrature(const Vector v, quadrature quad);
void quadrature_free(quadrature *quad);

#ifdef __cplusplus
}
#endif

#endif
