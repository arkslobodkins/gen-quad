#ifndef NODE_ELIMINATION_H
#define NODE_ELIMINATION_H

#include "GENERAL_QUADRATURE.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PASSED true
#define FAILED false

void NodeElimination(const_quadrature *quad_initial, quadrature *quad_final, elim_history *history);

#ifdef __cplusplus
}
#endif

#endif
