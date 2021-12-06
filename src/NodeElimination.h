#ifndef NODE_ELIMINATION_H
#define NODE_ELIMINATION_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

#define PASSED true
#define FAILED false

void NodeElimination(const quadrature *quad_initial, quadrature *quad_final, history *hist);

#ifdef __cplusplus
}
#endif

#endif
