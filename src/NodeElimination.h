#ifndef NODE_ELIMINATION_H
#define NODE_ELIMINATION_H

#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif


/***************************************************************************************************
 * A new node elimination scheme that eliminates one node at a time and computes
 * the initial guess for constrained Newton's method. Subsequently, Newton's method is called
 * to obtain quadrature rule with fewer nodes. The procedure is repeated
 * until no more nodes can be eliminated.
 ***************************************************************************************************/
void NodeElimination(const quadrature *quad_initial, quadrature *quad_final, history *hist);

#ifdef __cplusplus
}
#endif

#endif
