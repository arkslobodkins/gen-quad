#ifndef OUTPUT_H
#define OUTPUT_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void Output(int num_nodes_initial, double res, const_quadrature *quad_final, const elim_history history);
void DumpCubatureRule(const_quadrature *quad);

#ifdef __cplusplus
}
#endif

#endif
