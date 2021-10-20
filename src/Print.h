#ifndef PRINT_H
#define PRINT_H

#include "GENERAL_QUADRATURE.h"

/**********************************************
\* Functions for printing/debugging purposes \*
**********************************************/

#ifdef __cplusplus
extern "C" {
#endif

void PrintNodes(const_quadrature *q, const char *name);
void PrintNodesAndWeights(const_quadrature *q, const char *name);
void PrintNodeAndWeight(int id, const_quadrature *q, const char *name);
void PrintNodeInfo(int iters, double error_norm, const_quadrature *q, const char *name);
void PrintElimInfo(int dim, int num_nodes, int opt, double opt_factor);

void PrintInt(int x, const char *name);
void PrintBool(bool x, const char *name);
void PrintFloat(float x, const char *name);
void PrintDouble(double x, const char *name);
void Print(const char *x);
void PRINT_ERR(const char *x, int line_num, const char *file_name);

#ifdef __cplusplus
}
#endif

#endif
