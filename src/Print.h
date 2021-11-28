#ifndef PRINT_H
#define PRINT_H

#include "GENERAL_QUADRATURE.h"
#include <stdio.h>

/**********************************************
\* Functions for printing/debugging purposes \*
**********************************************/

#ifdef __cplusplus
extern "C" {
#endif

void PrintNodes(const quadrature *q, const char *name);
void PrintNodesAndWeights(const quadrature *q, const char *name);
void PrintNodesAndWeightsToFile(const quadrature *q, const char *name, FILE *file);
void PrintNodeAndWeight(int id, const quadrature *q, const char *name);
void PrintQuadLSQInfo(int iters, double error_norm, const quadrature *q, const char *name);
void PrintElimInfo(int dim, int num_nodes, int opt, double opt_factor);

void PrintInt(int x, const char *name);
void PrintBool(bool x, const char *name);
void PrintFloat(float x, const char *name);
void PrintDouble(double x, const char *name);
void Print(const char *x);
void PRINT_ERR(const char *x, int line_num, const char *file_name);
void PRINT_WARN(const char *x, int line_num, const char *file_name);
const char *ERR_STRING(int return_val);

void PrintHistElem(void *data);
#ifdef __cplusplus
}
#endif

#endif
