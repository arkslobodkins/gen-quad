/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef OUTPUT_H
#define OUTPUT_H

#include "GENERAL_QUADRATURE.h"
#include "Quadrature.h"

#ifdef __cplusplus
extern "C" {
#endif

// Prints history and quadrature parameters to a file.
void HistoryToFile(const quadrature *q, int arr_size, history **hist_arr);

// Prints final quadrature to a file
void QuadratureToFile(const quadrature *quad);

// Prints how many nodes contain coordinates < 0.1 and > 0.9, and how many coordinates per node
void BoundaryCubeStats(const quadrature *quad);

#ifdef __cplusplus
}
#endif

#endif
