#ifndef OUTPUT_H
#define OUTPUT_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

void HistoryToFile(const quadrature *q, int arr_size, history **hist_arr);
void QuadratureToFile(const quadrature *quad);
void BoundaryCubeStats(const quadrature *quad);

#ifdef __cplusplus
}
#endif

#endif
