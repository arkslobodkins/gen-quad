#ifndef QUAD_DRIVER_H
#define QUAD_DRIVER_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

// Receives dimension sizes and degree of precision as command line arguments.
// Performs tests to test appropriateness of the input arguments. If inputs are
// valid, it proceeds to recursive initial guess procedure and Node Elimination algorithm.
void QuadDriver(int argc, char **argv);

#ifdef __cplusplus
}
#endif

#endif
