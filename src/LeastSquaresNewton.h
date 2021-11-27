#ifndef LEAST_SQUARES_NEWTONPLASMA_H
#define LEAST_SQUARES_NEWTONPLASMA_H

#include "GENERAL_QUADRATURE.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SOL_FOUND true
#define SOL_NOT_FOUND false


#define CONSTR_NOT_NEEDED 1

#define CONSTR_SUCCESS  0
#define CONSTR_FAIL -1
#define CONSTR_INV_INPUT -2
#define CONSTR_UNEXPECTED -3

bool LeastSquaresNewton(LibraryType LibraryType, const bool_enum FLAG_CONSTR, const INT_8 *basis, quadrature *q_orig, int *its);


#ifdef __cplusplus
}
#endif

#endif


