/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_dcabs1.c, normal z -> c, Mon Nov 22 19:22:22 2021
 *
 **/

#include "core_blas.h"

#include <math.h>

/***************************************************************************//**
 *
 * @ingroup core_cabs1
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 *******************************************************************************
 *
 * @retval Complex 1-norm absolute value: abs(real(alpha)) + abs(imag(alpha)).
 *
 *******************************************************************************
 *
 * @sa core_scabs1
 *
 ******************************************************************************/
float core_scabs1(plasma_complex32_t alpha)
{
    return fabsf(creal(alpha)) + fabsf(cimag(alpha));
}
