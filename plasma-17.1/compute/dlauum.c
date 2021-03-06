/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/zlauum.c, normal z -> d, Thu Mar 10 18:57:24 2022
 *
 **/

#include "plasma.h"
#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

/***************************************************************************//**
 *
 * @ingroup plasma_lauum
 *
 *  Computes the product U * U^T or L^T * L, where the triangular
 *  factor U or L is stored in the upper or lower triangular part of
 *  the array A.
 *
 *  If uplo = 'U' or 'u' then the upper triangle of the result is stored,
 *  overwriting the factor U in A.
 *  If uplo = 'L' or 'l' then the lower triangle of the result is stored,
 *  overwriting the factor L in A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaUpper: Upper triangle of A is stored;
 *          - PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] n
 *          The order of the matrix A. n >= 0.
 *
 * @param[in,out] pA
 *          On entry, the triangular factor U or L.
 *          On exit, if UPLO = 'U', the upper triangle of A is
 *          overwritten with the upper triangle of the product U * U^T;
 *          if UPLO = 'L', the lower triangle of A is overwritten with
 *          the lower triangle of the product L^T * L.
 *          The diagonal is assumed to be real with no imaginary part.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,n).
 *
 *******************************************************************************
 *
 * @retval PlasmaSuccess successful exit.
 * @retval < 0 if -i, the i-th argument had an illegal value.
 *
 *******************************************************************************
 *
 * @sa plasma_clauum
 * @sa plasma_dlauum
 * @sa plasma_slauum
 *
 ******************************************************************************/
int plasma_dlauum(plasma_enum_t uplo, int n,
                  double *pA, int lda)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA not initialized");
        return PlasmaErrorNotInitialized;
    }

    // Check input arguments.
    if (uplo != PlasmaUpper && uplo != PlasmaLower) {
        plasma_error("illegal value of uplo");
        return -1;
    }
    if (n < 0) {
        plasma_error("illegal value of n");
        return -2;
    }
    if (lda < imax(1, n)) {
        plasma_error("illegal value of lda");
        return -4;
    }
    // quick return
    if (imax(n, 0) == 0)
        return PlasmaSuccess;

    // Set tiling parameters.
    int nb = plasma->nb;

    // Create tile matrix.
    plasma_desc_t A;
    int retval;
    retval = plasma_desc_general_create(PlasmaRealDouble, nb, nb,
                                        n, n, 0, 0, n, n, &A);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_desc_general_create() failed");
        return retval;
    }
    // Create sequence.
    plasma_sequence_t *sequence = NULL;
    retval = plasma_sequence_create(&sequence);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_sequence_create() failed");
        return retval;
    }

    // Initialize request.
    plasma_request_t request = PlasmaRequestInitializer;

    // Asynchronous block.
    #pragma omp parallel
    #pragma omp master
    {
        // Translate to tile layout.
        plasma_omp_dge2desc(pA, lda, A, sequence, &request);

        // Call the tile async function.
        plasma_omp_dlauum(uplo, A, sequence, &request);

        // Translate back to LAPACK layout.
        plasma_omp_ddesc2ge(A, pA, lda, sequence, &request);
    }
        // Implicit synchronization.

    // Free matrix A in tile layout.
    plasma_desc_destroy(&A);

    // Return status.
    int status = sequence->status;
    plasma_sequence_destroy(sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup plasma_lauum
 *  Computes the product U * U^T or L^T * L, where the
 *  triangular factor U or L is stored in the upper or lower triangular part of
 *  the array A.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaUpper: Upper triangle of A is stored;
 *          - PlasmaLower: Lower triangle of A is stored.
 *
 *
 * @param[in] A
 *          Descriptor of matrix A.
 *
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).  Check
 *          the sequence->status for errors.
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 * @retval void
 *          Errors are returned by setting sequence->status and
 *          request->status to error values.  The sequence->status and
 *          request->status should never be set to PlasmaSuccess (the
 *          initial values) since another async call may be setting a
 *          failure value at the same time.
 *
 *******************************************************************************
 *
 * @sa plasma_dlauum
 * @sa plasma_omp_dlauum
 * @sa plasma_omp_dlauum
 * @sa plasma_omp_clauum
 * @sa plasma_omp_slauum
 *
 ******************************************************************************/
void plasma_omp_dlauum(plasma_enum_t uplo, plasma_desc_t A,
                       plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA not initialized");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }

    // Check input arguments.
    if ((uplo != PlasmaUpper) &&
        (uplo != PlasmaLower)) {
        plasma_error("illegal value of uplo");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (plasma_desc_check(A) != PlasmaSuccess) {
        plasma_error("invalid A");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (sequence == NULL) {
        plasma_error("NULL sequence");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (request == NULL) {
        plasma_error("NULL request");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }


    // Quick return
    if (A.n == 0)
        return;

    // Call the parallel function.
    plasma_pdlauum(uplo, A, sequence, request);
}
