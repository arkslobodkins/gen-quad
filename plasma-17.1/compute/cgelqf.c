/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/zgelqf.c, normal z -> c, Thu Mar 10 18:58:09 2022
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
 * @ingroup plasma_gelqf
 *
 *  Computes tile LQ factorization of a complex m-by-n matrix A.
 *  The factorization has the form
 *    \f[ A = L \times Q \f],
 *  where L is a lower trapezoidal with positive diagonal and Q is a matrix with
 *  orthonormal rows.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrix A. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A. n >= 0.
 *
 * @param[in,out] pA
 *          On entry, pointer to the m-by-n matrix A.
 *          On exit, the elements on and below the diagonal of the array
 *          contain the m-by-min(m,n) lower trapezoidal matrix L (L is lower
 *          triangular if M <= N); the elements above the diagonal represent
 *          the unitary matrix Q as a product of elementary reflectors, stored
 *          by tiles.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,m).
 *
 * @param[out] T
 *          On exit, auxiliary factorization data, required by plasma_cgelqs
 *          to solve the system of equations.
 *          Matrix of T is allocated inside this function and needs to be
 *          destroyed by plasma_desc_destroy.
 *
 *******************************************************************************
 *
 * @retval PlasmaSuccess successful exit
 * @retval < 0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa plasma_omp_cgelqf
 * @sa plasma_cgelqf
 * @sa plasma_dgelqf
 * @sa plasma_sgelqf
 * @sa plasma_cgelqs
 *
 ******************************************************************************/
int plasma_cgelqf(int m, int n,
                  plasma_complex32_t *pA, int lda,
                  plasma_desc_t *T)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA not initialized");
        return PlasmaErrorNotInitialized;
    }

    // Check input arguments.
    if (m < 0) {
        plasma_error("illegal value of m");
        return -1;
    }
    if (n < 0) {
        plasma_error("illegal value of n");
        return -2;
    }
    if (lda < imax(1, m)) {
        plasma_error("illegal value of lda");
        return -4;
    }

    // quick return
    if (imin(m, n) == 0)
        return PlasmaSuccess;

    // Set tiling parameters.
    int ib = plasma->ib;
    int nb = plasma->nb;
    int householder_mode = plasma->householder_mode;

    // Create tile matrix.
    plasma_desc_t A;
    int retval;
    retval = plasma_desc_general_create(PlasmaComplexFloat, nb, nb,
                                        m, n, 0, 0, m, n, &A);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_desc_general_create() failed");
        return retval;
    }

    // Prepare descriptor T.
    retval = plasma_descT_create(A, ib, householder_mode, T);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_descT_create() failed");
        return retval;
    }

    // Allocate workspace.
    plasma_workspace_t work;
    size_t lwork = nb + ib*nb;  // gelqt: tau + work
    retval = plasma_workspace_create(&work, lwork, PlasmaComplexFloat);
    if (retval != PlasmaSuccess) {
        plasma_error("plasma_workspace_create() failed");
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

    // asynchronous block
    #pragma omp parallel
    #pragma omp master
    {
        // Translate to tile layout.
        plasma_omp_cge2desc(pA, lda, A, sequence, &request);

        // Call the tile async function.
        plasma_omp_cgelqf(A, *T, work, sequence, &request);

        // Translate back to LAPACK layout.
        plasma_omp_cdesc2ge(A, pA, lda, sequence, &request);
    }
    // implicit synchronization

    plasma_workspace_destroy(&work);

    // Free matrix A in tile layout.
    plasma_desc_destroy(&A);

    // Return status.
    int status = sequence->status;
    plasma_sequence_destroy(sequence);
    return status;
}

/***************************************************************************//**
 *
 * @ingroup plasma_gelqf
 *
 *  Computes the tile LQ factorization of a matrix.
 *  Non-blocking tile version of plasma_cgelqf().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in,out] A
 *          Descriptor of matrix A.
 *          A is stored in the tile layout.
 *
 * @param[out] T
 *          Descriptor of matrix T.
 *          On exit, auxiliary factorization data, required by plasma_cgelqs to
 *          solve the system of equations.
 *
 * @param[in] work
 *          Workspace for the auxiliary arrays needed by some coreblas kernels.
 *          For LQ factorization, contains preallocated space for tau and work
 *          arrays. Allocated by the plasma_workspace_create function.
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
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
 * @sa plasma_cgelqf
 * @sa plasma_omp_cgelqf
 * @sa plasma_omp_dgelqf
 * @sa plasma_omp_sgelqf
 * @sa plasma_omp_cgelqs
 *
 ******************************************************************************/
void plasma_omp_cgelqf(plasma_desc_t A, plasma_desc_t T,
                       plasma_workspace_t work,
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
    if (plasma_desc_check(A) != PlasmaSuccess) {
        plasma_error("invalid A");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (plasma_desc_check(T) != PlasmaSuccess) {
        plasma_error("invalid T");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (sequence == NULL) {
        plasma_fatal_error("NULL sequence");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }
    if (request == NULL) {
        plasma_fatal_error("NULL request");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }

    // quick return
    if (imin(A.m, A.n) == 0)
        return;

    // Call the parallel function.
    if (plasma->householder_mode == PlasmaTreeHouseholder) {
        plasma_pcgelqfrh(A, T, work, sequence, request);
    }
    else {
        plasma_pcgelqf(A, T, work, sequence, request);
    }
}
