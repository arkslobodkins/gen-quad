/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/zgbtrf.c, normal z -> c, Thu Mar 10 18:57:20 2022
 *
 **/

#include "plasma.h"
#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"

#include "mkl_lapacke.h"

/***************************************************************************//**
 *
 * @ingroup plasma_gbtrf
 *
 * Computes an LU factorization of a real m-by-n band matrix A
 * using partial pivoting with row interchanges.
 *
 *******************************************************************************
 *
 * @param[in] m
 *          The number of rows of the matrix A. n >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix A. n >= 0.
 *
 * @param[in] kl
 *          The number of subdiagonals within the band of A. kl >= 0.
 *
 * @param[in] ku
 *          The number of superdiagonals within the band of A. ku >= 0.
 *
 * @param[in,out] AB
 *          Details of the LU factorization of the band matrix A, as
 *          computed by plasma_cgbtrf.
 *
 * @param[in] ldab
 *          The leading dimension of the array AB.
 *
 * @param[out] ipiv
 *          The pivot indices; for 1 <= i <= min(m,n), row i of the
 *          matrix was interchanged with row ipiv(i).
 *
 ******************************************************************************/
int plasma_cgbtrf(int m, int n, int kl, int ku,
                  plasma_complex32_t *pAB, int ldab, int *ipiv)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA not initialized");
        return PlasmaErrorNotInitialized;
    }

    if (m < 0) {
        plasma_error("illegal value of m");
        return -1;
    }
    if (n < 0) {
        plasma_error("illegal value of n");
        return -2;
    }
    if (kl < 0) {
        plasma_error("illegal value of kl");
        return -3;
    }
    if (ku < 0) {
        plasma_error("illegal value of ku");
        return -4;
    }
    if (ldab < imax(1, 1+kl+ku)) {
        plasma_error("illegal value of ldab");
        return -6;
    }

    // quick return

    // Set tiling parameters.
    int nb = plasma->nb;

    // Initialize barrier.
    int num_panel_threads = plasma->num_panel_threads;
    plasma_barrier_init(&plasma->barrier, num_panel_threads);

    // Create tile matrix.
    plasma_desc_t AB;
    int tku = (ku+kl+nb-1)/nb; // number of tiles in upper band (not including diagonal)
    int tkl = (kl+nb-1)/nb;    // number of tiles in lower band (not including diagonal)
    int lm = (tku+tkl+1)*nb;   // since we use cgetrf on panel, we pivot back within panel.
                               // this could fill the last tile of the panel,
                               // and we need extra NB space on the bottom
    int retval;
    retval = plasma_desc_general_band_create(PlasmaComplexFloat, PlasmaGeneral,
                                             nb, nb, lm, n, 0, 0, m, n, kl, ku, &AB);
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

    #pragma omp parallel
    #pragma omp master
    {
        // Translate to tile layout.
        plasma_omp_cpb2desc(pAB, ldab, AB, sequence, &request);
    }

    #pragma omp parallel
    #pragma omp master
    {
        // Call the tile async function.
        plasma_omp_cgbtrf(AB, ipiv, sequence, &request);
    }

    #pragma omp parallel
    #pragma omp master
    {
        // Translate back to LAPACK layout.
        plasma_omp_cdesc2pb(AB, pAB, ldab, sequence, &request);
    }

    // Free matrix A in tile layout.
    plasma_desc_destroy(&AB);

    // Return status.
    int status = sequence->status;
    plasma_sequence_destroy(sequence);
    return status;
}

/***************************************************************************//**
 *
 * Computes an LU factorization of a real m-by-n band matrix A
 * using partial pivoting with row interchanges.
 * Non-blocking tile version of plasma_cgbsv().
 * Operates on matrices stored by tiles.
 * All matrices are passed through descriptors.
 * All dimensions are taken from the descriptors.
 * Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in,out] AB
 *          Descriptor of matrix A.
 *
 * @param[out] ipiv
 *          The pivot indices; for 1 <= i <= min(m,n), row i of the
 *          matrix was interchanged with row ipiv(i).
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
 ******************************************************************************/
void plasma_omp_cgbtrf(plasma_desc_t AB, int *ipiv,
                       plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_fatal_error("PLASMA not initialized");
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        return;
    }

    // Check input arguments.
    if (plasma_desc_check(AB) != PlasmaSuccess) {
        plasma_request_fail(sequence, request, PlasmaErrorIllegalValue);
        plasma_error("invalid AB");
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

    // Call the parallel function.
    plasma_pcgbtrf(AB, ipiv, sequence, request);
}
