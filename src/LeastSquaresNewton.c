/* Arkadijs Slobodkins
 * SMU Mathematics
 * August 2021
 */


#include "LeastSquaresNewton.h"
#include "ConstrainedOptimization.h"

#include "GetFunction.h"
#include "GetJacobian.h"
#include "Quadrature.h"
#include "Print.h"
#include "Conditional_Debug.h"
#include "get_time.h"
#include "LINALG.h"
#include "Basis.h"
#include <omp.h>

#include <stdlib.h>
#include <string.h>
#include <assert.h>


static int CheckForStop(int INFO, double errorNorm, double errorNormPrev, Vector least_sq_sol);

// LeastSquaresNewton
// Receives initial quadrature guess. Primarily solves
// underdetermined systems of equations in the least
// squares sense. Returns success if algorithm converged
// and all nodes are inside of the domain and if all nodes are positive.
double CONSTR_TIME = 0.0;
double LSQ_TIME = 0.0;
#define MAX_ELIM_WEIGHTS 10
bool LeastSquaresNewton(const bool_enum CONSTR_OPT, quadrature *q_orig, int *its)
{
   assert(q_orig->num_nodes >= 1);

   int k         = q_orig->num_nodes;
   int numFuncs  = q_orig->basis->numFuncs;
   int dim       = q_orig->dim;
   bool SOL_FLAG = SOL_NOT_FOUND;

   int itsLoc   = 0;
   int maxiter  = 25;
   double q_tol = QUAD_TOL; // 10^(-14);

   // initialize LAPACK
   int nrows  = numFuncs;
   int ncols = (dim+1)*k;
   int LEAD_DIM = MAX(nrows, ncols);
   int SMALL_DIM = MIN(nrows, ncols);
   int LWORK = 5*ncols;
   CMatrix JACOBIAN = CMatrix_init(nrows, ncols);
   Vector LEAST_SQ_SOL = Vector_init(LEAD_DIM);
   Vector RHS          = Vector_init(nrows);
   double *WORK        = (double *)malloc(LWORK*sizeof(double));

   quadrature *q_prev = quadrature_make_full_copy(q_orig);
   quadrature *q_next = quadrature_make_full_copy(q_orig);

   int elim_weights = 0;
   double errorNorm = 1.0, errorNormPrev = 1.0, errorNormUpdate = 1.0;
   ConstrVectData cVectData = ConstrVectDataInit();

   int (*leastsquares_ptr)(CMatrix A, Vector RHS_TO_X);
   void (*getFunc_ptr)(quadrature *q, Vector f);
   void (*getFunctionAndJacobian_ptr)(quadrature *q, Vector f, CMatrix JACOBIAN);
#ifdef _OPENMP
   QuadAllocBasisOmp(q_next, omp_get_max_threads());
   QuadAllocBasisOmp(q_prev, omp_get_max_threads());
   AllocVectorOmpData(&RHS);
   getFunc_ptr = &GetFunctionOmp;
   getFunctionAndJacobian_ptr = &GetFunctionAndJacobianOmp;
   if(PLASMA_CONDITION()) leastsquares_ptr = DGELS_PLASMA;
   else                   leastsquares_ptr = DGELS_LAPACK;
#else
   leastsquares_ptr = &DGELS_LAPACK;
   getFunc_ptr      = &GetFunction;
   getFunctionAndJacobian_ptr = &GetFunctionAndJacobian;
#endif

   // return if input is a satisfactory quadrature
   {
      getFunc_ptr(q_prev, RHS);
      double eNorm = V_InfNorm(RHS);
      if( (eNorm < q_tol) && (QuadInConstraint(q_prev) == true) ) {
         SOL_FLAG = SOL_FOUND;
         *its = 0;
         goto FREERETURN;
      }
   }


   while( (itsLoc < maxiter) && (errorNormUpdate > q_tol) )
   {

      if(CONSTR_OPT == ON  && cVectData.ACTIVE == ON  && cVectData.N_OR_W == WEIGHT)
      {
         if( (k == 1) || (++elim_weights > MAX_ELIM_WEIGHTS) ) {
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
         ncols = (dim+1)*--k;
         SMALL_DIM = MIN(nrows, ncols);
         LEAD_DIM  = MAX(nrows, ncols);
         CMatrix_realloc(nrows, ncols, &JACOBIAN);
         Vector_realloc(LEAD_DIM, &LEAST_SQ_SOL);
         quadrature_remove_element(cVectData.boundaryNodeId, q_next);
         quadrature_remove_element(cVectData.boundaryNodeId, q_prev);
      }

      getFunctionAndJacobian_ptr(q_prev, RHS, JACOBIAN); // computes RHS at essentially zero cost
      for(int i = 0; i < SMALL_DIM; ++i)        LEAST_SQ_SOL.id[i] = RHS.id[i];
      for(int i = SMALL_DIM; i < LEAD_DIM; ++i) LEAST_SQ_SOL.id[i] = 0.0;

      double start_time = get_cur_time();
      int INFO = leastsquares_ptr(JACOBIAN, LEAST_SQ_SOL);
      LSQ_TIME += get_cur_time() - start_time;

      for(int i = 0; i < ncols; ++i)
         q_next->z.id[i] = q_prev->z.id[i] - LEAST_SQ_SOL.id[i];

      errorNormPrev = errorNorm;
      getFunc_ptr(q_next, RHS);
      errorNorm = V_InfNorm(RHS);
      int check_values = CheckForStop(INFO, errorNorm, errorNormPrev, LEAST_SQ_SOL);
      if(check_values != GQ_SUCCESS) {
         SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }

      if(CONSTR_OPT == ON)
      {
         // Constrained optimization part 1: project onto the boundary
         int P_FLAG = ConstrainedProjection(q_prev, q_next);
         if(P_FLAG != CONSTR_SUCCESS) {
            COND_PRINT_ERR(P_FLAG);
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }

         // Constrained optimization part 2: shorten onto the boundary
         ConstrOptData *data = ConstrainedOptimizationInit(q_next);
         int C_FLAG = ConstrainedOptimization(data, q_prev, q_next, &cVectData);
         ConstrainedOptimizationFree(data);
         if(C_FLAG < 0) {
            PRINT_ERR("Constrained optimization failed", __LINE__, __FILE__);
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
      }
      else
      {
         if(!QuadInConstraint(q_next) && itsLoc > 10) {
            SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
      }
      getFunc_ptr(q_next, RHS);
      // if statement since norm might change after constrained optimization is applied
      if(CONSTR_OPT == ON) errorNormUpdate = V_InfNorm(RHS);
      else errorNormUpdate = errorNorm;

      quadrature_assign(q_next, q_prev);
      ++itsLoc;
   }


   if( (itsLoc < maxiter)
         && (isnan(errorNormUpdate) == 0)
         && (isinf(errorNormUpdate) == 0)
         && (errorNormUpdate <= q_tol)
         && QuadInConstraint(q_next) )
   {
      quadrature_assign_resize(q_next, q_orig);
      SOL_FLAG = SOL_FOUND;
   }
   else SOL_FLAG = SOL_NOT_FOUND;


FREERETURN:
#ifdef _OPENMP
   QuadFreeBasisOmp(q_next, omp_get_max_threads());
   QuadFreeBasisOmp(q_prev, omp_get_max_threads());
   FreeVectorOmpData(RHS);
#endif
   CMatrix_free(JACOBIAN);
   Vector_free(LEAST_SQ_SOL);
   free(WORK);
   Vector_free(RHS);
   quadrature_free(q_prev);
   quadrature_free(q_next);

   *its = itsLoc;
   return SOL_FLAG;
}// end LeastSquaresNewton



static int CheckForStop(int INFO, double errorNorm, double errorNormPrev, Vector least_sq_sol)
{
   if(errorNorm > errorNormPrev+2)     // fail if method is not converging
      return NOT_CONVERGE;

   if(V_InfNorm(least_sq_sol) > 100)   // fail if solution is exploding
      return DIVERGE_ERR;

   else if(INFO != 0)                  // fail if LAPACK routine has failed
      return LAPACK_ERR;

   else if(V_CheckInf(least_sq_sol))
      return INF_VAL;

   else if(V_CheckNan(least_sq_sol))
      return NAN_VAL;

   else return GQ_SUCCESS;
}

