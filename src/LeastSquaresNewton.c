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

typedef struct
{
   bool succeeded;
   int errors[4];
} errInfo;


static double ComputePenalty(quadrature *q_prev, quadrature *q_next);
static errInfo CheckForStop(int INFO, double errorNorm, double errorNormPrev, Vector least_sq_sol);

double CONSTR_TIME = 0.0;
double LSQ_TIME = 0.0;
#define MAX_ELIM_WEIGHTS 10

LSQ_out LeastSquaresNewton(const bool_enum CONSTR_OPT, quadrature *q_orig)
{
   assert(q_orig->num_nodes >= 1);
   LSQ_out lsq_out = {0, SOL_NOT_FOUND};

   int dim       = q_orig->dim;
   int num_nodes = q_orig->num_nodes;
   int numFuncs  = q_orig->basis->numFuncs;

   int maxiter  = 25;
   double q_tol = QUAD_TOL; // 10^(-14);

   int nrows    = numFuncs;
   int ncols    = (dim+1)*num_nodes;
   int LEAD_DIM = MAX(nrows, ncols);
   CMatrix JACOBIAN    = CMatrix_init(nrows, ncols);
   Vector LEAST_SQ_SOL = Vector_init(LEAD_DIM);
   Vector RHS          = Vector_init(nrows);

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
   getFunc_ptr(q_prev, RHS);
   double eNorm = V_InfNorm(RHS);
   if( (eNorm < q_tol) && (QuadInConstraint(q_prev) == true) )
   {
      lsq_out.SOL_FLAG = SOL_FOUND;
      lsq_out.its = 0;
      goto FREERETURN;
   }


   while( (lsq_out.its < maxiter) && (errorNormUpdate > q_tol) )
   {

      if( (CONSTR_OPT == ON) && (cVectData.N_OR_W == WEIGHT) )
      {
         if( q_next->w[cVectData.boundaryNodeId] > POW_DOUBLE(10, -12) )
            PRINT_ERR("weight was not shortened to 0", __LINE__, __FILE__);

         if( (num_nodes == 1) || (++elim_weights > MAX_ELIM_WEIGHTS) )
         {
            lsq_out.SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
         else
         {
            ncols = (dim+1) * --num_nodes;
            LEAD_DIM = MAX(nrows, ncols);
            CMatrix_realloc(nrows, ncols, &JACOBIAN);
            Vector_realloc(LEAD_DIM, &LEAST_SQ_SOL);
            quadrature_remove_element(cVectData.boundaryNodeId, q_next);
            quadrature_remove_element(cVectData.boundaryNodeId, q_prev);
         }
      }

      getFunctionAndJacobian_ptr(q_prev, RHS, JACOBIAN); // computes RHS at essentially zero cost
      for(int i = 0; i < LEAST_SQ_SOL.len; ++i) LEAST_SQ_SOL.id[i] = 0.0;
      for(int i = 0; i < RHS.len; ++i)          LEAST_SQ_SOL.id[i] = RHS.id[i];

      double start_time = get_cur_time();
      int INFO = leastsquares_ptr(JACOBIAN, LEAST_SQ_SOL);
      LSQ_TIME += get_cur_time() - start_time;


      for(int i = 0; i < ncols; ++i)
         q_next->z.id[i] = q_prev->z.id[i] - LEAST_SQ_SOL.id[i];

      double alpha = 1.0;
      if(CONSTR_OPT == OFF && lsq_out.its < 5)
         alpha = ComputePenalty(q_prev, q_next);
      if(alpha < 1.0)
         for(int i = 0; i < ncols; ++i)
            q_next->z.id[i] = q_prev->z.id[i] - alpha * LEAST_SQ_SOL.id[i];

      errorNormPrev = errorNorm;
      getFunc_ptr(q_next, RHS);
      errorNorm = V_InfNorm(RHS);
      errInfo check_values = CheckForStop(INFO, errorNorm, errorNormPrev, LEAST_SQ_SOL);
      if(check_values.succeeded != true)
      {
         lsq_out.SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }

      if(!QuadInConstraint(q_next) && lsq_out.its >= 7) {
         lsq_out.SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }

      if(CONSTR_OPT == ON)
      {
         // Constrained optimization part 1: project onto the boundary
         int P_FLAG = ConstrainedProjection(q_prev, q_next);
         if(P_FLAG < 0) {
            COND_PRINT_ERR(CONSTR_ERROR);
            lsq_out.SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }

         // Constrained optimization part 2: shorten onto the boundary
         ConstrOptData *data = ConstrainedOptimizationInit(q_next);
         int C_FLAG = ConstrainedOptimization(data, q_prev, q_next, &cVectData);
         ConstrainedOptimizationFree(data);
         if(C_FLAG < 0) {
            COND_PRINT_ERR(CONSTR_ERROR);
            lsq_out.SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }

      }

      getFunc_ptr(q_next, RHS);
      // if statement since norm might change after constrained optimization is applied
      if(CONSTR_OPT == ON) errorNormUpdate = V_InfNorm(RHS);
      else errorNormUpdate = errorNorm;

      quadrature_assign(q_next, q_prev);
      ++lsq_out.its;
   }


   if( (lsq_out.its < maxiter)
    && (isnan(errorNormUpdate) == 0)
    && (isinf(errorNormUpdate) == 0)
    && (errorNormUpdate <= q_tol)
    && QuadInConstraint(q_next) )
   {
      quadrature_assign_resize(q_next, q_orig);
      lsq_out.SOL_FLAG = SOL_FOUND;
   }
   else lsq_out.SOL_FLAG = SOL_NOT_FOUND;


FREERETURN:
#ifdef _OPENMP
   QuadFreeBasisOmp(q_next, omp_get_max_threads());
   QuadFreeBasisOmp(q_prev, omp_get_max_threads());
   FreeVectorOmpData(RHS);
#endif
   CMatrix_free(JACOBIAN);
   Vector_free(LEAST_SQ_SOL);
   Vector_free(RHS);
   quadrature_free(q_prev);
   quadrature_free(q_next);

   return lsq_out;
}// end LeastSquaresNewton

#undef MAX_ELIM_WEIGHTS


static double ComputePenalty(quadrature *q_prev, quadrature *q_next)
{
   if(!QuadInConstraint(q_prev)) return 1.0;
   if(!QuadInConstraint(q_next)) return 1.0;

   double distNext = QuadMinDistFromTheBoundary(q_next);
   double distPrev = QuadMinDistFromTheBoundary(q_prev);

   if(distNext > distPrev) return 1.0;
   else                    return distNext / distPrev;
}


static errInfo CheckForStop(int INFO, double errorNorm, double errorNormPrev, Vector least_sq_sol)
{
   errInfo info;
   info.succeeded = true;
   info.errors[0] = GQ_SUCCESS;
   info.errors[1] = GQ_SUCCESS;
   info.errors[2] = GQ_SUCCESS;
   info.errors[3] = GQ_SUCCESS;

   if(INFO != 0) {                       // fail if LAPACK routine has failed
      info.succeeded = false;
      info.errors[0] = LAPACK_ERR;
   }
   if(errorNorm > errorNormPrev+4) {     // fail if method is not converging
      info.succeeded = false;
      info.errors[1] =  NOT_CONVERGE;
   }
   if(V_InfNorm(least_sq_sol) > 100) {   // fail if solution is exploding
      info.succeeded = false;
      info.errors[2] = DIVERGE_ERR;
   }

   if(V_CheckInf(least_sq_sol)) {
      info.succeeded = false;
      info.errors[3] = INF_VAL;
   }
   if(V_CheckNan(least_sq_sol)) {
      info.succeeded = false;
      info.errors[3] = NAN_VAL;
   }

   return info;
}

