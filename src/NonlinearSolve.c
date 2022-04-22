/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#include "NonlinearSolve.h"
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


static double ComputePenalty(const quadrature *q_prev, const quadrature *q_next);
static errInfo CheckForStop(int INFO, double errorNorm, double errorNormPrev, const Vector least_sq_sol);

double LSQ_TIME = 0.0;
#define MAX_ELIM_WEIGHTS 10

LSQ_out LeastSquaresNewton(const bool_enum CONSTR_OPT, quadrature *q_orig)
{
   assert(q_orig->num_nodes >= 1);
   LSQ_out lsq_out = {0, SOL_NOT_FOUND};
   if(q_orig->isFullyInitialized != GQ_TRUE) {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return lsq_out;
   }
   if(V_IsUninitialized(q_orig->z))
      PRINT_ERR(STR_INV_INPUT, __LINE__, __FILE__);
   if(V_InfNorm(q_orig->z) >= QUAD_HUGE) {
      PRINT_ERR(STR_QUAD_HUGE_ERR, __LINE__, __FILE__);
      return lsq_out;
   }

   int dim       = q_orig->dim;
   int numFuncs  = q_orig->basis->numFuncs;
   int num_nodes = q_orig->num_nodes;

   int maxiter  = 25;
   double q_tol = QUAD_TOL; // 10^(-14);

   int nrows    = numFuncs;
   int ncols    = (dim+1)*num_nodes;
   CMatrix JACOBIAN = CMatrix_init(nrows, ncols);
   Vector dz        = Vector_init(ncols);
   Vector FPrev     = Vector_init(nrows);
   Vector FNext     = Vector_init(nrows);

   quadrature *q_prev = quadrature_make_full_copy(q_orig);
   quadrature *q_next = quadrature_make_full_copy(q_orig);
   Vector *z_prev = &q_prev->z;
   Vector *z_next = &q_next->z;

   int elim_weights = 0;
   double errorNorm = 1.0, errorNormPrev = 1.0;
   ConstrVectData cVectData = ConstrVectDataInit();

   int (*leastsquares_ptr)(CMatrix A, Vector b, Vector x);
   void (*getFunc_ptr)(quadrature *q, Vector f);
   void (*getFunctionAndJacobian_ptr)(quadrature *q, Vector f, CMatrix JACOBIAN);
#ifdef _OPENMP
   bool is_alloc_omp = false;
   if(OMP_CONDITION(q_orig->deg, dim))
   {
      QuadAllocBasisOmp(q_next, omp_get_max_threads());
      QuadAllocBasisOmp(q_prev, omp_get_max_threads());
      AllocVectorOmpData(&FPrev, omp_get_max_threads());
      AllocVectorOmpData(&FNext, omp_get_max_threads());
      is_alloc_omp = true;

      leastsquares_ptr = DGELS_PLASMA;
      getFunc_ptr = &GetFunctionOmp;
      getFunctionAndJacobian_ptr = &GetFunctionAndJacobianOmp;
   }
   else
   {
      leastsquares_ptr = DGELS_LAPACK;
      getFunc_ptr = &GetFunction;
      getFunctionAndJacobian_ptr = &GetFunctionAndJacobian;
   }
#else
   leastsquares_ptr = &DGELS_LAPACK;
   getFunc_ptr      = &GetFunction;
   getFunctionAndJacobian_ptr = &GetFunctionAndJacobian;
#endif


   // return if input is a satisfactory quadrature
   getFunc_ptr(q_prev, FPrev);
   double errorNormUpdate = V_InfNorm(FPrev);
   if( (errorNormUpdate <= q_tol) && (QuadInConstraint(q_next) == true) )
   {
      lsq_out.SOL_FLAG = SOL_FOUND;
      lsq_out.its = 0;
      goto FREERETURN;
   }


   while( (lsq_out.its < maxiter) && (errorNormUpdate > q_tol) )                                     // while infinity norm > q_tol
   {
      if( (CONSTR_OPT == ON) && (cVectData.N_OR_W == WEIGHT) )                                       // constrained optimization : remove 1 node from  q_prev and q_next, preserve other nodes,
                                                                                                     // remove (dim+1) columns from JACOBIAN, remove (dim+1) entries from dz vector
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
            CMatrix_realloc(nrows, ncols, &JACOBIAN);
            Vector_realloc(ncols, &dz);
            quadrature_remove_element(cVectData.boundaryNodeId, q_prev);
            quadrature_remove_element(cVectData.boundaryNodeId, q_next);
         }
      }

      getFunctionAndJacobian_ptr(q_prev, FPrev, JACOBIAN);                                           // F(q_prev)

      double start_time = get_cur_time();
      int INFO = leastsquares_ptr(JACOBIAN, FPrev, dz);                                              // dz = J(q_prev) \ F(q_prev) (pseudoinverse)
      LSQ_TIME += get_cur_time() - start_time;

      if(INFO != 0)                                                                                  // \\\\\ return if could not solve the system /////
      {
         lsq_out.SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }

      VScale(-1.0, dz);
      VAdd(*z_prev, dz, *z_next);                                                                    // z_next = z_prev + dz

      errorNormPrev = errorNormUpdate;
      getFunc_ptr(q_next, FNext);
      errorNorm = V_InfNorm(FNext);
      errInfo check_values = CheckForStop(INFO, errorNorm, errorNormPrev, dz);                       // \\\\\ return if certain criteria is not satisfied /////
      if(check_values.succeeded != true) {
         lsq_out.SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }

      if(CONSTR_OPT == OFF && !QuadInConstraint(q_next)) {                                           // \\\\\ return if not inside the boundary in unconstrained mode /////
         lsq_out.SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }

      if(CONSTR_OPT == OFF && lsq_out.its < 5) {                                                     // z_next = z_prev + α * dz
         double alpha = ComputePenalty(q_prev, q_next);
         if(alpha < 1.0)
            VAddScale(1.0, *z_prev, alpha, dz, *z_next);
      }

      if(CONSTR_OPT == ON)                                                                           // project and shorten in constrained mode, update q_next === z_next if necessary
      {
         // Constrained optimization part 1: project onto the boundary
         int P_FLAG = ConstrainedProjection(q_prev, q_next);                                         // \\\\\\ return if failed constrained projection /////
         if(P_FLAG < 0) {
            COND_PRINT_ERR(CONSTR_ERROR);
            lsq_out.SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }

         // Constrained optimization part 2: shorten onto the boundary
         ConstrOptData *data = ConstrainedOptimizationInit(q_next);                                  // \\\\\ return if did not shorten to the boundary /////
         int C_FLAG = ConstrainedOptimization(data, q_prev, q_next, &cVectData);
         ConstrainedOptimizationFree(data);
         if(C_FLAG < 0) {
            COND_PRINT_ERR(CONSTR_ERROR);
            lsq_out.SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }
      }

      getFunc_ptr(q_next, FNext);                                                                    // F(q_next)
      // if statement since norm might change after constrained optimization is applied
      if(CONSTR_OPT == ON) errorNormUpdate = V_InfNorm(FNext);
      else errorNormUpdate = errorNorm;

      Vector_Assign(*z_next, *z_prev);                                                               // z_prev = z_next
      ++lsq_out.its;
   }


   if( (lsq_out.its <= maxiter)
    && (isnan(errorNormUpdate) == false)
    && (isinf(errorNormUpdate) == false)
    && (errorNormUpdate <= q_tol)
    && QuadInConstraint(q_next) )
   {
      quadrature_assign_resize(q_next, q_orig);
      lsq_out.SOL_FLAG = SOL_FOUND;
   }
   else lsq_out.SOL_FLAG = SOL_NOT_FOUND;


FREERETURN:
//   if(lsq_out.SOL_FLAG == SOL_NOT_FOUND)
//   {
//      PrintInt(lsq_out.its, "its_when_failed");
//      PrintDouble(errorNormUpdate, "errorNormUpdate");
//      PrintBool(QuadInDomain(q_next), "InDomain");
//   }
#ifdef _OPENMP
   if(is_alloc_omp)
   {
      FreeVectorOmpData(FPrev);
      FreeVectorOmpData(FNext);
      QuadFreeBasisOmp(q_prev);
      QuadFreeBasisOmp(q_next);
   }
#endif
   CMatrix_free(JACOBIAN);
   Vector_free(dz);
   Vector_free(FPrev);
   Vector_free(FNext);
   quadrature_free(q_prev);
   quadrature_free(q_next);

   return lsq_out;
}// end LeastSquaresNewton




LSQ_out LevenbergMarquardt(const bool_enum CONSTR_OPT, quadrature *q_orig)
{
   assert(q_orig->num_nodes >= 1);
   LSQ_out lsq_out = {0, SOL_NOT_FOUND};
   if(q_orig->isFullyInitialized != GQ_TRUE) {
      PRINT_ERR(STR_QUAD_NOT_FULL_INIT, __LINE__, __FILE__);
      return lsq_out;
   }
   if(V_IsUninitialized(q_orig->z))
      PRINT_ERR(STR_INV_INPUT, __LINE__, __FILE__);
   if(V_InfNorm(q_orig->z) >= QUAD_HUGE) {
      PRINT_ERR(STR_QUAD_HUGE_ERR, __LINE__, __FILE__);
      return lsq_out;
   }

   const int dim       = q_orig->dim;
   const int num_nodes = q_orig->num_nodes;
   const int numFuncs  = q_orig->basis->numFuncs;

   const int maxiter  = 75;
   const double q_tol = QUAD_TOL; // 10^(-14);

   const int nrows    = numFuncs;
   const int ncols    = (dim+1)*num_nodes;

   CMatrix JACOBIAN     = CMatrix_init(nrows, ncols);
   CMatrix JACOBIAN_TR  = CMatrix_init(ncols, nrows);
   CMatrix JT_J_lmd     = CMatrix_init(ncols, ncols);

   double alpha_lvmr = 0.;
   if(CONSTR_OPT == OFF)     alpha_lvmr = 1.;                    // more local initial guess, consider other values
   else if(CONSTR_OPT == ON) alpha_lvmr = 100.;                  // more global search when solution is harder to find, consider other values

   Vector FPrev     = Vector_init(nrows);
   Vector FCur      = Vector_init(nrows);
   Vector FNext     = Vector_init(nrows);
   Vector LevMarRHS = Vector_init(ncols);
   Vector dz        = Vector_init(ncols);

   quadrature *q_prev = quadrature_make_full_copy(q_orig);
   quadrature *q_cur  = quadrature_make_full_copy(q_orig);
   quadrature *q_next = quadrature_make_full_copy(q_orig);
   Vector *z_prev = &q_prev->z;
   Vector *z_cur  = &q_cur->z;
   Vector *z_next = &q_next->z;

   int (*leastsquares_ptr)(CMatrix A, Vector b, Vector x);
   void (*getFunc_ptr)(quadrature *q, Vector f);
   void (*getFunctionAndJacobian_ptr)(quadrature *q, Vector f, CMatrix JACOBIAN);
   int (*dgemm_ptr)(CMatrix A, CMatrix B, CMatrix C);

#ifdef _OPENMP
   bool is_alloc_omp = false;
   if(OMP_CONDITION(q_orig->deg, dim))
   {
      QuadAllocBasisOmp(q_prev, omp_get_max_threads());
      QuadAllocBasisOmp(q_cur, omp_get_max_threads());
      QuadAllocBasisOmp(q_next, omp_get_max_threads());
      AllocVectorOmpData(&FPrev, omp_get_max_threads());
      AllocVectorOmpData(&FCur, omp_get_max_threads());
      AllocVectorOmpData(&FNext, omp_get_max_threads());
      is_alloc_omp = true;

      leastsquares_ptr           = DGELS_PLASMA;
      getFunc_ptr                = &GetFunctionOmp;
      getFunctionAndJacobian_ptr = &GetFunctionAndJacobianOmp;
      dgemm_ptr                  = &DGEMM_PLASMA;
   }
   else
   {
      leastsquares_ptr           = DGELS_LAPACK;
      getFunc_ptr                = &GetFunction;
      getFunctionAndJacobian_ptr = &GetFunctionAndJacobian;
      dgemm_ptr                  = &DGEMM_LAPACK;
   }
#else
   leastsquares_ptr           = &DGELS_LAPACK;
   getFunc_ptr                = &GetFunction;
   getFunctionAndJacobian_ptr = &GetFunctionAndJacobian;
   dgemm_ptr                  = &DGEMM_LAPACK;
#endif

   // return if input is a satisfactory quadrature
   getFunc_ptr(q_prev, FPrev);
   double errorNormUpdate = V_InfNorm(FPrev);
   if( (errorNormUpdate <= q_tol) && (QuadInConstraint(q_prev) == true) )
   {
      lsq_out.SOL_FLAG = SOL_FOUND;
      lsq_out.its = 0;
      goto FREERETURN;
   }


   while( (lsq_out.its < maxiter) && (errorNormUpdate > q_tol) )  // while infinity norm > q_tol
   {
      double alpha_lvmr_down;
      if( errorNormUpdate > POW_DOUBLE(10, -8) )
         alpha_lvmr_down = alpha_lvmr / 3.0;                                                         // consider dividing by 5
      else alpha_lvmr_down = alpha_lvmr / 5.;

      getFunctionAndJacobian_ptr(q_prev, FPrev, JACOBIAN);                                           // J(xp)                            m x n
      CMatrix_Assign_Transpose(JACOBIAN, JACOBIAN_TR);                                               // JT(xp)                           n x m

      double start_time = get_cur_time();                                                            // include dgemm in LSQ_TIME
      int INFO = 0;
      if( errorNormUpdate > POW_DOUBLE(10, -10) )
      {

         dgemm_ptr(JACOBIAN_TR, JACOBIAN, JT_J_lmd);                                                 // JT_J_lmd = JT * J                n x n
         for(int i = 0; i < JT_J_lmd.rows; ++i)
            JT_J_lmd.cid[i][i] += alpha_lvmr * JT_J_lmd.cid[i][i];                                   // JT_J_lmd += λ*D                  n x n
         CMatVecAccurate(JACOBIAN_TR, FPrev, LevMarRHS);                                             // LevMarRHS = JT * F               n x m * m -> n
         VScale(-1.0, LevMarRHS);
         INFO = leastsquares_ptr(JT_J_lmd, LevMarRHS, dz);                                           // dz = -(JT * J + λD) \ (JT * F)   n x n * n -> n
      }
      else                                                                                           // apply pseudoinverse directly when error is small enough
      {
         INFO = leastsquares_ptr(JACOBIAN, FPrev, dz);
         VScale(-1.0, dz);

      }
      LSQ_TIME += get_cur_time() - start_time;

      if(INFO != 0)                                                                                  // \\\\\ return if could not solve the system /////
      {
         lsq_out.SOL_FLAG = SOL_NOT_FOUND;
         goto FREERETURN;
      }

      VAdd(*z_prev, dz, *z_cur);                                                                     // z_cur = z_prev + dz -> q_cur = q_prev + dz

      // |||||||||| start region for updating q_next === z_next and alpha_lvmr ||||||||||
      getFunc_ptr(q_cur, FCur);                                                                      // F(q_cur) === F(z_cur)
      if( ( V_InfNorm(FCur) > V_InfNorm(FPrev) ) || (INFO != GQ_SUCCESS) ) {                         // if F(q_cur) > F(q_prev) : alpha *= 2, keep z_next = z_prev
         alpha_lvmr = MIN(alpha_lvmr*2.0, POW_DOUBLE(10.0, 7));
      }
      else {                                                                                         // else if F(q_cur) < F(q_prev) : alpha /= 3, update z_next
         alpha_lvmr = MAX(alpha_lvmr_down, POW_DOUBLE(10.0, -9));
         Vector_Assign(*z_cur, *z_next);                                                             // z_next = z_cur, z_next might be updated below

         if( CONSTR_OPT == OFF && !QuadInConstraint(q_next) )                                        // \\\\\ return if not inside the boundary in UNCONSTRAINED mode /////
         {
            lsq_out.SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;

         }
         else if( CONSTR_OPT == ON && !QuadInConstraint(q_next) && lsq_out.its >= 30 ) {             // \\\\\ return if not inside the boundary after 30 iterations in CONSTRAINED mode /////
            lsq_out.SOL_FLAG = SOL_NOT_FOUND;
            goto FREERETURN;
         }

         if(CONSTR_OPT == ON)
         {
            if(!QuadInConstraint(q_next))
            {
               ConstrVectData cVectData = ConstrVectDataInit();
               ConstrOptData *data = ConstrainedOptimizationInit(q_next);
               int C_FLAG = ConstrainedOptimization(data, q_prev, q_next, &cVectData);               // q_next -> q_next_b === z_next -> z_next_b
               ConstrainedOptimizationFree(data);
               if(C_FLAG < 0) {                                                                      // \\\\\ return if did not shorten to the boundary /////
                  COND_PRINT_ERR(CONSTR_ERROR);
                  lsq_out.SOL_FLAG = SOL_NOT_FOUND;
                  goto FREERETURN;
               }
               else if(C_FLAG == CONSTR_SUCCESS)
               {
                  Vector dz_short = Vector_init(ncols);
                  VAddScale(1.0, *z_next, -1.0, *z_prev, dz_short);                                  // dz_short_b = z_next_b - z_prev
                  VAddScale(1.0, *z_prev, 0.5, dz_short, *z_next);                                   // z_next = z_prev + 0.5 * dz_short_b
                  Vector_free(dz_short);
               }
            }
         }

         if(lsq_out.its < 25) {                                                                      // add boundary penalty
            double alpha = ComputePenalty(q_prev, q_next);                                           //
            if(alpha < 1.0)                                                                          // alpha >= 0
               VAddScale(1.0, *z_prev, alpha, dz, *z_next);                                          // z_next = z_prev + alpha * dz
         }
      }
      // |||||||||| end region for updating q_next === z_next and alpha_lvmr ||||||||||

      getFunc_ptr(q_next, FNext);                                                                    // F(q_next) === F(z_next)
      errorNormUpdate = V_InfNorm(FNext);

      Vector_Assign(*z_next, *z_prev);                                                               // z_prev = z_next
      ++lsq_out.its;
   }


   if( (lsq_out.its <= maxiter)
    && (isnan(errorNormUpdate) == false)
    && (isinf(errorNormUpdate) == false)
    && (errorNormUpdate <= q_tol)
    && QuadInConstraint(q_next) )
   {
      quadrature_assign_resize(q_next, q_orig);
      lsq_out.SOL_FLAG = SOL_FOUND;
   }
   else lsq_out.SOL_FLAG = SOL_NOT_FOUND;


FREERETURN:
//   if(lsq_out.SOL_FLAG == SOL_NOT_FOUND)
//   {
//      PrintInt(lsq_out.its, "its_when_failed");
//      PrintDouble(errorNormUpdate, "errorNormUpdate");
//      PrintBool(QuadInDomain(q_next), "InDomain");
//   }
#ifdef _OPENMP
   if(is_alloc_omp)
   {
      FreeVectorOmpData(FPrev);
      FreeVectorOmpData(FCur);
      FreeVectorOmpData(FNext);
      QuadFreeBasisOmp(q_prev);
      QuadFreeBasisOmp(q_cur);
      QuadFreeBasisOmp(q_next);
   }
#endif
   CMatrix_free(JACOBIAN);
   CMatrix_free(JACOBIAN_TR);
   CMatrix_free(JT_J_lmd);

   Vector_free(FPrev);
   Vector_free(FCur);
   Vector_free(FNext);
   Vector_free(LevMarRHS);
   Vector_free(dz);

   quadrature_free(q_prev);
   quadrature_free(q_cur);
   quadrature_free(q_next);

   return lsq_out;
}// end LeastSquaresNewton

#undef MAX_ELIM_WEIGHTS


static errInfo CheckForStop(int INFO, double errorNorm, double errorNormPrev, const Vector least_sq_sol)
{
   errInfo info;
   info.succeeded = true;
   info.errors[0] = GQ_SUCCESS;
   info.errors[1] = GQ_SUCCESS;
   info.errors[2] = GQ_SUCCESS;
   info.errors[3] = GQ_SUCCESS;

   if(INFO != 0) {                       // fail if LAPACK/PLASMA routine did not return success
      info.succeeded = false;
      info.errors[0] = LAPACK_ERR;
   }
   if(errorNorm > 10*errorNormPrev) {    // fail if method is not converging
      info.succeeded = false;
      info.errors[1] =  NOT_CONVERGE;
   }
   if(V_InfNorm(least_sq_sol) > 10) {    // fail if solution is exploding
      info.succeeded = false;
      info.errors[2] = DIVERGE_ERR;
   }

   if(V_CheckInf(least_sq_sol)) {
      info.succeeded = false;
      info.errors[3] = INF_VAL;
      PRINT_ERR(STR_INF_VAL, __LINE__, __FILE__);
   }
   if(V_CheckNan(least_sq_sol)) {
      info.succeeded = false;
      info.errors[3] = NAN_VAL;
      PRINT_ERR(STR_NAN_VAL, __LINE__, __FILE__);
   }

   return info;
}


static double ComputePenalty(const quadrature *q_prev, const quadrature *q_next)
{
   if(!QuadInConstraint(q_prev)) return 1.0;
   if(!QuadInConstraint(q_next)) return 1.0;

   double distNext = QuadMinDistFromTheBoundary(q_next);
   double distPrev = QuadMinDistFromTheBoundary(q_prev);
#ifdef QUAD_DEBUG_ON
   char errString[100] = STR_BOUND_ERROR;
   if(distNext < 0 || distPrev < 0)
      PRINT_ERR(strcat(errString, " , most likely a bug in QuadMinDistFromTheBoundary"), __LINE__, __FILE__);
#endif

   if(distNext > distPrev) return 1.0;
   else                    return distNext / distPrev;
}


