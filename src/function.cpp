// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/function.hpp"

#include "../include/util.hpp"

namespace gquad {

double EvalFunction::function_time_total = 0;
double EvalJacobian::jacobian_time_total = 0;

////////////////////////////////////////////////////////////////////////////////////////////////////
// basis_array is initialized for either serial or parallel computation
EvalFunction::EvalFunction(Basis& basis, basis_type _type)
    : bases{},
      f(basis.size()),
      basis_array(basis),
      type{_type} {
   if(type == orthogonal) {
      functions_fptr = &Basis::orthog_basis;
      integrals_fptr = &Basis::orthog_integrals;
   } else if(type == monomial) {
      functions_fptr = &Basis::monomial_basis;
      integrals_fptr = &Basis::monomial_integrals;
   }
}

Array1D& EvalFunction::operator()(const QuadDomain& q) {
   util::timer t;
   eval_function(q);
   function_time_total += t.wall_time();
   return f;
}

double function_residual(const QuadDomain& q, Basis& basis, basis_type type) {
   EvalFunction eval_function{basis, type};
   auto z = eval_function(q);
   return norm_inf(z);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
EvalJacobian::EvalJacobian(Basis& basis, basis_type _type) : basis_array(basis), type{_type} {
   if(type == orthogonal) {
      functions_fptr = &Basis::orthog_basis;
      derivatives_fptr = &Basis::orthog_der;
   } else if(type == monomial) {
      functions_fptr = &Basis::monomial_basis;
      derivatives_fptr = &Basis::monomial_der;
   }
}

void EvalJacobian::operator()(const QuadDomain& q, Matrix2D& Jacobian) {
   util::timer t;
   eval_jacobian(q, Jacobian);
   jacobian_time_total += t.wall_time();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// calls either serial or parallel algorithm for computing function
void EvalFunction::eval_function(const QuadDomain& q) {
   auto& b0 = basis_array[0];
   GEN_QUAD_ASSERT_DEBUG(q.deg() == b0->deg());
   GEN_QUAD_ASSERT_DEBUG(q.dim() == b0->dim());

   bases.resize(f.size(), q.num_nodes());

   if(basis_array.compute_type == BasisArray::serial) {
      eval_function_serial(q);
   }
#ifdef _OPENMP
   else if(basis_array.compute_type == BasisArray::parallel) {
      eval_function_omp(q);
   }
#endif
}

void EvalFunction::eval_function_serial(const QuadDomain& q) {
   auto& b = basis_array[0];
   const auto& integrals = ((*b).*integrals_fptr)();
   f = -1. * integrals;

   for(gq_int k = 0; k < q.num_nodes(); ++k) {
      bases.col(k) = ((*b).*functions_fptr)(q.node(k));
   }

   for(gq_int i = 0; i < f.size(); ++i) {
      f[i] += stable_dot_prod(bases.row(i), q.weights());
   }
}

#ifdef _OPENMP
void EvalFunction::eval_function_omp(const QuadDomain& q) {
   gq_int nthreads = gq_int(basis_array.size());
   auto& b = basis_array[0];

   const auto& integrals = ((*b).*integrals_fptr)();
   f = -1. * integrals;

   Array2D bases(f.size(), q.num_nodes());
#pragma omp parallel default(none) shared(bases, q, f) num_threads(nthreads)
   {
      auto& bLoc = basis_array[gq_int(omp_get_thread_num())];
#pragma omp for schedule(static)
      for(gq_int k = 0; k < q.num_nodes(); ++k) {
         bases.col(k) = ((*bLoc).*functions_fptr)(q.node(k));
      }

#pragma omp for schedule(static)
      for(gq_int i = 0; i < f.size(); ++i) {
         f[i] += stable_dot_prod(bases.row(i), q.weights());
      }
   }  // end omp parallel
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
// calls either serial or parallel algorithm for computing Jacobian
void EvalJacobian::eval_jacobian(const QuadDomain& q, Matrix2D& J) {
   auto& b0 = basis_array[0];
   GEN_QUAD_ASSERT_DEBUG(q.deg() == b0->deg());
   GEN_QUAD_ASSERT_DEBUG(q.dim() == b0->dim());

   GEN_QUAD_ASSERT_DEBUG(J.cols() == q.size());
   GEN_QUAD_ASSERT_DEBUG(J.rows() == b0->size());

   if(basis_array.compute_type == BasisArray::serial) {
      eval_jacobian_serial(q, J);
   }
#ifdef _OPENMP
   else if(basis_array.compute_type == BasisArray::parallel) {
      eval_jacobian_omp(q, J);
   }
#endif
}

void EvalJacobian::eval_jacobian_serial(const QuadDomain& q, Matrix2D& J) {
   auto& b = basis_array[0];
   for(gq_int j = 0; j < q.num_nodes(); ++j) {
      const auto& functions = ((*b).*functions_fptr)(q.node(j));
      J.col(j) = functions;
   }

   gq_int dim{q.dim()};
   gq_int num_nodes{q.num_nodes()};
   for(gq_int j = 0; j < q.num_nodes(); ++j) {
      const auto& derivatives = ((*b).*derivatives_fptr)(q.node(j)) * q.w(j);
      J(all, seqN(num_nodes + j * dim, dim)) = derivatives.transpose();
   }
}

#ifdef _OPENMP
void EvalJacobian::eval_jacobian_omp(const QuadDomain& q, Matrix2D& J) {
   gq_int dim{q.dim()};
   gq_int num_nodes{q.num_nodes()};
   gq_int nthreads{basis_array.size()};

#pragma omp parallel default(shared) num_threads(nthreads)
   {
      auto& bLoc = basis_array[gq_int(omp_get_thread_num())];
#pragma omp for schedule(static)
      for(gq_int j = 0; j < num_nodes; ++j) {
         const auto& functions = ((*bLoc).*functions_fptr)(q.node(j));
         J.col(j) = functions;
      }

#pragma omp for schedule(static)
      for(gq_int j = 0; j < num_nodes; ++j) {
         const auto& derivatives = ((*bLoc).*derivatives_fptr)(q.node(j)) * q.w(j);
         J(all, seqN(num_nodes + j * dim, dim)) = derivatives.transpose();
      }
   }
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////
BasisArray::BasisArray(Basis& basis) {
#ifdef _OPENMP
   if(omp_condition(basis.deg(), basis.dim())) {
      compute_type = parallel;
      for(gq_int i = 0; i < omp_get_max_threads(); ++i) {
         components.push_back(basis.clone());
      }
   } else {
      compute_type = serial;
      components.push_back(basis.clone());
   }
#else
   compute_type = serial;
   components.push_back(basis.clone());
#endif
}

}  // namespace gquad

#endif

