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
EvalFunction::EvalFunction(Basis& basis, basis_type btype)
    : bases_{},
      f_(basis.size()),
      basis_array_(basis),
      btype_{btype} {
   if(btype_ == orthogonal) {
      functions_ptr_ = &Basis::orthog_basis;
      integrals_fptr_ = &Basis::orthog_integrals;
   } else if(btype_ == monomial) {
      functions_ptr_ = &Basis::monomial_basis;
      integrals_fptr_ = &Basis::monomial_integrals;
   }
}


Array1D& EvalFunction::operator()(const QuadDomain& q) {
   util::timer t;
   this->eval_function(q);
   function_time_total += t.wall_time();
   return f_;
}


double function_residual(const QuadDomain& q, Basis& basis, basis_type btype) {
   return norm_inf(EvalFunction{basis, btype}(q));
}


////////////////////////////////////////////////////////////////////////////////////////////////////
EvalJacobian::EvalJacobian(Basis& basis, basis_type btype) : basis_array_(basis), btype_{btype} {
   if(btype_ == orthogonal) {
      functions_ptr_ = &Basis::orthog_basis;
      derivatives_fptr_ = &Basis::orthog_der;
   } else if(btype_ == monomial) {
      functions_ptr_ = &Basis::monomial_basis;
      derivatives_fptr_ = &Basis::monomial_der;
   }
}


void EvalJacobian::operator()(const QuadDomain& q, Matrix2D& JT) {
   util::timer t;
   this->eval_jacobian(q, JT);
   jacobian_time_total += t.wall_time();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// calls either serial or parallel algorithm for computing function
void EvalFunction::eval_function(const QuadDomain& q) {
   auto& b0 = basis_array_[0];
   GEN_QUAD_ASSERT_DEBUG(q.deg() == b0->deg());
   GEN_QUAD_ASSERT_DEBUG(q.dim() == b0->dim());

   bases_.resize(f_.size(), q.num_nodes());

   if(basis_array_.compute_type_ == BasisArray::serial) {
      this->eval_function_serial(q);
   }
#ifdef _OPENMP
   else if(basis_array_.compute_type_ == BasisArray::parallel) {
      this->eval_function_omp(q);
   }
#endif
}


void EvalFunction::eval_function_serial(const QuadDomain& q) {
   auto& b = basis_array_[0];
   const auto& integrals = ((*b).*integrals_fptr_)();
   f_ = -1. * integrals;

   for(gq_int k = 0; k < q.num_nodes(); ++k) {
      bases_.col(k) = ((*b).*functions_ptr_)(q.node(k));
   }

   for(gq_int i = 0; i < f_.size(); ++i) {
      f_[i] += stable_dot_prod(bases_.row(i), q.weights());
   }
}


#ifdef _OPENMP
void EvalFunction::eval_function_omp(const QuadDomain& q) {
   gq_int nthreads = gq_int(basis_array_.size());
   auto& b = basis_array_[0];

   const auto& integrals = ((*b).*integrals_fptr_)();
   f_ = -1. * integrals;

   Array2D bases_(f_.size(), q.num_nodes());
#pragma omp parallel default(none) shared(bases_, q, f_) num_threads(nthreads)
   {
      auto& bLoc = basis_array_[gq_int(omp_get_thread_num())];
#pragma omp for schedule(static)
      for(gq_int k = 0; k < q.num_nodes(); ++k) {
         bases_.col(k) = ((*bLoc).*functions_ptr_)(q.node(k));
      }

#pragma omp for schedule(static)
      for(gq_int i = 0; i < f_.size(); ++i) {
         f_[i] += stable_dot_prod(bases_.row(i), q.weights());
      }
   }  // end omp parallel
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////
// calls either serial or parallel algorithm for computing Jacobian
void EvalJacobian::eval_jacobian(const QuadDomain& q, Matrix2D& JT) {
   auto& b0 = basis_array_[0];
   GEN_QUAD_ASSERT_DEBUG(q.deg() == b0->deg());
   GEN_QUAD_ASSERT_DEBUG(q.dim() == b0->dim());

   GEN_QUAD_ASSERT_DEBUG(JT.cols() == q.size());
   GEN_QUAD_ASSERT_DEBUG(JT.rows() == b0->size());

   if(basis_array_.compute_type_ == BasisArray::serial) {
      this->eval_jacobian_serial(q, JT);
   }
#ifdef _OPENMP
   else if(basis_array_.compute_type_ == BasisArray::parallel) {
      this->eval_jacobian_omp(q, JT);
   }
#endif
}


void EvalJacobian::eval_jacobian_serial(const QuadDomain& q, Matrix2D& JT) {
   auto& b = basis_array_[0];
   for(gq_int j = 0; j < q.num_nodes(); ++j) {
      const auto& functions = ((*b).*functions_ptr_)(q.node(j));
      JT.col(j) = functions;
   }

   gq_int dim{q.dim()};
   gq_int num_nodes{q.num_nodes()};
   for(gq_int j = 0; j < q.num_nodes(); ++j) {
      const auto& derivatives = ((*b).*derivatives_fptr_)(q.node(j)) * q.w(j);
      JT(all, seqN(num_nodes + j * dim, dim)) = derivatives.transpose();
   }
}


#ifdef _OPENMP
void EvalJacobian::eval_jacobian_omp(const QuadDomain& q, Matrix2D& JT) {
   gq_int dim{q.dim()};
   gq_int num_nodes{q.num_nodes()};
   gq_int nthreads{basis_array_.size()};

#pragma omp parallel default(shared) num_threads(nthreads)
   {
      auto& bLoc = basis_array_[gq_int(omp_get_thread_num())];
#pragma omp for schedule(static)
      for(gq_int j = 0; j < num_nodes; ++j) {
         const auto& functions = ((*bLoc).*functions_ptr_)(q.node(j));
         JT.col(j) = functions;
      }

#pragma omp for schedule(static)
      for(gq_int j = 0; j < num_nodes; ++j) {
         const auto& derivatives = ((*bLoc).*derivatives_fptr_)(q.node(j)) * q.w(j);
         JT(all, seqN(num_nodes + j * dim, dim)) = derivatives.transpose();
      }
   }
}
#endif


////////////////////////////////////////////////////////////////////////////////////////////////////
BasisArray::BasisArray(Basis& basis) {
#ifdef _OPENMP
   if(omp_condition(basis.deg(), basis.dim())) {
      compute_type_ = parallel;
      for(gq_int i = 0; i < omp_get_max_threads(); ++i) {
         components_.push_back(basis.clone());
      }
   } else {
      compute_type_ = serial;
      components_.push_back(basis.clone());
   }
#else
   compute_type_ = serial;
   components_.push_back(basis.clone());
#endif
}


}  // namespace gquad

#endif

