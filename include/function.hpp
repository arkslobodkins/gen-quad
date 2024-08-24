// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

// This header contains definition of EvalFunction and EvalJacobian function objects,
// which implement functionality for both orthogonal and monomial versions of basis.
// If basis type is not specified, orthogonal basis is computed by default.
// If problem size is sufficiently large, computation is parallelized using OpenMP.
#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "basis.hpp"
#include "quadrature.hpp"

namespace gquad {


struct FunctionArray {
private:
   StdVector<Array1D> components_;
   Array1D root_;

public:
   enum { serial, parallel } compute_type_;

   explicit FunctionArray(Basis& basis);

   FunctionArray(const FunctionArray&) = delete;
   FunctionArray& operator=(const FunctionArray&) = delete;

   gq_int size() const {
      return gq_int(components_.size());
   }

   auto& operator[](gq_int i) {
      return components_.at(gq_int(i));
   }

   Array1D& data() {
      return root_;
   }
};


class BasisArray {
private:
   StdVector<std::unique_ptr<Basis>> components_;

public:
   enum { serial, parallel } compute_type_;

   explicit BasisArray(Basis& basis);

   BasisArray(const BasisArray&) = delete;
   BasisArray& operator=(const BasisArray&) = delete;

   gq_int size() const {
      return gq_int(components_.size());
   }

   auto& operator[](gq_int i) {
      return components_.at(gq_int(i));
   }
};


enum basis_type { monomial, orthogonal };


class EvalFunction {
public:
   static double function_time_total;

   explicit EvalFunction(Basis& basis, basis_type btype);

   EvalFunction(const EvalFunction&) = delete;
   EvalFunction& operator=(const EvalFunction&) = delete;

   // Operator that computes function for the Least Squares Newton.
   Array1D& operator()(const QuadDomain& q);

private:
   using FunctionsFptr = const Array1D& (Basis::*)(const Array1D& x);
   using IntegralsFptr = const Array1D& (Basis::*)();

   Array2D bases_;
   Array1D f_;
   BasisArray basis_array_;
   basis_type btype_;
   FunctionsFptr functions_ptr_{nullptr};
   IntegralsFptr integrals_fptr_{nullptr};

   void eval_function(const QuadDomain& q);
   void eval_function_serial(const QuadDomain& q);
#ifdef _OPENMP
   void eval_function_omp(const QuadDomain& q);
#endif
};


double function_residual(const QuadDomain& q, Basis& basis, basis_type btype);


class EvalJacobian {
public:
   static double jacobian_time_total;

   explicit EvalJacobian(Basis& basis, basis_type _type);

   EvalJacobian(const EvalFunction&) = delete;
   EvalJacobian& operator=(const EvalJacobian&) = delete;

   // Operator that computes transpose of the Jacobian for the Least Squares Newton.
   // Since Matrix is stored in row-major format, computing transpose
   // rather than Jacobian is more efficient.
   void operator()(const QuadDomain& q, Matrix2D& JT);

private:
   using FunctionsFptr = const Array1D& (Basis::*)(const Array1D& x);
   using DerivativesFptr = const Array2D& (Basis::*)(const Array1D& x);

   BasisArray basis_array_;
   basis_type btype_;
   FunctionsFptr functions_ptr_{nullptr};
   DerivativesFptr derivatives_fptr_{nullptr};

   void eval_jacobian(const QuadDomain& q, Matrix2D& JT);
   void eval_jacobian_serial(const QuadDomain& q, Matrix2D& JT);
#ifdef _OPENMP
   void eval_jacobian_omp(const QuadDomain& q, Matrix2D& JT);
#endif
};


}  // namespace gquad

#endif

