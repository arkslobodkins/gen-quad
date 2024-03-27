// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#define EIGEN_NO_STATIC_ASSERT
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#ifdef GEN_QUAD_USE_MKL
#define EIGEN_USE_MKL_ALL
#endif

#include "../eigen-3.4.0/Eigen/Dense"
#include "exceptions.hpp"
#include "headers.hpp"

namespace gquad {

constexpr double pi = 3.141'592'653'589'793'2;

using gq_int = long int;

template <gq_int N>
using StaticArrayInt1D = Eigen::Array<int, N, 1, Eigen::ColMajor>;
template <gq_int N>
using StaticArray1D = Eigen::Array<double, N, 1, Eigen::ColMajor>;
template <gq_int M, gq_int N>
using StaticArray2D = Eigen::Array<double, M, N, Eigen::RowMajor>;

using ArrayInt1D = Eigen::Array<int, Eigen::Dynamic, 1, Eigen::ColMajor>;
using Array1D = Eigen::Array<double, Eigen::Dynamic, 1, Eigen::ColMajor>;
using Array1DLong = Eigen::Array<long double, Eigen::Dynamic, 1, Eigen::ColMajor>;
using Vector1D = Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::ColMajor>;

using ArrayInt2D = Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Array2D = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Matrix2D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

using Eigen::all;
using Eigen::last;
using Eigen::seq;
using Eigen::seqN;

using Eigen::Infinity;

class Array3D {
public:
   Array3D(gq_int i, gq_int j, gq_int k) : dims{i, j, k}, arr(i * j * k) {
      GEN_QUAD_ASSERT_DEBUG(i >= 0);
      GEN_QUAD_ASSERT_DEBUG(j >= 0);
      GEN_QUAD_ASSERT_DEBUG(k >= 0);
   }

   auto operator()(gq_int i, gq_int j) {
      GEN_QUAD_ASSERT_DEBUG(i > -1 && i < dims[0]);
      GEN_QUAD_ASSERT_DEBUG(j > -1 && j < dims[1]);
      return arr(seqN(i * dims[1] * dims[2] + j * dims[2], dims[2]));
   }

   double operator()(gq_int i, gq_int j, gq_int k) {
      GEN_QUAD_ASSERT_DEBUG(i > -1 && i < dims[0]);
      GEN_QUAD_ASSERT_DEBUG(j > -1 && j < dims[1]);
      GEN_QUAD_ASSERT_DEBUG(k > -1 && k < dims[2]);
      return arr[i * dims[1] * dims[2] + j * dims[2] + k];
   }

private:
   std::array<gq_int, 3> dims;
   Array1D arr;
};

constexpr double quad_tol{1.e-14};
constexpr double quad_huge{10.};

// avoid warnings on conversions from signed to unsigned and vice versa
template <typename T>
class StdVector : private std::vector<T> {
public:
   explicit StdVector() : std::vector<T>{} {
   }

   StdVector(std::initializer_list<T> list) : std::vector<T>(list) {
   }

   explicit StdVector(gq_int i) : std::vector<T>(static_cast<vsize_t>(i)) {
      GEN_QUAD_ASSERT_DEBUG(i > -1);
   }

   explicit StdVector(gq_int i, const T& val) : std::vector<T>(static_cast<vsize_t>(i), val) {
      GEN_QUAD_ASSERT_DEBUG(i > 0);
   }

   T& operator[](gq_int i) {
      GEN_QUAD_ASSERT_DEBUG(i > -1);
      return std::vector<T>::operator[](static_cast<vsize_t>(i));
   }

   const T& operator[](gq_int i) const {
      GEN_QUAD_ASSERT_DEBUG(i > -1);
      return std::vector<T>::operator[](static_cast<vsize_t>(i));
   }

   T& at(gq_int i) {
      GEN_QUAD_ASSERT_DEBUG(i > -1);
      return std::vector<T>::at(static_cast<vsize_t>(i));
   }

   const T& at(gq_int i) const {
      GEN_QUAD_ASSERT_DEBUG(i > -1);
      return std::vector<T>::at(static_cast<vsize_t>(i));
   }

   gq_int size() const {
      return static_cast<gq_int>(std::vector<T>::size());
   }

   void resize(gq_int i) {
      GEN_QUAD_ASSERT_DEBUG(i > -1);
      std::vector<T>::resize(static_cast<vsize_t>(i));
   }

   void reserve(gq_int i) {
      GEN_QUAD_ASSERT_DEBUG(i > -1);
      std::vector<T>::reserve(static_cast<vsize_t>(i));
   }

   using std::vector<T>::push_back;
   using std::vector<T>::insert;
   using std::vector<T>::begin;
   using std::vector<T>::end;
   using std::vector<T>::front;
   using std::vector<T>::back;
   using std::vector<T>::empty;

private:
   using vsize_t = typename std::vector<T>::size_type;
};

StdVector<gq_int> sequence(gq_int first, gq_int last, gq_int stride = 1);

template <typename T>
gq_int max_index(const T& z) {
   GEN_QUAD_ASSERT_DEBUG(z.size() > 0);
   gq_int max_ind = 0;
   for(gq_int i = 1; i < z.size(); ++i) {
      if(z(i) > z(max_ind)) {
         max_ind = i;
      }
   }

   return max_ind;
}

template <typename EigenArrayType>
double norm_inf(const EigenArrayType& x) {
   GEN_QUAD_ASSERT_DEBUG(x.size() > 0);
   return x.matrix().template lpNorm<Infinity>();
}

template <typename EigenArrayType>
double norm2(const EigenArrayType& x) {
   GEN_QUAD_ASSERT_DEBUG(x.size() > 0);
   return x.matrix().template lpNorm<2>();
}

template <typename EigenArrayType1, typename EigenArrayType2>
double dot_prod(const EigenArrayType1& x, const EigenArrayType2& y) {
   GEN_QUAD_ASSERT_DEBUG(x.size() > 0);
   GEN_QUAD_ASSERT_DEBUG(x.size() == y.size());
   if(x.rows() == y.rows()) {
      return (x.array() * y.array()).sum();
   } else {
      return (x.array() * y.array().transpose()).sum();
   }
}

template <typename EigenType>
double stable_sum(const EigenType& A) {
   GEN_QUAD_ASSERT_DEBUG(A.size() > 0);

   double s{};
   double c{};
   for(gq_int i = 0; i < A.size(); ++i) {
      volatile double t = s + A[i];
      if(std::abs(s) >= std::abs(A[i])) {
         volatile double z = s - t;
         c += z + A[i];
      } else {
         volatile double z = A[i] - t;
         c += z + s;
      }
      s = t;
   }
   return s + c;
}

template <typename T>
auto two_sum(T x, T y) {
   volatile T r = x + y;
   volatile T z = r - x;
   T s = ((x - (r - z)) + (y - z));
   return std::pair<T, T>{r, s};
}

template <typename T>
auto two_prod(T v1, T v2) {
   auto r = v1 * v2;
   auto s = std::fma<T>(v1, v2, -r);
   return std::pair<T, T>{r, s};
}

template <typename EigenType1, typename EigenType2>
auto stable_dot_prod(const EigenType1& A, const EigenType2& B) {
   using T = typename EigenType1::value_type;
   T p, s, q;

   std::tie(p, s) = two_prod(A[0], B[0]);
   for(gq_int i = 1; i < A.size(); ++i) {
      T h, r;
      std::tie(h, r) = two_prod(A[i], B[i]);
      std::tie(p, q) = two_sum(p, h);
      s = s + (q + r);
   }

   return p + s;
}

struct SearchWidth {
   explicit SearchWidth(gq_int width) : width{width} {
      GEN_QUAD_ASSERT_ALWAYS(width > 0 && "search width must be greater or equal to 1");
   }

   gq_int width;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// structures used to output history to a file

struct ElimData {
   gq_int nodes_tot{-1};
   gq_int success_node{-1};
   gq_int success_its{-1};
   gq_int num_solutions{-1};
   gq_int num_fails{-1};
};

struct History {
   StdVector<ElimData> hlist{};
   gq_int dim{-1};
   gq_int deg{-1};
   gq_int nodes_initial{-1};
   gq_int nodes_final{-1};
   gq_int nodes_optimal{-1};
   gq_int num_funcs{-1};
   gq_int total_elims{-1};
   double res{-1.};
   double efficiency{-1.};
   std::string domain_name{};
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// determines whether the problem should be parallelized or not
#ifdef _OPENMP
bool omp_condition(gq_int deg, gq_int dim);
#endif

// Establishes reasonable upper bounds for the maximum
// degree for a given dimension. Anything above these
// upper bounds is assumed to be either too time consuming
// or outside the scope of this work.
// In addition, it prevents the user from accidentally entering extremely high
// degrees and dimensions.
// Reasonable upper bounds are subject to change once symmetry is employed.
bool is_reasonable_degree_and_dimension(gq_int deg, gq_int dim);

}  // namespace gquad

#endif

