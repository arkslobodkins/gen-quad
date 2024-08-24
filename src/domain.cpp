// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/domain.hpp"

namespace gquad {


static IdealPolytope::Constraints CubeOrLineConstraints(gq_int dim);


////////////////////////////////////////////////////////////////////////////////////////////////////
IdealPolytope::IdealPolytope(Constraints c) : Polytope{}, A_(std::move(c.first)), b_(std::move(c.second)) {
}


// Applies to all polytopes determined by linear inequality constraints
bool IdealPolytope::in_domain(const Array1D& x) const {
   for(gq_int i = 0; i < A_.rows(); ++i) {
      if(dot_prod(A_.row(i), x) > b_[i]) {
         return false;
      }
   }
   return true;
}


double IdealPolytope::dist_from_boundary(const Array1D& x) const {
   auto dist = A_ * x.matrix() - b_;
   auto max = dist.maxCoeff();
   if(max > 0.) {
      return std::numeric_limits<double>::infinity();
   }
   return max;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
IdealPolytope::Constraints Interval::ConstructConstraints() {
   return CubeOrLineConstraints(1);
}


Interval::Interval() : IdealPolytope(ConstructConstraints()) {
}


////////////////////////////////////////////////////////////////////////////////////////////////////
IdealPolytope::Constraints Cube::ConstructConstraints(gq_int dim) {
   return CubeOrLineConstraints(dim);
}


Cube::Cube(gq_int dim) : IdealPolytope{ConstructConstraints(dim)}, dim_{dim} {
   GEN_QUAD_ASSERT_DEBUG(dim_ > 1);
}


Cube& Cube::operator=(const Cube& C) {
   GEN_QUAD_ASSERT_DEBUG(dim_ == C.dim_);
   return *this;
}


void Cube::reinit(gq_int dim) {
   dim_ = dim;
   Cube c(dim);
   A_ = c.A_;
   b_ = c.b_;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
IdealPolytope::Constraints Simplex::ConstructConstraints(gq_int dim) {
   Matrix2D A(dim + 1, dim);
   A(0, 0) = 1.;
   for(gq_int i = 1; i < A.rows() - 1; ++i) {
      A(i, i - 1) = -1.;
      A(i, i) = +1.;
   }
   A(A.size() - 1) = -1.;

   Vector1D b(dim + 1);
   b[0] = +1.;
   return {A, b};
}


Simplex::Simplex(gq_int dim) : IdealPolytope{ConstructConstraints(dim)}, dim_{dim} {
   GEN_QUAD_ASSERT_DEBUG(dim_ > 1);
}


Simplex& Simplex::operator=(const Simplex& S) {
   GEN_QUAD_ASSERT_DEBUG(dim_ == S.dim_);
   return *this;
}


bool Simplex::in_domain(const Array1D& x) const {
   if(x[0] > 1. || x(last) < 0.) {
      return false;
   }

   for(gq_int i = 1; i < x.size(); ++i) {
      if(x[i] > x[i - 1]) {
         return false;
      }
   }

   return true;
}


void Simplex::reinit(gq_int dim) {
   dim_ = dim;
   Simplex s(dim);
   A_ = s.A_;
   b_ = s.b_;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
IdealPolytope::Constraints CubeSimplex::ConstructConstraints(gq_int dim1, gq_int dim2) {
   Constraints cube_or_line = CubeOrLineConstraints(dim1);
   Simplex S{dim2};
   const auto& C_A = cube_or_line.first;
   const auto& C_b = cube_or_line.second;
   const auto& S_A = S.get_constraint_matrix();
   const auto& S_b = S.get_constraint_vector();

   gq_int c_rows = C_A.rows();
   gq_int c_cols = C_A.cols();
   gq_int s_rows = S_A.rows();
   gq_int s_cols = S_A.cols();
   gq_int rows = c_rows + s_rows;
   gq_int cols = c_cols + s_cols;

   Matrix2D A(rows, cols);
   A(seqN(0, c_rows), seqN(0, c_cols)) = C_A;
   A(seqN(c_rows, s_rows), seqN(c_cols, s_cols)) = S_A;

   Vector1D b(rows);
   b(seqN(0, c_rows)) = C_b;
   b(seqN(c_rows, s_rows)) = S_b;

   return {A, b};
}


CubeSimplex::CubeSimplex(std::array<gq_int, 2> dims)
    : IdealPolytope{ConstructConstraints(dims[0], dims[1])},
      dims_(dims) {
   GEN_QUAD_ASSERT_DEBUG(dims_[0] > 0);
   GEN_QUAD_ASSERT_DEBUG(dims_[1] > 1);
}


CubeSimplex& CubeSimplex::operator=(const CubeSimplex& cs) {
   GEN_QUAD_ASSERT_DEBUG(dims_ == cs.dims_);
   return *this;
}


bool CubeSimplex::in_domain(const Array1D& x) const {
   gq_int d1 = dim1(), d2 = dim2();

   for(gq_int i = 0; i < d1; ++i) {
      if(dot_prod(A_.row(i), x) > b_[i]) {
         return false;
      }
   }

   if(x[d1] > 1. || x(last) < 0.) {
      return false;
   }

   for(gq_int i = 1; i < d2; ++i) {
      if(x[d1 + i] > x[d1 + i - 1]) {
         return false;
      }
   }

   return true;
}


void CubeSimplex::reinit(gq_int dim1, gq_int dim2) {
   dims_ = {dim1, dim2};
   CubeSimplex cs(dims_);
   A_ = cs.A_;
   b_ = cs.b_;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
IdealPolytope::Constraints SimplexSimplex::ConstructConstraints(gq_int dim1, gq_int dim2) {
   Simplex S1{dim1};
   Simplex S2{dim2};
   const auto& S1_A = S1.get_constraint_matrix();
   const auto& S1_b = S1.get_constraint_vector();
   const auto& S2_A = S2.get_constraint_matrix();
   const auto& S2_b = S2.get_constraint_vector();

   gq_int s1_rows = S1_A.rows();
   gq_int s1_cols = S1_A.cols();
   gq_int s2_rows = S2_A.rows();
   gq_int s2_cols = S2_A.cols();
   gq_int rows = s1_rows + s2_rows;
   gq_int cols = s1_cols + s2_cols;

   Matrix2D A(rows, cols);
   A(seqN(0, s1_rows), seqN(0, s1_cols)) = S1_A;
   A(seqN(s1_rows, s2_rows), seqN(s1_cols, s2_cols)) = S2_A;

   Vector1D b(rows);
   b(seqN(0, s1_rows)) = S1_b;
   b(seqN(s1_rows, s2_rows)) = S2_b;

   return {A, b};
}


SimplexSimplex::SimplexSimplex(std::array<gq_int, 2> dims)
    : IdealPolytope{ConstructConstraints(dims[0], dims[1])},
      dims_(dims) {
   GEN_QUAD_ASSERT_DEBUG(dims_[0] > 1);
   GEN_QUAD_ASSERT_DEBUG(dims_[1] > 1);
}


SimplexSimplex& SimplexSimplex::operator=(const SimplexSimplex& ss) {
   GEN_QUAD_ASSERT_DEBUG(dims_ == ss.dims_);
   return *this;
}


bool SimplexSimplex::in_domain(const Array1D& x) const {
   gq_int d1 = dim1(), d2 = dim2();

   if(x[0] > 1. || x[d1 - 1] < 0.) {
      return false;
   }

   for(gq_int i = 1; i < d1; ++i) {
      if(x[i] > x[i - 1]) {
         return false;
      }
   }

   if(x[d1] > 1. || x(last) < 0.) {
      return false;
   }

   for(gq_int i = 1; i < d2; ++i) {
      if(x[d1 + i] > x[d1 + i - 1]) {
         return false;
      }
   }

   return true;
}


void SimplexSimplex::reinit(gq_int dim1, gq_int dim2) {
   dims_ = {dim1, dim2};
   SimplexSimplex ss(dims_);
   A_ = ss.A_;
   b_ = ss.b_;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
IdealPolytope::Constraints Pyramid3D::ConstructConstraints() {
   Matrix2D A(6, 3);

   // [-1 +0 +0] x >= 0
   // [+1 +0 +0] x <= 1
   // [+0 -1 +0] y >= 0
   // [-1 +1 +0] y <= x
   // [+0 +0 -1] z >= 0
   // [-1 +0 +1] z <= x

   A(0, 0) = -1.;
   A(1, 0) = +1.;
   A(2, 1) = -1.;
   A(3, 0) = -1.;
   A(3, 1) = +1.;
   A(4, 2) = -1.;
   A(5, 0) = -1.;
   A(5, 2) = +1.;

   Vector1D b(6);
   b[1] = +1.;
   return {A, b};
}


Pyramid3D::Pyramid3D() : IdealPolytope{ConstructConstraints()} {
}


bool Pyramid3D::in_domain(const Array1D& x) const {
   if(x[1] < 0. || x[2] < 0. || x[0] > 1. || x[1] > x[0] || x[2] > x[0]) {
      return false;
   }
   return true;
}


static IdealPolytope::Constraints CubeOrLineConstraints(gq_int dim) {
   Matrix2D A(2 * dim, dim);
   Vector1D b(2 * dim);

   for(gq_int j = 0; j < A.cols(); ++j) {
      A(2 * j, j) = -1.;
      A(2 * j + 1, j) = +1.;
   }
   for(gq_int i = 0; i < b.size(); i += 2) {
      b[i + 1] = +1.;
   }
   return {A, b};
}


}  // namespace gquad

#endif

