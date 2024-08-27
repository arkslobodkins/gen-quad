// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/basis.hpp"

#include "../include/domain.hpp"
#include "../include/math_func.hpp"
#include "../include/quad_tensor.hpp"
#include "../include/quadrature.hpp"

namespace gquad {


template <typename EigenType>
static inline void Monomial(double x, EigenType&& p);


template <typename EigenType>
static inline void Legendre(double x, EigenType&& p);


template <typename EigenType>
static void Jacobi(gq_int alpha, double x, EigenType&& p);


static void StandardMonomialFunctions(const BasisTable& index_table, const Array1D& x, Array1D& functions);
static void StandardMonomialDerivatives(const BasisTable& index_table, const Array1D& x,
                                        Array2D& derivatives);


////////////////////////////////////////////////////////////////////////////////////////////////////////
Basis::Basis(gq_int deg, gq_int dim, gq_int num_funcs)
    : deg_{deg},
      dim_{dim},
      num_funcs_{num_funcs},
      functions_(num_funcs),
      derivatives_(dim, num_funcs),
      integrals_(num_funcs) {
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
PolytopeBasis::PolytopeBasis(gq_int deg, gq_int dim, gq_int num_funcs) : Basis(deg, dim, num_funcs) {
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cube Basis
CubeBasis::CubeBasis(gq_int deg, gq_int dim)
    : PolytopeBasis(deg, dim, CubeBasisSize(deg, dim)),
      index_table_(deg, dim),
      fdiff_{} {
   GEN_QUAD_ASSERT_DEBUG(deg_ >= 1);
   GEN_QUAD_ASSERT_DEBUG(dim_ >= 2);

   fdiff_[0].resize(num_funcs_);
   fdiff_[1].resize(num_funcs_);
   fdiff_[2].resize(num_funcs_);
   fdiff_[3].resize(num_funcs_);
   StandardBasisIndices(index_table_);
}


const Array1D& CubeBasis::orthog_basis(const Array1D& x) {
   this->orthog_basis_internal(x, functions_);
   return functions_;
}


void CubeBasis::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   Array2D legendre(dim_, deg_ + 1);
   for(gq_int d = 0; d < dim_; ++d) {
      Legendre(x[d], legendre.row(d));
   }

   for(gq_int k = 0; k < num_funcs_; ++k) {
      buff[k] = legendre(0, index_table_(0, k));
   }
   for(gq_int d = 1; d < dim_; ++d) {
      for(gq_int k = 0; k < num_funcs_; ++k) {
         buff[k] *= legendre(d, index_table_(d, k));
      }
   }
}


const Array2D& CubeBasis::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives_;
}


const Array1D& CubeBasis::orthog_integrals() {
   integrals_ = 0.;
   integrals_[0] = 1.;
   return integrals_;
}


const Array1D& CubeBasis::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table_, x, functions_);
   return functions_;
}


const Array2D& CubeBasis::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table_, x, derivatives_);
   return derivatives_;
}


const Array1D& CubeBasis::monomial_integrals() {
   for(gq_int k = 0; k < num_funcs_; ++k) {
      // use integer instead of double and perform multiplication instead of division
      // to avoid error associated with intermediary divisions
      gq_int prod = 1;
      for(gq_int d = 0; d < dim_; ++d) {
         prod *= (index_table_(d, k) + 1);
      }
      integrals_[k] = 1. / prod;
   }
   return integrals_;
}


std::unique_ptr<Basis> CubeBasis::clone() const {
   return std::make_unique<CubeBasis>(*this);
}


gq_int CubeBasisSize(gq_int deg, gq_int dim) {
   GEN_QUAD_ASSERT_DEBUG(deg >= 1);
   GEN_QUAD_ASSERT_DEBUG(dim >= 2);
   return StandardBasisSize(deg, dim);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Simplex Basis
SimplexBasis::SimplexBasis(gq_int deg, gq_int dim)
    : PolytopeBasis(deg, dim, SimplexBasisSize(deg, dim)),
      power_table_(deg, dim),
      index_table_(deg, dim),
      fdiff_{} {
   GEN_QUAD_ASSERT_DEBUG(deg_ >= 1);
   GEN_QUAD_ASSERT_DEBUG(dim_ >= 2);

   fdiff_[0].resize(num_funcs_);
   fdiff_[1].resize(num_funcs_);
   fdiff_[2].resize(num_funcs_);
   fdiff_[3].resize(num_funcs_);

   StandardBasisIndices(index_table_);
   for(gq_int d = 1; d < dim_; ++d) {
      for(gq_int i = 0; i < d; ++i) {
         for(gq_int k = 0; k < num_funcs_; ++k) {
            power_table_(d, k) += index_table_(i, k);
         }
      }
   }
}


const Array1D& SimplexBasis::orthog_basis(const Array1D& x) {
   this->orthog_basis_internal(x, functions_);
   return functions_;
}


const Array2D& SimplexBasis::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives_;
}


const Array1D& SimplexBasis::orthog_integrals() {
   integrals_ = 0.;
   integrals_[0] = 1. / math::factorial(dim_);
   return integrals_;
}


const Array1D& SimplexBasis::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table_, x, functions_);
   return functions_;
}


const Array2D& SimplexBasis::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table_, x, derivatives_);
   return derivatives_;
}


const Array1D& SimplexBasis::monomial_integrals() {
   for(gq_int k = 0; k < num_funcs_; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < dim_; ++d) {
         gq_int s = index_table_.array()(seq(d, dim_ - 1), k).sum();
         prod *= (s + dim_ - d);
      }
      integrals_[k] = 1. / prod;
   }
   return integrals_;
}


void SimplexBasis::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   Array1D legendre(deg_ + 1);
   Legendre((x[dim_ - 2] - x[dim_ - 1]) / x[dim_ - 2], legendre);  // mapped (2y-x)/x to (x-y)/x

   ArrayLong1D xCoord(dim_ - 1);
   ArrayLong1D jCoord(dim_ - 1);
   for(gq_int d = 0; d < dim_ - 2; ++d) {
      xCoord[d] = x[dim_ - d - 2] / x[dim_ - d - 3];
      jCoord[d] = 1.L - 2.L * x[dim_ - d - 2] / x[dim_ - d - 3];
   }
   xCoord[dim_ - 2] = x[0];
   jCoord[dim_ - 2] = 1.L - 2.L * x[0];

   Array3D jacobi(dim_, deg_ + 1, deg_ + 1);
   for(gq_int d = 1; d < dim_; ++d) {
      for(gq_int j = 0; j < deg_ + 1; ++j) {
         Jacobi(2 * j + d, jCoord[d - 1], jacobi(d, j));
      }
   }

   Array2D map_table(dim_, deg_ + 1);
   for(gq_int d = 1; d < map_table.rows(); ++d) {
      for(gq_int j = 0; j < map_table.cols(); ++j) {
         map_table(d, j) = std::pow<double>(xCoord[d - 1], j);
      }
   }

   for(gq_int k = 0; k < num_funcs_; ++k) {
      buff[k] = legendre[index_table_(0, k)];
   }

   for(gq_int d = 1; d < dim_; ++d) {
      for(gq_int k = 0; k < num_funcs_; ++k) {
         buff[k] *= jacobi(d, power_table_(d, k), index_table_(d, k)) * map_table(d, power_table_(d, k));
      }
   }
}


std::unique_ptr<Basis> SimplexBasis::clone() const {
   return std::make_unique<SimplexBasis>(*this);
}


gq_int SimplexBasisSize(gq_int deg, gq_int dim) {
   GEN_QUAD_ASSERT_DEBUG(deg >= 1);
   GEN_QUAD_ASSERT_DEBUG(dim >= 2);
   return StandardBasisSize(deg, dim);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// CubeSimplex Basis
CubeSimplexBasis::CubeSimplexBasis(gq_int deg, gq_int dim1, gq_int dim2)
    : PolytopeBasis(deg, dim1 + dim2, CubeSimplexBasisSize(deg, dim1 + dim2)),
      dims_{dim1, dim2},
      power_table_second_(deg, dim2, CubeSimplexBasisSize(deg, dim1 + dim2)),
      index_table_(deg, dim1 + dim2),
      fdiff_{} {
   GEN_QUAD_ASSERT_DEBUG(deg_ >= 1);
   GEN_QUAD_ASSERT_DEBUG(dims_[0] >= 1);
   GEN_QUAD_ASSERT_DEBUG(dims_[1] >= 2);

   fdiff_[0].resize(num_funcs_);
   fdiff_[1].resize(num_funcs_);
   fdiff_[2].resize(num_funcs_);
   fdiff_[3].resize(num_funcs_);

   StandardBasisIndices(index_table_);
   for(gq_int d = 1; d < dims_[1]; ++d) {
      for(gq_int i = 0; i < d; ++i) {
         for(gq_int k = 0; k < num_funcs_; ++k) {
            power_table_second_(d, k) += index_table_(dims_[0] + i, k);
         }
      }
   }
}


const Array1D& CubeSimplexBasis::orthog_basis(const Array1D& x) {
   this->orthog_basis_internal(x, functions_);
   return functions_;
}


const Array2D& CubeSimplexBasis::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives_;
}


const Array1D& CubeSimplexBasis::orthog_integrals() {
   integrals_ = 0.;
   integrals_[0] = 1. / math::factorial(dims_[1]);
   return integrals_;
}


const Array1D& CubeSimplexBasis::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table_, x, functions_);
   return functions_;
}


const Array2D& CubeSimplexBasis::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table_, x, derivatives_);
   return derivatives_;
}


const Array1D& CubeSimplexBasis::monomial_integrals() {
   // integrals due to cube
   for(gq_int k = 0; k < num_funcs_; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < dims_[0]; ++d) {
         prod *= (index_table_(d, k) + 1);
      }
      integrals_[k] = 1. / prod;
   }

   // integrals due to simplex
   for(gq_int k = 0; k < num_funcs_; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < dims_[1]; ++d) {
         gq_int s = index_table_.array()(seq(dims_[0] + d, dim_ - 1), k).sum();
         prod *= (s + dims_[1] - d);
      }
      integrals_[k] /= prod;
   }

   return integrals_;
}


void CubeSimplexBasis::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   this->orthog_basis_polytopic_two_internal(x, buff);

   Array1D legendre(deg_ + 1);
   for(gq_int d = 0; d < dims_[0]; ++d) {
      Legendre(x[d], legendre);
      for(gq_int k = 0; k < num_funcs_; ++k) {
         buff[k] *= legendre[index_table_(d, k)];
      }
   }
}


void CubeSimplexBasis::orthog_basis_polytopic_two_internal(const Array1D& x, Array1D& buff) {
   gq_int dim1{dims_[0]};
   gq_int dim2{dims_[1]};
   gq_int dim_two{dim1 + dim2};

   Array1D legendre(deg_ + 1);
   // mapped (2y-x)/x to (x-y)/x
   Legendre((x[dim_two - 2] - x[dim_two - 1]) / x[dim_two - 2], legendre);

   ArrayLong1D xCoord(dim2 - 1);
   for(gq_int d = 0; d < dim2 - 2; ++d) {
      xCoord[d] = x[dim_two - d - 2] / x[dim_two - d - 3];
   }
   xCoord[dim2 - 2] = x[dim1];

   ArrayLong1D jCoord(dim2 - 1);
   for(gq_int d = 1; d < dim2 - 1; ++d) {
      jCoord[d - 1] = 1.L - 2.L * x[dim_two - d - 1] / x[dim_two - d - 2];
   }
   jCoord[dim2 - 2] = 1.L - 2.L * x[dim1];

   Array3D jacobi(dim2, deg_ + 1, deg_ + 1);
   for(gq_int d = 1; d < dim2; ++d) {
      for(gq_int j = 0; j < deg_ + 1; ++j) {
         Jacobi(2 * j + d, jCoord[d - 1], jacobi(d, j));
      }
   }

   Array2D map_table(dim2, deg_ + 1);
   for(gq_int d = 1; d < map_table.rows(); ++d) {
      for(gq_int j = 0; j < map_table.cols(); ++j) {
         map_table(d, j) = std::pow<double>(xCoord[d - 1], j);
      }
   }

   for(gq_int k = 0; k < num_funcs_; ++k) {
      buff[k] = legendre[index_table_(dim1, k)];
   }

   for(gq_int d = 1; d < dim2; ++d) {
      for(gq_int k = 0; k < num_funcs_; ++k) {
         buff[k] *= jacobi(d, power_table_second_(d, k), index_table_(dim1 + d, k))
                  * map_table(d, power_table_second_(d, k));
      }
   }
}


std::unique_ptr<Basis> CubeSimplexBasis::clone() const {
   return std::make_unique<CubeSimplexBasis>(*this);
}


gq_int CubeSimplexBasisSize(gq_int deg, gq_int dim) {
   GEN_QUAD_ASSERT_DEBUG(deg >= 1);
   GEN_QUAD_ASSERT_DEBUG(dim >= 3);
   return StandardBasisSize(deg, dim);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// SimplexSimplex Basis
SimplexSimplexBasis::SimplexSimplexBasis(gq_int deg, gq_int dim1, gq_int dim2)
    : PolytopeBasis(deg, dim1 + dim2, SimplexSimplexBasisSize(deg, dim1 + dim2)),
      dims_{dim1, dim2},
      power_table_first_(deg, dim1, SimplexSimplexBasisSize(deg, dim1 + dim2)),
      power_table_second_(deg, dim2, SimplexSimplexBasisSize(deg, dim1 + dim2)),
      index_table_(deg, dim1 + dim2),
      poly_buff1_{},
      poly_buff2_{},
      fdiff_{} {
   GEN_QUAD_ASSERT_DEBUG(deg_ >= 1);
   GEN_QUAD_ASSERT_DEBUG(dims_[0] >= 2);
   GEN_QUAD_ASSERT_DEBUG(dims_[1] >= 2);

   poly_buff1_.resize(num_funcs_);
   poly_buff2_.resize(num_funcs_);
   fdiff_[0].resize(num_funcs_);
   fdiff_[1].resize(num_funcs_);
   fdiff_[2].resize(num_funcs_);
   fdiff_[3].resize(num_funcs_);

   StandardBasisIndices(index_table_);

   for(gq_int d = 1; d < dims_[0]; ++d) {
      for(gq_int i = 0; i < d; ++i) {
         for(gq_int k = 0; k < num_funcs_; ++k) {
            power_table_first_(d, k) += index_table_(i, k);
         }
      }
   }

   for(gq_int d = 1; d < dims_[1]; ++d) {
      for(gq_int i = 0; i < d; ++i) {
         for(gq_int k = 0; k < num_funcs_; ++k) {
            power_table_second_(d, k) += index_table_(dims_[0] + i, k);
         }
      }
   }
}


const Array1D& SimplexSimplexBasis::orthog_basis(const Array1D& x) {
   this->orthog_basis_internal(x, functions_);
   return functions_;
}


const Array2D& SimplexSimplexBasis::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives_;
}


const Array1D& SimplexSimplexBasis::orthog_integrals() {
   integrals_ = 0.;
   integrals_[0] = 1. / (math::factorial(dims_[0]) * math::factorial(dims_[1]));
   return integrals_;
}


const Array1D& SimplexSimplexBasis::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table_, x, functions_);
   return functions_;
}


const Array2D& SimplexSimplexBasis::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table_, x, derivatives_);
   return derivatives_;
}


const Array1D& SimplexSimplexBasis::monomial_integrals() {
   // integrals due to simplex1
   for(gq_int k = 0; k < num_funcs_; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < dims_[0]; ++d) {
         gq_int s = index_table_.array()(seq(d, dims_[0] - 1), k).sum();
         prod *= (s + dims_[0] - d);
      }
      integrals_[k] = 1. / prod;
   }

   // integrals due to simplex2
   for(gq_int k = 0; k < num_funcs_; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < dims_[1]; ++d) {
         gq_int s = index_table_.array()(seq(dims_[0] + d, dim_ - 1), k).sum();
         prod *= (s + dims_[1] - d);
      }
      integrals_[k] /= prod;
   }

   return integrals_;
}


void SimplexSimplexBasis::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   this->orthog_basis_polytopic_one_internal(x, poly_buff1_);
   this->orthog_basis_polytopic_two_internal(x, poly_buff2_);

   buff = 1.;
   buff *= poly_buff1_;
   buff *= poly_buff2_;
}


void SimplexSimplexBasis::orthog_basis_polytopic_one_internal(const Array1D& x, Array1D& buff) {
   gq_int dim1{dims_[0]};

   Array1D legendre(deg_ + 1);
   // mapped (2y-x)/x to (x-y)/x
   Legendre((x[dim1 - 2] - x[dim1 - 1]) / x[dim1 - 2], legendre);

   ArrayLong1D xCoord(dim1 - 1);
   for(gq_int d = 0; d < dim1 - 2; ++d) {
      xCoord[d] = x[dim1 - d - 2] / x[dim1 - d - 3];
   }
   xCoord[dim1 - 2] = x[0];

   ArrayLong1D jCoord(dim1 - 1);
   for(gq_int d = 1; d < dim1 - 1; ++d) {
      jCoord[d - 1] = 1.L - 2.L * x[dim1 - d - 1] / x[dim1 - d - 2];
   }
   jCoord[dim1 - 2] = 1.L - 2.L * x[0];

   Array3D jacobi(dim1, deg_ + 1, deg_ + 1);
   for(gq_int d = 1; d < dim1; ++d) {
      for(gq_int j = 0; j < deg_ + 1; ++j) {
         Jacobi(2 * j + d, jCoord[d - 1], jacobi(d, j));
      }
   }

   Array2D map_table(dim1, deg_ + 1);
   for(gq_int d = 1; d < map_table.rows(); ++d) {
      for(gq_int j = 0; j < map_table.cols(); ++j) {
         map_table(d, j) = std::pow<double>(xCoord[d - 1], j);
      }
   }

   for(gq_int k = 0; k < num_funcs_; ++k) {
      buff[k] = legendre[index_table_(0, k)];
   }

   for(gq_int d = 1; d < dim1; ++d) {
      for(gq_int k = 0; k < num_funcs_; ++k) {
         buff[k] *= jacobi(d, power_table_first_(d, k), index_table_(d, k))
                  * map_table(d, power_table_first_(d, k));
      }
   }
}


void SimplexSimplexBasis::orthog_basis_polytopic_two_internal(const Array1D& x, Array1D& buff) {
   gq_int dim1{dims_[0]};
   gq_int dim2{dims_[1]};
   gq_int dim_two{dim1 + dim2};

   Array1D legendre(deg_ + 1);
   Legendre((x[dim_two - 2] - x[dim_two - 1]) / x[dim_two - 2], legendre);  // mapped (2y-x)/x to (x-y)/x

   ArrayLong1D xCoord(dim2 - 1);
   for(gq_int d = 0; d < dim2 - 2; ++d) {
      xCoord[d] = x[dim_two - d - 2] / x[dim_two - d - 3];
   }
   xCoord[dim2 - 2] = x[dim1];

   ArrayLong1D jCoord(dim2 - 1);
   for(gq_int d = 1; d < dim2 - 1; ++d) {
      jCoord[d - 1] = 1.L - 2.L * x[dim_two - d - 1] / x[dim_two - d - 2];
   }
   jCoord[dim2 - 2] = 1.L - 2.L * x[dim1];

   Array3D jacobi(dim2, deg_ + 1, deg_ + 1);
   for(gq_int d = 1; d < dim2; ++d) {
      for(gq_int j = 0; j < deg_ + 1; ++j) {
         Jacobi(2 * j + d, jCoord[d - 1], jacobi(d, j));
      }
   }

   Array2D map_table(dim2, deg_ + 1);
   for(gq_int d = 1; d < map_table.rows(); ++d) {
      for(gq_int j = 0; j < map_table.cols(); ++j) {
         map_table(d, j) = std::pow<double>(xCoord[d - 1], j);
      }
   }

   for(gq_int k = 0; k < num_funcs_; ++k) {
      buff[k] = legendre[index_table_(dim1, k)];
   }

   for(gq_int d = 1; d < dim2; ++d) {
      for(gq_int k = 0; k < num_funcs_; ++k) {
         buff[k] *= jacobi(d, power_table_second_(d, k), index_table_(dim1 + d, k))
                  * map_table(d, power_table_second_(d, k));
      }
   }
}


std::unique_ptr<Basis> SimplexSimplexBasis::clone() const {
   return std::make_unique<SimplexSimplexBasis>(*this);
}


gq_int SimplexSimplexBasisSize(gq_int deg, gq_int dim) {
   GEN_QUAD_ASSERT_DEBUG(deg >= 1);
   GEN_QUAD_ASSERT_DEBUG(dim >= 4);
   return StandardBasisSize(deg, dim);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Pyramid3D Basis
PyramidBasis3D::PyramidBasis3D(gq_int deg)
    : PolytopeBasis(deg, 3, Pyramid3DBasisSize(deg)),
      index_table_(deg, 3),
      fdiff_{} {
   GEN_QUAD_ASSERT_DEBUG(deg_ >= 1);

   fdiff_[0].resize(num_funcs_);
   fdiff_[1].resize(num_funcs_);
   fdiff_[2].resize(num_funcs_);
   fdiff_[3].resize(num_funcs_);
   StandardBasisIndices(index_table_);
}


const Array1D& PyramidBasis3D::orthog_basis(const Array1D& x) {
   this->orthog_basis_internal(x, functions_);
   return functions_;
}


const Array2D& PyramidBasis3D::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives_;
}


const Array1D& PyramidBasis3D::orthog_integrals() {
   integrals_ = 0.;
   integrals_[0] = 1. / 3.;
   return integrals_;
}


const Array1D& PyramidBasis3D::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table_, x, functions_);
   return functions_;
}


const Array2D& PyramidBasis3D::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table_, x, derivatives_);
   return derivatives_;
}


const Array1D& PyramidBasis3D::monomial_integrals() {
   for(gq_int k = 0; k < num_funcs_; ++k) {
      integrals_[k] = 1.
                    / ((index_table_(1, k) + 1) * (index_table_(2, k) + 1)
                       * (index_table_(0, k) + index_table_(1, k) + index_table_(2, k) + 3));
   }
   return integrals_;
}


void PyramidBasis3D::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   Array2D legendre(2, deg_ + 1);
   Legendre(1.L - (x[1] / x[0]), legendre.row(0));  // mapped 2(y/x)-1 to 1-y/x
   Legendre(1.L - (x[2] / x[0]), legendre.row(1));  // mapped 2(z/x)-1 to 1-z/x

   Array2D jacobi(deg_ + 1, deg_ + 1);
   for(gq_int j = 0; j < deg_ + 1; ++j) {
      Jacobi(2 * j + 2, 1.L - 2.L * x[0], jacobi.row(j));
   }

   Array1D xpower(deg_ + 1);
   Monomial(x[0], xpower);

   gq_int count = 0;
   for(gq_int k = 0; k < deg_ + 1; ++k) {
      for(gq_int j = 0; j < deg_ + 1 - k; ++j) {
         for(gq_int i = 0; i < deg_ + 1 - j - k; ++i) {
            buff[count++] = xpower[j + k] * jacobi(j + k, i) * legendre(0, j) * legendre(1, k);
         }
      }
   }
}


std::unique_ptr<Basis> PyramidBasis3D::clone() const {
   return std::make_unique<PyramidBasis3D>(*this);
}


gq_int Pyramid3DBasisSize(gq_int deg) {
   GEN_QUAD_ASSERT_DEBUG(deg >= 1);
   return StandardBasisSize(deg, 3);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
// Omega2D Basis
OmegaBasis2D::OmegaBasis2D(Omega2D omega, gq_int deg)
    : Basis(deg, 2, Omega2DBasisSize(deg)),
      quad_simplex_(new QuadSimplex{GaussTensorSimplex(2 * deg + 1, 2)}),
      quad_omega_(new QuadOmega2D{CreateOmegaComposite(omega, GaussTensorSimplex(deg, 2))}),
      omega_(std::move(omega)),
      index_table_(deg, 2),
      fdiff_{},
      A_{},
      B_{},
      C_{},
      D_{},
      phi_pos_{},
      phi0_{} {
   GEN_QUAD_ASSERT_DEBUG(deg_ >= 1);

   fdiff_[0].resize(num_funcs_);
   fdiff_[1].resize(num_funcs_);
   fdiff_[2].resize(num_funcs_);
   fdiff_[3].resize(num_funcs_);

   phi_pos_.resize(deg_ + 2);
   phi_pos_[0] = 0;
   for(gq_int n = 1; n <= deg_ + 1; ++n) {
      phi_pos_[n] = (n * (n + 1)) / 2;
   }

   phi0_ = std::pow<double>(omega_.area(), -0.5);

   A_.resize(deg_);
   B_.resize(deg_);
   C_.resize(deg_);
   D_.resize(deg_);
   for(gq_int n = 0; n < A_.size(); ++n) {
      A_[n].resize((n + 1), (n + 2));
      B_[n].resize((n + 1), (n + 2));
      C_[n].resize((n + 1), (n + 2));
      if(n > 0) {
         D_[n].resize(n, (n + 2));
      }
      this->recurrence_coeffs(n);
   }

   StandardBasisIndices(index_table_);
}


OmegaBasis2D::OmegaBasis2D(const OmegaBasis2D& b)
    : Basis(b),
      quad_simplex_(new QuadSimplex{*b.quad_simplex_}),
      quad_omega_(new QuadOmega2D{*b.quad_omega_}),
      omega_(b.omega_),
      index_table_(b.index_table_) {
   for(gq_int i = 0; i < 4; ++i) {
      fdiff_[i].resize(b.fdiff_[i].size());
   }
   A_ = b.A_;
   B_ = b.B_;
   C_ = b.C_;
   D_ = b.D_;
   phi_pos_ = b.phi_pos_;
   phi0_ = b.phi0_;
}


OmegaBasis2D::~OmegaBasis2D() {
   delete quad_simplex_;
   delete quad_omega_;
}


const Array1D& OmegaBasis2D::orthog_basis(const Array1D& x) {
   this->orthog_basis_internal(x, functions_);
   return functions_;
}


const Array2D& OmegaBasis2D::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives_;
}


const Array1D& OmegaBasis2D::orthog_integrals() {
   integrals_ = 0.;
   integrals_[0] = 1. / phi0_;
   return integrals_;
}


const Array1D& OmegaBasis2D::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table_, x, functions_);
   return functions_;
}


const Array2D& OmegaBasis2D::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table_, x, derivatives_);
   return derivatives_;
}


const Array1D& OmegaBasis2D::monomial_integrals() {
   for(gq_int k = 0; k < num_funcs_; ++k) {
      auto monomial_k = [k, this](const auto& x) {
         return std::pow<double>(x[0], index_table_(0, k)) * std::pow<double>(x[1], index_table_(1, k));
      };
      integrals_[k] = quad_approx(*quad_omega_, monomial_k);
   }

   return integrals_;
}


void OmegaBasis2D::recurrence_coeffs(gq_int n) {
   Matrix2D ALoc(2 * n + 1, 2 * n + 2);
   Matrix2D MLoc(2 * n + 2, 2 * n + 2);

   Array1D phi((n + 1) * (n + 2) / 2);
   double* phin = &phi[phi_pos_[n]];

   StaticArray1D<2> x;
   double* L = &x[0];
   if(n > 0) {
      L = &phi[1];
   }

   for(gq_int nt = 0; nt < omega_.triang.size(); ++nt) {
      for(gq_int np = 0; np < quad_simplex_->num_nodes(); ++np) {
         const auto w = quad_simplex_->w(np) * omega_.triang[nt].jacobian();
         map_from_unit(omega_.triang[nt], quad_simplex_->node(np), x);

         this->orthog_basis_internal_p(n, x, phi);

         for(gq_int l = 0; l <= n; ++l) {
            StaticArrayInt1D<2> row_m{l, l + n + 1};
            StaticArray1D<2> chi{L[0] * phin[l], L[1] * phin[l]};

            StaticArrayInt1D<2> row_a;
            StaticArray1D<2> phi_loc;
            gq_int sz;
            if(l < n) {
               double* phin1 = &phi[phi_pos_[n - 1]];
               phi_loc[0] = phin1[l];
               phi_loc[1] = phin[l];
               row_a[0] = l;
               row_a[1] = l + n;
               sz = 2;
            } else {
               phi_loc[0] = phin[l];
               row_a[0] = l + n;
               sz = 1;
            }

            for(gq_int lp = 0; lp <= n; ++lp) {
               StaticArray1D<2> chi_p{L[0] * phin[lp], L[1] * phin[lp]};
               StaticArrayInt1D<2> col{lp, lp + n + 1};

               for(gq_int i = 0; i < 2; ++i) {
                  for(gq_int ip = 0; ip < 2; ++ip) {
                     MLoc(row_m[i], col[ip]) += chi[i] * chi_p[ip] * w;
                  }
               }

               for(gq_int i = 0; i < sz; ++i) {
                  for(gq_int ip = 0; ip < 2; ++ip) {
                     ALoc(row_a[i], col[ip]) += phi_loc[i] * chi_p[ip] * w;
                  }
               }
            }
         }
      }
   }

   MLoc -= ALoc.transpose() * ALoc;
   Eigen::SelfAdjointEigenSolver<Matrix2D> eig_solve{MLoc};
   auto eig_val = eig_solve.eigenvalues();
   auto eig_vec = eig_solve.eigenvectors();

   Matrix2D PLoc(2 * A_[n].rows(), A_[n].cols());
   for(gq_int k = 0; k < PLoc.cols(); ++k) {
      PLoc.col(k) = eig_vec.col(2 * n + 1 - k) * std::pow(eig_val[2 * n + 1 - k], -0.5);
   }

   Matrix2D QLoc = -(ALoc * PLoc);

   A_[n] = PLoc(seqN(0, A_[n].rows()), all);
   B_[n] = PLoc(seqN(n + 1, B_[n].rows()), all);

   if(n > 0) {
      D_[n] = QLoc(seqN(0, D_[n].rows()), all);
      C_[n] = QLoc(seqN(n, C_[n].rows()), all);
   } else {
      C_[n] = QLoc;
   }
}


void OmegaBasis2D::test_orthogonal(gq_int p, bool verbose) {
   GEN_QUAD_ASSERT_DEBUG(p <= deg_);

   gq_int N = ((p + 2) * (p + 1)) / 2;
   Array1D x(2);
   Array1D phi(N);
   Matrix2D I(N, N);

   for(gq_int nt = 0; nt < omega_.triang.size(); ++nt) {
      for(gq_int np = 0; np < quad_simplex_->num_nodes(); ++np) {
         const auto w = quad_simplex_->w(np) * omega_.triang[nt].jacobian();
         map_from_unit(omega_.triang[nt], quad_simplex_->node(np), x);

         this->orthog_basis_internal_p(p, x, phi);
         for(gq_int i = 0; i < N; ++i) {
            for(gq_int j = 0; j < N; ++j) {
               I(i, j) += phi[i] * phi[j] * w;
            }
         }
      }
   }

   if(verbose) {
      std::cout << I << std::endl;
   }

   Eigen::SelfAdjointEigenSolver<Matrix2D> eig_solve{I};
   auto eig_val = eig_solve.eigenvalues();
   double cnd_num = eig_val[N - 1] / eig_val[0];
   std::cout << "Condition number of mass matrix using orthogonal polynomials: " << cnd_num << std::endl;
}


void OmegaBasis2D::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   return this->orthog_basis_internal_p(deg_, x, buff);
}


void OmegaBasis2D::orthog_basis_internal_p(gq_int p, const Array1D& x, Array1D& buff) {
   GEN_QUAD_ASSERT_DEBUG(p <= deg_);
   buff[0] = phi0_;
   if(p == 0) {
      return;
   }

   buff[1] = (A_[0](0, 0) * x[0] + B_[0](0, 0) * x[1] + C_[0](0, 0)) * buff[0];
   buff[2] = (A_[0](0, 1) * x[0] + B_[0](0, 1) * x[1] + C_[0](0, 1)) * buff[0];

   for(gq_int n = 1; n < p; ++n) {
      const auto yp = buff(seqN(phi_pos_[n - 1], D_[n].rows()));
      const auto yc = buff(seqN(phi_pos_[n], A_[n].rows()));
      auto yn = buff(seqN(phi_pos_[n + 1], A_[n].cols()));
      yn = ((A_[n].transpose() * buff[1] + B_[n].transpose() * buff[2] + C_[n].transpose()) * yc.matrix())
               .array();
      yn += (D_[n].transpose() * yp.matrix()).array();
   }
}


std::unique_ptr<Basis> OmegaBasis2D::clone() const {
   return std::make_unique<OmegaBasis2D>(*this);
}


gq_int Omega2DBasisSize(gq_int deg) {
   GEN_QUAD_ASSERT_DEBUG(deg >= 1);
   return StandardBasisSize(deg, 2);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
BasisTable::BasisTable(gq_int deg, gq_int dim)
    : deg_{deg},
      dim_{dim},
      num_elem_{StandardBasisSize(deg, dim)},
      data_(dim_, num_elem_) {
   GEN_QUAD_ASSERT_DEBUG(deg_ >= 0);
   GEN_QUAD_ASSERT_DEBUG(dim_ >= 1);
}


BasisTable::BasisTable(gq_int deg, gq_int dim, gq_int num_elem)
    : deg_{deg},
      dim_{dim},
      num_elem_{num_elem},
      data_(dim_, num_elem_) {
   GEN_QUAD_ASSERT_DEBUG(deg_ >= 0);
   GEN_QUAD_ASSERT_DEBUG(dim_ >= 1);
   GEN_QUAD_ASSERT_DEBUG(num_elem_ >= 1);
}


// print transpose
std::ostream& operator<<(std::ostream& os, const BasisTable& table) {
   return os << table.array().transpose();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
namespace internal {


template <typename BasisType>
void fourth_order_difference(BasisType& basis, const Array1D& x) {
   const gq_int dim{basis.dim()};
   const gq_int num_funcs{basis.size()};
   constexpr double h{5.e-5};

   Array1D x_backw2(dim);
   Array1D x_backw1(dim);
   Array1D x_forw1(dim);
   Array1D x_forw2(dim);

   for(gq_int d = 0; d < dim; ++d) {
      x_backw1 = x;
      x_backw2 = x;
      x_forw1 = x;
      x_forw2 = x;

      x_backw2[d] -= 2. * h;
      x_backw1[d] -= h;
      x_forw1[d] += h;
      x_forw2[d] += 2. * h;

      basis.orthog_basis_internal(x_backw2, basis.fdiff_[0]);
      basis.orthog_basis_internal(x_backw1, basis.fdiff_[1]);
      basis.orthog_basis_internal(x_forw1, basis.fdiff_[2]);
      basis.orthog_basis_internal(x_forw2, basis.fdiff_[3]);
      for(gq_int k = 0; k < num_funcs; ++k) {
         basis.derivatives_(d, k) = (1. * basis.fdiff_[0][k] - 8. * basis.fdiff_[1][k]
                                     + 8. * basis.fdiff_[2][k] - 1. * basis.fdiff_[3][k])
                                  / (12. * h);
      }
   }
}


}  // namespace internal


// computes fact(deg+dim)/( fact(deg) * fact(dim) )
// special care is taken to avoid integer overflow
gq_int StandardBasisSize(gq_int deg, gq_int dim) {
   gq_int product = 1;
   for(gq_int dd = deg + 1, dimc = 1; dd <= deg + dim; ++dd, ++dimc) {
      product = product * dd / dimc;  // can prove that (product * dd) is divisible by dimc
   }
   return product;
}


void StandardBasisIndices(BasisTable& table) {
   gq_int deg{table.deg()};
   gq_int dim{table.dim()};

   if(dim == 1) {
      for(gq_int i = 0; i <= deg; ++i) {
         table(0, i) = i;
      }
      return;
   }

   // compute basis indices using nested recursion if dimension >= 2
   gq_int elem_count = 0;
   for(gq_int deg_cur = 0; deg_cur <= deg; ++deg_cur) {
      BasisTable table_cur(deg - deg_cur, dim - 1);
      StandardBasisIndices(table_cur);

      for(gq_int d = 0; d < dim - 1; ++d) {
         for(gq_int k = 0; k < table_cur.num_elem(); ++k) {
            table(d, elem_count + k) = table_cur(d, k);
         }
      }
      for(gq_int k = 0; k < table_cur.num_elem(); ++k) {
         table(dim - 1, elem_count + k) = deg_cur;
      }

      elem_count += table_cur.num_elem();
   }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
template <typename EigenType>
static inline void Monomial(double x, EigenType&& p) {
   p[0] = 1.;
   for(gq_int i = 1; i < p.size(); ++i) {
      p[i] = std::pow<double>(x, i);
   }
}


template <typename EigenType>
static inline void Legendre(double x, EigenType&& p) {
   long double x_map = 2.L * x - 1.L;

   long double fac1, fac2;
   long double prev = 1.L, cur = x_map, next;
   p[0] = double(prev);
   p[1] = double(cur);

   for(gq_int k = 1; k < p.size() - 1; ++k) {
      fac1 = (2.L * k + 1.L) / (k + 1.L);
      fac2 = k / (k + 1.L);
      next = fac1 * x_map * cur - fac2 * prev;
      p[k + 1] = double(next);

      prev = cur;
      cur = next;
   }
}


template <typename EigenType>
static void Jacobi(gq_int alpha, double x, EigenType&& p) {
   long double prev = 1.L;
   long double cur = 0.5L * (x - 1.L) * (alpha + 2.L) + alpha + 1.L;
   long double next;

   p[0] = double(prev);
   p[1] = double(cur);

   long double fac1Part1, fac1Part2, fac2, fac3, c1, c2, c3;
   gq_int k = 0;
   while(++k < p.size() - 1) {
      fac1Part1 = (alpha + 2.L * k + 1.L) * (alpha + 2.L * k + 2.L) * (alpha + 2.L * k);
      fac1Part2 = (alpha + 2.L * k + 1.L) * (alpha * alpha);
      fac2 = -2.L * (alpha + k) * k * (alpha + 2.L * k + 2.L);
      fac3 = 1.L / (2.L * (k + 1.L) * (alpha + k + 1.L) * (alpha + 2.L * k));
      c1 = fac3 * fac1Part1;
      c2 = fac3 * fac1Part2;
      c3 = fac3 * fac2;

      next = (c1 * x + c2) * cur + c3 * prev;
      p[k + 1] = double(next);

      prev = cur;
      cur = next;
   }
}


static void StandardMonomialFunctions(const BasisTable& index_table, const Array1D& x, Array1D& functions) {
   gq_int deg{index_table.deg()};
   gq_int dim{index_table.dim()};
   gq_int num_funcs{index_table.num_elem()};

   Array2D monomial(dim, deg + 1);
   for(gq_int d = 0; d < dim; ++d) {
      Monomial(x[d], monomial.row(d));
   }

   functions = 1.;
   for(gq_int d = 0; d < dim; ++d) {
      for(gq_int k = 0; k < num_funcs; ++k) {
         functions[k] *= monomial(d, index_table(d, k));
      }
   }
}


static void StandardMonomialDerivatives(const BasisTable& index_table, const Array1D& x,
                                        Array2D& derivatives) {
   gq_int dim{index_table.dim()};
   gq_int num_funcs{index_table.num_elem()};
   constexpr double h{5.e-5};

   Array1D fdiff[4];
   for(auto& fd : fdiff) {
      fd.resize(num_funcs);
   }

   Array1D x_backw2(dim);
   Array1D x_backw1(dim);
   Array1D x_forw1(dim);
   Array1D x_forw2(dim);

   for(gq_int d = 0; d < dim; ++d) {
      x_backw1 = x;
      x_backw2 = x;
      x_forw1 = x;
      x_forw2 = x;

      x_backw2[d] -= 2. * h;
      x_backw1[d] -= h;
      x_forw1[d] += h;
      x_forw2[d] += 2. * h;

      StandardMonomialFunctions(index_table, x_backw2, fdiff[0]);
      StandardMonomialFunctions(index_table, x_backw1, fdiff[1]);
      StandardMonomialFunctions(index_table, x_forw1, fdiff[2]);
      StandardMonomialFunctions(index_table, x_forw2, fdiff[3]);
      for(gq_int k = 0; k < num_funcs; ++k) {
         derivatives(d, k)
             = (1. * fdiff[0][k] - 8. * fdiff[1][k] + 8. * fdiff[2][k] - 1. * fdiff[3][k]) / (12. * h);
      }
   }
}


}  // namespace gquad

#endif

