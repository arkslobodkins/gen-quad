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
    : m_deg{deg},
      m_dim{dim},
      m_num_funcs{num_funcs},
      functions(num_funcs),
      derivatives(dim, num_funcs),
      integrals(num_funcs) {
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
PolytopeBasis::PolytopeBasis(gq_int deg, gq_int dim, gq_int num_funcs) : Basis(deg, dim, num_funcs) {
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cube Basis
CubeBasis::CubeBasis(gq_int deg, gq_int dim)
    : PolytopeBasis(deg, dim, CubeBasisSize(deg, dim)),
      index_table(deg, dim),
      fdiff{} {
   GEN_QUAD_ASSERT_DEBUG(m_deg >= 1);
   GEN_QUAD_ASSERT_DEBUG(m_dim >= 2);

   fdiff[0].resize(m_num_funcs);
   fdiff[1].resize(m_num_funcs);
   fdiff[2].resize(m_num_funcs);
   fdiff[3].resize(m_num_funcs);
   StandardBasisIndices(index_table);
}

const Array1D& CubeBasis::orthog_basis(const Array1D& x) {
   orthog_basis_internal(x, functions);
   return functions;
}

void CubeBasis::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   Array2D legendre(m_dim, m_deg + 1);
   for(gq_int d = 0; d < m_dim; ++d) {
      Legendre(x[d], legendre.row(d));
   }

   for(gq_int k = 0; k < m_num_funcs; ++k) {
      buff[k] = legendre(0, index_table(0, k));
   }
   for(gq_int d = 1; d < m_dim; ++d) {
      for(gq_int k = 0; k < m_num_funcs; ++k) {
         buff[k] *= legendre(d, index_table(d, k));
      }
   }
}

const Array2D& CubeBasis::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives;
}

const Array1D& CubeBasis::orthog_integrals() {
   integrals = 0.;
   integrals[0] = 1.;
   return integrals;
}

const Array1D& CubeBasis::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table, x, functions);
   return functions;
}

const Array2D& CubeBasis::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table, x, derivatives);
   return derivatives;
}

const Array1D& CubeBasis::monomial_integrals() {
   for(gq_int k = 0; k < m_num_funcs; ++k) {
      gq_int prod = 1;  // use integer instead of double and perform multiplication instead of division
                        // to avoid error associated with intermediary divisions
      for(gq_int d = 0; d < m_dim; ++d) {
         prod *= (index_table(d, k) + 1);
      }
      integrals[k] = 1. / prod;
   }
   return integrals;
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
      power_table(deg, dim),
      index_table(deg, dim),
      fdiff{} {
   GEN_QUAD_ASSERT_DEBUG(m_deg >= 1);
   GEN_QUAD_ASSERT_DEBUG(m_dim >= 2);

   fdiff[0].resize(m_num_funcs);
   fdiff[1].resize(m_num_funcs);
   fdiff[2].resize(m_num_funcs);
   fdiff[3].resize(m_num_funcs);

   StandardBasisIndices(index_table);
   for(gq_int d = 1; d < m_dim; ++d) {
      for(gq_int i = 0; i < d; ++i) {
         for(gq_int k = 0; k < m_num_funcs; ++k) {
            power_table(d, k) += index_table(i, k);
         }
      }
   }
}

const Array1D& SimplexBasis::orthog_basis(const Array1D& x) {
   orthog_basis_internal(x, functions);
   return functions;
}

const Array2D& SimplexBasis::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives;
}

const Array1D& SimplexBasis::orthog_integrals() {
   integrals = 0.;
   integrals[0] = 1. / math::factorial(m_dim);
   return integrals;
}

const Array1D& SimplexBasis::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table, x, functions);
   return functions;
}

const Array2D& SimplexBasis::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table, x, derivatives);
   return derivatives;
}

const Array1D& SimplexBasis::monomial_integrals() {
   for(gq_int k = 0; k < m_num_funcs; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < m_dim; ++d) {
         gq_int s = index_table.array()(seq(d, m_dim - 1), k).sum();
         prod *= (s + m_dim - d);
      }
      integrals[k] = 1. / prod;
   }
   return integrals;
}

void SimplexBasis::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   Array1D legendre(m_deg + 1);
   Legendre((x[m_dim - 2] - x[m_dim - 1]) / x[m_dim - 2], legendre);  // mapped (2y-x)/x to (x-y)/x

   Array1DLong xCoord(m_dim - 1);
   Array1DLong jCoord(m_dim - 1);
   for(gq_int d = 0; d < m_dim - 2; ++d) {
      xCoord[d] = x[m_dim - d - 2] / x[m_dim - d - 3];
      jCoord[d] = 1.L - 2.L * x[m_dim - d - 2] / x[m_dim - d - 3];
   }
   xCoord[m_dim - 2] = x[0];
   jCoord[m_dim - 2] = 1.L - 2.L * x[0];

   Array3D jacobi(m_dim, m_deg + 1, m_deg + 1);
   for(gq_int d = 1; d < m_dim; ++d) {
      for(gq_int j = 0; j < m_deg + 1; ++j) {
         Jacobi(2 * j + d, jCoord[d - 1], jacobi(d, j));
      }
   }

   Array2D map_table(m_dim, m_deg + 1);
   for(gq_int d = 1; d < map_table.rows(); ++d) {
      for(gq_int j = 0; j < map_table.cols(); ++j) {
         map_table(d, j) = std::pow<double>(xCoord[d - 1], j);
      }
   }

   for(gq_int k = 0; k < m_num_funcs; ++k) {
      buff[k] = legendre[index_table(0, k)];
   }

   for(gq_int d = 1; d < m_dim; ++d) {
      for(gq_int k = 0; k < m_num_funcs; ++k) {
         buff[k] *= jacobi(d, power_table(d, k), index_table(d, k)) * map_table(d, power_table(d, k));
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
      m_dims{dim1, dim2},
      power_table_second(deg, dim2, CubeSimplexBasisSize(deg, dim1 + dim2)),
      index_table(deg, dim1 + dim2),
      fdiff{} {
   GEN_QUAD_ASSERT_DEBUG(m_deg >= 1);
   GEN_QUAD_ASSERT_DEBUG(m_dims[0] >= 1);
   GEN_QUAD_ASSERT_DEBUG(m_dims[1] >= 2);

   fdiff[0].resize(m_num_funcs);
   fdiff[1].resize(m_num_funcs);
   fdiff[2].resize(m_num_funcs);
   fdiff[3].resize(m_num_funcs);

   StandardBasisIndices(index_table);
   for(gq_int d = 1; d < m_dims[1]; ++d) {
      for(gq_int i = 0; i < d; ++i) {
         for(gq_int k = 0; k < m_num_funcs; ++k) {
            power_table_second(d, k) += index_table(m_dims[0] + i, k);
         }
      }
   }
}

const Array1D& CubeSimplexBasis::orthog_basis(const Array1D& x) {
   orthog_basis_internal(x, functions);
   return functions;
}

const Array2D& CubeSimplexBasis::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives;
}

const Array1D& CubeSimplexBasis::orthog_integrals() {
   integrals = 0.;
   integrals[0] = 1. / math::factorial(m_dims[1]);
   return integrals;
}

const Array1D& CubeSimplexBasis::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table, x, functions);
   return functions;
}

const Array2D& CubeSimplexBasis::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table, x, derivatives);
   return derivatives;
}

const Array1D& CubeSimplexBasis::monomial_integrals() {
   // integrals due to cube
   for(gq_int k = 0; k < m_num_funcs; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < m_dims[0]; ++d) {
         prod *= (index_table(d, k) + 1);
      }
      integrals[k] = 1. / prod;
   }

   // integrals due to simplex
   for(gq_int k = 0; k < m_num_funcs; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < m_dims[1]; ++d) {
         gq_int s = index_table.array()(seq(m_dims[0] + d, m_dim - 1), k).sum();
         prod *= (s + m_dims[1] - d);
      }
      integrals[k] /= prod;
   }

   return integrals;
}

void CubeSimplexBasis::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   orthog_basis_polytopic_two_internal(x, buff);

   Array1D legendre(m_deg + 1);
   for(gq_int d = 0; d < m_dims[0]; ++d) {
      Legendre(x[d], legendre);
      for(gq_int k = 0; k < m_num_funcs; ++k) {
         buff[k] *= legendre[index_table(d, k)];
      }
   }
}

void CubeSimplexBasis::orthog_basis_polytopic_two_internal(const Array1D& x, Array1D& buff) {
   gq_int dim1{m_dims[0]};
   gq_int dim2{m_dims[1]};
   gq_int dim_two{dim1 + dim2};

   Array1D legendre(m_deg + 1);
   Legendre((x[dim_two - 2] - x[dim_two - 1]) / x[dim_two - 2], legendre);  // mapped (2y-x)/x to (x-y)/x

   Array1DLong xCoord(dim2 - 1);
   for(gq_int d = 0; d < dim2 - 2; ++d) {
      xCoord[d] = x[dim_two - d - 2] / x[dim_two - d - 3];
   }
   xCoord[dim2 - 2] = x[dim1];

   Array1DLong jCoord(dim2 - 1);
   for(gq_int d = 1; d < dim2 - 1; ++d) {
      jCoord[d - 1] = 1.L - 2.L * x[dim_two - d - 1] / x[dim_two - d - 2];
   }
   jCoord[dim2 - 2] = 1.L - 2.L * x[dim1];

   Array3D jacobi(dim2, m_deg + 1, m_deg + 1);
   for(gq_int d = 1; d < dim2; ++d) {
      for(gq_int j = 0; j < m_deg + 1; ++j) {
         Jacobi(2 * j + d, jCoord[d - 1], jacobi(d, j));
      }
   }

   Array2D map_table(dim2, m_deg + 1);
   for(gq_int d = 1; d < map_table.rows(); ++d) {
      for(gq_int j = 0; j < map_table.cols(); ++j) {
         map_table(d, j) = std::pow<double>(xCoord[d - 1], j);
      }
   }

   for(gq_int k = 0; k < m_num_funcs; ++k) {
      buff[k] = legendre[index_table(dim1, k)];
   }

   for(gq_int d = 1; d < dim2; ++d) {
      for(gq_int k = 0; k < m_num_funcs; ++k) {
         buff[k] *= jacobi(d, power_table_second(d, k), index_table(dim1 + d, k))
                  * map_table(d, power_table_second(d, k));
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
      m_dims{dim1, dim2},
      power_table_first(deg, dim1, SimplexSimplexBasisSize(deg, dim1 + dim2)),
      power_table_second(deg, dim2, SimplexSimplexBasisSize(deg, dim1 + dim2)),
      index_table(deg, dim1 + dim2),
      poly_buff1{},
      poly_buff2{},
      fdiff{} {
   GEN_QUAD_ASSERT_DEBUG(m_deg >= 1);
   GEN_QUAD_ASSERT_DEBUG(m_dims[0] >= 2);
   GEN_QUAD_ASSERT_DEBUG(m_dims[1] >= 2);

   poly_buff1.resize(m_num_funcs);
   poly_buff2.resize(m_num_funcs);
   fdiff[0].resize(m_num_funcs);
   fdiff[1].resize(m_num_funcs);
   fdiff[2].resize(m_num_funcs);
   fdiff[3].resize(m_num_funcs);

   StandardBasisIndices(index_table);

   for(gq_int d = 1; d < m_dims[0]; ++d) {
      for(gq_int i = 0; i < d; ++i) {
         for(gq_int k = 0; k < m_num_funcs; ++k) {
            power_table_first(d, k) += index_table(i, k);
         }
      }
   }

   for(gq_int d = 1; d < m_dims[1]; ++d) {
      for(gq_int i = 0; i < d; ++i) {
         for(gq_int k = 0; k < m_num_funcs; ++k) {
            power_table_second(d, k) += index_table(m_dims[0] + i, k);
         }
      }
   }
}

const Array1D& SimplexSimplexBasis::orthog_basis(const Array1D& x) {
   orthog_basis_internal(x, functions);
   return functions;
}

const Array2D& SimplexSimplexBasis::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives;
}

const Array1D& SimplexSimplexBasis::orthog_integrals() {
   integrals = 0.;
   integrals[0] = 1. / (math::factorial(m_dims[0]) * math::factorial(m_dims[1]));
   return integrals;
}

const Array1D& SimplexSimplexBasis::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table, x, functions);
   return functions;
}

const Array2D& SimplexSimplexBasis::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table, x, derivatives);
   return derivatives;
}

const Array1D& SimplexSimplexBasis::monomial_integrals() {
   // integrals due to simplex1
   for(gq_int k = 0; k < m_num_funcs; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < m_dims[0]; ++d) {
         gq_int s = index_table.array()(seq(d, m_dims[0] - 1), k).sum();
         prod *= (s + m_dims[0] - d);
      }
      integrals[k] = 1. / prod;
   }

   // integrals due to simplex2
   for(gq_int k = 0; k < m_num_funcs; ++k) {
      gq_int prod = 1;
      for(gq_int d = 0; d < m_dims[1]; ++d) {
         gq_int s = index_table.array()(seq(m_dims[0] + d, m_dim - 1), k).sum();
         prod *= (s + m_dims[1] - d);
      }
      integrals[k] /= prod;
   }

   return integrals;
}

void SimplexSimplexBasis::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   orthog_basis_polytopic_one_internal(x, poly_buff1);
   orthog_basis_polytopic_two_internal(x, poly_buff2);

   buff = 1.;
   buff *= poly_buff1;
   buff *= poly_buff2;
}

void SimplexSimplexBasis::orthog_basis_polytopic_one_internal(const Array1D& x, Array1D& buff) {
   gq_int dim1{m_dims[0]};

   Array1D legendre(m_deg + 1);
   Legendre((x[dim1 - 2] - x[dim1 - 1]) / x[dim1 - 2], legendre);  // mapped (2y-x)/x to (x-y)/x

   Array1DLong xCoord(dim1 - 1);
   for(gq_int d = 0; d < dim1 - 2; ++d) {
      xCoord[d] = x[dim1 - d - 2] / x[dim1 - d - 3];
   }
   xCoord[dim1 - 2] = x[0];

   Array1DLong jCoord(dim1 - 1);
   for(gq_int d = 1; d < dim1 - 1; ++d) {
      jCoord[d - 1] = 1.L - 2.L * x[dim1 - d - 1] / x[dim1 - d - 2];
   }
   jCoord[dim1 - 2] = 1.L - 2.L * x[0];

   Array3D jacobi(dim1, m_deg + 1, m_deg + 1);
   for(gq_int d = 1; d < dim1; ++d) {
      for(gq_int j = 0; j < m_deg + 1; ++j) {
         Jacobi(2 * j + d, jCoord[d - 1], jacobi(d, j));
      }
   }

   Array2D map_table(dim1, m_deg + 1);
   for(gq_int d = 1; d < map_table.rows(); ++d) {
      for(gq_int j = 0; j < map_table.cols(); ++j) {
         map_table(d, j) = std::pow<double>(xCoord[d - 1], j);
      }
   }

   for(gq_int k = 0; k < m_num_funcs; ++k) {
      buff[k] = legendre[index_table(0, k)];
   }

   for(gq_int d = 1; d < dim1; ++d) {
      for(gq_int k = 0; k < m_num_funcs; ++k) {
         buff[k]
             *= jacobi(d, power_table_first(d, k), index_table(d, k)) * map_table(d, power_table_first(d, k));
      }
   }
}

void SimplexSimplexBasis::orthog_basis_polytopic_two_internal(const Array1D& x, Array1D& buff) {
   gq_int dim1{m_dims[0]};
   gq_int dim2{m_dims[1]};
   gq_int dim_two{dim1 + dim2};

   Array1D legendre(m_deg + 1);
   Legendre((x[dim_two - 2] - x[dim_two - 1]) / x[dim_two - 2], legendre);  // mapped (2y-x)/x to (x-y)/x

   Array1DLong xCoord(dim2 - 1);
   for(gq_int d = 0; d < dim2 - 2; ++d) {
      xCoord[d] = x[dim_two - d - 2] / x[dim_two - d - 3];
   }
   xCoord[dim2 - 2] = x[dim1];

   Array1DLong jCoord(dim2 - 1);
   for(gq_int d = 1; d < dim2 - 1; ++d) {
      jCoord[d - 1] = 1.L - 2.L * x[dim_two - d - 1] / x[dim_two - d - 2];
   }
   jCoord[dim2 - 2] = 1.L - 2.L * x[dim1];

   Array3D jacobi(dim2, m_deg + 1, m_deg + 1);
   for(gq_int d = 1; d < dim2; ++d) {
      for(gq_int j = 0; j < m_deg + 1; ++j) {
         Jacobi(2 * j + d, jCoord[d - 1], jacobi(d, j));
      }
   }

   Array2D map_table(dim2, m_deg + 1);
   for(gq_int d = 1; d < map_table.rows(); ++d) {
      for(gq_int j = 0; j < map_table.cols(); ++j) {
         map_table(d, j) = std::pow<double>(xCoord[d - 1], j);
      }
   }

   for(gq_int k = 0; k < m_num_funcs; ++k) {
      buff[k] = legendre[index_table(dim1, k)];
   }

   for(gq_int d = 1; d < dim2; ++d) {
      for(gq_int k = 0; k < m_num_funcs; ++k) {
         buff[k] *= jacobi(d, power_table_second(d, k), index_table(dim1 + d, k))
                  * map_table(d, power_table_second(d, k));
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
      index_table(deg, 3),
      fdiff{} {
   GEN_QUAD_ASSERT_DEBUG(m_deg >= 1);

   fdiff[0].resize(m_num_funcs);
   fdiff[1].resize(m_num_funcs);
   fdiff[2].resize(m_num_funcs);
   fdiff[3].resize(m_num_funcs);
   StandardBasisIndices(index_table);
}

const Array1D& PyramidBasis3D::orthog_basis(const Array1D& x) {
   orthog_basis_internal(x, functions);
   return functions;
}

const Array2D& PyramidBasis3D::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives;
}

const Array1D& PyramidBasis3D::orthog_integrals() {
   integrals = 0.;
   integrals[0] = 1. / 3.;
   return integrals;
}

const Array1D& PyramidBasis3D::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table, x, functions);
   return functions;
}

const Array2D& PyramidBasis3D::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table, x, derivatives);
   return derivatives;
}

const Array1D& PyramidBasis3D::monomial_integrals() {
   for(gq_int k = 0; k < m_num_funcs; ++k) {
      integrals[k] = 1.
                   / ((index_table(1, k) + 1) * (index_table(2, k) + 1)
                      * (index_table(0, k) + index_table(1, k) + index_table(2, k) + 3));
   }
   return integrals;
}

void PyramidBasis3D::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   Array2D legendre(2, m_deg + 1);
   Legendre(1.L - (x[1] / x[0]), legendre.row(0));  // mapped 2(y/x)-1 to 1-y/x
   Legendre(1.L - (x[2] / x[0]), legendre.row(1));  // mapped 2(z/x)-1 to 1-z/x

   Array2D jacobi(m_deg + 1, m_deg + 1);
   for(gq_int j = 0; j < m_deg + 1; ++j) {
      Jacobi(2 * j + 2, 1.L - 2.L * x[0], jacobi.row(j));
   }

   Array1D xpower(m_deg + 1);
   Monomial(x[0], xpower);

   gq_int count = 0;
   for(gq_int k = 0; k < m_deg + 1; ++k) {
      for(gq_int j = 0; j < m_deg + 1 - k; ++j) {
         for(gq_int i = 0; i < m_deg + 1 - j - k; ++i) {
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
OmegaBasis2D::OmegaBasis2D(Omega2D omega_, gq_int deg)
    : Basis(deg, 2, Omega2DBasisSize(deg)),
      qts(new QuadSimplex{GaussTensorSimplex(2 * deg + 1, 2)}),
      qt_omega(new QuadOmega2D{CreateOmegaComposite(omega_, GaussTensorSimplex(deg, 2))}),
      omega(std::move(omega_)),
      index_table(deg, 2),
      fdiff{},
      A{},
      B{},
      C{},
      D{},
      phi_pos{},
      phi0{} {
   GEN_QUAD_ASSERT_DEBUG(m_deg >= 1);

   fdiff[0].resize(m_num_funcs);
   fdiff[1].resize(m_num_funcs);
   fdiff[2].resize(m_num_funcs);
   fdiff[3].resize(m_num_funcs);

   phi_pos.resize(m_deg + 2);
   phi_pos[0] = 0;
   for(gq_int n = 1; n <= m_deg + 1; ++n) {
      phi_pos[n] = (n * (n + 1)) / 2;
   }

   phi0 = std::pow<double>(omega.area(), -0.5);

   A.resize(m_deg);
   B.resize(m_deg);
   C.resize(m_deg);
   D.resize(m_deg);
   for(gq_int n = 0; n < A.size(); ++n) {
      A[n].resize((n + 1), (n + 2));
      B[n].resize((n + 1), (n + 2));
      C[n].resize((n + 1), (n + 2));
      if(n > 0) {
         D[n].resize(n, (n + 2));
      }
      recurrence_coeffs(n);
   }

   StandardBasisIndices(index_table);
}

OmegaBasis2D::OmegaBasis2D(const OmegaBasis2D& b)
    : Basis(b),
      qts(new QuadSimplex{*b.qts}),
      qt_omega(new QuadOmega2D{*b.qt_omega}),
      omega(b.omega),
      index_table(b.index_table) {
   for(gq_int i = 0; i < 4; ++i) {
      fdiff[i].resize(b.fdiff[i].size());
   }
   A = b.A;
   B = b.B;
   C = b.C;
   D = b.D;
   phi_pos = b.phi_pos;
   phi0 = b.phi0;
}

OmegaBasis2D::~OmegaBasis2D() {
   delete qts;
   delete qt_omega;
}

const Array1D& OmegaBasis2D::orthog_basis(const Array1D& x) {
   orthog_basis_internal(x, functions);
   return functions;
}

const Array2D& OmegaBasis2D::orthog_der(const Array1D& x) {
   internal::fourth_order_difference(*this, x);
   return derivatives;
}

const Array1D& OmegaBasis2D::orthog_integrals() {
   integrals = 0.;
   integrals[0] = 1. / phi0;
   return integrals;
}

const Array1D& OmegaBasis2D::monomial_basis(const Array1D& x) {
   StandardMonomialFunctions(index_table, x, functions);
   return functions;
}

const Array2D& OmegaBasis2D::monomial_der(const Array1D& x) {
   StandardMonomialDerivatives(index_table, x, derivatives);
   return derivatives;
}

const Array1D& OmegaBasis2D::monomial_integrals() {
   for(gq_int k = 0; k < m_num_funcs; ++k) {
      auto monomial_k = [k, this](const auto& x) {
         return std::pow<double>(x[0], index_table(0, k)) * std::pow<double>(x[1], index_table(1, k));
      };
      integrals[k] = quad_approx(*qt_omega, monomial_k);
   }

   return integrals;
}

void OmegaBasis2D::recurrence_coeffs(gq_int n) {
   Matrix2D ALoc(2 * n + 1, 2 * n + 2);
   Matrix2D MLoc(2 * n + 2, 2 * n + 2);

   Array1D phi((n + 1) * (n + 2) / 2);
   double* phin = &phi[phi_pos[n]];

   StaticArray1D<2> x;
   double* L = &x[0];
   if(n > 0) {
      L = &phi[1];
   }

   for(gq_int nt = 0; nt < omega.triang.size(); ++nt) {
      for(gq_int np = 0; np < qts->num_nodes(); ++np) {
         const auto w = qts->w(np) * omega.triang[nt].jacobian();
         map_from_unit(omega.triang[nt], qts->node(np), x);

         orthog_basis_internal_p(n, x, phi);

         for(gq_int l = 0; l <= n; ++l) {
            StaticArrayInt1D<2> row_m{l, l + n + 1};
            StaticArray1D<2> chi{L[0] * phin[l], L[1] * phin[l]};

            StaticArrayInt1D<2> row_a;
            StaticArray1D<2> phi_loc;
            gq_int sz;
            if(l < n) {
               double* phin1 = &phi[phi_pos[n - 1]];
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

   Matrix2D PLoc(2 * A[n].rows(), A[n].cols());
   for(gq_int k = 0; k < PLoc.cols(); ++k) {
      PLoc.col(k) = eig_vec.col(2 * n + 1 - k) * std::pow(eig_val[2 * n + 1 - k], -0.5);
   }

   Matrix2D QLoc = -(ALoc * PLoc);

   A[n] = PLoc(seqN(0, A[n].rows()), all);
   B[n] = PLoc(seqN(n + 1, B[n].rows()), all);

   if(n > 0) {
      D[n] = QLoc(seqN(0, D[n].rows()), all);
      C[n] = QLoc(seqN(n, C[n].rows()), all);
   } else {
      C[n] = QLoc;
   }
}

void OmegaBasis2D::TestOrthogonal(gq_int p, bool verbose) {
   GEN_QUAD_ASSERT_DEBUG(p <= m_deg);

   gq_int N = ((p + 2) * (p + 1)) / 2;
   Array1D x(2);
   Array1D phi(N);
   Matrix2D I(N, N);

   for(gq_int nt = 0; nt < omega.triang.size(); ++nt) {
      for(gq_int np = 0; np < qts->num_nodes(); ++np) {
         const auto w = qts->w(np) * omega.triang[nt].jacobian();
         map_from_unit(omega.triang[nt], qts->node(np), x);

         orthog_basis_internal_p(p, x, phi);

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
   double cond = eig_val[N - 1] / eig_val[0];
   std::cout << "Condition number of mass matrix using orthogonal polynomials: " << cond << std::endl;
}

void OmegaBasis2D::orthog_basis_internal(const Array1D& x, Array1D& buff) {
   return orthog_basis_internal_p(m_deg, x, buff);
}

void OmegaBasis2D::orthog_basis_internal_p(gq_int p, const Array1D& x, Array1D& buff) {
   GEN_QUAD_ASSERT_DEBUG(p <= m_deg);
   buff[0] = phi0;
   if(p == 0) {
      return;
   }

   buff[1] = (A[0](0, 0) * x[0] + B[0](0, 0) * x[1] + C[0](0, 0)) * buff[0];
   buff[2] = (A[0](0, 1) * x[0] + B[0](0, 1) * x[1] + C[0](0, 1)) * buff[0];

   for(gq_int n = 1; n < p; ++n) {
      const auto yp = buff(seqN(phi_pos[n - 1], D[n].rows()));
      const auto yc = buff(seqN(phi_pos[n], A[n].rows()));
      auto yn = buff(seqN(phi_pos[n + 1], A[n].cols()));
      yn = ((A[n].transpose() * buff[1] + B[n].transpose() * buff[2] + C[n].transpose()) * yc.matrix())
               .array();
      yn += (D[n].transpose() * yp.matrix()).array();
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
    : m_deg{deg},
      m_dim{dim},
      m_num_elem{StandardBasisSize(deg, dim)},
      data(m_dim, m_num_elem) {
   GEN_QUAD_ASSERT_DEBUG(m_deg >= 0);
   GEN_QUAD_ASSERT_DEBUG(m_dim >= 1);
}

BasisTable::BasisTable(gq_int deg, gq_int dim, gq_int num_elem)
    : m_deg{deg},
      m_dim{dim},
      m_num_elem{num_elem},
      data(m_dim, m_num_elem) {
   GEN_QUAD_ASSERT_DEBUG(m_deg >= 0);
   GEN_QUAD_ASSERT_DEBUG(m_dim >= 1);
   GEN_QUAD_ASSERT_DEBUG(m_num_elem >= 1);
}

BasisTable::BasisTable(const BasisTable& t) : BasisTable{t.deg(), t.dim(), t.num_elem()} {
   data = t.data;
}

std::ostream& operator<<(std::ostream& os, const BasisTable& t) {
   for(gq_int i = 0; i < t.num_elem(); ++i) {
      for(gq_int j = 0; j < t.dim(); ++j) {
         os << t(j, i) << "  ";
      }
      os << "\n";
   }
   return os;
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

      basis.orthog_basis_internal(x_backw2, basis.fdiff[0]);
      basis.orthog_basis_internal(x_backw1, basis.fdiff[1]);
      basis.orthog_basis_internal(x_forw1, basis.fdiff[2]);
      basis.orthog_basis_internal(x_forw2, basis.fdiff[3]);
      for(gq_int k = 0; k < num_funcs; ++k) {
         basis.derivatives(d, k) = (1. * basis.fdiff[0][k] - 8. * basis.fdiff[1][k] + 8. * basis.fdiff[2][k]
                                    - 1. * basis.fdiff[3][k])
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
   return static_cast<gq_int>(product);
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
   for(auto& v : fdiff) {
      v.resize(num_funcs);
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

