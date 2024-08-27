// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022
// Classes for computing orthogonal and polynomial bases, their derivatives, and their integrals.

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "gen_quad.hpp"
#include "mesh.hpp"

namespace gquad {


// BasisTable class is a matrix of type int. Its purpose is to store
// indexes that are used throughout the entire duration of Basis objects,
// thus avoiding recomputing them when Basis routines are called.
class BasisTable {
public:
   // number of rows will be set to the number of basis functions
   explicit BasisTable(gq_int deg, gq_int dim);

   // number of rows will be set to num_elem
   explicit BasisTable(gq_int deg, gq_int dim, gq_int num_elem);

   BasisTable(const BasisTable&) = default;
   BasisTable& operator=(const BasisTable) = delete;

   int& operator()(gq_int d, gq_int j) {
      return data_(d, j);
   }

   const int& operator()(gq_int d, gq_int j) const {
      return data_(d, j);
   }

   gq_int deg() const {
      return deg_;
   }

   gq_int dim() const {
      return dim_;
   }

   gq_int num_elem() const {
      return num_elem_;
   }

   const auto& array() const {
      return data_;
   }

private:
   const gq_int deg_;
   const gq_int dim_;
   const gq_int num_elem_;

   ArrayInt2D data_;  // dim x num_elem table
};


gq_int StandardBasisSize(gq_int deg, gq_int dim);
void StandardBasisIndices(BasisTable& table);


std::ostream& operator<<(std::ostream& os, const BasisTable& table);


namespace internal {
template <typename BasisType>
void fourth_order_difference(BasisType& basis, const Array1D& x);
}


// Abstract base class for all bases.
// Functions, their derivatives, and their integrals are
// implemented for both orthogonal and monomial polynomials.

// Any basis contains information about its degree, dimension, and number of basis functions,
// as well as arrays for computing basis functions, derivatives, and integrals.
// Currently they are initialized in this base class(subject to change in the future).
class Basis {
public:
   virtual ~Basis() = default;
   virtual std::unique_ptr<Basis> clone() const = 0;

   virtual const Array1D& orthog_basis(const Array1D& x) = 0;
   virtual const Array2D& orthog_der(const Array1D& x) = 0;
   virtual const Array1D& orthog_integrals() = 0;

   virtual const Array1D& monomial_basis(const Array1D& x) = 0;
   virtual const Array2D& monomial_der(const Array1D& x) = 0;
   virtual const Array1D& monomial_integrals() = 0;

   gq_int deg() const {
      return deg_;
   }

   gq_int dim() const {
      return dim_;
   }

   gq_int size() const {
      return num_funcs_;
   }

protected:
   const gq_int deg_;
   const gq_int dim_;
   const gq_int num_funcs_;

   Array1D functions_;
   Array2D derivatives_;
   Array1D integrals_;

   Basis(gq_int deg, gq_int dim, gq_int num_funcs);
   Basis(const Basis&) = default;
   Basis& operator=(const Basis&) = delete;
};


// Abstract base class for all polytopes.
class PolytopeBasis : public Basis {
public:
   virtual ~PolytopeBasis() = default;

protected:
   PolytopeBasis(gq_int deg, gq_int dim, gq_int num_funcs);
   PolytopeBasis(const PolytopeBasis&) = default;
   PolytopeBasis& operator=(const PolytopeBasis&) = delete;
};


class CubeBasis : public PolytopeBasis {
public:
   explicit CubeBasis(gq_int deg, gq_int dim);
   CubeBasis(const CubeBasis&) = default;
   std::unique_ptr<Basis> clone() const override;

   CubeBasis& operator=(const CubeBasis&) = delete;
   ~CubeBasis() override = default;

   const Array1D& orthog_basis(const Array1D& x) final;
   const Array2D& orthog_der(const Array1D& x) final;
   const Array1D& orthog_integrals() final;

   const Array1D& monomial_basis(const Array1D& x) final;
   const Array2D& monomial_der(const Array1D& x) final;
   const Array1D& monomial_integrals() final;

private:
   BasisTable index_table_;
   Array1D fdiff_[4];

   void orthog_basis_internal(const Array1D& x, Array1D& buff);

   template <typename BasisType>
   friend void internal::fourth_order_difference(BasisType& basis, const Array1D& x);
};


gq_int CubeBasisSize(gq_int deg, gq_int dim);


// Besides precomputing index_table, all Bases that contain Simplex
// also precompute power_table, which are used to efficiently compute
// orthogonal Simplex functions.
class SimplexBasis : public PolytopeBasis {
public:
   explicit SimplexBasis(gq_int deg, gq_int dim);
   SimplexBasis(const SimplexBasis&) = default;
   std::unique_ptr<Basis> clone() const override;

   SimplexBasis& operator=(const SimplexBasis&) = delete;
   ~SimplexBasis() override = default;

   const Array1D& orthog_basis(const Array1D& x) final;
   const Array2D& orthog_der(const Array1D& x) final;
   const Array1D& orthog_integrals() final;

   const Array1D& monomial_basis(const Array1D& x) final;
   const Array2D& monomial_der(const Array1D& x) final;
   const Array1D& monomial_integrals() final;

private:
   BasisTable power_table_;
   BasisTable index_table_;
   Array1D fdiff_[4];

   void orthog_basis_internal(const Array1D& x, Array1D& buff);

   template <typename BasisType>
   friend void internal::fourth_order_difference(BasisType& basis, const Array1D& x);
};


gq_int SimplexBasisSize(gq_int deg, gq_int dim);


class CubeSimplexBasis : public PolytopeBasis {
public:
   explicit CubeSimplexBasis(gq_int deg, gq_int dim1, gq_int dim2);
   CubeSimplexBasis(const CubeSimplexBasis&) = default;
   std::unique_ptr<Basis> clone() const override;

   CubeSimplexBasis& operator=(const CubeSimplexBasis&) = delete;
   ~CubeSimplexBasis() override = default;

   const Array1D& orthog_basis(const Array1D& x) final;
   const Array2D& orthog_der(const Array1D& x) final;
   const Array1D& orthog_integrals() final;

   const Array1D& monomial_basis(const Array1D& x) final;
   const Array2D& monomial_der(const Array1D& x) final;
   const Array1D& monomial_integrals() final;

private:
   std::array<gq_int, 2> dims_;
   BasisTable power_table_second_;
   BasisTable index_table_;
   Array1D fdiff_[4];

   void orthog_basis_internal(const Array1D& x, Array1D& buff);
   void orthog_basis_polytopic_two_internal(const Array1D& x, Array1D& buff);

   template <typename BasisType>
   friend void internal::fourth_order_difference(BasisType& basis, const Array1D& x);
};


gq_int CubeSimplexBasisSize(gq_int deg, gq_int dim);


class SimplexSimplexBasis : public PolytopeBasis {
public:
   explicit SimplexSimplexBasis(gq_int deg, gq_int dim1, gq_int dim2);
   SimplexSimplexBasis(const SimplexSimplexBasis&) = default;
   std::unique_ptr<Basis> clone() const override;

   SimplexSimplexBasis& operator=(const SimplexSimplexBasis&) = delete;
   ~SimplexSimplexBasis() override = default;

   const Array1D& orthog_basis(const Array1D& x) final;
   const Array2D& orthog_der(const Array1D& x) final;
   const Array1D& orthog_integrals() final;

   const Array1D& monomial_basis(const Array1D& x) final;
   const Array2D& monomial_der(const Array1D& x) final;
   const Array1D& monomial_integrals() final;

private:
   std::array<gq_int, 2> dims_;
   BasisTable power_table_first_;
   BasisTable power_table_second_;
   BasisTable index_table_;

   Array1D poly_buff1_;
   Array1D poly_buff2_;
   Array1D fdiff_[4];

   void orthog_basis_internal(const Array1D& x, Array1D& buff);
   void orthog_basis_polytopic_one_internal(const Array1D& x, Array1D& buff);
   void orthog_basis_polytopic_two_internal(const Array1D& x, Array1D& buff);

   template <typename BasisType>
   friend void internal::fourth_order_difference(BasisType& basis, const Array1D& x);
};


gq_int SimplexSimplexBasisSize(gq_int deg, gq_int dim);


class PyramidBasis3D : public PolytopeBasis {
public:
   explicit PyramidBasis3D(gq_int deg);
   PyramidBasis3D(const PyramidBasis3D&) = default;
   std::unique_ptr<Basis> clone() const override;

   PyramidBasis3D& operator=(const PyramidBasis3D&) = delete;
   ~PyramidBasis3D() override = default;

   const Array1D& orthog_basis(const Array1D& x) final;
   const Array2D& orthog_der(const Array1D& x) final;
   const Array1D& orthog_integrals() final;

   const Array1D& monomial_basis(const Array1D& x) final;
   const Array2D& monomial_der(const Array1D& x) final;
   const Array1D& monomial_integrals() final;

private:
   BasisTable index_table_;
   Array1D fdiff_[4];

   void orthog_basis_internal(const Array1D& x, Array1D& buff);

   template <typename BasisType>
   friend void internal::fourth_order_difference(BasisType& basis, const Array1D& x);
};


gq_int Pyramid3DBasisSize(gq_int deg);


// forward declare and include quadrature header in basis.cpp
class QuadSimplex;
class QuadOmega2D;


class OmegaBasis2D : public Basis {
public:
   explicit OmegaBasis2D(Omega2D omega, gq_int deg);
   OmegaBasis2D(const OmegaBasis2D&);
   std::unique_ptr<Basis> clone() const override;

   OmegaBasis2D& operator=(const PyramidBasis3D&) = delete;
   ~OmegaBasis2D() override;

   const Array1D& orthog_basis(const Array1D& x) final;
   const Array2D& orthog_der(const Array1D& x) final;
   const Array1D& orthog_integrals() final;

   const Array1D& monomial_basis(const Array1D& x) final;
   const Array2D& monomial_der(const Array1D& x) final;
   const Array1D& monomial_integrals() final;

   void test_orthogonal(gq_int max_deg, bool verbose = false);

private:
   // types are incomplete so pointers are used
   QuadSimplex* quad_simplex_;
   QuadOmega2D* quad_omega_;
   const Omega2D omega_;
   BasisTable index_table_;
   Array1D fdiff_[4];
   StdVector<Matrix2D> A_, B_, C_, D_;
   ArrayInt1D phi_pos_;
   double phi0_;

   void recurrence_coeffs(gq_int n);
   void orthog_basis_internal(const Array1D& x, Array1D& buff);
   void orthog_basis_internal_p(gq_int p, const Array1D& x, Array1D& buff);

   template <typename BasisType>
   friend void internal::fourth_order_difference(BasisType& basis, const Array1D& x);
};


gq_int Omega2DBasisSize(gq_int deg);


}  // namespace gquad

#endif

