// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "basis.hpp"
#include "domain.hpp"
#include "gen_quad.hpp"
#include "mesh.hpp"

namespace gquad {


class QuadArray {
private:
   Array1D array_;
   gq_int deg_;
   gq_int dim_;
   gq_int num_nodes_;

public:
   explicit QuadArray(gq_int deg, gq_int dim, gq_int num_nodes);
   QuadArray(const QuadArray& q);
   QuadArray(QuadArray&& q);
   QuadArray& operator=(const QuadArray& q) &;
   QuadArray& operator=(QuadArray&& q) &;

   inline double& w(gq_int i) &;
   inline const double& w(gq_int i) const&;
   inline auto node(gq_int i) &;
   inline const auto node(gq_int i) const&;

   // node_v returns a vector that has semantics of a matrix rather than array
   inline auto node_v(gq_int i) &;
   inline const auto node_v(gq_int i) const&;

   auto weights() & {
      return array_(seq(0, num_nodes_ - 1));
   }

   const auto weights() const& {
      return array_(seq(0, num_nodes_ - 1));
   }

   auto nodes() & {
      return array_(seq(num_nodes_, last));
   }

   const auto nodes() const& {
      return array_(seq(num_nodes_, last));
   }

   double* w_ptr() & {
      return &array_[0];
   }

   const double* w_ptr() const& {
      return &array_[0];
   }

   double* x_ptr() & {
      return &array_[num_nodes()];
   }

   const double* x_ptr() const& {
      return &array_[num_nodes()];
   }

   double& operator[](gq_int i) & {
      return array_[i];
   }

   const double& operator[](gq_int i) const& {
      return array_[i];
   }

   Array1D& array() & {
      return array_;
   }

   const Array1D& array() const& {
      return array_;
   }

   gq_int deg() const {
      return deg_;
   }

   gq_int dim() const {
      return dim_;
   }

   gq_int num_nodes() const {
      return num_nodes_;
   }

   gq_int size() const {
      return array_.size();
   }

   void resize(gq_int new_num_nodes);
   void resize_and_assign(const QuadArray& q);
   void reinit(gq_int new_dim, gq_int new_num_nodes);
};


std::ostream& operator<<(std::ostream& os, const QuadArray& q);


inline double& QuadArray::w(gq_int i) & {
   GEN_QUAD_ASSERT_DEBUG(i > -1 && i < num_nodes_);
   return array_[i];
}


inline const double& QuadArray::w(gq_int i) const& {
   GEN_QUAD_ASSERT_DEBUG(i > -1 && i < num_nodes_);
   return array_[i];
}


inline auto QuadArray::node(gq_int i) & {
   GEN_QUAD_ASSERT_DEBUG(i > -1 && i < num_nodes_);
   gq_int offset = num_nodes_ + i * dim_;
   return array_(seqN(offset, dim_));
}


inline const auto QuadArray::node(gq_int i) const& {
   GEN_QUAD_ASSERT_DEBUG(i > -1 && i < num_nodes_);
   gq_int offset = num_nodes_ + i * dim_;
   return array_(seqN(offset, dim_));
}


inline auto QuadArray::node_v(gq_int i) & {
   return node(i).matrix();
}


inline const auto QuadArray::node_v(gq_int i) const& {
   return node(i).matrix();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract class that is used to handle quadrature domains that are derived from it.
class QuadDomain : protected QuadArray {
public:
   virtual std::unique_ptr<QuadDomain> clone_quad_domain() const = 0;
   virtual ~QuadDomain() = default;

   virtual const Domain& get_domain() const = 0;
   virtual std::unique_ptr<Basis> create_basis() const;

   virtual double relative_exponential_residual() const = 0;
   virtual std::string quad_file_name() const = 0;
   virtual std::string quad_dims_file_name() const = 0;

   using QuadArray::node;
   using QuadArray::node_v;
   using QuadArray::nodes;
   using QuadArray::w;
   using QuadArray::w_ptr;
   using QuadArray::weights;
   using QuadArray::x_ptr;
   using QuadArray::operator[];
   using QuadArray::array;

   using QuadArray::deg;
   using QuadArray::dim;
   using QuadArray::num_nodes;
   using QuadArray::size;

   using QuadArray::resize;

   virtual gq_int num_nodes_optimal() const = 0;
   double efficiency() const;

   QuadDomain& resize_and_assign(const QuadDomain& q);

protected:
   using QA = QuadArray;

   explicit QuadDomain(gq_int deg, gq_int dim, gq_int num_nodes);
   QuadDomain(const QuadDomain&) = default;
   QuadDomain(QuadDomain&&) = default;
   QuadDomain& operator=(const QuadDomain&) = default;
   QuadDomain& operator=(QuadDomain&&) = default;

   friend std::ostream& operator<<(std::ostream& os, const QuadDomain& q);
};


template <typename F>
double quad_approx(const QuadDomain& q, F f) {
   Array1D eval(q.num_nodes());
   for(gq_int i = 0; i < q.num_nodes(); ++i) {
      eval[i] = f(q.node(i));
   }

   return stable_dot_prod(eval, q.weights());
}


double quad_exp_approx(const QuadDomain& q);


gq_int mult_nodes(const QuadDomain& q1, const QuadDomain& q2);


bool in_domain_elem(const QuadDomain& q, gq_int i);
bool in_domain(const QuadDomain& q);
bool in_constraint(const QuadDomain& q);
bool pos_weights(const QuadDomain& q);


std::ostream& operator<<(std::ostream& os, const QuadDomain& q);


class QuadPolytope : public QuadDomain {
public:
   virtual std::unique_ptr<QuadPolytope> clone_quad_polytope() const = 0;
   virtual const Polytope& get_polytope() const = 0;
   virtual ~QuadPolytope() = default;

protected:
   explicit QuadPolytope(gq_int deg, gq_int dim, gq_int num_nodes);
   QuadPolytope(const QuadPolytope&) = default;
   QuadPolytope(QuadPolytope&&) = default;
   QuadPolytope& operator=(const QuadPolytope&) = default;
   QuadPolytope& operator=(QuadPolytope&&) = default;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class QuadConvexPolytope : public QuadPolytope {
public:
   virtual ~QuadConvexPolytope() = default;
   virtual const ConvexPolytope& get_convex_polytope() const = 0;

protected:
   explicit QuadConvexPolytope(gq_int deg, gq_int dim, gq_int num_nodes);
   QuadConvexPolytope(const QuadConvexPolytope&) = default;
   QuadConvexPolytope(QuadConvexPolytope&&) = default;
   QuadConvexPolytope& operator=(const QuadConvexPolytope&) = default;
   QuadConvexPolytope& operator=(QuadConvexPolytope&&) = default;
};


double dist_from_boundary_min(const QuadPolytope& q);
double dist_from_constr_min(const QuadPolytope& q);


////////////////////////////////////////////////////////////////////////////////////////////////////
class QuadInterval : public QuadConvexPolytope {
public:
   explicit QuadInterval(gq_int deg, gq_int num_nodes = 0);

   std::unique_ptr<QuadDomain> clone_quad_domain() const override;
   std::unique_ptr<QuadPolytope> clone_quad_polytope() const override;

   const Domain& get_domain() const override;
   const Polytope& get_polytope() const override;
   const ConvexPolytope& get_convex_polytope() const override;

   double relative_exponential_residual() const override;
   std::string quad_file_name() const override;
   std::string quad_dims_file_name() const override;

   gq_int num_nodes_optimal() const override;

private:
   Interval interval;
};


QuadInterval quadrature_gauss_legendre(gq_int deg);
QuadInterval quadrature_gauss_jacobi(gq_int deg, gq_int alpha, gq_int beta);


////////////////////////////////////////////////////////////////////////////////////////////////////
class QuadCube : public QuadConvexPolytope {
public:
   explicit QuadCube(gq_int deg, gq_int dim, gq_int num_nodes = 0);

   std::unique_ptr<QuadDomain> clone_quad_domain() const override;
   std::unique_ptr<QuadPolytope> clone_quad_polytope() const override;

   void reinit(gq_int dim, gq_int num_nodes);
   void reinit_copy(const QuadCube& qc);

   const Domain& get_domain() const override;
   const Polytope& get_polytope() const override;
   const ConvexPolytope& get_convex_polytope() const override;
   std::unique_ptr<Basis> create_basis() const override;

   double relative_exponential_residual() const override;
   std::string quad_file_name() const override;
   std::string quad_dims_file_name() const override;

   gq_int num_nodes_optimal() const override;

private:
   Cube cube_;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class QuadSimplex : public QuadConvexPolytope {
public:
   QuadSimplex(gq_int deg, gq_int dim, gq_int num_nodes = 0);

   std::unique_ptr<QuadDomain> clone_quad_domain() const override;
   std::unique_ptr<QuadPolytope> clone_quad_polytope() const override;

   void reinit(gq_int dim, gq_int num_nodes);
   void reinit_copy(const QuadSimplex& qs);

   const Domain& get_domain() const override;
   const Polytope& get_polytope() const override;
   const ConvexPolytope& get_convex_polytope() const override;
   std::unique_ptr<Basis> create_basis() const override;

   double relative_exponential_residual() const override;
   std::string quad_file_name() const override;
   std::string quad_dims_file_name() const override;

   gq_int num_nodes_optimal() const override;

private:
   Simplex simplex_;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class QuadCubeSimplex : public QuadConvexPolytope {
public:
   QuadCubeSimplex(gq_int deg, gq_int dim1, gq_int dim2, gq_int num_nodes = 0);

   std::unique_ptr<QuadDomain> clone_quad_domain() const override;
   std::unique_ptr<QuadPolytope> clone_quad_polytope() const override;

   void reinit(gq_int dim1, gq_int dim2, gq_int num_nodes);
   void reinit_copy(const QuadCubeSimplex& qcs);

   const Domain& get_domain() const override;
   const Polytope& get_polytope() const override;
   const ConvexPolytope& get_convex_polytope() const override;
   std::unique_ptr<Basis> create_basis() const override;

   double relative_exponential_residual() const override;
   std::string quad_file_name() const override;
   std::string quad_dims_file_name() const override;

   gq_int dim1() const {
      return cubesimplex_.dim1();
   }

   gq_int dim2() const {
      return cubesimplex_.dim2();
   }

   gq_int num_nodes_optimal() const override;

private:
   CubeSimplex cubesimplex_;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class QuadSimplexSimplex : public QuadConvexPolytope {
public:
   QuadSimplexSimplex(gq_int deg, gq_int dim1, gq_int dim2, gq_int num_nodes = 0);

   std::unique_ptr<QuadDomain> clone_quad_domain() const override;
   std::unique_ptr<QuadPolytope> clone_quad_polytope() const override;

   void reinit(gq_int dim1, gq_int dim2, gq_int num_nodes);
   void reinit_copy(const QuadSimplexSimplex& qss);

   const Domain& get_domain() const override;
   const Polytope& get_polytope() const override;
   const ConvexPolytope& get_convex_polytope() const override;
   std::unique_ptr<Basis> create_basis() const override;

   double relative_exponential_residual() const override;
   std::string quad_file_name() const override;
   std::string quad_dims_file_name() const override;

   gq_int dim1() const {
      return simplexsimplex_.dim1();
   }

   gq_int dim2() const {
      return simplexsimplex_.dim2();
   }

   gq_int num_nodes_optimal() const override;

private:
   SimplexSimplex simplexsimplex_;
};


// remaining classes do not have reinit and reinit_copy functions, since their
// dimension never changes
////////////////////////////////////////////////////////////////////////////////////////////////////
class QuadPyramid3D : public QuadConvexPolytope {
public:
   explicit QuadPyramid3D(gq_int deg, gq_int num_nodes = 0);

   std::unique_ptr<QuadDomain> clone_quad_domain() const override;
   std::unique_ptr<QuadPolytope> clone_quad_polytope() const override;

   const Domain& get_domain() const override;
   const Polytope& get_polytope() const override;
   const ConvexPolytope& get_convex_polytope() const override;
   std::unique_ptr<Basis> create_basis() const override;

   double relative_exponential_residual() const override;
   std::string quad_file_name() const override;
   std::string quad_dims_file_name() const override;

   gq_int num_nodes_optimal() const override;

private:
   Pyramid3D pyramid3D_;
};


class QuadOmega2D : public QuadPolytope {
public:
   explicit QuadOmega2D(Omega2D omega, gq_int deg, gq_int num_nodes = 0);

   std::unique_ptr<QuadDomain> clone_quad_domain() const override;
   std::unique_ptr<QuadPolytope> clone_quad_polytope() const override;

   const Domain& get_domain() const override;
   const Polytope& get_polytope() const override;
   const Omega2D& get_domain_self() const;

   std::unique_ptr<Basis> create_basis() const override;
   std::unique_ptr<OmegaBasis2D> create_basis_self() const;

   double relative_exponential_residual() const override;
   std::string quad_file_name() const override;
   std::string quad_dims_file_name() const override;

   gq_int num_nodes_optimal() const override;

private:
   const Omega2D omega_;
};


QuadOmega2D CreateOmegaComposite(Omega2D omega, const QuadSimplex& qs);


}  // namespace gquad

#endif

