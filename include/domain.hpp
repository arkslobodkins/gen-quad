// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "gen_quad.hpp"

namespace gquad {


////////////////////////////////////////////////////////////////////////////////////////////////////
class Domain {
public:
   virtual ~Domain() = default;

   virtual bool in_domain(const Array1D& x) const = 0;
   virtual std::string domain_name() const = 0;

protected:
   explicit Domain() = default;
   Domain(const Domain&) = default;
   Domain& operator=(const Domain&) = delete;
};


class Polytope : public Domain {
public:
   virtual ~Polytope() = default;

   virtual double dist_from_boundary(const Array1D& x) const = 0;

protected:
   explicit Polytope() = default;
   Polytope(const Polytope&) = default;
   Polytope& operator=(const Polytope&) = delete;
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class IdealPolytope : public Polytope {
public:
   using Constraints = std::pair<Matrix2D, Vector1D>;  // Linear inequality constraints
   virtual ~IdealPolytope() = default;

   bool in_domain(const Array1D& x) const override;
   double dist_from_boundary(const Array1D& x) const override;

   const Matrix2D& get_constraint_matrix() const {
      return A_;
   }

   const Vector1D& get_constraint_vector() const {
      return b_;
   }

protected:
   Matrix2D A_;
   Vector1D b_;

   explicit IdealPolytope(Constraints C);
   IdealPolytope(const IdealPolytope&) = default;
   IdealPolytope& operator=(const IdealPolytope&) = delete;
};


class Interval : public IdealPolytope {
public:
   explicit Interval();
   Interval(const Interval&) = default;

   Interval& operator=(const Interval&) {
      return *this;
   }

   std::string domain_name() const override {
      return "interval";
   }

   static gq_int dim() {
      return 1;
   }

private:
   Constraints ConstructConstraints();
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class Cube : public IdealPolytope {
public:
   explicit Cube(gq_int dim);
   Cube(const Cube&) = default;
   Cube& operator=(const Cube&);

   std::string domain_name() const override {
      return "cube";
   }

   gq_int dim() const {
      return dim_;
   }

   void reinit(gq_int dim);

private:
   gq_int dim_;
   Constraints ConstructConstraints(gq_int dim);
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class Simplex : public IdealPolytope {
public:
   explicit Simplex(gq_int dim);
   Simplex(const Simplex&) = default;
   Simplex& operator=(const Simplex&);

   // implements in_domain differently to completely avoid roundoff errors
   bool in_domain(const Array1D& x) const override;

   std::string domain_name() const override {
      return "simplex";
   }

   gq_int dim() const {
      return dim_;
   }

   void reinit(gq_int dim);

private:
   gq_int dim_;
   Constraints ConstructConstraints(gq_int dim);
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class CubeSimplex : public IdealPolytope {
public:
   explicit CubeSimplex(std::array<gq_int, 2> dims);
   CubeSimplex(const CubeSimplex&) = default;
   CubeSimplex& operator=(const CubeSimplex&);

   // implements in_domain differently to completely avoid roundoff errors
   bool in_domain(const Array1D& x) const override;

   std::string domain_name() const override {
      return "cubesimplex";
   }

   gq_int dim() const {
      return dims_[0] + dims_[1];
   }

   gq_int dim1() const {
      return dims_[0];
   }

   gq_int dim2() const {
      return dims_[1];
   }

   void reinit(gq_int dim1, gq_int dim2);

private:
   std::array<gq_int, 2> dims_;
   Constraints ConstructConstraints(gq_int dim1, gq_int dim2);
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class SimplexSimplex : public IdealPolytope {
public:
   explicit SimplexSimplex(std::array<gq_int, 2> dims);
   SimplexSimplex(const SimplexSimplex&) = default;
   SimplexSimplex& operator=(const SimplexSimplex&);

   // implements in_domain differently to completely avoid roundoff errors
   bool in_domain(const Array1D& x) const override;

   std::string domain_name() const override {
      return "simplexsimplex";
   }

   gq_int dim() const {
      return dims_[0] + dims_[1];
   }

   gq_int dim1() const {
      return dims_[0];
   }

   gq_int dim2() const {
      return dims_[1];
   }

   void reinit(gq_int dim1, gq_int dim2);

private:
   std::array<gq_int, 2> dims_;
   Constraints ConstructConstraints(gq_int dim1, gq_int dim2);
};


////////////////////////////////////////////////////////////////////////////////////////////////////
class Pyramid3D : public IdealPolytope {
public:
   explicit Pyramid3D();
   Pyramid3D(const Pyramid3D&) = default;

   Pyramid3D& operator=(const Pyramid3D&) {
      return *this;
   };

   // implements in_domain differently to completely avoid roundoff errors
   bool in_domain(const Array1D& x) const override;

   std::string domain_name() const override {
      return "pyramid";
   }

   static gq_int dim() {
      return 3;
   }

private:
   Constraints ConstructConstraints();
};


}  // namespace gquad

#endif

