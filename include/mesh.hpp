#pragma once

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "domain.hpp"
#include "gen_quad.hpp"

namespace gquad {

template <gq_int dimension>
using Point = StaticArray1D<dimension>;

struct Vertex {
   Vertex() : x{} {
   }

   explicit Vertex(double x1, double x2) : x{x1, x2} {
   }

   bool operator==(const Vertex& v) const;

   Point<2> x;
};

std::ostream& operator<<(std::ostream& os, const Vertex& vertex);

struct Edge {
   Edge() : vert{} {
   }

   explicit Edge(const Vertex& v1, const Vertex& v2) : vert{v1, v2} {
   }

   bool operator==(const Edge& e) const;

   std::array<Vertex, 2> vert;

   struct EdgeParams {
      StaticArray1D<2> tangent;
      StaticArray1D<2> normal;
      double length;
   };

   EdgeParams getparams() const;
};

std::ostream& operator<<(std::ostream& os, const Edge& edge);

struct Triangle {
   Triangle() : vert{}, edg{} {
   }

   explicit Triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3)
       : vert{v1, v2, v3},
         edg{Edge{v1, v2}, Edge{v2, v3}, Edge{v3, v1}} {
   }

   bool operator==(const Triangle& t) const;
   double jacobian() const;

   std::array<Vertex, 3> vert;
   std::array<Edge, 3> edg;
};

std::ostream& operator<<(std::ostream& os, const Triangle& triangle);

struct Omega2D : public Polytope {
   explicit Omega2D(std::vector<Vertex> vertexes_, std::vector<Edge> edges_, std::vector<Triangle> triangles_,
                    std::vector<bool> vertexes_on_boundary_, std::vector<bool> edges_on_boundary_,
                    std::string name_);
   Omega2D(const Omega2D&) = default;
   Omega2D& operator=(const Omega2D&);  // asserts that all members are equal

   const std::vector<Vertex> vert;
   const std::vector<Edge> edg;
   const std::vector<Triangle> triang;
   const std::vector<bool> vert_on_boundary;
   const std::vector<bool> edg_on_boundary;
   const Point<2> lower_left;
   const Point<2> upper_right;

   bool in_domain(const Array1D& x) const override;
   std::tuple<bool, gq_int, Point<2>> in_domain_info(const Array1D& x) const;
   double dist_from_boundary(const Array1D& x) const override;

   gq_int dim() const {
      return 2;
   }

   std::string domain_name() const override {
      return name;
   }

   double area() const;

private:
   std::string name;
};

Omega2D ReadOmega(const std::string& file_vert, const std::string& file_edg, const std::string& file_triang);
Omega2D CreateSquare();
Omega2D CreatePentagon();
Omega2D CreateHexagon();
Omega2D CreateIrregTrapezoid();
Omega2D CreateIrreg5();

std::ostream& operator<<(std::ostream& os, const Omega2D& omega);
std::pair<Point<2>, Point<2>> get_bounding_box_2D(const std::vector<Vertex>& v);

template <typename EigT1, typename EigT2>
void map_to_unit(const Triangle& t, const EigT1& x, EigT2&& xmap);

template <typename EigT1, typename EigT2>
void map_from_unit(const Triangle& t, const EigT1& x, EigT2&& xmap);

template <typename EigT1, typename EigT2>
void map_to_unit(const Triangle& t, const EigT1& x, EigT2&& xmap) {
   StaticArray1D<2> r;
   StaticArray2D<2, 2> A;

   for(gq_int d = 0; d < 2; ++d) {
      r[d] = x[d] - t.vert[0].x[d];
      A(d, 0) = t.vert[1].x[d] - t.vert[0].x[d];
      A(d, 1) = t.vert[2].x[d] - t.vert[1].x[d];
   }
   double det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
   xmap[0] = (A(1, 1) * r[0] - A(0, 1) * r[1]) / det;
   xmap[1] = (A(0, 0) * r[1] - A(1, 0) * r[0]) / det;
}

template <typename EigT1, typename EigT2>
void map_from_unit(const Triangle& t, const EigT1& x, EigT2&& xmap) {
   StaticArray2D<2, 2> A;

   for(gq_int d = 0; d < 2; ++d) {
      A(0, d) = t.vert[1].x[d] - t.vert[0].x[d];
      A(1, d) = t.vert[2].x[d] - t.vert[1].x[d];
   }
   xmap[0] = t.vert[0].x[0] + A(0, 0) * x[0] + A(1, 0) * x[1];
   xmap[1] = t.vert[0].x[1] + A(0, 1) * x[0] + A(1, 1) * x[1];
}

}  // namespace gquad

#endif

