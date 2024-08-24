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
   Vertex() = default;
   Vertex(const Vertex&) = default;
   Vertex& operator=(const Vertex&) = default;

   explicit Vertex(double x1, double x2);
   bool operator==(const Vertex& v) const;

   Point<2> x;
};


std::ostream& operator<<(std::ostream& os, const Vertex& v);


struct Edge {
   struct EdgeParams {
      StaticArray1D<2> tangent;
      StaticArray1D<2> normal;
      double length;
   };

   Edge() = default;
   Edge(const Edge&) = default;
   Edge& operator=(const Edge&) = default;

   explicit Edge(const Vertex& v1, const Vertex& v2);
   bool operator==(const Edge& e) const;
   EdgeParams getparams() const;

   std::array<Vertex, 2> vert;
};


std::ostream& operator<<(std::ostream& os, const Edge& edge);


struct Triangle {
   Triangle() = default;
   Triangle(const Triangle&) = default;
   Triangle& operator=(const Triangle&) = default;

   explicit Triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3);
   bool operator==(const Triangle& t) const;
   double jacobian() const;

   std::array<Vertex, 3> vert;
   std::array<Edge, 3> edg;
};


std::ostream& operator<<(std::ostream& os, const Triangle& triangle);


struct Omega2D : public Polytope {
   explicit Omega2D(StdVector<Vertex> vertexes, StdVector<Edge> edges, StdVector<Triangle> triangles,
                    StdVector<bool> vertexes_on_boundary, StdVector<bool> edges_on_boundary,
                    std::string name);

   Omega2D(const Omega2D&) = default;
   Omega2D& operator=(const Omega2D&);  // asserts that all members are equal

   bool in_domain(const Array1D& x) const override;
   std::tuple<bool, gq_int, Point<2>> in_domain_info(const Array1D& x) const;
   double dist_from_boundary(const Array1D& x) const override;

   gq_int dim() const {
      return 2;
   }

   std::string domain_name() const override {
      return s;
   }

   double area() const;

   const StdVector<Vertex> vert;
   const StdVector<Edge> edg;
   const StdVector<Triangle> triang;
   const StdVector<bool> vert_on_boundary;
   const StdVector<bool> edg_on_boundary;
   const Point<2> lower_left;
   const Point<2> upper_right;

private:
   std::string s;
};


Omega2D ReadOmega(const std::string& file_vert, const std::string& file_edg, const std::string& file_triang);
Omega2D CreateSquare();
Omega2D CreatePentagon();
Omega2D CreateHexagon();
Omega2D CreateIrregTrapezoid();
Omega2D CreateIrreg5();


std::ostream& operator<<(std::ostream& os, const Omega2D& omega);
std::pair<Point<2>, Point<2>> get_bounding_box_2D(const StdVector<Vertex>& v);


template <typename EigT1, typename EigT2>
void map_to_unit(const Triangle& t, const EigT1& x, EigT2&& xmap) {
   StaticArray1D<2> R;
   StaticArray2D<2, 2> A;

   for(gq_int d = 0; d < 2; ++d) {
      R[d] = x[d] - t.vert[0].x[d];
      A(d, 0) = t.vert[1].x[d] - t.vert[0].x[d];
      A(d, 1) = t.vert[2].x[d] - t.vert[1].x[d];
   }
   double det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
   xmap[0] = (A(1, 1) * R[0] - A(0, 1) * R[1]) / det;
   xmap[1] = (A(0, 0) * R[1] - A(1, 0) * R[0]) / det;
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

