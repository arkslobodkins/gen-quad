// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2023

#include <limits>
#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/domain.hpp"
#include "../include/mesh.hpp"

namespace gquad {

Vertex::Vertex(double x1, double x2) : x{x1, x2} {
}

Edge::Edge(const Vertex& v1, const Vertex& v2) : vert{v1, v2} {
}

Triangle::Triangle(const Vertex& v1, const Vertex& v2, const Vertex& v3)
    : vert{v1, v2, v3},
      edg{Edge{v1, v2}, Edge{v2, v3}, Edge{v3, v1}} {
}

static std::istream& operator>>(std::istream& is, std::pair<Vertex, bool>& vb) {
   is >> vb.first.x[0];
   is >> vb.first.x[1];
   is >> vb.second;
   return is;
}

static std::istream& operator>>(std::istream& is, std::tuple<gq_int, gq_int, bool>& ib) {
   is >> std::get<0>(ib);
   is >> std::get<1>(ib);
   is >> std::get<2>(ib);
   return is;
}

static std::istream& operator>>(std::istream& is, std::tuple<gq_int, gq_int, gq_int>& iii) {
   is >> std::get<0>(iii);
   is >> std::get<1>(iii);
   is >> std::get<2>(iii);
   return is;
}

template <typename T>
static bool equal_members(const std::vector<T>& v1, const std::vector<T>& v2) {
   GEN_QUAD_ASSERT_DEBUG(v1.size() == v2.size());
   for(std::size_t i = 0; i < v1.size(); ++i) {
      if(!(v1[i] == v2[i])) {
         return false;
      }
   }
   return true;
}

std::pair<Point<2>, Point<2>> get_bounding_box_2D(const std::vector<Vertex>& v) {
   GEN_QUAD_ASSERT_DEBUG(!v.empty());
   Point<2> ll{v[0].x[0], v[0].x[1]};
   Point<2> uu = ll;

   for(std::size_t i = 1; i < v.size(); ++i) {
      if(v[i].x[0] < ll[0]) {
         ll[0] = v[i].x[0];
      }
      if(v[i].x[1] < ll[1]) {
         ll[1] = v[i].x[1];
      }

      if(v[i].x[0] > uu[0]) {
         uu[0] = v[i].x[0];
      }
      if(v[i].x[1] > uu[1]) {
         uu[1] = v[i].x[1];
      }
   }
   return {ll, uu};
}

std::ostream& operator<<(std::ostream& os, const Vertex& vertex) {
   char s1{}, s2{};
   if(vertex.x[0] >= 0) {
      s1 = ' ';
   }
   if(vertex.x[1] >= 0) {
      s2 = ' ';
   }

   std::cout << std::fixed;
   std::cout << std::setprecision(4);
   std::cout << "(" << s1 << vertex.x[0] << "  " << s2 << vertex.x[1] << ")";
   return os;
}

std::ostream& operator<<(std::ostream& os, const Edge& edge) {
   std::cout << std::fixed;
   std::cout << std::setprecision(4);
   std::cout << "{ " << edge.vert[0] << ",  " << edge.vert[1] << " }";
   return os;
}

std::ostream& operator<<(std::ostream& os, const Triangle& triangle) {
   std::cout << std::fixed;
   std::cout << std::setprecision(4);
   std::cout << "{ " << triangle.vert[0] << ",  " << triangle.vert[1] << ",  " << triangle.vert[2] << " }";
   return os;
}

Edge::EdgeParams Edge::getparams() const {
   EdgeParams ep;
   const auto& v0 = vert[0];
   const auto& v1 = vert[1];
   ep.tangent[0] = v1.x[0] - v0.x[0];
   ep.tangent[1] = v1.x[1] - v0.x[1];
   ep.length = std::sqrt(dot_prod(ep.tangent, ep.tangent));
   ep.tangent[0] /= ep.length;
   ep.tangent[1] /= ep.length;
   ep.normal[0] = ep.tangent[1];
   ep.normal[1] = -ep.tangent[0];
   return {ep};
}

bool Vertex::operator==(const Vertex& v) const {
   return (x[0] == v.x[0]) && (x[1] == v.x[1]);
}

bool Edge::operator==(const Edge& e) const {
   return (vert[0] == e.vert[0]) && (vert[1] == e.vert[1]);
}

bool Triangle::operator==(const Triangle& t) const {
   return (vert[0] == t.vert[0]) && (vert[1] == t.vert[1]) && (vert[2] == t.vert[2]);
}

double Triangle::jacobian() const {
   StaticArray1D<2> a1, a2;
   for(gq_int d = 0; d < 2; ++d) {
      a1[d] = vert[1].x[d] - vert[0].x[d];
      a2[d] = vert[2].x[d] - vert[1].x[d];
   }

   return std::fabs(a1[0] * a2[1] - a1[1] * a2[0]);
}

Omega2D::Omega2D(std::vector<Vertex> vert_, std::vector<Edge> edg_, std::vector<Triangle> triang_,
                 std::vector<bool> vert_on_boundary_, std::vector<bool> edg_on_boundary_, std::string name_)
    : Polytope{},
      vert{std::move(vert_)},
      edg{std::move(edg_)},
      triang{std::move(triang_)},
      vert_on_boundary{std::move(vert_on_boundary_)},
      edg_on_boundary{std::move(edg_on_boundary_)},
      lower_left{get_bounding_box_2D(vert).first},
      upper_right{get_bounding_box_2D(vert).second},
      name{std::move(name_)} {
   GEN_QUAD_ASSERT_DEBUG(!vert.empty());
   GEN_QUAD_ASSERT_DEBUG(!edg.empty());
   GEN_QUAD_ASSERT_DEBUG(!triang.empty());

   GEN_QUAD_ASSERT_DEBUG(vert.size() == vert_on_boundary.size());
   GEN_QUAD_ASSERT_DEBUG(edg.size() == edg_on_boundary.size());
}

Omega2D& Omega2D::operator=(const Omega2D& omega) {
   GEN_QUAD_ASSERT_DEBUG(vert.size() == omega.vert.size());
   GEN_QUAD_ASSERT_DEBUG(edg.size() == omega.vert.size());
   GEN_QUAD_ASSERT_DEBUG(triang.size() == omega.vert.size());

   GEN_QUAD_ASSERT_DEBUG(equal_members(vert, omega.vert));
   GEN_QUAD_ASSERT_DEBUG(equal_members(edg, omega.edg));
   GEN_QUAD_ASSERT_DEBUG(equal_members(triang, omega.triang));

   GEN_QUAD_ASSERT_DEBUG(equal_members(vert_on_boundary, omega.vert_on_boundary));
   GEN_QUAD_ASSERT_DEBUG(equal_members(edg_on_boundary, omega.edg_on_boundary));

   return *this;
}

bool Omega2D::in_domain(const Array1D& x) const {
   Simplex s(2);
   Array1D xmap(2);

   for(std::size_t nt = 0; nt < this->triang.size(); ++nt) {
      map_to_unit(this->triang[nt], x, xmap);
      if(s.in_domain(xmap)) {
         return true;
      }
   }

   return false;
}

std::tuple<bool, gq_int, Point<2>> Omega2D::in_domain_info(const Array1D& x) const {
   Simplex s(2);
   Array1D xmap(2);

   for(std::size_t nt = 0; nt < this->triang.size(); ++nt) {
      map_to_unit(this->triang[nt], x, xmap);
      if(s.in_domain(xmap)) {
         return {true, gq_int(nt), {xmap[0], xmap[1]}};
      }
   }

   Array1D z(2);
   z[0] = std::numeric_limits<double>::infinity();
   z[1] = std::numeric_limits<double>::infinity();
   return {false, -1, z};
}

double Omega2D::dist_from_boundary(const Array1D& x) const {
   double dist = std::numeric_limits<double>::infinity();
   for(std::size_t k = 0; k < edg.size(); ++k) {
      const auto params = edg[k].getparams();

      if(edg_on_boundary[k]) {
         const auto& v0 = edg[k].vert[0];
         StaticArray1D<2> r{x[0] - v0.x[0], x[1] - v0.x[1]};
         double alf = r[0] * params.tangent[0] + r[1] * params.tangent[1];

         double rr;
         if(alf < -1.e-10 * params.length) {
            rr = std::sqrt(r[0] * r[0] + r[1] * r[1]);
         } else if(alf > 1.0000000001 * params.length) {
            const auto& v1 = edg[k].vert[1];
            r[0] = x[0] - v1.x[0];
            r[1] = x[1] - v1.x[1];
            rr = std::sqrt(r[0] * r[0] + r[1] * r[1]);
         } else {
            double bet = r[0] * params.normal[0] + r[1] * params.normal[1];
            rr = std::fabs(bet);
         }
         if(rr < dist) {
            dist = rr;
         }
      }
   }
   return -dist;
}

double Omega2D::area() const {
   double a = 0.;
   for(std::size_t n = 0; n < triang.size(); ++n) {
      a += triang[n].jacobian();
   }
   return a *= 0.5;
}

Omega2D ReadOmega(const std::string& file_vert, const std::string& file_edg, const std::string& file_triang,
                  const std::string domain_name) {
   std::vector<Vertex> vert;
   std::vector<Edge> edg;
   std::vector<Triangle> triang;
   std::vector<bool> vert_bound;
   std::vector<bool> edg_bound;

   std::ifstream ifs_vert{file_vert};
   std::ifstream ifs_edg{file_edg};
   std::ifstream ifs_triang{file_triang};
   if(!ifs_vert) {
      GEN_QUAD_ASSERT_ALWAYS_MSG(false, "file input failed for vertexes");
   }
   if(!ifs_edg) {
      GEN_QUAD_ASSERT_ALWAYS_MSG(false, "file input failed for edges");
   }
   if(!ifs_triang) {
      GEN_QUAD_ASSERT_ALWAYS_MSG(false, "file input failed for triangles");
   }

   std::pair<Vertex, bool> vb;
   while(ifs_vert >> vb) {
      vert.push_back(vb.first);
      vert_bound.push_back(vb.second);
   }
   if(!ifs_vert.eof()) {
      GEN_QUAD_ASSERT_ALWAYS_MSG(false, "vertexes did not reach end of file successfully");
   }

   std::tuple<gq_int, gq_int, bool> ib;
   while(ifs_edg >> ib) {
      Edge e{vert[std::get<0>(ib)], vert[std::get<1>(ib)]};
      edg.push_back(e);
      edg_bound.push_back(std::get<2>(ib));
   }
   if(!ifs_edg.eof()) {
      GEN_QUAD_ASSERT_ALWAYS_MSG(false, "edges did not reach end of file successfully");
   }

   std::tuple<gq_int, gq_int, gq_int> iii;
   while(ifs_triang >> iii) {
      Triangle t{vert[std::get<0>(iii)], vert[std::get<1>(iii)], vert[std::get<2>(iii)]};
      triang.push_back(t);
   }
   if(!ifs_triang.eof()) {
      GEN_QUAD_ASSERT_ALWAYS_MSG(false, "triangles did not reach end of file successfully");
   }

   return Omega2D{std::move(vert),
                  std::move(edg),
                  std::move(triang),
                  std::move(vert_bound),
                  std::move(edg_bound),
                  domain_name};
}

Omega2D CreateSquare() {
   std::vector<bool> vb{1, 1, 1, 1};
   std::vector<Vertex> v{Vertex{0., 0.}, Vertex{1., 0.}, Vertex{1., 1.}, Vertex{0., 1.}};

   std::vector<bool> eb{1, 1, 1, 1, 0};
   std::vector<Edge> e{
       Edge{v[0], v[1]}, Edge{v[1], v[2]}, Edge{v[2], v[3]}, Edge{v[3], v[0]}, Edge{v[0], v[2]}};

   std::vector<Triangle> t{Triangle{v[0], v[1], v[2]}, Triangle{v[2], v[3], v[0]}};

   return Omega2D{std::move(v), std::move(e), std::move(t), std::move(vb), std::move(eb), "square"};
}

Omega2D CreatePentagon() {
   std::vector<bool> vb{1, 1, 1, 1, 1};
   std::vector<Vertex> v{Vertex{0., 1.},
                         Vertex{0.951056516295154, 0.309016994374947},
                         Vertex{0.587785252292473, -0.809016994374947},
                         Vertex{-0.587785252292473, -0.809016994374948},
                         Vertex{-0.951056516295154, 0.309016994374947}};

   std::vector<bool> eb{0, 0, 1, 1, 1, 1, 1};
   std::vector<Edge> e{Edge{v[4], v[1]},
                       Edge{v[4], v[2]},
                       Edge{v[0], v[1]},
                       Edge{v[1], v[2]},
                       Edge{v[2], v[3]},
                       Edge{v[3], v[4]},
                       Edge{v[4], v[0]}};

   std::vector<Triangle> t{
       Triangle{v[4], v[0], v[1]}, Triangle{v[4], v[1], v[2]}, Triangle{v[4], v[2], v[3]}};

   return Omega2D{std::move(v), std::move(e), std::move(t), std::move(vb), std::move(eb), "pentagon"};
}

Omega2D CreateHexagon() {
   std::vector<bool> vb{1, 1, 1, 1, 1, 1};
   double sqt32 = std::sqrt(3.) / 2.;
   std::vector<Vertex> v{Vertex{1., 0.},
                         Vertex{0.5, sqt32},
                         Vertex{-0.5, sqt32},
                         Vertex{-1., 0.},
                         Vertex{-0.5, -sqt32},
                         Vertex{0.5, -sqt32}};

   std::vector<bool> eb{0, 0, 0, 1, 1, 1, 1, 1, 1};
   std::vector<Edge> e{Edge{v[2], v[4]},
                       Edge{v[4], v[1]},
                       Edge{v[1], v[5]},
                       Edge{v[1], v[0]},
                       Edge{v[2], v[1]},
                       Edge{v[3], v[2]},
                       Edge{v[4], v[3]},
                       Edge{v[5], v[4]},
                       Edge{v[0], v[5]}};

   std::vector<Triangle> t{Triangle{v[2], v[3], v[4]},
                           Triangle{v[2], v[4], v[1]},
                           Triangle{v[1], v[4], v[5]},
                           Triangle{v[5], v[0], v[1]}};

   return Omega2D{std::move(v), std::move(e), std::move(t), std::move(vb), std::move(eb), "hexagon"};
}

Omega2D CreateIrregTrapezoid() {
   std::vector<bool> vb{1, 1, 1, 1};
   std::vector<Vertex> v{Vertex{-1., 0.}, Vertex{-0.5, 1.}, Vertex{0.75, 0.75}, Vertex{1., 0.}};

   std::vector<bool> eb{0, 1, 1, 1, 1};
   std::vector<Edge> e{
       Edge{v[0], v[2]}, Edge{v[0], v[1]}, Edge{v[1], v[2]}, Edge{v[2], v[3]}, Edge{v[3], v[0]}};

   std::vector<Triangle> t{Triangle{v[0], v[1], v[2]}, Triangle{v[2], v[3], v[0]}};
   return Omega2D{std::move(v), std::move(e), std::move(t), std::move(vb), std::move(eb), "irreg_trapezoid"};
}

Omega2D CreateIrreg5() {
   std::vector<bool> vb{1, 1, 1, 1, 1};
   std::vector<Vertex> v{Vertex{0., 0.},
                         Vertex{0., 1.},
                         Vertex{0.951056516295154, 0.309016994374947},
                         Vertex{0.587785252292473, -0.809016994374947},
                         Vertex{-0.587785252292473, -0.809016994374947}};

   std::vector<bool> eb{0, 0, 1, 1, 1, 1};
   std::vector<Edge> e{Edge{v[0], v[2]},
                       Edge{v[0], v[3]},
                       Edge{v[0], v[4]},
                       Edge{v[2], v[1]},
                       Edge{v[3], v[2]},
                       Edge{v[4], v[3]}};

   std::vector<Triangle> t{
       Triangle{v[0], v[1], v[2]}, Triangle{v[0], v[2], v[3]}, Triangle{v[0], v[3], v[4]}};
   return Omega2D{std::move(v), std::move(e), std::move(t), std::move(vb), std::move(eb), "irreg5"};
}

std::ostream& operator<<(std::ostream& os, const Omega2D& omega) {
   std::cout << "vertexes:" << "\n\n";
   for(std::size_t i = 0; i < omega.vert.size(); ++i) {
      std::cout << omega.vert[i] << ",  " << omega.vert_on_boundary[i] << '\n';
   }

   std::cout << "\nedges:" << "\n\n";
   for(std::size_t i = 0; i < omega.edg.size(); ++i) {
      std::cout << omega.edg[i] << ",  " << omega.edg_on_boundary[i] << '\n';
   }

   std::cout << "\ntriangles:" << "\n\n";
   for(std::size_t i = 0; i < omega.triang.size(); ++i) {
      std::cout << omega.triang[i] << '\n';
   }

   std::cout << "\nlower left bound:" << "\n";
   std::cout << omega.lower_left[0] << "  " << omega.lower_left[1];

   std::cout << "\n\nupper right bound:" << "\n";
   std::cout << omega.upper_right[0] << "  " << omega.upper_right[1];

   return os;
}

}  // namespace gquad

#endif
