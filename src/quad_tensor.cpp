// Arkadijs Slobodkins
// Johannes Tausch
// SMU Mathematics
// Copyright 2022

#if __cplusplus < 201402L
#error requires c++14 or higher
#else

#include "../include/quad_tensor.hpp"

namespace gquad {


static void assert_tensor_common(const QuadDomain& q1, const QuadDomain& q2, const QuadDomain& q3) {
   GEN_QUAD_ASSERT_DEBUG(q1.deg() == q2.deg() && q2.deg() == q3.deg());
   GEN_QUAD_ASSERT_DEBUG(q1.dim() + q2.dim() == q3.dim());
   GEN_QUAD_ASSERT_DEBUG(q1.num_nodes() * q2.num_nodes() == q3.num_nodes());
}


QuadCube GaussTensorCube(gq_int deg, gq_int dim) {
   GEN_QUAD_ASSERT_DEBUG(deg > 0);
   GEN_QUAD_ASSERT_DEBUG(dim > 1);

   QuadInterval qi = quadrature_gauss_legendre(deg);
   QuadCube qc(deg, 2, mult_nodes(qi, qi));
   AddLineCube(qi, qi, qc);

   QuadCube qc_next = qc;

   for(gq_int d = 3; d <= dim; ++d) {
      qc_next.reinit(d, mult_nodes(qi, qc));
      AddLineCube(qi, qc, qc_next);
      qc.reinit_copy(qc_next);
   }

   return qc_next;
}


QuadSimplex GaussTensorSimplex(gq_int deg, gq_int dim) {
   GEN_QUAD_ASSERT_DEBUG(deg > 0);
   GEN_QUAD_ASSERT_DEBUG(dim > 1);

   QuadInterval qgl = quadrature_gauss_legendre(deg);
   QuadInterval qgj = quadrature_gauss_jacobi(deg, 0, 1);
   QuadSimplex qs(deg, 2, mult_nodes(qgj, qgl));
   AddLineSimplex(qgj, qgl, qs);

   QuadSimplex qs_next = qs;

   for(gq_int d = 3; d <= dim; ++d) {
      QuadInterval qgj_next = quadrature_gauss_jacobi(deg, 0, d - 1);
      qs_next.reinit(d, mult_nodes(qgj_next, qs));
      AddLineSimplex(qgj_next, qs, qs_next);

      qs.reinit_copy(qs_next);
   }

   return qs_next;
}


QuadCubeSimplex GaussTensorCubeSimplex(gq_int deg, gq_int dim1, gq_int dim2) {
   GEN_QUAD_ASSERT_DEBUG(deg > 0);
   GEN_QUAD_ASSERT_DEBUG(dim1 > 0);
   GEN_QUAD_ASSERT_DEBUG(dim2 > 1);

   QuadSimplex qs = GaussTensorSimplex(deg, dim2);

   if(dim1 == 1) {
      QuadInterval qi = quadrature_gauss_legendre(deg);
      QuadCubeSimplex qcs(deg, dim1, dim2, mult_nodes(qi, qs));
      AddLineCube(qi, qs, qcs);
      return qcs;
   } else {
      QuadCube qc = GaussTensorCube(deg, dim1);
      QuadCubeSimplex qcs(deg, dim1, dim2, mult_nodes(qc, qs));
      MixedTensor(qc, qs, qcs);
      return qcs;
   }
}


QuadSimplexSimplex GaussTensorSimplexSimplex(gq_int deg, gq_int dim1, gq_int dim2) {
   GEN_QUAD_ASSERT_DEBUG(deg > 0);
   GEN_QUAD_ASSERT_DEBUG(dim1 > 1);
   GEN_QUAD_ASSERT_DEBUG(dim2 > 1);

   auto qs1 = GaussTensorSimplex(deg, dim1);
   auto qs2 = GaussTensorSimplex(deg, dim2);
   QuadSimplexSimplex qss(deg, dim1, dim2, mult_nodes(qs1, qs2));
   MixedTensor(qs1, qs2, qss);

   return qss;
}


void NodesTensor2D(const QuadInterval& q1, const QuadInterval& q2, QuadDomain& q3) {
   assert_tensor_common(q1, q2, q3);

   gq_int n1 = q1.num_nodes();
   gq_int n2 = q2.num_nodes();
   const double* x1 = q1.x_ptr();
   const double* x2 = q2.x_ptr();
   double* x3 = q3.x_ptr();

   for(gq_int i = 0; i < n1; ++i) {
      for(gq_int j = 0; j < n2; ++j) {
         x3[2 * (n2 * i + j)] = x1[i];
         x3[2 * (n2 * i + j) + 1] = x2[j];
      }
   }
}


void WeightsTensor2D(const QuadInterval& q1, const QuadInterval& q2, QuadDomain& q3) {
   assert_tensor_common(q1, q2, q3);

   gq_int n1 = q1.num_nodes();
   gq_int n2 = q2.num_nodes();
   for(gq_int i = 0; i < n1; ++i) {
      for(gq_int j = 0; j < n2; ++j) {
         q3.w(i * n2 + j) = q1.w(i) * q2.w(j);
      }
   }
}


void Tensor2D(const QuadInterval& q1, const QuadInterval& q2, QuadDomain& q3) {
   NodesTensor2D(q1, q2, q3);
   WeightsTensor2D(q1, q2, q3);
}


static void WeightsTensor(const QuadDomain& q1, const QuadDomain& q2, QuadDomain& q3) {
   assert_tensor_common(q1, q2, q3);

   gq_int n1 = q1.num_nodes();
   gq_int n2 = q2.num_nodes();
   for(gq_int i = 0; i < n1; ++i) {
      for(gq_int j = 0; j < n2; ++j) {
         q3.w(i * n2 + j) = q1.w(i) * q2.w(j);
      }
   }
}


void MixedTensor(const QuadDomain& q1, const QuadDomain& q2, QuadDomain& q3) {
   assert_tensor_common(q1, q2, q3);

   const gq_int dim1 = q1.dim(), dim2 = q2.dim();
   const gq_int n1 = q1.num_nodes(), n2 = q2.num_nodes();

   for(gq_int i = 0; i < n1; ++i) {
      for(gq_int j = 0; j < n2; ++j) {
         for(gq_int d1 = 0; d1 < dim1; ++d1) {
            q3.node(i * n2 + j)(d1) = q1.node(i)(d1);
         }
         for(gq_int d2 = 0; d2 < dim2; ++d2) {
            q3.node(i * n2 + j)(dim1 + d2) = q2.node(j)(d2);
         }
      }
   }
   WeightsTensor(q1, q2, q3);
}


void AddLineCube(const QuadInterval& q1, const QuadDomain& q2, QuadDomain& q3) {
   assert_tensor_common(q1, q2, q3);

   gq_int n1 = q1.num_nodes();
   gq_int n2 = q2.num_nodes();
   gq_int dim2 = q2.dim();

   for(gq_int i = 0; i < n1; ++i) {
      for(gq_int j = 0; j < n2; ++j) {
         q3.node(i * n2 + j)(0) = q1.node(i)(0);
      }
   }

   for(gq_int i = 0; i < n1; ++i) {
      for(gq_int j = 0; j < n2; ++j) {
         for(gq_int d = 0; d < dim2; ++d) {
            q3.node(i * n2 + j)(d + 1) = q2.node(j)(d);
         }
      }
   }

   WeightsTensor(q1, q2, q3);
}


void AddLineSimplex(const QuadInterval& q1, const QuadDomain& q2, QuadDomain& q3) {
   assert_tensor_common(q1, q2, q3);

   AddLineCube(q1, q2, q3);
   // apply one level of Duffy Transformation
   for(gq_int i = 0; i < q3.num_nodes(); ++i) {
      for(gq_int d = 1; d < q3.dim(); ++d) {
         q3.node(i)(d) *= q3.node(i)(0);
      }
   }
}


void AddLinePyramid(const QuadInterval& q1, const QuadDomain& q2, QuadDomain& q3) {
   assert_tensor_common(q1, q2, q3);
   assert(q3.dim() == 3);

   AddLineCube(q1, q2, q3);
   for(gq_int i = 0; i < q3.num_nodes(); ++i) {
      q3.node(i)(1) *= q3.node(i)(0);  // yz
      q3.node(i)(2) *= q3.node(i)(0);  // xz
   }
}


}  // namespace gquad

#endif

