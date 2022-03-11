#ifndef COMPUTE_DOMAIN_TYPE_H
#define COMPUTE_DOMAIN_TYPE_H

#include "GENERAL_QUADRATURE.h"

#ifdef __cplusplus
extern "C" {
#endif

// All routines print history and quadrature rules to results directory.

// returns Gauss-Legendre quadrature on [0, 1]
void ComputeInterval(int degree);


// runs Node Elimination but starts with expensive
// generalized tensor product of G-L quadrature.
// Mostly for comparison/testing purposes.
void ComputeCubeFull(int deg, int dim);


// runs Node Elimination but starts with expensive
// generalized tensor product of G-J quadrature.
// Mostly for comparison/testing purposes.
void ComputeSimplexFull(int deg, int dim);


// Runs Node Elimination algorithm using recursive initial guess.
void ComputeCube(int degree, int dim);


// Runs Node Elimination algorithm using recursive initial guess.
void ComputeSimplex(int degree, int dim);


// Runs Node Elimination algorithm using recursive initial guess for
// dim-1 cube and dim-2 simplex, but takes a tensor product of (cube x simplex).
void ComputeCubeSimplexTensor(int deg, int dim1, int dim2);


// Runs Node Elimination algorithm using recursive initial guess.
void ComputeCubeSimplex(int degree, int dim1, int dim2);


// Runs Node Elimination algorithm using recursive initial guess for
// dim-1 simplex and dim-2 simplex, but takes a tensor product of(simplex1 x simplex2).
void ComputeSimplexSimplex(int degree, int dim1, int dim2);

#ifdef __cplusplus
}
#endif

#endif
