/* Copyright 1985 Pierre Asselin.
 *
 * This package may be used and distributed freely, but may not be sold.
 * The present notice must be retained in all copies.
 * If incorporated in a commercial product,
 *  1) the package's origin and availability must be acknowledged
 *     prominently in the commercial product's documentation;
 *  2) the source code and documentation of the package must be
 *     made available to purchasers on request and at no extra cost.
 **/

/* Translation of module Jacobi (export).
 **/

// October 15 2021 :
// Function definitions and declarations have been changed from old-style
// to new-style by Arkadijs Slobodkins for C++ compatibility.

int Jacobi(int n, double alpha, double beta, double abscis[], double weight[]);
int Radau_Jacobi(int n, double alpha, double beta, double abscis[], double weight[], double *leftw);
int Lobatto_Jacobi(int n, double alpha, double beta, double abscis[], double weight[], double *leftw, double *rightw);

/* ARGUMENT LISTS:
 * Jacobi(n, alpha, beta, abscis, weight)
 * int n;
 * double alpha, beta;
 * double abscis[], weight[];
 *
 * Radau_Jacobi(n, alpha, beta, abscis, weight, leftw)
 * int n;
 * double alpha, beta;
 * double abscis[], weight[], *leftw;
 *
 * Lobatto_Jacobi(n, alpha, beta, abscis, weight, leftw, rightw)
 * int n;
 * double alpha, beta;
 * double abscis[], weight[], *leftw, *rightw;
 **/
