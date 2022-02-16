#ifndef VECTOR_H
#define VECTOR_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

// Level-1 BLAS
void dcopy_(int *N, double *DX, int *INCX, double *DY, int *INCY);
void dscal_(int *n, double *c, double *x, int *incx);
void daxpy_(int *n, double *c, double *x, int *incx, double *y, int *incy);
double dnrm2_(int *n, double *x, int *incx);
int idamax_(int *n, double *x, int *incx);

typedef struct
{
   double *id;
   double **ompId;
   int len;
} Vector;

typedef struct
{
   int min_index;
   double min_value;
} VMin;

Vector Vector_init(int n);
void Vector_realloc(int n, Vector *V);
void Vector_Assign(Vector v1, Vector v2);
void Vector_free(Vector V);

#ifdef _OPENMP
void AllocVectorOmpData(Vector *v);
void FreeVectorOmpData(Vector v);
#endif

void VPrint(Vector V);
double VDot(Vector a, Vector b);
void VSetToOne(Vector v);
void VectorScale(double c, Vector V);
void Vector_daxpy(double a, Vector x, Vector y);
void double_daxpy(int len, double a, double *x, double *y);
void VectorAddScale(double c1, Vector V1, double c2, Vector V2, Vector V3);
VMin VectorMin(Vector v);
void VectorRemoveElement(int index, Vector *z);
double V_ScaledTwoNorm(Vector z);
double V_TwoNorm(Vector z);
double V_InfNorm(Vector z);
bool V_CheckNan(Vector z);
bool V_CheckInf(Vector z);
void double_dcopy(int n, double *x, double *y);

#ifdef __cplusplus
}
#endif

#endif
