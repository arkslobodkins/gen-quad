#ifndef VECTOR_H
#define VECTOR_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
   int *id;
   int len;
} VectorInt;

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

#define StaticVectorInit(vecName, __len)  \
   Vector vecName;                        \
   double __ ##vecName[__len];            \
   vecName.len = __len;                   \
   vecName.id = __##vecName;

VectorInt VectorInt_init(int n);
void VectorInt_realloc(int n, VectorInt *V);
void VectorInt_Assign(VectorInt v1, VectorInt v2);
void VectorInt_free(VectorInt V);

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
