/* Arkadijs Slobodkins
 * SMU Mathematics
 * Copyright 2022
 */


#ifndef VECTOR_H
#define VECTOR_H

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>

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

typedef struct
{
   int max_index;
   double max_value;
} VMax;

#define StaticVectorInit(size, vecName)           \
   Vector vecName;                                \
   double __##vecName[size];                      \
   memset(__##vecName, 0, (size)*sizeof(double)); \
   vecName.len = size;                            \
   vecName.id = __##vecName;

VectorInt VectorInt_init(int n);
void VectorInt_realloc(int n, VectorInt *V);
void VectorInt_Assign(VectorInt v1, VectorInt v2);
void VectorInt_free(VectorInt V);

Vector Vector_init(int n);
Vector Vector_uninitialized(int n);
void Vector_realloc(int n, Vector *V);
void Vector_Assign(Vector v1, Vector v2);
void Vector_free(Vector V);

#ifdef _OPENMP
void AllocVectorOmpData(Vector *v);
void FreeVectorOmpData(Vector v);
#endif

void VPrint(Vector V);
double VDot(Vector a, Vector b);
void VSetToZero(Vector z);
void VSetToOne(Vector v);
void VScale(double c, Vector V);
void Vector_daxpy(double a, Vector x, Vector y);
void VAdd(Vector V1, Vector V2, Vector V3);
void VAddScale(double c1, Vector V1, double c2, Vector V2, Vector V3);
void VRemoveElement(int index, Vector *z);
VMin VectorMin(Vector v);
VMax VectorMax(Vector v);

double V_ScaledTwoNorm(Vector z);
double V_TwoNorm(Vector z);
double V_InfNorm(Vector z);
bool V_CheckNan(Vector z);
bool V_CheckInf(Vector z);
bool V_CheckInfOrNan(Vector z);
bool V_IsUninitialized(Vector v);

void double_daxpy(int len, double a, double *x, double *y);
void double_dcopy(int n, double *x, double *y);


static inline void set_elem(Vector v, int i, double x)
{
   assert(i >= 0 && i <= v.len-1);
   v.id[i] = x;
}

static inline double get_elem(Vector v, int i)
{
   assert(i >= 0 && i <= v.len-1);
   return v.id[i];
}

static inline double DDot(int len, double *a, double *b)
{
   double d = 0.0;
   for(int i = 0; i < len; ++i)
      d += a[i] * b[i];
   return d;
}

static inline Vector DToVec(int len, double *x)
{
   Vector v;
   v.ompId = NULL;
   v.len = len;
   v.id = x;
   return v;
}


#ifdef __cplusplus
}
#endif

#endif
