#ifndef VECTOR_H
#define VECTOR_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
   double *id;
   int len;
} Vector;

typedef struct
{
   int min_index;
   double min_value;
} VMin;

Vector Vector_init(int n);
void Vector_realloc(int n, Vector *V);
void Vector_Assign(const Vector v1, Vector v2);
void Vector_free(Vector V);

void VPrint(const Vector V);
double VDot(const Vector a, const Vector b);
void VectorAddScale(double c1, const Vector V1, double c2, const Vector V2, Vector V3);
VMin VectorMin(const Vector v);
void VectorRemoveElement(int index, Vector *z);
double V_ScaledTwoNorm(const Vector z);
double V_TwoNorm(const Vector z);
double V_InfNorm(const Vector z);
bool V_CheckNan(const Vector z);
bool V_CheckInf(const Vector z);

int IntPower(int x, int power);

#ifdef __cplusplus
}
#endif

#endif
