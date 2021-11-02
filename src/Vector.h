#ifndef VECTOR_H
#define VECTOR_H

#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
   bool *id;
   int len;
} Vector_bool;

typedef struct
{
   int *id;
   int len;
} Vector_int;

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

Vector_bool Vector_bool_init(int n);
void Vector_bool_realloc(int n, Vector_bool *V);
void Vector_bool_free(Vector_bool V);

Vector_int Vector_int_init(int n);
void Vector_int_realloc(int n, Vector_int *V);
void Vector_int_free(Vector_int V);

Vector Vector_init(int n);
void Vector_realloc(int n, Vector *V);
void Vector_Assign(const Vector v1, Vector v2);
void Vector_free(Vector V);

void Vector_Print(const Vector V);
double VDot(const Vector a, const Vector b);
void VectorAddScale(double c1, const Vector V1, double c2, const Vector V2, Vector V3);
VMin VectorMin(const Vector v);

void VRemoveElement(int index, Vector *z);
void VIntRemoveElement(int index, Vector_int *z);

double V_ScaledTwoNorm(const Vector z);
double V_TwoNorm(const Vector z);
double V_InfNorm(const Vector z);
bool V_CheckInfAndNan(const Vector z);


int IntPower(int x, int power);

#ifdef __cplusplus
}
#endif

#endif
