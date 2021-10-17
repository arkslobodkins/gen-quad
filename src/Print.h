#ifndef PRINT_H
#define PRINT_H

#include "GENERAL_QUADRATURE.h"


#ifdef __cplusplus
extern "C" {

#define PRINT
#define PRINT_ERR

}

#else

/*********************************************************
\* MACROS and functions for printing/debugging purposes \*
*********************************************************/

#define PRINT(x, str) _Generic( (x),    \
int:      print_int,                    \
float:    print_float,                  \
double:   print_double,                 \
bool:     print_bool,                   \
char *:   print_string                  \
)(x, str)

#define PRINT_ERR(x, str, line_num, file_name) _Generic( (x),        \
int:      print_int_err,                                             \
float:    print_float_err,                                           \
double:   print_double_err,                                          \
bool:     print_bool_err,                                            \
char *:   print_string_err                                           \
)(x, str, line_num, file_name)


#endif


#ifdef __cplusplus
extern "C" {
#endif

void PrintNodes(const_quadrature *q, const char *name);
void PrintNodesAndWeights(const_quadrature *q, const char *name);
void PrintNodeAndWeight(int id, const_quadrature *q, const char *name);
void PrintNodeInfo(int iters, double error_norm, const _DomainFuncs dom_funcs, const_quadrature *q, const char *name);
void PrintElimInfo(int dim, int num_nodes, int opt);

void print_int(int x, const char *name);
void print_bool(bool x, const char *name);
void print_float(float x, const char *name);
void print_double(double x, const char *name);
void print_string(const char *x, int format);

void print_int_err(int x, const char *name, int line_num, const char *file_name);
void print_bool_err(bool x, const char *name, int line_num, const char *file_name);
void print_float_err(float x, const char *name, int line_num, const char *file_name);
void print_double_err(double x, const char *name, int line_num, const char *file_name);
void print_string_err(const char *x, int format, int line_num, const char *file_name);

#ifdef __cplusplus
}
#endif

#endif
