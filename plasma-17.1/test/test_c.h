/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_z.h, normal z -> c, Tue Feb  8 19:15:03 2022
 *
 **/
#ifndef TEST_C_H
#define TEST_C_H

#include "test.h"

//==============================================================================
// test routines
//==============================================================================
void test_scamax(param_value_t param[], char *info);
void test_cgbsv(param_value_t param[], char *info);
void test_cgbtrf(param_value_t param[], char *info);
void test_cgeadd(param_value_t param[], char *info);
void test_cgelqf(param_value_t param[], char *info);
void test_cgelqs(param_value_t param[], char *info);
void test_cgels(param_value_t param[], char *info);
void test_cgemm(param_value_t param[], char *info);
void test_cgeqrf(param_value_t param[], char *info);
void test_cgeqrs(param_value_t param[], char *info);
void test_cgesv(param_value_t param[], char *info);
void test_cgetrf(param_value_t param[], char *info);
void test_cgetri(param_value_t param[], char *info);
void test_cgetri_aux(param_value_t param[], char *info);
void test_cgetrs(param_value_t param[], char *info);
void test_chemm(param_value_t param[], char *info);
void test_cher2k(param_value_t param[], char *info);
void test_cherk(param_value_t param[], char *info);
void test_clacpy(param_value_t param[], char *info);
void test_clag2z(param_value_t param[], char *info);
void test_clange(param_value_t param[], char *info);
void test_clanhe(param_value_t param[], char *info);
void test_clansy(param_value_t param[], char *info);
void test_clantr(param_value_t param[], char *info);
void test_clascl(param_value_t param[], char *info);
void test_claset(param_value_t param[], char *info);
void test_claswp(param_value_t param[], char *info);
void test_clauum(param_value_t param[], char *info);
void test_cpbsv(param_value_t param[], char *info);
void test_cpbtrf(param_value_t param[], char *info);
void test_cposv(param_value_t param[], char *info);
void test_cpotrf(param_value_t param[], char *info);
void test_cpotri(param_value_t param[], char *info);
void test_cpotrs(param_value_t param[], char *info);
void test_csymm(param_value_t param[], char *info);
void test_csyr2k(param_value_t param[], char *info);
void test_csyrk(param_value_t param[], char *info);
void test_ctradd(param_value_t param[], char *info);
void test_ctrmm(param_value_t param[], char *info);
void test_ctrsm(param_value_t param[], char *info);
void test_ctrtri(param_value_t param[], char *info);
void test_cunmlq(param_value_t param[], char *info);
void test_cunmqr(param_value_t param[], char *info);

#endif // TEST_C_H
