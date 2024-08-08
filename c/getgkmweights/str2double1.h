/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * str2double1.h
 *
 * Code generation for function 'str2double1'
 *
 */

#ifndef STR2DOUBLE1_H
#define STR2DOUBLE1_H

/* Include files */
#include "getgkmweights_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void readfloat(emxArray_char_T *s1, int *idx, const emxArray_char_T *s, int *k,
               int n, bool *isimag, bool *b_finite, double *nfv,
               bool *foundsign, bool *success);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (str2double1.h) */
