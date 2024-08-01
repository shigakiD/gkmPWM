/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rot90.h
 *
 * Code generation for function 'rot90'
 *
 */

#ifndef ROT90_H
#define ROT90_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_rot90(const double A[16], double B[16]);

void c_rot90(const double A[16], double B[16]);

void d_rot90(const emxArray_real_T *A, emxArray_real_T *B);

void rot90(const emxArray_real_T *A, emxArray_real_T *B);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (rot90.h) */
