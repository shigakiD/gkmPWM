/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * nullAssignment.h
 *
 * Code generation for function 'nullAssignment'
 *
 */

#ifndef NULLASSIGNMENT_H
#define NULLASSIGNMENT_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_nullAssignment(double x_data[], int x_size[2], const int idx[2]);

void c_nullAssignment(double x_data[], int x_size[2], const int idx[2]);

void nullAssignment(emxArray_real_T *x, int idx);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (nullAssignment.h) */
