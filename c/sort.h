/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sort.h
 *
 * Code generation for function 'sort'
 *
 */

#ifndef SORT_H
#define SORT_H

/* Include files */
#include "gkmPWMlasso3_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_sort(emxArray_real_T *x, emxArray_int32_T *idx);

void c_sort(double x_data[], const int x_size[2], int idx_data[],
            int idx_size[2]);

void sort(emxArray_real_T *x, emxArray_int32_T *idx);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (sort.h) */
