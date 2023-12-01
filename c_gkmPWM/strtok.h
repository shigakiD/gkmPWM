/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * strtok.h
 *
 * Code generation for function 'strtok'
 *
 */

#ifndef STRTOK_H
#define STRTOK_H

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
void b_strtok(const emxArray_char_T *x, emxArray_char_T *token,
              emxArray_char_T *remain);

void c_strtok(const emxArray_char_T *x, emxArray_char_T *token,
              emxArray_char_T *remain);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (strtok.h) */
