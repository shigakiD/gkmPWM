/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sortedInsertion.h
 *
 * Code generation for function 'sortedInsertion'
 *
 */

#ifndef SORTEDINSERTION_H
#define SORTEDINSERTION_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void sortedInsertion(double x, int ix, emxArray_real_T *b, int *nb, int blen,
                     emxArray_int32_T *idx);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (sortedInsertion.h) */
