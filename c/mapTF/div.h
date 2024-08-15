/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * div.h
 *
 * Code generation for function 'div'
 *
 */

#ifndef DIV_H
#define DIV_H

/* Include files */
#include "mapTF_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_binary_expand_op(emxArray_real_T *ff, const emxArray_real_T *kmat,
                        const emxArray_cell_wrap_2 *seqindmat, int b_I, int i,
                        const emxArray_real_T *minnorm,
                        const emxArray_int32_T *p1, const emxArray_real_T *f);

void binary_expand_op(emxArray_real_T *ff, const emxArray_real_T *kmat,
                      const emxArray_cell_wrap_2 *seqindmat, int b_I, int i,
                      const emxArray_real_T *minnorm, const emxArray_real_T *p2,
                      const emxArray_real_T *f);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (div.h) */
