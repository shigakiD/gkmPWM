/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eml_setop.h
 *
 * Code generation for function 'eml_setop'
 *
 */

#ifndef EML_SETOP_H
#define EML_SETOP_H

/* Include files */
#include "mapTF_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
                emxArray_real_T *c, emxArray_int32_T *ia, int *ib_size);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (eml_setop.h) */