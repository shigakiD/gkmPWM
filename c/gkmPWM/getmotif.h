/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getmotif.h
 *
 * Code generation for function 'getmotif'
 *
 */

#ifndef GETMOTIF_H
#define GETMOTIF_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void getmotif(const emxArray_char_T *filename, const emxArray_real_T *m,
              emxArray_cell_wrap_4 *mat, emxArray_cell_wrap_0 *names);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (getmotif.h) */
