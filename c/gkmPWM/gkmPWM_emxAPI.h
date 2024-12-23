/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWM_emxAPI.h
 *
 * Code generation for function 'gkmPWM_emxAPI'
 *
 */

#ifndef GKMPWM_EMXAPI_H
#define GKMPWM_EMXAPI_H

/* Include files */
#include "gkmPWM_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern emxArray_char_T *emxCreateND_char_T(int numDimensions, const int *size);

extern emxArray_char_T *emxCreateWrapperND_char_T(char *data, int numDimensions,
                                                  const int *size);

extern emxArray_char_T *emxCreateWrapper_char_T(char *data, int rows, int cols);

extern emxArray_char_T *emxCreate_char_T(int rows, int cols);

extern void emxDestroyArray_char_T(emxArray_char_T *emxArray);

extern void emxInitArray_char_T(emxArray_char_T **pEmxArray, int numDimensions);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (gkmPWM_emxAPI.h) */
