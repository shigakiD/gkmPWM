/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * combineVectorElements.c
 *
 * Code generation for function 'combineVectorElements'
 *
 */

/* Include files */
#include "combineVectorElements.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void combineVectorElements(const emxArray_boolean_T *x, emxArray_int32_T *y)
{
  unsigned int sz[2];
  int j;
  int k;
  int vlen;
  int vstride;
  int xoffset;
  int *y_data;
  const bool *x_data;
  x_data = x->data;
  vlen = x->size[1];
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    for (j = 0; j < 2; j++) {
      sz[j] = (unsigned int)x->size[j];
    }
    j = y->size[0];
    y->size[0] = (int)sz[0];
    emxEnsureCapacity_int32_T(y, j);
    y_data = y->data;
    xoffset = (int)sz[0];
    for (j = 0; j < xoffset; j++) {
      y_data[j] = 0;
    }
  } else {
    vstride = x->size[0];
    j = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_int32_T(y, j);
    y_data = y->data;
    for (j = 0; j < vstride; j++) {
      y_data[j] = x_data[j];
    }
    for (k = 2; k <= vlen; k++) {
      xoffset = (k - 1) * vstride;
      for (j = 0; j < vstride; j++) {
        y_data[j] += x_data[xoffset + j];
      }
    }
  }
}

/* End of code generation (combineVectorElements.c) */
