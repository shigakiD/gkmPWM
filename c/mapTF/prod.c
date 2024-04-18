/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * prod.c
 *
 * Code generation for function 'prod'
 *
 */

/* Include files */
#include "prod.h"
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void prod(const emxArray_real_T *x, emxArray_real_T *y)
{
  const double *x_data;
  double *y_data;
  int j;
  int k;
  int vlen;
  int vstride;
  int xoffset;
  x_data = x->data;
  vlen = x->size[1];
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    j = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real_T(y, j);
    y_data = y->data;
    xoffset = x->size[0];
    for (j = 0; j < xoffset; j++) {
      y_data[j] = 1.0;
    }
  } else {
    vstride = x->size[0];
    j = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real_T(y, j);
    y_data = y->data;
    for (j = 0; j < vstride; j++) {
      y_data[j] = x_data[j];
    }
    for (k = 2; k <= vlen; k++) {
      xoffset = (k - 1) * vstride;
      for (j = 0; j < vstride; j++) {
        y_data[j] *= x_data[xoffset + j];
      }
    }
  }
}

/* End of code generation (prod.c) */
