/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repelem.c
 *
 * Code generation for function 'repelem'
 *
 */

/* Include files */
#include "repelem.h"
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void repelem(const double x[16], double varargin_1, emxArray_real_T *y)
{
  double *y_data;
  int idx;
  int j;
  int k;
  int ni;
  idx = y->size[0];
  y->size[0] = (int)varargin_1 << 4;
  emxEnsureCapacity_real_T(y, idx);
  y_data = y->data;
  idx = -1;
  ni = (int)varargin_1;
  for (k = 0; k < 16; k++) {
    for (j = 0; j < ni; j++) {
      y_data[(idx + j) + 1] = x[k];
    }
    idx += (int)varargin_1;
  }
}

/* End of code generation (repelem.c) */
