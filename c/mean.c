/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mean.c
 *
 * Code generation for function 'mean'
 *
 */

/* Include files */
#include "mean.h"
#include "blockedSummation.h"
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void mean(const emxArray_real_T *x, double y_data[], int y_size[2])
{
  emxArray_real_T *y;
  double *b_y_data;
  int i;
  int loop_ub;
  emxInit_real_T(&y, 2);
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real_T(y, i);
    b_y_data = y->data;
    loop_ub = x->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_y_data[i] = 0.0;
    }
  } else {
    colMajorFlatIter(x, x->size[0], y);
    b_y_data = y->data;
  }
  y_size[0] = 1;
  y_size[1] = y->size[1];
  loop_ub = y->size[1];
  for (i = 0; i < loop_ub; i++) {
    y_data[i] = b_y_data[i] / (double)x->size[0];
  }
  emxFree_real_T(&y);
}

/* End of code generation (mean.c) */
