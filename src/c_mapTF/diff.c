/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * diff.c
 *
 * Code generation for function 'diff'
 *
 */

/* Include files */
#include "diff.h"
#include "mapTF2_ls_emxutil.h"
#include "mapTF2_ls_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void diff(const emxArray_boolean_T *x, emxArray_real_T *y)
{
  double *y_data;
  int dimSize;
  int i;
  int m;
  int tmp1;
  int work_data;
  const bool *x_data;
  x_data = x->data;
  dimSize = x->size[1];
  if (x->size[1] == 0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    work_data = x->size[1] - 1;
    if (work_data > 1) {
      work_data = 1;
    }
    if (work_data < 1) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1] - 1;
      emxEnsureCapacity_real_T(y, i);
      y_data = y->data;
      if (x->size[1] - 1 != 0) {
        work_data = x_data[0];
        for (m = 2; m <= dimSize; m++) {
          tmp1 = x_data[m - 1];
          i = tmp1;
          tmp1 -= work_data;
          work_data = i;
          y_data[m - 2] = tmp1;
        }
      }
    }
  }
}

/* End of code generation (diff.c) */
