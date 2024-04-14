/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * div.c
 *
 * Code generation for function 'div'
 *
 */

/* Include files */
#include "div.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include <string.h>

/* Function Definitions */
void v_binary_expand_op(emxArray_real_T *b, const emxArray_real_T *b_b,
                        const emxArray_real_T *fit)
{
  emxArray_real_T *c_b;
  const double *b_data;
  const double *fit_data;
  double *b_b_data;
  double *c_b_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  fit_data = fit->data;
  b_data = b_b->data;
  b_b_data = b->data;
  emxInit_real_T(&c_b, 1);
  i = c_b->size[0];
  if (fit->size[0] == 1) {
    if (b->size[0] == 1) {
      c_b->size[0] = b_b->size[0];
    } else {
      c_b->size[0] = b->size[0];
    }
  } else {
    c_b->size[0] = fit->size[0];
  }
  emxEnsureCapacity_real_T(c_b, i);
  c_b_data = c_b->data;
  stride_0_0 = (b_b->size[0] != 1);
  stride_1_0 = (b->size[0] != 1);
  stride_2_0 = (fit->size[0] != 1);
  if (fit->size[0] == 1) {
    if (b->size[0] == 1) {
      loop_ub = b_b->size[0];
    } else {
      loop_ub = b->size[0];
    }
  } else {
    loop_ub = fit->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    c_b_data[i] = (b_data[i * stride_0_0] - b_b_data[i * stride_1_0]) /
                  (fit_data[i * stride_2_0] + 1.0);
  }
  i = b->size[0];
  b->size[0] = c_b->size[0];
  emxEnsureCapacity_real_T(b, i);
  b_b_data = b->data;
  loop_ub = c_b->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_b_data[i] = c_b_data[i];
  }
  emxFree_real_T(&c_b);
}

/* End of code generation (div.c) */
