/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * find.c
 *
 * Code generation for function 'find'
 *
 */

/* Include files */
#include "find.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include <string.h>

/* Function Definitions */
void b_binary_expand_op(emxArray_int32_T *idx, const emxArray_real_T *indc,
                        const emxArray_real_T *BY)
{
  emxArray_boolean_T *b_indc;
  const double *BY_data;
  const double *indc_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  bool *b_indc_data;
  BY_data = BY->data;
  indc_data = indc->data;
  emxInit_boolean_T(&b_indc, 1);
  i = b_indc->size[0];
  if (BY->size[0] == 1) {
    b_indc->size[0] = indc->size[0];
  } else {
    b_indc->size[0] = BY->size[0];
  }
  emxEnsureCapacity_boolean_T(b_indc, i);
  b_indc_data = b_indc->data;
  stride_0_0 = (indc->size[0] != 1);
  stride_1_0 = (BY->size[0] != 1);
  if (BY->size[0] == 1) {
    loop_ub = indc->size[0];
  } else {
    loop_ub = BY->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_indc_data[i] =
        (indc_data[i * stride_0_0] + BY_data[i * stride_1_0] != 0.0);
  }
  b_eml_find(b_indc, idx);
  emxFree_boolean_T(&b_indc);
}

/*
 *
 */
void b_eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i)
{
  int idx;
  int ii;
  int nx;
  int *i_data;
  const bool *x_data;
  bool exitg1;
  x_data = x->data;
  nx = x->size[0];
  idx = 0;
  ii = i->size[0];
  i->size[0] = x->size[0];
  emxEnsureCapacity_int32_T(i, ii);
  i_data = i->data;
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx - 1)) {
    if (x_data[ii]) {
      idx++;
      i_data[idx - 1] = ii + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }
  if (x->size[0] == 1) {
    if (idx == 0) {
      i->size[0] = 0;
    }
  } else {
    ii = i->size[0];
    if (1 > idx) {
      i->size[0] = 0;
    } else {
      i->size[0] = idx;
    }
    emxEnsureCapacity_int32_T(i, ii);
  }
}

void e_binary_expand_op(emxArray_int32_T *idx, const emxArray_real_T *loc,
                        const emxArray_real_T *Z, double varargin_4)
{
  emxArray_boolean_T *b_loc;
  const double *Z_data;
  const double *loc_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  bool *b_loc_data;
  Z_data = Z->data;
  loc_data = loc->data;
  emxInit_boolean_T(&b_loc, 1);
  i = b_loc->size[0];
  if (Z->size[0] == 1) {
    b_loc->size[0] = loc->size[0];
  } else {
    b_loc->size[0] = Z->size[0];
  }
  emxEnsureCapacity_boolean_T(b_loc, i);
  b_loc_data = b_loc->data;
  stride_0_0 = (loc->size[0] != 1);
  stride_1_0 = (Z->size[0] != 1);
  if (Z->size[0] == 1) {
    loop_ub = loc->size[0];
  } else {
    loop_ub = Z->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_loc_data[i] =
        (loc_data[i * stride_0_0] / Z_data[i * stride_1_0] >= varargin_4);
  }
  b_eml_find(b_loc, idx);
  emxFree_boolean_T(&b_loc);
}

/*
 *
 */
void eml_find(const emxArray_boolean_T *x, emxArray_int32_T *i)
{
  int idx;
  int ii;
  int nx;
  int *i_data;
  const bool *x_data;
  bool exitg1;
  x_data = x->data;
  nx = x->size[1];
  idx = 0;
  ii = i->size[0] * i->size[1];
  i->size[0] = 1;
  i->size[1] = x->size[1];
  emxEnsureCapacity_int32_T(i, ii);
  i_data = i->data;
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx - 1)) {
    if (x_data[ii]) {
      idx++;
      i_data[idx - 1] = ii + 1;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        ii++;
      }
    } else {
      ii++;
    }
  }
  if (x->size[1] == 1) {
    if (idx == 0) {
      i->size[0] = 1;
      i->size[1] = 0;
    }
  } else {
    ii = i->size[0] * i->size[1];
    if (1 > idx) {
      i->size[1] = 0;
    } else {
      i->size[1] = idx;
    }
    emxEnsureCapacity_int32_T(i, ii);
  }
}

/* End of code generation (find.c) */
