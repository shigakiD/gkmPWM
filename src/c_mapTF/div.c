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
#include "mapTF2_ls_emxutil.h"
#include "mapTF2_ls_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void binary_expand_op(emxArray_real_T *f, const emxArray_real_T *c,
                      const emxArray_cell_wrap_0 *seqindmat, int b_I, int i,
                      const emxArray_real_T *minnorm,
                      const emxArray_real_T *maxnorm)
{
  const cell_wrap_0 *seqindmat_data;
  const double *c_data;
  const double *maxnorm_data;
  const double *minnorm_data;
  double *f_data;
  int b_i;
  int b_seqindmat;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  maxnorm_data = maxnorm->data;
  minnorm_data = minnorm->data;
  seqindmat_data = seqindmat->data;
  c_data = c->data;
  b_seqindmat = (int)seqindmat_data[b_I].f1->data[i];
  b_i = c->size[0];
  stride_0_0 = f->size[0];
  if (maxnorm->size[0] == 1) {
    if (minnorm->size[0] == 1) {
      f->size[0] = b_i;
    } else {
      f->size[0] = minnorm->size[0];
    }
  } else {
    f->size[0] = maxnorm->size[0];
  }
  emxEnsureCapacity_real_T(f, stride_0_0);
  f_data = f->data;
  stride_0_0 = (b_i != 1);
  stride_1_0 = (minnorm->size[0] != 1);
  stride_2_0 = (maxnorm->size[0] != 1);
  if (maxnorm->size[0] == 1) {
    if (minnorm->size[0] == 1) {
      loop_ub = b_i;
    } else {
      loop_ub = minnorm->size[0];
    }
  } else {
    loop_ub = maxnorm->size[0];
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    f_data[b_i] = (c_data[b_i * stride_0_0 + c->size[0] * (b_seqindmat - 1)] -
                   minnorm_data[b_i * stride_1_0]) /
                  maxnorm_data[b_i * stride_2_0];
  }
}

/* End of code generation (div.c) */
