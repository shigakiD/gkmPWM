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
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
void b_binary_expand_op(emxArray_real_T *ff, const emxArray_real_T *kmat,
                        const emxArray_cell_wrap_2 *seqindmat, int b_I, int i,
                        const emxArray_real_T *minnorm,
                        const emxArray_int32_T *p1, const emxArray_real_T *f)
{
  const cell_wrap_2 *seqindmat_data;
  const double *f_data;
  const double *kmat_data;
  const double *minnorm_data;
  double *ff_data;
  const int *p1_data;
  int b_i;
  int b_seqindmat;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  f_data = f->data;
  p1_data = p1->data;
  minnorm_data = minnorm->data;
  seqindmat_data = seqindmat->data;
  kmat_data = kmat->data;
  b_seqindmat = (int)seqindmat_data[b_I].f1->data[i];
  b_i = kmat->size[0];
  stride_0_0 = ff->size[0];
  if (p1->size[0] == 1) {
    ff->size[0] = b_i;
  } else {
    ff->size[0] = p1->size[0];
  }
  emxEnsureCapacity_real_T(ff, stride_0_0);
  ff_data = ff->data;
  stride_0_0 = (b_i != 1);
  stride_1_0 = (p1->size[0] != 1);
  stride_2_0 = (p1->size[0] != 1);
  if (p1->size[0] == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = p1->size[0];
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    ff_data[b_i] =
        (kmat_data[b_i * stride_0_0 + kmat->size[0] * (b_seqindmat - 1)] -
         minnorm_data[p1_data[b_i * stride_1_0] - 1]) /
        f_data[p1_data[b_i * stride_2_0] - 1];
  }
}

void binary_expand_op(emxArray_real_T *ff, const emxArray_real_T *kmat,
                      const emxArray_cell_wrap_2 *seqindmat, int b_I, int i,
                      const emxArray_real_T *minnorm, const emxArray_real_T *p2,
                      const emxArray_real_T *f)
{
  const cell_wrap_2 *seqindmat_data;
  const double *f_data;
  const double *kmat_data;
  const double *minnorm_data;
  const double *p2_data;
  double *ff_data;
  int b_i;
  int b_seqindmat;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  f_data = f->data;
  p2_data = p2->data;
  minnorm_data = minnorm->data;
  seqindmat_data = seqindmat->data;
  kmat_data = kmat->data;
  b_seqindmat =
      (int)seqindmat_data[b_I].f1->data[i + seqindmat_data[b_I].f1->size[0]];
  b_i = kmat->size[0];
  stride_0_0 = ff->size[0];
  if (p2->size[0] == 1) {
    ff->size[0] = b_i;
  } else {
    ff->size[0] = p2->size[0];
  }
  emxEnsureCapacity_real_T(ff, stride_0_0);
  ff_data = ff->data;
  stride_0_0 = (b_i != 1);
  stride_1_0 = (p2->size[0] != 1);
  stride_2_0 = (p2->size[0] != 1);
  if (p2->size[0] == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = p2->size[0];
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    ff_data[b_i] =
        (kmat_data[b_i * stride_0_0 + kmat->size[0] * (b_seqindmat - 1)] -
         minnorm_data[(int)p2_data[b_i * stride_1_0] - 1]) /
        f_data[(int)p2_data[b_i * stride_2_0] - 1];
  }
}

/* End of code generation (div.c) */
