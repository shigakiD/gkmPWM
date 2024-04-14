/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * bsxfun.c
 *
 * Code generation for function 'bsxfun'
 *
 */

/* Include files */
#include "bsxfun.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_bsxfun(const emxArray_real_T *a, const emxArray_real_T *b,
              emxArray_real_T *c)
{
  const double *a_data;
  const double *b_data;
  double *c_data;
  int acoef;
  int csz_idx_0;
  int i;
  int u1;
  b_data = b->data;
  a_data = a->data;
  acoef = b->size[0];
  u1 = a->size[0];
  if (acoef <= u1) {
    u1 = acoef;
  }
  if (a->size[0] == 1) {
    csz_idx_0 = b->size[0];
  } else if (a->size[0] == b->size[0]) {
    csz_idx_0 = a->size[0];
  } else {
    csz_idx_0 = u1;
  }
  i = c->size[0];
  acoef = b->size[0];
  u1 = a->size[0];
  if (acoef <= u1) {
    u1 = acoef;
  }
  if (a->size[0] == 1) {
    c->size[0] = b->size[0];
  } else if (a->size[0] == b->size[0]) {
    c->size[0] = a->size[0];
  } else {
    c->size[0] = u1;
  }
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  if (csz_idx_0 != 0) {
    acoef = (a->size[0] != 1);
    i = csz_idx_0 - 1;
    for (u1 = 0; u1 <= i; u1++) {
      c_data[u1] = a_data[acoef * u1] - b_data[u1];
    }
  }
}

/*
 *
 */
void bsxfun(const emxArray_real_T *a, const emxArray_real_T *b,
            emxArray_real_T *c)
{
  const double *a_data;
  const double *b_data;
  double *c_data;
  int acoef;
  int b_k;
  int bcoef;
  int csz_idx_1;
  int i;
  int k;
  int u0;
  int u1;
  b_data = b->data;
  a_data = a->data;
  u0 = b->size[1];
  u1 = a->size[1];
  if (u0 <= u1) {
    u1 = u0;
  }
  if (b->size[1] == 1) {
    csz_idx_1 = a->size[1];
  } else if (a->size[1] == 1) {
    csz_idx_1 = b->size[1];
  } else if (a->size[1] == b->size[1]) {
    csz_idx_1 = a->size[1];
  } else {
    csz_idx_1 = u1;
  }
  i = c->size[0] * c->size[1];
  c->size[0] = a->size[0];
  u0 = b->size[1];
  u1 = a->size[1];
  if (u0 <= u1) {
    u1 = u0;
  }
  if (b->size[1] == 1) {
    c->size[1] = a->size[1];
  } else if (a->size[1] == 1) {
    c->size[1] = b->size[1];
  } else if (a->size[1] == b->size[1]) {
    c->size[1] = a->size[1];
  } else {
    c->size[1] = u1;
  }
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  if (csz_idx_1 != 0) {
    acoef = (a->size[1] != 1);
    bcoef = (b->size[1] != 1);
    i = csz_idx_1 - 1;
    for (k = 0; k <= i; k++) {
      u0 = acoef * k;
      u1 = bcoef * k;
      csz_idx_1 = c->size[0] - 1;
      for (b_k = 0; b_k <= csz_idx_1; b_k++) {
        c_data[b_k + c->size[0] * k] =
            a_data[b_k + a->size[0] * u0] - b_data[u1];
      }
    }
  }
}

/*
 *
 */
void c_bsxfun(const emxArray_real_T *a, emxArray_real_T *c)
{
  const double *a_data;
  double *c_data;
  int acoef;
  int b_acoef;
  int b_k;
  int i;
  int i1;
  int k;
  int varargin_2;
  a_data = a->data;
  i = c->size[0] * c->size[1];
  c->size[0] = a->size[0];
  c->size[1] = a->size[1];
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  if ((a->size[0] != 0) && (a->size[1] != 0)) {
    acoef = (a->size[1] != 1);
    i = a->size[1] - 1;
    for (k = 0; k <= i; k++) {
      varargin_2 = acoef * k;
      b_acoef = (a->size[0] != 1);
      i1 = c->size[0] - 1;
      for (b_k = 0; b_k <= i1; b_k++) {
        c_data[b_k + c->size[0] * k] =
            a_data[b_acoef * b_k + a->size[0] * varargin_2];
      }
    }
  }
}

/*
 *
 */
void d_bsxfun(const emxArray_real_T *a, const emxArray_real_T *b,
              emxArray_real_T *c)
{
  const double *a_data;
  const double *b_data;
  double *c_data;
  int acoef;
  int b_k;
  int bcoef;
  int i;
  int i1;
  int k;
  int u1;
  b_data = b->data;
  a_data = a->data;
  i = c->size[0] * c->size[1];
  bcoef = b->size[0];
  u1 = a->size[0];
  if (bcoef <= u1) {
    u1 = bcoef;
  }
  if (a->size[0] == 1) {
    c->size[0] = b->size[0];
  } else if (a->size[0] == b->size[0]) {
    c->size[0] = a->size[0];
  } else {
    c->size[0] = u1;
  }
  c->size[1] = b->size[1];
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  bcoef = b->size[0];
  u1 = a->size[0];
  if (bcoef <= u1) {
    u1 = bcoef;
  }
  if (a->size[0] == 1) {
    u1 = b->size[0];
  } else if (a->size[0] == b->size[0]) {
    u1 = a->size[0];
  }
  if ((u1 != 0) && (b->size[1] != 0)) {
    bcoef = (b->size[1] != 1);
    i = b->size[1] - 1;
    for (k = 0; k <= i; k++) {
      u1 = bcoef * k;
      acoef = (a->size[0] != 1);
      i1 = c->size[0] - 1;
      for (b_k = 0; b_k <= i1; b_k++) {
        c_data[b_k + c->size[0] * k] =
            a_data[acoef * b_k] - b_data[b_k + b->size[0] * u1];
      }
    }
  }
}

/* End of code generation (bsxfun.c) */
