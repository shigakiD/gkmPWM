/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * genIndex.c
 *
 * Code generation for function 'genIndex'
 *
 */

/* Include files */
#include "genIndex.h"
#include "find.h"
#include "flip.h"
#include "fliplr.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include "minOrMax.h"
#include "nchoosek.h"
#include "sort.h"
#include "sum.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void b_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               unsigned int j, const emxArray_real_T *vec);

static void c_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               int i, const emxArray_real_T *vec);

static void d_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *C,
                               const emxArray_int32_T *f,
                               const emxArray_real_T *d, int i);

static void e_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               const emxArray_int32_T *ind,
                               const emxArray_int32_T *f, int i,
                               const emxArray_real_T *C);

static void h_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               const emxArray_int32_T *ind,
                               const emxArray_real_T *d, int i,
                               const emxArray_real_T *C);

/* Function Definitions */
static void b_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               unsigned int j, const emxArray_real_T *vec)
{
  const double *c_data;
  const double *vec_data;
  double *x_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  vec_data = vec->data;
  c_data = c->data;
  i = c->size[1];
  stride_0_1 = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (vec->size[1] == 1) {
    x->size[1] = i;
  } else {
    x->size[1] = vec->size[1];
  }
  emxEnsureCapacity_real_T(x, stride_0_1);
  x_data = x->data;
  stride_0_1 = (i != 1);
  stride_1_1 = (vec->size[1] != 1);
  if (vec->size[1] == 1) {
    loop_ub = i;
  } else {
    loop_ub = vec->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    x_data[i] = c_data[((int)j + c->size[0] * (i * stride_0_1)) - 1] -
                vec_data[i * stride_1_1];
  }
}

static void c_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               int i, const emxArray_real_T *vec)
{
  const double *c_data;
  const double *vec_data;
  double *x_data;
  int b_i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  vec_data = vec->data;
  c_data = c->data;
  b_i = c->size[1];
  stride_0_1 = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (vec->size[1] == 1) {
    x->size[1] = b_i;
  } else {
    x->size[1] = vec->size[1];
  }
  emxEnsureCapacity_real_T(x, stride_0_1);
  x_data = x->data;
  stride_0_1 = (b_i != 1);
  stride_1_1 = (vec->size[1] != 1);
  if (vec->size[1] == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = vec->size[1];
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    x_data[b_i] = c_data[i + c->size[0] * (b_i * stride_0_1)] -
                  vec_data[b_i * stride_1_1];
  }
}

static void d_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *C,
                               const emxArray_int32_T *f,
                               const emxArray_real_T *d, int i)
{
  emxArray_real_T *b_C;
  const double *C_data;
  const double *d_data;
  double *b_C_data;
  double *x_data;
  const int *f_data;
  int b_f;
  int b_i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  d_data = d->data;
  f_data = f->data;
  C_data = C->data;
  x_data = x->data;
  emxInit_real_T(&b_C, 2);
  b_f = f_data[(int)d_data[i] - 1];
  b_i = C->size[1];
  stride_0_1 = b_C->size[0] * b_C->size[1];
  b_C->size[0] = 1;
  if (x->size[1] == 1) {
    b_C->size[1] = b_i;
  } else {
    b_C->size[1] = x->size[1];
  }
  emxEnsureCapacity_real_T(b_C, stride_0_1);
  b_C_data = b_C->data;
  stride_0_1 = (b_i != 1);
  stride_1_1 = (x->size[1] != 1);
  if (x->size[1] == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = x->size[1];
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    b_C_data[b_i] =
        (C_data[(b_f + C->size[0] * (b_i * stride_0_1)) - 1] + 1.0) -
        x_data[b_i * stride_1_1];
  }
  b_i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = b_C->size[1];
  emxEnsureCapacity_real_T(x, b_i);
  x_data = x->data;
  loop_ub = b_C->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    x_data[b_i] = b_C_data[b_i];
  }
  emxFree_real_T(&b_C);
}

static void e_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               const emxArray_int32_T *ind,
                               const emxArray_int32_T *f, int i,
                               const emxArray_real_T *C)
{
  const double *C_data;
  const double *c_data;
  double *x_data;
  const int *f_data;
  const int *ind_data;
  int b_i;
  int b_ind;
  int i1;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  C_data = C->data;
  f_data = f->data;
  ind_data = ind->data;
  c_data = c->data;
  b_ind = ind_data[f_data[i] - 1];
  b_i = c->size[1];
  i1 = C->size[1];
  loop_ub = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (i1 == 1) {
    x->size[1] = b_i;
  } else {
    x->size[1] = i1;
  }
  emxEnsureCapacity_real_T(x, loop_ub);
  x_data = x->data;
  stride_0_1 = (b_i != 1);
  stride_1_1 = (i1 != 1);
  if (i1 == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = i1;
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    x_data[b_i] = c_data[(b_ind + c->size[0] * (b_i * stride_0_1)) - 1] -
                  C_data[i + C->size[0] * (b_i * stride_1_1)];
  }
}

static void h_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               const emxArray_int32_T *ind,
                               const emxArray_real_T *d, int i,
                               const emxArray_real_T *C)
{
  const double *C_data;
  const double *c_data;
  const double *d_data;
  double *x_data;
  const int *ind_data;
  int b_i;
  int b_ind;
  int i1;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  C_data = C->data;
  d_data = d->data;
  ind_data = ind->data;
  c_data = c->data;
  b_ind = ind_data[(int)d_data[i] - 1];
  b_i = c->size[1];
  i1 = C->size[1];
  loop_ub = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (i1 == 1) {
    x->size[1] = b_i;
  } else {
    x->size[1] = i1;
  }
  emxEnsureCapacity_real_T(x, loop_ub);
  x_data = x->data;
  stride_0_1 = (b_i != 1);
  stride_1_1 = (i1 != 1);
  if (i1 == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = i1;
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    x_data[b_i] = c_data[(b_ind + c->size[0] * (b_i * stride_0_1)) - 1] -
                  C_data[i + C->size[0] * (b_i * stride_1_1)];
  }
}

/*
 * function [c,C,I,ind,mat,rcnum] = genIndex(l,k)
 */
void b_genIndex(double l, double k, emxArray_real_T *c)
{
  emxArray_boolean_T *b_d;
  emxArray_boolean_T *b_e;
  emxArray_int32_T *e;
  emxArray_int32_T *f;
  emxArray_int32_T *iidx;
  emxArray_int32_T *ind;
  emxArray_int32_T *r;
  emxArray_real_T *C;
  emxArray_real_T *d;
  emxArray_real_T *first_occ;
  emxArray_real_T *mat;
  emxArray_real_T *varargin_2;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  double a;
  double xtmp;
  double *C_data;
  double *c_data;
  double *d_data;
  double *first_occ_data;
  double *vec_data;
  double *x_data;
  int b_i;
  unsigned int b_j;
  int i;
  int i1;
  int j;
  int j2;
  int n;
  int nd2;
  int nx;
  int *e_data;
  int *f_data;
  int *iidx_data;
  bool empty_non_axis_sizes;
  bool *b_e_data;
  /* l = k-mer length */
  /* k = # of ungapped positions */
  /*  Alternative to combnk */
  /*  c = combnk(1:l,k); */
  /* 'genIndex:7' c = flip(nchoosek(1:l,k)); */
  emxInit_real_T(&first_occ, 2);
  if (l < 1.0) {
    first_occ->size[0] = 1;
    first_occ->size[1] = 0;
  } else {
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    j2 = (int)floor(l - 1.0);
    first_occ->size[1] = j2 + 1;
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    for (i = 0; i <= j2; i++) {
      first_occ_data[i] = (double)i + 1.0;
    }
  }
  nchoosek(first_occ, k, c);
  flip(c);
  c_data = c->data;
  /* 'genIndex:9' d = ones(length(c),1); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    n = -1;
  } else if (c->size[0] > c->size[1]) {
    n = c->size[0] - 1;
  } else {
    n = c->size[1] - 1;
  }
  emxInit_real_T(&d, 1);
  i = d->size[0];
  d->size[0] = n + 1;
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  for (i = 0; i <= n; i++) {
    d_data[i] = 1.0;
  }
  emxInit_int32_T(&e, 2);
  /*  Alternative to continual growth */
  /*  e = []; */
  /* 'genIndex:13' e = ones(1,length(d)-1) * -1; */
  i = e->size[0] * e->size[1];
  e->size[0] = 1;
  e->size[1] = n;
  emxEnsureCapacity_int32_T(e, i);
  e_data = e->data;
  for (i = 0; i < n; i++) {
    e_data[i] = -1;
  }
  /* 'genIndex:15' a = 1; */
  a = 1.0;
  /* 'genIndex:16' for i = 1:length(d)-1 */
  emxInit_real_T(&vec, 2);
  emxInit_real_T(&x, 2);
  for (b_i = 0; b_i < n; b_i++) {
    /* 'genIndex:17' vec = l+1-fliplr(c(i,:)); */
    j2 = c->size[1];
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    vec->size[1] = c->size[1];
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    for (i = 0; i < j2; i++) {
      vec_data[i] = c_data[b_i + c->size[0] * i];
    }
    nd2 = c->size[1] >> 1;
    for (nx = 0; nx < nd2; nx++) {
      j2 = (c->size[1] - nx) - 1;
      xtmp = vec_data[nx];
      vec_data[nx] = vec_data[j2];
      vec_data[j2] = xtmp;
    }
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    j2 = vec->size[1] - 1;
    for (i = 0; i <= j2; i++) {
      vec_data[i] = (l + 1.0) - vec_data[i];
    }
    /* 'genIndex:18' if d(i) ~= 0 */
    if (d_data[b_i] != 0.0) {
      /* 'genIndex:19' if sum(abs(c(i,:)-vec))==0 */
      j2 = c->size[1];
      if (c->size[1] == vec->size[1]) {
        i = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = c->size[1];
        emxEnsureCapacity_real_T(x, i);
        x_data = x->data;
        for (i = 0; i < j2; i++) {
          x_data[i] = c_data[b_i + c->size[0] * i] - vec_data[i];
        }
      } else {
        c_binary_expand_op(x, c, b_i, vec);
        x_data = x->data;
      }
      nx = x->size[1];
      i = first_occ->size[0] * first_occ->size[1];
      first_occ->size[0] = 1;
      first_occ->size[1] = x->size[1];
      emxEnsureCapacity_real_T(first_occ, i);
      first_occ_data = first_occ->data;
      for (nd2 = 0; nd2 < nx; nd2++) {
        first_occ_data[nd2] = fabs(x_data[nd2]);
      }
      if (sum(first_occ) == 0.0) {
        /* 'genIndex:20' e(a) = i; */
        e_data[(int)a - 1] = b_i + 1;
        /* 'genIndex:21' d(i) = 0; */
        d_data[b_i] = 0.0;
        /* 'genIndex:22' a = a+1; */
        a++;
      } else {
        /* 'genIndex:23' else */
        /* 'genIndex:24' for j = i+1:length(d) */
        i = d->size[0] - b_i;
        for (j = 0; j <= i - 2; j++) {
          b_j = ((unsigned int)b_i + j) + 2U;
          /* 'genIndex:25' if sum(abs(c(j,:)-vec))==0 */
          j2 = c->size[1];
          if (c->size[1] == vec->size[1]) {
            i1 = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = c->size[1];
            emxEnsureCapacity_real_T(x, i1);
            x_data = x->data;
            for (i1 = 0; i1 < j2; i1++) {
              x_data[i1] =
                  c_data[((int)b_j + c->size[0] * i1) - 1] - vec_data[i1];
            }
          } else {
            b_binary_expand_op(x, c, b_j, vec);
            x_data = x->data;
          }
          nx = x->size[1];
          i1 = first_occ->size[0] * first_occ->size[1];
          first_occ->size[0] = 1;
          first_occ->size[1] = x->size[1];
          emxEnsureCapacity_real_T(first_occ, i1);
          first_occ_data = first_occ->data;
          for (nd2 = 0; nd2 < nx; nd2++) {
            first_occ_data[nd2] = fabs(x_data[nd2]);
          }
          if (sum(first_occ) == 0.0) {
            /* 'genIndex:26' d(j) = 0; */
            d_data[(int)b_j - 1] = 0.0;
          }
        }
      }
    }
  }
  emxInit_boolean_T(&b_e, 2);
  /*  Remove all trailing -1 */
  /* 'genIndex:33' first_occ = find(e==-1); */
  i = b_e->size[0] * b_e->size[1];
  b_e->size[0] = 1;
  b_e->size[1] = e->size[1];
  emxEnsureCapacity_boolean_T(b_e, i);
  b_e_data = b_e->data;
  j2 = e->size[1];
  for (i = 0; i < j2; i++) {
    b_e_data[i] = (e_data[i] == -1);
  }
  emxInit_int32_T(&r, 2);
  eml_find(b_e, r);
  e_data = r->data;
  i = first_occ->size[0] * first_occ->size[1];
  first_occ->size[0] = 1;
  first_occ->size[1] = r->size[1];
  emxEnsureCapacity_real_T(first_occ, i);
  first_occ_data = first_occ->data;
  j2 = r->size[1];
  emxFree_boolean_T(&b_e);
  for (i = 0; i < j2; i++) {
    first_occ_data[i] = e_data[i];
  }
  emxFree_int32_T(&r);
  emxInit_boolean_T(&b_d, 1);
  /* 'genIndex:34' e = e(1:first_occ(1)-1); */
  if (1.0 > (double)(int)first_occ_data[0] - 1.0) {
    nx = 0;
  } else {
    nx = (int)first_occ_data[0] - 1;
  }
  i = e->size[0] * e->size[1];
  if (1.0 > (double)(int)first_occ_data[0] - 1.0) {
    e->size[1] = 0;
  } else {
    e->size[1] = (int)first_occ_data[0] - 1;
  }
  emxEnsureCapacity_int32_T(e, i);
  e_data = e->data;
  /* 'genIndex:36' rcnum = a-1; */
  /* 'genIndex:37' f = find(d==1); */
  i = b_d->size[0];
  b_d->size[0] = d->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  j2 = d->size[0];
  for (i = 0; i < j2; i++) {
    b_e_data[i] = (d_data[i] == 1.0);
  }
  emxInit_int32_T(&f, 1);
  emxInit_int32_T(&iidx, 1);
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  j2 = iidx->size[0];
  for (i = 0; i < j2; i++) {
    f_data[i] = iidx_data[i];
  }
  /* 'genIndex:38' c=[c(f,:);c(e,:)]; */
  i = b_d->size[0];
  b_d->size[0] = d->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  j2 = d->size[0];
  for (i = 0; i < j2; i++) {
    b_e_data[i] = (d_data[i] == 1.0);
  }
  emxInit_real_T(&mat, 2);
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  j2 = c->size[1];
  i = mat->size[0] * mat->size[1];
  mat->size[0] = iidx->size[0];
  mat->size[1] = c->size[1];
  emxEnsureCapacity_real_T(mat, i);
  vec_data = mat->data;
  for (i = 0; i < j2; i++) {
    nd2 = iidx->size[0];
    for (i1 = 0; i1 < nd2; i1++) {
      vec_data[i1 + mat->size[0] * i] =
          c_data[(iidx_data[i1] + c->size[0] * i) - 1];
    }
  }
  emxInit_real_T(&varargin_2, 2);
  j2 = c->size[1];
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = nx;
  varargin_2->size[1] = c->size[1];
  emxEnsureCapacity_real_T(varargin_2, i);
  first_occ_data = varargin_2->data;
  for (i = 0; i < j2; i++) {
    for (i1 = 0; i1 < nx; i1++) {
      first_occ_data[i1 + varargin_2->size[0] * i] =
          c_data[(e_data[i1] + c->size[0] * i) - 1];
    }
  }
  emxFree_int32_T(&e);
  if ((f->size[0] != 0) && (c->size[1] != 0)) {
    nd2 = c->size[1];
  } else if ((nx != 0) && (c->size[1] != 0)) {
    nd2 = c->size[1];
  } else {
    if (c->size[1] > 0) {
      nd2 = c->size[1];
    } else {
      nd2 = 0;
    }
    if (c->size[1] > nd2) {
      nd2 = c->size[1];
    }
  }
  empty_non_axis_sizes = (nd2 == 0);
  if (empty_non_axis_sizes || ((f->size[0] != 0) && (c->size[1] != 0))) {
    j2 = f->size[0];
  } else {
    j2 = 0;
  }
  if ((!empty_non_axis_sizes) && ((nx == 0) || (c->size[1] == 0))) {
    nx = 0;
  }
  i = c->size[0] * c->size[1];
  c->size[0] = j2 + nx;
  c->size[1] = nd2;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  for (i = 0; i < nd2; i++) {
    for (i1 = 0; i1 < j2; i1++) {
      c_data[i1 + c->size[0] * i] = vec_data[i1 + j2 * i];
    }
  }
  for (i = 0; i < nd2; i++) {
    for (i1 = 0; i1 < nx; i1++) {
      c_data[(i1 + j2) + c->size[0] * i] = first_occ_data[i1 + nx * i];
    }
  }
  emxFree_real_T(&varargin_2);
  emxInit_real_T(&C, 2);
  /* 'genIndex:39' C = c; */
  i = C->size[0] * C->size[1];
  C->size[0] = c->size[0];
  C->size[1] = c->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  j2 = c->size[0] * c->size[1];
  for (i = 0; i < j2; i++) {
    C_data[i] = c_data[i];
  }
  /* 'genIndex:40' for i = 1:length(c) */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    n = 0;
  } else {
    nd2 = c->size[0];
    n = c->size[1];
    if (nd2 >= n) {
      n = nd2;
    }
  }
  for (b_i = 0; b_i < n; b_i++) {
    /* 'genIndex:41' c2 = l+1-fliplr(c(i,:)); */
    j2 = c->size[1];
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    first_occ->size[1] = c->size[1];
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    for (i = 0; i < j2; i++) {
      first_occ_data[i] = c_data[b_i + c->size[0] * i];
    }
    nd2 = c->size[1] >> 1;
    for (nx = 0; nx < nd2; nx++) {
      j2 = (c->size[1] - nx) - 1;
      xtmp = first_occ_data[nx];
      first_occ_data[nx] = first_occ_data[j2];
      first_occ_data[j2] = xtmp;
    }
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    j2 = first_occ->size[1] - 1;
    for (i = 0; i <= j2; i++) {
      first_occ_data[i] = (l + 1.0) - first_occ_data[i];
    }
    /* 'genIndex:42' if sum(0.5.^c(i,:)) < sum(0.5.^c2) */
    j2 = c->size[1];
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    vec->size[1] = c->size[1];
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    for (i = 0; i < j2; i++) {
      a = c_data[b_i + c->size[0] * i];
      vec_data[i] = pow(0.5, a);
    }
    i = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = first_occ->size[1];
    emxEnsureCapacity_real_T(x, i);
    x_data = x->data;
    j2 = first_occ->size[1];
    for (i = 0; i < j2; i++) {
      a = first_occ_data[i];
      x_data[i] = pow(0.5, a);
    }
    if (sum(vec) < sum(x)) {
      /* 'genIndex:43' C(i,:) = c2; */
      j2 = first_occ->size[1];
      for (i = 0; i < j2; i++) {
        C_data[b_i + C->size[0] * i] = first_occ_data[i];
      }
    }
  }
  emxFree_real_T(&vec);
  /* 'genIndex:46' c=C; */
  i = c->size[0] * c->size[1];
  c->size[0] = C->size[0];
  c->size[1] = C->size[1];
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  j2 = C->size[0] * C->size[1];
  for (i = 0; i < j2; i++) {
    c_data[i] = C_data[i];
  }
  /* 'genIndex:47' [~,ind] = sort(sum(0.5.^C,2),'descend'); */
  i = mat->size[0] * mat->size[1];
  mat->size[0] = C->size[0];
  mat->size[1] = C->size[1];
  emxEnsureCapacity_real_T(mat, i);
  vec_data = mat->data;
  j2 = C->size[0] * C->size[1];
  for (i = 0; i < j2; i++) {
    a = C_data[i];
    vec_data[i] = pow(0.5, a);
  }
  emxInit_int32_T(&ind, 1);
  b_sum(mat, d);
  sort(d, iidx);
  iidx_data = iidx->data;
  i = ind->size[0];
  ind->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(ind, i);
  e_data = ind->data;
  j2 = iidx->size[0];
  for (i = 0; i < j2; i++) {
    e_data[i] = iidx_data[i];
  }
  /* 'genIndex:48' C = C(ind,:); */
  nd2 = C->size[1] - 1;
  i = mat->size[0] * mat->size[1];
  mat->size[0] = ind->size[0];
  mat->size[1] = C->size[1];
  emxEnsureCapacity_real_T(mat, i);
  vec_data = mat->data;
  for (i = 0; i <= nd2; i++) {
    j2 = ind->size[0];
    for (i1 = 0; i1 < j2; i1++) {
      vec_data[i1 + mat->size[0] * i] =
          C_data[(e_data[i1] + C->size[0] * i) - 1];
    }
  }
  i = C->size[0] * C->size[1];
  C->size[0] = mat->size[0];
  C->size[1] = mat->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  j2 = mat->size[0] * mat->size[1];
  for (i = 0; i < j2; i++) {
    C_data[i] = vec_data[i];
  }
  /* 'genIndex:49' f = find(C(:,1)==1); */
  j2 = C->size[0];
  i = b_d->size[0];
  b_d->size[0] = C->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  for (i = 0; i < j2; i++) {
    b_e_data[i] = (C_data[i] == 1.0);
  }
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  j2 = iidx->size[0];
  for (i = 0; i < j2; i++) {
    f_data[i] = iidx_data[i];
  }
  /* 'genIndex:50' S = length(f); */
  /* 'genIndex:51' ff = find(C(f,end)~=l); */
  i = b_d->size[0];
  b_d->size[0] = f->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  j2 = f->size[0];
  for (i = 0; i < j2; i++) {
    b_e_data[i] =
        (C_data[(f_data[i] + C->size[0] * (C->size[1] - 1)) - 1] != l);
  }
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = d->size[0];
  d->size[0] = iidx->size[0];
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  j2 = iidx->size[0];
  emxFree_boolean_T(&b_d);
  for (i = 0; i < j2; i++) {
    d_data[i] = iidx_data[i];
  }
  /* 'genIndex:52' for i = 1:length(ff)-1 */
  i = d->size[0];
  for (b_i = 0; b_i <= i - 2; b_i++) {
    /* 'genIndex:53' for j = i+1:length(ff) */
    i1 = d->size[0] - b_i;
    for (j = 0; j <= i1 - 2; j++) {
      b_j = ((unsigned int)b_i + j) + 2U;
      /* 'genIndex:54' if sum(abs(C(f(ff(i)),:)+1-fliplr(l+1-C(f(ff(j)),:))))==0
       */
      j2 = C->size[1];
      nd2 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = C->size[1];
      emxEnsureCapacity_real_T(x, nd2);
      x_data = x->data;
      for (nd2 = 0; nd2 < j2; nd2++) {
        x_data[nd2] =
            (l + 1.0) -
            C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] + C->size[0] * nd2) -
                   1];
      }
      n = x->size[1] - 1;
      nd2 = x->size[1] >> 1;
      for (nx = 0; nx < nd2; nx++) {
        j2 = n - nx;
        xtmp = x_data[nx];
        x_data[nx] = x_data[j2];
        x_data[j2] = xtmp;
      }
      j2 = C->size[1];
      if (C->size[1] == x->size[1]) {
        nd2 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = C->size[1];
        emxEnsureCapacity_real_T(x, nd2);
        x_data = x->data;
        for (nd2 = 0; nd2 < j2; nd2++) {
          x_data[nd2] =
              (C_data[(f_data[(int)d_data[b_i] - 1] + C->size[0] * nd2) - 1] +
               1.0) -
              x_data[nd2];
        }
      } else {
        d_binary_expand_op(x, C, f, d, b_i);
        x_data = x->data;
      }
      nx = x->size[1];
      nd2 = first_occ->size[0] * first_occ->size[1];
      first_occ->size[0] = 1;
      first_occ->size[1] = x->size[1];
      emxEnsureCapacity_real_T(first_occ, nd2);
      first_occ_data = first_occ->data;
      for (nd2 = 0; nd2 < nx; nd2++) {
        first_occ_data[nd2] = fabs(x_data[nd2]);
      }
      if (sum(first_occ) == 0.0) {
        /* 'genIndex:55' C(f(ff(j)),:) = fliplr(l+1-C(f(ff(j)),:)); */
        j2 = C->size[1];
        nd2 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = C->size[1];
        emxEnsureCapacity_real_T(x, nd2);
        x_data = x->data;
        for (nd2 = 0; nd2 < j2; nd2++) {
          x_data[nd2] =
              (l + 1.0) - C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] +
                                  C->size[0] * nd2) -
                                 1];
        }
        n = x->size[1] - 1;
        nd2 = x->size[1] >> 1;
        for (nx = 0; nx < nd2; nx++) {
          j2 = n - nx;
          xtmp = x_data[nx];
          x_data[nx] = x_data[j2];
          x_data[j2] = xtmp;
        }
        j2 = x->size[1];
        for (nd2 = 0; nd2 < j2; nd2++) {
          C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] + C->size[0] * nd2) -
                 1] = x_data[nd2];
        }
      }
    }
  }
  emxFree_int32_T(&f);
  /* 'genIndex:59' [~,ind2] = sort(sum(0.5.^C,2),'descend'); */
  i = mat->size[0] * mat->size[1];
  mat->size[0] = C->size[0];
  mat->size[1] = C->size[1];
  emxEnsureCapacity_real_T(mat, i);
  vec_data = mat->data;
  j2 = C->size[0] * C->size[1];
  for (i = 0; i < j2; i++) {
    a = C_data[i];
    vec_data[i] = pow(0.5, a);
  }
  b_sum(mat, d);
  sort(d, iidx);
  iidx_data = iidx->data;
  i = d->size[0];
  d->size[0] = iidx->size[0];
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  j2 = iidx->size[0];
  for (i = 0; i < j2; i++) {
    d_data[i] = iidx_data[i];
  }
  emxFree_int32_T(&iidx);
  /* 'genIndex:60' C = C(ind2,:); */
  nd2 = C->size[1] - 1;
  i = mat->size[0] * mat->size[1];
  mat->size[0] = d->size[0];
  mat->size[1] = C->size[1];
  emxEnsureCapacity_real_T(mat, i);
  vec_data = mat->data;
  for (i = 0; i <= nd2; i++) {
    j2 = d->size[0];
    for (i1 = 0; i1 < j2; i1++) {
      vec_data[i1 + mat->size[0] * i] =
          C_data[((int)d_data[i1] + C->size[0] * i) - 1];
    }
  }
  i = C->size[0] * C->size[1];
  C->size[0] = mat->size[0];
  C->size[1] = mat->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  j2 = mat->size[0] * mat->size[1];
  for (i = 0; i < j2; i++) {
    C_data[i] = vec_data[i];
  }
  emxFree_real_T(&mat);
  /* 'genIndex:61' ind = ind(ind2); */
  /* 'genIndex:62' S = sum(C(:,1)==1); */
  /* 'genIndex:63' mat = zeros(S, max(c(:,1))-1); */
  /* 'genIndex:64' for i = 1:S */
  /* 'genIndex:80' mat = [(1:S)' mat]; */
  /* 'genIndex:81' I = zeros(length(C),1); */
  /* 'genIndex:82' I(1)=2; */
  /* 'genIndex:83' for i=2:length(I) */
  /* 'genIndex:93' for i = 1:length(c) */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    n = 0;
  } else {
    nd2 = c->size[0];
    n = c->size[1];
    if (nd2 >= n) {
      n = nd2;
    }
  }
  for (b_i = 0; b_i < n; b_i++) {
    /* 'genIndex:94' if sum(abs(c(ind(i),:)-C(i,:)))~=0 */
    j2 = c->size[1];
    if (c->size[1] == C->size[1]) {
      i = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = c->size[1];
      emxEnsureCapacity_real_T(x, i);
      x_data = x->data;
      for (i = 0; i < j2; i++) {
        x_data[i] =
            c_data[(e_data[(int)d_data[b_i] - 1] + c->size[0] * i) - 1] -
            C_data[b_i + C->size[0] * i];
      }
    } else {
      h_binary_expand_op(x, c, ind, d, b_i, C);
      x_data = x->data;
    }
    nx = x->size[1];
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    first_occ->size[1] = x->size[1];
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    for (nd2 = 0; nd2 < nx; nd2++) {
      first_occ_data[nd2] = fabs(x_data[nd2]);
    }
    if (sum(first_occ) != 0.0) {
      /* 'genIndex:95' c(ind(i),:) = l+1-fliplr(c(ind(i),:)); */
      j2 = c->size[1];
      i = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = c->size[1];
      emxEnsureCapacity_real_T(x, i);
      x_data = x->data;
      for (i = 0; i < j2; i++) {
        x_data[i] = c_data[(e_data[(int)d_data[b_i] - 1] + c->size[0] * i) - 1];
      }
      nd2 = c->size[1] >> 1;
      for (nx = 0; nx < nd2; nx++) {
        j2 = (c->size[1] - nx) - 1;
        xtmp = x_data[nx];
        x_data[nx] = x_data[j2];
        x_data[j2] = xtmp;
      }
      j2 = x->size[1];
      for (i = 0; i < j2; i++) {
        c_data[(e_data[(int)d_data[b_i] - 1] + c->size[0] * i) - 1] =
            (l + 1.0) - x_data[i];
      }
    }
  }
  emxFree_real_T(&x);
  emxFree_real_T(&C);
  emxFree_int32_T(&ind);
  emxFree_real_T(&first_occ);
  emxFree_real_T(&d);
}

/*
 * function [c,C,I,ind,mat,rcnum] = genIndex(l,k)
 */
void genIndex(double l, double k, emxArray_real_T *c, emxArray_real_T *C,
              emxArray_real_T *b_I, emxArray_real_T *ind, emxArray_real_T *mat,
              double *rcnum)
{
  emxArray_boolean_T *b_d;
  emxArray_boolean_T *b_e;
  emxArray_int32_T *b_ind;
  emxArray_int32_T *e;
  emxArray_int32_T *f;
  emxArray_int32_T *iidx;
  emxArray_int32_T *r;
  emxArray_real_T *d;
  emxArray_real_T *first_occ;
  emxArray_real_T *varargin_1;
  emxArray_real_T *varargin_2;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  double a;
  double b_l;
  double *C_data;
  double *c_data;
  double *d_data;
  double *first_occ_data;
  double *vec_data;
  double *x_data;
  int b_i;
  unsigned int b_j;
  int exitg1;
  int i;
  int i1;
  int input_sizes_idx_0;
  int j;
  int loop_ub;
  int nx;
  int nz;
  int result;
  int u1;
  int *e_data;
  int *f_data;
  int *iidx_data;
  signed char input_sizes_idx_1;
  bool empty_non_axis_sizes;
  bool *b_e_data;
  /* l = k-mer length */
  /* k = # of ungapped positions */
  /*  Alternative to combnk */
  /*  c = combnk(1:l,k); */
  /* 'genIndex:7' c = flip(nchoosek(1:l,k)); */
  emxInit_real_T(&first_occ, 2);
  if (l < 1.0) {
    first_occ->size[0] = 1;
    first_occ->size[1] = 0;
  } else {
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    loop_ub = (int)floor(l - 1.0);
    first_occ->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    for (i = 0; i <= loop_ub; i++) {
      first_occ_data[i] = (double)i + 1.0;
    }
  }
  nchoosek(first_occ, k, c);
  flip(c);
  c_data = c->data;
  /* 'genIndex:9' d = ones(length(c),1); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    input_sizes_idx_0 = -1;
  } else if (c->size[0] > c->size[1]) {
    input_sizes_idx_0 = c->size[0] - 1;
  } else {
    input_sizes_idx_0 = c->size[1] - 1;
  }
  emxInit_real_T(&d, 1);
  i = d->size[0];
  d->size[0] = input_sizes_idx_0 + 1;
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  for (i = 0; i <= input_sizes_idx_0; i++) {
    d_data[i] = 1.0;
  }
  emxInit_int32_T(&e, 2);
  /*  Alternative to continual growth */
  /*  e = []; */
  /* 'genIndex:13' e = ones(1,length(d)-1) * -1; */
  i = e->size[0] * e->size[1];
  e->size[0] = 1;
  e->size[1] = input_sizes_idx_0;
  emxEnsureCapacity_int32_T(e, i);
  e_data = e->data;
  for (i = 0; i < input_sizes_idx_0; i++) {
    e_data[i] = -1;
  }
  /* 'genIndex:15' a = 1; */
  a = 1.0;
  /* 'genIndex:16' for i = 1:length(d)-1 */
  emxInit_real_T(&vec, 2);
  emxInit_real_T(&x, 2);
  for (b_i = 0; b_i < input_sizes_idx_0; b_i++) {
    /* 'genIndex:17' vec = l+1-fliplr(c(i,:)); */
    loop_ub = c->size[1];
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    vec->size[1] = c->size[1];
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    for (i = 0; i < loop_ub; i++) {
      vec_data[i] = c_data[b_i + c->size[0] * i];
    }
    fliplr(vec);
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    loop_ub = vec->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      vec_data[i] = (l + 1.0) - vec_data[i];
    }
    /* 'genIndex:18' if d(i) ~= 0 */
    if (d_data[b_i] != 0.0) {
      /* 'genIndex:19' if sum(abs(c(i,:)-vec))==0 */
      loop_ub = c->size[1];
      if (c->size[1] == vec->size[1]) {
        i = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = c->size[1];
        emxEnsureCapacity_real_T(x, i);
        x_data = x->data;
        for (i = 0; i < loop_ub; i++) {
          x_data[i] = c_data[b_i + c->size[0] * i] - vec_data[i];
        }
      } else {
        c_binary_expand_op(x, c, b_i, vec);
        x_data = x->data;
      }
      nx = x->size[1];
      i = first_occ->size[0] * first_occ->size[1];
      first_occ->size[0] = 1;
      first_occ->size[1] = x->size[1];
      emxEnsureCapacity_real_T(first_occ, i);
      first_occ_data = first_occ->data;
      for (loop_ub = 0; loop_ub < nx; loop_ub++) {
        first_occ_data[loop_ub] = fabs(x_data[loop_ub]);
      }
      if (sum(first_occ) == 0.0) {
        /* 'genIndex:20' e(a) = i; */
        e_data[(int)a - 1] = b_i + 1;
        /* 'genIndex:21' d(i) = 0; */
        d_data[b_i] = 0.0;
        /* 'genIndex:22' a = a+1; */
        a++;
      } else {
        /* 'genIndex:23' else */
        /* 'genIndex:24' for j = i+1:length(d) */
        i = d->size[0] - b_i;
        for (j = 0; j <= i - 2; j++) {
          b_j = ((unsigned int)b_i + j) + 2U;
          /* 'genIndex:25' if sum(abs(c(j,:)-vec))==0 */
          loop_ub = c->size[1];
          if (c->size[1] == vec->size[1]) {
            i1 = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = c->size[1];
            emxEnsureCapacity_real_T(x, i1);
            x_data = x->data;
            for (i1 = 0; i1 < loop_ub; i1++) {
              x_data[i1] =
                  c_data[((int)b_j + c->size[0] * i1) - 1] - vec_data[i1];
            }
          } else {
            b_binary_expand_op(x, c, b_j, vec);
            x_data = x->data;
          }
          nx = x->size[1];
          i1 = first_occ->size[0] * first_occ->size[1];
          first_occ->size[0] = 1;
          first_occ->size[1] = x->size[1];
          emxEnsureCapacity_real_T(first_occ, i1);
          first_occ_data = first_occ->data;
          for (loop_ub = 0; loop_ub < nx; loop_ub++) {
            first_occ_data[loop_ub] = fabs(x_data[loop_ub]);
          }
          if (sum(first_occ) == 0.0) {
            /* 'genIndex:26' d(j) = 0; */
            d_data[(int)b_j - 1] = 0.0;
          }
        }
      }
    }
  }
  emxInit_boolean_T(&b_e, 2);
  /*  Remove all trailing -1 */
  /* 'genIndex:33' first_occ = find(e==-1); */
  i = b_e->size[0] * b_e->size[1];
  b_e->size[0] = 1;
  b_e->size[1] = e->size[1];
  emxEnsureCapacity_boolean_T(b_e, i);
  b_e_data = b_e->data;
  loop_ub = e->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (e_data[i] == -1);
  }
  emxInit_int32_T(&r, 2);
  eml_find(b_e, r);
  e_data = r->data;
  i = first_occ->size[0] * first_occ->size[1];
  first_occ->size[0] = 1;
  first_occ->size[1] = r->size[1];
  emxEnsureCapacity_real_T(first_occ, i);
  first_occ_data = first_occ->data;
  loop_ub = r->size[1];
  emxFree_boolean_T(&b_e);
  for (i = 0; i < loop_ub; i++) {
    first_occ_data[i] = e_data[i];
  }
  emxFree_int32_T(&r);
  emxInit_boolean_T(&b_d, 1);
  /* 'genIndex:34' e = e(1:first_occ(1)-1); */
  if (1.0 > (double)(int)first_occ_data[0] - 1.0) {
    nx = 0;
  } else {
    nx = (int)first_occ_data[0] - 1;
  }
  i = e->size[0] * e->size[1];
  if (1.0 > (double)(int)first_occ_data[0] - 1.0) {
    e->size[1] = 0;
  } else {
    e->size[1] = (int)first_occ_data[0] - 1;
  }
  emxEnsureCapacity_int32_T(e, i);
  e_data = e->data;
  /* 'genIndex:36' rcnum = a-1; */
  *rcnum = a - 1.0;
  /* 'genIndex:37' f = find(d==1); */
  i = b_d->size[0];
  b_d->size[0] = d->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  loop_ub = d->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (d_data[i] == 1.0);
  }
  emxInit_int32_T(&f, 1);
  emxInit_int32_T(&iidx, 1);
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = iidx_data[i];
  }
  /* 'genIndex:38' c=[c(f,:);c(e,:)]; */
  i = b_d->size[0];
  b_d->size[0] = d->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  loop_ub = d->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (d_data[i] == 1.0);
  }
  emxInit_real_T(&varargin_1, 2);
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  loop_ub = c->size[1];
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = iidx->size[0];
  varargin_1->size[1] = c->size[1];
  emxEnsureCapacity_real_T(varargin_1, i);
  x_data = varargin_1->data;
  for (i = 0; i < loop_ub; i++) {
    input_sizes_idx_0 = iidx->size[0];
    for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
      x_data[i1 + varargin_1->size[0] * i] =
          c_data[(iidx_data[i1] + c->size[0] * i) - 1];
    }
  }
  emxInit_real_T(&varargin_2, 2);
  loop_ub = c->size[1];
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = nx;
  varargin_2->size[1] = c->size[1];
  emxEnsureCapacity_real_T(varargin_2, i);
  vec_data = varargin_2->data;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < nx; i1++) {
      vec_data[i1 + varargin_2->size[0] * i] =
          c_data[(e_data[i1] + c->size[0] * i) - 1];
    }
  }
  emxFree_int32_T(&e);
  if ((f->size[0] != 0) && (c->size[1] != 0)) {
    result = c->size[1];
  } else if ((nx != 0) && (c->size[1] != 0)) {
    result = c->size[1];
  } else {
    if (c->size[1] > 0) {
      result = c->size[1];
    } else {
      result = 0;
    }
    if (c->size[1] > result) {
      result = c->size[1];
    }
  }
  empty_non_axis_sizes = (result == 0);
  if (empty_non_axis_sizes || ((f->size[0] != 0) && (c->size[1] != 0))) {
    input_sizes_idx_0 = f->size[0];
  } else {
    input_sizes_idx_0 = 0;
  }
  if ((!empty_non_axis_sizes) && ((nx == 0) || (c->size[1] == 0))) {
    nx = 0;
  }
  i = c->size[0] * c->size[1];
  c->size[0] = input_sizes_idx_0 + nx;
  c->size[1] = result;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  for (i = 0; i < result; i++) {
    for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
      c_data[i1 + c->size[0] * i] = x_data[i1 + input_sizes_idx_0 * i];
    }
  }
  for (i = 0; i < result; i++) {
    for (i1 = 0; i1 < nx; i1++) {
      c_data[(i1 + input_sizes_idx_0) + c->size[0] * i] = vec_data[i1 + nx * i];
    }
  }
  emxFree_real_T(&varargin_2);
  /* 'genIndex:39' C = c; */
  i = C->size[0] * C->size[1];
  C->size[0] = c->size[0];
  C->size[1] = c->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  loop_ub = c->size[0] * c->size[1];
  for (i = 0; i < loop_ub; i++) {
    C_data[i] = c_data[i];
  }
  /* 'genIndex:40' for i = 1:length(c) */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    result = 0;
  } else {
    input_sizes_idx_0 = c->size[0];
    result = c->size[1];
    if (input_sizes_idx_0 >= result) {
      result = input_sizes_idx_0;
    }
  }
  for (b_i = 0; b_i < result; b_i++) {
    /* 'genIndex:41' c2 = l+1-fliplr(c(i,:)); */
    loop_ub = c->size[1];
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    first_occ->size[1] = c->size[1];
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    for (i = 0; i < loop_ub; i++) {
      first_occ_data[i] = c_data[b_i + c->size[0] * i];
    }
    fliplr(first_occ);
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    loop_ub = first_occ->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      first_occ_data[i] = (l + 1.0) - first_occ_data[i];
    }
    /* 'genIndex:42' if sum(0.5.^c(i,:)) < sum(0.5.^c2) */
    loop_ub = c->size[1];
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    vec->size[1] = c->size[1];
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    for (i = 0; i < loop_ub; i++) {
      a = c_data[b_i + c->size[0] * i];
      vec_data[i] = pow(0.5, a);
    }
    i = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = first_occ->size[1];
    emxEnsureCapacity_real_T(x, i);
    x_data = x->data;
    loop_ub = first_occ->size[1];
    for (i = 0; i < loop_ub; i++) {
      a = first_occ_data[i];
      x_data[i] = pow(0.5, a);
    }
    if (sum(vec) < sum(x)) {
      /* 'genIndex:43' C(i,:) = c2; */
      loop_ub = first_occ->size[1];
      for (i = 0; i < loop_ub; i++) {
        C_data[b_i + C->size[0] * i] = first_occ_data[i];
      }
    }
  }
  /* 'genIndex:46' c=C; */
  i = c->size[0] * c->size[1];
  c->size[0] = C->size[0];
  c->size[1] = C->size[1];
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  loop_ub = C->size[0] * C->size[1];
  for (i = 0; i < loop_ub; i++) {
    c_data[i] = C_data[i];
  }
  /* 'genIndex:47' [~,ind] = sort(sum(0.5.^C,2),'descend'); */
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = C->size[0];
  varargin_1->size[1] = C->size[1];
  emxEnsureCapacity_real_T(varargin_1, i);
  x_data = varargin_1->data;
  loop_ub = C->size[0] * C->size[1];
  for (i = 0; i < loop_ub; i++) {
    a = C_data[i];
    x_data[i] = pow(0.5, a);
  }
  emxInit_int32_T(&b_ind, 1);
  b_sum(varargin_1, d);
  sort(d, iidx);
  iidx_data = iidx->data;
  i = b_ind->size[0];
  b_ind->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(b_ind, i);
  e_data = b_ind->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    e_data[i] = iidx_data[i];
  }
  /* 'genIndex:48' C = C(ind,:); */
  input_sizes_idx_0 = C->size[1] - 1;
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = b_ind->size[0];
  varargin_1->size[1] = C->size[1];
  emxEnsureCapacity_real_T(varargin_1, i);
  x_data = varargin_1->data;
  for (i = 0; i <= input_sizes_idx_0; i++) {
    loop_ub = b_ind->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      x_data[i1 + varargin_1->size[0] * i] =
          C_data[(e_data[i1] + C->size[0] * i) - 1];
    }
  }
  i = C->size[0] * C->size[1];
  C->size[0] = varargin_1->size[0];
  C->size[1] = varargin_1->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  loop_ub = varargin_1->size[0] * varargin_1->size[1];
  for (i = 0; i < loop_ub; i++) {
    C_data[i] = x_data[i];
  }
  /* 'genIndex:49' f = find(C(:,1)==1); */
  loop_ub = C->size[0];
  i = b_d->size[0];
  b_d->size[0] = C->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (C_data[i] == 1.0);
  }
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = iidx_data[i];
  }
  /* 'genIndex:50' S = length(f); */
  /* 'genIndex:51' ff = find(C(f,end)~=l); */
  i = b_d->size[0];
  b_d->size[0] = f->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] =
        (C_data[(f_data[i] + C->size[0] * (C->size[1] - 1)) - 1] != l);
  }
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = d->size[0];
  d->size[0] = iidx->size[0];
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_data[i] = iidx_data[i];
  }
  /* 'genIndex:52' for i = 1:length(ff)-1 */
  i = d->size[0];
  for (b_i = 0; b_i <= i - 2; b_i++) {
    /* 'genIndex:53' for j = i+1:length(ff) */
    i1 = d->size[0] - b_i;
    if (0 <= i1 - 2) {
      b_l = l + 1.0;
    }
    for (j = 0; j <= i1 - 2; j++) {
      b_j = ((unsigned int)b_i + j) + 2U;
      /* 'genIndex:54' if sum(abs(C(f(ff(i)),:)+1-fliplr(l+1-C(f(ff(j)),:))))==0
       */
      loop_ub = C->size[1];
      input_sizes_idx_0 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = C->size[1];
      emxEnsureCapacity_real_T(x, input_sizes_idx_0);
      x_data = x->data;
      for (input_sizes_idx_0 = 0; input_sizes_idx_0 < loop_ub;
           input_sizes_idx_0++) {
        x_data[input_sizes_idx_0] =
            b_l - C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] +
                          C->size[0] * input_sizes_idx_0) -
                         1];
      }
      fliplr(x);
      loop_ub = C->size[1];
      if (C->size[1] == x->size[1]) {
        input_sizes_idx_0 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = C->size[1];
        emxEnsureCapacity_real_T(x, input_sizes_idx_0);
        x_data = x->data;
        for (input_sizes_idx_0 = 0; input_sizes_idx_0 < loop_ub;
             input_sizes_idx_0++) {
          x_data[input_sizes_idx_0] = (C_data[(f_data[(int)d_data[b_i] - 1] +
                                               C->size[0] * input_sizes_idx_0) -
                                              1] +
                                       1.0) -
                                      x_data[input_sizes_idx_0];
        }
      } else {
        d_binary_expand_op(x, C, f, d, b_i);
        x_data = x->data;
      }
      nx = x->size[1];
      input_sizes_idx_0 = first_occ->size[0] * first_occ->size[1];
      first_occ->size[0] = 1;
      first_occ->size[1] = x->size[1];
      emxEnsureCapacity_real_T(first_occ, input_sizes_idx_0);
      first_occ_data = first_occ->data;
      for (loop_ub = 0; loop_ub < nx; loop_ub++) {
        first_occ_data[loop_ub] = fabs(x_data[loop_ub]);
      }
      if (sum(first_occ) == 0.0) {
        /* 'genIndex:55' C(f(ff(j)),:) = fliplr(l+1-C(f(ff(j)),:)); */
        loop_ub = C->size[1];
        input_sizes_idx_0 = vec->size[0] * vec->size[1];
        vec->size[0] = 1;
        vec->size[1] = C->size[1];
        emxEnsureCapacity_real_T(vec, input_sizes_idx_0);
        vec_data = vec->data;
        for (input_sizes_idx_0 = 0; input_sizes_idx_0 < loop_ub;
             input_sizes_idx_0++) {
          vec_data[input_sizes_idx_0] =
              (l + 1.0) - C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] +
                                  C->size[0] * input_sizes_idx_0) -
                                 1];
        }
        fliplr(vec);
        vec_data = vec->data;
        loop_ub = vec->size[1];
        for (input_sizes_idx_0 = 0; input_sizes_idx_0 < loop_ub;
             input_sizes_idx_0++) {
          C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] +
                  C->size[0] * input_sizes_idx_0) -
                 1] = vec_data[input_sizes_idx_0];
        }
      }
    }
  }
  /* 'genIndex:59' [~,ind2] = sort(sum(0.5.^C,2),'descend'); */
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = C->size[0];
  varargin_1->size[1] = C->size[1];
  emxEnsureCapacity_real_T(varargin_1, i);
  x_data = varargin_1->data;
  loop_ub = C->size[0] * C->size[1];
  for (i = 0; i < loop_ub; i++) {
    a = C_data[i];
    x_data[i] = pow(0.5, a);
  }
  b_sum(varargin_1, d);
  sort(d, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = iidx_data[i];
  }
  emxFree_int32_T(&iidx);
  /* 'genIndex:60' C = C(ind2,:); */
  input_sizes_idx_0 = C->size[1] - 1;
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = f->size[0];
  varargin_1->size[1] = C->size[1];
  emxEnsureCapacity_real_T(varargin_1, i);
  x_data = varargin_1->data;
  for (i = 0; i <= input_sizes_idx_0; i++) {
    loop_ub = f->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      x_data[i1 + varargin_1->size[0] * i] =
          C_data[(f_data[i1] + C->size[0] * i) - 1];
    }
  }
  i = C->size[0] * C->size[1];
  C->size[0] = varargin_1->size[0];
  C->size[1] = varargin_1->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  loop_ub = varargin_1->size[0] * varargin_1->size[1];
  for (i = 0; i < loop_ub; i++) {
    C_data[i] = x_data[i];
  }
  /* 'genIndex:61' ind = ind(ind2); */
  i = ind->size[0];
  ind->size[0] = f->size[0];
  emxEnsureCapacity_real_T(ind, i);
  vec_data = ind->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = e_data[f_data[i] - 1];
  }
  /* 'genIndex:62' S = sum(C(:,1)==1); */
  loop_ub = C->size[0];
  i = b_d->size[0];
  b_d->size[0] = C->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (C_data[i] == 1.0);
  }
  input_sizes_idx_0 = b_d->size[0];
  if (b_d->size[0] == 0) {
    nz = 0;
  } else {
    nz = b_e_data[0];
    for (loop_ub = 2; loop_ub <= input_sizes_idx_0; loop_ub++) {
      nz += b_e_data[loop_ub - 1];
    }
  }
  emxFree_boolean_T(&b_d);
  /* 'genIndex:63' mat = zeros(S, max(c(:,1))-1); */
  loop_ub = c->size[0];
  i = d->size[0];
  d->size[0] = c->size[0];
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  for (i = 0; i < loop_ub; i++) {
    d_data[i] = c_data[i];
  }
  a = maximum(d);
  i = mat->size[0] * mat->size[1];
  mat->size[0] = nz;
  mat->size[1] = (int)(a - 1.0);
  emxEnsureCapacity_real_T(mat, i);
  vec_data = mat->data;
  loop_ub = nz * (int)(a - 1.0);
  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = 0.0;
  }
  /* 'genIndex:64' for i = 1:S */
  for (b_i = 0; b_i < nz; b_i++) {
    /* 'genIndex:65' if C(i,end) ~= l */
    if (C_data[b_i + C->size[0] * (C->size[1] - 1)] != l) {
      /* 'genIndex:66' for j = 2:max(c(:,1)) */
      loop_ub = c->size[0];
      i = d->size[0];
      d->size[0] = c->size[0];
      emxEnsureCapacity_real_T(d, i);
      d_data = d->data;
      for (i = 0; i < loop_ub; i++) {
        d_data[i] = c_data[i];
      }
      i = (int)(maximum(d) + -1.0);
      if (0 <= i - 1) {
        if ((c->size[0] == 0) || (c->size[1] == 0)) {
          u1 = 0;
        } else {
          input_sizes_idx_0 = c->size[0];
          u1 = c->size[1];
          if (input_sizes_idx_0 >= u1) {
            u1 = input_sizes_idx_0;
          }
        }
      }
      for (j = 0; j < i; j++) {
        /* 'genIndex:67' a = S+1; */
        a = (double)nz + 1.0;
        /* 'genIndex:68' while sum(abs(C(i,:)+j-1-C(a,:))) ~= 0 */
        do {
          exitg1 = 0;
          loop_ub = C->size[1];
          i1 = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = C->size[1];
          emxEnsureCapacity_real_T(x, i1);
          x_data = x->data;
          for (i1 = 0; i1 < loop_ub; i1++) {
            x_data[i1] =
                ((C_data[b_i + C->size[0] * i1] + ((double)j + 2.0)) - 1.0) -
                C_data[((int)a + C->size[0] * i1) - 1];
          }
          nx = x->size[1];
          i1 = first_occ->size[0] * first_occ->size[1];
          first_occ->size[0] = 1;
          first_occ->size[1] = x->size[1];
          emxEnsureCapacity_real_T(first_occ, i1);
          first_occ_data = first_occ->data;
          for (loop_ub = 0; loop_ub < nx; loop_ub++) {
            first_occ_data[loop_ub] = fabs(x_data[loop_ub]);
          }
          if (sum(first_occ) != 0.0) {
            /* 'genIndex:69' a = a + 1; */
            a++;
            /* 'genIndex:70' if a > length(c) */
            if ((c->size[0] == 0) || (c->size[1] == 0)) {
              result = 0;
            } else {
              input_sizes_idx_0 = c->size[0];
              result = c->size[1];
              if (input_sizes_idx_0 >= result) {
                result = input_sizes_idx_0;
              }
            }
            if ((unsigned int)a > (unsigned int)result) {
              exitg1 = 1;
            }
          } else {
            exitg1 = 1;
          }
        } while (exitg1 == 0);
        /* 'genIndex:74' if a <= length(c) */
        if ((unsigned int)a <= (unsigned int)u1) {
          /* 'genIndex:75' mat(i,j-1) = a; */
          vec_data[b_i + mat->size[0] * j] = a;
        }
      }
    }
  }
  /* 'genIndex:80' mat = [(1:S)' mat]; */
  if (nz < 1) {
    first_occ->size[0] = 1;
    first_occ->size[1] = 0;
  } else {
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    first_occ->size[1] = nz;
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    loop_ub = nz - 1;
    for (i = 0; i <= loop_ub; i++) {
      first_occ_data[i] = (double)i + 1.0;
    }
  }
  i = d->size[0];
  d->size[0] = first_occ->size[1];
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  loop_ub = first_occ->size[1];
  for (i = 0; i < loop_ub; i++) {
    d_data[i] = first_occ_data[i];
  }
  if (d->size[0] != 0) {
    result = d->size[0];
  } else if ((mat->size[0] != 0) && (mat->size[1] != 0)) {
    result = mat->size[0];
  } else {
    result = 0;
    if (mat->size[0] > 0) {
      result = mat->size[0];
    }
  }
  empty_non_axis_sizes = (result == 0);
  if (empty_non_axis_sizes || (d->size[0] != 0)) {
    input_sizes_idx_1 = 1;
  } else {
    input_sizes_idx_1 = 0;
  }
  if (empty_non_axis_sizes || ((mat->size[0] != 0) && (mat->size[1] != 0))) {
    input_sizes_idx_0 = mat->size[1];
  } else {
    input_sizes_idx_0 = 0;
  }
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = result;
  varargin_1->size[1] = input_sizes_idx_1 + input_sizes_idx_0;
  emxEnsureCapacity_real_T(varargin_1, i);
  x_data = varargin_1->data;
  loop_ub = input_sizes_idx_1;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < result; i1++) {
      x_data[i1] = d_data[i1];
    }
  }
  emxFree_real_T(&d);
  for (i = 0; i < input_sizes_idx_0; i++) {
    for (i1 = 0; i1 < result; i1++) {
      x_data[i1 + varargin_1->size[0] * (i + input_sizes_idx_1)] =
          vec_data[i1 + result * i];
    }
  }
  i = mat->size[0] * mat->size[1];
  mat->size[0] = varargin_1->size[0];
  mat->size[1] = varargin_1->size[1];
  emxEnsureCapacity_real_T(mat, i);
  vec_data = mat->data;
  loop_ub = varargin_1->size[0] * varargin_1->size[1];
  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = x_data[i];
  }
  emxFree_real_T(&varargin_1);
  /* 'genIndex:81' I = zeros(length(C),1); */
  if ((C->size[0] == 0) || (C->size[1] == 0)) {
    result = 0;
  } else {
    input_sizes_idx_0 = C->size[0];
    result = C->size[1];
    if (input_sizes_idx_0 >= result) {
      result = input_sizes_idx_0;
    }
  }
  i = b_I->size[0];
  b_I->size[0] = result;
  emxEnsureCapacity_real_T(b_I, i);
  vec_data = b_I->data;
  for (i = 0; i < result; i++) {
    vec_data[i] = 0.0;
  }
  /* 'genIndex:82' I(1)=2; */
  vec_data[0] = 2.0;
  /* 'genIndex:83' for i=2:length(I) */
  for (b_i = 0; b_i <= result - 2; b_i++) {
    /* 'genIndex:84' a = 1; */
    /* 'genIndex:85' while C(i-1,a)==C(i,a) */
    for (a = 1.0; C_data[b_i + C->size[0] * ((int)a - 1)] ==
                  C_data[(b_i + C->size[0] * ((int)a - 1)) + 1];
         a++) {
      /* 'genIndex:86' a = a+1; */
    }
    /* 'genIndex:88' I(i) = a; */
    vec_data[b_i + 1] = a;
    /* 'genIndex:89' if I(i) < 2 */
    if (a < 2.0) {
      /* 'genIndex:90' I(i)=2; */
      vec_data[b_i + 1] = 2.0;
    }
  }
  /* 'genIndex:93' for i = 1:length(c) */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    result = 0;
  } else {
    input_sizes_idx_0 = c->size[0];
    result = c->size[1];
    if (input_sizes_idx_0 >= result) {
      result = input_sizes_idx_0;
    }
  }
  for (b_i = 0; b_i < result; b_i++) {
    /* 'genIndex:94' if sum(abs(c(ind(i),:)-C(i,:)))~=0 */
    loop_ub = c->size[1];
    if (c->size[1] == C->size[1]) {
      i = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = c->size[1];
      emxEnsureCapacity_real_T(x, i);
      x_data = x->data;
      for (i = 0; i < loop_ub; i++) {
        x_data[i] = c_data[(e_data[f_data[b_i] - 1] + c->size[0] * i) - 1] -
                    C_data[b_i + C->size[0] * i];
      }
    } else {
      e_binary_expand_op(x, c, b_ind, f, b_i, C);
      x_data = x->data;
    }
    nx = x->size[1];
    i = first_occ->size[0] * first_occ->size[1];
    first_occ->size[0] = 1;
    first_occ->size[1] = x->size[1];
    emxEnsureCapacity_real_T(first_occ, i);
    first_occ_data = first_occ->data;
    for (loop_ub = 0; loop_ub < nx; loop_ub++) {
      first_occ_data[loop_ub] = fabs(x_data[loop_ub]);
    }
    if (sum(first_occ) != 0.0) {
      /* 'genIndex:95' c(ind(i),:) = l+1-fliplr(c(ind(i),:)); */
      loop_ub = c->size[1];
      i = vec->size[0] * vec->size[1];
      vec->size[0] = 1;
      vec->size[1] = c->size[1];
      emxEnsureCapacity_real_T(vec, i);
      vec_data = vec->data;
      for (i = 0; i < loop_ub; i++) {
        vec_data[i] = c_data[(e_data[f_data[b_i] - 1] + c->size[0] * i) - 1];
      }
      fliplr(vec);
      vec_data = vec->data;
      loop_ub = vec->size[1];
      for (i = 0; i < loop_ub; i++) {
        c_data[(e_data[f_data[b_i] - 1] + c->size[0] * i) - 1] =
            (l + 1.0) - vec_data[i];
      }
    }
  }
  emxFree_real_T(&x);
  emxFree_int32_T(&b_ind);
  emxFree_int32_T(&f);
  emxFree_real_T(&first_occ);
  emxFree_real_T(&vec);
}

/* End of code generation (genIndex.c) */