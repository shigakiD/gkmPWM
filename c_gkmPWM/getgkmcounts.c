/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getgkmcounts.c
 *
 * Code generation for function 'getgkmcounts'
 *
 */

/* Include files */
#include "getgkmcounts.h"
#include "blockedSummation.h"
#include "feof.h"
#include "fgetl.h"
#include "fgets.h"
#include "fileManager.h"
#include "frewind.h"
#include "fseek.h"
#include "ftell.h"
#include "genIndex.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include "repmat.h"
#include "rot90.h"
#include "str2double.h"
#include "strcmp.h"
#include "strip.h"
#include "strtok.h"
#include "sum.h"
#include "cblas.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Function Declarations */
static void encodekmers(double l, double k, emxArray_real_T *mat);

static void encodekmers_norc(double l, double k, emxArray_real_T *mat);

static void f_binary_expand_op(emxArray_real_T *gkmc, double a,
                               const emxArray_real_T *temp);

static void g_binary_expand_op(emxArray_real_T *temp, int i1, int i2, int i3);

/* Function Definitions */
/*
 * function mat = encodekmers(l,k)
 */
static void encodekmers(double l, double k, emxArray_real_T *mat)
{
  emxArray_real_T *a;
  emxArray_real_T *c;
  emxArray_real_T *seqvec;
  emxArray_real_T *seqvec2;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  emxArray_real_T *y;
  double c_tmp;
  double xtmp;
  double *a_data;
  double *c_data;
  double *mat_data;
  double *seqvec2_data;
  double *seqvec_data;
  double *vec_data;
  double *y_data;
  int b_j1;
  int i;
  int j2;
  int nd2;
  int nx;
  unsigned int u;
  int u1;
  emxInit_real_T(&c, 2);
  /* 'getgkmcounts:261' c = genIndex(l,k); */
  b_genIndex(l, k, c);
  c_data = c->data;
  /* 'getgkmcounts:262' lcnum = length(c); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    u1 = 0;
  } else {
    nx = c->size[0];
    u1 = c->size[1];
    if (nx >= u1) {
      u1 = nx;
    }
  }
  emxInit_real_T(&seqvec, 2);
  /* 'getgkmcounts:263' seqvec = zeros(4^l, l); */
  c_tmp = pow(4.0, l);
  j2 = seqvec->size[0] * seqvec->size[1];
  seqvec->size[0] = (int)c_tmp;
  b_j1 = (int)l;
  seqvec->size[1] = (int)l;
  emxEnsureCapacity_real_T(seqvec, j2);
  seqvec_data = seqvec->data;
  nx = (int)c_tmp * (int)l;
  for (j2 = 0; j2 < nx; j2++) {
    seqvec_data[j2] = 0.0;
  }
  /* 'getgkmcounts:264' vec = (1:4^l)'-1; */
  emxInit_real_T(&y, 2);
  y_data = y->data;
  if (c_tmp < 1.0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    j2 = y->size[0] * y->size[1];
    y->size[0] = 1;
    nx = (int)floor(c_tmp - 1.0);
    y->size[1] = nx + 1;
    emxEnsureCapacity_real_T(y, j2);
    y_data = y->data;
    for (j2 = 0; j2 <= nx; j2++) {
      y_data[j2] = (double)j2 + 1.0;
    }
  }
  emxInit_real_T(&vec, 1);
  j2 = vec->size[0];
  vec->size[0] = y->size[1];
  emxEnsureCapacity_real_T(vec, j2);
  vec_data = vec->data;
  nx = y->size[1];
  for (j2 = 0; j2 < nx; j2++) {
    vec_data[j2] = y_data[j2] - 1.0;
  }
  /* 'getgkmcounts:265' for i = 1:l */
  emxInit_real_T(&x, 1);
  for (i = 0; i < b_j1; i++) {
    /* 'getgkmcounts:266' seqvec(:,i) = mod(floor(vec/4^(i-1)), 4); */
    xtmp = pow(4.0, ((double)i + 1.0) - 1.0);
    nx = vec->size[0];
    j2 = x->size[0];
    x->size[0] = vec->size[0];
    emxEnsureCapacity_real_T(x, j2);
    y_data = x->data;
    for (j2 = 0; j2 < nx; j2++) {
      y_data[j2] = vec_data[j2] / xtmp;
    }
    nx = x->size[0];
    for (nd2 = 0; nd2 < nx; nd2++) {
      y_data[nd2] = floor(y_data[nd2]);
    }
    nx = x->size[0];
    for (j2 = 0; j2 < nx; j2++) {
      xtmp = y_data[j2];
      if (xtmp == 0.0) {
        xtmp = 0.0;
      } else {
        xtmp = fmod(xtmp, 4.0);
        if (xtmp == 0.0) {
          xtmp = 0.0;
        }
      }
      seqvec_data[j2 + seqvec->size[0] * i] = xtmp;
    }
  }
  emxInit_real_T(&seqvec2, 2);
  /* 'getgkmcounts:268' seqvec2 = 3-fliplr(seqvec); */
  j2 = seqvec2->size[0] * seqvec2->size[1];
  seqvec2->size[0] = seqvec->size[0];
  seqvec2->size[1] = seqvec->size[1];
  emxEnsureCapacity_real_T(seqvec2, j2);
  seqvec2_data = seqvec2->data;
  nx = seqvec->size[0] * seqvec->size[1];
  for (j2 = 0; j2 < nx; j2++) {
    seqvec2_data[j2] = seqvec_data[j2];
  }
  nx = seqvec->size[0];
  nd2 = seqvec->size[1] >> 1;
  for (b_j1 = 0; b_j1 < nd2; b_j1++) {
    j2 = (seqvec->size[1] - b_j1) - 1;
    for (i = 0; i < nx; i++) {
      xtmp = seqvec2_data[i + seqvec2->size[0] * b_j1];
      seqvec2_data[i + seqvec2->size[0] * b_j1] =
          seqvec2_data[i + seqvec2->size[0] * j2];
      seqvec2_data[i + seqvec2->size[0] * j2] = xtmp;
    }
  }
  nx = seqvec2->size[0] * seqvec2->size[1];
  for (j2 = 0; j2 < nx; j2++) {
    seqvec2_data[j2] = 3.0 - seqvec2_data[j2];
  }
  /* 'getgkmcounts:269' mat = zeros(4^l, 2*lcnum); */
  j2 = mat->size[0] * mat->size[1];
  mat->size[0] = (int)c_tmp;
  mat->size[1] = 2 * u1;
  emxEnsureCapacity_real_T(mat, j2);
  mat_data = mat->data;
  nx = (int)c_tmp * (2 * u1);
  for (j2 = 0; j2 < nx; j2++) {
    mat_data[j2] = 0.0;
  }
  /* 'getgkmcounts:270' pow = 4.^(0:k-1)'; */
  if (k - 1.0 < 0.0) {
    y->size[1] = 0;
  } else {
    j2 = y->size[0] * y->size[1];
    y->size[0] = 1;
    nx = (int)floor(k - 1.0);
    y->size[1] = nx + 1;
    emxEnsureCapacity_real_T(y, j2);
    y_data = y->data;
    for (j2 = 0; j2 <= nx; j2++) {
      y_data[j2] = j2;
    }
  }
  j2 = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity_real_T(y, j2);
  y_data = y->data;
  nx = y->size[1] - 1;
  for (j2 = 0; j2 <= nx; j2++) {
    xtmp = y_data[j2];
    y_data[j2] = pow(4.0, xtmp);
  }
  j2 = vec->size[0];
  vec->size[0] = y->size[1];
  emxEnsureCapacity_real_T(vec, j2);
  vec_data = vec->data;
  nx = y->size[1];
  for (j2 = 0; j2 < nx; j2++) {
    vec_data[j2] = y_data[j2];
  }
  emxFree_real_T(&y);
  /* 'getgkmcounts:271' for i = 1:lcnum */
  emxInit_real_T(&a, 2);
  for (i = 0; i < u1; i++) {
    /* 'getgkmcounts:272' mat(:,i) = seqvec(:,c(i,:))*pow+4^k*(i-1)+1; */
    nx = seqvec->size[0];
    nd2 = c->size[1];
    j2 = a->size[0] * a->size[1];
    a->size[0] = seqvec->size[0];
    a->size[1] = c->size[1];
    emxEnsureCapacity_real_T(a, j2);
    a_data = a->data;
    for (j2 = 0; j2 < nd2; j2++) {
      for (b_j1 = 0; b_j1 < nx; b_j1++) {
        a_data[b_j1 + a->size[0] * j2] =
            seqvec_data[b_j1 + seqvec->size[0] *
                                   ((int)c_data[i + c->size[0] * j2] - 1)];
      }
    }
    nx = seqvec->size[0];
    if ((seqvec->size[0] == 0) || (c->size[1] == 0) || (vec->size[0] == 0)) {
      j2 = x->size[0];
      x->size[0] = seqvec->size[0];
      emxEnsureCapacity_real_T(x, j2);
      y_data = x->data;
      for (j2 = 0; j2 < nx; j2++) {
        y_data[j2] = 0.0;
      }
    } else {
      j2 = x->size[0];
      x->size[0] = seqvec->size[0];
      emxEnsureCapacity_real_T(x, j2);
      y_data = x->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)seqvec->size[0], (blasint)1, (blasint)c->size[1],
                  1.0, &a_data[0], (blasint)seqvec->size[0], &vec_data[0],
                  (blasint)vec->size[0], 0.0, &y_data[0],
                  (blasint)seqvec->size[0]);
    }
    c_tmp = pow(4.0, k);
    xtmp = c_tmp * (((double)i + 1.0) - 1.0);
    nx = x->size[0];
    for (j2 = 0; j2 < nx; j2++) {
      mat_data[j2 + mat->size[0] * i] = (y_data[j2] + xtmp) + 1.0;
    }
    /* 'getgkmcounts:273' mat(:,i+lcnum) =
     * seqvec2(:,c(i,:))*pow+4^k*(lcnum+i-1)+1; */
    u = (unsigned int)i + u1;
    nx = seqvec2->size[0];
    nd2 = c->size[1];
    j2 = a->size[0] * a->size[1];
    a->size[0] = seqvec2->size[0];
    a->size[1] = c->size[1];
    emxEnsureCapacity_real_T(a, j2);
    a_data = a->data;
    for (j2 = 0; j2 < nd2; j2++) {
      for (b_j1 = 0; b_j1 < nx; b_j1++) {
        a_data[b_j1 + a->size[0] * j2] =
            seqvec2_data[b_j1 + seqvec2->size[0] *
                                    ((int)c_data[i + c->size[0] * j2] - 1)];
      }
    }
    nx = seqvec2->size[0];
    if ((seqvec2->size[0] == 0) || (c->size[1] == 0) || (vec->size[0] == 0)) {
      j2 = x->size[0];
      x->size[0] = seqvec2->size[0];
      emxEnsureCapacity_real_T(x, j2);
      y_data = x->data;
      for (j2 = 0; j2 < nx; j2++) {
        y_data[j2] = 0.0;
      }
    } else {
      j2 = x->size[0];
      x->size[0] = seqvec2->size[0];
      emxEnsureCapacity_real_T(x, j2);
      y_data = x->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)seqvec2->size[0], (blasint)1, (blasint)c->size[1],
                  1.0, &a_data[0], (blasint)seqvec2->size[0], &vec_data[0],
                  (blasint)vec->size[0], 0.0, &y_data[0],
                  (blasint)seqvec2->size[0]);
    }
    xtmp = c_tmp * ((double)(u + 1U) - 1.0);
    nx = x->size[0];
    for (j2 = 0; j2 < nx; j2++) {
      mat_data[j2 + mat->size[0] * (int)u] = (y_data[j2] + xtmp) + 1.0;
    }
  }
  emxFree_real_T(&x);
  emxFree_real_T(&a);
  emxFree_real_T(&seqvec2);
  emxFree_real_T(&vec);
  emxFree_real_T(&seqvec);
  emxFree_real_T(&c);
}

/*
 * function mat = encodekmers_norc(l,k)
 */
static void encodekmers_norc(double l, double k, emxArray_real_T *mat)
{
  emxArray_real_T *a;
  emxArray_real_T *c;
  emxArray_real_T *seqvec;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  emxArray_real_T *y;
  double b_c;
  double c_tmp;
  double *a_data;
  double *c_data;
  double *mat_data;
  double *seqvec_data;
  double *vec_data;
  double *y_data;
  int b_i;
  int b_k;
  int i;
  int i1;
  int nx;
  int u1;
  emxInit_real_T(&c, 2);
  /* 'getgkmcounts:277' c = genIndex(l,k); */
  b_genIndex(l, k, c);
  c_data = c->data;
  /* 'getgkmcounts:278' lcnum = length(c); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    u1 = 0;
  } else {
    nx = c->size[0];
    u1 = c->size[1];
    if (nx >= u1) {
      u1 = nx;
    }
  }
  emxInit_real_T(&seqvec, 2);
  /* 'getgkmcounts:279' seqvec = zeros(4^l, l); */
  c_tmp = pow(4.0, l);
  i = seqvec->size[0] * seqvec->size[1];
  seqvec->size[0] = (int)c_tmp;
  i1 = (int)l;
  seqvec->size[1] = (int)l;
  emxEnsureCapacity_real_T(seqvec, i);
  seqvec_data = seqvec->data;
  nx = (int)c_tmp * (int)l;
  for (i = 0; i < nx; i++) {
    seqvec_data[i] = 0.0;
  }
  /* 'getgkmcounts:280' vec = (1:4^l)'-1; */
  emxInit_real_T(&y, 2);
  y_data = y->data;
  if (c_tmp < 1.0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    nx = (int)floor(c_tmp - 1.0);
    y->size[1] = nx + 1;
    emxEnsureCapacity_real_T(y, i);
    y_data = y->data;
    for (i = 0; i <= nx; i++) {
      y_data[i] = (double)i + 1.0;
    }
  }
  emxInit_real_T(&vec, 1);
  i = vec->size[0];
  vec->size[0] = y->size[1];
  emxEnsureCapacity_real_T(vec, i);
  vec_data = vec->data;
  nx = y->size[1];
  for (i = 0; i < nx; i++) {
    vec_data[i] = y_data[i] - 1.0;
  }
  /* 'getgkmcounts:281' for i = 1:l */
  emxInit_real_T(&x, 1);
  for (b_i = 0; b_i < i1; b_i++) {
    /* 'getgkmcounts:282' seqvec(:,i) = mod(floor(vec/4^(i-1)), 4); */
    b_c = pow(4.0, ((double)b_i + 1.0) - 1.0);
    nx = vec->size[0];
    i = x->size[0];
    x->size[0] = vec->size[0];
    emxEnsureCapacity_real_T(x, i);
    y_data = x->data;
    for (i = 0; i < nx; i++) {
      y_data[i] = vec_data[i] / b_c;
    }
    nx = x->size[0];
    for (b_k = 0; b_k < nx; b_k++) {
      y_data[b_k] = floor(y_data[b_k]);
    }
    nx = x->size[0];
    for (i = 0; i < nx; i++) {
      b_c = y_data[i];
      if (b_c == 0.0) {
        b_c = 0.0;
      } else {
        b_c = fmod(b_c, 4.0);
        if (b_c == 0.0) {
          b_c = 0.0;
        }
      }
      seqvec_data[i + seqvec->size[0] * b_i] = b_c;
    }
  }
  emxFree_real_T(&x);
  /* 'getgkmcounts:284' mat = zeros(4^l, lcnum); */
  i = mat->size[0] * mat->size[1];
  mat->size[0] = (int)c_tmp;
  mat->size[1] = u1;
  emxEnsureCapacity_real_T(mat, i);
  mat_data = mat->data;
  nx = (int)c_tmp * u1;
  for (i = 0; i < nx; i++) {
    mat_data[i] = 0.0;
  }
  /* 'getgkmcounts:285' pow = 4.^(0:k-1)'; */
  if (k - 1.0 < 0.0) {
    y->size[1] = 0;
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    nx = (int)floor(k - 1.0);
    y->size[1] = nx + 1;
    emxEnsureCapacity_real_T(y, i);
    y_data = y->data;
    for (i = 0; i <= nx; i++) {
      y_data[i] = i;
    }
  }
  i = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity_real_T(y, i);
  y_data = y->data;
  nx = y->size[1] - 1;
  for (i = 0; i <= nx; i++) {
    b_c = y_data[i];
    y_data[i] = pow(4.0, b_c);
  }
  /* 'getgkmcounts:286' for i = 1:lcnum */
  emxInit_real_T(&a, 2);
  for (b_i = 0; b_i < u1; b_i++) {
    /* 'getgkmcounts:287' mat(:,i) = seqvec(:,c(i,:))*pow+4^k*(i-1)+1; */
    nx = seqvec->size[0];
    b_k = c->size[1];
    i = a->size[0] * a->size[1];
    a->size[0] = seqvec->size[0];
    a->size[1] = c->size[1];
    emxEnsureCapacity_real_T(a, i);
    a_data = a->data;
    for (i = 0; i < b_k; i++) {
      for (i1 = 0; i1 < nx; i1++) {
        a_data[i1 + a->size[0] * i] =
            seqvec_data[i1 + seqvec->size[0] *
                                 ((int)c_data[b_i + c->size[0] * i] - 1)];
      }
    }
    if ((c->size[1] == 0) || (y->size[1] == 0)) {
      nx = seqvec->size[0];
      i = vec->size[0];
      vec->size[0] = seqvec->size[0];
      emxEnsureCapacity_real_T(vec, i);
      vec_data = vec->data;
      for (i = 0; i < nx; i++) {
        vec_data[i] = 0.0;
      }
    } else {
      i = vec->size[0];
      vec->size[0] = seqvec->size[0];
      emxEnsureCapacity_real_T(vec, i);
      vec_data = vec->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                  (blasint)seqvec->size[0], (blasint)1, (blasint)c->size[1],
                  1.0, &a_data[0], (blasint)seqvec->size[0], &y_data[0],
                  (blasint)1, 0.0, &vec_data[0], (blasint)seqvec->size[0]);
    }
    b_c = pow(4.0, k) * (((double)b_i + 1.0) - 1.0);
    nx = vec->size[0];
    for (i = 0; i < nx; i++) {
      mat_data[i + mat->size[0] * b_i] = (vec_data[i] + b_c) + 1.0;
    }
  }
  emxFree_real_T(&a);
  emxFree_real_T(&y);
  emxFree_real_T(&vec);
  emxFree_real_T(&seqvec);
  emxFree_real_T(&c);
}

static void f_binary_expand_op(emxArray_real_T *gkmc, double a,
                               const emxArray_real_T *temp)
{
  emxArray_real_T *b_gkmc;
  const double *temp_data;
  double *b_gkmc_data;
  double *gkmc_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  temp_data = temp->data;
  gkmc_data = gkmc->data;
  emxInit_real_T(&b_gkmc, 1);
  i = b_gkmc->size[0];
  if (temp->size[0] == 1) {
    b_gkmc->size[0] = gkmc->size[0];
  } else {
    b_gkmc->size[0] = temp->size[0];
  }
  emxEnsureCapacity_real_T(b_gkmc, i);
  b_gkmc_data = b_gkmc->data;
  stride_0_0 = (gkmc->size[0] != 1);
  stride_1_0 = (temp->size[0] != 1);
  if (temp->size[0] == 1) {
    loop_ub = gkmc->size[0];
  } else {
    loop_ub = temp->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_gkmc_data[i] = gkmc_data[i * stride_0_0] + a * temp_data[i * stride_1_0];
  }
  i = gkmc->size[0];
  gkmc->size[0] = b_gkmc->size[0];
  emxEnsureCapacity_real_T(gkmc, i);
  gkmc_data = gkmc->data;
  loop_ub = b_gkmc->size[0];
  for (i = 0; i < loop_ub; i++) {
    gkmc_data[i] = b_gkmc_data[i];
  }
  emxFree_real_T(&b_gkmc);
}

static void g_binary_expand_op(emxArray_real_T *temp, int i1, int i2, int i3)
{
  emxArray_real_T *b_temp;
  double *b_temp_data;
  double *temp_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  temp_data = temp->data;
  emxInit_real_T(&b_temp, 1);
  i = b_temp->size[0];
  if ((i3 - i2) + 1 == 1) {
    b_temp->size[0] = i1 + 1;
  } else {
    b_temp->size[0] = (i3 - i2) + 1;
  }
  emxEnsureCapacity_real_T(b_temp, i);
  b_temp_data = b_temp->data;
  stride_0_0 = (i1 + 1 != 1);
  stride_1_0 = ((i3 - i2) + 1 != 1);
  if ((i3 - i2) + 1 == 1) {
    loop_ub = i1 + 1;
  } else {
    loop_ub = (i3 - i2) + 1;
  }
  for (i = 0; i < loop_ub; i++) {
    b_temp_data[i] = temp_data[i * stride_0_0] + temp_data[i2 + i * stride_1_0];
  }
  i = temp->size[0];
  temp->size[0] = b_temp->size[0];
  emxEnsureCapacity_real_T(temp, i);
  temp_data = temp->data;
  loop_ub = b_temp->size[0];
  for (i = 0; i < loop_ub; i++) {
    temp_data[i] = b_temp_data[i];
  }
  emxFree_real_T(&b_temp);
}

/*
 * function [gkmc, GCpos, GCneg, mat, mat2] = getgkmcounts(filename,l,k, lk,RC)
 */
void getgkmcounts(const emxArray_char_T *filename, double l, double k,
                  double RC, emxArray_real_T *gkmc, double *GCpos,
                  double *GCneg, double mat[16], double mat2[16])
{
  static const char cv1[12] = {'_', 's', 'v', 'a', 'l', 'p',
                               'h', 'a', '.', 'o', 'u', 't'};
  static const char cv2[10] = {'.', 'm', 'o', 'd', 'e',
                               'l', '.', 't', 'x', 't'};
  static const char cv[9] = {'_', 's', 'v', 's', 'e', 'q', '.', 'f', 'a'};
  cell_wrap_0 r;
  cell_wrap_0 r1;
  cell_wrap_10 *sequences_data;
  emxArray_boolean_T *b_x;
  emxArray_boolean_T *c_x;
  emxArray_cell_wrap_10 *sequences;
  emxArray_char_T *b_fileid;
  emxArray_char_T *c_fileid;
  emxArray_char_T *cur_alpha;
  emxArray_char_T *cur_line;
  emxArray_char_T *d_fileid;
  emxArray_char_T *filename_tmp;
  emxArray_char_T *line;
  emxArray_int32_T *r2;
  emxArray_int32_T *r4;
  emxArray_real_T *A;
  emxArray_real_T *Ker;
  emxArray_real_T *a__2;
  emxArray_real_T *a__5;
  emxArray_real_T *alpha;
  emxArray_real_T *c_y;
  emxArray_real_T *comb;
  emxArray_real_T *ss;
  emxArray_real_T *temp;
  emxArray_real_T *x;
  emxArray_uint32_T *len;
  creal_T dc;
  double dv[16];
  double dv2[16];
  double dv1[4];
  double a;
  double alphasum;
  double c_tmp;
  double curr_pos;
  double d;
  double idx;
  double lcnum2;
  double per;
  double rcnum;
  double *A_data;
  double *Ker_data;
  double *alpha_data;
  double *ss_data;
  double *temp_data;
  double *x_data;
  double *y_data;
  int b_i;
  int b_k;
  int b_loop_ub;
  int b_y;
  int c_loop_ub;
  int d_y;
  int e_y;
  int exitg1;
  int i;
  int i1;
  int i2;
  int loop_ub;
  int loop_ub_tmp;
  int nz;
  int u1;
  int ver;
  int vlen;
  unsigned int *len_data;
  int *r3;
  const char *filename_data;
  signed char fileid;
  char *filename_tmp_data;
  char *line_data;
  bool exitg2;
  bool guard1 = false;
  bool y;
  bool *b_x_data;
  filename_data = filename->data;
  emxInit_char_T(&filename_tmp, 2);
  /* gets the gapped kmer counts using the alphas to weight each support vector
   */
  /*  filename is the name of the support vector sequence fa ('_svseq.fa') */
  /*  l and k are the parameters for the length of the gapped kmer (l) and the
   * number of ungapped positions (k) */
  /*  Alternative to isfile */
  /*  if isfile([filename '_svseq.fa']) && isfile([filename '_svalpha.out']) */
  /* 'getgkmcounts:7' if fopen([filename '_svseq.fa']) ~= -1 && fopen([filename
   * '_svalpha.out']) ~= -1 */
  i = filename_tmp->size[0] * filename_tmp->size[1];
  filename_tmp->size[0] = 1;
  filename_tmp->size[1] = filename->size[1] + 9;
  emxEnsureCapacity_char_T(filename_tmp, i);
  filename_tmp_data = filename_tmp->data;
  loop_ub = filename->size[1];
  for (i = 0; i < loop_ub; i++) {
    filename_tmp_data[i] = filename_data[i];
  }
  for (i = 0; i < 9; i++) {
    filename_tmp_data[i + filename->size[1]] = cv[i];
  }
  if (b_strcmp(filename_tmp)) {
    vlen = 0;
  } else {
    fileid = cfopen(filename_tmp, "rb");
    vlen = fileid;
  }
  emxInit_real_T(&alpha, 1);
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&line, 2);
  emxInit_char_T(&cur_alpha, 2);
  emxInit_cell_wrap_10(&sequences);
  emxInitStruct_cell_wrap_0(&r);
  emxInitStruct_cell_wrap_0(&r1);
  emxInit_char_T(&b_fileid, 2);
  emxInit_char_T(&c_fileid, 2);
  emxInit_char_T(&d_fileid, 2);
  guard1 = false;
  if (vlen != -1) {
    i = line->size[0] * line->size[1];
    line->size[0] = 1;
    line->size[1] = filename->size[1] + 12;
    emxEnsureCapacity_char_T(line, i);
    line_data = line->data;
    loop_ub = filename->size[1];
    for (i = 0; i < loop_ub; i++) {
      line_data[i] = filename_data[i];
    }
    for (i = 0; i < 12; i++) {
      line_data[i + filename->size[1]] = cv1[i];
    }
    if (b_strcmp(line)) {
      vlen = 0;
    } else {
      fileid = cfopen(line, "rb");
      vlen = fileid;
    }
    if (vlen != -1) {
      /* 'getgkmcounts:8' filenameseq = [filename '_svseq.fa']; */
      /* 'getgkmcounts:9' filenamealpha = [filename '_svalpha.out']; */
      /*  Alternative to dlmread */
      /*  alpha = dlmread([filename '_svalpha.out'],'\t',0,1); */
      /* 'getgkmcounts:13' fp = fopen(filenamealpha, 'r'); */
      fileid = cfopen(line, "rb");
      /* 'getgkmcounts:14' idx=0; */
      idx = 0.0;
      /* 'getgkmcounts:15' while ~feof(fp) */
      do {
        exitg1 = 0;
        d = b_feof(fileid);
        if (d == 0.0) {
          /* 'getgkmcounts:16' idx=idx+1; */
          idx++;
          /* 'getgkmcounts:17' fgetl(fp); */
          b_fgets(fileid, b_fileid);
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
      /* 'getgkmcounts:19' alpha = zeros(idx, 1); */
      loop_ub_tmp = (int)idx;
      i = alpha->size[0];
      alpha->size[0] = (int)idx;
      emxEnsureCapacity_real_T(alpha, i);
      alpha_data = alpha->data;
      for (i = 0; i < loop_ub_tmp; i++) {
        alpha_data[i] = 0.0;
      }
      /*  pre-allocate alpha vector */
      /* 'getgkmcounts:20' frewind(fp); */
      frewind(fileid);
      /*  Back to the start */
      /* 'getgkmcounts:21' idx=0; */
      idx = 0.0;
      /* 'getgkmcounts:22' while ~feof(fp) */
      do {
        exitg1 = 0;
        d = b_feof(fileid);
        if (d == 0.0) {
          /* 'getgkmcounts:23' idx = idx+1; */
          idx++;
          /* 'getgkmcounts:24' cur_line = fgetl(fp); */
          fgetl(fileid, cur_line);
          /* 'getgkmcounts:25' [~, cur_alpha] = strtok(cur_line, char(9)); */
          b_strtok(cur_line, line, cur_alpha);
          /* 'getgkmcounts:26' alpha(idx, 1) = real(str2double(cur_alpha)); */
          dc = str2double(cur_alpha);
          alpha_data[(int)idx - 1] = dc.re;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
      /* 'getgkmcounts:28' fclose(fp); */
      cfclose(fileid);
      /*  Alternative to importdata */
      /*  sequences = importdata([filename '_svseq.fa']); */
      /* 'getgkmcounts:32' fp = fopen(filenameseq, 'r'); */
      fileid = cfopen(filename_tmp, "rb");
      /* 'getgkmcounts:33' idx=0; */
      idx = 0.0;
      /* 'getgkmcounts:34' while ~feof(fp) */
      do {
        exitg1 = 0;
        d = b_feof(fileid);
        if (d == 0.0) {
          /* 'getgkmcounts:35' idx=idx+1; */
          idx++;
          /* 'getgkmcounts:36' fgetl(fp); */
          b_fgets(fileid, d_fileid);
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
      /* 'getgkmcounts:38' sequences = cell(idx,1); */
      vlen = (int)idx;
      i = sequences->size[0];
      sequences->size[0] = (int)idx;
      emxEnsureCapacity_cell_wrap_10(sequences, i);
      sequences_data = sequences->data;
      for (i = 0; i < vlen; i++) {
        sequences_data[i].f1[0].f1->size[1] = 0;
        sequences_data[i].f1[0].f1->size[0] = 1;
      }
      /* 'getgkmcounts:39' sequences = coder.nullcopy(sequences); */
      /* 'getgkmcounts:40' for cur_idx=1:idx */
      if (0 <= (int)idx - 1) {
        r.f1->size[0] = 1;
        r.f1->size[1] = 0;
      }
      for (loop_ub = 0; loop_ub < vlen; loop_ub++) {
        /* 'getgkmcounts:41' sequences{cur_idx} = cellstr(""); */
        emxCopyStruct_cell_wrap_0(&sequences_data[loop_ub].f1[0], &r);
      }
      /* 'getgkmcounts:43' frewind(fp); */
      frewind(fileid);
      /* 'getgkmcounts:44' for cur_idx=1:idx */
      for (loop_ub = 0; loop_ub < vlen; loop_ub++) {
        /* 'getgkmcounts:45' sequences{cur_idx,1} = cellstr(string(fgetl(fp)));
         */
        fgetl(fileid, r1.f1);
        emxCopyStruct_cell_wrap_0(&sequences_data[loop_ub].f1[0], &r1);
      }
      /* 'getgkmcounts:47' fclose(fp); */
      cfclose(fileid);
      /* 'getgkmcounts:49' ver = 0; */
      ver = 0;
      /*  Alternative to isfile */
      /*  elseif isfile([filename '.model.txt']) */
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }
  if (guard1) {
    i = filename_tmp->size[0] * filename_tmp->size[1];
    filename_tmp->size[0] = 1;
    filename_tmp->size[1] = filename->size[1] + 10;
    emxEnsureCapacity_char_T(filename_tmp, i);
    filename_tmp_data = filename_tmp->data;
    loop_ub = filename->size[1];
    for (i = 0; i < loop_ub; i++) {
      filename_tmp_data[i] = filename_data[i];
    }
    for (i = 0; i < 10; i++) {
      filename_tmp_data[i + filename->size[1]] = cv2[i];
    }
    if (b_strcmp(filename_tmp)) {
      vlen = 0;
    } else {
      fileid = cfopen(filename_tmp, "rb");
      vlen = fileid;
    }
    if (vlen != -1) {
      /* 'getgkmcounts:53' elseif fopen([filename '.model.txt']) ~= -1 */
      /* 'getgkmcounts:54' fid = fopen([filename '.model.txt'], 'r'); */
      fileid = cfopen(filename_tmp, "rb");
      /* 'getgkmcounts:55' a = 1; */
      /* 'getgkmcounts:56' while a == 1 */
      do {
        /* 'getgkmcounts:57' line = fgetl(fid); */
        fgetl(fileid, line);
        /* 'getgkmcounts:58' if strcmp('SV', line) */
      } while (!c_strcmp(line));
      /* 'getgkmcounts:59' a = 0; */
      /*  Alternative to textscan */
      /*  X = textscan(fid,'%f %s\n'); */
      /*  sequences = X{2}; */
      /*  alpha = X{1}; */
      /* 'getgkmcounts:67' curr_pos = ftell(fid); */
      curr_pos = b_ftell(fileid);
      /* 'getgkmcounts:68' idx=0; */
      idx = 0.0;
      /* 'getgkmcounts:69' while ~feof(fid) */
      do {
        exitg1 = 0;
        d = b_feof(fileid);
        if (d == 0.0) {
          /* 'getgkmcounts:70' idx=idx+1; */
          idx++;
          /* 'getgkmcounts:71' fgetl(fid); */
          b_fgets(fileid, c_fileid);
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
      /* 'getgkmcounts:73' fseek(fid, curr_pos, 'bof'); */
      b_fseek(fileid, curr_pos);
      /* 'getgkmcounts:75' alpha = zeros(idx, 1); */
      loop_ub_tmp = (int)idx;
      i = alpha->size[0];
      alpha->size[0] = (int)idx;
      emxEnsureCapacity_real_T(alpha, i);
      alpha_data = alpha->data;
      for (i = 0; i < loop_ub_tmp; i++) {
        alpha_data[i] = 0.0;
      }
      /*  pre-allocate alpha vector */
      /* 'getgkmcounts:76' sequences = cell(idx,1); */
      i = sequences->size[0];
      sequences->size[0] = (int)idx;
      emxEnsureCapacity_cell_wrap_10(sequences, i);
      sequences_data = sequences->data;
      for (i = 0; i < loop_ub_tmp; i++) {
        sequences_data[i].f1[0].f1->size[1] = 0;
        sequences_data[i].f1[0].f1->size[0] = 1;
      }
      /*  pre-allocate seq cell arr */
      /* 'getgkmcounts:77' sequences = coder.nullcopy(sequences); */
      /* 'getgkmcounts:78' for cur_idx=1:idx */
      loop_ub = 0;
      exitg2 = false;
      while ((!exitg2) && (loop_ub <= (int)idx - 1)) {
        /* 'getgkmcounts:79' cur_line = fgetl(fid); */
        fgetl(fileid, cur_line);
        /* 'getgkmcounts:80' if cur_line == -1 */
        y = (cur_line->size[1] != 0);
        if (y) {
          y = (0 > cur_line->size[1] - 1);
        }
        if (y) {
          exitg2 = true;
        } else {
          /* 'getgkmcounts:83' [cur_alpha, cur_seq] = strtok(cur_line, ' '); */
          c_strtok(cur_line, cur_alpha, filename_tmp);
          /* 'getgkmcounts:84' alpha(cur_idx,1) = real(str2double(cur_alpha));
           */
          dc = str2double(cur_alpha);
          alpha_data[loop_ub] = dc.re;
          /* 'getgkmcounts:85' sequences{cur_idx} =
           * cellstr(string(strip(cur_seq))); */
          strip(filename_tmp, r.f1);
          emxCopyStruct_cell_wrap_0(&sequences_data[loop_ub].f1[0], &r);
          loop_ub++;
        }
      }
      /* 'getgkmcounts:87' fclose(fid); */
      cfclose(fileid);
      /* 'getgkmcounts:89' ver = 1; */
      ver = 1;
    } else {
      /* 'getgkmcounts:90' else */
      /* 'getgkmcounts:91' error(['Needs ' filename '_svseq.fa and ' filename
       * '_svalpha.out or ' filename '.model.txt']); */
      /*  Dummy Initialization to Appease Matlab Coder */
      /* 'getgkmcounts:93' idx = 1; */
      /* 'getgkmcounts:94' alpha = zeros(idx,1); */
      i = alpha->size[0];
      alpha->size[0] = 1;
      emxEnsureCapacity_real_T(alpha, i);
      alpha_data = alpha->data;
      alpha_data[0] = 0.0;
      /* 'getgkmcounts:95' sequences = cell(idx,1); */
      /* 'getgkmcounts:96' sequences = coder.nullcopy(sequences); */
      i = sequences->size[0];
      sequences->size[0] = 1;
      emxEnsureCapacity_cell_wrap_10(sequences, i);
      sequences_data = sequences->data;
      /* 'getgkmcounts:97' for cur_idx=1:idx */
      /* 'getgkmcounts:98' sequences{cur_idx} = cellstr(""); */
      r.f1->size[0] = 1;
      r.f1->size[1] = 0;
      emxCopyStruct_cell_wrap_0(&sequences_data[0].f1[0], &r);
      /* 'getgkmcounts:100' ver = 0; */
      ver = 0;
    }
  }
  emxFree_char_T(&d_fileid);
  emxFree_char_T(&c_fileid);
  emxFree_char_T(&b_fileid);
  emxFreeStruct_cell_wrap_0(&r1);
  emxFreeStruct_cell_wrap_0(&r);
  emxFree_char_T(&cur_alpha);
  emxFree_char_T(&cur_line);
  /* 'getgkmcounts:103' if RC */
  emxInit_real_T(&x, 2);
  if (RC != 0.0) {
    /* 'getgkmcounts:104' x = encodekmers(l, k); */
    encodekmers(l, k, x);
    x_data = x->data;
  } else {
    /* 'getgkmcounts:105' else */
    /* 'getgkmcounts:106' x = encodekmers_norc(l, k); */
    encodekmers_norc(l, k, x);
    x_data = x->data;
  }
  emxInit_real_T(&a__2, 2);
  emxInit_real_T(&a__5, 2);
  emxInit_real_T(&Ker, 1);
  emxInit_real_T(&temp, 1);
  emxInit_real_T(&comb, 2);
  /* 'getgkmcounts:109' [comb,~,~,~,~,rcnum] = genIndex(l,k); */
  genIndex(l, k, comb, a__2, Ker, temp, a__5, &rcnum);
  /* 'getgkmcounts:110' lcnum = length(comb); */
  emxFree_real_T(&a__5);
  emxFree_real_T(&a__2);
  if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
    u1 = 0;
  } else {
    vlen = comb->size[0];
    u1 = comb->size[1];
    if (vlen >= u1) {
      u1 = vlen;
    }
  }
  emxFree_real_T(&comb);
  /* 'getgkmcounts:111' n = length(alpha); */
  /*  disp(['# of support vectors: ' num2str(n)]) */
  /* 'getgkmcounts:113' fprintf("# of support vectors: %d\n", int32(n)); */
  printf("# of support vectors: %d\n", alpha->size[0]);
  fflush(stdout);
  /* 'getgkmcounts:115' alphasum = sum(alpha(alpha>0)); */
  loop_ub = alpha->size[0] - 1;
  vlen = 0;
  for (b_i = 0; b_i <= loop_ub; b_i++) {
    if (alpha_data[b_i] > 0.0) {
      vlen++;
    }
  }
  emxInit_int32_T(&r2, 1);
  i = r2->size[0];
  r2->size[0] = vlen;
  emxEnsureCapacity_int32_T(r2, i);
  r3 = r2->data;
  vlen = 0;
  for (b_i = 0; b_i <= loop_ub; b_i++) {
    if (alpha_data[b_i] > 0.0) {
      r3[vlen] = b_i + 1;
      vlen++;
    }
  }
  i = Ker->size[0];
  Ker->size[0] = r2->size[0];
  emxEnsureCapacity_real_T(Ker, i);
  Ker_data = Ker->data;
  loop_ub = r2->size[0];
  for (i = 0; i < loop_ub; i++) {
    Ker_data[i] = alpha_data[r3[i] - 1];
  }
  alphasum = blockedSummation(Ker, r2->size[0]);
  /* 'getgkmcounts:116' gkmc = zeros(4^k*lcnum,1); */
  c_tmp = pow(4.0, k);
  loop_ub_tmp = (int)(c_tmp * (double)u1);
  i = gkmc->size[0];
  gkmc->size[0] = loop_ub_tmp;
  emxEnsureCapacity_real_T(gkmc, i);
  Ker_data = gkmc->data;
  emxFree_int32_T(&r2);
  for (i = 0; i < loop_ub_tmp; i++) {
    Ker_data[i] = 0.0;
  }
  emxInit_boolean_T(&b_x, 1);
  /* 'getgkmcounts:117' GCpos = 0; */
  *GCpos = 0.0;
  /* 'getgkmcounts:118' GCneg = 0; */
  *GCneg = 0.0;
  /* 'getgkmcounts:119' np = sum(alpha>0); */
  i = b_x->size[0];
  b_x->size[0] = alpha->size[0];
  emxEnsureCapacity_boolean_T(b_x, i);
  b_x_data = b_x->data;
  loop_ub = alpha->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_x_data[i] = (alpha_data[i] > 0.0);
  }
  vlen = b_x->size[0];
  if (b_x->size[0] == 0) {
    nz = 0;
  } else {
    nz = b_x_data[0];
    for (b_k = 2; b_k <= vlen; b_k++) {
      nz += b_x_data[b_k - 1];
    }
  }
  vlen = b_x->size[0];
  if (b_x->size[0] == 0) {
    b_y = 0;
  } else {
    b_y = b_x_data[0];
    for (b_k = 2; b_k <= vlen; b_k++) {
      b_y += b_x_data[b_k - 1];
    }
  }
  emxFree_boolean_T(&b_x);
  emxInit_uint32_T(&len);
  /* 'getgkmcounts:120' nn = sum(alpha<0); */
  /* 'getgkmcounts:121' len = zeros(n,1); */
  i = len->size[0];
  len->size[0] = alpha->size[0];
  emxEnsureCapacity_uint32_T(len, i);
  len_data = len->data;
  /* 'getgkmcounts:122' mat = zeros(4); */
  /* 'getgkmcounts:123' mat2 = zeros(4); */
  memset(&mat[0], 0, 16U * sizeof(double));
  memset(&mat2[0], 0, 16U * sizeof(double));
  /* 'getgkmcounts:124' pow = 4.^(0:l-1).'; */
  emxInit_real_T(&c_y, 2);
  if (l - 1.0 < 0.0) {
    c_y->size[1] = 0;
  } else {
    i = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    loop_ub = (int)floor(l - 1.0);
    c_y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(c_y, i);
    y_data = c_y->data;
    for (i = 0; i <= loop_ub; i++) {
      y_data[i] = i;
    }
  }
  i = c_y->size[0] * c_y->size[1];
  c_y->size[0] = 1;
  emxEnsureCapacity_real_T(c_y, i);
  y_data = c_y->data;
  loop_ub = c_y->size[1] - 1;
  for (i = 0; i <= loop_ub; i++) {
    curr_pos = y_data[i];
    y_data[i] = pow(4.0, curr_pos);
  }
  /* 'getgkmcounts:125' n4 = 4^(l-1); */
  idx = pow(4.0, l - 1.0);
  /* 'getgkmcounts:126' slen = numel(sequences); */
  /* 'getgkmcounts:127' per = 10; */
  per = 10.0;
  /* 'getgkmcounts:128' a = 2; */
  a = 2.0;
  /* 'getgkmcounts:129' s = ''; */
  line->size[0] = 1;
  line->size[1] = 0;
  /* 'getgkmcounts:130' l2 = l-1; */
  /* 'getgkmcounts:131' lcnum2 = lcnum-rcnum; */
  lcnum2 = (double)u1 - rcnum;
  /* 'getgkmcounts:132' Ker = zeros(l+1,1); */
  /* 'getgkmcounts:133' M2 = 0; */
  /* 'getgkmcounts:134' for i = 0:(l-k) */
  /* 'getgkmcounts:137' for i = 1:n */
  i = alpha->size[0];
  if (0 <= alpha->size[0] - 1) {
    d_y = (int)floor((double)alpha->size[0] / 10.0);
    if (1.0 > l) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = (int)l;
    }
    i1 = x->size[1];
    c_loop_ub = x->size[1];
  }
  emxInit_real_T(&ss, 2);
  emxInit_int32_T(&r4, 2);
  emxInit_real_T(&A, 2);
  emxInit_boolean_T(&c_x, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'getgkmcounts:138' if mod(i, floor(n/10))==0 */
    vlen = b_i + 1;
    if (d_y != 0) {
      vlen = (int)fmod((double)b_i + 1.0, d_y);
    }
    if (vlen == 0) {
      /* 'getgkmcounts:139' fprintf('%d...', int32(per)); */
      printf("%d...", (int)per);
      fflush(stdout);
      /* 'getgkmcounts:140' per = per+10; */
      per += 10.0;
    }
    /* 'getgkmcounts:142' if ver == 0 */
    if (ver == 0) {
      /* 'getgkmcounts:143' cur_seq = char(string(sequences{a})); */
      i2 = filename_tmp->size[0] * filename_tmp->size[1];
      filename_tmp->size[0] = 1;
      filename_tmp->size[1] = sequences_data[(int)a - 1].f1[0].f1->size[1];
      emxEnsureCapacity_char_T(filename_tmp, i2);
      filename_tmp_data = filename_tmp->data;
      loop_ub = sequences_data[(int)a - 1].f1[0].f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        filename_tmp_data[i2] = sequences_data[(int)a - 1].f1[0].f1->data[i2];
      }
      /* 'getgkmcounts:144' while ~strcmp('>', cur_seq(1)) */
      exitg2 = false;
      while ((!exitg2) && ('>' != filename_tmp_data[0])) {
        /* 'getgkmcounts:145' s = [s cur_seq]; */
        i2 = line->size[1];
        loop_ub = filename_tmp->size[1];
        e_y = line->size[0] * line->size[1];
        line->size[1] += filename_tmp->size[1];
        emxEnsureCapacity_char_T(line, e_y);
        line_data = line->data;
        for (e_y = 0; e_y < loop_ub; e_y++) {
          line_data[i2 + e_y] = filename_tmp_data[e_y];
        }
        /* 'getgkmcounts:146' a = a+1; */
        a++;
        /* 'getgkmcounts:147' if a > slen */
        if (a > sequences->size[0]) {
          exitg2 = true;
        } else {
          /* 'getgkmcounts:150' cur_seq = char(string(sequences{a})); */
          i2 = filename_tmp->size[0] * filename_tmp->size[1];
          filename_tmp->size[0] = 1;
          e_y = (int)a - 1;
          filename_tmp->size[1] = sequences_data[e_y].f1[0].f1->size[1];
          emxEnsureCapacity_char_T(filename_tmp, i2);
          filename_tmp_data = filename_tmp->data;
          loop_ub = sequences_data[e_y].f1[0].f1->size[1];
          for (i2 = 0; i2 < loop_ub; i2++) {
            filename_tmp_data[i2] = sequences_data[e_y].f1[0].f1->data[i2];
          }
        }
      }
    } else {
      /* 'getgkmcounts:152' else */
      /* 'getgkmcounts:153' if ~isempty(sequences{i}) */
      /* 'getgkmcounts:154' s = char(string(sequences{i})); */
      i2 = line->size[0] * line->size[1];
      line->size[0] = 1;
      line->size[1] = sequences_data[b_i].f1[0].f1->size[1];
      emxEnsureCapacity_char_T(line, i2);
      line_data = line->data;
      loop_ub = sequences_data[b_i].f1[0].f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        line_data[i2] = sequences_data[b_i].f1[0].f1->data[i2];
      }
    }
    /* 'getgkmcounts:158' ss = letterconvert(s); */
    letterconvert(line, ss);
    ss_data = ss->data;
    /* 'getgkmcounts:159' if RC */
    if (RC != 0.0) {
      /* 'getgkmcounts:160' temp = zeros(4^k*(lcnum*2), 1); */
      loop_ub = (int)(c_tmp * ((double)u1 * 2.0));
      i2 = temp->size[0];
      temp->size[0] = loop_ub;
      emxEnsureCapacity_real_T(temp, i2);
      temp_data = temp->data;
      for (i2 = 0; i2 < loop_ub; i2++) {
        temp_data[i2] = 0.0;
      }
    } else {
      /* 'getgkmcounts:161' else */
      /* 'getgkmcounts:162' temp = zeros(4^k*(lcnum), 1); */
      i2 = temp->size[0];
      temp->size[0] = loop_ub_tmp;
      emxEnsureCapacity_real_T(temp, i2);
      temp_data = temp->data;
      for (i2 = 0; i2 < loop_ub_tmp; i2++) {
        temp_data[i2] = 0.0;
      }
    }
    /* 'getgkmcounts:164' vec = ss(1:l)*pow; */
    i2 = A->size[0] * A->size[1];
    A->size[0] = 1;
    A->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(A, i2);
    A_data = A->data;
    for (i2 = 0; i2 < b_loop_ub; i2++) {
      A_data[i2] = ss_data[i2];
    }
    if (b_loop_ub < 1) {
      curr_pos = 0.0;
    } else {
      curr_pos = cblas_ddot((blasint)b_loop_ub, &A_data[0], (blasint)1,
                            &y_data[0], (blasint)1);
    }
    /* 'getgkmcounts:165' en = x(vec+1, :); */
    /* 'getgkmcounts:166' temp(en) = temp(en) + 1; */
    i2 = r4->size[0] * r4->size[1];
    r4->size[0] = 1;
    r4->size[1] = i1;
    emxEnsureCapacity_int32_T(r4, i2);
    r3 = r4->data;
    for (i2 = 0; i2 < c_loop_ub; i2++) {
      r3[i2] = (int)x_data[((int)(curr_pos + 1.0) + x->size[0] * i2) - 1];
    }
    i2 = A->size[0] * A->size[1];
    A->size[0] = 1;
    A->size[1] = r4->size[1];
    emxEnsureCapacity_real_T(A, i2);
    A_data = A->data;
    loop_ub = r4->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      A_data[i2] = 1.0;
    }
    loop_ub = A->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      temp_data[r3[i2] - 1] = 1.0;
    }
    /* 'getgkmcounts:167' for j = 2:length(ss)-l2 */
    i2 = (int)(((double)ss->size[1] - (l - 1.0)) + -1.0);
    for (b_k = 0; b_k < i2; b_k++) {
      /* 'getgkmcounts:168' vec = (vec-ss(j-1))/4+n4*ss(j+l2); */
      curr_pos = (curr_pos - ss_data[b_k]) / 4.0 +
                 idx * ss_data[(int)(((double)b_k + 2.0) + (l - 1.0)) - 1];
      /* 'getgkmcounts:169' en = x(vec+1, :); */
      /* 'getgkmcounts:170' temp(en) = temp(en) + 1; */
      loop_ub = x->size[1];
      e_y = r4->size[0] * r4->size[1];
      r4->size[0] = 1;
      r4->size[1] = x->size[1];
      emxEnsureCapacity_int32_T(r4, e_y);
      r3 = r4->data;
      for (e_y = 0; e_y < loop_ub; e_y++) {
        r3[e_y] = (int)x_data[((int)(curr_pos + 1.0) + x->size[0] * e_y) - 1];
      }
      e_y = A->size[0] * A->size[1];
      A->size[0] = 1;
      A->size[1] = r4->size[1];
      emxEnsureCapacity_real_T(A, e_y);
      A_data = A->data;
      loop_ub = r4->size[1];
      for (e_y = 0; e_y < loop_ub; e_y++) {
        A_data[e_y] =
            temp_data[(int)x_data[((int)(curr_pos + 1.0) + x->size[0] * e_y) -
                                  1] -
                      1] +
            1.0;
      }
      loop_ub = A->size[1];
      for (e_y = 0; e_y < loop_ub; e_y++) {
        temp_data[r3[e_y] - 1] = A_data[e_y];
      }
    }
    /* 'getgkmcounts:172' if RC */
    if (RC != 0.0) {
      /* 'getgkmcounts:173' temp =
       * temp(1:4^k*lcnum)+temp(4^k*lcnum+1:4^k*2*lcnum); */
      d = pow(4.0, k) * (double)u1;
      if (1.0 > d) {
        loop_ub = 0;
      } else {
        loop_ub = (int)d;
      }
      curr_pos = c_tmp * 2.0 * (double)u1;
      if (d + 1.0 > curr_pos) {
        i2 = 0;
        e_y = 0;
      } else {
        i2 = (int)(d + 1.0) - 1;
        e_y = (int)curr_pos;
      }
      if (loop_ub == e_y - i2) {
        for (e_y = 0; e_y < loop_ub; e_y++) {
          temp_data[e_y] += temp_data[i2 + e_y];
        }
        i2 = temp->size[0];
        temp->size[0] = loop_ub;
        emxEnsureCapacity_real_T(temp, i2);
        temp_data = temp->data;
      } else {
        g_binary_expand_op(temp, loop_ub - 1, i2, e_y - 1);
        temp_data = temp->data;
      }
      /*          temp = sum(reshape(temp, [], 2), 2); */
    }
    /* 'getgkmcounts:176' if rcnum > 0 && RC */
    if ((rcnum > 0.0) && (RC != 0.0)) {
      /* 'getgkmcounts:177'
       * temp(4^k*lcnum2+1:4^k*lcnum)=temp(4^k*lcnum2+1:4^k*lcnum)/sqrt(2); */
      d = c_tmp * lcnum2 + 1.0;
      curr_pos = pow(4.0, k) * (double)u1;
      if (d > curr_pos) {
        i2 = 1;
        e_y = -1;
        vlen = 0;
      } else {
        i2 = (int)d;
        e_y = (int)d - 2;
        vlen = (int)curr_pos;
      }
      loop_ub = (vlen - e_y) - 1;
      vlen = A->size[0] * A->size[1];
      A->size[0] = 1;
      A->size[1] = loop_ub;
      emxEnsureCapacity_real_T(A, vlen);
      A_data = A->data;
      for (vlen = 0; vlen < loop_ub; vlen++) {
        A_data[vlen] = temp_data[(i2 + vlen) - 1] / 1.4142135623730951;
      }
      loop_ub = A->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        temp_data[(e_y + i2) + 1] = A_data[i2];
      }
    }
    /* 'getgkmcounts:179' if lk == 1 */
    /* 'getgkmcounts:180' gkmc = gkmc + alpha(i)/sqrt(temp'*temp)*temp; */
    if (temp->size[0] < 1) {
      curr_pos = 0.0;
    } else {
      curr_pos = cblas_ddot((blasint)temp->size[0], &temp_data[0], (blasint)1,
                            &temp_data[0], (blasint)1);
    }
    curr_pos = alpha_data[b_i] / sqrt(curr_pos);
    loop_ub = gkmc->size[0];
    if (gkmc->size[0] == temp->size[0]) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        Ker_data[i2] += curr_pos * temp_data[i2];
      }
    } else {
      f_binary_expand_op(gkmc, curr_pos, temp);
      Ker_data = gkmc->data;
    }
    /* 'getgkmcounts:181' len(i) = length(ss); */
    len_data[b_i] = (unsigned int)ss->size[1];
    /* 'getgkmcounts:217' if alpha(i) > 0 */
    if (alpha_data[b_i] > 0.0) {
      /* 'getgkmcounts:218' GCpos = GCpos + alpha(i)*(sum(ss==1)+sum(ss==2)); */
      i2 = c_x->size[0] * c_x->size[1];
      c_x->size[0] = 1;
      c_x->size[1] = ss->size[1];
      emxEnsureCapacity_boolean_T(c_x, i2);
      b_x_data = c_x->data;
      loop_ub = ss->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_x_data[i2] = (ss_data[i2] == 1.0);
      }
      vlen = c_x->size[1];
      if (c_x->size[1] == 0) {
        e_y = 0;
      } else {
        e_y = b_x_data[0];
        for (b_k = 2; b_k <= vlen; b_k++) {
          e_y += b_x_data[b_k - 1];
        }
      }
      i2 = c_x->size[0] * c_x->size[1];
      c_x->size[0] = 1;
      c_x->size[1] = ss->size[1];
      emxEnsureCapacity_boolean_T(c_x, i2);
      b_x_data = c_x->data;
      loop_ub = ss->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_x_data[i2] = (ss_data[i2] == 2.0);
      }
      vlen = c_x->size[1];
      if (c_x->size[1] == 0) {
        loop_ub = 0;
      } else {
        loop_ub = b_x_data[0];
        for (b_k = 2; b_k <= vlen; b_k++) {
          loop_ub += b_x_data[b_k - 1];
        }
      }
      *GCpos += alpha_data[b_i] * ((double)e_y + (double)loop_ub);
      /* 'getgkmcounts:219' ss = ss+1; */
      i2 = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity_real_T(ss, i2);
      ss_data = ss->data;
      loop_ub = ss->size[1] - 1;
      for (i2 = 0; i2 <= loop_ub; i2++) {
        ss_data[i2]++;
      }
      /* 'getgkmcounts:220' for j = 1:length(ss)-1 */
      i2 = ss->size[1];
      for (b_k = 0; b_k <= i2 - 2; b_k++) {
        /* 'getgkmcounts:221' mat2(ss(j),ss(j+1)) = mat2(ss(j),ss(j+1)) +
         * abs(alpha(i)); */
        vlen = ((int)ss_data[b_k] + (((int)ss_data[b_k + 1] - 1) << 2)) - 1;
        mat2[vlen] += fabs(alpha_data[b_i]);
      }
    } else {
      /* 'getgkmcounts:223' else */
      /* 'getgkmcounts:224' GCneg = GCneg + alpha(i)*(sum(ss==1)+sum(ss==2)); */
      i2 = c_x->size[0] * c_x->size[1];
      c_x->size[0] = 1;
      c_x->size[1] = ss->size[1];
      emxEnsureCapacity_boolean_T(c_x, i2);
      b_x_data = c_x->data;
      loop_ub = ss->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_x_data[i2] = (ss_data[i2] == 1.0);
      }
      vlen = c_x->size[1];
      if (c_x->size[1] == 0) {
        e_y = 0;
      } else {
        e_y = b_x_data[0];
        for (b_k = 2; b_k <= vlen; b_k++) {
          e_y += b_x_data[b_k - 1];
        }
      }
      i2 = c_x->size[0] * c_x->size[1];
      c_x->size[0] = 1;
      c_x->size[1] = ss->size[1];
      emxEnsureCapacity_boolean_T(c_x, i2);
      b_x_data = c_x->data;
      loop_ub = ss->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_x_data[i2] = (ss_data[i2] == 2.0);
      }
      vlen = c_x->size[1];
      if (c_x->size[1] == 0) {
        loop_ub = 0;
      } else {
        loop_ub = b_x_data[0];
        for (b_k = 2; b_k <= vlen; b_k++) {
          loop_ub += b_x_data[b_k - 1];
        }
      }
      *GCneg += alpha_data[b_i] * ((double)e_y + (double)loop_ub);
      /* 'getgkmcounts:225' ss = ss+1; */
      i2 = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity_real_T(ss, i2);
      ss_data = ss->data;
      loop_ub = ss->size[1] - 1;
      for (i2 = 0; i2 <= loop_ub; i2++) {
        ss_data[i2]++;
      }
      /* 'getgkmcounts:226' for j = 1:length(ss)-1 */
      i2 = ss->size[1];
      for (b_k = 0; b_k <= i2 - 2; b_k++) {
        /* 'getgkmcounts:227' mat(ss(j),ss(j+1)) = mat(ss(j),ss(j+1)) +
         * abs(alpha(i)); */
        vlen = ((int)ss_data[b_k] + (((int)ss_data[b_k + 1] - 1) << 2)) - 1;
        mat[vlen] += fabs(alpha_data[b_i]);
      }
    }
    /* 'getgkmcounts:230' s = ''; */
    line->size[0] = 1;
    line->size[1] = 0;
    /* 'getgkmcounts:231' a = a+1; */
    a++;
  }
  emxFree_boolean_T(&c_x);
  emxFree_real_T(&A);
  emxFree_char_T(&filename_tmp);
  emxFree_real_T(&c_y);
  emxFree_int32_T(&r4);
  emxFree_cell_wrap_10(&sequences);
  emxFree_real_T(&temp);
  emxFree_real_T(&ss);
  emxFree_real_T(&x);
  emxFree_char_T(&line);
  emxFree_real_T(&alpha);
  /* 'getgkmcounts:233' fprintf('\n') */
  printf("\n");
  fflush(stdout);
  /* 'getgkmcounts:234' GCpos = GCpos/alphasum/mean(len(1:np)); */
  if (1 > b_y) {
    loop_ub = 0;
  } else {
    loop_ub = nz;
  }
  i = Ker->size[0];
  Ker->size[0] = loop_ub;
  emxEnsureCapacity_real_T(Ker, i);
  Ker_data = Ker->data;
  for (i = 0; i < loop_ub; i++) {
    Ker_data[i] = len_data[i];
  }
  *GCpos =
      *GCpos / alphasum / (blockedSummation(Ker, loop_ub) / (double)loop_ub);
  /* 'getgkmcounts:235' GCneg = -1*GCneg/alphasum/mean(len(np+1:end)); */
  if ((double)nz + 1.0 > len->size[0]) {
    i = -1;
    i1 = -1;
  } else {
    i = b_y - 1;
    i1 = len->size[0] - 1;
  }
  loop_ub = i1 - i;
  i1 = Ker->size[0];
  Ker->size[0] = loop_ub;
  emxEnsureCapacity_real_T(Ker, i1);
  Ker_data = Ker->data;
  for (i1 = 0; i1 < loop_ub; i1++) {
    Ker_data[i1] = len_data[(i + i1) + 1];
  }
  emxFree_uint32_T(&len);
  *GCneg =
      -*GCneg / alphasum / (blockedSummation(Ker, loop_ub) / (double)loop_ub);
  /* 'getgkmcounts:236' if RC */
  emxFree_real_T(&Ker);
  if (RC != 0.0) {
    /* 'getgkmcounts:237' mat = (mat+rot90(rot90(mat,3)'))/2; */
    rot90(mat, dv);
    for (i = 0; i < 4; i++) {
      i1 = i << 2;
      dv2[i1] = dv[i];
      dv2[i1 + 1] = dv[i + 4];
      dv2[i1 + 2] = dv[i + 8];
      dv2[i1 + 3] = dv[i + 12];
    }
    b_rot90(dv2, dv);
    for (i = 0; i < 16; i++) {
      mat[i] = (mat[i] + dv[i]) / 2.0;
    }
    /* 'getgkmcounts:238' mat2 = (mat2+rot90(rot90(mat2,3)'))/2; */
    rot90(mat2, dv);
    for (i = 0; i < 4; i++) {
      i1 = i << 2;
      dv2[i1] = dv[i];
      dv2[i1 + 1] = dv[i + 4];
      dv2[i1 + 2] = dv[i + 8];
      dv2[i1 + 3] = dv[i + 12];
    }
    b_rot90(dv2, dv);
    for (i = 0; i < 16; i++) {
      mat2[i] = (mat2[i] + dv[i]) / 2.0;
    }
  }
  /* 'getgkmcounts:240' mat = mat./repmat(sum(mat,2),1,4); */
  c_sum(mat, dv1);
  repmat(dv1, dv);
  for (i = 0; i < 16; i++) {
    mat[i] /= dv[i];
  }
  /* 'getgkmcounts:241' mat2 = mat2./repmat(sum(mat2,2),1,4); */
  c_sum(mat2, dv1);
  repmat(dv1, dv);
  for (i = 0; i < 16; i++) {
    mat2[i] /= dv[i];
  }
}

/*
 * function en = letterconvert(s)
 */
void letterconvert(const emxArray_char_T *s, emxArray_real_T *en)
{
  double *en_data;
  int b_i;
  int i;
  const char *s_data;
  char c;
  s_data = s->data;
  /* 'getgkmcounts:246' l = length(s); */
  /* 'getgkmcounts:247' en = zeros(1,l); */
  i = en->size[0] * en->size[1];
  en->size[0] = 1;
  en->size[1] = s->size[1];
  emxEnsureCapacity_real_T(en, i);
  en_data = en->data;
  /* 'getgkmcounts:248' for i = 1:l */
  i = s->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'getgkmcounts:249' if strcmp(s(i),'A') || strcmp(s(i), 'a') */
    c = s_data[b_i];
    if ((!(c != 'A')) || (!(c != 'a'))) {
      /* 'getgkmcounts:250' en(i) = 0; */
      en_data[b_i] = 0.0;
    } else if ((!(c != 'C')) || (!(c != 'c'))) {
      /* 'getgkmcounts:251' elseif strcmp(s(i),'C') || strcmp(s(i),'c') */
      /* 'getgkmcounts:252' en(i) = 1; */
      en_data[b_i] = 1.0;
    } else if ((!(c != 'G')) || (!(c != 'g'))) {
      /* 'getgkmcounts:253' elseif strcmp(s(i),'G') || strcmp(s(i),'g') */
      /* 'getgkmcounts:254' en(i) = 2; */
      en_data[b_i] = 2.0;
    } else {
      /* 'getgkmcounts:255' else */
      /* 'getgkmcounts:256' en(i) = 3; */
      en_data[b_i] = 3.0;
    }
  }
}

/* End of code generation (getgkmcounts.c) */
