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
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include "mod.h"
#include "mtimes.h"
#include "nchoosek.h"
#include "print_processing.h"
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
static void encodekmers(double l, double k, const emxArray_real_T *c,
                        emxArray_real_T *mat);

static void encodekmers_norc(double l, double k, const emxArray_real_T *c,
                             emxArray_real_T *mat);

static void g_binary_expand_op(emxArray_real_T *gkmc, double a,
                               const emxArray_real_T *temp);

static void h_binary_expand_op(emxArray_real_T *y, const emxArray_real_T *r3,
                               const emxArray_real_T *M);

static void i_binary_expand_op(emxArray_real_T *temp, int i1, int i2, int i3);

static void letterconvert(const emxArray_char_T *s, emxArray_real_T *en);

static void plus(emxArray_real_T *M2, const emxArray_real_T *r1);

/* Function Definitions */
/*
 * function mat = encodekmers(l,k,c)
 */
static void encodekmers(double l, double k, const emxArray_real_T *c,
                        emxArray_real_T *mat)
{
  emxArray_real_T *a;
  emxArray_real_T *seqvec;
  emxArray_real_T *seqvec2;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  emxArray_real_T *y;
  const double *c_data;
  double c_tmp;
  double xtmp;
  double *a_data;
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
  c_data = c->data;
  /* 'getgkmcounts:270' lcnum = length(c); */
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
  /* 'getgkmcounts:271' seqvec = zeros(4^l, l); */
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
  /* 'getgkmcounts:272' vec = (1:4^l)'-1; */
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
  /* 'getgkmcounts:273' for i = 1:l */
  emxInit_real_T(&x, 1);
  for (i = 0; i < b_j1; i++) {
    /* 'getgkmcounts:274' seqvec(:,i) = mod(floor(vec/4^(i-1)), 4); */
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
      seqvec_data[j2 + seqvec->size[0] * i] = b_mod(xtmp, 4.0);
    }
  }
  emxInit_real_T(&seqvec2, 2);
  /* 'getgkmcounts:276' seqvec2 = 3-fliplr(seqvec); */
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
  /* 'getgkmcounts:277' mat = zeros(4^l, 2*lcnum); */
  j2 = mat->size[0] * mat->size[1];
  mat->size[0] = (int)c_tmp;
  mat->size[1] = 2 * u1;
  emxEnsureCapacity_real_T(mat, j2);
  mat_data = mat->data;
  nx = (int)c_tmp * (2 * u1);
  for (j2 = 0; j2 < nx; j2++) {
    mat_data[j2] = 0.0;
  }
  /* 'getgkmcounts:278' pow = 4.^(0:k-1)'; */
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
  /* 'getgkmcounts:279' for i = 1:lcnum */
  emxInit_real_T(&a, 2);
  for (i = 0; i < u1; i++) {
    /* 'getgkmcounts:280' mat(:,i) = seqvec(:,c(i,:))*pow+4^k*(i-1)+1; */
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
    /* 'getgkmcounts:281' mat(:,i+lcnum) =
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
}

/*
 * function mat = encodekmers_norc(l,k,c)
 */
static void encodekmers_norc(double l, double k, const emxArray_real_T *c,
                             emxArray_real_T *mat)
{
  emxArray_real_T *a;
  emxArray_real_T *seqvec;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  emxArray_real_T *y;
  const double *c_data;
  double b_c;
  double c_tmp;
  double *a_data;
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
  c_data = c->data;
  /* 'getgkmcounts:285' lcnum = length(c); */
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
  /* 'getgkmcounts:286' seqvec = zeros(4^l, l); */
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
  /* 'getgkmcounts:287' vec = (1:4^l)'-1; */
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
  /* 'getgkmcounts:288' for i = 1:l */
  emxInit_real_T(&x, 1);
  for (b_i = 0; b_i < i1; b_i++) {
    /* 'getgkmcounts:289' seqvec(:,i) = mod(floor(vec/4^(i-1)), 4); */
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
      seqvec_data[i + seqvec->size[0] * b_i] = b_mod(b_c, 4.0);
    }
  }
  emxFree_real_T(&x);
  /* 'getgkmcounts:291' mat = zeros(4^l, lcnum); */
  i = mat->size[0] * mat->size[1];
  mat->size[0] = (int)c_tmp;
  mat->size[1] = u1;
  emxEnsureCapacity_real_T(mat, i);
  mat_data = mat->data;
  nx = (int)c_tmp * u1;
  for (i = 0; i < nx; i++) {
    mat_data[i] = 0.0;
  }
  /* 'getgkmcounts:292' pow = 4.^(0:k-1)'; */
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
  /* 'getgkmcounts:293' for i = 1:lcnum */
  emxInit_real_T(&a, 2);
  for (b_i = 0; b_i < u1; b_i++) {
    /* 'getgkmcounts:294' mat(:,i) = seqvec(:,c(i,:))*pow+4^k*(i-1)+1; */
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
}

static void g_binary_expand_op(emxArray_real_T *gkmc, double a,
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

static void h_binary_expand_op(emxArray_real_T *y, const emxArray_real_T *r3,
                               const emxArray_real_T *M)
{
  emxArray_real_T *r1;
  const double *M_data;
  const double *b_r;
  double *r2;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  M_data = M->data;
  b_r = r3->data;
  emxInit_real_T(&r1, 2);
  i = r1->size[0] * r1->size[1];
  if (M->size[0] == 1) {
    r1->size[0] = r3->size[0];
  } else {
    r1->size[0] = M->size[0];
  }
  if (M->size[1] == 1) {
    r1->size[1] = r3->size[1];
  } else {
    r1->size[1] = M->size[1];
  }
  emxEnsureCapacity_real_T(r1, i);
  r2 = r1->data;
  stride_0_0 = (r3->size[0] != 1);
  stride_0_1 = (r3->size[1] != 1);
  stride_1_0 = (M->size[0] != 1);
  stride_1_1 = (M->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (M->size[1] == 1) {
    loop_ub = r3->size[1];
  } else {
    loop_ub = M->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    if (M->size[0] == 1) {
      b_loop_ub = r3->size[0];
    } else {
      b_loop_ub = M->size[0];
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      r2[i1 + r1->size[0] * i] = b_r[i1 * stride_0_0 + r3->size[0] * aux_0_1] +
                                 M_data[i1 * stride_1_0 + M->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  c_sum(r1, y);
  emxFree_real_T(&r1);
}

static void i_binary_expand_op(emxArray_real_T *temp, int i1, int i2, int i3)
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
 * function en = letterconvert(s)
 */
static void letterconvert(const emxArray_char_T *s, emxArray_real_T *en)
{
  double *en_data;
  int b_i;
  int i;
  const char *s_data;
  char c;
  s_data = s->data;
  /* 'getgkmcounts:255' l = length(s); */
  /* 'getgkmcounts:256' en = zeros(1,l); */
  i = en->size[0] * en->size[1];
  en->size[0] = 1;
  en->size[1] = s->size[1];
  emxEnsureCapacity_real_T(en, i);
  en_data = en->data;
  /* 'getgkmcounts:257' for i = 1:l */
  i = s->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'getgkmcounts:258' if strcmp(s(i),'A') || strcmp(s(i), 'a') */
    c = s_data[b_i];
    if ((!(c != 'A')) || (!(c != 'a'))) {
      /* 'getgkmcounts:259' en(i) = 0; */
      en_data[b_i] = 0.0;
    } else if ((!(c != 'C')) || (!(c != 'c'))) {
      /* 'getgkmcounts:260' elseif strcmp(s(i),'C') || strcmp(s(i),'c') */
      /* 'getgkmcounts:261' en(i) = 1; */
      en_data[b_i] = 1.0;
    } else if ((!(c != 'G')) || (!(c != 'g'))) {
      /* 'getgkmcounts:262' elseif strcmp(s(i),'G') || strcmp(s(i),'g') */
      /* 'getgkmcounts:263' en(i) = 2; */
      en_data[b_i] = 2.0;
    } else {
      /* 'getgkmcounts:264' else */
      /* 'getgkmcounts:265' en(i) = 3; */
      en_data[b_i] = 3.0;
    }
  }
}

static void plus(emxArray_real_T *M2, const emxArray_real_T *r1)
{
  emxArray_real_T *b_M2;
  const double *b_r;
  double *M2_data;
  double *b_M2_data;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  b_r = r1->data;
  M2_data = M2->data;
  emxInit_real_T(&b_M2, 2);
  i = b_M2->size[0] * b_M2->size[1];
  if (r1->size[0] == 1) {
    b_M2->size[0] = M2->size[0];
  } else {
    b_M2->size[0] = r1->size[0];
  }
  if (r1->size[1] == 1) {
    b_M2->size[1] = M2->size[1];
  } else {
    b_M2->size[1] = r1->size[1];
  }
  emxEnsureCapacity_real_T(b_M2, i);
  b_M2_data = b_M2->data;
  stride_0_0 = (M2->size[0] != 1);
  stride_0_1 = (M2->size[1] != 1);
  stride_1_0 = (r1->size[0] != 1);
  stride_1_1 = (r1->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (r1->size[1] == 1) {
    loop_ub = M2->size[1];
  } else {
    loop_ub = r1->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    if (r1->size[0] == 1) {
      b_loop_ub = M2->size[0];
    } else {
      b_loop_ub = r1->size[0];
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_M2_data[i1 + b_M2->size[0] * i] =
          M2_data[i1 * stride_0_0 + M2->size[0] * aux_0_1] +
          b_r[i1 * stride_1_0 + r1->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = M2->size[0] * M2->size[1];
  M2->size[0] = b_M2->size[0];
  M2->size[1] = b_M2->size[1];
  emxEnsureCapacity_real_T(M2, i);
  M2_data = M2->data;
  loop_ub = b_M2->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = b_M2->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      M2_data[i1 + M2->size[0] * i] = b_M2_data[i1 + b_M2->size[0] * i];
    }
  }
  emxFree_real_T(&b_M2);
}

/*
 * function [gkmc, GCpos, GCneg, mat, mat2] = getgkmcounts(filename,l,k, lk,RC,
 * comb, rcnum)
 */
void getgkmcounts(const emxArray_char_T *filename, double l, double k,
                  const double lk_data[], const int lk_size[2], bool RC,
                  const emxArray_real_T *comb, double rcnum,
                  emxArray_real_T *gkmc, double *GCpos, double *GCneg,
                  double mat[16], double mat2[16])
{
  static const char cv1[12] = {'_', 's', 'v', 'a', 'l', 'p',
                               'h', 'a', '.', 'o', 'u', 't'};
  static const char cv2[10] = {'.', 'm', 'o', 'd', 'e',
                               'l', '.', 't', 'x', 't'};
  static const char cv[9] = {'_', 's', 'v', 's', 'e', 'q', '.', 'f', 'a'};
  cell_wrap_0 C[8];
  cell_wrap_8 validatedHoleFilling[3];
  cell_wrap_8 b_r;
  cell_wrap_8 r1;
  cell_wrap_9 *sequences_data;
  emxArray_boolean_T *b_x;
  emxArray_boolean_T *c_x;
  emxArray_cell_wrap_9 *sequences;
  emxArray_char_T *b_fileid;
  emxArray_char_T *c_fileid;
  emxArray_char_T *cur_alpha;
  emxArray_char_T *d_fileid;
  emxArray_int32_T *r2;
  emxArray_int32_T *r6;
  emxArray_real_T *Ker;
  emxArray_real_T *M;
  emxArray_real_T *M2;
  emxArray_real_T *alpha;
  emxArray_real_T *b_Ker;
  emxArray_real_T *c_y;
  emxArray_real_T *r4;
  emxArray_real_T *r7;
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
  double ll;
  double loop_ub_tmp_tmp;
  double per;
  double *Ker_data;
  double *M2_data;
  double *M_data;
  double *alpha_data;
  double *b_Ker_data;
  double *gkmc_data;
  double *r5;
  double *r8;
  double *ss_data;
  double *temp_data;
  double *x_data;
  int LK;
  int b_i;
  int b_k;
  int b_loop_ub;
  int b_y;
  int c_loop_ub;
  int exitg1;
  int i;
  int i1;
  int i2;
  int ii;
  int j;
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
  bool exitg2;
  bool guard1 = false;
  bool y;
  bool *b_x_data;
  filename_data = filename->data;
  emxInitMatrix_cell_wrap_81(validatedHoleFilling);
  /* gets the gapped kmer counts using the alphas to weight each support vector
   */
  /*  filename is the name of the support vector sequence fa ('_svseq.fa') */
  /*  l and k are the parameters for the length of the gapped kmer (l) and the
   * number of ungapped positions (k) */
  /*  Alternative to isfile */
  /*  if isfile([filename '_svseq.fa']) && isfile([filename '_svalpha.out']) */
  /* 'getgkmcounts:7' if fopen([filename '_svseq.fa']) ~= -1 && fopen([filename
   * '_svalpha.out']) ~= -1 */
  i = validatedHoleFilling[1].f1->size[0] * validatedHoleFilling[1].f1->size[1];
  validatedHoleFilling[1].f1->size[0] = 1;
  validatedHoleFilling[1].f1->size[1] = filename->size[1] + 9;
  emxEnsureCapacity_char_T(validatedHoleFilling[1].f1, i);
  loop_ub = filename->size[1];
  for (i = 0; i < loop_ub; i++) {
    validatedHoleFilling[1].f1->data[i] = filename_data[i];
  }
  for (i = 0; i < 9; i++) {
    validatedHoleFilling[1].f1->data[i + filename->size[1]] = cv[i];
  }
  if (b_strcmp(validatedHoleFilling[1].f1)) {
    vlen = 0;
  } else {
    fileid = cfopen(validatedHoleFilling[1].f1, "rb");
    vlen = fileid;
  }
  emxInit_real_T(&alpha, 1);
  emxInit_char_T(&cur_alpha, 2);
  emxInit_cell_wrap_9(&sequences);
  emxInitStruct_cell_wrap_8(&b_r);
  emxInitStruct_cell_wrap_8(&r1);
  emxInit_char_T(&b_fileid, 2);
  emxInit_char_T(&c_fileid, 2);
  emxInit_char_T(&d_fileid, 2);
  guard1 = false;
  if (vlen != -1) {
    i = validatedHoleFilling[0].f1->size[0] *
        validatedHoleFilling[0].f1->size[1];
    validatedHoleFilling[0].f1->size[0] = 1;
    validatedHoleFilling[0].f1->size[1] = filename->size[1] + 12;
    emxEnsureCapacity_char_T(validatedHoleFilling[0].f1, i);
    loop_ub = filename->size[1];
    for (i = 0; i < loop_ub; i++) {
      validatedHoleFilling[0].f1->data[i] = filename_data[i];
    }
    for (i = 0; i < 12; i++) {
      validatedHoleFilling[0].f1->data[i + filename->size[1]] = cv1[i];
    }
    if (b_strcmp(validatedHoleFilling[0].f1)) {
      vlen = 0;
    } else {
      fileid = cfopen(validatedHoleFilling[0].f1, "rb");
      vlen = fileid;
    }
    if (vlen != -1) {
      /* 'getgkmcounts:8' filenameseq = [filename '_svseq.fa']; */
      /* 'getgkmcounts:9' filenamealpha = [filename '_svalpha.out']; */
      /*  Alternative to dlmread */
      /*  alpha = dlmread([filename '_svalpha.out'],'\t',0,1); */
      /* 'getgkmcounts:13' fp = fopen(filenamealpha, 'r'); */
      fileid = cfopen(validatedHoleFilling[0].f1, "rb");
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
          fgetl(fileid, validatedHoleFilling[2].f1);
          /* 'getgkmcounts:25' [~, cur_alpha] = strtok(cur_line, char(9)); */
          b_strtok(validatedHoleFilling[2].f1, validatedHoleFilling[0].f1,
                   cur_alpha);
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
      fileid = cfopen(validatedHoleFilling[1].f1, "rb");
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
      emxEnsureCapacity_cell_wrap_9(sequences, i);
      sequences_data = sequences->data;
      for (i = 0; i < vlen; i++) {
        sequences_data[i].f1[0].f1->size[1] = 0;
        sequences_data[i].f1[0].f1->size[0] = 1;
      }
      /* 'getgkmcounts:39' sequences = coder.nullcopy(sequences); */
      /* 'getgkmcounts:40' for cur_idx=1:idx */
      if (0 <= (int)idx - 1) {
        b_r.f1->size[0] = 1;
        b_r.f1->size[1] = 0;
      }
      for (loop_ub = 0; loop_ub < vlen; loop_ub++) {
        /* 'getgkmcounts:41' sequences{cur_idx} = cellstr(""); */
        emxCopyStruct_cell_wrap_8(&sequences_data[loop_ub].f1[0], &b_r);
      }
      /* 'getgkmcounts:43' frewind(fp); */
      frewind(fileid);
      /* 'getgkmcounts:44' for cur_idx=1:idx */
      for (loop_ub = 0; loop_ub < vlen; loop_ub++) {
        /* 'getgkmcounts:45' sequences{cur_idx,1} = cellstr(string(fgetl(fp)));
         */
        fgetl(fileid, r1.f1);
        emxCopyStruct_cell_wrap_8(&sequences_data[loop_ub].f1[0], &r1);
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
    i = validatedHoleFilling[1].f1->size[0] *
        validatedHoleFilling[1].f1->size[1];
    validatedHoleFilling[1].f1->size[0] = 1;
    validatedHoleFilling[1].f1->size[1] = filename->size[1] + 10;
    emxEnsureCapacity_char_T(validatedHoleFilling[1].f1, i);
    loop_ub = filename->size[1];
    for (i = 0; i < loop_ub; i++) {
      validatedHoleFilling[1].f1->data[i] = filename_data[i];
    }
    for (i = 0; i < 10; i++) {
      validatedHoleFilling[1].f1->data[i + filename->size[1]] = cv2[i];
    }
    if (b_strcmp(validatedHoleFilling[1].f1)) {
      vlen = 0;
    } else {
      fileid = cfopen(validatedHoleFilling[1].f1, "rb");
      vlen = fileid;
    }
    if (vlen != -1) {
      /* 'getgkmcounts:53' elseif fopen([filename '.model.txt']) ~= -1 */
      /* 'getgkmcounts:54' fid = fopen([filename '.model.txt'], 'r'); */
      fileid = cfopen(validatedHoleFilling[1].f1, "rb");
      /* 'getgkmcounts:55' a = 1; */
      /* 'getgkmcounts:56' while a == 1 */
      do {
        /* 'getgkmcounts:57' line = fgetl(fid); */
        fgetl(fileid, validatedHoleFilling[0].f1);
        /* 'getgkmcounts:58' if strcmp('SV', line) */
      } while (!c_strcmp(validatedHoleFilling[0].f1));
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
      emxEnsureCapacity_cell_wrap_9(sequences, i);
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
        fgetl(fileid, validatedHoleFilling[2].f1);
        /* 'getgkmcounts:80' if cur_line == -1 */
        y = (validatedHoleFilling[2].f1->size[1] != 0);
        if (y) {
          y = (0 > validatedHoleFilling[2].f1->size[1] - 1);
        }
        if (y) {
          exitg2 = true;
        } else {
          /* 'getgkmcounts:83' [cur_alpha, cur_seq] = strtok(cur_line, ' '); */
          c_strtok(validatedHoleFilling[2].f1, cur_alpha,
                   validatedHoleFilling[1].f1);
          /* 'getgkmcounts:84' alpha(cur_idx,1) = real(str2double(cur_alpha));
           */
          dc = str2double(cur_alpha);
          alpha_data[loop_ub] = dc.re;
          /* 'getgkmcounts:85' sequences{cur_idx} =
           * cellstr(string(strip(cur_seq))); */
          strip(validatedHoleFilling[1].f1, b_r.f1);
          emxCopyStruct_cell_wrap_8(&sequences_data[loop_ub].f1[0], &b_r);
          loop_ub++;
        }
      }
      /* 'getgkmcounts:87' fclose(fid); */
      cfclose(fileid);
      /* 'getgkmcounts:89' ver = 1; */
      ver = 1;
    } else {
      /* 'getgkmcounts:90' else */
      /* 'getgkmcounts:91' fprintf('ERROR: Cannot open %s_svseq.fa and
       * %s_svalpha.out or %s.model.txt\n', filename, filename, filename); */
      print_processing(filename, filename, filename, validatedHoleFilling);
      printf(
          "ERROR: Cannot open %s_svseq.fa and %s_svalpha.out or %s.model.txt\n",
          &validatedHoleFilling[0].f1->data[0],
          &validatedHoleFilling[1].f1->data[0],
          &validatedHoleFilling[2].f1->data[0]);
      fflush(stdout);
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
      emxEnsureCapacity_cell_wrap_9(sequences, i);
      sequences_data = sequences->data;
      /* 'getgkmcounts:97' for cur_idx=1:idx */
      /* 'getgkmcounts:98' sequences{cur_idx} = cellstr(""); */
      b_r.f1->size[0] = 1;
      b_r.f1->size[1] = 0;
      emxCopyStruct_cell_wrap_8(&sequences_data[0].f1[0], &b_r);
      /* 'getgkmcounts:100' ver = 0; */
      ver = 0;
    }
  }
  emxFree_char_T(&d_fileid);
  emxFree_char_T(&c_fileid);
  emxFree_char_T(&b_fileid);
  emxFreeStruct_cell_wrap_8(&r1);
  emxFreeStruct_cell_wrap_8(&b_r);
  emxFree_char_T(&cur_alpha);
  /* 'getgkmcounts:103' if RC */
  emxInit_real_T(&x, 2);
  if (RC) {
    /* 'getgkmcounts:104' x = encodekmers(l, k, comb); */
    encodekmers(l, k, comb, x);
    x_data = x->data;
  } else {
    /* 'getgkmcounts:105' else */
    /* 'getgkmcounts:106' x = encodekmers_norc(l, k, comb); */
    encodekmers_norc(l, k, comb, x);
    x_data = x->data;
  }
  emxInit_real_T(&Ker, 1);
  /* 'getgkmcounts:108' ll = -1; */
  ll = -1.0;
  /* 'getgkmcounts:109' kk = -1; */
  /* 'getgkmcounts:110' Ker = 0; */
  i = Ker->size[0];
  Ker->size[0] = 1;
  emxEnsureCapacity_real_T(Ker, i);
  Ker_data = Ker->data;
  Ker_data[0] = 0.0;
  /* 'getgkmcounts:111' if length(lk) > 1 */
  if (lk_size[1] > 1) {
    /* 'getgkmcounts:112' fprintf('Using mismatches to calculate norms\n'); */
    printf("Using mismatches to calculate norms\n");
    fflush(stdout);
    /* 'getgkmcounts:113' ll = lk(1); */
    ll = lk_data[0];
    /* 'getgkmcounts:114' kk = lk(2); */
    /* 'getgkmcounts:115' LK = 0; */
    LK = 0;
    /* 'getgkmcounts:116' Ker = zeros(ll+1,1); */
    loop_ub = (int)(lk_data[0] + 1.0);
    i = Ker->size[0];
    Ker->size[0] = loop_ub;
    emxEnsureCapacity_real_T(Ker, i);
    Ker_data = Ker->data;
    for (i = 0; i < loop_ub; i++) {
      Ker_data[i] = 0.0;
    }
    /* 'getgkmcounts:117' for i = 0:(ll-kk) */
    d = lk_data[0] - lk_data[1];
    i = (int)(d + 1.0);
    for (b_i = 0; b_i < i; b_i++) {
      /* 'getgkmcounts:118' Ker(ll-i+1) = nchoosek(ll-i,ll-kk-i); */
      curr_pos = ll - (double)b_i;
      Ker_data[(int)(curr_pos + 1.0) - 1] =
          b_nchoosek(curr_pos, d - (double)b_i);
    }
  } else {
    /* 'getgkmcounts:120' else */
    /* 'getgkmcounts:121' fprintf('Using gapped k-mer count vector to calculate
     * norms\n'); */
    printf("Using gapped k-mer count vector to calculate norms\n");
    fflush(stdout);
    /* 'getgkmcounts:122' LK = 1; */
    LK = 1;
  }
  /* 'getgkmcounts:124' lcnum = length(comb); */
  if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
    u1 = 0;
  } else {
    vlen = comb->size[0];
    u1 = comb->size[1];
    if (vlen >= u1) {
      u1 = vlen;
    }
  }
  /* 'getgkmcounts:125' n = length(alpha); */
  /*  disp(['# of support vectors: ' num2str(n)]) */
  /* 'getgkmcounts:127' fprintf("# of support vectors: %d\n", int32(n)); */
  printf("# of support vectors: %d\n", alpha->size[0]);
  fflush(stdout);
  /* 'getgkmcounts:128' alphasum = sum(alpha(alpha>0)); */
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
  emxInit_real_T(&temp, 1);
  i = temp->size[0];
  temp->size[0] = r2->size[0];
  emxEnsureCapacity_real_T(temp, i);
  temp_data = temp->data;
  loop_ub = r2->size[0];
  for (i = 0; i < loop_ub; i++) {
    temp_data[i] = alpha_data[r3[i] - 1];
  }
  alphasum = blockedSummation(temp, r2->size[0]);
  /* 'getgkmcounts:129' gkmc = zeros(4^k*lcnum,1); */
  c_tmp = pow(4.0, k);
  loop_ub_tmp_tmp = c_tmp * (double)u1;
  loop_ub_tmp = (int)loop_ub_tmp_tmp;
  i = gkmc->size[0];
  gkmc->size[0] = (int)loop_ub_tmp_tmp;
  emxEnsureCapacity_real_T(gkmc, i);
  gkmc_data = gkmc->data;
  emxFree_int32_T(&r2);
  for (i = 0; i < loop_ub_tmp; i++) {
    gkmc_data[i] = 0.0;
  }
  emxInit_boolean_T(&b_x, 1);
  /* 'getgkmcounts:130' GCpos = 0; */
  *GCpos = 0.0;
  /* 'getgkmcounts:131' GCneg = 0; */
  *GCneg = 0.0;
  /* 'getgkmcounts:132' np = sum(alpha>0); */
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
  /* 'getgkmcounts:133' nn = sum(alpha<0); */
  /* 'getgkmcounts:134' len = zeros(n,1); */
  i = len->size[0];
  len->size[0] = alpha->size[0];
  emxEnsureCapacity_uint32_T(len, i);
  len_data = len->data;
  loop_ub = alpha->size[0];
  for (i = 0; i < loop_ub; i++) {
    len_data[i] = 0U;
  }
  /* 'getgkmcounts:135' mat = zeros(4); */
  /* 'getgkmcounts:136' mat2 = zeros(4); */
  memset(&mat[0], 0, 16U * sizeof(double));
  memset(&mat2[0], 0, 16U * sizeof(double));
  /* 'getgkmcounts:137' pow = 4.^(0:l-1).'; */
  emxInit_real_T(&c_y, 2);
  M_data = c_y->data;
  if (l - 1.0 < 0.0) {
    c_y->size[0] = 1;
    c_y->size[1] = 0;
  } else {
    i = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    loop_ub = (int)floor(l - 1.0);
    c_y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(c_y, i);
    M_data = c_y->data;
    for (i = 0; i <= loop_ub; i++) {
      M_data[i] = i;
    }
  }
  emxInit_real_T(&r4, 2);
  i = r4->size[0] * r4->size[1];
  r4->size[0] = 1;
  r4->size[1] = c_y->size[1];
  emxEnsureCapacity_real_T(r4, i);
  r5 = r4->data;
  loop_ub = c_y->size[1];
  for (i = 0; i < loop_ub; i++) {
    curr_pos = M_data[i];
    r5[i] = pow(4.0, curr_pos);
  }
  emxInit_real_T(&M2, 2);
  /* 'getgkmcounts:138' n4 = 4^(l-1); */
  idx = pow(4.0, l - 1.0);
  /* 'getgkmcounts:139' slen = numel(sequences); */
  /* 'getgkmcounts:140' per = 10; */
  per = 10.0;
  /* 'getgkmcounts:141' a = 2; */
  a = 2.0;
  /* 'getgkmcounts:142' s = ''; */
  validatedHoleFilling[0].f1->size[0] = 1;
  validatedHoleFilling[0].f1->size[1] = 0;
  /* 'getgkmcounts:143' l2 = l-1; */
  /* 'getgkmcounts:144' lcnum2 = lcnum-rcnum; */
  lcnum2 = (double)u1 - rcnum;
  /* 'getgkmcounts:145' M2 = 0; */
  i = M2->size[0] * M2->size[1];
  M2->size[0] = 1;
  M2->size[1] = 1;
  emxEnsureCapacity_real_T(M2, i);
  M2_data = M2->data;
  M2_data[0] = 0.0;
  /* 'getgkmcounts:146' for i = 1:n */
  i = alpha->size[0];
  if (0 <= alpha->size[0] - 1) {
    if (1.0 > l) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = (int)l;
    }
    i1 = x->size[1];
    c_loop_ub = x->size[1];
  }
  emxInit_real_T(&ss, 2);
  emxInitMatrix_cell_wrap_0(C);
  emxInit_real_T(&M, 2);
  emxInit_int32_T(&r6, 2);
  emxInit_boolean_T(&c_x, 2);
  emxInit_real_T(&b_Ker, 2);
  emxInit_real_T(&r7, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'getgkmcounts:147' if mod(i, floor(n/10))==0 */
    if (b_mod((double)b_i + 1.0, floor((double)alpha->size[0] / 10.0)) == 0.0) {
      /* 'getgkmcounts:148' fprintf('%d...', int32(per)); */
      printf("%d...", (int)per);
      fflush(stdout);
      /* 'getgkmcounts:149' per = per+10; */
      per += 10.0;
    }
    /* 'getgkmcounts:151' if ver == 0 */
    if (ver == 0) {
      /* 'getgkmcounts:152' cur_seq = char(string(sequences{a})); */
      i2 = validatedHoleFilling[1].f1->size[0] *
           validatedHoleFilling[1].f1->size[1];
      validatedHoleFilling[1].f1->size[0] = 1;
      validatedHoleFilling[1].f1->size[1] =
          sequences_data[(int)a - 1].f1[0].f1->size[1];
      emxEnsureCapacity_char_T(validatedHoleFilling[1].f1, i2);
      loop_ub = sequences_data[(int)a - 1].f1[0].f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        validatedHoleFilling[1].f1->data[i2] =
            sequences_data[(int)a - 1].f1[0].f1->data[i2];
      }
      /* 'getgkmcounts:153' while ~strcmp('>', cur_seq(1)) */
      exitg2 = false;
      while ((!exitg2) && ('>' != validatedHoleFilling[1].f1->data[0])) {
        /* 'getgkmcounts:154' s = [s cur_seq]; */
        i2 = validatedHoleFilling[0].f1->size[1];
        b_k = validatedHoleFilling[0].f1->size[0] *
              validatedHoleFilling[0].f1->size[1];
        validatedHoleFilling[0].f1->size[0] = 1;
        validatedHoleFilling[0].f1->size[1] +=
            validatedHoleFilling[1].f1->size[1];
        emxEnsureCapacity_char_T(validatedHoleFilling[0].f1, b_k);
        loop_ub = validatedHoleFilling[1].f1->size[1];
        for (b_k = 0; b_k < loop_ub; b_k++) {
          validatedHoleFilling[0].f1->data[b_k + i2] =
              validatedHoleFilling[1].f1->data[b_k];
        }
        /* 'getgkmcounts:155' a = a+1; */
        a++;
        /* 'getgkmcounts:156' if a > slen */
        if (a > sequences->size[0]) {
          exitg2 = true;
        } else {
          /* 'getgkmcounts:159' cur_seq = char(string(sequences{a})); */
          i2 = validatedHoleFilling[1].f1->size[0] *
               validatedHoleFilling[1].f1->size[1];
          validatedHoleFilling[1].f1->size[0] = 1;
          b_k = (int)a - 1;
          validatedHoleFilling[1].f1->size[1] =
              sequences_data[b_k].f1[0].f1->size[1];
          emxEnsureCapacity_char_T(validatedHoleFilling[1].f1, i2);
          loop_ub = sequences_data[b_k].f1[0].f1->size[1];
          for (i2 = 0; i2 < loop_ub; i2++) {
            validatedHoleFilling[1].f1->data[i2] =
                sequences_data[b_k].f1[0].f1->data[i2];
          }
        }
      }
    } else {
      /* 'getgkmcounts:161' else */
      /* 'getgkmcounts:162' if ~isempty(sequences{i}) */
      /* 'getgkmcounts:163' s = char(string(sequences{i})); */
      i2 = validatedHoleFilling[0].f1->size[0] *
           validatedHoleFilling[0].f1->size[1];
      validatedHoleFilling[0].f1->size[0] = 1;
      validatedHoleFilling[0].f1->size[1] =
          sequences_data[b_i].f1[0].f1->size[1];
      emxEnsureCapacity_char_T(validatedHoleFilling[0].f1, i2);
      loop_ub = sequences_data[b_i].f1[0].f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        validatedHoleFilling[0].f1->data[i2] =
            sequences_data[b_i].f1[0].f1->data[i2];
      }
    }
    /* 'getgkmcounts:167' ss = letterconvert(s); */
    letterconvert(validatedHoleFilling[0].f1, ss);
    ss_data = ss->data;
    /* 'getgkmcounts:168' if RC */
    if (RC) {
      /* 'getgkmcounts:169' temp = zeros(4^k*(lcnum*2), 1); */
      loop_ub = (int)(c_tmp * ((double)u1 * 2.0));
      i2 = temp->size[0];
      temp->size[0] = loop_ub;
      emxEnsureCapacity_real_T(temp, i2);
      temp_data = temp->data;
      for (i2 = 0; i2 < loop_ub; i2++) {
        temp_data[i2] = 0.0;
      }
    } else {
      /* 'getgkmcounts:170' else */
      /* 'getgkmcounts:171' temp = zeros(4^k*(lcnum), 1); */
      i2 = temp->size[0];
      temp->size[0] = (int)loop_ub_tmp_tmp;
      emxEnsureCapacity_real_T(temp, i2);
      temp_data = temp->data;
      for (i2 = 0; i2 < loop_ub_tmp; i2++) {
        temp_data[i2] = 0.0;
      }
    }
    /* 'getgkmcounts:173' vec = ss(1:l)*pow; */
    i2 = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(c_y, i2);
    M_data = c_y->data;
    for (i2 = 0; i2 < b_loop_ub; i2++) {
      M_data[i2] = ss_data[i2];
    }
    if (b_loop_ub < 1) {
      curr_pos = 0.0;
    } else {
      curr_pos = cblas_ddot((blasint)b_loop_ub, &M_data[0], (blasint)1, &r5[0],
                            (blasint)1);
    }
    /* 'getgkmcounts:174' en = x(vec+1, :); */
    /* 'getgkmcounts:175' temp(en) = temp(en) + 1; */
    i2 = r6->size[0] * r6->size[1];
    r6->size[0] = 1;
    r6->size[1] = i1;
    emxEnsureCapacity_int32_T(r6, i2);
    r3 = r6->data;
    for (i2 = 0; i2 < c_loop_ub; i2++) {
      r3[i2] = (int)x_data[((int)(curr_pos + 1.0) + x->size[0] * i2) - 1];
    }
    i2 = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = r6->size[1];
    emxEnsureCapacity_real_T(c_y, i2);
    M_data = c_y->data;
    loop_ub = r6->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      M_data[i2] = 1.0;
    }
    loop_ub = c_y->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      temp_data[r3[i2] - 1] = 1.0;
    }
    /* 'getgkmcounts:176' for j = 2:length(ss)-l2 */
    i2 = (int)(((double)ss->size[1] - (l - 1.0)) + -1.0);
    for (j = 0; j < i2; j++) {
      /* 'getgkmcounts:177' vec = (vec-ss(j-1))/4+n4*ss(j+l2); */
      curr_pos = (curr_pos - ss_data[j]) / 4.0 +
                 idx * ss_data[(int)(((double)j + 2.0) + (l - 1.0)) - 1];
      /* 'getgkmcounts:178' en = x(vec+1, :); */
      /* 'getgkmcounts:179' temp(en) = temp(en) + 1; */
      loop_ub = x->size[1];
      b_k = r6->size[0] * r6->size[1];
      r6->size[0] = 1;
      r6->size[1] = x->size[1];
      emxEnsureCapacity_int32_T(r6, b_k);
      r3 = r6->data;
      for (b_k = 0; b_k < loop_ub; b_k++) {
        r3[b_k] = (int)x_data[((int)(curr_pos + 1.0) + x->size[0] * b_k) - 1];
      }
      b_k = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = r6->size[1];
      emxEnsureCapacity_real_T(c_y, b_k);
      M_data = c_y->data;
      loop_ub = r6->size[1];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        M_data[b_k] =
            temp_data[(int)x_data[((int)(curr_pos + 1.0) + x->size[0] * b_k) -
                                  1] -
                      1] +
            1.0;
      }
      loop_ub = c_y->size[1];
      for (b_k = 0; b_k < loop_ub; b_k++) {
        temp_data[r3[b_k] - 1] = M_data[b_k];
      }
    }
    /* 'getgkmcounts:181' if RC */
    if (RC) {
      /* 'getgkmcounts:182' temp =
       * temp(1:4^k*lcnum)+temp(4^k*lcnum+1:4^k*2*lcnum); */
      if (1.0 > loop_ub_tmp_tmp) {
        loop_ub = 0;
      } else {
        loop_ub = (int)loop_ub_tmp_tmp;
      }
      d = pow(4.0, k) * (double)u1 + 1.0;
      curr_pos = c_tmp * 2.0 * (double)u1;
      if (d > curr_pos) {
        i2 = 0;
        b_k = 0;
      } else {
        i2 = (int)d - 1;
        b_k = (int)curr_pos;
      }
      if (loop_ub == b_k - i2) {
        for (b_k = 0; b_k < loop_ub; b_k++) {
          temp_data[b_k] += temp_data[i2 + b_k];
        }
        i2 = temp->size[0];
        temp->size[0] = loop_ub;
        emxEnsureCapacity_real_T(temp, i2);
        temp_data = temp->data;
      } else {
        i_binary_expand_op(temp, loop_ub - 1, i2, b_k - 1);
        temp_data = temp->data;
      }
      /*          temp = sum(reshape(temp, [], 2), 2); */
    }
    /* 'getgkmcounts:185' if rcnum > 0 && RC */
    if ((rcnum > 0.0) && RC) {
      /* 'getgkmcounts:186'
       * temp(4^k*lcnum2+1:4^k*lcnum)=temp(4^k*lcnum2+1:4^k*lcnum)/sqrt(2); */
      d = c_tmp * lcnum2 + 1.0;
      curr_pos = pow(4.0, k) * (double)u1;
      if (d > curr_pos) {
        i2 = 1;
        b_k = -1;
        vlen = 0;
      } else {
        i2 = (int)d;
        b_k = (int)d - 2;
        vlen = (int)curr_pos;
      }
      loop_ub = (vlen - b_k) - 1;
      vlen = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = loop_ub;
      emxEnsureCapacity_real_T(c_y, vlen);
      M_data = c_y->data;
      for (vlen = 0; vlen < loop_ub; vlen++) {
        M_data[vlen] = temp_data[(i2 + vlen) - 1] / 1.4142135623730951;
      }
      loop_ub = c_y->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        temp_data[(b_k + i2) + 1] = M_data[i2];
      }
    }
    /* 'getgkmcounts:188' if LK == 1 */
    if (LK == 1) {
      /* 'getgkmcounts:189' gkmc = gkmc + alpha(i)/sqrt(temp'*temp)*temp; */
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
          gkmc_data[i2] += curr_pos * temp_data[i2];
        }
      } else {
        g_binary_expand_op(gkmc, curr_pos, temp);
        gkmc_data = gkmc->data;
      }
      /* 'getgkmcounts:190' len(i) = length(ss); */
      len_data[b_i] = (unsigned int)ss->size[1];
    } else {
      /* 'getgkmcounts:191' else */
      /* 'getgkmcounts:192' C = cell(4,2); */
      /* 'getgkmcounts:193' len(i) = length(ss); */
      len_data[b_i] = (unsigned int)ss->size[1];
      /* 'getgkmcounts:194' for ii = 1:4 */
      /* 'getgkmcounts:195' for j = 1:2 */
      /* 'getgkmcounts:196' C{ii,j} = zeros(len(i)-ll+1, ll); */
      vlen = (int)(((double)(unsigned int)ss->size[1] - ll) + 1.0);
      i2 = C[0].f1->size[0] * C[0].f1->size[1];
      C[0].f1->size[0] = vlen;
      C[0].f1->size[1] = (int)ll;
      emxEnsureCapacity_real_T(C[0].f1, i2);
      loop_ub = (int)(((double)(unsigned int)ss->size[1] - ll) + 1.0) * (int)ll;
      for (i2 = 0; i2 < loop_ub; i2++) {
        C[0].f1->data[i2] = 0.0;
      }
      /* 'getgkmcounts:196' C{ii,j} = zeros(len(i)-ll+1, ll); */
      curr_pos = ((double)len_data[b_i] - ll) + 1.0;
      i2 = C[4].f1->size[0] * C[4].f1->size[1];
      C[4].f1->size[0] = (int)curr_pos;
      C[4].f1->size[1] = (int)ll;
      emxEnsureCapacity_real_T(C[4].f1, i2);
      vlen = (int)curr_pos * (int)ll;
      for (i2 = 0; i2 < vlen; i2++) {
        C[4].f1->data[i2] = 0.0;
      }
      /* 'getgkmcounts:195' for j = 1:2 */
      /* 'getgkmcounts:196' C{ii,j} = zeros(len(i)-ll+1, ll); */
      i2 = C[1].f1->size[0] * C[1].f1->size[1];
      C[1].f1->size[0] = (int)curr_pos;
      C[1].f1->size[1] = (int)ll;
      emxEnsureCapacity_real_T(C[1].f1, i2);
      for (i2 = 0; i2 < vlen; i2++) {
        C[1].f1->data[i2] = 0.0;
      }
      /* 'getgkmcounts:196' C{ii,j} = zeros(len(i)-ll+1, ll); */
      i2 = C[5].f1->size[0] * C[5].f1->size[1];
      C[5].f1->size[0] = (int)curr_pos;
      C[5].f1->size[1] = (int)ll;
      emxEnsureCapacity_real_T(C[5].f1, i2);
      for (i2 = 0; i2 < vlen; i2++) {
        C[5].f1->data[i2] = 0.0;
      }
      /* 'getgkmcounts:195' for j = 1:2 */
      /* 'getgkmcounts:196' C{ii,j} = zeros(len(i)-ll+1, ll); */
      i2 = C[2].f1->size[0] * C[2].f1->size[1];
      C[2].f1->size[0] = (int)curr_pos;
      C[2].f1->size[1] = (int)ll;
      emxEnsureCapacity_real_T(C[2].f1, i2);
      for (i2 = 0; i2 < vlen; i2++) {
        C[2].f1->data[i2] = 0.0;
      }
      /* 'getgkmcounts:196' C{ii,j} = zeros(len(i)-ll+1, ll); */
      i2 = C[6].f1->size[0] * C[6].f1->size[1];
      C[6].f1->size[0] = (int)curr_pos;
      C[6].f1->size[1] = (int)ll;
      emxEnsureCapacity_real_T(C[6].f1, i2);
      for (i2 = 0; i2 < vlen; i2++) {
        C[6].f1->data[i2] = 0.0;
      }
      /* 'getgkmcounts:195' for j = 1:2 */
      /* 'getgkmcounts:196' C{ii,j} = zeros(len(i)-ll+1, ll); */
      i2 = C[3].f1->size[0] * C[3].f1->size[1];
      C[3].f1->size[0] = (int)curr_pos;
      C[3].f1->size[1] = (int)ll;
      emxEnsureCapacity_real_T(C[3].f1, i2);
      for (i2 = 0; i2 < vlen; i2++) {
        C[3].f1->data[i2] = 0.0;
      }
      /* 'getgkmcounts:196' C{ii,j} = zeros(len(i)-ll+1, ll); */
      i2 = C[7].f1->size[0] * C[7].f1->size[1];
      C[7].f1->size[0] = (int)curr_pos;
      C[7].f1->size[1] = (int)ll;
      emxEnsureCapacity_real_T(C[7].f1, i2);
      for (i2 = 0; i2 < vlen; i2++) {
        C[7].f1->data[i2] = 0.0;
      }
      /* 'getgkmcounts:199' for ii = 1:len(i) */
      i2 = (int)len_data[b_i];
      for (ii = 0; ii < i2; ii++) {
        /* 'getgkmcounts:200' for j = 0:ll-1 */
        b_k = (int)((ll - 1.0) + 1.0);
        for (j = 0; j < b_k; j++) {
          /* 'getgkmcounts:201' if ii - j <= len(i)-ll+1 && ii - j > 0 */
          vlen = ii - j;
          if ((vlen + 1 <= curr_pos) && (vlen + 1 > 0)) {
            /* 'getgkmcounts:202' C{ss(ii)+1,1}(ii-j,j+1) = 1; */
            C[(int)(ss_data[ii] + 1.0) - 1]
                .f1
                ->data[vlen + C[(int)(ss_data[ii] + 1.0) - 1].f1->size[0] * j] =
                1.0;
          }
        }
      }
      /* 'getgkmcounts:206' if RC */
      if (RC) {
        /* 'getgkmcounts:207' for ii = 0:3 */
        /* 'getgkmcounts:208' C{ii+1,2} = rot90(C{4-ii,1}, 2); */
        rot90(C[3].f1, C[4].f1);
        /* 'getgkmcounts:208' C{ii+1,2} = rot90(C{4-ii,1}, 2); */
        rot90(C[2].f1, C[5].f1);
        /* 'getgkmcounts:208' C{ii+1,2} = rot90(C{4-ii,1}, 2); */
        rot90(C[1].f1, C[6].f1);
        /* 'getgkmcounts:208' C{ii+1,2} = rot90(C{4-ii,1}, 2); */
        rot90(C[0].f1, C[7].f1);
        /* 'getgkmcounts:210' M2 = zeros(len(i)-ll+1); */
        i2 = M2->size[0] * M2->size[1];
        M2->size[0] = (int)curr_pos;
        M2->size[1] = (int)curr_pos;
        emxEnsureCapacity_real_T(M2, i2);
        M2_data = M2->data;
        loop_ub = (int)curr_pos * (int)curr_pos;
        for (i2 = 0; i2 < loop_ub; i2++) {
          M2_data[i2] = 0.0;
        }
      }
      /* 'getgkmcounts:212' M = zeros(len(i)-ll+1); */
      i2 = M->size[0] * M->size[1];
      M->size[0] = (int)curr_pos;
      M->size[1] = (int)curr_pos;
      emxEnsureCapacity_real_T(M, i2);
      M_data = M->data;
      loop_ub = (int)curr_pos * (int)curr_pos;
      for (i2 = 0; i2 < loop_ub; i2++) {
        M_data[i2] = 0.0;
      }
      /* 'getgkmcounts:213' for ii = 1:4 */
      for (ii = 0; ii < 4; ii++) {
        /* 'getgkmcounts:214' M = M + C{ii,1}*C{ii,1}.'; */
        mtimes(C[ii].f1, C[ii].f1, b_Ker);
        b_Ker_data = b_Ker->data;
        if ((M->size[0] == b_Ker->size[0]) && (M->size[1] == b_Ker->size[1])) {
          loop_ub = M->size[0] * M->size[1];
          for (i2 = 0; i2 < loop_ub; i2++) {
            M_data[i2] += b_Ker_data[i2];
          }
        } else {
          plus(M, b_Ker);
          M_data = M->data;
        }
        /* 'getgkmcounts:215' if RC */
        if (RC) {
          /* 'getgkmcounts:216' M2 = M2 + C{ii,1}*C{ii,2}.'; */
          mtimes(C[ii].f1, C[ii + 4].f1, b_Ker);
          b_Ker_data = b_Ker->data;
          if ((M2->size[0] == b_Ker->size[0]) &&
              (M2->size[1] == b_Ker->size[1])) {
            loop_ub = M2->size[0] * M2->size[1];
            for (i2 = 0; i2 < loop_ub; i2++) {
              M2_data[i2] += b_Ker_data[i2];
            }
          } else {
            plus(M2, b_Ker);
            M2_data = M2->data;
          }
        }
      }
      /* 'getgkmcounts:219' if RC */
      if (RC) {
        /* 'getgkmcounts:220' norm = sqrt(sum(sum(Ker(1 + M)+Ker(1 + M2)))); */
        i2 = b_Ker->size[0] * b_Ker->size[1];
        b_Ker->size[0] = M->size[0];
        b_Ker->size[1] = M->size[1];
        emxEnsureCapacity_real_T(b_Ker, i2);
        b_Ker_data = b_Ker->data;
        loop_ub = M->size[0] * M->size[1];
        for (i2 = 0; i2 < loop_ub; i2++) {
          b_Ker_data[i2] = Ker_data[(int)(M_data[i2] + 1.0) - 1];
        }
        i2 = M->size[0] * M->size[1];
        M->size[0] = M2->size[0];
        M->size[1] = M2->size[1];
        emxEnsureCapacity_real_T(M, i2);
        M_data = M->data;
        loop_ub = M2->size[0] * M2->size[1];
        for (i2 = 0; i2 < loop_ub; i2++) {
          M_data[i2] = Ker_data[(int)(M2_data[i2] + 1.0) - 1];
        }
        if ((b_Ker->size[0] == M->size[0]) && (b_Ker->size[1] == M->size[1])) {
          i2 = r7->size[0] * r7->size[1];
          r7->size[0] = b_Ker->size[0];
          r7->size[1] = b_Ker->size[1];
          emxEnsureCapacity_real_T(r7, i2);
          r8 = r7->data;
          loop_ub = b_Ker->size[0] * b_Ker->size[1];
          for (i2 = 0; i2 < loop_ub; i2++) {
            r8[i2] = b_Ker_data[i2] + M_data[i2];
          }
          c_sum(r7, c_y);
        } else {
          h_binary_expand_op(c_y, b_Ker, M);
        }
        curr_pos = sqrt(sum(c_y));
      } else {
        /* 'getgkmcounts:221' else */
        /* 'getgkmcounts:222' norm = sqrt(sum(sum(Ker(1 + M)))); */
        i2 = b_Ker->size[0] * b_Ker->size[1];
        b_Ker->size[0] = M->size[0];
        b_Ker->size[1] = M->size[1];
        emxEnsureCapacity_real_T(b_Ker, i2);
        b_Ker_data = b_Ker->data;
        loop_ub = M->size[0] * M->size[1];
        for (i2 = 0; i2 < loop_ub; i2++) {
          b_Ker_data[i2] = Ker_data[(int)(M_data[i2] + 1.0) - 1];
        }
        c_sum(b_Ker, c_y);
        curr_pos = sqrt(sum(c_y));
      }
      /* 'getgkmcounts:224' gkmc = gkmc + alpha(i)/norm*temp; */
      curr_pos = alpha_data[b_i] / curr_pos;
      loop_ub = gkmc->size[0];
      if (gkmc->size[0] == temp->size[0]) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          gkmc_data[i2] += curr_pos * temp_data[i2];
        }
      } else {
        g_binary_expand_op(gkmc, curr_pos, temp);
        gkmc_data = gkmc->data;
      }
    }
    /* 'getgkmcounts:226' if alpha(i) > 0 */
    if (alpha_data[b_i] > 0.0) {
      /* 'getgkmcounts:227' GCpos = GCpos + alpha(i)*(sum(ss==1)+sum(ss==2)); */
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
        ii = 0;
      } else {
        ii = b_x_data[0];
        for (b_k = 2; b_k <= vlen; b_k++) {
          ii += b_x_data[b_k - 1];
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
      *GCpos += alpha_data[b_i] * ((double)ii + (double)loop_ub);
      /* 'getgkmcounts:228' ss = ss+1; */
      i2 = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity_real_T(ss, i2);
      ss_data = ss->data;
      loop_ub = ss->size[1] - 1;
      for (i2 = 0; i2 <= loop_ub; i2++) {
        ss_data[i2]++;
      }
      /* 'getgkmcounts:229' for j = 1:length(ss)-1 */
      i2 = ss->size[1];
      for (j = 0; j <= i2 - 2; j++) {
        /* 'getgkmcounts:230' mat2(ss(j),ss(j+1)) = mat2(ss(j),ss(j+1)) +
         * abs(alpha(i)); */
        vlen = ((int)ss_data[j] + (((int)ss_data[j + 1] - 1) << 2)) - 1;
        mat2[vlen] += fabs(alpha_data[b_i]);
      }
    } else {
      /* 'getgkmcounts:232' else */
      /* 'getgkmcounts:233' GCneg = GCneg + alpha(i)*(sum(ss==1)+sum(ss==2)); */
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
        ii = 0;
      } else {
        ii = b_x_data[0];
        for (b_k = 2; b_k <= vlen; b_k++) {
          ii += b_x_data[b_k - 1];
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
      *GCneg += alpha_data[b_i] * ((double)ii + (double)loop_ub);
      /* 'getgkmcounts:234' ss = ss+1; */
      i2 = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity_real_T(ss, i2);
      ss_data = ss->data;
      loop_ub = ss->size[1] - 1;
      for (i2 = 0; i2 <= loop_ub; i2++) {
        ss_data[i2]++;
      }
      /* 'getgkmcounts:235' for j = 1:length(ss)-1 */
      i2 = ss->size[1];
      for (j = 0; j <= i2 - 2; j++) {
        /* 'getgkmcounts:236' mat(ss(j),ss(j+1)) = mat(ss(j),ss(j+1)) +
         * abs(alpha(i)); */
        vlen = ((int)ss_data[j] + (((int)ss_data[j + 1] - 1) << 2)) - 1;
        mat[vlen] += fabs(alpha_data[b_i]);
      }
    }
    /* 'getgkmcounts:239' s = ''; */
    validatedHoleFilling[0].f1->size[0] = 1;
    validatedHoleFilling[0].f1->size[1] = 0;
    /* 'getgkmcounts:240' a = a+1; */
    a++;
  }
  emxFree_real_T(&r7);
  emxFree_real_T(&b_Ker);
  emxFree_boolean_T(&c_x);
  emxFreeMatrix_cell_wrap_81(validatedHoleFilling);
  emxFree_real_T(&c_y);
  emxFree_int32_T(&r6);
  emxFree_real_T(&r4);
  emxFree_cell_wrap_9(&sequences);
  emxFree_real_T(&M);
  emxFreeMatrix_cell_wrap_0(C);
  emxFree_real_T(&ss);
  emxFree_real_T(&M2);
  emxFree_real_T(&Ker);
  emxFree_real_T(&x);
  emxFree_real_T(&alpha);
  /* 'getgkmcounts:242' fprintf('\n') */
  printf("\n");
  fflush(stdout);
  /* 'getgkmcounts:243' GCpos = GCpos/alphasum/mean(len(1:np)); */
  if (1 > b_y) {
    loop_ub = 0;
  } else {
    loop_ub = nz;
  }
  i = temp->size[0];
  temp->size[0] = loop_ub;
  emxEnsureCapacity_real_T(temp, i);
  temp_data = temp->data;
  for (i = 0; i < loop_ub; i++) {
    temp_data[i] = len_data[i];
  }
  *GCpos =
      *GCpos / alphasum / (blockedSummation(temp, loop_ub) / (double)loop_ub);
  /* 'getgkmcounts:244' GCneg = -1*GCneg/alphasum/mean(len(np+1:end)); */
  if ((double)nz + 1.0 > len->size[0]) {
    i = -1;
    i1 = -1;
  } else {
    i = b_y - 1;
    i1 = len->size[0] - 1;
  }
  loop_ub = i1 - i;
  i1 = temp->size[0];
  temp->size[0] = loop_ub;
  emxEnsureCapacity_real_T(temp, i1);
  temp_data = temp->data;
  for (i1 = 0; i1 < loop_ub; i1++) {
    temp_data[i1] = len_data[(i + i1) + 1];
  }
  emxFree_uint32_T(&len);
  *GCneg =
      -*GCneg / alphasum / (blockedSummation(temp, loop_ub) / (double)loop_ub);
  /* 'getgkmcounts:245' if RC */
  emxFree_real_T(&temp);
  if (RC) {
    /* 'getgkmcounts:246' mat = (mat+rot90(rot90(mat,3)'))/2; */
    b_rot90(mat, dv);
    for (i = 0; i < 4; i++) {
      i1 = i << 2;
      dv2[i1] = dv[i];
      dv2[i1 + 1] = dv[i + 4];
      dv2[i1 + 2] = dv[i + 8];
      dv2[i1 + 3] = dv[i + 12];
    }
    c_rot90(dv2, dv);
    for (i = 0; i < 16; i++) {
      mat[i] = (mat[i] + dv[i]) / 2.0;
    }
    /* 'getgkmcounts:247' mat2 = (mat2+rot90(rot90(mat2,3)'))/2; */
    b_rot90(mat2, dv);
    for (i = 0; i < 4; i++) {
      i1 = i << 2;
      dv2[i1] = dv[i];
      dv2[i1 + 1] = dv[i + 4];
      dv2[i1 + 2] = dv[i + 8];
      dv2[i1 + 3] = dv[i + 12];
    }
    c_rot90(dv2, dv);
    for (i = 0; i < 16; i++) {
      mat2[i] = (mat2[i] + dv[i]) / 2.0;
    }
  }
  /* 'getgkmcounts:249' mat = mat./repmat(sum(mat,2),1,4); */
  d_sum(mat, dv1);
  repmat(dv1, dv);
  for (i = 0; i < 16; i++) {
    mat[i] /= dv[i];
  }
  /* 'getgkmcounts:250' mat2 = mat2./repmat(sum(mat2,2),1,4); */
  d_sum(mat2, dv1);
  repmat(dv1, dv);
  for (i = 0; i < 16; i++) {
    mat2[i] /= dv[i];
  }
}

/* End of code generation (getgkmcounts.c) */
