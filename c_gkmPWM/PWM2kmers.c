/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PWM2kmers.c
 *
 * Code generation for function 'PWM2kmers'
 *
 */

/* Include files */
#include "PWM2kmers.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include "repmat.h"
#include "rot90.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void q_binary_expand_op(emxArray_cell_wrap_14 *KC,
                               const emxArray_real_T *ind,
                               const emxArray_real_T *x, int i2, int j,
                               const emxArray_real_T *r1,
                               const emxArray_cell_wrap_14 *ktree, double k,
                               const emxArray_cell_wrap_14 *ktree2);

static void r_binary_expand_op(emxArray_cell_wrap_14 *KC,
                               const emxArray_real_T *ind, int i2,
                               const emxArray_cell_wrap_14 *ktree, double k,
                               const emxArray_cell_wrap_14 *ktree2);

static void s_binary_expand_op(emxArray_real_T *kweig, double n);

/* Function Definitions */
static void q_binary_expand_op(emxArray_cell_wrap_14 *KC,
                               const emxArray_real_T *ind,
                               const emxArray_real_T *x, int i2, int j,
                               const emxArray_real_T *r1,
                               const emxArray_cell_wrap_14 *ktree, double k,
                               const emxArray_cell_wrap_14 *ktree2)
{
  const cell_wrap_14 *ktree2_data;
  const cell_wrap_14 *ktree_data;
  cell_wrap_14 *KC_data;
  const double *ind_data;
  const double *r;
  const double *x_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  ktree2_data = ktree2->data;
  ktree_data = ktree->data;
  r = r1->data;
  x_data = x->data;
  ind_data = ind->data;
  KC_data = KC->data;
  i = KC_data[(int)ind_data[(int)x_data[i2 + x->size[0] * j] - 1] - 1]
          .f1->size[0];
  if (ktree2_data[(int)k - 1].f1->size[0] == 1) {
    if (ktree_data[(int)k - 1].f1->size[0] == 1) {
      KC_data[(int)ind_data[(int)x_data[i2 + x->size[0] * j] - 1] - 1]
          .f1->size[0] = r1->size[0];
    } else {
      KC_data[(int)ind_data[(int)x_data[i2 + x->size[0] * j] - 1] - 1]
          .f1->size[0] = ktree_data[(int)k - 1].f1->size[0];
    }
  } else {
    KC_data[(int)ind_data[(int)x_data[i2 + x->size[0] * j] - 1] - 1]
        .f1->size[0] = ktree2_data[(int)k - 1].f1->size[0];
  }
  emxEnsureCapacity_real_T(
      KC_data[(int)ind_data[(int)x_data[i2 + x->size[0] * j] - 1] - 1].f1, i);
  stride_0_0 = (r1->size[0] != 1);
  stride_1_0 = (ktree_data[(int)k - 1].f1->size[0] != 1);
  stride_2_0 = (ktree2_data[(int)k - 1].f1->size[0] != 1);
  if (ktree2_data[(int)k - 1].f1->size[0] == 1) {
    if (ktree_data[(int)k - 1].f1->size[0] == 1) {
      loop_ub = r1->size[0];
    } else {
      loop_ub = ktree_data[(int)k - 1].f1->size[0];
    }
  } else {
    loop_ub = ktree2_data[(int)k - 1].f1->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    KC_data[(int)ind_data[(int)x_data[i2 + x->size[0] * j] - 1] - 1]
        .f1->data[i] =
        (r[i * stride_0_0] + ktree_data[(int)k - 1].f1->data[i * stride_1_0]) +
        ktree2_data[(int)k - 1].f1->data[i * stride_2_0];
  }
}

static void r_binary_expand_op(emxArray_cell_wrap_14 *KC,
                               const emxArray_real_T *ind, int i2,
                               const emxArray_cell_wrap_14 *ktree, double k,
                               const emxArray_cell_wrap_14 *ktree2)
{
  const cell_wrap_14 *ktree2_data;
  const cell_wrap_14 *ktree_data;
  cell_wrap_14 *KC_data;
  emxArray_real_T *b_KC;
  const double *ind_data;
  double *b_KC_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  ktree2_data = ktree2->data;
  ktree_data = ktree->data;
  ind_data = ind->data;
  KC_data = KC->data;
  emxInit_real_T(&b_KC, 1);
  i = b_KC->size[0];
  if (ktree2_data[(int)k - 1].f1->size[0] == 1) {
    if (ktree_data[(int)k - 1].f1->size[0] == 1) {
      b_KC->size[0] = KC_data[(int)ind_data[i2] - 1].f1->size[0];
    } else {
      b_KC->size[0] = ktree_data[(int)k - 1].f1->size[0];
    }
  } else {
    b_KC->size[0] = ktree2_data[(int)k - 1].f1->size[0];
  }
  emxEnsureCapacity_real_T(b_KC, i);
  b_KC_data = b_KC->data;
  stride_0_0 = (KC_data[(int)ind_data[i2] - 1].f1->size[0] != 1);
  stride_1_0 = (ktree_data[(int)k - 1].f1->size[0] != 1);
  stride_2_0 = (ktree2_data[(int)k - 1].f1->size[0] != 1);
  if (ktree2_data[(int)k - 1].f1->size[0] == 1) {
    if (ktree_data[(int)k - 1].f1->size[0] == 1) {
      loop_ub = KC_data[(int)ind_data[i2] - 1].f1->size[0];
    } else {
      loop_ub = ktree_data[(int)k - 1].f1->size[0];
    }
  } else {
    loop_ub = ktree2_data[(int)k - 1].f1->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_KC_data[i] = (KC_data[(int)ind_data[i2] - 1].f1->data[i * stride_0_0] +
                    ktree_data[(int)k - 1].f1->data[i * stride_1_0]) +
                   ktree2_data[(int)k - 1].f1->data[i * stride_2_0];
  }
  i = KC_data[(int)ind_data[i2] - 1].f1->size[0];
  KC_data[(int)ind_data[i2] - 1].f1->size[0] = b_KC->size[0];
  emxEnsureCapacity_real_T(KC_data[(int)ind_data[i2] - 1].f1, i);
  loop_ub = b_KC->size[0];
  for (i = 0; i < loop_ub; i++) {
    KC_data[(int)ind_data[i2] - 1].f1->data[i] = b_KC_data[i];
  }
  emxFree_real_T(&b_KC);
}

static void s_binary_expand_op(emxArray_real_T *kweig, double n)
{
  emxArray_real_T *b_kweig;
  double *b_kweig_data;
  double *kweig_data;
  int i;
  int loop_ub;
  int stride_0_0;
  kweig_data = kweig->data;
  emxInit_real_T(&b_kweig, 1);
  i = b_kweig->size[0];
  if ((int)n == 1) {
    b_kweig->size[0] = kweig->size[0];
  } else {
    b_kweig->size[0] = (int)n;
  }
  emxEnsureCapacity_real_T(b_kweig, i);
  b_kweig_data = b_kweig->data;
  stride_0_0 = (kweig->size[0] != 1);
  if ((int)n == 1) {
    loop_ub = kweig->size[0];
  } else {
    loop_ub = (int)n;
  }
  for (i = 0; i < loop_ub; i++) {
    b_kweig_data[i] = kweig_data[i * stride_0_0];
  }
  i = kweig->size[0];
  kweig->size[0] = b_kweig->size[0];
  emxEnsureCapacity_real_T(kweig, i);
  kweig_data = kweig->data;
  loop_ub = b_kweig->size[0];
  for (i = 0; i < loop_ub; i++) {
    kweig_data[i] = b_kweig_data[i];
  }
  emxFree_real_T(&b_kweig);
}

/*
 * function [kweig] = PWM2kmers(mat,negmat,c,s,ind,indloc,x,l,k,rcnum)
 */
void PWM2kmers(const emxArray_real_T *mat, const double negmat[16],
               const emxArray_real_T *c, const emxArray_real_T *s,
               const emxArray_real_T *ind, const emxArray_real_T *indloc,
               const emxArray_real_T *x, double l, double k, double rcnum,
               emxArray_real_T *kweig)
{
  cell_wrap_13 *p_data;
  cell_wrap_14 *KC_data;
  cell_wrap_14 *ktree2_data;
  cell_wrap_14 *ktree_data;
  emxArray_cell_wrap_13 *p;
  emxArray_cell_wrap_14 *KC;
  emxArray_cell_wrap_14 *ktree;
  emxArray_cell_wrap_14 *ktree2;
  emxArray_real_T *a;
  emxArray_real_T *b_ktree;
  emxArray_real_T *indvec;
  emxArray_real_T *mat2;
  emxArray_real_T *r;
  emxArray_real_T *repmatt;
  emxArray_real_T *sPWM;
  emxArray_real_T *sPWM2;
  emxArray_uint32_T *X;
  double b_p[16];
  double b_matt[4];
  const double *c_data;
  const double *ind_data;
  const double *indloc_data;
  const double *mat_data;
  const double *s_data;
  const double *x_data;
  double M;
  double b_i3;
  double c_tmp;
  double d;
  double d1;
  double matt;
  double n;
  double n_tmp;
  double *indvec_data;
  double *kweig_data;
  double *mat2_data;
  double *r1;
  double *repmatt_data;
  double *sPWM2_data;
  double *sPWM_data;
  int b_i;
  int b_i2;
  int i;
  int i1;
  int i2;
  int i3;
  int idx;
  int j;
  int loop_ub;
  int m;
  int ni;
  int rx;
  unsigned int *X_data;
  x_data = x->data;
  indloc_data = indloc->data;
  ind_data = ind->data;
  s_data = s->data;
  c_data = c->data;
  mat_data = mat->data;
  emxInit_cell_wrap_13(&p);
  /* makes gkm-pwm faster */
  /* 'PWM2kmers:3' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_13(p, i);
  p_data = p->data;
  /* 'PWM2kmers:4' p = coder.nullcopy(p); */
  /* 'PWM2kmers:5' p{1} = eye(4); */
  for (i = 0; i < 16; i++) {
    p_data[0].f1[i] = 0.0;
  }
  p_data[0].f1[0] = 1.0;
  p_data[0].f1[5] = 1.0;
  p_data[0].f1[10] = 1.0;
  p_data[0].f1[15] = 1.0;
  /* 'PWM2kmers:6' for i = 1:l-1 */
  i = (int)(l - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:7' p{i+1} = p{i}*negmat; */
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < 4; i2++) {
        ni = i2 << 2;
        b_p[i1 + ni] = ((p_data[b_i].f1[i1] * negmat[ni] +
                         p_data[b_i].f1[i1 + 4] * negmat[ni + 1]) +
                        p_data[b_i].f1[i1 + 8] * negmat[ni + 2]) +
                       p_data[b_i].f1[i1 + 12] * negmat[ni + 3];
      }
    }
    for (i1 = 0; i1 < 16; i1++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i1] = b_p[i1];
    }
  }
  /* 'PWM2kmers:9' m = length(x); */
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    m = 0;
  } else {
    idx = x->size[0];
    m = x->size[1];
    if (idx >= m) {
      m = idx;
    }
  }
  /* 'PWM2kmers:10' M = length(mat)-l; */
  idx = mat->size[0];
  if (idx < 4) {
    idx = 4;
  }
  if (mat->size[0] == 0) {
    idx = 0;
  }
  M = (double)idx - l;
  /* 'PWM2kmers:11' n = 4^k*length(c); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    ni = 0;
  } else {
    idx = c->size[0];
    ni = c->size[1];
    if (idx >= ni) {
      ni = idx;
    }
  }
  emxInit_real_T(&mat2, 2);
  n_tmp = pow(4.0, k);
  n = n_tmp * (double)ni;
  /* number of possible k-mers */
  /* 'PWM2kmers:12' mat2 = rot90(mat,2); */
  c_rot90(mat, mat2);
  mat2_data = mat2->data;
  /* 'PWM2kmers:13' kweig = zeros(n,1); */
  idx = (int)n;
  i = kweig->size[0];
  kweig->size[0] = (int)n;
  emxEnsureCapacity_real_T(kweig, i);
  kweig_data = kweig->data;
  for (i = 0; i < idx; i++) {
    kweig_data[i] = 0.0;
  }
  emxInit_cell_wrap_14(&ktree);
  /* 'PWM2kmers:14' kweig2= zeros(n,1); */
  /* 'PWM2kmers:15' ktree = cell(k,1); */
  idx = (int)k;
  i = ktree->size[0];
  ktree->size[0] = (int)k;
  emxEnsureCapacity_cell_wrap_14(ktree, i);
  ktree_data = ktree->data;
  for (i = 0; i < idx; i++) {
    ktree_data[i].f1->size[0] = 0;
  }
  emxInit_cell_wrap_14(&ktree2);
  /* 'PWM2kmers:16' ktree = coder.nullcopy(ktree); */
  /* 'PWM2kmers:17' ktree2 = cell(k,1); */
  i = ktree2->size[0];
  ktree2->size[0] = (int)k;
  emxEnsureCapacity_cell_wrap_14(ktree2, i);
  ktree2_data = ktree2->data;
  for (i = 0; i < idx; i++) {
    ktree2_data[i].f1->size[0] = 0;
  }
  emxInit_uint32_T(&X);
  /* 'PWM2kmers:18' ktree2 = coder.nullcopy(ktree2); */
  /* 'PWM2kmers:19' [rx,cx] = size(x); */
  rx = x->size[0];
  /* 'PWM2kmers:20' X = cx*ones(length(mat)-l+1,1); */
  idx = mat->size[0];
  if (idx < 4) {
    idx = 4;
  }
  if (mat->size[0] == 0) {
    idx = 0;
  }
  loop_ub = (int)(((double)idx - l) + 1.0);
  i = X->size[0];
  X->size[0] = loop_ub;
  emxEnsureCapacity_uint32_T(X, i);
  X_data = X->data;
  for (i = 0; i < loop_ub; i++) {
    X_data[i] = 1U;
  }
  loop_ub = X->size[0];
  for (i = 0; i < loop_ub; i++) {
    X_data[i] *= x->size[1];
  }
  /* 'PWM2kmers:21' KC = cell(length(c),1); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    ni = 0;
  } else {
    idx = c->size[0];
    ni = c->size[1];
    if (idx >= ni) {
      ni = idx;
    }
  }
  emxInit_cell_wrap_14(&KC);
  i = KC->size[0];
  KC->size[0] = ni;
  emxEnsureCapacity_cell_wrap_14(KC, i);
  KC_data = KC->data;
  for (i = 0; i < ni; i++) {
    KC_data[i].f1->size[0] = 0;
  }
  /* 'PWM2kmers:22' KC = coder.nullcopy(KC); */
  /* 'PWM2kmers:23' for i = 1:length(c) */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    ni = 0;
  } else {
    idx = c->size[0];
    ni = c->size[1];
    if (idx >= ni) {
      ni = idx;
    }
  }
  for (b_i = 0; b_i < ni; b_i++) {
    /* 'PWM2kmers:24' KC{i} = zeros(4^k,1); */
    idx = (int)n_tmp;
    i = KC_data[b_i].f1->size[0];
    KC_data[b_i].f1->size[0] = (int)n_tmp;
    emxEnsureCapacity_real_T(KC_data[b_i].f1, i);
    for (i = 0; i < idx; i++) {
      KC_data[b_i].f1->data[i] = 0.0;
    }
  }
  /* 'PWM2kmers:26' for i = 1:cx */
  i = x->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:27' X(i) = i; */
    X_data[b_i] = (unsigned int)(b_i + 1);
  }
  /* 'PWM2kmers:29' for i = 2:k */
  i = (int)(k + -1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:30' ktree{i} = zeros(4^i,1); */
    c_tmp = pow(4.0, (double)b_i + 2.0);
    idx = (int)pow(4.0, (double)b_i + 2.0);
    i1 = ktree_data[b_i + 1].f1->size[0];
    ktree_data[b_i + 1].f1->size[0] = (int)c_tmp;
    emxEnsureCapacity_real_T(ktree_data[b_i + 1].f1, i1);
    for (i1 = 0; i1 < idx; i1++) {
      ktree_data[b_i + 1].f1->data[i1] = 0.0;
    }
    /* 'PWM2kmers:31' ktree2{i} = zeros(4^i,1); */
    i1 = ktree2_data[b_i + 1].f1->size[0];
    ktree2_data[b_i + 1].f1->size[0] = (int)c_tmp;
    emxEnsureCapacity_real_T(ktree2_data[b_i + 1].f1, i1);
    for (i1 = 0; i1 < idx; i1++) {
      ktree2_data[b_i + 1].f1->data[i1] = 0.0;
    }
  }
  /* 'PWM2kmers:33' a = 0; */
  /* 'PWM2kmers:34' for i = 0:M */
  i = (int)(M + 1.0);
  emxInit_real_T(&indvec, 2);
  emxInit_real_T(&sPWM, 2);
  emxInit_real_T(&sPWM2, 2);
  emxInit_real_T(&repmatt, 1);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&r, 1);
  emxInit_real_T(&b_ktree, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:35' if i == M-1 */
    if (b_i == M - 1.0) {
      /* 'PWM2kmers:36' m = length(c); */
      if ((c->size[0] == 0) || (c->size[1] == 0)) {
        m = 0;
      } else {
        idx = c->size[0];
        m = c->size[1];
        if (idx >= m) {
          m = idx;
        }
      }
    }
    /* 'PWM2kmers:38' for i2 = 1:m */
    for (b_i2 = 0; b_i2 < m; b_i2++) {
      /* 'PWM2kmers:39' if ~(i == M-1 && i2 > rx && i2 ~= m) */
      if ((b_i != M - 1.0) || (b_i2 + 1 <= rx) || (b_i2 + 1 == m)) {
        /* 'PWM2kmers:40' indvec = c(i2,:)+i; */
        loop_ub = c->size[1];
        i1 = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i1);
        indvec_data = indvec->data;
        for (i1 = 0; i1 < loop_ub; i1++) {
          indvec_data[i1] = c_data[b_i2 + c->size[0] * i1] + (double)b_i;
        }
        /* 'PWM2kmers:41' loc = indloc(indvec); */
        /* 'PWM2kmers:42' sPWM = mat(indvec,:).'; */
        i1 = repmatt->size[0];
        repmatt->size[0] = indvec->size[1];
        emxEnsureCapacity_real_T(repmatt, i1);
        repmatt_data = repmatt->data;
        loop_ub = indvec->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          repmatt_data[i1] = indvec_data[i1];
        }
        i1 = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = 4;
        sPWM->size[1] = repmatt->size[0];
        emxEnsureCapacity_real_T(sPWM, i1);
        sPWM_data = sPWM->data;
        loop_ub = repmatt->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          idx = (int)repmatt_data[i1] - 1;
          sPWM_data[4 * i1] = mat_data[idx];
          sPWM_data[4 * i1 + 1] = mat_data[idx + mat->size[0]];
          sPWM_data[4 * i1 + 2] = mat_data[idx + mat->size[0] * 2];
          sPWM_data[4 * i1 + 3] = mat_data[idx + mat->size[0] * 3];
        }
        /* 'PWM2kmers:43' sPWM2 = mat2(indvec,:).'; */
        i1 = sPWM2->size[0] * sPWM2->size[1];
        sPWM2->size[0] = 4;
        sPWM2->size[1] = repmatt->size[0];
        emxEnsureCapacity_real_T(sPWM2, i1);
        sPWM2_data = sPWM2->data;
        loop_ub = repmatt->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          idx = (int)repmatt_data[i1] - 1;
          sPWM2_data[4 * i1] = mat2_data[idx];
          sPWM2_data[4 * i1 + 1] = mat2_data[idx + mat2->size[0]];
          sPWM2_data[4 * i1 + 2] = mat2_data[idx + mat2->size[0] * 2];
          sPWM2_data[4 * i1 + 3] = mat2_data[idx + mat2->size[0] * 3];
        }
        /* 'PWM2kmers:44' ktree{1} = sPWM(:,1); */
        i1 = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree_data[0].f1, i1);
        /* 'PWM2kmers:45' ktree2{1} = sPWM2(:,1); */
        i1 = ktree2_data[0].f1->size[0];
        ktree2_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree2_data[0].f1, i1);
        ktree_data[0].f1->data[0] = sPWM_data[0];
        ktree2_data[0].f1->data[0] = sPWM2_data[0];
        ktree_data[0].f1->data[1] = sPWM_data[1];
        ktree2_data[0].f1->data[1] = sPWM2_data[1];
        ktree_data[0].f1->data[2] = sPWM_data[2];
        ktree2_data[0].f1->data[2] = sPWM2_data[2];
        ktree_data[0].f1->data[3] = sPWM_data[3];
        ktree2_data[0].f1->data[3] = sPWM2_data[3];
        /* 'PWM2kmers:46' for i3 = s(i2):k */
        i1 = (int)(k + (1.0 - s_data[b_i2]));
        for (i3 = 0; i3 < i1; i3++) {
          b_i3 = s_data[b_i2] + (double)i3;
          /* 'PWM2kmers:47' if loc(i3)==0 */
          if (indloc_data[(int)indvec_data[(int)b_i3 - 1] - 1] == 0.0) {
            /* 'PWM2kmers:48' if loc(i3-1) == 1 */
            if (indloc_data[(int)indvec_data[(int)(b_i3 - 1.0) - 1] - 1] ==
                1.0) {
              /* 'PWM2kmers:49' matt = mat(1,:)*p{indvec(i3)-indvec(i3-1)+1}; */
              for (i2 = 0; i2 < 16; i2++) {
                b_p[i2] = p_data[(int)(unsigned int)indvec_data[(int)b_i3 - 1] -
                                 (int)(unsigned int)
                                     indvec_data[(int)(b_i3 - 1.0) - 1]]
                              .f1[i2];
              }
              /* 'PWM2kmers:50' for i4 = 1:4 */
              c_tmp = pow(4.0, b_i3 - 1.0);
              for (j = 0; j < 4; j++) {
                idx = j << 2;
                matt = ((mat_data[0] * b_p[idx] +
                         mat_data[mat->size[0]] * b_p[idx + 1]) +
                        mat_data[mat->size[0] * 2] * b_p[idx + 2]) +
                       mat_data[mat->size[0] * 3] * b_p[idx + 3];
                b_matt[j] = matt;
                /* 'PWM2kmers:51' ktree{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4)) =
                 * ktree{i3-1}*matt(i4); */
                d = (((double)j + 1.0) - 1.0) * c_tmp + 1.0;
                d1 = c_tmp * ((double)j + 1.0);
                if (d > d1) {
                  i2 = -1;
                  ni = 0;
                } else {
                  i2 = (int)d - 2;
                  ni = (int)d1;
                }
                idx = (ni - i2) - 1;
                ni = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = idx;
                emxEnsureCapacity_real_T(b_ktree, ni);
                repmatt_data = b_ktree->data;
                for (ni = 0; ni < idx; ni++) {
                  repmatt_data[ni] =
                      ktree_data[(int)(b_i3 - 1.0) - 1].f1->data[ni] * matt;
                }
                loop_ub = b_ktree->size[1];
                for (ni = 0; ni < loop_ub; ni++) {
                  ktree_data[(int)b_i3 - 1].f1->data[(i2 + ni) + 1] =
                      repmatt_data[ni];
                }
                /* 'PWM2kmers:52' ktree2{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4))
                 * = ktree2{i3-1}*matt(i4); */
                if (d > d1) {
                  i2 = -1;
                  ni = 0;
                } else {
                  i2 = (int)d - 2;
                  ni = (int)d1;
                }
                idx = (ni - i2) - 1;
                ni = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = idx;
                emxEnsureCapacity_real_T(b_ktree, ni);
                repmatt_data = b_ktree->data;
                for (ni = 0; ni < idx; ni++) {
                  repmatt_data[ni] =
                      ktree2_data[(int)(b_i3 - 1.0) - 1].f1->data[ni] *
                      b_matt[j];
                }
                loop_ub = b_ktree->size[1];
                for (ni = 0; ni < loop_ub; ni++) {
                  ktree2_data[(int)b_i3 - 1].f1->data[(i2 + ni) + 1] =
                      repmatt_data[ni];
                }
              }
              /*                          ktree{i3} = reshape(ktree{i3-1}*matt,
               * [], 1); */
              /*                          ktree2{i3} =
               * reshape(ktree2{i3-1}*matt, [], 1); */
            } else {
              /* 'PWM2kmers:56' else */
              /* 'PWM2kmers:57' matt = p{indvec(i3)-indvec(i3-1)+1}; */
              /* 'PWM2kmers:58' repmatt = repelem(matt(:), 4^(i3-2)); */
              c_tmp = pow(4.0, b_i3 - 2.0);
              i2 = repmatt->size[0];
              repmatt->size[0] = (int)c_tmp << 4;
              emxEnsureCapacity_real_T(repmatt, i2);
              repmatt_data = repmatt->data;
              idx = -1;
              ni = (int)c_tmp;
              for (loop_ub = 0; loop_ub < 16; loop_ub++) {
                for (j = 0; j < ni; j++) {
                  repmatt_data[(idx + j) + 1] =
                      p_data[(int)(unsigned int)indvec_data[(int)b_i3 - 1] -
                             (int)(unsigned int)
                                 indvec_data[(int)(b_i3 - 1.0) - 1]]
                          .f1[loop_ub];
                }
                idx += (int)c_tmp;
              }
              /* 'PWM2kmers:59' ktree{i3} = repmat(ktree{i3-1}, 4, 1).*repmatt;
               */
              c_repmat(ktree_data[(int)(b_i3 - 1.0) - 1].f1, r);
              r1 = r->data;
              if (r->size[0] == repmatt->size[0]) {
                i2 = ktree_data[(int)b_i3 - 1].f1->size[0];
                ktree_data[(int)b_i3 - 1].f1->size[0] = r->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)b_i3 - 1].f1, i2);
                loop_ub = r->size[0];
                for (i2 = 0; i2 < loop_ub; i2++) {
                  ktree_data[(int)b_i3 - 1].f1->data[i2] =
                      r1[i2] * repmatt_data[i2];
                }
              } else {
                p_binary_expand_op(ktree, b_i3, r, repmatt);
                ktree_data = ktree->data;
              }
              /* 'PWM2kmers:60' ktree2{i3} = repmat(ktree2{i3-1}, 4,
               * 1).*repmatt; */
              c_repmat(ktree2_data[(int)(b_i3 - 1.0) - 1].f1, r);
              r1 = r->data;
              if (r->size[0] == repmatt->size[0]) {
                i2 = ktree2_data[(int)b_i3 - 1].f1->size[0];
                ktree2_data[(int)b_i3 - 1].f1->size[0] = r->size[0];
                emxEnsureCapacity_real_T(ktree2_data[(int)b_i3 - 1].f1, i2);
                loop_ub = r->size[0];
                for (i2 = 0; i2 < loop_ub; i2++) {
                  ktree2_data[(int)b_i3 - 1].f1->data[i2] =
                      r1[i2] * repmatt_data[i2];
                }
              } else {
                p_binary_expand_op(ktree2, b_i3, r, repmatt);
                ktree2_data = ktree2->data;
              }
            }
          } else {
            /* 'PWM2kmers:62' else */
            /* 'PWM2kmers:63' a = ktree{i3-1}.*sPWM(:,i3).'; */
            i2 = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, i2);
            repmatt_data = a->data;
            loop_ub = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
            for (i2 = 0; i2 < 4; i2++) {
              for (ni = 0; ni < loop_ub; ni++) {
                repmatt_data[ni + a->size[0] * i2] =
                    ktree_data[(int)(b_i3 - 1.0) - 1].f1->data[ni] *
                    sPWM_data[i2 + 4 * ((int)b_i3 - 1)];
              }
            }
            /* 'PWM2kmers:64' ktree{i3} = a(:); */
            i2 = ktree_data[(int)b_i3 - 1].f1->size[0];
            ktree_data[(int)b_i3 - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree_data[(int)b_i3 - 1].f1, i2);
            loop_ub = a->size[0] << 2;
            for (i2 = 0; i2 < loop_ub; i2++) {
              ktree_data[(int)b_i3 - 1].f1->data[i2] = repmatt_data[i2];
            }
            /* 'PWM2kmers:65' a = ktree2{i3-1}.*sPWM2(:,i3).'; */
            i2 = a->size[0] * a->size[1];
            a->size[0] = ktree2_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, i2);
            repmatt_data = a->data;
            loop_ub = ktree2_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
            for (i2 = 0; i2 < 4; i2++) {
              for (ni = 0; ni < loop_ub; ni++) {
                repmatt_data[ni + a->size[0] * i2] =
                    ktree2_data[(int)(b_i3 - 1.0) - 1].f1->data[ni] *
                    sPWM2_data[i2 + 4 * ((int)b_i3 - 1)];
              }
            }
            /* 'PWM2kmers:66' ktree2{i3} = a(:); */
            i2 = ktree2_data[(int)b_i3 - 1].f1->size[0];
            ktree2_data[(int)b_i3 - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree2_data[(int)b_i3 - 1].f1, i2);
            loop_ub = a->size[0] << 2;
            for (i2 = 0; i2 < loop_ub; i2++) {
              ktree2_data[(int)b_i3 - 1].f1->data[i2] = repmatt_data[i2];
            }
          }
        }
        /* 'PWM2kmers:69' if i2 <= rx */
        if (b_i2 + 1 <= rx) {
          /* 'PWM2kmers:70' for j = 1:X(i+1) */
          i1 = (int)X_data[b_i];
          for (j = 0; j < i1; j++) {
            /* 'PWM2kmers:71' if x(i2,j) ~= 0 */
            if (x_data[b_i2 + x->size[0] * j] != 0.0) {
              /* 'PWM2kmers:72' KC{ind(x(i2,j))} = KC{ind(x(i2,j))} + ktree{k} +
               * ktree2{k}; */
              i2 = r->size[0];
              r->size[0] =
                  KC_data[(int)
                              ind_data[(int)x_data[b_i2 + x->size[0] * j] - 1] -
                          1]
                      .f1->size[0];
              emxEnsureCapacity_real_T(r, i2);
              r1 = r->data;
              loop_ub =
                  KC_data[(int)
                              ind_data[(int)x_data[b_i2 + x->size[0] * j] - 1] -
                          1]
                      .f1->size[0];
              for (i2 = 0; i2 < loop_ub; i2++) {
                r1[i2] =
                    KC_data[(int)ind_data[(int)x_data[b_i2 + x->size[0] * j] -
                                          1] -
                            1]
                        .f1->data[i2];
              }
              if (r->size[0] == 1) {
                i2 = ktree_data[(int)k - 1].f1->size[0];
              } else {
                i2 = r->size[0];
              }
              if ((r->size[0] == ktree_data[(int)k - 1].f1->size[0]) &&
                  (i2 == ktree2_data[(int)k - 1].f1->size[0])) {
                i2 = KC_data[(int)ind_data[(int)x_data[b_i2 + x->size[0] * j] -
                                           1] -
                             1]
                         .f1->size[0];
                KC_data[(int)ind_data[(int)x_data[b_i2 + x->size[0] * j] - 1] -
                        1]
                    .f1->size[0] = r->size[0];
                emxEnsureCapacity_real_T(
                    KC_data[(int)ind_data[(int)x_data[b_i2 + x->size[0] * j] -
                                          1] -
                            1]
                        .f1,
                    i2);
                loop_ub = r->size[0];
                for (i2 = 0; i2 < loop_ub; i2++) {
                  KC_data[(int)
                              ind_data[(int)x_data[b_i2 + x->size[0] * j] - 1] -
                          1]
                      .f1->data[i2] =
                      (r1[i2] + ktree_data[(int)k - 1].f1->data[i2]) +
                      ktree2_data[(int)k - 1].f1->data[i2];
                }
              } else {
                q_binary_expand_op(KC, ind, x, b_i2, j, r, ktree, k, ktree2);
                KC_data = KC->data;
              }
              /* kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) =
               * kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) +
               * ktree{k}+ktree2{k}; */
            }
          }
        } else {
          /* 'PWM2kmers:76' else */
          /* 'PWM2kmers:77' KC{ind(i2)} = KC{ind(i2)} + ktree{k} + ktree2{k}; */
          if (KC_data[(int)ind_data[b_i2] - 1].f1->size[0] == 1) {
            idx = ktree_data[(int)k - 1].f1->size[0];
          } else {
            idx = KC_data[(int)ind_data[b_i2] - 1].f1->size[0];
          }
          if ((KC_data[(int)ind_data[b_i2] - 1].f1->size[0] ==
               ktree_data[(int)k - 1].f1->size[0]) &&
              (idx == ktree2_data[(int)k - 1].f1->size[0])) {
            loop_ub = KC_data[(int)ind_data[b_i2] - 1].f1->size[0];
            for (i1 = 0; i1 < loop_ub; i1++) {
              KC_data[(int)ind_data[b_i2] - 1].f1->data[i1] =
                  (KC_data[(int)ind_data[b_i2] - 1].f1->data[i1] +
                   ktree_data[(int)k - 1].f1->data[i1]) +
                  ktree2_data[(int)k - 1].f1->data[i1];
            }
          } else {
            r_binary_expand_op(KC, ind, b_i2, ktree, k, ktree2);
            KC_data = KC->data;
          }
          /* kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) =
           * kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) + ktree{k}+ktree2{k}; */
        }
      }
    }
  }
  emxFree_real_T(&r);
  emxFree_real_T(&a);
  emxFree_cell_wrap_13(&p);
  emxFree_real_T(&repmatt);
  emxFree_real_T(&sPWM2);
  emxFree_real_T(&sPWM);
  emxFree_real_T(&indvec);
  emxFree_uint32_T(&X);
  emxFree_cell_wrap_14(&ktree2);
  emxFree_cell_wrap_14(&ktree);
  emxFree_real_T(&mat2);
  /* 'PWM2kmers:83' for i = 1:length(c) */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    ni = 0;
  } else {
    idx = c->size[0];
    ni = c->size[1];
    if (idx >= ni) {
      ni = idx;
    }
  }
  for (b_i = 0; b_i < ni; b_i++) {
    /* 'PWM2kmers:84' kweig((4^k*(i-1)+1):4^k*i) = KC{i}; */
    c_tmp = n_tmp * (((double)b_i + 1.0) - 1.0) + 1.0;
    d = n_tmp * ((double)b_i + 1.0);
    if (c_tmp > d) {
      i = -1;
      i1 = 0;
    } else {
      i = (int)c_tmp - 2;
      i1 = (int)d;
    }
    loop_ub = (i1 - i) - 1;
    for (i1 = 0; i1 < loop_ub; i1++) {
      kweig_data[(i + i1) + 1] = KC_data[b_i].f1->data[i1];
    }
  }
  emxFree_cell_wrap_14(&KC);
  /* 'PWM2kmers:86' alen = length(c)-rcnum; */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    ni = 0;
  } else {
    idx = c->size[0];
    ni = c->size[1];
    if (idx >= ni) {
      ni = idx;
    }
  }
  /* 'PWM2kmers:87' kweig(4^k*alen+1:end) = kweig(4^k*alen+1:end)/sqrt(2); */
  c_tmp = n_tmp * ((double)ni - rcnum) + 1.0;
  if (c_tmp > kweig->size[0]) {
    i = 1;
    i1 = -1;
    i2 = 0;
  } else {
    i = (int)c_tmp;
    i1 = (int)c_tmp - 2;
    i2 = kweig->size[0];
  }
  idx = (i2 - i1) - 1;
  i2 = b_ktree->size[0] * b_ktree->size[1];
  b_ktree->size[0] = 1;
  b_ktree->size[1] = idx;
  emxEnsureCapacity_real_T(b_ktree, i2);
  repmatt_data = b_ktree->data;
  for (i2 = 0; i2 < idx; i2++) {
    repmatt_data[i2] = kweig_data[(i + i2) - 1] / 1.4142135623730951;
  }
  loop_ub = b_ktree->size[1];
  for (i = 0; i < loop_ub; i++) {
    kweig_data[(i1 + i) + 1] = repmatt_data[i];
  }
  emxFree_real_T(&b_ktree);
  /* 'PWM2kmers:88' kweig = kweig+kweig2; */
  if (kweig->size[0] != (int)n) {
    s_binary_expand_op(kweig, n);
  }
}

void p_binary_expand_op(emxArray_cell_wrap_14 *ktree2, double i3,
                        const emxArray_real_T *r1,
                        const emxArray_real_T *repmatt)
{
  cell_wrap_14 *ktree2_data;
  const double *r;
  const double *repmatt_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  repmatt_data = repmatt->data;
  r = r1->data;
  ktree2_data = ktree2->data;
  i = ktree2_data[(int)i3 - 1].f1->size[0];
  if (repmatt->size[0] == 1) {
    ktree2_data[(int)i3 - 1].f1->size[0] = r1->size[0];
  } else {
    ktree2_data[(int)i3 - 1].f1->size[0] = repmatt->size[0];
  }
  emxEnsureCapacity_real_T(ktree2_data[(int)i3 - 1].f1, i);
  stride_0_0 = (r1->size[0] != 1);
  stride_1_0 = (repmatt->size[0] != 1);
  if (repmatt->size[0] == 1) {
    loop_ub = r1->size[0];
  } else {
    loop_ub = repmatt->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    ktree2_data[(int)i3 - 1].f1->data[i] =
        r[i * stride_0_0] * repmatt_data[i * stride_1_0];
  }
}

/* End of code generation (PWM2kmers.c) */
