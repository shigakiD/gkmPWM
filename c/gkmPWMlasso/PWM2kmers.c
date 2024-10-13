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
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include "rot90.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void k_binary_expand_op(emxArray_cell_wrap_3 *KC,
                               const emxArray_real_T *ind,
                               const emxArray_real_T *x, int i2, int j,
                               const emxArray_real_T *r1,
                               const emxArray_cell_wrap_3 *ktree, double k,
                               const emxArray_cell_wrap_3 *ktree2);

static void l_binary_expand_op(emxArray_cell_wrap_3 *KC,
                               const emxArray_real_T *ind, int i2,
                               const emxArray_cell_wrap_3 *ktree, double k,
                               const emxArray_cell_wrap_3 *ktree2);

static void m_binary_expand_op(emxArray_real_T *kweig, double n);

/* Function Definitions */
static void k_binary_expand_op(emxArray_cell_wrap_3 *KC,
                               const emxArray_real_T *ind,
                               const emxArray_real_T *x, int i2, int j,
                               const emxArray_real_T *r1,
                               const emxArray_cell_wrap_3 *ktree, double k,
                               const emxArray_cell_wrap_3 *ktree2)
{
  const cell_wrap_3 *ktree2_data;
  const cell_wrap_3 *ktree_data;
  cell_wrap_3 *KC_data;
  const double *b_r;
  const double *ind_data;
  const double *x_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  int stride_2_0;
  ktree2_data = ktree2->data;
  ktree_data = ktree->data;
  b_r = r1->data;
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
        .f1->data[i] = (b_r[i * stride_0_0] +
                        ktree_data[(int)k - 1].f1->data[i * stride_1_0]) +
                       ktree2_data[(int)k - 1].f1->data[i * stride_2_0];
  }
}

static void l_binary_expand_op(emxArray_cell_wrap_3 *KC,
                               const emxArray_real_T *ind, int i2,
                               const emxArray_cell_wrap_3 *ktree, double k,
                               const emxArray_cell_wrap_3 *ktree2)
{
  const cell_wrap_3 *ktree2_data;
  const cell_wrap_3 *ktree_data;
  cell_wrap_3 *KC_data;
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

static void m_binary_expand_op(emxArray_real_T *kweig, double n)
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
  cell_wrap_14 *p_data;
  cell_wrap_3 *KC_data;
  cell_wrap_3 *ktree2_data;
  cell_wrap_3 *ktree_data;
  emxArray_cell_wrap_14 *p;
  emxArray_cell_wrap_3 *KC;
  emxArray_cell_wrap_3 *ktree;
  emxArray_cell_wrap_3 *ktree2;
  emxArray_real_T *a;
  emxArray_real_T *b_ktree;
  emxArray_real_T *b_r;
  emxArray_real_T *indvec;
  emxArray_real_T *mat2;
  emxArray_real_T *repmatt;
  emxArray_real_T *sPWM;
  emxArray_real_T *sPWM2;
  emxArray_uint32_T *X;
  double b_p[16];
  double matt[4];
  const double *c_data;
  const double *ind_data;
  const double *indloc_data;
  const double *mat_data;
  const double *s_data;
  const double *x_data;
  double M;
  double alen;
  double c_i3;
  double d;
  double d1;
  double d2;
  double m;
  double n;
  double n_tmp_tmp;
  double *indvec_data;
  double *kweig_data;
  double *mat2_data;
  double *r1;
  double *repmatt_data;
  double *sPWM2_data;
  double *sPWM_data;
  int b_i;
  int b_i2;
  int b_i3;
  int b_k;
  int i;
  int i1;
  int i2;
  int i3;
  int idx;
  int j;
  int ni;
  int rx;
  unsigned int *X_data;
  x_data = x->data;
  indloc_data = indloc->data;
  ind_data = ind->data;
  s_data = s->data;
  c_data = c->data;
  mat_data = mat->data;
  emxInit_cell_wrap_14(&p);
  /* makes gkm-pwm faster */
  /* 'PWM2kmers:3' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_14(p, i);
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
        i3 = i2 << 2;
        b_p[i1 + i3] = ((p_data[b_i].f1[i1] * negmat[i3] +
                         p_data[b_i].f1[i1 + 4] * negmat[i3 + 1]) +
                        p_data[b_i].f1[i1 + 8] * negmat[i3 + 2]) +
                       p_data[b_i].f1[i1 + 12] * negmat[i3 + 3];
      }
    }
    for (i1 = 0; i1 < 16; i1++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i1] = b_p[i1];
    }
  }
  /* 'PWM2kmers:9' m = length(x); */
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    ni = 0;
  } else {
    idx = x->size[0];
    ni = x->size[1];
    if (idx >= ni) {
      ni = idx;
    }
  }
  m = ni;
  /* 'PWM2kmers:10' M = length(mat)-l; */
  if ((mat->size[0] == 0) || (mat->size[1] == 0)) {
    ni = 0;
  } else {
    idx = mat->size[0];
    ni = mat->size[1];
    if (idx >= ni) {
      ni = idx;
    }
  }
  emxInit_real_T(&mat2, 2);
  M = (double)ni - l;
  /* 'PWM2kmers:11' n = 4^k*numel(c)/k; */
  n_tmp_tmp = pow(4.0, k);
  n = n_tmp_tmp * (double)(c->size[0] * c->size[1]) / k;
  /* number of possible k-mers */
  /* 'PWM2kmers:12' mat2 = rot90(mat,2); */
  rot90(mat, mat2);
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
  emxInit_cell_wrap_3(&ktree, 1);
  /* 'PWM2kmers:14' kweig2= zeros(n,1); */
  /* 'PWM2kmers:15' ktree = cell(k,1); */
  idx = (int)k;
  i = ktree->size[0];
  ktree->size[0] = (int)k;
  emxEnsureCapacity_cell_wrap_32(ktree, i);
  ktree_data = ktree->data;
  for (i = 0; i < idx; i++) {
    ktree_data[i].f1->size[0] = 0;
  }
  emxInit_cell_wrap_3(&ktree2, 1);
  /* 'PWM2kmers:16' ktree = coder.nullcopy(ktree); */
  /* 'PWM2kmers:17' ktree2 = cell(k,1); */
  i = ktree2->size[0];
  ktree2->size[0] = (int)k;
  emxEnsureCapacity_cell_wrap_32(ktree2, i);
  ktree2_data = ktree2->data;
  for (i = 0; i < idx; i++) {
    ktree2_data[i].f1->size[0] = 0;
  }
  /* 'PWM2kmers:18' ktree2 = coder.nullcopy(ktree2); */
  /* 'PWM2kmers:19' [rx,cx] = size(x); */
  rx = x->size[0];
  /* 'PWM2kmers:20' X = cx*ones(length(mat)-l+1,1); */
  if ((mat->size[0] == 0) || (mat->size[1] == 0)) {
    ni = 0;
  } else {
    idx = mat->size[0];
    ni = mat->size[1];
    if (idx >= ni) {
      ni = idx;
    }
  }
  emxInit_uint32_T(&X);
  ni = (int)(((double)ni - l) + 1.0);
  i = X->size[0];
  X->size[0] = ni;
  emxEnsureCapacity_uint32_T(X, i);
  X_data = X->data;
  for (i = 0; i < ni; i++) {
    X_data[i] = 1U;
  }
  ni = X->size[0];
  for (i = 0; i < ni; i++) {
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
  emxInit_cell_wrap_3(&KC, 1);
  i = KC->size[0];
  KC->size[0] = ni;
  emxEnsureCapacity_cell_wrap_32(KC, i);
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
    idx = (int)n_tmp_tmp;
    i = KC_data[b_i].f1->size[0];
    KC_data[b_i].f1->size[0] = (int)n_tmp_tmp;
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
    alen = pow(4.0, (double)b_i + 2.0);
    idx = (int)pow(4.0, (double)b_i + 2.0);
    i1 = ktree_data[b_i + 1].f1->size[0];
    ktree_data[b_i + 1].f1->size[0] = (int)alen;
    emxEnsureCapacity_real_T(ktree_data[b_i + 1].f1, i1);
    for (i1 = 0; i1 < idx; i1++) {
      ktree_data[b_i + 1].f1->data[i1] = 0.0;
    }
    /* 'PWM2kmers:31' ktree2{i} = zeros(4^i,1); */
    i1 = ktree2_data[b_i + 1].f1->size[0];
    ktree2_data[b_i + 1].f1->size[0] = (int)alen;
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
  emxInit_real_T(&b_r, 1);
  emxInit_real_T(&b_ktree, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:35' if i == M-1 */
    if (b_i == M - 1.0) {
      /* 'PWM2kmers:36' m = numel(c)/k; */
      m = (double)(c->size[0] * c->size[1]) / k;
    }
    /* 'PWM2kmers:38' for i2 = 1:m */
    i1 = (int)m;
    for (b_i2 = 0; b_i2 < i1; b_i2++) {
      /* 'PWM2kmers:39' if ~(i == M-1 && i2 > rx && i2 ~= m) */
      if ((b_i != M - 1.0) || (b_i2 + 1U <= (unsigned int)rx) ||
          ((double)b_i2 + 1.0 == m)) {
        /* 'PWM2kmers:40' indvec = c(i2,:)+i; */
        ni = c->size[1];
        i2 = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i2);
        indvec_data = indvec->data;
        for (i2 = 0; i2 < ni; i2++) {
          indvec_data[i2] = c_data[b_i2 + c->size[0] * i2] + (double)b_i;
        }
        /* 'PWM2kmers:41' loc = indloc(indvec); */
        /* 'PWM2kmers:42' sPWM = mat(indvec,:).'; */
        i2 = repmatt->size[0];
        repmatt->size[0] = indvec->size[1];
        emxEnsureCapacity_real_T(repmatt, i2);
        repmatt_data = repmatt->data;
        ni = indvec->size[1];
        for (i2 = 0; i2 < ni; i2++) {
          repmatt_data[i2] = indvec_data[i2];
        }
        ni = mat->size[1];
        i2 = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = mat->size[1];
        sPWM->size[1] = repmatt->size[0];
        emxEnsureCapacity_real_T(sPWM, i2);
        sPWM_data = sPWM->data;
        idx = repmatt->size[0];
        for (i2 = 0; i2 < idx; i2++) {
          for (i3 = 0; i3 < ni; i3++) {
            sPWM_data[i3 + sPWM->size[0] * i2] =
                mat_data[((int)repmatt_data[i2] + mat->size[0] * i3) - 1];
          }
        }
        /* 'PWM2kmers:43' sPWM2 = mat2(indvec,:).'; */
        ni = mat2->size[1];
        i2 = sPWM2->size[0] * sPWM2->size[1];
        sPWM2->size[0] = mat2->size[1];
        sPWM2->size[1] = repmatt->size[0];
        emxEnsureCapacity_real_T(sPWM2, i2);
        sPWM2_data = sPWM2->data;
        idx = repmatt->size[0];
        for (i2 = 0; i2 < idx; i2++) {
          for (i3 = 0; i3 < ni; i3++) {
            sPWM2_data[i3 + sPWM2->size[0] * i2] =
                mat2_data[((int)repmatt_data[i2] + mat2->size[0] * i3) - 1];
          }
        }
        /* 'PWM2kmers:44' ktree{1} = sPWM(:,1); */
        ni = sPWM->size[0];
        i2 = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = sPWM->size[0];
        emxEnsureCapacity_real_T(ktree_data[0].f1, i2);
        for (i2 = 0; i2 < ni; i2++) {
          ktree_data[0].f1->data[i2] = sPWM_data[i2];
        }
        /* 'PWM2kmers:45' ktree2{1} = sPWM2(:,1); */
        ni = sPWM2->size[0];
        i2 = ktree2_data[0].f1->size[0];
        ktree2_data[0].f1->size[0] = sPWM2->size[0];
        emxEnsureCapacity_real_T(ktree2_data[0].f1, i2);
        for (i2 = 0; i2 < ni; i2++) {
          ktree2_data[0].f1->data[i2] = sPWM2_data[i2];
        }
        /* 'PWM2kmers:46' for i3 = s(i2):k */
        i2 = (int)(k + (1.0 - s_data[b_i2]));
        for (b_i3 = 0; b_i3 < i2; b_i3++) {
          c_i3 = s_data[b_i2] + (double)b_i3;
          /* 'PWM2kmers:47' if loc(i3)==0 */
          if (indloc_data[(int)indvec_data[(int)c_i3 - 1] - 1] == 0.0) {
            /* 'PWM2kmers:48' if loc(i3-1) == 1 */
            if (indloc_data[(int)indvec_data[(int)(c_i3 - 1.0) - 1] - 1] ==
                1.0) {
              /* 'PWM2kmers:49' matt = mat(1,:)*p{indvec(i3)-indvec(i3-1)+1}; */
              for (i3 = 0; i3 < 16; i3++) {
                b_p[i3] = p_data[(int)(unsigned int)indvec_data[(int)c_i3 - 1] -
                                 (int)(unsigned int)
                                     indvec_data[(int)(c_i3 - 1.0) - 1]]
                              .f1[i3];
              }
              /* 'PWM2kmers:50' for i4 = 1:4 */
              d = pow(4.0, c_i3 - 1.0);
              for (b_k = 0; b_k < 4; b_k++) {
                idx = b_k << 2;
                alen = ((mat_data[0] * b_p[idx] +
                         mat_data[mat->size[0]] * b_p[idx + 1]) +
                        mat_data[mat->size[0] * 2] * b_p[idx + 2]) +
                       mat_data[mat->size[0] * 3] * b_p[idx + 3];
                matt[b_k] = alen;
                /* 'PWM2kmers:51' ktree{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4)) =
                 * ktree{i3-1}*matt(i4); */
                d1 = (((double)b_k + 1.0) - 1.0) * d + 1.0;
                d2 = d * ((double)b_k + 1.0);
                if (d1 > d2) {
                  i3 = -1;
                  j = 0;
                } else {
                  i3 = (int)d1 - 2;
                  j = (int)d2;
                }
                idx = (j - i3) - 1;
                j = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = idx;
                emxEnsureCapacity_real_T(b_ktree, j);
                repmatt_data = b_ktree->data;
                for (j = 0; j < idx; j++) {
                  repmatt_data[j] =
                      ktree_data[(int)(c_i3 - 1.0) - 1].f1->data[j] * alen;
                }
                ni = b_ktree->size[1];
                for (j = 0; j < ni; j++) {
                  ktree_data[(int)c_i3 - 1].f1->data[(i3 + j) + 1] =
                      repmatt_data[j];
                }
                /* 'PWM2kmers:52' ktree2{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4))
                 * = ktree2{i3-1}*matt(i4); */
                d1 = (((double)b_k + 1.0) - 1.0) * pow(4.0, c_i3 - 1.0) + 1.0;
                d2 = pow(4.0, c_i3 - 1.0) * ((double)b_k + 1.0);
                if (d1 > d2) {
                  i3 = -1;
                  j = 0;
                } else {
                  i3 = (int)d1 - 2;
                  j = (int)d2;
                }
                idx = (j - i3) - 1;
                j = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = idx;
                emxEnsureCapacity_real_T(b_ktree, j);
                repmatt_data = b_ktree->data;
                for (j = 0; j < idx; j++) {
                  repmatt_data[j] =
                      ktree2_data[(int)(c_i3 - 1.0) - 1].f1->data[j] *
                      matt[b_k];
                }
                ni = b_ktree->size[1];
                for (j = 0; j < ni; j++) {
                  ktree2_data[(int)c_i3 - 1].f1->data[(i3 + j) + 1] =
                      repmatt_data[j];
                }
              }
            } else {
              /* 'PWM2kmers:54' else */
              /* 'PWM2kmers:55' matt = p{indvec(i3)-indvec(i3-1)+1}; */
              /* 'PWM2kmers:56' repmatt = repelem(matt(:), 4^(i3-2)); */
              alen = pow(4.0, c_i3 - 2.0);
              i3 = repmatt->size[0];
              repmatt->size[0] = (int)alen << 4;
              emxEnsureCapacity_real_T(repmatt, i3);
              repmatt_data = repmatt->data;
              idx = -1;
              ni = (int)alen;
              for (b_k = 0; b_k < 16; b_k++) {
                for (j = 0; j < ni; j++) {
                  repmatt_data[(idx + j) + 1] =
                      p_data[(int)(unsigned int)indvec_data[(int)c_i3 - 1] -
                             (int)(unsigned int)
                                 indvec_data[(int)(c_i3 - 1.0) - 1]]
                          .f1[b_k];
                }
                idx += (int)alen;
              }
              /* 'PWM2kmers:57' ktree{i3} = repmat(ktree{i3-1}, 4, 1).*repmatt;
               */
              i3 = b_r->size[0];
              b_r->size[0] = ktree_data[(int)(c_i3 - 1.0) - 1].f1->size[0] << 2;
              emxEnsureCapacity_real_T(b_r, i3);
              r1 = b_r->data;
              idx = ktree_data[(int)(c_i3 - 1.0) - 1].f1->size[0];
              for (j = 0; j < 4; j++) {
                ni = j * idx;
                for (b_k = 0; b_k < idx; b_k++) {
                  r1[ni + b_k] =
                      ktree_data[(int)(c_i3 - 1.0) - 1].f1->data[b_k];
                }
              }
              if (b_r->size[0] == repmatt->size[0]) {
                i3 = ktree_data[(int)c_i3 - 1].f1->size[0];
                ktree_data[(int)c_i3 - 1].f1->size[0] = b_r->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)c_i3 - 1].f1, i3);
                ni = b_r->size[0];
                for (i3 = 0; i3 < ni; i3++) {
                  ktree_data[(int)c_i3 - 1].f1->data[i3] =
                      r1[i3] * repmatt_data[i3];
                }
              } else {
                j_binary_expand_op(ktree, c_i3, b_r, repmatt);
                ktree_data = ktree->data;
              }
              /* 'PWM2kmers:58' ktree2{i3} = repmat(ktree2{i3-1}, 4,
               * 1).*repmatt; */
              i3 = b_r->size[0];
              b_r->size[0] = ktree2_data[(int)(c_i3 - 1.0) - 1].f1->size[0]
                             << 2;
              emxEnsureCapacity_real_T(b_r, i3);
              r1 = b_r->data;
              idx = ktree2_data[(int)(c_i3 - 1.0) - 1].f1->size[0];
              for (j = 0; j < 4; j++) {
                ni = j * idx;
                for (b_k = 0; b_k < idx; b_k++) {
                  r1[ni + b_k] =
                      ktree2_data[(int)(c_i3 - 1.0) - 1].f1->data[b_k];
                }
              }
              if (b_r->size[0] == repmatt->size[0]) {
                i3 = ktree2_data[(int)c_i3 - 1].f1->size[0];
                ktree2_data[(int)c_i3 - 1].f1->size[0] = b_r->size[0];
                emxEnsureCapacity_real_T(ktree2_data[(int)c_i3 - 1].f1, i3);
                ni = b_r->size[0];
                for (i3 = 0; i3 < ni; i3++) {
                  ktree2_data[(int)c_i3 - 1].f1->data[i3] =
                      r1[i3] * repmatt_data[i3];
                }
              } else {
                j_binary_expand_op(ktree2, c_i3, b_r, repmatt);
                ktree2_data = ktree2->data;
              }
            }
          } else {
            /* 'PWM2kmers:60' else */
            /* 'PWM2kmers:61' a = ktree{i3-1}.*sPWM(:,i3).'; */
            ni = sPWM->size[0];
            i3 = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(c_i3 - 1.0) - 1].f1->size[0];
            a->size[1] = sPWM->size[0];
            emxEnsureCapacity_real_T(a, i3);
            repmatt_data = a->data;
            for (i3 = 0; i3 < ni; i3++) {
              idx = ktree_data[(int)(c_i3 - 1.0) - 1].f1->size[0];
              for (j = 0; j < idx; j++) {
                repmatt_data[j + a->size[0] * i3] =
                    ktree_data[(int)(c_i3 - 1.0) - 1].f1->data[j] *
                    sPWM_data[i3 + sPWM->size[0] * ((int)c_i3 - 1)];
              }
            }
            /* 'PWM2kmers:62' ktree{i3} = a(:); */
            i3 = ktree_data[(int)c_i3 - 1].f1->size[0];
            ktree_data[(int)c_i3 - 1].f1->size[0] = a->size[0] * a->size[1];
            emxEnsureCapacity_real_T(ktree_data[(int)c_i3 - 1].f1, i3);
            ni = a->size[0] * a->size[1];
            for (i3 = 0; i3 < ni; i3++) {
              ktree_data[(int)c_i3 - 1].f1->data[i3] = repmatt_data[i3];
            }
            /* 'PWM2kmers:63' a = ktree2{i3-1}.*sPWM2(:,i3).'; */
            ni = sPWM2->size[0];
            i3 = a->size[0] * a->size[1];
            a->size[0] = ktree2_data[(int)(c_i3 - 1.0) - 1].f1->size[0];
            a->size[1] = sPWM2->size[0];
            emxEnsureCapacity_real_T(a, i3);
            repmatt_data = a->data;
            for (i3 = 0; i3 < ni; i3++) {
              idx = ktree2_data[(int)(c_i3 - 1.0) - 1].f1->size[0];
              for (j = 0; j < idx; j++) {
                repmatt_data[j + a->size[0] * i3] =
                    ktree2_data[(int)(c_i3 - 1.0) - 1].f1->data[j] *
                    sPWM2_data[i3 + sPWM2->size[0] * ((int)c_i3 - 1)];
              }
            }
            /* 'PWM2kmers:64' ktree2{i3} = a(:); */
            i3 = ktree2_data[(int)c_i3 - 1].f1->size[0];
            ktree2_data[(int)c_i3 - 1].f1->size[0] = a->size[0] * a->size[1];
            emxEnsureCapacity_real_T(ktree2_data[(int)c_i3 - 1].f1, i3);
            ni = a->size[0] * a->size[1];
            for (i3 = 0; i3 < ni; i3++) {
              ktree2_data[(int)c_i3 - 1].f1->data[i3] = repmatt_data[i3];
            }
          }
        }
        /* 'PWM2kmers:67' if i2 <= rx */
        if (b_i2 + 1U <= (unsigned int)rx) {
          /* 'PWM2kmers:68' for j = 1:X(i+1) */
          i2 = (int)X_data[b_i];
          for (j = 0; j < i2; j++) {
            /* 'PWM2kmers:69' if x(i2,j) ~= 0 */
            if (x_data[b_i2 + x->size[0] * j] != 0.0) {
              /* 'PWM2kmers:70' KC{ind(x(i2,j))} = KC{ind(x(i2,j))} + ktree{k} +
               * ktree2{k}; */
              i3 = b_r->size[0];
              b_r->size[0] =
                  KC_data[(int)
                              ind_data[(int)x_data[b_i2 + x->size[0] * j] - 1] -
                          1]
                      .f1->size[0];
              emxEnsureCapacity_real_T(b_r, i3);
              r1 = b_r->data;
              ni =
                  KC_data[(int)
                              ind_data[(int)x_data[b_i2 + x->size[0] * j] - 1] -
                          1]
                      .f1->size[0];
              for (i3 = 0; i3 < ni; i3++) {
                r1[i3] =
                    KC_data[(int)ind_data[(int)x_data[b_i2 + x->size[0] * j] -
                                          1] -
                            1]
                        .f1->data[i3];
              }
              if (b_r->size[0] == 1) {
                i3 = ktree_data[(int)k - 1].f1->size[0];
              } else {
                i3 = b_r->size[0];
              }
              if ((b_r->size[0] == ktree_data[(int)k - 1].f1->size[0]) &&
                  (i3 == ktree2_data[(int)k - 1].f1->size[0])) {
                i3 = KC_data[(int)ind_data[(int)x_data[b_i2 + x->size[0] * j] -
                                           1] -
                             1]
                         .f1->size[0];
                KC_data[(int)ind_data[(int)x_data[b_i2 + x->size[0] * j] - 1] -
                        1]
                    .f1->size[0] = b_r->size[0];
                emxEnsureCapacity_real_T(
                    KC_data[(int)ind_data[(int)x_data[b_i2 + x->size[0] * j] -
                                          1] -
                            1]
                        .f1,
                    i3);
                ni = b_r->size[0];
                for (i3 = 0; i3 < ni; i3++) {
                  KC_data[(int)
                              ind_data[(int)x_data[b_i2 + x->size[0] * j] - 1] -
                          1]
                      .f1->data[i3] =
                      (r1[i3] + ktree_data[(int)k - 1].f1->data[i3]) +
                      ktree2_data[(int)k - 1].f1->data[i3];
                }
              } else {
                k_binary_expand_op(KC, ind, x, b_i2, j, b_r, ktree, k, ktree2);
                KC_data = KC->data;
              }
              /* kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) =
               * kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) +
               * ktree{k}+ktree2{k}; */
            }
          }
        } else {
          /* 'PWM2kmers:74' else */
          /* 'PWM2kmers:75' KC{ind(i2)} = KC{ind(i2)} + ktree{k} + ktree2{k}; */
          if (KC_data[(int)ind_data[b_i2] - 1].f1->size[0] == 1) {
            idx = ktree_data[(int)k - 1].f1->size[0];
          } else {
            idx = KC_data[(int)ind_data[b_i2] - 1].f1->size[0];
          }
          if ((KC_data[(int)ind_data[b_i2] - 1].f1->size[0] ==
               ktree_data[(int)k - 1].f1->size[0]) &&
              (idx == ktree2_data[(int)k - 1].f1->size[0])) {
            ni = KC_data[(int)ind_data[b_i2] - 1].f1->size[0];
            for (i2 = 0; i2 < ni; i2++) {
              KC_data[(int)ind_data[b_i2] - 1].f1->data[i2] =
                  (KC_data[(int)ind_data[b_i2] - 1].f1->data[i2] +
                   ktree_data[(int)k - 1].f1->data[i2]) +
                  ktree2_data[(int)k - 1].f1->data[i2];
            }
          } else {
            l_binary_expand_op(KC, ind, b_i2, ktree, k, ktree2);
            KC_data = KC->data;
          }
          /* kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) =
           * kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) + ktree{k}+ktree2{k}; */
        }
      }
    }
  }
  emxFree_real_T(&b_r);
  emxFree_real_T(&a);
  emxFree_cell_wrap_14(&p);
  emxFree_real_T(&repmatt);
  emxFree_real_T(&sPWM2);
  emxFree_real_T(&sPWM);
  emxFree_real_T(&indvec);
  emxFree_uint32_T(&X);
  emxFree_cell_wrap_3(&ktree2);
  emxFree_cell_wrap_3(&ktree);
  emxFree_real_T(&mat2);
  /* 'PWM2kmers:81' for i = 1:numel(c)/k */
  i = (int)((double)(c->size[0] * c->size[1]) / k);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:82' kweig((4^k*(i-1)+1):4^k*i) = KC{i}; */
    d = n_tmp_tmp * (((double)b_i + 1.0) - 1.0) + 1.0;
    d1 = n_tmp_tmp * ((double)b_i + 1.0);
    if (d > d1) {
      i1 = -1;
      i2 = 0;
    } else {
      i1 = (int)d - 2;
      i2 = (int)d1;
    }
    ni = (i2 - i1) - 1;
    for (i2 = 0; i2 < ni; i2++) {
      kweig_data[(i1 + i2) + 1] = KC_data[b_i].f1->data[i2];
    }
  }
  emxFree_cell_wrap_3(&KC);
  /* 'PWM2kmers:84' alen = numel(c)/k-rcnum; */
  alen = (double)(c->size[0] * c->size[1]) / k - rcnum;
  /* 'PWM2kmers:85' kweig(4^k*alen+1:end) = kweig(4^k*alen+1:end)/sqrt(2); */
  d = n_tmp_tmp * alen + 1.0;
  if (d > kweig->size[0]) {
    i = 1;
  } else {
    i = (int)d;
  }
  d = pow(4.0, k) * alen + 1.0;
  if (d > kweig->size[0]) {
    i1 = -1;
    i2 = 0;
  } else {
    i1 = (int)d - 2;
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
  ni = b_ktree->size[1];
  for (i = 0; i < ni; i++) {
    kweig_data[(i1 + i) + 1] = repmatt_data[i];
  }
  emxFree_real_T(&b_ktree);
  /* 'PWM2kmers:86' kweig = kweig+kweig2; */
  if (kweig->size[0] != (int)n) {
    m_binary_expand_op(kweig, n);
  }
}

void j_binary_expand_op(emxArray_cell_wrap_3 *ktree2, double i3,
                        const emxArray_real_T *r1,
                        const emxArray_real_T *repmatt)
{
  cell_wrap_3 *ktree2_data;
  const double *b_r;
  const double *repmatt_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  repmatt_data = repmatt->data;
  b_r = r1->data;
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
        b_r[i * stride_0_0] * repmatt_data[i * stride_1_0];
  }
}

/* End of code generation (PWM2kmers.c) */
