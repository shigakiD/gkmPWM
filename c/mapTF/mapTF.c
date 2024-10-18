/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF.c
 *
 * Code generation for function 'mapTF'
 *
 */

/* Include files */
#include "mapTF.h"
#include "applyScalarFunction.h"
#include "blockedSummation.h"
#include "colon.h"
#include "div.h"
#include "eml_setop.h"
#include "feof.h"
#include "fgetl.h"
#include "fgets.h"
#include "fileManager.h"
#include "fileread.h"
#include "find.h"
#include "fliplr.h"
#include "fseek.h"
#include "ftell.h"
#include "genIndex.h"
#include "getmotif.h"
#include "mapTF_data.h"
#include "mapTF_emxutil.h"
#include "mapTF_initialize.h"
#include "mapTF_rtwutil.h"
#include "mapTF_types.h"
#include "minOrMax.h"
#include "mod.h"
#include "mtimes.h"
#include "print_processing.h"
#include "prod.h"
#include "repmat.h"
#include "rot90.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "spdiags.h"
#include "std.h"
#include "str2double.h"
#include "strip.h"
#include "strtok.h"
#include "sum.h"
#include "tic.h"
#include "toc.h"
#include "cblas.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Variable Definitions */
static const char cv[11] = { '_', 'm', 'o', 't', 'i', 'f', 's', '.', 'o', 'u',
  't' };

/* Function Declarations */
static void MAPTF(const emxArray_real_T *ss, emxArray_real_T *pwm_prob, double
                  l_svm, const emxArray_real_T *LEN, const emxArray_real_T
                  *LEN_2, const emxArray_real_T *shift, const emxArray_real_T
                  *gkmprob, const emxArray_cell_wrap_4 *names, double a,
                  emxArray_real_T *Lmat, emxArray_cell_wrap_4 *NAME);
static void PWM_corr(const emxArray_char_T *fn, const emxArray_cell_wrap_0 *VV,
                     const emxArray_cell_wrap_7 *NN, const emxArray_cell_wrap_5 *
                     LL, const emxArray_cell_wrap_4 *seq, double dsvmcut);
static double countGC(const emxArray_cell_wrap_3 *s);
static void d_binary_expand_op(emxArray_real_T *mat, const emxArray_real_T *D,
  int i, const emxArray_real_T *path, const emxArray_real_T *pwm_prob);
static void getMOTIF(const emxArray_char_T *fn, emxArray_cell_wrap_5 *mat,
                     emxArray_cell_wrap_4 *names, emxArray_real_T *len);
static void getdenovomotif(const emxArray_char_T *filename, emxArray_cell_wrap_6
  *mat, emxArray_real_T *w);
static void i_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *mat,
  const emxArray_real_T *r1);
static void letterconvert(const emxArray_char_T *s, emxArray_real_T *en);
static void minus(emxArray_real_T *f, const emxArray_real_T *minnorm);
static void plus(emxArray_real_T *pwm_prob, const emxArray_real_T *b);
static void ppmsim(emxArray_cell_wrap_6 *mot, const emxArray_real_T *lenvec,
                   double *ind, double *M);
static void process_motifs(const emxArray_char_T *dfn, const emxArray_char_T
  *lfn, const emxArray_char_T *memefn, const emxArray_char_T *ofn, double
  PWMcorrcut, double Negative);
static void scoreseqkmer(const emxArray_cell_wrap_5 *PWM2, const
  emxArray_cell_wrap_5 *lPWM2, const emxArray_real_T *Lmat, const
  emxArray_real_T *ss, const emxArray_cell_wrap_6 *Smat, double l_svm, const
  emxArray_real_T *dsvm, emxArray_real_T *varscore);
static void seq2pv(const emxArray_char_T *sfn, const emxArray_char_T *wfn,
                   double l_svm, double Negative, emxArray_cell_wrap_0 *P,
                   emxArray_cell_wrap_1 *V, emxArray_cell_wrap_2 *seqindmat,
                   emxArray_cell_wrap_3 *seqout, emxArray_cell_wrap_4 *seq);
static void trim_pwm(emxArray_cell_wrap_6 *p, emxArray_real_T *info,
                     emxArray_real_T *len);

/* Function Definitions */
/*
 * function [Lmat, NAME] = MAPTF(fn, ss, pwm_prob, l_svm, k_svm, LEN, LEN_2, shift, gkmprob, names, a, b)
 */
static void MAPTF(const emxArray_real_T *ss, emxArray_real_T *pwm_prob, double
                  l_svm, const emxArray_real_T *LEN, const emxArray_real_T
                  *LEN_2, const emxArray_real_T *shift, const emxArray_real_T
                  *gkmprob, const emxArray_cell_wrap_4 *names, double a,
                  emxArray_real_T *Lmat, emxArray_cell_wrap_4 *NAME)
{
  const cell_wrap_4 *names_data;
  cell_wrap_4 *NAME_data;
  emxArray_boolean_T *b_path;
  emxArray_int32_T *f;
  emxArray_int32_T *ind;
  emxArray_real_T *C;
  emxArray_real_T *C2;
  emxArray_real_T *D;
  emxArray_real_T *L2;
  emxArray_real_T *R;
  emxArray_real_T *b;
  emxArray_real_T *b_C;
  emxArray_real_T *b_L2;
  emxArray_real_T *c;
  emxArray_real_T *mat;
  emxArray_real_T *neg;
  emxArray_real_T *path;
  const double *LEN_2_data;
  const double *LEN_data;
  const double *gkmprob_data;
  const double *shift_data;
  double L;
  double b_pwm_prob;
  double name_len;
  double *C2_data;
  double *C_data;
  double *D_data;
  double *Lmat_data;
  double *R_data;
  double *mat_data;
  double *neg_data;
  double *path_data;
  double *pwm_prob_data;
  int b_i;
  int b_loop_ub;
  int i;
  int i1;
  int i2;
  int ibmat;
  int k;
  int loop_ub;
  int nx;
  int *ind_data;
  bool *b_path_data;
  names_data = names->data;
  gkmprob_data = gkmprob->data;
  shift_data = shift->data;
  LEN_2_data = LEN_2->data;
  LEN_data = LEN->data;
  pwm_prob_data = pwm_prob->data;
  emxInit_real_T(&mat, 2);

  /* 'mapTF:277' L = length(ss)-l_svm+1; */
  L = ((double)ss->size[1] - l_svm) + 1.0;

  /* 'mapTF:278' n = sum(LEN)+1; */
  name_len = blockedSummation(LEN, LEN->size[0]);

  /* 'mapTF:279' mat = zeros(n,L)-Inf; */
  i = mat->size[0] * mat->size[1];
  mat->size[0] = (int)(name_len + 1.0);
  loop_ub = (int)L;
  mat->size[1] = (int)L;
  emxEnsureCapacity_real_T(mat, i);
  mat_data = mat->data;
  nx = (int)(name_len + 1.0) * (int)L;
  for (i = 0; i < nx; i++) {
    mat_data[i] = rtMinusInf;
  }

  emxInit_int32_T(&ind, 2);

  /* 'mapTF:280' ind = zeros(n,L); */
  i = ind->size[0] * ind->size[1];
  ind->size[0] = (int)(name_len + 1.0);
  ind->size[1] = (int)L;
  emxEnsureCapacity_int32_T(ind, i);
  ind_data = ind->data;
  for (i = 0; i < nx; i++) {
    ind_data[i] = 0;
  }

  emxInit_real_T(&C, 1);

  /* 'mapTF:281' LEN = [0;LEN]; */
  i = C->size[0];
  C->size[0] = LEN->size[0] + 1;
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  C_data[0] = 0.0;
  b_loop_ub = LEN->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    C_data[i + 1] = LEN_data[i];
  }

  /* 'mapTF:282' C = cumsum(LEN); */
  if (C->size[0] != 1) {
    i = C->size[0];
    for (k = 0; k <= i - 2; k++) {
      C_data[k + 1] += C_data[k];
    }
  }

  /* 'mapTF:283' D = setdiff(1:C(end),C+1)'; */
  emxInit_real_T(&R, 2);
  R_data = R->data;
  if (rtIsNaN(C_data[C->size[0] - 1])) {
    i = R->size[0] * R->size[1];
    R->size[0] = 1;
    R->size[1] = 1;
    emxEnsureCapacity_real_T(R, i);
    R_data = R->data;
    R_data[0] = rtNaN;
  } else if (C_data[C->size[0] - 1] < 1.0) {
    R->size[0] = 1;
    R->size[1] = 0;
  } else if (rtIsInf(C_data[C->size[0] - 1]) && (1.0 == C_data[C->size[0] - 1]))
  {
    i = R->size[0] * R->size[1];
    R->size[0] = 1;
    R->size[1] = 1;
    emxEnsureCapacity_real_T(R, i);
    R_data = R->data;
    R_data[0] = rtNaN;
  } else {
    i = R->size[0] * R->size[1];
    R->size[0] = 1;
    R->size[1] = (int)floor(C_data[C->size[0] - 1] - 1.0) + 1;
    emxEnsureCapacity_real_T(R, i);
    R_data = R->data;
    b_loop_ub = (int)floor(C_data[C->size[0] - 1] - 1.0);
    for (i = 0; i <= b_loop_ub; i++) {
      R_data[i] = (double)i + 1.0;
    }
  }

  emxInit_real_T(&c, 2);
  emxInit_real_T(&b_C, 1);
  i = b_C->size[0];
  b_C->size[0] = C->size[0];
  emxEnsureCapacity_real_T(b_C, i);
  Lmat_data = b_C->data;
  b_loop_ub = C->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    Lmat_data[i] = C_data[i] + 1.0;
  }

  emxInit_real_T(&D, 1);
  emxInit_int32_T(&f, 1);
  do_vectors(R, b_C, c, f, &nx);
  Lmat_data = c->data;
  i = D->size[0];
  D->size[0] = c->size[1];
  emxEnsureCapacity_real_T(D, i);
  D_data = D->data;
  b_loop_ub = c->size[1];
  for (i = 0; i < b_loop_ub; i++) {
    D_data[i] = Lmat_data[i];
  }

  emxFree_real_T(&c);

  /* 'mapTF:284' C2 = [C(2:end);n]; */
  if (2 > C->size[0]) {
    i = 0;
    i1 = 0;
  } else {
    i = 1;
    i1 = C->size[0];
  }

  emxInit_real_T(&C2, 1);
  b_loop_ub = i1 - i;
  i1 = C2->size[0];
  C2->size[0] = b_loop_ub + 1;
  emxEnsureCapacity_real_T(C2, i1);
  C2_data = C2->data;
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    C2_data[i1] = C_data[i + i1];
  }

  emxInit_real_T(&neg, 1);
  C2_data[b_loop_ub] = name_len + 1.0;

  /* 'mapTF:285' pos = log(gkmprob); */
  /* 'mapTF:286' neg = log(1-gkmprob); */
  i = neg->size[0];
  neg->size[0] = gkmprob->size[0];
  emxEnsureCapacity_real_T(neg, i);
  neg_data = neg->data;
  b_loop_ub = gkmprob->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    neg_data[i] = 1.0 - gkmprob_data[i];
  }

  nx = neg->size[0];
  for (k = 0; k < nx; k++) {
    neg_data[k] = log(neg_data[k]);
  }

  emxInit_real_T(&path, 1);

  /* 'mapTF:287' pwm_prob = pwm_prob+repmat(pos',n-1,1); */
  i = path->size[0];
  path->size[0] = gkmprob->size[0];
  emxEnsureCapacity_real_T(path, i);
  path_data = path->data;
  b_loop_ub = gkmprob->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    path_data[i] = gkmprob_data[i];
  }

  nx = gkmprob->size[0];
  for (k = 0; k < nx; k++) {
    path_data[k] = log(path_data[k]);
  }

  emxInit_real_T(&b, 2);
  i = (int)((name_len + 1.0) - 1.0);
  i1 = b->size[0] * b->size[1];
  b->size[0] = (int)((name_len + 1.0) - 1.0);
  b->size[1] = path->size[0];
  emxEnsureCapacity_real_T(b, i1);
  Lmat_data = b->data;
  nx = path->size[0];
  for (b_loop_ub = 0; b_loop_ub < nx; b_loop_ub++) {
    ibmat = b_loop_ub * (int)((name_len + 1.0) - 1.0);
    for (k = 0; k < i; k++) {
      Lmat_data[ibmat + k] = path_data[b_loop_ub];
    }
  }

  if ((pwm_prob->size[0] == b->size[0]) && (pwm_prob->size[1] == b->size[1])) {
    b_loop_ub = pwm_prob->size[0] * pwm_prob->size[1];
    for (i = 0; i < b_loop_ub; i++) {
      pwm_prob_data[i] += Lmat_data[i];
    }
  } else {
    plus(pwm_prob, b);
    pwm_prob_data = pwm_prob->data;
  }

  /* 'mapTF:288' for i = 1:a */
  i = (int)a;
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:289' mat(C(i)+1,1) = pwm_prob(C(i)+1,1); */
    nx = (int)(C_data[b_i] + 1.0) - 1;
    mat_data[nx] = pwm_prob_data[nx];
  }

  /* 'mapTF:291' mat(n,1) = neg(1); */
  mat_data[(int)(name_len + 1.0) - 1] = neg_data[0];

  /* 'mapTF:292' ind(D,:) = repmat(D-1,1,L); */
  i1 = path->size[0];
  path->size[0] = D->size[0];
  emxEnsureCapacity_real_T(path, i1);
  path_data = path->data;
  b_loop_ub = D->size[0];
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    path_data[i1] = D_data[i1] - 1.0;
  }

  i1 = b->size[0] * b->size[1];
  b->size[0] = path->size[0];
  b->size[1] = (int)L;
  emxEnsureCapacity_real_T(b, i1);
  Lmat_data = b->data;
  nx = path->size[0];
  for (b_loop_ub = 0; b_loop_ub < loop_ub; b_loop_ub++) {
    ibmat = b_loop_ub * nx;
    for (k = 0; k < nx; k++) {
      Lmat_data[ibmat + k] = path_data[k];
    }
  }

  b_loop_ub = b->size[1];
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    nx = b->size[0];
    for (i2 = 0; i2 < nx; i2++) {
      ind_data[((int)D_data[i2] + ind->size[0] * i1) - 1] = (int)Lmat_data[i2 +
        b->size[0] * i1];
    }
  }

  emxFree_real_T(&b);

  /* 'mapTF:293' for i = 2:L */
  i1 = (int)(L + -1.0);
  for (b_i = 0; b_i < i1; b_i++) {
    /* 'mapTF:294' for j = 1:a */
    for (k = 0; k < i; k++) {
      /* 'mapTF:295' [mat(C(j)+1,i),ind(C(j)+1,i)] = max(mat(C2,i-1)+pwm_prob(C(j)+1,i)); */
      ibmat = (int)(C_data[k] + 1.0) - 1;
      b_pwm_prob = pwm_prob_data[ibmat + pwm_prob->size[0] * (b_i + 1)];
      i2 = b_C->size[0];
      b_C->size[0] = C2->size[0];
      emxEnsureCapacity_real_T(b_C, i2);
      Lmat_data = b_C->data;
      b_loop_ub = C2->size[0];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        Lmat_data[i2] = mat_data[((int)C2_data[i2] + mat->size[0] * b_i) - 1] +
          b_pwm_prob;
      }

      d_maximum(b_C, &b_pwm_prob, &nx);
      mat_data[((int)(C_data[k] + 1.0) + mat->size[0] * (b_i + 1)) - 1] =
        b_pwm_prob;
      ind_data[ibmat + ind->size[0] * (b_i + 1)] = nx;

      /* 'mapTF:296' ind(C(j)+1,i) = C2(ind(C(j)+1,i)); */
      ind_data[ibmat + ind->size[0] * (b_i + 1)] = (int)C2_data[ind_data[ibmat +
        ind->size[0] * (b_i + 1)] - 1];
    }

    /* 'mapTF:298' mat(D,i) = mat(D-1,i-1)+pwm_prob(D,i); */
    b_loop_ub = D->size[0];
    i2 = path->size[0];
    path->size[0] = D->size[0];
    emxEnsureCapacity_real_T(path, i2);
    path_data = path->data;
    for (i2 = 0; i2 < b_loop_ub; i2++) {
      path_data[i2] = mat_data[((int)D_data[i2] + mat->size[0] * b_i) - 2];
    }

    if (path->size[0] == D->size[0]) {
      b_loop_ub = path->size[0];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        nx = (int)D_data[i2] - 1;
        mat_data[nx + mat->size[0] * (b_i + 1)] = path_data[i2] +
          pwm_prob_data[nx + pwm_prob->size[0] * (b_i + 1)];
      }
    } else {
      d_binary_expand_op(mat, D, b_i, path, pwm_prob);
      mat_data = mat->data;
    }

    /* 'mapTF:299' [mat(n,i),ind(n,i)] = max(mat(C2,i-1)+neg(i)); */
    i2 = b_C->size[0];
    b_C->size[0] = C2->size[0];
    emxEnsureCapacity_real_T(b_C, i2);
    Lmat_data = b_C->data;
    b_loop_ub = C2->size[0];
    for (i2 = 0; i2 < b_loop_ub; i2++) {
      Lmat_data[i2] = mat_data[((int)C2_data[i2] + mat->size[0] * b_i) - 1] +
        neg_data[b_i + 1];
    }

    d_maximum(b_C, &b_pwm_prob, &nx);
    mat_data[((int)(name_len + 1.0) + mat->size[0] * (b_i + 1)) - 1] =
      b_pwm_prob;
    ind_data[((int)(name_len + 1.0) + ind->size[0] * (b_i + 1)) - 1] = nx;

    /* 'mapTF:300' ind(n,i) = C2(ind(n,i)); */
    ind_data[((int)(name_len + 1.0) + ind->size[0] * (b_i + 1)) - 1] = (int)
      C2_data[ind_data[((int)(name_len + 1.0) + ind->size[0] * (b_i + 1)) - 1] -
      1];
  }

  emxFree_real_T(&neg);
  emxFree_real_T(&C2);
  emxFree_real_T(&D);
  emxFree_real_T(&mat);

  /* 'mapTF:302' path = zeros(L,1); */
  i = path->size[0];
  path->size[0] = (int)L;
  emxEnsureCapacity_real_T(path, i);
  path_data = path->data;
  for (i = 0; i < loop_ub; i++) {
    path_data[i] = 0.0;
  }

  /* 'mapTF:303' path(end) = n; */
  path_data[(int)L - 1] = name_len + 1.0;

  /* 'mapTF:305' for i = L-1:-1:1 */
  i = (int)-((-1.0 - (L - 1.0)) + 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    name_len = (L - 1.0) + -(double)b_i;

    /* 'mapTF:306' path(i) = ind(path(i+1),i+1); */
    path_data[(int)name_len - 1] = ind_data[((int)path_data[(int)(unsigned int)
      name_len] + ind->size[0] * (int)(unsigned int)name_len) - 1];
  }

  emxFree_int32_T(&ind);
  emxInit_real_T(&L2, 2);
  C2_data = L2->data;

  /* 'mapTF:309' total_len = length(ss); */
  /* 'mapTF:310' for i = 1:a */
  /* 'mapTF:323' L2 = []; */
  L2->size[0] = 0;
  L2->size[1] = 4;

  /* 'mapTF:324' for i = 1:length(LEN_2) */
  i = LEN_2->size[0];
  emxInit_boolean_T(&b_path, 1);
  emxInit_real_T(&b_L2, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:325' f = find(path==C(i)+1); */
    loop_ub = path->size[0];
    i1 = b_path->size[0];
    b_path->size[0] = path->size[0];
    emxEnsureCapacity_boolean_T(b_path, i1);
    b_path_data = b_path->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_path_data[i1] = (path_data[i1] == C_data[b_i] + 1.0);
    }

    eml_find(b_path, f);
    ind_data = f->data;

    /* 'mapTF:326' if ~isempty(f) */
    if (f->size[0] != 0) {
      /* 'mapTF:327' for j = 1:length(f) */
      i1 = f->size[0];
      for (k = 0; k < i1; k++) {
        /* 'mapTF:328' if f(j)+shift(i)+LEN_2(i)< length(ss)-l_svm+1 && f(j)+shift(i) > l_svm-1 */
        i2 = ind_data[k];
        name_len = (double)i2 + shift_data[b_i];
        if ((name_len + LEN_2_data[b_i] < ((double)ss->size[1] - l_svm) + 1.0) &&
            (name_len > l_svm - 1.0)) {
          /* 'mapTF:329' [~,ind] = max(gkmprob(f(j):f(j)+2*shift(i)+LEN_2(i)-l_svm)); */
          name_len = (((double)i2 + 2.0 * shift_data[b_i]) + LEN_2_data[b_i]) -
            l_svm;
          if (i2 > name_len) {
            i2 = 0;
            nx = 0;
          } else {
            i2 = ind_data[k] - 1;
            nx = (int)name_len;
          }

          loop_ub = nx - i2;
          nx = b_C->size[0];
          b_C->size[0] = loop_ub;
          emxEnsureCapacity_real_T(b_C, nx);
          Lmat_data = b_C->data;
          for (nx = 0; nx < loop_ub; nx++) {
            Lmat_data[nx] = gkmprob_data[i2 + nx];
          }

          d_maximum(b_C, &name_len, &nx);

          /* 'mapTF:330' if f(j)-1+ind < 3 */
          i2 = ind_data[k];
          name_len = ((double)i2 - 1.0) + (double)nx;
          if (name_len < 3.0) {
            /* 'mapTF:331' R = 1:5; */
            nx = R->size[0] * R->size[1];
            R->size[0] = 1;
            R->size[1] = 5;
            emxEnsureCapacity_real_T(R, nx);
            R_data = R->data;
            for (nx = 0; nx < 5; nx++) {
              R_data[nx] = (double)nx + 1.0;
            }
          } else if (name_len > L - 2.0) {
            /* 'mapTF:332' elseif f(j)-1+ind > L-2 */
            /* 'mapTF:333' R = L-4:L; */
            if (L < L - 4.0) {
              R->size[0] = 1;
              R->size[1] = 0;
            } else if ((rtIsInf(L - 4.0) || rtIsInf(L)) && (L - 4.0 == L)) {
              nx = R->size[0] * R->size[1];
              R->size[0] = 1;
              R->size[1] = 1;
              emxEnsureCapacity_real_T(R, nx);
              R_data = R->data;
              R_data[0] = rtNaN;
            } else if (L - 4.0 == L - 4.0) {
              nx = R->size[0] * R->size[1];
              R->size[0] = 1;
              loop_ub = (int)floor(L - (L - 4.0));
              R->size[1] = loop_ub + 1;
              emxEnsureCapacity_real_T(R, nx);
              R_data = R->data;
              for (nx = 0; nx <= loop_ub; nx++) {
                R_data[nx] = (L - 4.0) + (double)nx;
              }
            } else {
              eml_float_colon(L - 4.0, L, R);
              R_data = R->data;
            }
          } else {
            /* 'mapTF:334' else */
            /* 'mapTF:335' R = f(j)+ind-3:f(j)+ind+1; */
            name_len = (double)i2 + (double)nx;
            if (name_len + 1.0 < name_len - 3.0) {
              R->size[0] = 1;
              R->size[1] = 0;
            } else if (name_len - 3.0 == name_len - 3.0) {
              nx = R->size[0] * R->size[1];
              R->size[0] = 1;
              loop_ub = (int)((name_len + 1.0) - (name_len - 3.0));
              R->size[1] = loop_ub + 1;
              emxEnsureCapacity_real_T(R, nx);
              R_data = R->data;
              for (nx = 0; nx <= loop_ub; nx++) {
                R_data[nx] = (name_len - 3.0) + (double)nx;
              }
            } else {
              eml_float_colon(name_len - 3.0, name_len + 1.0, R);
              R_data = R->data;
            }
          }

          /* 'mapTF:337' L2 = [L2; f(j)+shift(i) f(j)+shift(i)+LEN_2(i)-1 mean(gkmprob(R)) i]; */
          nx = b_C->size[0];
          b_C->size[0] = R->size[1];
          emxEnsureCapacity_real_T(b_C, nx);
          Lmat_data = b_C->data;
          loop_ub = R->size[1];
          for (nx = 0; nx < loop_ub; nx++) {
            Lmat_data[nx] = gkmprob_data[(int)R_data[nx] - 1];
          }

          loop_ub = L2->size[0];
          nx = b_L2->size[0] * b_L2->size[1];
          b_L2->size[0] = L2->size[0] + 1;
          b_L2->size[1] = 4;
          emxEnsureCapacity_real_T(b_L2, nx);
          Lmat_data = b_L2->data;
          for (nx = 0; nx < 4; nx++) {
            for (ibmat = 0; ibmat < loop_ub; ibmat++) {
              Lmat_data[ibmat + b_L2->size[0] * nx] = C2_data[ibmat + L2->size[0]
                * nx];
            }
          }

          name_len = (double)i2 + shift_data[b_i];
          Lmat_data[L2->size[0]] = name_len;
          Lmat_data[L2->size[0] + b_L2->size[0]] = (name_len + LEN_2_data[b_i])
            - 1.0;
          Lmat_data[L2->size[0] + b_L2->size[0] * 2] = blockedSummation(b_C,
            R->size[1]) / (double)R->size[1];
          Lmat_data[L2->size[0] + b_L2->size[0] * 3] = (double)b_i + 1.0;
          i2 = L2->size[0] * L2->size[1];
          L2->size[0] = b_L2->size[0];
          L2->size[1] = 4;
          emxEnsureCapacity_real_T(L2, i2);
          C2_data = L2->data;
          loop_ub = b_L2->size[0] * 4;
          for (i2 = 0; i2 < loop_ub; i2++) {
            C2_data[i2] = Lmat_data[i2];
          }
        }
      }
    }
  }

  emxFree_boolean_T(&b_path);
  emxFree_real_T(&b_C);
  emxFree_real_T(&R);
  emxFree_real_T(&C);

  /* 'mapTF:342' name_len = numel(L2)/4; */
  name_len = (double)(L2->size[0] << 2) / 4.0;

  /* 'mapTF:343' NAME = cell(name_len,1); */
  /* 'mapTF:344' Lmat = zeros(name_len,4); */
  i = (int)name_len;
  i1 = Lmat->size[0] * Lmat->size[1];
  Lmat->size[0] = (int)name_len;
  Lmat->size[1] = 4;
  emxEnsureCapacity_real_T(Lmat, i1);
  Lmat_data = Lmat->data;
  loop_ub = (int)name_len << 2;
  for (i1 = 0; i1 < loop_ub; i1++) {
    Lmat_data[i1] = 0.0;
  }

  /* 'mapTF:345' for i=1:name_len */
  i1 = NAME->size[0];
  NAME->size[0] = (int)name_len;
  emxEnsureCapacity_cell_wrap_4(NAME, i1);
  NAME_data = NAME->data;
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:346' NAME{i} = ''; */
    NAME_data[b_i].f1->size[0] = 1;
    NAME_data[b_i].f1->size[1] = 0;
  }

  /* 'mapTF:348' if ~isempty(L2) */
  if (L2->size[0] != 0) {
    /* 'mapTF:349' [~,b] = sort(L2(:,1)); */
    loop_ub = L2->size[0];
    i1 = path->size[0];
    path->size[0] = L2->size[0];
    emxEnsureCapacity_real_T(path, i1);
    path_data = path->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      path_data[i1] = C2_data[i1];
    }

    sort(path, f);
    ind_data = f->data;

    /* 'mapTF:350' L2 = L2(b,:); */
    i1 = b_L2->size[0] * b_L2->size[1];
    b_L2->size[0] = f->size[0];
    b_L2->size[1] = 4;
    emxEnsureCapacity_real_T(b_L2, i1);
    Lmat_data = b_L2->data;
    loop_ub = f->size[0];
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        Lmat_data[i2 + b_L2->size[0] * i1] = C2_data[(ind_data[i2] + L2->size[0]
          * i1) - 1];
      }
    }

    i1 = L2->size[0] * L2->size[1];
    L2->size[0] = b_L2->size[0];
    L2->size[1] = 4;
    emxEnsureCapacity_real_T(L2, i1);
    C2_data = L2->data;
    loop_ub = b_L2->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      C2_data[i1] = Lmat_data[i1];
    }

    /* 'mapTF:351' for i = 1:name_len */
    i1 = NAME->size[0];
    NAME->size[0] = (int)name_len;
    emxEnsureCapacity_cell_wrap_4(NAME, i1);
    NAME_data = NAME->data;
    i1 = Lmat->size[0] * Lmat->size[1];
    Lmat->size[0] = (int)name_len;
    Lmat->size[1] = 4;
    emxEnsureCapacity_real_T(Lmat, i1);
    Lmat_data = Lmat->data;
    for (b_i = 0; b_i < i; b_i++) {
      /* 'mapTF:352' NAME{i} = names{L2(i,4)}; */
      i1 = NAME_data[b_i].f1->size[0] * NAME_data[b_i].f1->size[1];
      NAME_data[b_i].f1->size[0] = 1;
      NAME_data[b_i].f1->size[1] = names_data[(int)C2_data[b_i + L2->size[0] * 3]
        - 1].f1->size[1];
      emxEnsureCapacity_char_T(NAME_data[b_i].f1, i1);
      loop_ub = names_data[(int)C2_data[b_i + L2->size[0] * 3] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        NAME_data[b_i].f1->data[i1] = names_data[(int)C2_data[b_i + L2->size[0] *
          3] - 1].f1->data[i1];
      }

      /* 'mapTF:353' Lmat(i,:) = [L2(i,4) L2(i,1) L2(i,2) L2(i,3)]; */
      Lmat_data[b_i] = C2_data[b_i + L2->size[0] * 3];
      Lmat_data[b_i + Lmat->size[0]] = C2_data[b_i];
      Lmat_data[b_i + Lmat->size[0] * 2] = C2_data[b_i + L2->size[0]];
      Lmat_data[b_i + Lmat->size[0] * 3] = C2_data[b_i + L2->size[0] * 2];
    }
  }

  emxFree_real_T(&b_L2);
  emxFree_real_T(&L2);
  emxFree_int32_T(&f);
  emxFree_real_T(&path);
}

/*
 * function PWM_corr(fn, VV, NN, LL, seq, dsvmcut)
 */
static void PWM_corr(const emxArray_char_T *fn, const emxArray_cell_wrap_0 *VV,
                     const emxArray_cell_wrap_7 *NN, const emxArray_cell_wrap_5 *
                     LL, const emxArray_cell_wrap_4 *seq, double dsvmcut)
{
  static const char b_cv[19] = { '_', 'T', 'F', 'B', 'S', '_', 'l', 'o', 'c',
    'a', 't', 'i', 'o', 'n', 's', '.', 'o', 'u', 't' };

  FILE* b_NULL;
  FILE* filestar;
  const cell_wrap_0 *VV_data;
  const cell_wrap_4 *seq_data;
  const cell_wrap_5 *LL_data;
  const cell_wrap_7 *NN_data;
  emxArray_char_T *s;
  emxArray_char_T *varargin_2;
  emxArray_char_T *varargin_8;
  double d;
  double d1;
  int b_i;
  int i;
  int i1;
  int i2;
  int i3;
  int j;
  int loop_ub;
  int loop_ub_tmp;
  unsigned int unnamed_idx_1;
  const char *fn_data;
  signed char fileid;
  char *s_data;
  char *varargin_2_data;
  char *varargin_8_data;
  bool autoflush;
  seq_data = seq->data;
  LL_data = LL->data;
  NN_data = NN->data;
  VV_data = VV->data;
  fn_data = fn->data;
  emxInit_char_T(&s, 2);

  /* 'mapTF:456' n = length(VV); */
  /* 'mapTF:457' fid1 = fopen([fn '_TFBS_locations.out'],'w'); */
  i = s->size[0] * s->size[1];
  s->size[0] = 1;
  s->size[1] = fn->size[1] + 19;
  emxEnsureCapacity_char_T(s, i);
  s_data = s->data;
  loop_ub = fn->size[1];
  for (i = 0; i < loop_ub; i++) {
    s_data[i] = fn_data[i];
  }

  for (i = 0; i < 19; i++) {
    s_data[i + fn->size[1]] = b_cv[i];
  }

  fileid = cfopen(s, "wb");

  /* 'mapTF:458' for j = 1:n */
  i = VV->size[0];
  emxInit_char_T(&varargin_2, 2);
  emxInit_char_T(&varargin_8, 2);
  for (j = 0; j < i; j++) {
    /* 'mapTF:459' NAME = NN{j}; */
    /* 'mapTF:460' Lmat = LL{j}; */
    /* 'mapTF:461' varscore = VV{j}; */
    /* 'mapTF:462' s = seq{j*2}; */
    i1 = s->size[0] * s->size[1];
    s->size[0] = 1;
    loop_ub_tmp = ((j + 1) << 1) - 1;
    loop_ub = seq_data[loop_ub_tmp].f1->size[1];
    s->size[1] = loop_ub;
    emxEnsureCapacity_char_T(s, i1);
    s_data = s->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      s_data[i1] = seq_data[loop_ub_tmp].f1->data[i1];
    }

    /* 'mapTF:463' L = length(NAME); */
    /* 'mapTF:464' for i = 1:L */
    i1 = NN_data[j].f1->size[0];
    for (b_i = 0; b_i < i1; b_i++) {
      /* 'mapTF:465' if varscore(i) > dsvmcut */
      if (VV_data[j].f1->data[b_i] > dsvmcut) {
        /* 'mapTF:466' fprintf(fid1, '%d\t%s\t%d\t%d\t%d\t%f\t%f\t%s\n', int32(j), NAME{i}, ... */
        /* 'mapTF:467'                 int32(Lmat(i,1)), int32(Lmat(i,2)), ... */
        /* 'mapTF:468'                 int32(Lmat(i,3)), Lmat(i,4), varscore(i), s(Lmat(i,2):Lmat(i,3))); */
        d = LL_data[j].f1->data[b_i + LL_data[j].f1->size[0]];
        d1 = LL_data[j].f1->data[b_i + LL_data[j].f1->size[0] * 2];
        if (d > d1) {
          loop_ub_tmp = -1;
          i2 = 0;
        } else {
          loop_ub_tmp = (int)d - 2;
          i2 = (int)d1;
        }

        i3 = varargin_2->size[0] * varargin_2->size[1];
        varargin_2->size[0] = 1;
        varargin_2->size[1] = NN_data[j].f1->data[b_i].f1->size[1] + 1;
        emxEnsureCapacity_char_T(varargin_2, i3);
        varargin_2_data = varargin_2->data;
        loop_ub = NN_data[j].f1->data[b_i].f1->size[1];
        for (i3 = 0; i3 < loop_ub; i3++) {
          varargin_2_data[i3] = NN_data[j].f1->data[b_i].f1->data[i3];
        }

        varargin_2_data[NN_data[j].f1->data[b_i].f1->size[1]] = '\x00';
        unnamed_idx_1 = (unsigned int)((i2 - loop_ub_tmp) - 1);
        i2 = varargin_8->size[0] * varargin_8->size[1];
        varargin_8->size[0] = 1;
        varargin_8->size[1] = (int)unnamed_idx_1 + 1;
        emxEnsureCapacity_char_T(varargin_8, i2);
        varargin_8_data = varargin_8->data;
        loop_ub = (int)unnamed_idx_1;
        for (i2 = 0; i2 < loop_ub; i2++) {
          varargin_8_data[i2] = s_data[(loop_ub_tmp + i2) + 1];
        }

        varargin_8_data[(int)unnamed_idx_1] = '\x00';
        b_NULL = NULL;
        getfilestar(fileid, &filestar, &autoflush);
        if (!(filestar == b_NULL)) {
          fprintf(filestar, "%d\t%s\t%d\t%d\t%d\t%f\t%f\t%s\n", j + 1,
                  &varargin_2_data[0], (int)rt_roundd_snf(LL_data[j].f1->
                   data[b_i]), (int)rt_roundd_snf(LL_data[j].f1->data[b_i +
                   LL_data[j].f1->size[0]]), (int)rt_roundd_snf(LL_data[j]
                   .f1->data[b_i + LL_data[j].f1->size[0] * 2]), LL_data[j]
                  .f1->data[b_i + LL_data[j].f1->size[0] * 3], VV_data[j]
                  .f1->data[b_i], &varargin_8_data[0]);
          if (autoflush) {
            fflush(filestar);
          }
        }
      }
    }
  }

  emxFree_char_T(&varargin_8);
  emxFree_char_T(&varargin_2);
  emxFree_char_T(&s);

  /* 'mapTF:472' fclose(fid1); */
  cfclose(fileid);
}

/*
 * function GC = countGC(s)
 */
static double countGC(const emxArray_cell_wrap_3 *s)
{
  const cell_wrap_3 *s_data;
  double GC;
  double L;
  double d;
  int b_i;
  int i;
  int i1;
  int j;
  s_data = s->data;

  /* 'mapTF:637' GC = 0; */
  GC = 0.0;

  /* 'mapTF:638' L = 0; */
  L = 0.0;

  /* 'mapTF:639' for i = 1:length(s) */
  i = s->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:640' S = s{i}; */
    /* 'mapTF:641' N = length(S); */
    /* 'mapTF:642' L = N + L; */
    i1 = s_data[b_i].f1->size[1];
    L += (double)s_data[b_i].f1->size[1];

    /* 'mapTF:643' for j = 1:N */
    for (j = 0; j < i1; j++) {
      /* 'mapTF:644' if S(j) == 2 || S(j) == 3 */
      d = s_data[b_i].f1->data[j];
      if ((d == 2.0) || (d == 3.0)) {
        /* 'mapTF:645' GC = GC + 1; */
        GC++;
      }
    }
  }

  /* 'mapTF:649' GC = GC/L; */
  GC /= L;
  return GC;
}

static void d_binary_expand_op(emxArray_real_T *mat, const emxArray_real_T *D,
  int i, const emxArray_real_T *path, const emxArray_real_T *pwm_prob)
{
  const double *D_data;
  const double *path_data;
  const double *pwm_prob_data;
  double *mat_data;
  int b_i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  pwm_prob_data = pwm_prob->data;
  path_data = path->data;
  D_data = D->data;
  mat_data = mat->data;
  stride_0_0 = (path->size[0] != 1);
  stride_1_0 = (D->size[0] != 1);
  if (D->size[0] == 1) {
    loop_ub = path->size[0];
  } else {
    loop_ub = D->size[0];
  }

  for (b_i = 0; b_i < loop_ub; b_i++) {
    mat_data[((int)D_data[b_i] + mat->size[0] * (i + 1)) - 1] = path_data[b_i *
      stride_0_0] + pwm_prob_data[((int)D_data[b_i * stride_1_0] +
      pwm_prob->size[0] * (i + 1)) - 1];
  }
}

/*
 * function [mat, names, len] = getMOTIF(fn)
 */
static void getMOTIF(const emxArray_char_T *fn, emxArray_cell_wrap_5 *mat,
                     emxArray_cell_wrap_4 *names, emxArray_real_T *len)
{
  static const char b[5] = { 'M', 'O', 'T', 'I', 'F' };

  static const char b_cv[3] = { 'a', 'l', 'l' };

  FILE* b_NULL;
  FILE* filestar;
  int st;
  int wherefrom;
  long position_t;
  cell_wrap_4 *names_data;
  cell_wrap_5 *mat_data;
  emxArray_char_T *b_remain;
  emxArray_char_T *line;
  emxArray_char_T *remain;
  emxArray_char_T *v1;
  emxArray_char_T *v2;
  emxArray_real_T *tmp;
  creal_T dc;
  creal_T dc1;
  creal_T dc2;
  creal_T dc3;
  double curr_pos;
  double idx;
  double *len_data;
  double *tmp_data;
  unsigned int b_i;
  int exitg1;
  int fid;
  int i;
  unsigned int line_count;
  int ret;
  char a[5];
  const char *fn_data;
  signed char fileid;
  char *line_data;
  bool b_bool;
  fn_data = fn->data;

  /* 'mapTF:407' fid = fopen(fn); */
  b_bool = false;
  if (fn->size[1] == 3) {
    ret = 0;
    do {
      exitg1 = 0;
      if (ret < 3) {
        if (fn_data[ret] != b_cv[ret]) {
          exitg1 = 1;
        } else {
          ret++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    fid = 0;
  } else {
    fileid = cfopen(fn, "rb");
    fid = fileid;
  }

  /* 'mapTF:408' if fid < 0 */
  if (fid < 0) {
    /* 'mapTF:409' fprintf('ERROR: Cannot open gkmPWM motif files\n'); */
    printf("ERROR: Cannot open gkmPWM motif files\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:411' curr_pos = ftell(fid); */
  getfilestar(fid, &filestar, &b_bool);
  if ((fid == 0) || (fid == 1) || (fid == 2)) {
    filestar = NULL;
  }

  if (filestar == NULL) {
    curr_pos = -1.0;
  } else {
    position_t = ftell(filestar);
    curr_pos = (double)position_t;
  }

  /* 'mapTF:412' idx=0; */
  idx = 0.0;

  /* 'mapTF:413' while ~feof(fid) */
  b_NULL = NULL;
  emxInit_char_T(&line, 2);
  do {
    exitg1 = 0;
    getfilestar(fid, &filestar, &b_bool);
    if (filestar == b_NULL) {
      ret = 0;
    } else {
      st = feof(filestar);
      ret = ((int)st != 0);
    }

    if (ret == 0) {
      /* 'mapTF:414' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF:415' if length(line) >= 5 && strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF:416' idx=idx+1; */
          idx++;
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /* 'mapTF:419' fseek(fid, curr_pos, 'bof'); */
  wherefrom = SEEK_SET;
  if ((!rtIsInf(curr_pos)) && (!rtIsNaN(curr_pos)) && (floor(curr_pos) ==
       curr_pos)) {
    getfilestar(fid, &filestar, &b_bool);
    if ((fid == 0) || (fid == 1) || (fid == 2)) {
      filestar = NULL;
    }

    if (!(filestar == NULL)) {
      fseek(filestar, (long int)curr_pos, wherefrom);
    }
  }

  /* 'mapTF:421' mat = cell(idx, 1); */
  ret = (int)idx;
  i = mat->size[0];
  mat->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_5(mat, i);
  mat_data = mat->data;

  /* 'mapTF:422' mat = coder.nullcopy(mat); */
  /* 'mapTF:423' names = cell(idx, 1); */
  i = names->size[0];
  names->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_4(names, i);
  names_data = names->data;
  for (i = 0; i < ret; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 4;
    names_data[i].f1->size[0] = 1;
    names_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:424' names = coder.nullcopy(names); */
  /* 'mapTF:425' len = zeros(idx, 1); */
  i = len->size[0];
  len->size[0] = (int)idx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;
  for (i = 0; i < ret; i++) {
    len_data[i] = 0.0;
  }

  /* 'mapTF:427' i=0; */
  b_i = 0U;

  /* 'mapTF:428' while ~feof(fid) */
  b_NULL = NULL;
  emxInit_real_T(&tmp, 2);
  emxInit_char_T(&v1, 2);
  emxInit_char_T(&remain, 2);
  emxInit_char_T(&v2, 2);
  emxInit_char_T(&b_remain, 2);
  do {
    exitg1 = 0;
    getfilestar(fid, &filestar, &b_bool);
    if (filestar == b_NULL) {
      ret = 0;
    } else {
      st = feof(filestar);
      ret = ((int)st != 0);
    }

    if (ret == 0) {
      /* 'mapTF:429' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF:430' if length(line) < 5 || ~strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF:433' i = i + 1; */
          b_i++;

          /* 'mapTF:434' names{i} = fgetl(fid); */
          fgetl(fid, names_data[(int)b_i - 1].f1);

          /* 'mapTF:435' len(i) = real(str2double(fgetl(fid))); */
          fgetl(fid, line);
          dc = str2double(line);
          len_data[(int)b_i - 1] = dc.re;

          /* 'mapTF:436' tmp = zeros(len(i), 4); */
          i = (int)len_data[(int)b_i - 1];
          ret = tmp->size[0] * tmp->size[1];
          tmp->size[0] = i;
          tmp->size[1] = 4;
          emxEnsureCapacity_real_T(tmp, ret);
          tmp_data = tmp->data;
          ret = i << 2;
          for (i = 0; i < ret; i++) {
            tmp_data[i] = 0.0;
          }

          /* 'mapTF:437' line_count = 1; */
          line_count = 1U;

          /* 'mapTF:438' line = fgetl(fid); */
          fgetl(fid, line);

          /* 'mapTF:439' while ~isempty(line) */
          while (line->size[1] != 0) {
            /* 'mapTF:440' [v1, remain] = strtok(line); */
            b_strtok(line, v1, remain);

            /* 'mapTF:441' [v2, remain] = strtok(remain); */
            b_strtok(remain, v2, b_remain);

            /* 'mapTF:442' [v3, v4] = strtok(remain); */
            b_strtok(b_remain, line, remain);

            /* 'mapTF:443' vals = [real(str2double(v1)),  */
            /* 'mapTF:444'                 real(str2double(v2)),  */
            /* 'mapTF:445'                 real(str2double(v3)),  */
            /* 'mapTF:446'                 real(str2double(v4))]'; */
            dc = str2double(v1);
            dc1 = str2double(v2);
            dc2 = str2double(line);
            dc3 = str2double(remain);
            tmp_data[(int)line_count - 1] = dc.re;
            tmp_data[((int)line_count + tmp->size[0]) - 1] = dc1.re;
            tmp_data[((int)line_count + tmp->size[0] * 2) - 1] = dc2.re;
            tmp_data[((int)line_count + tmp->size[0] * 3) - 1] = dc3.re;

            /* 'mapTF:447' tmp(line_count,:) = vals; */
            /* 'mapTF:448' line = fgetl(fid); */
            fgetl(fid, line);

            /* 'mapTF:449' line_count = line_count + 1; */
            line_count++;
          }

          /* 'mapTF:451' mat{i} = tmp; */
          i = mat_data[(int)b_i - 1].f1->size[0] * mat_data[(int)b_i - 1]
            .f1->size[1];
          mat_data[(int)b_i - 1].f1->size[0] = tmp->size[0];
          mat_data[(int)b_i - 1].f1->size[1] = 4;
          emxEnsureCapacity_real_T(mat_data[(int)b_i - 1].f1, i);
          ret = tmp->size[0] * 4;
          for (i = 0; i < ret; i++) {
            mat_data[(int)b_i - 1].f1->data[i] = tmp_data[i];
          }
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_remain);
  emxFree_char_T(&v2);
  emxFree_char_T(&remain);
  emxFree_char_T(&v1);
  emxFree_real_T(&tmp);
  emxFree_char_T(&line);

  /* 'mapTF:453' fclose(fid); */
  cfclose(fid);
}

/*
 * function [mat,w] = getdenovomotif(filename)
 */
static void getdenovomotif(const emxArray_char_T *filename, emxArray_cell_wrap_6
  *mat, emxArray_real_T *w)
{
  static const char cv1[6] = { 'w', 'e', 'i', 'g', 'h', 't' };

  static const char b[5] = { 'M', 'O', 'T', 'I', 'F' };

  static const char b_cv[3] = { 'a', 'l', 'l' };

  FILE* b_NULL;
  FILE* filestar;
  int st;
  int wherefrom;
  long position_t;
  cell_wrap_6 *mat_data;
  emxArray_char_T *a__12;
  emxArray_char_T *a__9;
  emxArray_char_T *b_remain;
  emxArray_char_T *curr_alphabet;
  emxArray_char_T *curr_w;
  emxArray_char_T *line;
  emxArray_char_T *remain;
  emxArray_int32_T *match_out;
  emxArray_int32_T *matches;
  emxArray_real_T *tmp;
  creal_T dc;
  creal_T dc1;
  creal_T dc2;
  creal_T dc3;
  double vals[4];
  double curr_pos;
  double idx;
  double *tmp_data;
  double *w_data;
  int out_size[2];
  unsigned int b_i;
  int c_i;
  int exitg1;
  int fid;
  int i;
  int j;
  unsigned int line_count;
  int ret;
  int *match_out_data;
  int *matches_data;
  char a[5];
  const char *filename_data;
  signed char fileid;
  char *line_data;
  bool b_bool;
  filename_data = filename->data;

  /*  filename is the meme file that contains the motifs; */
  /*  n is the nth motif in the file */
  /* 'mapTF:934' fid = fopen(filename); */
  b_bool = false;
  if (filename->size[1] == 3) {
    ret = 0;
    do {
      exitg1 = 0;
      if (ret < 3) {
        if (filename_data[ret] != b_cv[ret]) {
          exitg1 = 1;
        } else {
          ret++;
        }
      } else {
        b_bool = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  if (b_bool) {
    fid = 0;
  } else {
    fileid = cfopen(filename, "rb");
    fid = fileid;
  }

  /* 'mapTF:935' if fid < 0 */
  if (fid < 0) {
    /* 'mapTF:936' fprintf('ERROR: Cannot open gkmPWM motif files\n'); */
    printf("ERROR: Cannot open gkmPWM motif files\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:938' curr_pos = ftell(fid); */
  getfilestar(fid, &filestar, &b_bool);
  if ((fid == 0) || (fid == 1) || (fid == 2)) {
    filestar = NULL;
  }

  if (filestar == NULL) {
    curr_pos = -1.0;
  } else {
    position_t = ftell(filestar);
    curr_pos = (double)position_t;
  }

  /* 'mapTF:939' idx=0; */
  idx = 0.0;

  /* 'mapTF:940' while ~feof(fid) */
  b_NULL = NULL;
  emxInit_char_T(&line, 2);
  do {
    exitg1 = 0;
    getfilestar(fid, &filestar, &b_bool);
    if (filestar == b_NULL) {
      ret = 0;
    } else {
      st = feof(filestar);
      ret = ((int)st != 0);
    }

    if (ret == 0) {
      /* 'mapTF:941' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF:942' if length(line) >= 5 && strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF:943' idx=idx+1; */
          idx++;
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /* 'mapTF:946' fseek(fid, curr_pos, 'bof'); */
  wherefrom = SEEK_SET;
  if ((!rtIsInf(curr_pos)) && (!rtIsNaN(curr_pos)) && (floor(curr_pos) ==
       curr_pos)) {
    getfilestar(fid, &filestar, &b_bool);
    if ((fid == 0) || (fid == 1) || (fid == 2)) {
      filestar = NULL;
    }

    if (!(filestar == NULL)) {
      fseek(filestar, (long int)curr_pos, wherefrom);
    }
  }

  /* 'mapTF:948' mat = cell(idx, 1); */
  ret = (int)idx;
  i = mat->size[0];
  mat->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_6(mat, i);
  mat_data = mat->data;
  for (i = 0; i < ret; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:949' mat = coder.nullcopy(mat); */
  /* 'mapTF:950' w = zeros(idx, 1); */
  i = w->size[0];
  w->size[0] = (int)idx;
  emxEnsureCapacity_real_T(w, i);
  w_data = w->data;
  for (i = 0; i < ret; i++) {
    w_data[i] = 0.0;
  }

  /* 'mapTF:952' i = 0; */
  b_i = 0U;

  /* 'mapTF:953' while ~feof(fid) */
  b_NULL = NULL;
  emxInit_char_T(&a__9, 2);
  emxInit_char_T(&a__12, 2);
  emxInit_real_T(&tmp, 2);
  emxInit_char_T(&remain, 2);
  emxInit_char_T(&curr_w, 2);
  emxInit_char_T(&b_remain, 2);
  emxInit_char_T(&curr_alphabet, 2);
  emxInit_int32_T(&match_out, 2);
  emxInit_int32_T(&matches, 2);
  do {
    exitg1 = 0;
    getfilestar(fid, &filestar, &b_bool);
    if (filestar == b_NULL) {
      ret = 0;
    } else {
      st = feof(filestar);
      ret = ((int)st != 0);
    }

    if (ret == 0) {
      /* 'mapTF:954' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF:955' if length(line) < 5 || ~strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF:958' i = i + 1; */
          b_i++;

          /* 'mapTF:959' line = fgetl(fid); */
          fgetl(fid, line);
          line_data = line->data;

          /* 'mapTF:961' [~, remain] = strtok(line); */
          b_strtok(line, a__9, remain);

          /* 'mapTF:962' [curr_w, remain] = strtok(remain); */
          b_strtok(remain, curr_w, b_remain);

          /* 'mapTF:963' [~, remain] = strtok(remain); */
          b_strtok(b_remain, a__9, remain);

          /* 'mapTF:964' [curr_alphabet, remain] = strtok(remain); */
          b_strtok(remain, curr_alphabet, b_remain);

          /* 'mapTF:965' [~, remain] = strtok(remain); */
          b_strtok(b_remain, a__9, remain);

          /* 'mapTF:966' [curr_length, ~] = strtok(remain); */
          b_strtok(remain, a__9, a__12);

          /* 'mapTF:967' w(i) = real(str2double(curr_w)); */
          dc = str2double(curr_w);
          w_data[(int)b_i - 1] = dc.re;

          /* 'mapTF:968' curr_alphabet = real(str2double(curr_alphabet)); */
          dc = str2double(curr_alphabet);

          /* 'mapTF:969' curr_length = real(str2double(curr_length)); */
          dc1 = str2double(a__9);

          /* 'mapTF:970' tmp = zeros(curr_length, curr_alphabet); */
          i = tmp->size[0] * tmp->size[1];
          tmp->size[0] = (int)dc1.re;
          tmp->size[1] = (int)dc.re;
          emxEnsureCapacity_real_T(tmp, i);
          tmp_data = tmp->data;
          ret = (int)dc1.re * (int)dc.re;
          for (i = 0; i < ret; i++) {
            tmp_data[i] = 0.0;
          }

          /* 'mapTF:971' line_count = 1; */
          line_count = 1U;

          /* 'mapTF:972' while ~isempty(line) */
          while (line->size[1] != 0) {
            /* 'mapTF:973' if strfind(line, 'weight') */
            i = matches->size[0] * matches->size[1];
            matches->size[0] = 1;
            matches->size[1] = line->size[1];
            emxEnsureCapacity_int32_T(matches, i);
            matches_data = matches->data;
            ret = 0;
            i = line->size[1];
            for (c_i = 0; c_i <= i - 6; c_i++) {
              j = 1;
              while ((j <= 6) && (line_data[(c_i + j) - 1] == cv1[j - 1])) {
                j++;
              }

              if (j > 6) {
                matches_data[ret] = c_i + 1;
                ret++;
              }
            }

            i = match_out->size[0] * match_out->size[1];
            match_out->size[0] = 1;
            match_out->size[1] = ret;
            emxEnsureCapacity_int32_T(match_out, i);
            match_out_data = match_out->data;
            for (c_i = 0; c_i < ret; c_i++) {
              match_out_data[c_i] = matches_data[c_i];
            }

            for (i = 0; i < 2; i++) {
              out_size[i] = match_out->size[i];
            }

            if (out_size[1] != 0) {
              /* 'mapTF:974' line = fgetl(fid); */
              fgetl(fid, line);
            }

            /* 'mapTF:976' [v1, remain] = strtok(line); */
            b_strtok(line, curr_w, remain);

            /* 'mapTF:977' [v2, remain] = strtok(remain); */
            b_strtok(remain, curr_alphabet, b_remain);

            /* 'mapTF:978' [v3, v4] = strtok(remain); */
            b_strtok(b_remain, a__9, a__12);

            /* 'mapTF:979' vals = [real(str2double(v1)),  */
            /* 'mapTF:980'                 real(str2double(v2)),  */
            /* 'mapTF:981'                 real(str2double(v3)),  */
            /* 'mapTF:982'                 real(str2double(v4))]'; */
            dc = str2double(curr_w);
            dc1 = str2double(curr_alphabet);
            dc2 = str2double(a__9);
            dc3 = str2double(a__12);
            vals[0] = dc.re;
            vals[1] = dc1.re;
            vals[2] = dc2.re;
            vals[3] = dc3.re;

            /* 'mapTF:983' tmp(line_count,:) = vals; */
            ret = tmp->size[1];
            for (i = 0; i < ret; i++) {
              tmp_data[((int)line_count + tmp->size[0] * i) - 1] = vals[i];
            }

            /* 'mapTF:984' line = fgetl(fid); */
            fgetl(fid, line);
            line_data = line->data;

            /* 'mapTF:985' line_count = line_count + 1; */
            line_count++;
          }

          /* 'mapTF:987' mat{i} = tmp; */
          i = mat_data[(int)b_i - 1].f1->size[0] * mat_data[(int)b_i - 1]
            .f1->size[1];
          mat_data[(int)b_i - 1].f1->size[0] = tmp->size[0];
          mat_data[(int)b_i - 1].f1->size[1] = tmp->size[1];
          emxEnsureCapacity_real_T(mat_data[(int)b_i - 1].f1, i);
          ret = tmp->size[0] * tmp->size[1];
          for (i = 0; i < ret; i++) {
            mat_data[(int)b_i - 1].f1->data[i] = tmp_data[i];
          }
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_int32_T(&matches);
  emxFree_int32_T(&match_out);
  emxFree_char_T(&curr_alphabet);
  emxFree_char_T(&b_remain);
  emxFree_char_T(&curr_w);
  emxFree_char_T(&remain);
  emxFree_real_T(&tmp);
  emxFree_char_T(&a__12);
  emxFree_char_T(&a__9);
  emxFree_char_T(&line);

  /* 'mapTF:989' fclose(fid); */
  cfclose(fid);
}

static void i_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *mat,
  const emxArray_real_T *r1)
{
  emxArray_real_T *b_mat;
  const double *mat_data;
  const double *r;
  double *b_mat_data;
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
  r = r1->data;
  mat_data = mat->data;
  emxInit_real_T(&b_mat, 2);
  i = b_mat->size[0] * b_mat->size[1];
  if (r1->size[0] == 1) {
    b_mat->size[0] = mat->size[0];
  } else {
    b_mat->size[0] = r1->size[0];
  }

  if (r1->size[1] == 1) {
    b_mat->size[1] = mat->size[1];
  } else {
    b_mat->size[1] = r1->size[1];
  }

  emxEnsureCapacity_real_T(b_mat, i);
  b_mat_data = b_mat->data;
  stride_0_0 = (mat->size[0] != 1);
  stride_0_1 = (mat->size[1] != 1);
  stride_1_0 = (r1->size[0] != 1);
  stride_1_1 = (r1->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (r1->size[1] == 1) {
    loop_ub = mat->size[1];
  } else {
    loop_ub = r1->size[1];
  }

  for (i = 0; i < loop_ub; i++) {
    if (r1->size[0] == 1) {
      b_loop_ub = mat->size[0];
    } else {
      b_loop_ub = r1->size[0];
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_mat_data[i1 + b_mat->size[0] * i] = mat_data[i1 * stride_0_0 + mat->
        size[0] * aux_0_1] * r[i1 * stride_1_0 + r1->size[0] * aux_1_1] /
        0.69314718055994529;
    }

    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  sum(b_mat, vec);
  emxFree_real_T(&b_mat);
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

  /* 'mapTF:389' l = length(s); */
  /* 'mapTF:390' en = zeros(1,l); */
  i = en->size[0] * en->size[1];
  en->size[0] = 1;
  en->size[1] = s->size[1];
  emxEnsureCapacity_real_T(en, i);
  en_data = en->data;

  /* 'mapTF:391' for i = 1:l */
  i = s->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:392' if strcmp(s(i),'A') || strcmp(s(i), 'a') */
    c = s_data[b_i];
    if ((!(c != 'A')) || (!(c != 'a'))) {
      /* 'mapTF:393' en(i) = 0; */
      en_data[b_i] = 0.0;
    } else if ((!(c != 'C')) || (!(c != 'c'))) {
      /* 'mapTF:394' elseif strcmp(s(i),'C') || strcmp(s(i),'c') */
      /* 'mapTF:395' en(i) = 1; */
      en_data[b_i] = 1.0;
    } else if ((!(c != 'G')) || (!(c != 'g'))) {
      /* 'mapTF:396' elseif strcmp(s(i),'G') || strcmp(s(i),'g') */
      /* 'mapTF:397' en(i) = 2; */
      en_data[b_i] = 2.0;
    } else {
      /* 'mapTF:398' else */
      /* 'mapTF:399' en(i) = 3; */
      en_data[b_i] = 3.0;
    }
  }
}

static void minus(emxArray_real_T *f, const emxArray_real_T *minnorm)
{
  emxArray_real_T *b_f;
  const double *minnorm_data;
  double *b_f_data;
  double *f_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  minnorm_data = minnorm->data;
  f_data = f->data;
  emxInit_real_T(&b_f, 1);
  i = b_f->size[0];
  if (minnorm->size[0] == 1) {
    b_f->size[0] = f->size[0];
  } else {
    b_f->size[0] = minnorm->size[0];
  }

  emxEnsureCapacity_real_T(b_f, i);
  b_f_data = b_f->data;
  stride_0_0 = (f->size[0] != 1);
  stride_1_0 = (minnorm->size[0] != 1);
  if (minnorm->size[0] == 1) {
    loop_ub = f->size[0];
  } else {
    loop_ub = minnorm->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_f_data[i] = f_data[i * stride_0_0] - minnorm_data[i * stride_1_0];
  }

  i = f->size[0];
  f->size[0] = b_f->size[0];
  emxEnsureCapacity_real_T(f, i);
  f_data = f->data;
  loop_ub = b_f->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = b_f_data[i];
  }

  emxFree_real_T(&b_f);
}

static void plus(emxArray_real_T *pwm_prob, const emxArray_real_T *b)
{
  emxArray_real_T *b_pwm_prob;
  const double *b_data;
  double *b_pwm_prob_data;
  double *pwm_prob_data;
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
  b_data = b->data;
  pwm_prob_data = pwm_prob->data;
  emxInit_real_T(&b_pwm_prob, 2);
  i = b_pwm_prob->size[0] * b_pwm_prob->size[1];
  if (b->size[0] == 1) {
    b_pwm_prob->size[0] = pwm_prob->size[0];
  } else {
    b_pwm_prob->size[0] = b->size[0];
  }

  if (b->size[1] == 1) {
    b_pwm_prob->size[1] = pwm_prob->size[1];
  } else {
    b_pwm_prob->size[1] = b->size[1];
  }

  emxEnsureCapacity_real_T(b_pwm_prob, i);
  b_pwm_prob_data = b_pwm_prob->data;
  stride_0_0 = (pwm_prob->size[0] != 1);
  stride_0_1 = (pwm_prob->size[1] != 1);
  stride_1_0 = (b->size[0] != 1);
  stride_1_1 = (b->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (b->size[1] == 1) {
    loop_ub = pwm_prob->size[1];
  } else {
    loop_ub = b->size[1];
  }

  for (i = 0; i < loop_ub; i++) {
    if (b->size[0] == 1) {
      b_loop_ub = pwm_prob->size[0];
    } else {
      b_loop_ub = b->size[0];
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_pwm_prob_data[i1 + b_pwm_prob->size[0] * i] = pwm_prob_data[i1 *
        stride_0_0 + pwm_prob->size[0] * aux_0_1] + b_data[i1 * stride_1_0 +
        b->size[0] * aux_1_1];
    }

    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  i = pwm_prob->size[0] * pwm_prob->size[1];
  pwm_prob->size[0] = b_pwm_prob->size[0];
  pwm_prob->size[1] = b_pwm_prob->size[1];
  emxEnsureCapacity_real_T(pwm_prob, i);
  pwm_prob_data = pwm_prob->data;
  loop_ub = b_pwm_prob->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = b_pwm_prob->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      pwm_prob_data[i1 + pwm_prob->size[0] * i] = b_pwm_prob_data[i1 +
        b_pwm_prob->size[0] * i];
    }
  }

  emxFree_real_T(&b_pwm_prob);
}

/*
 * function [ind, M] = ppmsim(mot,lenvec)
 */
static void ppmsim(emxArray_cell_wrap_6 *mot, const emxArray_real_T *lenvec,
                   double *ind, double *M)
{
  cell_wrap_6 *mot_data;
  emxArray_real_T *A;
  emxArray_real_T *diag_a;
  emxArray_real_T *diag_b;
  emxArray_real_T *diag_overall;
  emxArray_real_T *mat;
  emxArray_real_T *rmat;
  double MM;
  double *diag_overall_data;
  double *mat_data;
  double *rmat_data;
  int b_i;
  int i;
  int i1;
  int i2;
  int loop_ub;
  mot_data = mot->data;

  /* 'mapTF:907' n = length(lenvec)-1; */
  /* 'mapTF:908' simmat = ones(n-1,1); */
  /* 'mapTF:909' for i = 1:n+1 */
  i = lenvec->size[0];
  emxInit_real_T(&diag_a, 2);
  emxInit_real_T(&A, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:910' mot{i} = mot{i}-1/4; */
    loop_ub = mot_data[b_i].f1->size[0] * mot_data[b_i].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      mot_data[b_i].f1->data[i1] -= 0.25;
    }

    /* 'mapTF:911' mot{i} = mot{i}/sqrt(sum(sum(mot{i}.^2))); */
    i1 = A->size[0] * A->size[1];
    A->size[0] = mot_data[b_i].f1->size[0];
    A->size[1] = mot_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(A, i1);
    mat_data = A->data;
    loop_ub = mot_data[b_i].f1->size[0] * mot_data[b_i].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      MM = mot_data[b_i].f1->data[i1];
      mat_data[i1] = MM * MM;
    }

    b_sum(A, diag_a);
    MM = sqrt(c_sum(diag_a));
    loop_ub = mot_data[b_i].f1->size[0] * mot_data[b_i].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      mot_data[b_i].f1->data[i1] /= MM;
    }
  }

  /* 'mapTF:913' M = 0; */
  *M = 0.0;

  /* 'mapTF:914' ind = 1; */
  *ind = 1.0;

  /* 'mapTF:915' for j = 2:n+1 */
  i = lenvec->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&rmat, 2);
  emxInit_real_T(&diag_b, 2);
  emxInit_real_T(&diag_overall, 1);
  for (b_i = 0; b_i <= i - 2; b_i++) {
    /* 'mapTF:916' mat = mot{1}*mot{j}'; */
    if ((mot_data[0].f1->size[0] == 0) || (mot_data[0].f1->size[1] == 0) ||
        (mot_data[b_i + 1].f1->size[0] == 0) || (mot_data[b_i + 1].f1->size[1] ==
         0)) {
      i1 = mat->size[0] * mat->size[1];
      mat->size[0] = mot_data[0].f1->size[0];
      mat->size[1] = mot_data[b_i + 1].f1->size[0];
      emxEnsureCapacity_real_T(mat, i1);
      mat_data = mat->data;
      loop_ub = mot_data[0].f1->size[0] * mot_data[b_i + 1].f1->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        mat_data[i1] = 0.0;
      }
    } else {
      i1 = mat->size[0] * mat->size[1];
      mat->size[0] = mot_data[0].f1->size[0];
      mat->size[1] = mot_data[b_i + 1].f1->size[0];
      emxEnsureCapacity_real_T(mat, i1);
      mat_data = mat->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)mot_data[0].
                  f1->size[0], (blasint)mot_data[b_i + 1].f1->size[0], (blasint)
                  mot_data[0].f1->size[1], 1.0, &mot_data[0].f1->data[0],
                  (blasint)mot_data[0].f1->size[0], &mot_data[b_i + 1].f1->data
                  [0], (blasint)mot_data[b_i + 1].f1->size[0], 0.0, &mat_data[0],
                  (blasint)mot_data[0].f1->size[0]);
    }

    /* 'mapTF:917' rmat = rot90(mot{1},2)*mot{j}'; */
    rot90(mot_data[0].f1, A);
    mat_data = A->data;
    if ((A->size[0] == 0) || (A->size[1] == 0) || (mot_data[b_i + 1].f1->size[0]
         == 0) || (mot_data[b_i + 1].f1->size[1] == 0)) {
      i1 = rmat->size[0] * rmat->size[1];
      rmat->size[0] = A->size[0];
      rmat->size[1] = mot_data[b_i + 1].f1->size[0];
      emxEnsureCapacity_real_T(rmat, i1);
      rmat_data = rmat->data;
      loop_ub = A->size[0] * mot_data[b_i + 1].f1->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        rmat_data[i1] = 0.0;
      }
    } else {
      i1 = rmat->size[0] * rmat->size[1];
      rmat->size[0] = A->size[0];
      rmat->size[1] = mot_data[b_i + 1].f1->size[0];
      emxEnsureCapacity_real_T(rmat, i1);
      rmat_data = rmat->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)A->size[0],
                  (blasint)mot_data[b_i + 1].f1->size[0], (blasint)A->size[1],
                  1.0, &mat_data[0], (blasint)A->size[0], &mot_data[b_i + 1].
                  f1->data[0], (blasint)mot_data[b_i + 1].f1->size[0], 0.0,
                  &rmat_data[0], (blasint)A->size[0]);
    }

    /*  MM = max([sum(spdiags(mat)) sum(spdiags(rmat))]); */
    /* 'mapTF:919' diag_a = sum(spdiags(mat)); */
    spdiags(mat, A);
    b_sum(A, diag_a);
    mat_data = diag_a->data;

    /* 'mapTF:920' diag_b = sum(spdiags(rmat)); */
    spdiags(rmat, A);
    b_sum(A, diag_b);
    rmat_data = diag_b->data;

    /* 'mapTF:921' diag_overall = zeros(length(diag_a)+length(diag_b),1); */
    i1 = diag_overall->size[0];
    diag_overall->size[0] = (int)((unsigned int)diag_a->size[1] + diag_b->size[1]);
    emxEnsureCapacity_real_T(diag_overall, i1);
    diag_overall_data = diag_overall->data;
    loop_ub = (int)((unsigned int)diag_a->size[1] + diag_b->size[1]);
    for (i1 = 0; i1 < loop_ub; i1++) {
      diag_overall_data[i1] = 0.0;
    }

    /* 'mapTF:922' diag_overall(1:length(diag_a)) = diag_a; */
    loop_ub = diag_a->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      diag_overall_data[i1] = mat_data[i1];
    }

    /* 'mapTF:923' diag_overall(length(diag_a)+1:end) = diag_b; */
    if (diag_a->size[1] + 1U > (unsigned int)diag_overall->size[0]) {
      i1 = 0;
    } else {
      i1 = diag_a->size[1];
    }

    loop_ub = diag_b->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      diag_overall_data[i1 + i2] = rmat_data[i2];
    }

    /* 'mapTF:924' MM = max(diag_overall); */
    MM = maximum(diag_overall);

    /* 'mapTF:925' if MM > M */
    if (MM > *M) {
      /* 'mapTF:926' M = MM; */
      *M = MM;

      /* 'mapTF:927' ind = j-1; */
      *ind = ((double)b_i + 2.0) - 1.0;
    }
  }

  emxFree_real_T(&A);
  emxFree_real_T(&diag_overall);
  emxFree_real_T(&diag_b);
  emxFree_real_T(&diag_a);
  emxFree_real_T(&rmat);
  emxFree_real_T(&mat);
}

/*
 * function process_motifs(dfn, lfn, memefn, ofn, PWMcorrcut, Negative)
 */
static void process_motifs(const emxArray_char_T *dfn, const emxArray_char_T
  *lfn, const emxArray_char_T *memefn, const emxArray_char_T *ofn, double
  PWMcorrcut, double Negative)
{
  static const char b_cv[5] = { 'M', 'O', 'T', 'I', 'F' };

  FILE* b_NULL;
  FILE* c_NULL;
  FILE* d_NULL;
  FILE* filestar;
  cell_wrap_4 *names_data;
  cell_wrap_6 r1;
  cell_wrap_6 *PWM2_data;
  cell_wrap_6 *P_data;
  cell_wrap_6 *cur_PWM_data;
  cell_wrap_6 *cur_PWM_tmp_data;
  emxArray_boolean_T *b_clus;
  emxArray_cell_wrap_4 *names;
  emxArray_cell_wrap_6 *P;
  emxArray_cell_wrap_6 *PWM2;
  emxArray_cell_wrap_6 *b_cur_PWM_tmp;
  emxArray_cell_wrap_6 *b_new_PWM2;
  emxArray_cell_wrap_6 *c_cur_PWM_tmp;
  emxArray_cell_wrap_6 *cur_PWM;
  emxArray_cell_wrap_6 *cur_PWM_tmp;
  emxArray_cell_wrap_6 *d_cur_PWM_tmp;
  emxArray_cell_wrap_6 *new_PWM2;
  emxArray_char_T *b_fileid;
  emxArray_char_T *b_remain;
  emxArray_char_T *c_fileid;
  emxArray_char_T *cur_line;
  emxArray_char_T *d_fileid;
  emxArray_char_T *e_fileid;
  emxArray_char_T *f_fileid;
  emxArray_char_T *p1;
  emxArray_char_T *p2;
  emxArray_char_T *p4;
  emxArray_char_T *remain;
  emxArray_int32_T *match_out;
  emxArray_int32_T *matches;
  emxArray_int32_T *r;
  emxArray_real_T *I_2;
  emxArray_real_T *LEN_2;
  emxArray_real_T *clus;
  emxArray_real_T *f;
  emxArray_real_T *lasso_weight;
  emxArray_real_T *uid;
  emxArray_real_T *vec;
  emxArray_real_T *vec2;
  emxArray_real_T *w;
  emxArray_real_T *y;
  emxArray_real_T *zscore;
  creal_T dc;
  double validatedHoleFilling[4];
  double LEN_2_new_tmp[2];
  double a;
  double curr_pos;
  double idx;
  double x;
  double *I_2_data;
  double *clus_data;
  double *f_data;
  double *lasso_weight_data;
  double *uid_data;
  double *vec2_data;
  double *vec_data;
  double *w_data;
  double *zscore_data;
  int x_size[2];
  int N;
  unsigned int b_cur_idx;
  int b_i;
  int cur_idx;
  int exitg1;
  int i;
  int i1;
  int loop_ub;
  int match_idx;
  int u0;
  int unnamed_idx_0;
  int *match_out_data;
  int *matches_data;
  const char *ofn_data;
  signed char fileid;
  char *cur_line_data;
  bool exitg2;
  bool hal;
  bool *b_clus_data;
  ofn_data = ofn->data;
  emxInit_real_T(&f, 1);
  emxInit_cell_wrap_6(&cur_PWM_tmp, 1);

  /* dfn: file name for denovo motifs */
  /* lfn: file name for lasso motifs */
  /* memefn: file name for the meme input for gkmPWMlasso */
  /* ofn: output filename */
  /* 'mapTF:656' a = 1; */
  a = 1.0;

  /* 'mapTF:657' LEN = zeros(1,1); */
  /* 'mapTF:658' shift = zeros(1,1); */
  /* 'mapTF:659' [p,w] = getdenovomotif(dfn); */
  getdenovomotif(dfn, cur_PWM_tmp, f);
  f_data = f->data;
  cur_PWM_tmp_data = cur_PWM_tmp->data;

  /* 'mapTF:660' N = numel(w); */
  N = f->size[0];

  /*  fid = fopen(lfn,'r'); */
  /*  X = textscan(fid,'%f\t%f\t%s\t%f\t%f\t%f\n','delimiter', '\t', 'headerlines', 4); */
  /*  fclose(fid); */
  /* 'mapTF:665' fid = fopen(lfn, 'r'); */
  fileid = cfopen(lfn, "rb");

  /* 'mapTF:666' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF:667' fprintf("ERROR: gkmPWMlasso output file cannot be opened.\n") */
    printf("ERROR: gkmPWMlasso output file cannot be opened.\n");
    fflush(stdout);
    exit(1);
  }

  emxInit_char_T(&b_fileid, 2);
  emxInit_char_T(&c_fileid, 2);
  emxInit_char_T(&d_fileid, 2);
  emxInit_char_T(&e_fileid, 2);

  /* 'mapTF:670' fgetl(fid); */
  b_fgets(fileid, b_fileid);

  /* 'mapTF:671' fgetl(fid); */
  b_fgets(fileid, c_fileid);

  /* 'mapTF:672' fgetl(fid); */
  b_fgets(fileid, d_fileid);

  /* 'mapTF:673' fgetl(fid); */
  b_fgets(fileid, e_fileid);

  /* 'mapTF:674' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF:675' idx=0; */
  idx = 0.0;

  /* 'mapTF:676' while ~feof(fid) */
  emxFree_char_T(&e_fileid);
  emxFree_char_T(&d_fileid);
  emxFree_char_T(&c_fileid);
  emxFree_char_T(&b_fileid);
  emxInit_char_T(&f_fileid, 2);
  do {
    exitg1 = 0;
    x = b_feof(fileid);
    if (!(x != 0.0)) {
      /* 'mapTF:677' idx=idx+1; */
      idx++;

      /* 'mapTF:678' fgetl(fid); */
      b_fgets(fileid, f_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&f_fileid);
  emxInit_real_T(&clus, 1);

  /* 'mapTF:680' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF:681' clus = zeros(idx, 1); */
  match_idx = (int)idx;
  i = clus->size[0];
  clus->size[0] = (int)idx;
  emxEnsureCapacity_real_T(clus, i);
  clus_data = clus->data;
  for (i = 0; i < match_idx; i++) {
    clus_data[i] = 0.0;
  }

  emxInit_real_T(&uid, 1);

  /* 'mapTF:682' uid = zeros(idx, 1); */
  i = uid->size[0];
  uid->size[0] = (int)idx;
  emxEnsureCapacity_real_T(uid, i);
  uid_data = uid->data;
  for (i = 0; i < match_idx; i++) {
    uid_data[i] = 0.0;
  }

  emxInit_real_T(&lasso_weight, 1);

  /* 'mapTF:683' lasso_weight = zeros(idx, 1); */
  i = lasso_weight->size[0];
  lasso_weight->size[0] = (int)idx;
  emxEnsureCapacity_real_T(lasso_weight, i);
  lasso_weight_data = lasso_weight->data;
  for (i = 0; i < match_idx; i++) {
    lasso_weight_data[i] = 0.0;
  }

  emxInit_real_T(&zscore, 1);

  /* 'mapTF:684' zscore = zeros(idx, 1); */
  i = zscore->size[0];
  zscore->size[0] = (int)idx;
  emxEnsureCapacity_real_T(zscore, i);
  zscore_data = zscore->data;
  for (i = 0; i < match_idx; i++) {
    zscore_data[i] = 0.0;
  }

  /* 'mapTF:686' for cur_idx=1:idx */
  cur_idx = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&p1, 2);
  emxInit_char_T(&remain, 2);
  emxInit_char_T(&p2, 2);
  emxInit_char_T(&b_remain, 2);
  emxInit_char_T(&p4, 2);
  exitg2 = false;
  while ((!exitg2) && (cur_idx <= (int)idx - 1)) {
    /* 'mapTF:687' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF:688' if cur_line == -1 */
    for (i = 0; i < 2; i++) {
      x_size[i] = cur_line->size[i];
    }

    hal = (x_size[1] != 0);
    if (hal) {
      hal = (0 > x_size[1] - 1);
    }

    if (hal) {
      exitg2 = true;
    } else {
      /* 'mapTF:691' [p1, remain] = strtok(cur_line, char(9)); */
      c_strtok(cur_line, p1, remain);

      /* 'mapTF:692' [p2, remain] = strtok(remain, char(9)); */
      c_strtok(remain, p2, b_remain);

      /* 'mapTF:693' [p3, remain] = strtok(remain, char(9)); */
      c_strtok(b_remain, cur_line, remain);

      /* 'mapTF:694' [p4, remain] = strtok(remain, char(9)); */
      c_strtok(remain, p4, b_remain);

      /* 'mapTF:695' [p5, p6] = strtok(remain, char(9)); */
      c_strtok(b_remain, cur_line, remain);

      /* 'mapTF:696' clus(cur_idx, 1) = real(str2double(p1)); */
      dc = str2double(p1);
      clus_data[cur_idx] = dc.re;

      /* 'mapTF:697' uid(cur_idx, 1) = real(str2double(p2)); */
      dc = str2double(p2);
      uid_data[cur_idx] = dc.re;

      /* 'mapTF:698' lasso_weight(cur_idx, 1) = real(str2double(p4)); */
      dc = str2double(p4);
      lasso_weight_data[cur_idx] = dc.re;

      /* 'mapTF:699' zscore(cur_idx, 1) = real(str2double(p5)); */
      dc = str2double(cur_line);
      zscore_data[cur_idx] = dc.re;
      cur_idx++;
    }
  }

  emxFree_char_T(&p4);
  emxFree_char_T(&b_remain);
  emxFree_char_T(&p2);
  emxFree_char_T(&remain);
  emxFree_char_T(&p1);
  emxInit_real_T(&w, 1);
  emxInit_real_T(&vec, 1);
  emxInit_real_T(&vec2, 1);

  /* 'mapTF:701' fclose(fid); */
  cfclose(fileid);

  /*  n = X{1}(end); */
  /* 'mapTF:705' n = clus(end); */
  /* 'mapTF:706' vec = zeros(n,1); */
  i = vec->size[0];
  vec->size[0] = (int)clus_data[clus->size[0] - 1];
  emxEnsureCapacity_real_T(vec, i);
  vec_data = vec->data;

  /* 'mapTF:707' vec2 = zeros(n,1); */
  i = vec2->size[0];
  vec2->size[0] = (int)clus_data[clus->size[0] - 1];
  emxEnsureCapacity_real_T(vec2, i);
  vec2_data = vec2->data;

  /* 'mapTF:708' w = [w; zeros(n,1)]; */
  i = w->size[0];
  w->size[0] = f->size[0] + (int)clus_data[clus->size[0] - 1];
  emxEnsureCapacity_real_T(w, i);
  w_data = w->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    w_data[i] = f_data[i];
  }

  loop_ub = (int)clus_data[clus->size[0] - 1];
  for (i = 0; i < loop_ub; i++) {
    w_data[i + f->size[0]] = 0.0;
  }

  /* 'mapTF:709' for i = 1:n */
  i = (int)clus_data[clus->size[0] - 1];
  emxInit_int32_T(&r, 1);
  emxInit_boolean_T(&b_clus, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /*      f = find(X{1}==i); */
    /*      vec(i) = X{2}(f(1)); */
    /*      vec2(i) = X{5}(f(1)); */
    /*      w(i+N) = X{4}(f(1)); */
    /* 'mapTF:714' f = find(clus == i); */
    loop_ub = clus->size[0];
    i1 = b_clus->size[0];
    b_clus->size[0] = clus->size[0];
    emxEnsureCapacity_boolean_T(b_clus, i1);
    b_clus_data = b_clus->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_clus_data[i1] = (clus_data[i1] == (double)b_i + 1.0);
    }

    eml_find(b_clus, r);
    matches_data = r->data;
    i1 = f->size[0];
    f->size[0] = r->size[0];
    emxEnsureCapacity_real_T(f, i1);
    f_data = f->data;
    loop_ub = r->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[i1] = matches_data[i1];
    }

    /* 'mapTF:715' vec(i) = uid(f(1)); */
    vec_data[b_i] = uid_data[(int)f_data[0] - 1];

    /* 'mapTF:716' vec2(i) = zscore(f(1)); */
    vec2_data[b_i] = zscore_data[(int)f_data[0] - 1];

    /* 'mapTF:717' w(i+N) = lasso_weight(f(1)); */
    w_data[(int)((unsigned int)b_i + N)] = lasso_weight_data[(int)f_data[0] - 1];
  }

  emxFree_boolean_T(&b_clus);
  emxFree_int32_T(&r);
  emxFree_real_T(&zscore);

  /* 'mapTF:719' if Negative */
  if (Negative != 0.0) {
    /* 'mapTF:720' w = -1 * w; */
    loop_ub = w->size[0];
    for (i = 0; i < loop_ub; i++) {
      w_data[i] = -w_data[i];
    }

    /* 'mapTF:721' vec2 = -1 * vec2; */
    loop_ub = vec2->size[0];
    for (i = 0; i < loop_ub; i++) {
      vec2_data[i] = -vec2_data[i];
    }
  }

  emxInit_cell_wrap_6(&P, 1);
  emxInit_cell_wrap_6(&cur_PWM, 1);

  /* 'mapTF:723' P = getmotif(memefn,vec); */
  getmotif(memefn, vec, P);
  P_data = P->data;

  /*  [pp, info, len] = trim_pwm([p;P],0.25); */
  /* 'mapTF:726' denovo_len = length(p); */
  /* 'mapTF:727' database_len = length(P); */
  /* 'mapTF:728' cur_PWM = cell(denovo_len+database_len,1); */
  unnamed_idx_0 = (int)((unsigned int)cur_PWM_tmp->size[0] + P->size[0]);
  i = cur_PWM->size[0];
  cur_PWM->size[0] = (int)((unsigned int)cur_PWM_tmp->size[0] + P->size[0]);
  emxEnsureCapacity_cell_wrap_6(cur_PWM, i);
  cur_PWM_data = cur_PWM->data;
  emxFree_real_T(&vec);
  for (i = 0; i < unnamed_idx_0; i++) {
    cur_PWM_data[i].f1->size[0] = 0;
    cur_PWM_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:729' cur_PWM = coder.nullcopy(cur_PWM); */
  /* 'mapTF:730' for cur_idx=1:denovo_len */
  i = cur_PWM_tmp->size[0];
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'mapTF:731' cur_PWM{cur_idx} = p{cur_idx}; */
    i1 = cur_PWM_data[cur_idx].f1->size[0] * cur_PWM_data[cur_idx].f1->size[1];
    cur_PWM_data[cur_idx].f1->size[0] = cur_PWM_tmp_data[cur_idx].f1->size[0];
    cur_PWM_data[cur_idx].f1->size[1] = cur_PWM_tmp_data[cur_idx].f1->size[1];
    emxEnsureCapacity_real_T(cur_PWM_data[cur_idx].f1, i1);
    loop_ub = cur_PWM_tmp_data[cur_idx].f1->size[0] * cur_PWM_tmp_data[cur_idx].
      f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      cur_PWM_data[cur_idx].f1->data[i1] = cur_PWM_tmp_data[cur_idx].f1->data[i1];
    }
  }

  /* 'mapTF:733' for cur_idx=denovo_len+1:denovo_len+database_len */
  i = P->size[0];
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    b_cur_idx = ((unsigned int)cur_PWM_tmp->size[0] + cur_idx) + 1U;

    /* 'mapTF:734' cur_PWM{cur_idx} = P{cur_idx-denovo_len}; */
    i1 = cur_PWM_data[(int)b_cur_idx - 1].f1->size[0] * cur_PWM_data[(int)
      b_cur_idx - 1].f1->size[1];
    cur_PWM_data[(int)b_cur_idx - 1].f1->size[0] = P_data[(int)((double)
      b_cur_idx - (double)cur_PWM_tmp->size[0]) - 1].f1->size[0];
    cur_PWM_data[(int)b_cur_idx - 1].f1->size[1] = P_data[(int)((double)
      b_cur_idx - (double)cur_PWM_tmp->size[0]) - 1].f1->size[1];
    emxEnsureCapacity_real_T(cur_PWM_data[(int)b_cur_idx - 1].f1, i1);
    loop_ub = P_data[(int)((double)b_cur_idx - (double)cur_PWM_tmp->size[0]) - 1]
      .f1->size[0] * P_data[(int)((double)b_cur_idx - (double)cur_PWM_tmp->size
      [0]) - 1].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      cur_PWM_data[(int)b_cur_idx - 1].f1->data[i1] = P_data[(int)((double)
        b_cur_idx - (double)cur_PWM_tmp->size[0]) - 1].f1->data[i1];
    }
  }

  emxFree_cell_wrap_6(&P);
  emxInitStruct_cell_wrap_6(&r1);

  /* 'mapTF:736' [pp, info, len] = trim_pwm(cur_PWM,0.25); */
  trim_pwm(cur_PWM, lasso_weight, uid);
  uid_data = uid->data;
  lasso_weight_data = lasso_weight->data;
  cur_PWM_data = cur_PWM->data;

  /* 'mapTF:738' coder.varsize("PWM2"); */
  /* 'mapTF:739' PWM2 = {pp{1}}; */
  i = r1.f1->size[0] * r1.f1->size[1];
  r1.f1->size[0] = cur_PWM_data[0].f1->size[0];
  r1.f1->size[1] = cur_PWM_data[0].f1->size[1];
  emxEnsureCapacity_real_T(r1.f1, i);
  loop_ub = cur_PWM_data[0].f1->size[0] * cur_PWM_data[0].f1->size[1];
  for (i = 0; i < loop_ub; i++) {
    r1.f1->data[i] = cur_PWM_data[0].f1->data[i];
  }

  emxInit_cell_wrap_6(&PWM2, 2);
  emxInit_real_T(&LEN_2, 2);
  emxInit_real_T(&I_2, 2);
  i = PWM2->size[0] * PWM2->size[1];
  PWM2->size[0] = 1;
  PWM2->size[1] = 1;
  emxEnsureCapacity_cell_wrap_6(PWM2, i);
  PWM2_data = PWM2->data;
  emxCopyStruct_cell_wrap_6(&PWM2_data[0], &r1);

  /* 'mapTF:740' coder.varsize("LEN_2"); */
  /* 'mapTF:741' LEN_2 = zeros(1,1); */
  i = LEN_2->size[0] * LEN_2->size[1];
  LEN_2->size[0] = 1;
  LEN_2->size[1] = 1;
  emxEnsureCapacity_real_T(LEN_2, i);
  vec_data = LEN_2->data;
  vec_data[0] = 0.0;

  /* 'mapTF:742' coder.varsize("I_2"); */
  /* 'mapTF:743' I_2 = zeros(1,1); */
  i = I_2->size[0] * I_2->size[1];
  I_2->size[0] = 1;
  I_2->size[1] = 1;
  emxEnsureCapacity_real_T(I_2, i);
  I_2_data = I_2->data;
  I_2_data[0] = 0.0;

  /* 'mapTF:744' hal = true; */
  hal = true;

  /* 'mapTF:745' for ii = 1:length(w) */
  i = w->size[0];
  emxFreeStruct_cell_wrap_6(&r1);
  emxInit_cell_wrap_6(&b_cur_PWM_tmp, 1);
  emxInit_cell_wrap_6(&new_PWM2, 1);
  emxInit_cell_wrap_6(&c_cur_PWM_tmp, 1);
  emxInit_cell_wrap_6(&b_new_PWM2, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:746' if ii > N */
    if (b_i + 1 > N) {
      /* 'mapTF:747' hal = true; */
      hal = true;

      /*  [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]); */
      /* 'mapTF:749' PWM2_len = length(PWM2); */
      u0 = PWM2->size[0];
      if (u0 < 1) {
        u0 = 1;
      }

      /* 'mapTF:750' cur_PWM_tmp = cell(PWM2_len+1,1); */
      unnamed_idx_0 = u0 + 1;
      i1 = c_cur_PWM_tmp->size[0];
      c_cur_PWM_tmp->size[0] = u0 + 1;
      emxEnsureCapacity_cell_wrap_6(c_cur_PWM_tmp, i1);
      cur_PWM_tmp_data = c_cur_PWM_tmp->data;
      for (i1 = 0; i1 < unnamed_idx_0; i1++) {
        cur_PWM_tmp_data[i1].f1->size[0] = 0;
        cur_PWM_tmp_data[i1].f1->size[1] = 0;
      }

      /* 'mapTF:751' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
      i1 = cur_PWM_tmp->size[0];
      cur_PWM_tmp->size[0] = c_cur_PWM_tmp->size[0];
      emxEnsureCapacity_cell_wrap_6(cur_PWM_tmp, i1);
      cur_PWM_tmp_data = cur_PWM_tmp->data;

      /* 'mapTF:752' cur_PWM_tmp{1} = pp{ii}; */
      i1 = cur_PWM_tmp_data[0].f1->size[0] * cur_PWM_tmp_data[0].f1->size[1];
      cur_PWM_tmp_data[0].f1->size[0] = cur_PWM_data[b_i].f1->size[0];
      cur_PWM_tmp_data[0].f1->size[1] = cur_PWM_data[b_i].f1->size[1];
      emxEnsureCapacity_real_T(cur_PWM_tmp_data[0].f1, i1);
      loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cur_PWM_tmp_data[0].f1->data[i1] = cur_PWM_data[b_i].f1->data[i1];
      }

      /* 'mapTF:753' for cur_idx=1:PWM2_len */
      for (cur_idx = 0; cur_idx < u0; cur_idx++) {
        /* 'mapTF:754' cur_PWM_tmp{cur_idx+1} = PWM2{cur_idx}; */
        i1 = cur_PWM_tmp_data[cur_idx + 1].f1->size[0] *
          cur_PWM_tmp_data[cur_idx + 1].f1->size[1];
        cur_PWM_tmp_data[cur_idx + 1].f1->size[0] = PWM2_data[cur_idx].f1->size
          [0];
        cur_PWM_tmp_data[cur_idx + 1].f1->size[1] = PWM2_data[cur_idx].f1->size
          [1];
        emxEnsureCapacity_real_T(cur_PWM_tmp_data[cur_idx + 1].f1, i1);
        loop_ub = PWM2_data[cur_idx].f1->size[0] * PWM2_data[cur_idx].f1->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          cur_PWM_tmp_data[cur_idx + 1].f1->data[i1] = PWM2_data[cur_idx]
            .f1->data[i1];
        }
      }

      /* 'mapTF:756' cur_length = zeros(length(LEN_2)+1,1); */
      u0 = LEN_2->size[0];
      if (u0 < 1) {
        u0 = 1;
      }

      i1 = f->size[0];
      f->size[0] = u0 + 1;
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      for (i1 = 0; i1 <= u0; i1++) {
        f_data[i1] = 0.0;
      }

      /* 'mapTF:757' cur_length(1) = len(ii); */
      f_data[0] = uid_data[b_i];

      /* 'mapTF:758' cur_length(2:end) = LEN_2; */
      for (i1 = 0; i1 < u0; i1++) {
        f_data[i1 + 1] = vec_data[i1];
      }

      /* 'mapTF:759' [~,cor] = ppmsim(cur_PWM_tmp, cur_length); */
      ppmsim(cur_PWM_tmp, f, &curr_pos, &idx);

      /* 'mapTF:761' if cor > PWMcorrcut || vec2(ii-N) < 1.5 */
      if ((idx > PWMcorrcut) || (vec2_data[b_i - N] < 1.5)) {
        /* 'mapTF:762' hal = false; */
        hal = false;
      }
    } else if ((b_i + 1 > 1) && (a > 2.0)) {
      /* 'mapTF:764' elseif ii > 1 && a > 2 */
      /* 'mapTF:765' hal = true; */
      hal = true;

      /*  [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]); */
      /* 'mapTF:767' PWM2_len = length(PWM2); */
      u0 = PWM2->size[0];
      if (u0 < 1) {
        u0 = 1;
      }

      /* 'mapTF:768' cur_PWM_tmp = cell(PWM2_len+1,1); */
      unnamed_idx_0 = u0 + 1;
      i1 = b_cur_PWM_tmp->size[0];
      b_cur_PWM_tmp->size[0] = u0 + 1;
      emxEnsureCapacity_cell_wrap_6(b_cur_PWM_tmp, i1);
      P_data = b_cur_PWM_tmp->data;
      for (i1 = 0; i1 < unnamed_idx_0; i1++) {
        P_data[i1].f1->size[0] = 0;
        P_data[i1].f1->size[1] = 0;
      }

      /* 'mapTF:769' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
      i1 = cur_PWM_tmp->size[0];
      cur_PWM_tmp->size[0] = b_cur_PWM_tmp->size[0];
      emxEnsureCapacity_cell_wrap_6(cur_PWM_tmp, i1);
      cur_PWM_tmp_data = cur_PWM_tmp->data;

      /* 'mapTF:770' cur_PWM_tmp{1} = pp{ii}; */
      i1 = cur_PWM_tmp_data[0].f1->size[0] * cur_PWM_tmp_data[0].f1->size[1];
      cur_PWM_tmp_data[0].f1->size[0] = cur_PWM_data[b_i].f1->size[0];
      cur_PWM_tmp_data[0].f1->size[1] = cur_PWM_data[b_i].f1->size[1];
      emxEnsureCapacity_real_T(cur_PWM_tmp_data[0].f1, i1);
      loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cur_PWM_tmp_data[0].f1->data[i1] = cur_PWM_data[b_i].f1->data[i1];
      }

      /* 'mapTF:771' for cur_idx=1:PWM2_len */
      for (cur_idx = 0; cur_idx < u0; cur_idx++) {
        /* 'mapTF:772' cur_PWM_tmp{cur_idx+1} = PWM2{cur_idx}; */
        i1 = cur_PWM_tmp_data[cur_idx + 1].f1->size[0] *
          cur_PWM_tmp_data[cur_idx + 1].f1->size[1];
        cur_PWM_tmp_data[cur_idx + 1].f1->size[0] = PWM2_data[cur_idx].f1->size
          [0];
        cur_PWM_tmp_data[cur_idx + 1].f1->size[1] = PWM2_data[cur_idx].f1->size
          [1];
        emxEnsureCapacity_real_T(cur_PWM_tmp_data[cur_idx + 1].f1, i1);
        loop_ub = PWM2_data[cur_idx].f1->size[0] * PWM2_data[cur_idx].f1->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          cur_PWM_tmp_data[cur_idx + 1].f1->data[i1] = PWM2_data[cur_idx]
            .f1->data[i1];
        }
      }

      /* 'mapTF:774' cur_length = zeros(length(LEN_2)+1,1); */
      u0 = LEN_2->size[0];
      if (u0 < 1) {
        u0 = 1;
      }

      i1 = f->size[0];
      f->size[0] = u0 + 1;
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      for (i1 = 0; i1 <= u0; i1++) {
        f_data[i1] = 0.0;
      }

      /* 'mapTF:775' cur_length(1) = len(ii); */
      f_data[0] = uid_data[b_i];

      /* 'mapTF:776' cur_length(2:end) = LEN_2; */
      for (i1 = 0; i1 < u0; i1++) {
        f_data[i1 + 1] = vec_data[i1];
      }

      /* 'mapTF:777' [~,cor] = ppmsim(cur_PWM_tmp, cur_length); */
      ppmsim(cur_PWM_tmp, f, &curr_pos, &idx);

      /* 'mapTF:779' if cor > PWMcorrcut */
      if (idx > PWMcorrcut) {
        /* 'mapTF:780' hal = false; */
        hal = false;
      }
    }

    /* 'mapTF:783' if hal && w(ii) > 0 && len(ii) >= 6 */
    if (hal && (w_data[b_i] > 0.0) && (uid_data[b_i] >= 6.0)) {
      /* 'mapTF:784' if len(ii) > 10 */
      if (uid_data[b_i] > 10.0) {
        /* 'mapTF:785' if info(ii)/len(ii) > 0.7 */
        if (lasso_weight_data[b_i] / uid_data[b_i] > 0.7) {
          /*  PWM2{a} = pp{ii}; */
          /*  LEN_2(a) = len(ii); */
          /*  I_2(a) = info(ii); */
          /*  a = a+1; */
          /*  PWM2{a} = rot90(PWM2{a-1},2); */
          /*  LEN_2(a) = LEN_2(a-1); */
          /*  I_2(a) = I_2(a-1); */
          /*  a = a+1; */
          /* 'mapTF:795' new_len = a+1; */
          /* 'mapTF:796' new_PWM2 = cell(new_len,1); */
          match_idx = (int)(a + 1.0);
          i1 = b_new_PWM2->size[0];
          b_new_PWM2->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_cell_wrap_6(b_new_PWM2, i1);
          P_data = b_new_PWM2->data;
          for (i1 = 0; i1 < match_idx; i1++) {
            P_data[i1].f1->size[0] = 0;
            P_data[i1].f1->size[1] = 0;
          }

          /* 'mapTF:797' new_PWM2 = coder.nullcopy(new_PWM2); */
          i1 = cur_PWM_tmp->size[0];
          cur_PWM_tmp->size[0] = b_new_PWM2->size[0];
          emxEnsureCapacity_cell_wrap_6(cur_PWM_tmp, i1);
          cur_PWM_tmp_data = cur_PWM_tmp->data;

          /* 'mapTF:798' for idx = 1:new_len */
          for (unnamed_idx_0 = 0; unnamed_idx_0 < match_idx; unnamed_idx_0++) {
            /* 'mapTF:799' if idx <= a-1 */
            if ((double)unnamed_idx_0 + 1.0 <= a - 1.0) {
              /* 'mapTF:800' new_PWM2{idx} = PWM2{idx}; */
              i1 = cur_PWM_tmp_data[unnamed_idx_0].f1->size[0] *
                cur_PWM_tmp_data[unnamed_idx_0].f1->size[1];
              cur_PWM_tmp_data[unnamed_idx_0].f1->size[0] =
                PWM2_data[unnamed_idx_0].f1->size[0];
              cur_PWM_tmp_data[unnamed_idx_0].f1->size[1] =
                PWM2_data[unnamed_idx_0].f1->size[1];
              emxEnsureCapacity_real_T(cur_PWM_tmp_data[unnamed_idx_0].f1, i1);
              loop_ub = PWM2_data[unnamed_idx_0].f1->size[0] *
                PWM2_data[unnamed_idx_0].f1->size[1];
              for (i1 = 0; i1 < loop_ub; i1++) {
                cur_PWM_tmp_data[unnamed_idx_0].f1->data[i1] =
                  PWM2_data[unnamed_idx_0].f1->data[i1];
              }
            } else if ((double)unnamed_idx_0 + 1.0 == a) {
              /* 'mapTF:801' elseif idx == a */
              /* 'mapTF:802' new_PWM2{idx} = pp{ii}; */
              i1 = cur_PWM_tmp_data[unnamed_idx_0].f1->size[0] *
                cur_PWM_tmp_data[unnamed_idx_0].f1->size[1];
              cur_PWM_tmp_data[unnamed_idx_0].f1->size[0] = cur_PWM_data[b_i].
                f1->size[0];
              cur_PWM_tmp_data[unnamed_idx_0].f1->size[1] = cur_PWM_data[b_i].
                f1->size[1];
              emxEnsureCapacity_real_T(cur_PWM_tmp_data[unnamed_idx_0].f1, i1);
              loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i]
                .f1->size[1];
              for (i1 = 0; i1 < loop_ub; i1++) {
                cur_PWM_tmp_data[unnamed_idx_0].f1->data[i1] = cur_PWM_data[b_i]
                  .f1->data[i1];
              }
            } else if ((double)unnamed_idx_0 + 1.0 == a + 1.0) {
              /* 'mapTF:803' elseif idx == a+1 */
              /* 'mapTF:804' new_PWM2{a+1} = rot90(pp{ii},2); */
              rot90(cur_PWM_data[b_i].f1, cur_PWM_tmp_data[(int)(a + 1.0) - 1].
                    f1);
            }
          }

          /* 'mapTF:807' PWM2_len = length(new_PWM2); */
          /* 'mapTF:808' PWM2 = cell(PWM2_len, 1); */
          x_size[0] = cur_PWM_tmp->size[0];
          i1 = PWM2->size[0] * PWM2->size[1];
          PWM2->size[0] = cur_PWM_tmp->size[0];
          PWM2->size[1] = 1;
          emxEnsureCapacity_cell_wrap_6(PWM2, i1);
          PWM2_data = PWM2->data;
          for (i1 = 0; i1 < x_size[0]; i1++) {
            PWM2_data[i1].f1->size[0] = 0;
            PWM2_data[i1].f1->size[1] = 0;
          }

          /* 'mapTF:809' PWM2 = coder.nullcopy(PWM2); */
          /* 'mapTF:810' for cur_idx = 1:PWM2_len */
          i1 = cur_PWM_tmp->size[0];
          for (cur_idx = 0; cur_idx < i1; cur_idx++) {
            /* 'mapTF:811' PWM2{cur_idx} = new_PWM2{cur_idx}; */
            unnamed_idx_0 = PWM2_data[cur_idx].f1->size[0] * PWM2_data[cur_idx].
              f1->size[1];
            PWM2_data[cur_idx].f1->size[0] = cur_PWM_tmp_data[cur_idx].f1->size
              [0];
            PWM2_data[cur_idx].f1->size[1] = cur_PWM_tmp_data[cur_idx].f1->size
              [1];
            emxEnsureCapacity_real_T(PWM2_data[cur_idx].f1, unnamed_idx_0);
            loop_ub = cur_PWM_tmp_data[cur_idx].f1->size[0] *
              cur_PWM_tmp_data[cur_idx].f1->size[1];
            for (unnamed_idx_0 = 0; unnamed_idx_0 < loop_ub; unnamed_idx_0++) {
              PWM2_data[cur_idx].f1->data[unnamed_idx_0] =
                cur_PWM_tmp_data[cur_idx].f1->data[unnamed_idx_0];
            }
          }

          /* 'mapTF:813' LEN_2_new = zeros(a+1,1); */
          i1 = f->size[0];
          f->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_real_T(f, i1);
          f_data = f->data;
          for (i1 = 0; i1 < match_idx; i1++) {
            f_data[i1] = 0.0;
          }

          /* 'mapTF:814' LEN_2_new(1:a-1) = LEN_2(1:a-1); */
          if (1.0 > a - 1.0) {
            loop_ub = 0;
          } else {
            loop_ub = (int)(a - 1.0);
          }

          for (i1 = 0; i1 < loop_ub; i1++) {
            f_data[i1] = vec_data[i1];
          }

          /* 'mapTF:815' LEN_2_new(a:a+1) = [len(ii), len(ii)]; */
          for (i1 = 0; i1 < 2; i1++) {
            LEN_2_new_tmp[i1] = a + (double)i1;
          }

          f_data[(int)LEN_2_new_tmp[0] - 1] = uid_data[b_i];
          f_data[(int)LEN_2_new_tmp[1] - 1] = uid_data[b_i];

          /* 'mapTF:817' I_2_new = zeros(a+1,1); */
          i1 = clus->size[0];
          clus->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_real_T(clus, i1);
          clus_data = clus->data;
          for (i1 = 0; i1 < match_idx; i1++) {
            clus_data[i1] = 0.0;
          }

          /* 'mapTF:818' I_2_new(1:a-1) = I_2(1:a-1); */
          if (1.0 > a - 1.0) {
            loop_ub = 0;
          } else {
            loop_ub = (int)(a - 1.0);
          }

          for (i1 = 0; i1 < loop_ub; i1++) {
            clus_data[i1] = I_2_data[i1];
          }

          /* 'mapTF:819' I_2_new(a:a+1) = [info(ii), info(ii)]; */
          clus_data[(int)LEN_2_new_tmp[0] - 1] = lasso_weight_data[b_i];
          clus_data[(int)LEN_2_new_tmp[1] - 1] = lasso_weight_data[b_i];

          /* 'mapTF:821' LEN_2 = LEN_2_new; */
          i1 = LEN_2->size[0] * LEN_2->size[1];
          LEN_2->size[0] = f->size[0];
          LEN_2->size[1] = 1;
          emxEnsureCapacity_real_T(LEN_2, i1);
          vec_data = LEN_2->data;
          loop_ub = f->size[0];
          for (i1 = 0; i1 < loop_ub; i1++) {
            vec_data[i1] = f_data[i1];
          }

          /* 'mapTF:822' I_2 = I_2_new; */
          i1 = I_2->size[0] * I_2->size[1];
          I_2->size[0] = clus->size[0];
          I_2->size[1] = 1;
          emxEnsureCapacity_real_T(I_2, i1);
          I_2_data = I_2->data;
          loop_ub = clus->size[0];
          for (i1 = 0; i1 < loop_ub; i1++) {
            I_2_data[i1] = clus_data[i1];
          }

          /* 'mapTF:823' a = a+1; */
          a++;

          /* 'mapTF:824' a = a+1; */
          a++;
        }
      } else if ((lasso_weight_data[b_i] > 6.0) || (lasso_weight_data[b_i] /
                  uid_data[b_i] > 1.0)) {
        /* 'mapTF:826' elseif info(ii) > 6 || info(ii)/len(ii) > 1 */
        /*  PWM2{a} = pp{ii}; */
        /*  LEN_2(a) = len(ii); */
        /*  I_2(a) = info(ii); */
        /*  a = a+1; */
        /*  PWM2{a} = rot90(PWM2{a-1},2); */
        /*  LEN_2(a) = LEN_2(a-1); */
        /*  I_2(a) = I_2(a-1); */
        /*  a = a+1; */
        /* 'mapTF:835' new_len = a+1; */
        /* 'mapTF:836' new_PWM2 = cell(new_len,1); */
        match_idx = (int)(a + 1.0);
        i1 = new_PWM2->size[0];
        new_PWM2->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_cell_wrap_6(new_PWM2, i1);
        P_data = new_PWM2->data;
        for (i1 = 0; i1 < match_idx; i1++) {
          P_data[i1].f1->size[0] = 0;
          P_data[i1].f1->size[1] = 0;
        }

        /* 'mapTF:837' new_PWM2 = coder.nullcopy(new_PWM2); */
        i1 = cur_PWM_tmp->size[0];
        cur_PWM_tmp->size[0] = new_PWM2->size[0];
        emxEnsureCapacity_cell_wrap_6(cur_PWM_tmp, i1);
        cur_PWM_tmp_data = cur_PWM_tmp->data;

        /* 'mapTF:838' for idx = 1:new_len */
        for (unnamed_idx_0 = 0; unnamed_idx_0 < match_idx; unnamed_idx_0++) {
          /* 'mapTF:839' if idx <= a-1 */
          if ((double)unnamed_idx_0 + 1.0 <= a - 1.0) {
            /* 'mapTF:840' new_PWM2{idx} = PWM2{idx}; */
            i1 = cur_PWM_tmp_data[unnamed_idx_0].f1->size[0] *
              cur_PWM_tmp_data[unnamed_idx_0].f1->size[1];
            cur_PWM_tmp_data[unnamed_idx_0].f1->size[0] =
              PWM2_data[unnamed_idx_0].f1->size[0];
            cur_PWM_tmp_data[unnamed_idx_0].f1->size[1] =
              PWM2_data[unnamed_idx_0].f1->size[1];
            emxEnsureCapacity_real_T(cur_PWM_tmp_data[unnamed_idx_0].f1, i1);
            loop_ub = PWM2_data[unnamed_idx_0].f1->size[0] *
              PWM2_data[unnamed_idx_0].f1->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              cur_PWM_tmp_data[unnamed_idx_0].f1->data[i1] =
                PWM2_data[unnamed_idx_0].f1->data[i1];
            }
          } else if ((double)unnamed_idx_0 + 1.0 == a) {
            /* 'mapTF:841' elseif idx == a */
            /* 'mapTF:842' new_PWM2{idx} = pp{ii}; */
            i1 = cur_PWM_tmp_data[unnamed_idx_0].f1->size[0] *
              cur_PWM_tmp_data[unnamed_idx_0].f1->size[1];
            cur_PWM_tmp_data[unnamed_idx_0].f1->size[0] = cur_PWM_data[b_i]
              .f1->size[0];
            cur_PWM_tmp_data[unnamed_idx_0].f1->size[1] = cur_PWM_data[b_i]
              .f1->size[1];
            emxEnsureCapacity_real_T(cur_PWM_tmp_data[unnamed_idx_0].f1, i1);
            loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i].f1->
              size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              cur_PWM_tmp_data[unnamed_idx_0].f1->data[i1] = cur_PWM_data[b_i].
                f1->data[i1];
            }
          } else if ((double)unnamed_idx_0 + 1.0 == a + 1.0) {
            /* 'mapTF:843' elseif idx == a+1 */
            /* 'mapTF:844' new_PWM2{a+1} = rot90(pp{ii},2); */
            rot90(cur_PWM_data[b_i].f1, cur_PWM_tmp_data[(int)(a + 1.0) - 1].f1);
          }
        }

        /* 'mapTF:847' PWM2_len = length(new_PWM2); */
        /* 'mapTF:848' PWM2 = cell(PWM2_len, 1); */
        x_size[0] = cur_PWM_tmp->size[0];
        i1 = PWM2->size[0] * PWM2->size[1];
        PWM2->size[0] = cur_PWM_tmp->size[0];
        PWM2->size[1] = 1;
        emxEnsureCapacity_cell_wrap_6(PWM2, i1);
        PWM2_data = PWM2->data;
        for (i1 = 0; i1 < x_size[0]; i1++) {
          PWM2_data[i1].f1->size[0] = 0;
          PWM2_data[i1].f1->size[1] = 0;
        }

        /* 'mapTF:849' PWM2 = coder.nullcopy(PWM2); */
        /* 'mapTF:850' for cur_idx = 1:PWM2_len */
        i1 = cur_PWM_tmp->size[0];
        for (cur_idx = 0; cur_idx < i1; cur_idx++) {
          /* 'mapTF:851' PWM2{cur_idx} = new_PWM2{cur_idx}; */
          unnamed_idx_0 = PWM2_data[cur_idx].f1->size[0] * PWM2_data[cur_idx].
            f1->size[1];
          PWM2_data[cur_idx].f1->size[0] = cur_PWM_tmp_data[cur_idx].f1->size[0];
          PWM2_data[cur_idx].f1->size[1] = cur_PWM_tmp_data[cur_idx].f1->size[1];
          emxEnsureCapacity_real_T(PWM2_data[cur_idx].f1, unnamed_idx_0);
          loop_ub = cur_PWM_tmp_data[cur_idx].f1->size[0] *
            cur_PWM_tmp_data[cur_idx].f1->size[1];
          for (unnamed_idx_0 = 0; unnamed_idx_0 < loop_ub; unnamed_idx_0++) {
            PWM2_data[cur_idx].f1->data[unnamed_idx_0] =
              cur_PWM_tmp_data[cur_idx].f1->data[unnamed_idx_0];
          }
        }

        /* 'mapTF:853' LEN_2_new = zeros(a+1,1); */
        i1 = f->size[0];
        f->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        for (i1 = 0; i1 < match_idx; i1++) {
          f_data[i1] = 0.0;
        }

        /* 'mapTF:854' LEN_2_new(1:a-1) = LEN_2(1:a-1); */
        if (1.0 > a - 1.0) {
          loop_ub = 0;
        } else {
          loop_ub = (int)(a - 1.0);
        }

        for (i1 = 0; i1 < loop_ub; i1++) {
          f_data[i1] = vec_data[i1];
        }

        /* 'mapTF:855' LEN_2_new(a:a+1) = [len(ii), len(ii)]; */
        for (i1 = 0; i1 < 2; i1++) {
          LEN_2_new_tmp[i1] = a + (double)i1;
        }

        f_data[(int)LEN_2_new_tmp[0] - 1] = uid_data[b_i];
        f_data[(int)LEN_2_new_tmp[1] - 1] = uid_data[b_i];

        /* 'mapTF:857' I_2_new = zeros(a+1,1); */
        i1 = clus->size[0];
        clus->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_real_T(clus, i1);
        clus_data = clus->data;
        for (i1 = 0; i1 < match_idx; i1++) {
          clus_data[i1] = 0.0;
        }

        /* 'mapTF:858' I_2_new(1:a-1) = I_2(1:a-1); */
        if (1.0 > a - 1.0) {
          loop_ub = 0;
        } else {
          loop_ub = (int)(a - 1.0);
        }

        for (i1 = 0; i1 < loop_ub; i1++) {
          clus_data[i1] = I_2_data[i1];
        }

        /* 'mapTF:859' I_2_new(a:a+1) = [info(ii), info(ii)]; */
        clus_data[(int)LEN_2_new_tmp[0] - 1] = lasso_weight_data[b_i];
        clus_data[(int)LEN_2_new_tmp[1] - 1] = lasso_weight_data[b_i];

        /* 'mapTF:861' LEN_2 = LEN_2_new; */
        i1 = LEN_2->size[0] * LEN_2->size[1];
        LEN_2->size[0] = f->size[0];
        LEN_2->size[1] = 1;
        emxEnsureCapacity_real_T(LEN_2, i1);
        vec_data = LEN_2->data;
        loop_ub = f->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          vec_data[i1] = f_data[i1];
        }

        /* 'mapTF:862' I_2 = I_2_new; */
        i1 = I_2->size[0] * I_2->size[1];
        I_2->size[0] = clus->size[0];
        I_2->size[1] = 1;
        emxEnsureCapacity_real_T(I_2, i1);
        I_2_data = I_2->data;
        loop_ub = clus->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          I_2_data[i1] = clus_data[i1];
        }

        /* 'mapTF:863' a = a+1; */
        a++;

        /* 'mapTF:864' a = a+1; */
        a++;
      }
    }
  }

  emxFree_cell_wrap_6(&b_new_PWM2);
  emxFree_cell_wrap_6(&c_cur_PWM_tmp);
  emxFree_cell_wrap_6(&new_PWM2);
  emxFree_cell_wrap_6(&cur_PWM);
  emxFree_real_T(&vec2);
  emxFree_real_T(&uid);
  emxFree_real_T(&w);

  /* 'mapTF:869' num2 = length(strfind(fileread(memefn),'MOTIF')); */
  fileread(memefn, cur_line);
  cur_line_data = cur_line->data;
  if (cur_line->size[1] == 0) {
    for (b_i = 0; b_i < 2; b_i++) {
      x_size[b_i] = 1 - b_i;
    }
  } else {
    emxInit_int32_T(&match_out, 2);
    emxInit_int32_T(&matches, 2);
    i = matches->size[0] * matches->size[1];
    matches->size[0] = 1;
    matches->size[1] = cur_line->size[1];
    emxEnsureCapacity_int32_T(matches, i);
    matches_data = matches->data;
    match_idx = 0;
    i = cur_line->size[1];
    for (b_i = 0; b_i <= i - 5; b_i++) {
      N = 1;
      while ((N <= 5) && (cur_line_data[(b_i + N) - 1] == b_cv[N - 1])) {
        N++;
      }

      if (N > 5) {
        matches_data[match_idx] = b_i + 1;
        match_idx++;
      }
    }

    i = match_out->size[0] * match_out->size[1];
    match_out->size[0] = 1;
    match_out->size[1] = match_idx;
    emxEnsureCapacity_int32_T(match_out, i);
    match_out_data = match_out->data;
    for (b_i = 0; b_i < match_idx; b_i++) {
      match_out_data[b_i] = matches_data[b_i];
    }

    emxFree_int32_T(&matches);
    for (i = 0; i < 2; i++) {
      x_size[i] = match_out->size[i];
    }

    emxFree_int32_T(&match_out);
  }

  /* 'mapTF:870' [p,names] = getmotif(memefn,1:num2); */
  emxInit_real_T(&y, 2);
  if (x_size[1] < 1) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x_size[1];
    emxEnsureCapacity_real_T(y, i);
    zscore_data = y->data;
    loop_ub = x_size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      zscore_data[i] = (double)i + 1.0;
    }
  }

  emxInit_cell_wrap_4(&names);
  b_getmotif(memefn, y, cur_PWM_tmp, names);
  names_data = names->data;

  /* 'mapTF:871' [p,info,lenvec] = trim_pwm(p,0.25); */
  trim_pwm(cur_PWM_tmp, lasso_weight, f);
  f_data = f->data;
  cur_PWM_tmp_data = cur_PWM_tmp->data;

  /* 'mapTF:872' fid = fopen([ofn '_motifs.out'], 'w'); */
  i = cur_line->size[0] * cur_line->size[1];
  cur_line->size[0] = 1;
  cur_line->size[1] = ofn->size[1] + 11;
  emxEnsureCapacity_char_T(cur_line, i);
  cur_line_data = cur_line->data;
  loop_ub = ofn->size[1];
  emxFree_real_T(&y);
  emxFree_real_T(&lasso_weight);
  for (i = 0; i < loop_ub; i++) {
    cur_line_data[i] = ofn_data[i];
  }

  for (i = 0; i < 11; i++) {
    cur_line_data[i + ofn->size[1]] = cv[i];
  }

  fileid = cfopen(cur_line, "wb");

  /* 'mapTF:873' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF:874' fprintf("ERROR: Cannot create the combined motif file\n"); */
    printf("ERROR: Cannot create the combined motif file\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:876' a = 1; */
  a = 1.0;

  /* 'mapTF:877' for i = 1:length(PWM2) */
  u0 = PWM2->size[0];
  if (u0 < 1) {
    u0 = 1;
  }

  unnamed_idx_0 = cur_PWM_tmp->size[0] + 1;
  i = cur_PWM_tmp->size[0];
  loop_ub = f->size[0];
  emxInit_cell_wrap_6(&d_cur_PWM_tmp, 1);
  for (b_i = 0; b_i < u0; b_i++) {
    /*  [ind, r] = ppmsim([PWM2{i};p], [LEN_2(i);lenvec]); */
    /* 'mapTF:879' p_len = length(p); */
    /* 'mapTF:880' cur_PWM_tmp = cell(p_len+1,1); */
    i1 = d_cur_PWM_tmp->size[0];
    d_cur_PWM_tmp->size[0] = unnamed_idx_0;
    emxEnsureCapacity_cell_wrap_6(d_cur_PWM_tmp, i1);
    P_data = d_cur_PWM_tmp->data;
    for (i1 = 0; i1 < unnamed_idx_0; i1++) {
      P_data[i1].f1->size[0] = 0;
      P_data[i1].f1->size[1] = 0;
    }

    /* 'mapTF:881' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
    i1 = b_cur_PWM_tmp->size[0];
    b_cur_PWM_tmp->size[0] = d_cur_PWM_tmp->size[0];
    emxEnsureCapacity_cell_wrap_6(b_cur_PWM_tmp, i1);
    P_data = b_cur_PWM_tmp->data;

    /* 'mapTF:882' cur_PWM_tmp{1} = PWM2{i}; */
    i1 = P_data[0].f1->size[0] * P_data[0].f1->size[1];
    P_data[0].f1->size[0] = PWM2_data[b_i].f1->size[0];
    P_data[0].f1->size[1] = PWM2_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(P_data[0].f1, i1);
    match_idx = PWM2_data[b_i].f1->size[0] * PWM2_data[b_i].f1->size[1];
    for (i1 = 0; i1 < match_idx; i1++) {
      P_data[0].f1->data[i1] = PWM2_data[b_i].f1->data[i1];
    }

    /* 'mapTF:883' for cur_idx=1:p_len */
    for (cur_idx = 0; cur_idx < i; cur_idx++) {
      /* 'mapTF:884' cur_PWM_tmp{cur_idx+1} = p{cur_idx}; */
      i1 = P_data[cur_idx + 1].f1->size[0] * P_data[cur_idx + 1].f1->size[1];
      P_data[cur_idx + 1].f1->size[0] = cur_PWM_tmp_data[cur_idx].f1->size[0];
      P_data[cur_idx + 1].f1->size[1] = cur_PWM_tmp_data[cur_idx].f1->size[1];
      emxEnsureCapacity_real_T(P_data[cur_idx + 1].f1, i1);
      match_idx = cur_PWM_tmp_data[cur_idx].f1->size[0] *
        cur_PWM_tmp_data[cur_idx].f1->size[1];
      for (i1 = 0; i1 < match_idx; i1++) {
        P_data[cur_idx + 1].f1->data[i1] = cur_PWM_tmp_data[cur_idx].f1->data[i1];
      }
    }

    /* 'mapTF:886' [ind, r] = ppmsim(cur_PWM_tmp, [LEN_2(i);lenvec]); */
    i1 = clus->size[0];
    clus->size[0] = f->size[0] + 1;
    emxEnsureCapacity_real_T(clus, i1);
    clus_data = clus->data;
    clus_data[0] = vec_data[b_i];
    for (i1 = 0; i1 < loop_ub; i1++) {
      clus_data[i1 + 1] = f_data[i1];
    }

    ppmsim(b_cur_PWM_tmp, clus, &curr_pos, &idx);

    /* 'mapTF:888' if r > 0.80 */
    if (idx > 0.8) {
      /* 'mapTF:889' fprintf(fid,'MOTIF %d\n%s\n%d\n', int32(a), names{ind}, int32(LEN_2(i))); */
      i1 = cur_line->size[0] * cur_line->size[1];
      cur_line->size[0] = 1;
      cur_line->size[1] = names_data[(int)curr_pos - 1].f1->size[1] + 1;
      emxEnsureCapacity_char_T(cur_line, i1);
      cur_line_data = cur_line->data;
      match_idx = names_data[(int)curr_pos - 1].f1->size[1];
      for (i1 = 0; i1 < match_idx; i1++) {
        cur_line_data[i1] = names_data[(int)curr_pos - 1].f1->data[i1];
      }

      cur_line_data[names_data[(int)curr_pos - 1].f1->size[1]] = '\x00';
      b_NULL = NULL;
      getfilestar(fileid, &filestar, &hal);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "MOTIF %d\n%s\n%d\n", (int)a, &cur_line_data[0], (int)
                rt_roundd_snf(vec_data[b_i]));
        if (hal) {
          fflush(filestar);
        }
      }

      /* 'mapTF:890' a = a+1; */
      a++;

      /* 'mapTF:891' for j = 1:LEN_2(i) */
      i1 = (int)vec_data[b_i];
      if (0 <= i1 - 1) {
        d_NULL = NULL;
      }

      for (N = 0; N < i1; N++) {
        /* 'mapTF:892' fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4)); */
        print_processing(PWM2_data[b_i].f1->data[N], PWM2_data[b_i].f1->data[N +
                         PWM2_data[b_i].f1->size[0]], PWM2_data[b_i].f1->data[N
                         + PWM2_data[b_i].f1->size[0] * 2], PWM2_data[b_i]
                         .f1->data[N + PWM2_data[b_i].f1->size[0] * 3],
                         validatedHoleFilling);
        getfilestar(fileid, &filestar, &hal);
        if (!(filestar == d_NULL)) {
          fprintf(filestar, "%0.3f %0.3f %0.3f %0.3f\n", validatedHoleFilling[0],
                  validatedHoleFilling[1], validatedHoleFilling[2],
                  validatedHoleFilling[3]);
          if (hal) {
            fflush(filestar);
          }
        }
      }

      /* 'mapTF:894' fprintf(fid, '\n'); */
      b_NULL = NULL;
      getfilestar(fileid, &filestar, &hal);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "\n");
        if (hal) {
          fflush(filestar);
        }
      }
    } else {
      curr_pos = vec_data[b_i];
      if (I_2_data[b_i] / curr_pos > 1.0) {
        /* 'mapTF:895' elseif I_2(i)/LEN_2(i) > 1 */
        /* 'mapTF:896' fprintf(fid,'MOTIF %d\n%s\n%d\n', int32(a), consen(PWM2{i}, LEN_2(i)), int32(LEN_2(i))); */
        b_NULL = NULL;
        getfilestar(fileid, &filestar, &hal);
        if (!(filestar == b_NULL)) {
          fprintf(filestar, "MOTIF %d\n%s\n%d\n", (int)a, "", (int)rt_roundd_snf
                  (curr_pos));
          if (hal) {
            fflush(filestar);
          }
        }

        /* 'mapTF:897' a = a+1; */
        a++;

        /* 'mapTF:898' for j = 1:LEN_2(i) */
        i1 = (int)curr_pos;
        if (0 <= (int)curr_pos - 1) {
          c_NULL = NULL;
        }

        for (N = 0; N < i1; N++) {
          /* 'mapTF:899' fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4)); */
          print_processing(PWM2_data[b_i].f1->data[N], PWM2_data[b_i].f1->data[N
                           + PWM2_data[b_i].f1->size[0]], PWM2_data[b_i]
                           .f1->data[N + PWM2_data[b_i].f1->size[0] * 2],
                           PWM2_data[b_i].f1->data[N + PWM2_data[b_i].f1->size[0]
                           * 3], validatedHoleFilling);
          getfilestar(fileid, &filestar, &hal);
          if (!(filestar == c_NULL)) {
            fprintf(filestar, "%0.3f %0.3f %0.3f %0.3f\n", validatedHoleFilling
                    [0], validatedHoleFilling[1], validatedHoleFilling[2],
                    validatedHoleFilling[3]);
            if (hal) {
              fflush(filestar);
            }
          }
        }

        /* 'mapTF:901' fprintf(fid, '\n'); */
        b_NULL = NULL;
        getfilestar(fileid, &filestar, &hal);
        if (!(filestar == b_NULL)) {
          fprintf(filestar, "\n");
          if (hal) {
            fflush(filestar);
          }
        }
      }
    }
  }

  emxFree_cell_wrap_6(&d_cur_PWM_tmp);
  emxFree_cell_wrap_6(&cur_PWM_tmp);
  emxFree_cell_wrap_4(&names);
  emxFree_cell_wrap_6(&b_cur_PWM_tmp);
  emxFree_real_T(&I_2);
  emxFree_real_T(&LEN_2);
  emxFree_cell_wrap_6(&PWM2);
  emxFree_real_T(&f);
  emxFree_char_T(&cur_line);
  emxFree_real_T(&clus);

  /* 'mapTF:904' fclose(fid); */
  cfclose(fileid);
}

/*
 * function varscore = scoreseqkmer(PWM2, lPWM2, Lmat, ss, Smat, l_svm, k_svm, ofn, dsvm)
 */
static void scoreseqkmer(const emxArray_cell_wrap_5 *PWM2, const
  emxArray_cell_wrap_5 *lPWM2, const emxArray_real_T *Lmat, const
  emxArray_real_T *ss, const emxArray_cell_wrap_6 *Smat, double l_svm, const
  emxArray_real_T *dsvm, emxArray_real_T *varscore)
{
  static const signed char varc[12] = { 2, 1, 1, 1, 3, 3, 2, 2, 4, 4, 4, 3 };

  const cell_wrap_5 *PWM2_data;
  const cell_wrap_5 *lPWM2_data;
  const cell_wrap_6 *Smat_data;
  emxArray_real_T *B;
  emxArray_real_T *a;
  emxArray_real_T *b_dsvm;
  emxArray_real_T *ind;
  emxArray_real_T *matscore;
  emxArray_real_T *scores;
  emxArray_real_T *y;
  const double *Lmat_data;
  const double *dsvm_data;
  const double *ss_data;
  double L;
  double M;
  double d;
  double evec;
  double n;
  double *B_data;
  double *a_data;
  double *ind_data;
  double *matscore_data;
  double *scores_data;
  double *varscore_data;
  double *y_data;
  int b_i;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int ii;
  int j;
  int k;
  int loop_ub;
  int loop_ub_tmp;
  int nx;
  int unnamed_idx_1;
  dsvm_data = dsvm->data;
  Smat_data = Smat->data;
  ss_data = ss->data;
  Lmat_data = Lmat->data;
  lPWM2_data = lPWM2->data;
  PWM2_data = PWM2->data;

  /* 'mapTF:358' varc = [2 3 4; 1 3 4;1 2 4;1 2 3]; */
  /* 'mapTF:359' O = ones(1,l_svm-1); */
  /* 'mapTF:360' n = numel(Lmat)/4; */
  n = (double)(Lmat->size[0] << 2) / 4.0;

  /* 'mapTF:361' varscore = zeros(n,1); */
  i = varscore->size[0];
  i1 = (int)n;
  varscore->size[0] = (int)n;
  emxEnsureCapacity_real_T(varscore, i);
  varscore_data = varscore->data;

  /* 'mapTF:362' for i = 1:n */
  if (0 <= (int)n - 1) {
    unnamed_idx_1 = (int)(l_svm - 1.0);
    loop_ub = (int)(l_svm - 1.0);
  }

  emxInit_real_T(&ind, 2);
  emxInit_real_T(&matscore, 2);
  emxInit_real_T(&scores, 1);
  emxInit_real_T(&y, 1);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&B, 2);
  emxInit_real_T(&b_dsvm, 2);
  for (b_i = 0; b_i < i1; b_i++) {
    /* 'mapTF:363' M = Lmat(i,1); */
    M = Lmat_data[b_i];

    /* 'mapTF:364' ind = [O ss(Lmat(i,2):Lmat(i,3)) O]; */
    d = Lmat_data[b_i + Lmat->size[0]];
    n = Lmat_data[b_i + Lmat->size[0] * 2];
    if (d > n) {
      i = -1;
      i2 = -1;
    } else {
      i = (int)d - 2;
      i2 = (int)n - 1;
    }

    i3 = ind->size[0] * ind->size[1];
    ind->size[0] = 1;
    ind->size[1] = (((int)(l_svm - 1.0) + i2) - i) + unnamed_idx_1;
    emxEnsureCapacity_real_T(ind, i3);
    ind_data = ind->data;
    nx = (int)(l_svm - 1.0);
    for (i3 = 0; i3 < nx; i3++) {
      ind_data[i3] = 1.0;
    }

    nx = i2 - i;
    for (i3 = 0; i3 < nx; i3++) {
      ind_data[i3 + (int)(l_svm - 1.0)] = ss_data[(i + i3) + 1];
    }

    for (i3 = 0; i3 < loop_ub; i3++) {
      ind_data[((i3 + (int)(l_svm - 1.0)) + i2) - i] = 1.0;
    }

    /* 'mapTF:365' DSVM = dsvm(Lmat(i,2):Lmat(i,3),:); */
    d = Lmat_data[b_i + Lmat->size[0]];
    n = Lmat_data[b_i + Lmat->size[0] * 2];
    if (d > n) {
      i = -1;
      i2 = -1;
    } else {
      i = (int)d - 2;
      i2 = (int)n - 1;
    }

    /* 'mapTF:366' L = Lmat(i,3)-Lmat(i,2)+1; */
    L = (n - d) + 1.0;

    /* 'mapTF:367' matscore = zeros(L,3); */
    i3 = (int)((n - d) + 1.0);
    i4 = matscore->size[0] * matscore->size[1];
    matscore->size[0] = (int)L;
    matscore->size[1] = 3;
    emxEnsureCapacity_real_T(matscore, i4);
    matscore_data = matscore->data;
    nx = (int)L * 3;
    for (i4 = 0; i4 < nx; i4++) {
      matscore_data[i4] = 0.0;
    }

    /* 'mapTF:368' for ii = 1:L */
    if (0 <= (int)L - 1) {
      loop_ub_tmp = (int)(2.0 * l_svm - 1.0);
      i5 = (int)l_svm;
    }

    for (ii = 0; ii < i3; ii++) {
      /* 'mapTF:369' Lind = l_svm-1+ii; */
      n = (l_svm - 1.0) + ((double)ii + 1.0);

      /* 'mapTF:370' scores = zeros(2*l_svm-1,1); */
      i4 = scores->size[0];
      scores->size[0] = (int)(2.0 * l_svm - 1.0);
      emxEnsureCapacity_real_T(scores, i4);
      scores_data = scores->data;
      for (i4 = 0; i4 < loop_ub_tmp; i4++) {
        scores_data[i4] = 0.0;
      }

      /* 'mapTF:371' for j = 1:2*l_svm-1 */
      for (j = 0; j < loop_ub_tmp; j++) {
        /* 'mapTF:372' scores(j) = lPWM2{M}(Lind-l_svm+j,ind(Lind-l_svm+j)); */
        nx = (int)((n - l_svm) + ((double)j + 1.0));
        scores_data[j] = lPWM2_data[(int)M - 1].f1->data[(nx + lPWM2_data[(int)M
          - 1].f1->size[0] * ((int)ind_data[nx - 1] - 1)) - 1];
      }

      /* 'mapTF:374' evec = 0; */
      evec = 0.0;

      /* 'mapTF:375' scores(l_svm) = 0; */
      scores_data[(int)l_svm - 1] = 0.0;

      /* 'mapTF:376' for j = 1:l_svm */
      for (j = 0; j < i5; j++) {
        /* 'mapTF:377' evec = evec+sum(exp(Smat{l_svm-j+1}*scores(j:l_svm+j-1))); */
        d = (l_svm + ((double)j + 1.0)) - 1.0;
        if ((double)j + 1.0 > d) {
          i4 = 0;
          i6 = 0;
        } else {
          i4 = j;
          i6 = (int)d;
        }

        i7 = a->size[0] * a->size[1];
        k = (int)((l_svm - ((double)j + 1.0)) + 1.0) - 1;
        a->size[0] = Smat_data[k].f1->size[0];
        a->size[1] = Smat_data[k].f1->size[1];
        emxEnsureCapacity_real_T(a, i7);
        a_data = a->data;
        nx = Smat_data[k].f1->size[0] * Smat_data[k].f1->size[1];
        for (i7 = 0; i7 < nx; i7++) {
          a_data[i7] = Smat_data[k].f1->data[i7];
        }

        i7 = B->size[0] * B->size[1];
        B->size[0] = 1;
        nx = i6 - i4;
        B->size[1] = nx;
        emxEnsureCapacity_real_T(B, i7);
        B_data = B->data;
        for (i6 = 0; i6 < nx; i6++) {
          B_data[i6] = scores_data[i4 + i6];
        }

        if ((Smat_data[k].f1->size[0] == 0) || (Smat_data[k].f1->size[1] == 0) ||
            (nx == 0)) {
          i4 = y->size[0];
          y->size[0] = Smat_data[k].f1->size[0];
          emxEnsureCapacity_real_T(y, i4);
          y_data = y->data;
          nx = Smat_data[k].f1->size[0];
          for (i4 = 0; i4 < nx; i4++) {
            y_data[i4] = 0.0;
          }
        } else {
          i4 = y->size[0];
          y->size[0] = Smat_data[k].f1->size[0];
          emxEnsureCapacity_real_T(y, i4);
          y_data = y->data;
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)
                      Smat_data[(int)((l_svm - ((double)j + 1.0)) + 1.0) - 1].
                      f1->size[0], (blasint)1, (blasint)Smat_data[(int)((l_svm -
            ((double)j + 1.0)) + 1.0) - 1].f1->size[1], 1.0, &a_data[0],
                      (blasint)Smat_data[(int)((l_svm - ((double)j + 1.0)) + 1.0)
                      - 1].f1->size[0], &B_data[0], (blasint)1, 0.0, &y_data[0],
                      (blasint)Smat_data[(int)((l_svm - ((double)j + 1.0)) + 1.0)
                      - 1].f1->size[0]);
        }

        nx = y->size[0];
        for (k = 0; k < nx; k++) {
          y_data[k] = exp(y_data[k]);
        }

        evec += blockedSummation(y, y->size[0]);
      }

      /* 'mapTF:379' for iii = 1:3 */
      for (k = 0; k < 3; k++) {
        /* 'mapTF:380' V = PWM2{M}(Lind,varc(ind(Lind),iii))-PWM2{M}(Lind,ind(Lind)); */
        /* 'mapTF:381' matscore(ii,iii) = evec*V; */
        nx = (int)ind_data[(int)n - 1];
        matscore_data[ii + matscore->size[0] * k] = evec * (PWM2_data[(int)M - 1]
          .f1->data[((int)n + PWM2_data[(int)M - 1].f1->size[0] * (varc[(nx + (k
          << 2)) - 1] - 1)) - 1] - PWM2_data[(int)M - 1].f1->data[((int)n +
          PWM2_data[(int)M - 1].f1->size[0] * (nx - 1)) - 1]);
      }
    }

    /* 'mapTF:384' varscore(i) = ip(matscore(:), DSVM(:)); */
    i3 = scores->size[0];
    scores->size[0] = matscore->size[0] * 3;
    emxEnsureCapacity_real_T(scores, i3);
    scores_data = scores->data;
    nx = matscore->size[0] * 3;
    for (i3 = 0; i3 < nx; i3++) {
      scores_data[i3] = matscore_data[i3];
    }

    nx = i2 - i;
    i3 = b_dsvm->size[0] * b_dsvm->size[1];
    b_dsvm->size[0] = nx;
    b_dsvm->size[1] = 3;
    emxEnsureCapacity_real_T(b_dsvm, i3);
    ind_data = b_dsvm->data;
    for (i3 = 0; i3 < 3; i3++) {
      for (i4 = 0; i4 < nx; i4++) {
        ind_data[i4 + b_dsvm->size[0] * i3] = dsvm_data[((i + i4) + dsvm->size[0]
          * i3) + 1];
      }
    }

    nx *= 3;
    i3 = y->size[0];
    y->size[0] = nx;
    emxEnsureCapacity_real_T(y, i3);
    y_data = y->data;
    for (i3 = 0; i3 < nx; i3++) {
      y_data[i3] = ind_data[i3];
    }

    /* 'mapTF:404' c = x'*y/sqrt(x'*x)/sqrt(y'*y); */
    if (matscore->size[0] * 3 < 1) {
      n = 0.0;
      L = 0.0;
    } else {
      n = cblas_ddot((blasint)(matscore->size[0] * 3), &scores_data[0], (blasint)
                     1, &y_data[0], (blasint)1);
      L = cblas_ddot((blasint)(matscore->size[0] * 3), &scores_data[0], (blasint)
                     1, &scores_data[0], (blasint)1);
    }

    if (y->size[0] < 1) {
      evec = 0.0;
    } else {
      evec = cblas_ddot((blasint)((i2 - i) * 3), &y_data[0], (blasint)1,
                        &y_data[0], (blasint)1);
    }

    L = sqrt(L);
    evec = sqrt(evec);
    varscore_data[b_i] = n / L / evec;
  }

  emxFree_real_T(&b_dsvm);
  emxFree_real_T(&B);
  emxFree_real_T(&a);
  emxFree_real_T(&y);
  emxFree_real_T(&scores);
  emxFree_real_T(&matscore);
  emxFree_real_T(&ind);
}

/*
 * function [P,V,seqindmat,seqout,seq] = seq2pv(sfn, wfn, l_svm, Negative)
 */
static void seq2pv(const emxArray_char_T *sfn, const emxArray_char_T *wfn,
                   double l_svm, double Negative, emxArray_cell_wrap_0 *P,
                   emxArray_cell_wrap_1 *V, emxArray_cell_wrap_2 *seqindmat,
                   emxArray_cell_wrap_3 *seqout, emxArray_cell_wrap_4 *seq)
{
  static const signed char mat[12] = { 1, 0, 0, 0, 2, 2, 1, 1, 3, 3, 3, 2 };

  cell_wrap_0 *P_data;
  cell_wrap_1 *V_data;
  cell_wrap_2 *seqindmat_data;
  cell_wrap_3 *seqout_data;
  cell_wrap_4 *seq_data;
  cell_wrap_4 *sequences_data;
  emxArray_cell_wrap_4 *sequences;
  emxArray_char_T *b_fileid;
  emxArray_char_T *c_fileid;
  emxArray_char_T *cur_alpha;
  emxArray_char_T *cur_line;
  emxArray_char_T *cur_seq;
  emxArray_real_T *A;
  emxArray_real_T *RR;
  emxArray_real_T *S;
  emxArray_real_T *alpha;
  emxArray_real_T *b_pow;
  emxArray_real_T *b_w;
  emxArray_real_T *p;
  emxArray_real_T *pow2;
  emxArray_real_T *rs;
  emxArray_real_T *ss;
  emxArray_real_T *v;
  emxArray_real_T *w;
  creal_T dc;
  double b_j[2];
  double b_varargin_1_tmp[2];
  double varargin_1_tmp[2];
  double I2;
  double b_I;
  double curr_pos;
  double idx;
  double n;
  double *A_data;
  double *S_data;
  double *alpha_data;
  double *p_data;
  double *pow2_data;
  double *pow_data;
  double *rs_data;
  double *ss_data;
  double *v_data;
  double *w_data;
  int x_size[2];
  int b_i;
  int b_loop_ub;
  int c_loop_ub;
  int d_loop_ub;
  int exitg1;
  int i;
  int i1;
  int j;
  int l;
  int loop_ub;
  int m;
  int md2;
  int pow2_tmp;
  int unnamed_idx_0_tmp_tmp;
  int varargin_2;
  signed char fileid;
  bool exitg2;
  bool y;

  /* sfn: fasta file */
  /* wfn: kmer weight file */
  /* ofn: output prefix */
  /* 'mapTF:478' fid = fopen(wfn, 'r'); */
  fileid = cfopen(wfn, "rb");

  /* 'mapTF:479' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF:480' fprintf("ERROR: Weight file cannot be opened.\n") */
    printf("ERROR: Weight file cannot be opened.\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:483' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF:484' idx=0; */
  idx = 0.0;

  /* 'mapTF:485' while ~feof(fid) */
  emxInit_char_T(&b_fileid, 2);
  do {
    exitg1 = 0;
    I2 = b_feof(fileid);
    if (!(I2 != 0.0)) {
      /* 'mapTF:486' idx=idx+1; */
      idx++;

      /* 'mapTF:487' fgetl(fid); */
      b_fgets(fileid, b_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_fileid);
  emxInit_cell_wrap_4(&sequences);

  /* 'mapTF:489' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF:490' sequences = cell(idx, 1); */
  m = (int)idx;
  i = sequences->size[0];
  sequences->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_4(sequences, i);
  sequences_data = sequences->data;
  for (i = 0; i < m; i++) {
    sequences_data[i].f1->size[0] = 1;
    sequences_data[i].f1->size[1] = 0;
  }

  emxInit_real_T(&alpha, 1);

  /* 'mapTF:491' sequences = coder.nullcopy(sequences); */
  /* 'mapTF:492' alpha = zeros(idx, 1); */
  i = alpha->size[0];
  alpha->size[0] = (int)idx;
  emxEnsureCapacity_real_T(alpha, i);
  alpha_data = alpha->data;
  for (i = 0; i < m; i++) {
    alpha_data[i] = 0.0;
  }

  /* 'mapTF:493' for cur_idx=1:idx */
  m = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&cur_seq, 2);
  emxInit_char_T(&cur_alpha, 2);
  exitg2 = false;
  while ((!exitg2) && (m <= (int)idx - 1)) {
    /* 'mapTF:494' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF:495' if cur_line == -1 */
    for (i = 0; i < 2; i++) {
      x_size[i] = cur_line->size[i];
    }

    y = (x_size[1] != 0);
    if (y) {
      y = (0 > x_size[1] - 1);
    }

    if (y) {
      exitg2 = true;
    } else {
      /* 'mapTF:498' [cur_seq, cur_alpha] = strtok(cur_line, char(9)); */
      c_strtok(cur_line, cur_seq, cur_alpha);

      /* 'mapTF:499' alpha(cur_idx,1) = real(str2double(cur_alpha)); */
      dc = str2double(cur_alpha);
      alpha_data[m] = dc.re;

      /* 'mapTF:500' sequences{cur_idx} = (strip(cur_seq)); */
      strip(cur_seq, sequences_data[m].f1);
      m++;
    }
  }

  emxFree_char_T(&cur_alpha);
  emxFree_char_T(&cur_seq);

  /* 'mapTF:502' fclose(fid); */
  cfclose(fileid);

  /* 'mapTF:504' l = length(sequences{1}); */
  varargin_2 = sequences_data[0].f1->size[1];
  l = sequences_data[0].f1->size[1];

  /* 'mapTF:505' if l ~= l_svm */
  if (sequences_data[0].f1->size[1] != l_svm) {
    /* 'mapTF:506' fprintf("ERROR: L must be the same as the length of k-mer in the weight file\n"); */
    printf("ERROR: L must be the same as the length of k-mer in the weight file\n");
    fflush(stdout);
    exit(1);
  }

  emxInit_real_T(&w, 1);

  /* 'mapTF:509' w = zeros(4^l,1); */
  md2 = (int)rt_powd_snf(4.0, sequences_data[0].f1->size[1]);
  i = w->size[0];
  w->size[0] = (int)rt_powd_snf(4.0, sequences_data[0].f1->size[1]);
  emxEnsureCapacity_real_T(w, i);
  w_data = w->data;
  for (i = 0; i < md2; i++) {
    w_data[i] = 0.0;
  }

  /* 'mapTF:510' pow = (4.^(0:(l-1)))'; */
  emxInit_real_T(&rs, 2);
  if (sequences_data[0].f1->size[1] - 1 < 0) {
    rs->size[1] = 0;
  } else {
    i = rs->size[0] * rs->size[1];
    rs->size[0] = 1;
    rs->size[1] = sequences_data[0].f1->size[1];
    emxEnsureCapacity_real_T(rs, i);
    rs_data = rs->data;
    md2 = sequences_data[0].f1->size[1] - 1;
    for (i = 0; i <= md2; i++) {
      rs_data[i] = i;
    }
  }

  i = rs->size[0] * rs->size[1];
  rs->size[0] = 1;
  emxEnsureCapacity_real_T(rs, i);
  rs_data = rs->data;
  md2 = rs->size[1] - 1;
  for (i = 0; i <= md2; i++) {
    curr_pos = rs_data[i];
    rs_data[i] = rt_powd_snf(4.0, curr_pos);
  }

  emxInit_real_T(&b_pow, 1);
  i = b_pow->size[0];
  b_pow->size[0] = rs->size[1];
  emxEnsureCapacity_real_T(b_pow, i);
  pow_data = b_pow->data;
  md2 = rs->size[1];
  for (i = 0; i < md2; i++) {
    pow_data[i] = rs_data[i];
  }

  emxInit_real_T(&pow2, 1);

  /* 'mapTF:511' pow2 = flipud(pow); */
  i = pow2->size[0];
  pow2->size[0] = b_pow->size[0];
  emxEnsureCapacity_real_T(pow2, i);
  pow2_data = pow2->data;
  md2 = b_pow->size[0];
  for (i = 0; i < md2; i++) {
    pow2_data[i] = pow_data[i];
  }

  m = b_pow->size[0] - 1;
  md2 = b_pow->size[0] >> 1;
  for (b_i = 0; b_i < md2; b_i++) {
    curr_pos = pow2_data[b_i];
    pow2_tmp = m - b_i;
    pow2_data[b_i] = pow2_data[pow2_tmp];
    pow2_data[pow2_tmp] = curr_pos;
  }

  /* 'mapTF:512' fprintf('calculating indices\n'); */
  printf("calculating indices\n");
  fflush(stdout);

  /* 'mapTF:513' for i = 1:numel(alpha) */
  i = alpha->size[0];
  emxInit_real_T(&ss, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:514' ss = letterconvert(sequences{i}); */
    letterconvert(sequences_data[b_i].f1, ss);
    ss_data = ss->data;

    /* 'mapTF:515' rs = 3-fliplr(ss); */
    /* 'mapTF:516' w(ss*pow+1) = alpha(i); */
    if (ss->size[1] < 1) {
      curr_pos = 0.0;
    } else {
      curr_pos = cblas_ddot((blasint)ss->size[1], &ss_data[0], (blasint)1,
                            &pow_data[0], (blasint)1);
    }

    w_data[(int)(curr_pos + 1.0) - 1] = alpha_data[b_i];

    /* 'mapTF:517' w(rs*pow+1) = alpha(i); */
    fliplr(ss);
    i1 = ss->size[0] * ss->size[1];
    ss->size[0] = 1;
    emxEnsureCapacity_real_T(ss, i1);
    ss_data = ss->data;
    md2 = ss->size[1] - 1;
    for (i1 = 0; i1 <= md2; i1++) {
      ss_data[i1] = 3.0 - ss_data[i1];
    }

    if (ss->size[1] < 1) {
      curr_pos = 0.0;
    } else {
      curr_pos = cblas_ddot((blasint)ss->size[1], &ss_data[0], (blasint)1,
                            &pow_data[0], (blasint)1);
    }

    w_data[(int)(curr_pos + 1.0) - 1] = alpha_data[b_i];
  }

  /* 'mapTF:519' if Negative */
  if (Negative != 0.0) {
    /* 'mapTF:520' w = -1 * w; */
    md2 = w->size[0];
    for (i = 0; i < md2; i++) {
      w_data[i] = -w_data[i];
    }

    /* 'mapTF:521' m = -1 * mean(alpha); */
    curr_pos = -(blockedSummation(alpha, alpha->size[0]) / (double)alpha->size[0]);
  } else {
    /* 'mapTF:522' else */
    /* 'mapTF:523' m = mean(alpha); */
    curr_pos = blockedSummation(alpha, alpha->size[0]) / (double)alpha->size[0];
  }

  emxInit_real_T(&b_w, 1);

  /* 'mapTF:525' s = std(alpha); */
  idx = b_std(alpha);

  /* 'mapTF:526' W = (1/2)*(1+erf((w-m)/s/sqrt(2))); */
  i = b_w->size[0];
  b_w->size[0] = w->size[0];
  emxEnsureCapacity_real_T(b_w, i);
  alpha_data = b_w->data;
  md2 = w->size[0];
  for (i = 0; i < md2; i++) {
    alpha_data[i] = (w_data[i] - curr_pos) / idx / 1.4142135623730951;
  }

  applyScalarFunction(b_w, alpha);
  alpha_data = alpha->data;
  md2 = alpha->size[0];
  for (i = 0; i < md2; i++) {
    alpha_data[i] = 0.5 * (alpha_data[i] + 1.0);
  }

  /* 'mapTF:527' W(W>0.99) = 0.99; */
  m = alpha->size[0];
  for (b_i = 0; b_i < m; b_i++) {
    if (alpha_data[b_i] > 0.99) {
      alpha_data[b_i] = 0.99;
    }
  }

  /* 'mapTF:528' W(W<0.01) = 0.01; */
  m = alpha->size[0];
  for (b_i = 0; b_i < m; b_i++) {
    if (alpha_data[b_i] < 0.01) {
      alpha_data[b_i] = 0.01;
    }
  }

  /*  seq = importdata(sfn); */
  /* 'mapTF:531' fid = fopen(sfn, 'r'); */
  fileid = cfopen(sfn, "rb");

  /* 'mapTF:532' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF:533' fprintf("ERROR: Sequence file (.fa or .fasta) cannot be opened.\n") */
    printf("ERROR: Sequence file (.fa or .fasta) cannot be opened.\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:536' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF:537' idx=0; */
  idx = 0.0;

  /* 'mapTF:538' while ~feof(fid) */
  emxInit_char_T(&c_fileid, 2);
  do {
    exitg1 = 0;
    I2 = b_feof(fileid);
    if (!(I2 != 0.0)) {
      /* 'mapTF:539' idx=idx+1; */
      idx++;

      /* 'mapTF:540' fgetl(fid); */
      b_fgets(fileid, c_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&c_fileid);

  /* 'mapTF:542' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF:543' seq = cell(idx, 1); */
  m = (int)idx;
  i = seq->size[0];
  seq->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_4(seq, i);
  seq_data = seq->data;
  for (i = 0; i < m; i++) {
    seq_data[i].f1->size[0] = 1;
    seq_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:544' seq = coder.nullcopy(seq); */
  /* 'mapTF:545' for cur_idx=1:idx */
  m = 0;
  exitg2 = false;
  while ((!exitg2) && (m <= (int)idx - 1)) {
    /* 'mapTF:546' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF:547' if cur_line == -1 */
    for (i = 0; i < 2; i++) {
      x_size[i] = cur_line->size[i];
    }

    y = (x_size[1] != 0);
    if (y) {
      y = (0 > x_size[1] - 1);
    }

    if (y) {
      exitg2 = true;
    } else {
      /* 'mapTF:550' seq{cur_idx} = (strip(cur_line)); */
      strip(cur_line, seq_data[m].f1);
      m++;
    }
  }

  emxFree_char_T(&cur_line);

  /* 'mapTF:552' fclose(fid); */
  cfclose(fileid);

  /* 'mapTF:554' n = length(seq)/2; */
  n = (double)seq->size[0] / 2.0;

  /* 'mapTF:555' seqout = cell(n,1); */
  /* 'mapTF:556' seqindmat = cell(n,1); */
  unnamed_idx_0_tmp_tmp = (int)n;
  i = seqindmat->size[0];
  seqindmat->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_2(seqindmat, i);
  seqindmat_data = seqindmat->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    seqindmat_data[i].f1->size[0] = 0;
    seqindmat_data[i].f1->size[1] = 2;
  }

  /* 'mapTF:557' seqindmat = coder.nullcopy(seqindmat); */
  /* 'mapTF:558' fprintf('converting kmers to probabilities\n'); */
  printf("converting kmers to probabilities\n");
  fflush(stdout);

  /* 'mapTF:559' P = cell(n,1); */
  /* 'mapTF:560' for i = 1:n */
  i = seqout->size[0];
  seqout->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_3(seqout, i);
  seqout_data = seqout->data;
  i = P->size[0];
  P->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_0(P, i);
  P_data = P->data;
  if (0 <= (int)n - 1) {
    if (1 > sequences_data[0].f1->size[1]) {
      loop_ub = 0;
    } else {
      loop_ub = sequences_data[0].f1->size[1];
    }

    if (1 > sequences_data[0].f1->size[1]) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = sequences_data[0].f1->size[1];
    }
  }

  emxInit_real_T(&p, 1);
  emxInit_real_T(&A, 2);
  for (b_i = 0; b_i < unnamed_idx_0_tmp_tmp; b_i++) {
    /* 'mapTF:561' if mod(i,1000)==0 */
    if (fmod((double)b_i + 1.0, 1000.0) == 0.0) {
      /* 'mapTF:562' fprintf('%d sequences converted\n', int32(i)); */
      printf("%d sequences converted\n", b_i + 1);
      fflush(stdout);
    }

    /* 'mapTF:564' L = length(seq{2*i})-l+1; */
    for (i = 0; i < 2; i++) {
      x_size[i] = seq_data[((b_i + 1) << 1) - 1].f1->size[i];
    }

    idx = (double)(x_size[1] - varargin_2) + 1.0;

    /* 'mapTF:565' seqindmat{i} = zeros(L,2); */
    m = (int)idx;
    i = seqindmat_data[b_i].f1->size[0] * seqindmat_data[b_i].f1->size[1];
    seqindmat_data[b_i].f1->size[0] = (int)idx;
    seqindmat_data[b_i].f1->size[1] = 2;
    emxEnsureCapacity_real_T(seqindmat_data[b_i].f1, i);
    md2 = (int)idx << 1;
    for (i = 0; i < md2; i++) {
      seqindmat_data[b_i].f1->data[i] = 0.0;
    }

    /* 'mapTF:566' ss = letterconvert(seq{2*i}); */
    letterconvert(seq_data[((b_i + 1) << 1) - 1].f1, ss);
    ss_data = ss->data;

    /* 'mapTF:567' rs = 3-ss; */
    i = rs->size[0] * rs->size[1];
    rs->size[0] = 1;
    rs->size[1] = ss->size[1];
    emxEnsureCapacity_real_T(rs, i);
    rs_data = rs->data;
    md2 = ss->size[1];

    /* 'mapTF:568' seqout{i} = ss+1; */
    i = seqout_data[b_i].f1->size[0] * seqout_data[b_i].f1->size[1];
    seqout_data[b_i].f1->size[0] = 1;
    seqout_data[b_i].f1->size[1] = ss->size[1];
    emxEnsureCapacity_real_T(seqout_data[b_i].f1, i);
    for (i = 0; i < md2; i++) {
      curr_pos = ss_data[i];
      rs_data[i] = 3.0 - curr_pos;
      seqout_data[b_i].f1->data[i] = curr_pos + 1.0;
    }

    /* 'mapTF:569' p = zeros(L,1); */
    i = p->size[0];
    p->size[0] = (int)idx;
    emxEnsureCapacity_real_T(p, i);
    p_data = p->data;
    for (i = 0; i < m; i++) {
      p_data[i] = 0.0;
    }

    /* 'mapTF:570' I = ss(1:l)*pow; */
    i = A->size[0] * A->size[1];
    A->size[0] = 1;
    A->size[1] = loop_ub;
    emxEnsureCapacity_real_T(A, i);
    A_data = A->data;
    for (i = 0; i < loop_ub; i++) {
      A_data[i] = ss_data[i];
    }

    if (loop_ub < 1) {
      b_I = 0.0;
    } else {
      b_I = cblas_ddot((blasint)loop_ub, &A_data[0], (blasint)1, &pow_data[0],
                       (blasint)1);
    }

    /* 'mapTF:571' I2 = rs(1:l)*pow2; */
    i = A->size[0] * A->size[1];
    A->size[0] = 1;
    A->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(A, i);
    A_data = A->data;
    for (i = 0; i < b_loop_ub; i++) {
      A_data[i] = rs_data[i];
    }

    if (b_loop_ub < 1) {
      I2 = 0.0;
    } else {
      I2 = cblas_ddot((blasint)b_loop_ub, &A_data[0], (blasint)1, &pow2_data[0],
                      (blasint)1);
    }

    /* 'mapTF:572' seqindmat{i}(1,1) = I+1; */
    seqindmat_data[b_i].f1->data[0] = b_I + 1.0;

    /* 'mapTF:573' seqindmat{i}(1,2) = I2+1; */
    seqindmat_data[b_i].f1->data[seqindmat_data[b_i].f1->size[0]] = I2 + 1.0;

    /* 'mapTF:574' p(1) = W(I+1); */
    p_data[0] = alpha_data[(int)(b_I + 1.0) - 1];

    /* 'mapTF:575' for j = 2:L */
    i = (int)(idx + -1.0);
    for (j = 0; j < i; j++) {
      /* 'mapTF:576' I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1); */
      curr_pos = rt_powd_snf(4.0, (double)varargin_2 - 1.0);
      m = (int)((unsigned int)j + varargin_2);
      b_I = (b_I - ss_data[j]) / 4.0 + curr_pos * ss_data[m];

      /* I2 = (I2-rs(j-1))/4+4^(l-1)*rs(j+l-1); */
      /* 'mapTF:578' I2 = (I2-rs(j-1)*4^(l-1))*4+rs(j+l-1); */
      I2 = (I2 - rs_data[j] * curr_pos) * 4.0 + rs_data[m];

      /* 'mapTF:579' seqindmat{i}(j,1) = I+1; */
      seqindmat_data[b_i].f1->data[j + 1] = b_I + 1.0;

      /* 'mapTF:580' seqindmat{i}(j,2) = I2+1; */
      seqindmat_data[b_i].f1->data[(j + seqindmat_data[b_i].f1->size[0]) + 1] =
        I2 + 1.0;

      /* 'mapTF:581' p(j) = W(I+1); */
      p_data[j + 1] = alpha_data[(int)(b_I + 1.0) - 1];
    }

    /* seqindmat{i}(:,2) = flipud(seqindmat{i}(:,2)); */
    /* 'mapTF:584' P{i} = p; */
    i = P_data[b_i].f1->size[0];
    P_data[b_i].f1->size[0] = p->size[0];
    emxEnsureCapacity_real_T(P_data[b_i].f1, i);
    md2 = p->size[0];
    for (i = 0; i < md2; i++) {
      P_data[b_i].f1->data[i] = p_data[i];
    }
  }

  emxFree_real_T(&pow2);
  emxFree_real_T(&alpha);

  /* 'mapTF:587' fprintf('Running dsvm\n'); */
  printf("Running dsvm\n");
  fflush(stdout);

  /* 'mapTF:588' mat = [1 2 3;0 2 3;0 1 3;0 1 2]; */
  /* 'mapTF:589' O = ones(1,l); */
  /* 'mapTF:590' V = cell(n,1); */
  /* 'mapTF:591' for i = 1:n */
  i = V->size[0];
  V->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_1(V, i);
  V_data = V->data;
  if (0 <= (int)n - 1) {
    if (1 > sequences_data[0].f1->size[1]) {
      c_loop_ub = 0;
    } else {
      c_loop_ub = sequences_data[0].f1->size[1];
    }
  }

  emxFree_cell_wrap_4(&sequences);
  emxInit_real_T(&v, 2);
  emxInit_real_T(&RR, 2);
  pow2_data = RR->data;
  emxInit_real_T(&S, 2);
  for (b_i = 0; b_i < unnamed_idx_0_tmp_tmp; b_i++) {
    /* 'mapTF:592' if mod(i,1000)==0 */
    if (fmod((double)b_i + 1.0, 1000.0) == 0.0) {
      /* 'mapTF:593' fprintf('%d sequences converted\n', int32(i)); */
      printf("%d sequences converted\n", b_i + 1);
      fflush(stdout);
    }

    /* 'mapTF:595' L = length(seq{2*i}); */
    m = ((b_i + 1) << 1) - 1;
    pow2_tmp = seq_data[m].f1->size[1];

    /* 'mapTF:596' ss = seqout{i}-1; */
    i = ss->size[0] * ss->size[1];
    ss->size[0] = 1;
    ss->size[1] = seqout_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(ss, i);
    ss_data = ss->data;
    md2 = seqout_data[b_i].f1->size[1];
    for (i = 0; i < md2; i++) {
      ss_data[i] = seqout_data[b_i].f1->data[i] - 1.0;
    }

    /* 'mapTF:597' p = zeros(L+l-1,1); */
    md2 = (int)((double)((unsigned int)seq_data[m].f1->size[1] + varargin_2) -
                1.0);
    i = p->size[0];
    p->size[0] = md2;
    emxEnsureCapacity_real_T(p, i);
    p_data = p->data;
    for (i = 0; i < md2; i++) {
      p_data[i] = 0.0;
    }

    /* 'mapTF:598' I = ss(1:l)*pow; */
    i = A->size[0] * A->size[1];
    A->size[0] = 1;
    A->size[1] = c_loop_ub;
    emxEnsureCapacity_real_T(A, i);
    A_data = A->data;
    for (i = 0; i < c_loop_ub; i++) {
      A_data[i] = ss_data[i];
    }

    x_size[1] = c_loop_ub;
    if (c_loop_ub < 1) {
      b_I = 0.0;
    } else {
      b_I = cblas_ddot((blasint)c_loop_ub, &A_data[0], (blasint)1, &pow_data[0],
                       (blasint)1);
    }

    /* 'mapTF:599' p(1) = w(I+1); */
    p_data[0] = w_data[(int)(b_I + 1.0) - 1];

    /* 'mapTF:600' for j = 2:L-l+1 */
    i = seq_data[m].f1->size[1] - l;
    for (j = 0; j < i; j++) {
      /* 'mapTF:601' I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1); */
      b_I = (b_I - ss_data[j]) / 4.0 + rt_powd_snf(4.0, (double)varargin_2 - 1.0)
        * ss_data[(int)((unsigned int)j + varargin_2)];

      /* 'mapTF:602' p(j) = w(I+1); */
      p_data[j + 1] = w_data[(int)(b_I + 1.0) - 1];
    }

    /* 'mapTF:604' v = zeros(L,3); */
    i = seq_data[m].f1->size[1];
    i1 = v->size[0] * v->size[1];
    v->size[0] = i;
    v->size[1] = 3;
    emxEnsureCapacity_real_T(v, i1);
    v_data = v->data;
    md2 = i * 3;
    for (i = 0; i < md2; i++) {
      v_data[i] = 0.0;
    }

    /* 'mapTF:605' for j = 1:L */
    i = seq_data[m].f1->size[1];
    if (0 <= i - 1) {
      varargin_1_tmp[0] = 1.0;
      b_j[1] = ((double)pow2_tmp + (double)l) - 1.0;
      if (1 > l) {
        d_loop_ub = 0;
      } else {
        d_loop_ub = varargin_2;
      }

      x_size[1] = d_loop_ub;
    }

    for (j = 0; j < i; j++) {
      /* 'mapTF:606' R = max([1 j-l+1]):min([j+l-1 L]); */
      b_varargin_1_tmp[0] = 1.0;
      b_varargin_1_tmp[1] = (double)((j - l) + 1) + 1.0;
      curr_pos = b_maximum(b_varargin_1_tmp);
      b_varargin_1_tmp[0] = (unsigned int)j + varargin_2;
      b_varargin_1_tmp[1] = (unsigned int)pow2_tmp;
      I2 = minimum(b_varargin_1_tmp);
      if (rtIsNaN(curr_pos) || rtIsNaN(I2)) {
        i1 = rs->size[0] * rs->size[1];
        rs->size[0] = 1;
        rs->size[1] = 1;
        emxEnsureCapacity_real_T(rs, i1);
        rs_data = rs->data;
        rs_data[0] = rtNaN;
      } else if (I2 < curr_pos) {
        rs->size[0] = 1;
        rs->size[1] = 0;
      } else if ((rtIsInf(curr_pos) || rtIsInf(I2)) && (curr_pos == I2)) {
        i1 = rs->size[0] * rs->size[1];
        rs->size[0] = 1;
        rs->size[1] = 1;
        emxEnsureCapacity_real_T(rs, i1);
        rs_data = rs->data;
        rs_data[0] = rtNaN;
      } else if (floor(curr_pos) == curr_pos) {
        i1 = rs->size[0] * rs->size[1];
        rs->size[0] = 1;
        md2 = (int)floor(I2 - curr_pos);
        rs->size[1] = md2 + 1;
        emxEnsureCapacity_real_T(rs, i1);
        rs_data = rs->data;
        for (i1 = 0; i1 <= md2; i1++) {
          rs_data[i1] = curr_pos + (double)i1;
        }
      } else {
        eml_float_colon(curr_pos, I2, rs);
        rs_data = rs->data;
      }

      /* 'mapTF:607' RR = max([1 j-l+1]):min([j L+l-1]); */
      varargin_1_tmp[1] = (double)((j - varargin_2) + 1) + 1.0;
      idx = b_maximum(varargin_1_tmp);
      b_j[0] = (double)j + 1.0;
      curr_pos = minimum(b_j);
      if (rtIsNaN(idx) || rtIsNaN(curr_pos)) {
        i1 = RR->size[0] * RR->size[1];
        RR->size[0] = 1;
        RR->size[1] = 1;
        emxEnsureCapacity_real_T(RR, i1);
        pow2_data = RR->data;
        pow2_data[0] = rtNaN;
      } else if (curr_pos < idx) {
        RR->size[0] = 1;
        RR->size[1] = 0;
      } else if ((rtIsInf(idx) || rtIsInf(curr_pos)) && (idx == curr_pos)) {
        i1 = RR->size[0] * RR->size[1];
        RR->size[0] = 1;
        RR->size[1] = 1;
        emxEnsureCapacity_real_T(RR, i1);
        pow2_data = RR->data;
        pow2_data[0] = rtNaN;
      } else if (floor(idx) == idx) {
        i1 = RR->size[0] * RR->size[1];
        RR->size[0] = 1;
        md2 = (int)floor(curr_pos - idx);
        RR->size[1] = md2 + 1;
        emxEnsureCapacity_real_T(RR, i1);
        pow2_data = RR->data;
        for (i1 = 0; i1 <= md2; i1++) {
          pow2_data[i1] = idx + (double)i1;
        }
      } else {
        eml_float_colon(idx, curr_pos, RR);
        pow2_data = RR->data;
      }

      /* 'mapTF:608' S = ss(R); */
      i1 = S->size[0] * S->size[1];
      S->size[0] = 1;
      S->size[1] = rs->size[1];
      emxEnsureCapacity_real_T(S, i1);
      S_data = S->data;
      md2 = rs->size[1];
      for (i1 = 0; i1 < md2; i1++) {
        S_data[i1] = ss_data[(int)rs_data[i1] - 1];
      }

      /* 'mapTF:609' ref = sum(p(RR)); */
      i1 = b_w->size[0];
      b_w->size[0] = RR->size[1];
      emxEnsureCapacity_real_T(b_w, i1);
      alpha_data = b_w->data;
      md2 = RR->size[1];
      for (i1 = 0; i1 < md2; i1++) {
        alpha_data[i1] = p_data[(int)pow2_data[i1] - 1];
      }

      curr_pos = blockedSummation(b_w, RR->size[1]);

      /* 'mapTF:610' if max([1 j-l+1]) == 1 */
      if (idx == 1.0) {
        /* 'mapTF:611' cen = j; */
        m = j + 1;

        /* 'mapTF:612' cen2 = j; */
        md2 = j + 1;
      } else if (I2 == pow2_tmp) {
        /* 'mapTF:613' elseif min([j+l-1 L]) == L */
        /* 'mapTF:614' cen = l; */
        m = varargin_2;

        /* 'mapTF:615' cen2 = L-j+1; */
        md2 = pow2_tmp - j;
      } else {
        /* 'mapTF:616' else */
        /* 'mapTF:617' cen = l; */
        m = varargin_2;

        /* 'mapTF:618' cen2 = l; */
        md2 = varargin_2;
      }

      /* 'mapTF:620' for ii = 1:3 */
      for (b_loop_ub = 0; b_loop_ub < 3; b_loop_ub++) {
        /* 'mapTF:621' S(cen) = mat(ss(j)+1,ii); */
        S_data[m - 1] = mat[((int)(ss_data[j] + 1.0) + (b_loop_ub << 2)) - 1];

        /* 'mapTF:622' I = S(1:l)*pow; */
        i1 = A->size[0] * A->size[1];
        A->size[0] = 1;
        A->size[1] = d_loop_ub;
        emxEnsureCapacity_real_T(A, i1);
        A_data = A->data;
        for (i1 = 0; i1 < d_loop_ub; i1++) {
          A_data[i1] = S_data[i1];
        }

        if (x_size[1] < 1) {
          b_I = 0.0;
        } else {
          b_I = cblas_ddot((blasint)d_loop_ub, &A_data[0], (blasint)1,
                           &pow_data[0], (blasint)1);
        }

        /* 'mapTF:623' v(j,ii) = w(I+1); */
        v_data[j + v->size[0] * b_loop_ub] = w_data[(int)(b_I + 1.0) - 1];

        /* 'mapTF:624' if length(RR) > 1 */
        if (RR->size[1] > 1) {
          /* 'mapTF:625' for jj = 2:cen2 */
          for (loop_ub = 0; loop_ub <= md2 - 2; loop_ub++) {
            /* 'mapTF:626' I = (I-S(jj-1))/4+4^(l-1)*S(jj+l-1); */
            b_I = (b_I - S_data[loop_ub]) / 4.0 + rt_powd_snf(4.0, (double)
              varargin_2 - 1.0) * S_data[(int)((unsigned int)loop_ub +
              varargin_2)];

            /* 'mapTF:627' v(j,ii)=v(j,ii)+w(I+1); */
            v_data[j + v->size[0] * b_loop_ub] += w_data[(int)(b_I + 1.0) - 1];
          }
        }

        /* 'mapTF:630' v(j,ii) = v(j,ii)-ref; */
        v_data[j + v->size[0] * b_loop_ub] -= curr_pos;
      }
    }

    /* 'mapTF:633' V{i} = v; */
    i = V_data[b_i].f1->size[0] * V_data[b_i].f1->size[1];
    V_data[b_i].f1->size[0] = v->size[0];
    V_data[b_i].f1->size[1] = 3;
    emxEnsureCapacity_real_T(V_data[b_i].f1, i);
    md2 = v->size[0] * 3;
    for (i = 0; i < md2; i++) {
      V_data[b_i].f1->data[i] = v_data[i];
    }
  }

  emxFree_real_T(&b_w);
  emxFree_real_T(&A);
  emxFree_real_T(&S);
  emxFree_real_T(&RR);
  emxFree_real_T(&v);
  emxFree_real_T(&p);
  emxFree_real_T(&rs);
  emxFree_real_T(&ss);
  emxFree_real_T(&b_pow);
  emxFree_real_T(&w);
}

/*
 * function [pp, info, len] = trim_pwm(p,cut)
 */
static void trim_pwm(emxArray_cell_wrap_6 *p, emxArray_real_T *info,
                     emxArray_real_T *len)
{
  cell_wrap_6 *p_data;
  emxArray_real_T *b_mat;
  emxArray_real_T *mat;
  emxArray_real_T *r;
  emxArray_real_T *vec;
  double y;
  double *b_mat_data;
  double *info_data;
  double *len_data;
  double *mat_data;
  double *r1;
  int b_i;
  int c_i;
  int exitg1;
  int i;
  int idx;
  int j;
  int nrows;
  int nx;
  int nxin;
  bool guard1 = false;
  p_data = p->data;

  /* 'mapTF:992' l = length(p); */
  /* 'mapTF:993' info = zeros(l, 1); */
  nx = p->size[0];
  i = info->size[0];
  info->size[0] = nx;
  emxEnsureCapacity_real_T(info, i);
  info_data = info->data;

  /* 'mapTF:994' len = zeros(l,1); */
  i = len->size[0];
  len->size[0] = nx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;

  /* 'mapTF:995' for i = 1:l */
  i = p->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&vec, 1);
  emxInit_real_T(&r, 2);
  emxInit_real_T(&b_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:996' mat = p{i}+(p{i}==0); */
    idx = mat->size[0] * mat->size[1];
    mat->size[0] = p_data[b_i].f1->size[0];
    mat->size[1] = p_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(mat, idx);
    mat_data = mat->data;
    nx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
    for (idx = 0; idx < nx; idx++) {
      mat_data[idx] = p_data[b_i].f1->data[idx] + (double)(p_data[b_i].f1->
        data[idx] == 0.0);
    }

    /* 'mapTF:997' vec = 2+sum(mat.*log(mat)/log(2),2); */
    idx = r->size[0] * r->size[1];
    r->size[0] = mat->size[0];
    r->size[1] = mat->size[1];
    emxEnsureCapacity_real_T(r, idx);
    r1 = r->data;
    nx = mat->size[0] * mat->size[1];
    for (idx = 0; idx < nx; idx++) {
      r1[idx] = mat_data[idx];
    }

    nx = mat->size[0] * mat->size[1];
    for (nrows = 0; nrows < nx; nrows++) {
      r1[nrows] = log(r1[nrows]);
    }

    if ((mat->size[0] == r->size[0]) && (mat->size[1] == r->size[1])) {
      idx = b_mat->size[0] * b_mat->size[1];
      b_mat->size[0] = mat->size[0];
      b_mat->size[1] = mat->size[1];
      emxEnsureCapacity_real_T(b_mat, idx);
      b_mat_data = b_mat->data;
      nx = mat->size[0] * mat->size[1];
      for (idx = 0; idx < nx; idx++) {
        b_mat_data[idx] = mat_data[idx] * r1[idx] / 0.69314718055994529;
      }

      sum(b_mat, vec);
      mat_data = vec->data;
    } else {
      i_binary_expand_op(vec, mat, r);
      mat_data = vec->data;
    }

    nx = vec->size[0];
    for (idx = 0; idx < nx; idx++) {
      mat_data[idx] += 2.0;
    }

    /* 'mapTF:998' while (vec(1) < cut || mean(vec(1:3)) < cut || mean(vec(2:4)) < cut) && length(vec) > 4 */
    do {
      exitg1 = 0;
      guard1 = false;
      if (mat_data[0] < 0.25) {
        guard1 = true;
      } else {
        y = mat_data[0];
        for (nrows = 0; nrows < 2; nrows++) {
          y += mat_data[nrows + 1];
        }

        if (y / 3.0 < 0.25) {
          guard1 = true;
        } else {
          y = mat_data[1];
          for (nrows = 0; nrows < 2; nrows++) {
            y += mat_data[nrows + 2];
          }

          if (y / 3.0 < 0.25) {
            guard1 = true;
          } else {
            exitg1 = 1;
          }
        }
      }

      if (guard1) {
        if (vec->size[0] > 4) {
          /* 'mapTF:999' p{i}(1,:) = []; */
          nx = p_data[b_i].f1->size[0] - 2;
          nxin = p_data[b_i].f1->size[1];
          nrows = p_data[b_i].f1->size[0] - 1;
          for (j = 0; j < nxin; j++) {
            for (c_i = 0; c_i < nrows; c_i++) {
              p_data[b_i].f1->data[c_i + p_data[b_i].f1->size[0] * j] =
                p_data[b_i].f1->data[(c_i + p_data[b_i].f1->size[0] * j) + 1];
            }
          }

          if (1 > nrows) {
            nx = 0;
          } else {
            nx++;
          }

          nxin = p_data[b_i].f1->size[1] - 1;
          for (idx = 0; idx <= nxin; idx++) {
            for (nrows = 0; nrows < nx; nrows++) {
              p_data[b_i].f1->data[nrows + nx * idx] = p_data[b_i].f1->
                data[nrows + p_data[b_i].f1->size[0] * idx];
            }
          }

          idx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
          p_data[b_i].f1->size[0] = nx;
          p_data[b_i].f1->size[1] = nxin + 1;
          emxEnsureCapacity_real_T(p_data[b_i].f1, idx);

          /* 'mapTF:1000' vec(1) = []; */
          nxin = vec->size[0];
          nx = vec->size[0];
          for (nrows = 0; nrows <= nx - 2; nrows++) {
            mat_data[nrows] = mat_data[nrows + 1];
          }

          idx = vec->size[0];
          vec->size[0] = nxin - 1;
          emxEnsureCapacity_real_T(vec, idx);
          mat_data = vec->data;
        } else {
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);

    /* 'mapTF:1002' while (vec(end) < cut || mean(vec(end-2:end)) < cut || mean(vec(end-3:end-1)) < cut) && length(vec) > 4 */
    do {
      exitg1 = 0;
      guard1 = false;
      if (mat_data[vec->size[0] - 1] < 0.25) {
        guard1 = true;
      } else {
        y = mat_data[vec->size[0] - 3];
        for (nrows = 0; nrows < 2; nrows++) {
          y += mat_data[(nrows + vec->size[0]) - 2];
        }

        if (y / 3.0 < 0.25) {
          guard1 = true;
        } else {
          y = mat_data[vec->size[0] - 4];
          for (nrows = 0; nrows < 2; nrows++) {
            y += mat_data[(nrows + vec->size[0]) - 3];
          }

          if (y / 3.0 < 0.25) {
            guard1 = true;
          } else {
            exitg1 = 1;
          }
        }
      }

      if (guard1) {
        if (vec->size[0] > 4) {
          /* 'mapTF:1003' vec(end) = []; */
          idx = vec->size[0];
          nxin = vec->size[0];
          nx = vec->size[0] - 1;
          for (nrows = idx; nrows <= nx; nrows++) {
            mat_data[nrows - 1] = mat_data[nrows];
          }

          idx = vec->size[0];
          vec->size[0] = nxin - 1;
          emxEnsureCapacity_real_T(vec, idx);
          mat_data = vec->data;

          /* 'mapTF:1004' p{i}(end,:) = []; */
          idx = p_data[b_i].f1->size[0];
          nx = p_data[b_i].f1->size[0] - 2;
          nxin = p_data[b_i].f1->size[1];
          nrows = p_data[b_i].f1->size[0] - 1;
          for (j = 0; j < nxin; j++) {
            for (c_i = idx; c_i <= nrows; c_i++) {
              p_data[b_i].f1->data[(c_i + p_data[b_i].f1->size[0] * j) - 1] =
                p_data[b_i].f1->data[c_i + p_data[b_i].f1->size[0] * j];
            }
          }

          if (1 > nrows) {
            nx = 0;
          } else {
            nx++;
          }

          nxin = p_data[b_i].f1->size[1] - 1;
          for (idx = 0; idx <= nxin; idx++) {
            for (nrows = 0; nrows < nx; nrows++) {
              p_data[b_i].f1->data[nrows + nx * idx] = p_data[b_i].f1->
                data[nrows + p_data[b_i].f1->size[0] * idx];
            }
          }

          idx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
          p_data[b_i].f1->size[0] = nx;
          p_data[b_i].f1->size[1] = nxin + 1;
          emxEnsureCapacity_real_T(p_data[b_i].f1, idx);
        } else {
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);

    /* 'mapTF:1006' info(i) = sum(vec); */
    info_data[b_i] = blockedSummation(vec, vec->size[0]);

    /* 'mapTF:1007' [len(i), ~] = size(p{i}); */
    len_data[b_i] = p_data[b_i].f1->size[0];
  }

  emxFree_real_T(&b_mat);
  emxFree_real_T(&r);
  emxFree_real_T(&vec);
  emxFree_real_T(&mat);

  /* 'mapTF:1009' pp = p; */
}

/*
 * function mapTF(varargin)
 */
void mapTF(const emxArray_char_T *varargin_1, const emxArray_char_T *varargin_2,
           const emxArray_char_T *varargin_3, const emxArray_char_T *varargin_4,
           const emxArray_char_T *varargin_5, const emxArray_char_T *varargin_6,
           double varargin_7, double varargin_8, double varargin_9, double
           varargin_10, double varargin_11, double varargin_12, double
           varargin_13)
{
  static const char b_cv[9] = { '_', 'n', 'e', 'g', 'a', 't', 'i', 'v', 'e' };

  cell_wrap_0 *P_data;
  cell_wrap_0 *VV_data;
  cell_wrap_1 *V_data;
  cell_wrap_2 *seqindmat_data;
  cell_wrap_3 *ss_data;
  cell_wrap_5 *PWM2_data;
  cell_wrap_5 *PWM_data;
  cell_wrap_5 *lPWM2_data;
  cell_wrap_5 *pwm_data;
  cell_wrap_6 *Smat_data;
  cell_wrap_7 *NN_data;
  emxArray_boolean_T *b_lab;
  emxArray_cell_wrap_0 *P;
  emxArray_cell_wrap_0 *VV;
  emxArray_cell_wrap_1 *V;
  emxArray_cell_wrap_2 *seqindmat;
  emxArray_cell_wrap_3 *ss;
  emxArray_cell_wrap_4 *names;
  emxArray_cell_wrap_4 *seq;
  emxArray_cell_wrap_5 *LL;
  emxArray_cell_wrap_5 *PWM;
  emxArray_cell_wrap_5 *PWM2;
  emxArray_cell_wrap_5 *lPWM2;
  emxArray_cell_wrap_5 *lpwm;
  emxArray_cell_wrap_5 *p;
  emxArray_cell_wrap_5 *pwm;
  emxArray_cell_wrap_6 *Smat;
  emxArray_cell_wrap_7 *NN;
  emxArray_char_T *b_ofn;
  emxArray_char_T *ofn;
  emxArray_int32_T *p1;
  emxArray_int32_T *r2;
  emxArray_real_T *LEN;
  emxArray_real_T *LEN_2;
  emxArray_real_T *b_lpwm;
  emxArray_real_T *c;
  emxArray_real_T *f;
  emxArray_real_T *ff;
  emxArray_real_T *kmat;
  emxArray_real_T *lab;
  emxArray_real_T *len;
  emxArray_real_T *minnorm;
  emxArray_real_T *p2;
  emxArray_real_T *pwm_prob;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  emxArray_real_T *seqmat;
  emxArray_real_T *shift;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  emxArray_real_T *y;
  double GCmat[4];
  double b_varargin_7[2];
  double GC;
  double GCmat_tmp;
  double b;
  double *LEN_2_data;
  double *LEN_data;
  double *c_data;
  double *f_data;
  double *ff_data;
  double *kmat_data;
  double *lab_data;
  double *p2_data;
  double *seqmat_data;
  double *shift_data;
  int b_i;
  int b_loop_ub;
  unsigned int c_i;
  int i;
  int i1;
  int i2;
  int i3;
  int j;
  int loop_ub;
  int md2;
  int nx;
  int pl;
  int sizes_idx_0;
  int varargin_7_idx_0;
  int *p1_data;
  int *r3;
  const char *varargin_6_data;
  char *b_ofn_data;
  char *ofn_data;
  bool empty_non_axis_sizes;
  bool *b_lab_data;
  if (!isInitialized_mapTF) {
    mapTF_initialize();
  }

  varargin_6_data = varargin_6->data;
  emxInit_char_T(&ofn, 2);

  /*  mapTF maps the TFBS motifs found from gkmPWM and gkmPWMlasso to regions at */
  /*      base-pair resolution */
  /*   */
  /*      mapTF(seqfile, wfile, gkmPWMmemefile, gkmPWMlassofile, memefile, outputprefix...) */
  /*   */
  /*      Requires outputs of gkmPWM and gkmPWMlasso for a given model. Any combination */
  /*      of parameters for gkmPWM and gkmPWMlasso will work. */
  /*   */
  /*      Positional Parameters (Required): */
  /*   */
  /*      seqfile         The set of sequences to which the motifs will be mapped */
  /*                      (fasta format) */
  /*      wfile           Two column tab delimited file of kmer weights.  See the */
  /*                      README.md for detail on how to generate this file. */
  /*      gkmPWMmemefile  The gkmPWM meme output file */
  /*      gkmPWMlassofile The gkmPWMlasso output file */
  /*      memefile        The collection of PWMs in meme format */
  /*      outputprefix    The prefix of the output files.   */
  /*   */
  /*      Name Value Pair Parameters (Optional): */
  /*   */
  /*      'l'             The full length of the gapped k-mer.  This NEEDS to be  */
  /*                      the same as the l in the gkmSVM model (default: 11) */
  /*      'k'             The number of ungapped positions of the gapped k-mer. */
  /*                      This NEEDS to be the same as the k in the gkmSVM model */
  /*                      (default: 7) */
  /*      'KmerFrac'      Set the fraction of the total number of gapped k-mers to */
  /*                      use with mapTF.  This reduces the memory and runtime */
  /*                      needed.  If the total number of gapped k-mers is too high */
  /*                      with the given combination of (l,k,KmerFrac), KmerFrac will */
  /*                      be automatically set to a lower value to create a more  */
  /*                      workable number of gapped k-mers */
  /*      'LS'            Speeds up mapTF for large numbers of sequences by saving  */
  /*                      the PWM probabilities for each k-mer.  Needs around 15GB */
  /*                      of RAM.  Only set this if you have enough RAM.   */
  /*                      (default: false) */
  /*  */
  /*      'PWMcorrcut'    The correlation cutoff to remove redundant motifs in the  */
  /*                      gkmPWM and gkmPWMlasso list.  Motif selection is prioritized */
  /*                      by the regression weight.  (default: 0.80) */
  /*  */
  /*      'dSVMcorrcut'   The correlation cutoff to remove motifs calls that do not fit */
  /*                      the deltaSVM scores for all variants in the predicted TFBS. */
  /*                      (default: 0.60) */
  /*  */
  /*   */
  /*      Outputs 2 files named: */
  /*      outputprefix_TFBS_locations.out */
  /*          Contains the location of the motifs from gkmPWM and gkmPWMlasso.   */
  /*          See the README.md for details on the format of the output */
  /*      outputprefix_motifs.out */
  /*          A file containing the PWMs that were mapped.  This is NOT in meme */
  /*          format and is made to work with the python script "mapTF_profile.py" */
  /*  */
  /*      Example (files in the example_files directory): */
  /*      mapTF('GM12878.fa', 'GM12878_weights.out','GM12878_10_6_0_15_denovo.meme',... */
  /*          'GM12878_10_6_30_gkmPWMlasso.out','combined_db_v4.meme', 'GM12878',... */
  /*          'l', 11, 'k', 7,'KmerFrac', 1) */
  /*          Outputs GM12878_TFBS_locations.out and GM12878_motif.out */
  /* 'mapTF:61' if nargin < 6 */
  /* 'mapTF:65' fn = varargin{1}; */
  /* 'mapTF:66' wfn = varargin{2}; */
  /* 'mapTF:67' mfn1 = varargin{3}; */
  /* 'mapTF:68' mfn2 = varargin{4}; */
  /* 'mapTF:69' memefn = varargin{5}; */
  /* 'mapTF:70' ofn = varargin{6}; */
  i = ofn->size[0] * ofn->size[1];
  ofn->size[0] = 1;
  ofn->size[1] = varargin_6->size[1];
  emxEnsureCapacity_char_T(ofn, i);
  ofn_data = ofn->data;
  loop_ub = varargin_6->size[1];
  for (i = 0; i < loop_ub; i++) {
    ofn_data[i] = varargin_6_data[i];
  }

  /* 'mapTF:71' l_svm = varargin{7}; */
  /* 'mapTF:72' k_svm = varargin{8}; */
  /* 'mapTF:73' nfrac = varargin{9}; */
  /* 'mapTF:74' LS = varargin{10}; */
  /* 'mapTF:75' PWMcorrcut = varargin{11}; */
  /* 'mapTF:76' dsvmcut = varargin{12}; */
  /* 'mapTF:77' Negative = varargin{13}; */
  /* 'mapTF:79' if Negative */
  if (varargin_13 != 0.0) {
    /* 'mapTF:80' ofn = [ofn '_negative']; */
    i = ofn->size[0] * ofn->size[1];
    ofn->size[0] = 1;
    ofn->size[1] = varargin_6->size[1] + 9;
    emxEnsureCapacity_char_T(ofn, i);
    ofn_data = ofn->data;
    loop_ub = varargin_6->size[1];
    for (i = 0; i < loop_ub; i++) {
      ofn_data[i] = varargin_6_data[i];
    }

    for (i = 0; i < 9; i++) {
      ofn_data[i + varargin_6->size[1]] = b_cv[i];
    }
  }

  emxInit_cell_wrap_0(&P);
  emxInit_cell_wrap_1(&V);
  emxInit_cell_wrap_2(&seqindmat);
  emxInit_cell_wrap_3(&ss);
  emxInit_cell_wrap_4(&seq);
  emxInit_char_T(&b_ofn, 2);

  /* 'mapTF:83' fprintf('Processing Motifs\n'); */
  printf("Processing Motifs\n");
  fflush(stdout);

  /* 'mapTF:84' process_motifs(mfn1, mfn2, memefn, ofn, PWMcorrcut, Negative); */
  process_motifs(varargin_3, varargin_4, varargin_5, ofn, varargin_11,
                 varargin_13);

  /* 'mapTF:85' mfn = [ofn '_motifs.out']; */
  /* 'mapTF:86' [P,V, seqindmat, ss, seq] = seq2pv(fn, wfn,l_svm, Negative); */
  seq2pv(varargin_1, varargin_2, varargin_7, varargin_13, P, V, seqindmat, ss,
         seq);
  ss_data = ss->data;
  seqindmat_data = seqindmat->data;
  V_data = V->data;
  P_data = P->data;

  /* 'mapTF:87' GC = countGC(ss); */
  GC = countGC(ss);

  /* 'mapTF:88' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  GCmat_tmp = 0.5 - GC / 2.0;
  GCmat[0] = GCmat_tmp;
  GCmat[1] = GC / 2.0;
  GCmat[2] = GC / 2.0;
  GCmat[3] = GCmat_tmp;

  /* 'mapTF:89' b = 1; */
  b = 1.0;

  /* 'mapTF:90' [p,names,len] = getMOTIF(mfn); */
  i = b_ofn->size[0] * b_ofn->size[1];
  b_ofn->size[0] = 1;
  b_ofn->size[1] = ofn->size[1] + 11;
  emxEnsureCapacity_char_T(b_ofn, i);
  b_ofn_data = b_ofn->data;
  loop_ub = ofn->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_ofn_data[i] = ofn_data[i];
  }

  for (i = 0; i < 11; i++) {
    b_ofn_data[i + ofn->size[1]] = cv[i];
  }

  emxInit_cell_wrap_5(&PWM);
  emxInit_cell_wrap_5(&PWM2);
  emxInit_cell_wrap_5(&lPWM2);
  emxInit_cell_wrap_5(&p);
  emxInit_cell_wrap_4(&names);
  emxInit_real_T(&len, 1);
  getMOTIF(b_ofn, p, names, len);
  c_data = len->data;
  pwm_data = p->data;

  /* 'mapTF:91' a = numel(len); */
  /* 'mapTF:92' PWM = cell(a,1); */
  nx = len->size[0];
  i = PWM->size[0];
  PWM->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_5(PWM, i);
  PWM_data = PWM->data;

  /* 'mapTF:93' PWM = coder.nullcopy(PWM); */
  /* 'mapTF:94' PWM2 = cell(a,1); */
  i = PWM2->size[0];
  PWM2->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_5(PWM2, i);
  PWM2_data = PWM2->data;

  /* 'mapTF:95' PWM2 = coder.nullcopy(PWM2); */
  /* 'mapTF:96' lPWM2 = cell(a,1); */
  i = lPWM2->size[0];
  lPWM2->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_5(lPWM2, i);
  lPWM2_data = lPWM2->data;
  emxFree_char_T(&b_ofn);
  for (i = 0; i < nx; i++) {
    PWM_data[i].f1->size[0] = 0;
    PWM_data[i].f1->size[1] = 4;
    PWM2_data[i].f1->size[0] = 0;
    PWM2_data[i].f1->size[1] = 4;
    lPWM2_data[i].f1->size[0] = 0;
    lPWM2_data[i].f1->size[1] = 4;
  }

  emxInit_real_T(&LEN, 1);

  /* 'mapTF:97' lPWM2 = coder.nullcopy(lPWM2); */
  /* 'mapTF:99' LEN = zeros(a,1); */
  i = LEN->size[0];
  LEN->size[0] = len->size[0];
  emxEnsureCapacity_real_T(LEN, i);
  LEN_data = LEN->data;
  loop_ub = len->size[0];
  for (i = 0; i < loop_ub; i++) {
    LEN_data[i] = 0.0;
  }

  emxInit_real_T(&LEN_2, 1);
  emxInit_real_T(&shift, 1);

  /* 'mapTF:100' LEN_2 = zeros(a,1); */
  /* 'mapTF:101' shift = zeros(a,1); */
  /* 'mapTF:103' for i = 1:a */
  i = len->size[0];
  i1 = shift->size[0];
  shift->size[0] = len->size[0];
  emxEnsureCapacity_real_T(shift, i1);
  shift_data = shift->data;
  i1 = LEN_2->size[0];
  LEN_2->size[0] = len->size[0];
  emxEnsureCapacity_real_T(LEN_2, i1);
  LEN_2_data = LEN_2->data;
  if (0 <= len->size[0] - 1) {
    b_varargin_7[1] = 4.0;
    varargin_7_idx_0 = (int)(varargin_7 - 1.0);
    sizes_idx_0 = (int)(varargin_7 - 1.0);
    b_loop_ub = (int)(varargin_7 - 1.0);
    pl = (int)(varargin_7 - 1.0);
  }

  emxInit_real_T(&r, 2);
  emxInit_real_T(&r1, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:104' shift(i) = max([l_svm-len(i) 4]); */
    b_varargin_7[0] = varargin_7 - c_data[b_i];
    GC = b_maximum(b_varargin_7);
    shift_data[b_i] = GC;

    /* 'mapTF:105' PWM{i} = [repmat(GCmat,shift(i), 1); p{i} ;repmat(GCmat,shift(i), 1)]; */
    repmat(GCmat, GC, r);
    ff_data = r->data;
    repmat(GCmat, GC, r1);
    f_data = r1->data;
    loop_ub = pwm_data[b_i].f1->size[0];
    i1 = PWM_data[b_i].f1->size[0] * PWM_data[b_i].f1->size[1];
    PWM_data[b_i].f1->size[0] = (r->size[0] + pwm_data[b_i].f1->size[0]) +
      r1->size[0];
    PWM_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(PWM_data[b_i].f1, i1);

    /* 'mapTF:106' PWM2{i} = [ones(l_svm-1,4)/4; p{i} ;ones(l_svm-1,4)/4]; */
    i1 = PWM2_data[b_i].f1->size[0] * PWM2_data[b_i].f1->size[1];
    PWM2_data[b_i].f1->size[0] = (varargin_7_idx_0 + pwm_data[b_i].f1->size[0])
      + sizes_idx_0;
    PWM2_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(PWM2_data[b_i].f1, i1);
    nx = r->size[0];
    md2 = r1->size[0];
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < nx; i2++) {
        PWM_data[b_i].f1->data[i2 + PWM_data[b_i].f1->size[0] * i1] = ff_data[i2
          + r->size[0] * i1];
      }

      for (i2 = 0; i2 < loop_ub; i2++) {
        PWM_data[b_i].f1->data[(i2 + r->size[0]) + PWM_data[b_i].f1->size[0] *
          i1] = pwm_data[b_i].f1->data[i2 + pwm_data[b_i].f1->size[0] * i1];
      }

      for (i2 = 0; i2 < md2; i2++) {
        PWM_data[b_i].f1->data[((i2 + r->size[0]) + pwm_data[b_i].f1->size[0]) +
          PWM_data[b_i].f1->size[0] * i1] = f_data[i2 + r1->size[0] * i1];
      }

      for (i2 = 0; i2 < b_loop_ub; i2++) {
        PWM2_data[b_i].f1->data[i2 + PWM2_data[b_i].f1->size[0] * i1] = 0.25;
      }
    }

    loop_ub = pwm_data[b_i].f1->size[0];
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        PWM2_data[b_i].f1->data[(i2 + varargin_7_idx_0) + PWM2_data[b_i]
          .f1->size[0] * i1] = pwm_data[b_i].f1->data[i2 + pwm_data[b_i]
          .f1->size[0] * i1];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < pl; i2++) {
        PWM2_data[b_i].f1->data[((i2 + varargin_7_idx_0) + pwm_data[b_i]
          .f1->size[0]) + PWM2_data[b_i].f1->size[0] * i1] = 0.25;
      }
    }

    /* 'mapTF:107' lPWM2{i} = log((PWM2{i}+10^-10)/(1+4*10^-10)); */
    i1 = lPWM2_data[b_i].f1->size[0] * lPWM2_data[b_i].f1->size[1];
    lPWM2_data[b_i].f1->size[0] = PWM2_data[b_i].f1->size[0];
    lPWM2_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(lPWM2_data[b_i].f1, i1);
    loop_ub = PWM2_data[b_i].f1->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      lPWM2_data[b_i].f1->data[i1] = (PWM2_data[b_i].f1->data[i1] + 1.0E-10) /
        1.0000000004;
    }

    nx = lPWM2_data[b_i].f1->size[0] << 2;
    for (loop_ub = 0; loop_ub < nx; loop_ub++) {
      lPWM2_data[b_i].f1->data[loop_ub] = log(lPWM2_data[b_i].f1->data[loop_ub]);
    }

    /* 'mapTF:108' LEN_2(i) = len(i); */
    LEN_2_data[b_i] = c_data[b_i];

    /* 'mapTF:109' LEN(i) = len(i)+2*shift(i)-l_svm+1; */
    GC = ((c_data[b_i] + 2.0 * shift_data[b_i]) - varargin_7) + 1.0;
    LEN_data[b_i] = GC;

    /* 'mapTF:110' for j = 1:LEN(i) */
    i1 = (int)GC;
    for (j = 0; j < i1; j++) {
      /* 'mapTF:111' b = b+1; */
      b++;
    }
  }

  emxFree_real_T(&r1);
  emxFree_real_T(&r);
  emxFree_cell_wrap_5(&p);
  emxInit_cell_wrap_5(&pwm);
  emxInit_cell_wrap_5(&lpwm);

  /* 'mapTF:115' pwm = cell(b,1); */
  nx = (int)b;
  i = pwm->size[0];
  pwm->size[0] = (int)b;
  emxEnsureCapacity_cell_wrap_5(pwm, i);
  pwm_data = pwm->data;

  /* 'mapTF:116' pwm = coder.nullcopy(pwm); */
  /* 'mapTF:117' lpwm = cell(b,1); */
  i = lpwm->size[0];
  lpwm->size[0] = (int)b;
  emxEnsureCapacity_cell_wrap_5(lpwm, i);
  PWM2_data = lpwm->data;
  for (i = 0; i < nx; i++) {
    pwm_data[i].f1->size[0] = 0;
    pwm_data[i].f1->size[1] = 4;
    PWM2_data[i].f1->size[0] = 0;
    PWM2_data[i].f1->size[1] = 4;
  }

  emxInit_real_T(&lab, 1);

  /* 'mapTF:118' lpwm = coder.nullcopy(lpwm); */
  /* 'mapTF:119' lab = zeros(b,1); */
  i = lab->size[0];
  lab->size[0] = (int)b;
  emxEnsureCapacity_real_T(lab, i);
  lab_data = lab->data;
  for (i = 0; i < nx; i++) {
    lab_data[i] = 0.0;
  }

  /* 'mapTF:121' b = 1; */
  b = 1.0;

  /* 'mapTF:122' for i = 1:a */
  i = len->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:123' for j = 1:LEN(i) */
    i1 = (int)LEN_data[b_i];
    for (j = 0; j < i1; j++) {
      /* 'mapTF:124' pwm{b} = PWM{i}(j:j+l_svm-1,:); */
      GC = (((double)j + 1.0) + varargin_7) - 1.0;
      if ((double)j + 1.0 > GC) {
        i2 = 0;
        i3 = 0;
      } else {
        i2 = j;
        i3 = (int)GC;
      }

      loop_ub = i3 - i2;
      i3 = (int)(b + (double)j) - 1;
      b_loop_ub = pwm_data[i3].f1->size[0] * pwm_data[i3].f1->size[1];
      pwm_data[(int)(b + (double)j) - 1].f1->size[0] = loop_ub;
      pwm_data[(int)(b + (double)j) - 1].f1->size[1] = 4;
      emxEnsureCapacity_real_T(pwm_data[(int)(b + (double)j) - 1].f1, b_loop_ub);
      for (b_loop_ub = 0; b_loop_ub < 4; b_loop_ub++) {
        for (md2 = 0; md2 < loop_ub; md2++) {
          pwm_data[(int)(b + (double)j) - 1].f1->data[md2 + pwm_data[i3]
            .f1->size[0] * b_loop_ub] = PWM_data[b_i].f1->data[(i2 + md2) +
            PWM_data[b_i].f1->size[0] * b_loop_ub];
        }
      }

      /* 'mapTF:125' lpwm{b} = log((pwm{b}+10^-10)/(1+4*10^-10)); */
      i2 = PWM2_data[i3].f1->size[0] * PWM2_data[i3].f1->size[1];
      PWM2_data[(int)(b + (double)j) - 1].f1->size[0] = pwm_data[i3].f1->size[0];
      PWM2_data[(int)(b + (double)j) - 1].f1->size[1] = 4;
      emxEnsureCapacity_real_T(PWM2_data[(int)(b + (double)j) - 1].f1, i2);
      loop_ub = pwm_data[i3].f1->size[0] * 4;
      for (i2 = 0; i2 < loop_ub; i2++) {
        PWM2_data[(int)(b + (double)j) - 1].f1->data[i2] = (pwm_data[i3]
          .f1->data[i2] + 1.0E-10) / 1.0000000004;
      }

      nx = PWM2_data[i3].f1->size[0] << 2;
      for (loop_ub = 0; loop_ub < nx; loop_ub++) {
        PWM2_data[(int)(b + (double)j) - 1].f1->data[loop_ub] = log(PWM2_data[i3]
          .f1->data[loop_ub]);
      }

      /* 'mapTF:126' lab(b) = i; */
      lab_data[i3] = (double)b_i + 1.0;

      /* 'mapTF:127' b = b+1; */
    }

    b += (double)(i1 - 1) + 1.0;
  }

  emxFree_cell_wrap_5(&pwm);
  emxFree_cell_wrap_5(&PWM);

  /* 'mapTF:132' lab = lab(1:(b-1)); */
  i = lab->size[0];
  if (1.0 > b - 1.0) {
    lab->size[0] = 0;
  } else {
    lab->size[0] = (int)((unsigned int)b - 1U);
  }

  emxEnsureCapacity_real_T(lab, i);
  lab_data = lab->data;

  /* 'mapTF:133' lab2 = lab; */
  /* 'mapTF:134' f = (1:(b-1))'; */
  emxInit_real_T(&y, 2);
  c_data = y->data;
  if (b - 1.0 < 1.0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)((b - 1.0) - 1.0) + 1;
    emxEnsureCapacity_real_T(y, i);
    c_data = y->data;
    loop_ub = (int)((b - 1.0) - 1.0);
    for (i = 0; i <= loop_ub; i++) {
      c_data[i] = (double)i + 1.0;
    }
  }

  emxInit_real_T(&f, 1);
  i = f->size[0];
  f->size[0] = y->size[1];
  emxEnsureCapacity_real_T(f, i);
  f_data = f->data;
  loop_ub = y->size[1];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = c_data[i];
  }

  /* 'mapTF:135' for i = 2:2:a */
  i = (int)((double)len->size[0] / 2.0);
  emxInit_real_T(&ff, 1);
  emxInit_real_T(&minnorm, 1);
  emxInit_int32_T(&r2, 1);
  emxInit_boolean_T(&b_lab, 1);
  for (b_i = 0; b_i < i; b_i++) {
    c_i = (b_i << 1) + 2U;

    /* 'mapTF:136' ff = find(lab2==i); */
    loop_ub = lab->size[0];
    i1 = b_lab->size[0];
    b_lab->size[0] = lab->size[0];
    emxEnsureCapacity_boolean_T(b_lab, i1);
    b_lab_data = b_lab->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_lab_data[i1] = ((unsigned int)lab_data[i1] == c_i);
    }

    eml_find(b_lab, r2);
    r3 = r2->data;
    i1 = ff->size[0];
    ff->size[0] = r2->size[0];
    emxEnsureCapacity_real_T(ff, i1);
    ff_data = ff->data;
    loop_ub = r2->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      ff_data[i1] = r3[i1];
    }

    /* 'mapTF:137' f(ff) = flipud(f(ff)); */
    i1 = minnorm->size[0];
    minnorm->size[0] = ff->size[0];
    emxEnsureCapacity_real_T(minnorm, i1);
    LEN_data = minnorm->data;
    loop_ub = ff->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      LEN_data[i1] = f_data[(int)ff_data[i1] - 1];
    }

    nx = ff->size[0] - 1;
    md2 = ff->size[0] >> 1;
    for (b_loop_ub = 0; b_loop_ub < md2; b_loop_ub++) {
      GC = LEN_data[b_loop_ub];
      varargin_7_idx_0 = nx - b_loop_ub;
      LEN_data[b_loop_ub] = LEN_data[varargin_7_idx_0];
      LEN_data[varargin_7_idx_0] = GC;
    }

    loop_ub = minnorm->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[(int)ff_data[i1] - 1] = LEN_data[i1];
    }
  }

  /* 'mapTF:139' p1 = find(mod(lab,2)==1); */
  i = minnorm->size[0];
  minnorm->size[0] = lab->size[0];
  emxEnsureCapacity_real_T(minnorm, i);
  LEN_data = minnorm->data;
  loop_ub = lab->size[0];
  for (i = 0; i < loop_ub; i++) {
    GC = lab_data[i];
    if (rtIsNaN(GC) || rtIsInf(GC)) {
      GCmat_tmp = rtNaN;
    } else if (GC == 0.0) {
      GCmat_tmp = 0.0;
    } else {
      GCmat_tmp = fmod(GC, 2.0);
      if (GCmat_tmp == 0.0) {
        GCmat_tmp = 0.0;
      } else if (GC < 0.0) {
        GCmat_tmp += 2.0;
      }
    }

    LEN_data[i] = GCmat_tmp;
  }

  i = b_lab->size[0];
  b_lab->size[0] = minnorm->size[0];
  emxEnsureCapacity_boolean_T(b_lab, i);
  b_lab_data = b_lab->data;
  loop_ub = minnorm->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_lab_data[i] = (LEN_data[i] == 1.0);
  }

  emxInit_int32_T(&p1, 1);
  eml_find(b_lab, r2);
  r3 = r2->data;
  i = p1->size[0];
  p1->size[0] = r2->size[0];
  emxEnsureCapacity_int32_T(p1, i);
  p1_data = p1->data;
  loop_ub = r2->size[0];
  for (i = 0; i < loop_ub; i++) {
    p1_data[i] = r3[i];
  }

  /* 'mapTF:140' ff = find(mod(lab,2)==0); */
  /* 'mapTF:141' p2 = f(ff); */
  i = minnorm->size[0];
  minnorm->size[0] = lab->size[0];
  emxEnsureCapacity_real_T(minnorm, i);
  LEN_data = minnorm->data;
  loop_ub = lab->size[0];
  for (i = 0; i < loop_ub; i++) {
    GC = lab_data[i];
    if (rtIsNaN(GC) || rtIsInf(GC)) {
      GCmat_tmp = rtNaN;
    } else if (GC == 0.0) {
      GCmat_tmp = 0.0;
    } else {
      GCmat_tmp = fmod(GC, 2.0);
      if (GCmat_tmp == 0.0) {
        GCmat_tmp = 0.0;
      } else if (GC < 0.0) {
        GCmat_tmp += 2.0;
      }
    }

    LEN_data[i] = GCmat_tmp;
  }

  i = b_lab->size[0];
  b_lab->size[0] = minnorm->size[0];
  emxEnsureCapacity_boolean_T(b_lab, i);
  b_lab_data = b_lab->data;
  loop_ub = minnorm->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_lab_data[i] = (LEN_data[i] == 0.0);
  }

  emxInit_real_T(&p2, 1);
  eml_find(b_lab, r2);
  r3 = r2->data;
  i = p2->size[0];
  p2->size[0] = r2->size[0];
  emxEnsureCapacity_real_T(p2, i);
  p2_data = p2->data;
  loop_ub = r2->size[0];
  for (i = 0; i < loop_ub; i++) {
    p2_data[i] = f_data[r3[i] - 1];
  }

  emxInit_real_T(&c, 2);
  emxInit_real_T(&seqmat, 2);
  emxInit_real_T(&kmat, 2);

  /* 'mapTF:142' pl = length(p1); */
  pl = p1->size[0] - 1;

  /* 'mapTF:145' [c,~,~,~,~,rcnum] = genIndex(l_svm,k_svm,nfrac); */
  genIndex(varargin_7, varargin_8, varargin_9, seqmat, kmat, ff, lab, c, &GC);
  seqmat_data = seqmat->data;

  /* 'mapTF:146' c2 = c(1:numel(c)/k_svm-rcnum,:); */
  GC = (double)(seqmat->size[0] * seqmat->size[1]) / varargin_8 - GC;
  if (1.0 > GC) {
    loop_ub = 0;
  } else {
    loop_ub = (int)GC;
  }

  /* 'mapTF:147' c = [c;l_svm+1-fliplr(c2)]; */
  b_loop_ub = seqmat->size[1];
  i = kmat->size[0] * kmat->size[1];
  kmat->size[0] = loop_ub;
  kmat->size[1] = seqmat->size[1];
  emxEnsureCapacity_real_T(kmat, i);
  kmat_data = kmat->data;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      kmat_data[i1 + kmat->size[0] * i] = seqmat_data[i1 + seqmat->size[0] * i];
    }
  }

  b_fliplr(kmat);
  kmat_data = kmat->data;
  loop_ub = kmat->size[0] * kmat->size[1];
  for (i = 0; i < loop_ub; i++) {
    kmat_data[i] = (varargin_7 + 1.0) - kmat_data[i];
  }

  if ((seqmat->size[0] != 0) && (seqmat->size[1] != 0)) {
    nx = seqmat->size[1];
  } else if ((kmat->size[0] != 0) && (kmat->size[1] != 0)) {
    nx = kmat->size[1];
  } else {
    nx = seqmat->size[1];
    if (kmat->size[1] > seqmat->size[1]) {
      nx = kmat->size[1];
    }
  }

  empty_non_axis_sizes = (nx == 0);
  if (empty_non_axis_sizes || ((seqmat->size[0] != 0) && (seqmat->size[1] != 0)))
  {
    md2 = seqmat->size[0];
  } else {
    md2 = 0;
  }

  if (empty_non_axis_sizes || ((kmat->size[0] != 0) && (kmat->size[1] != 0))) {
    sizes_idx_0 = kmat->size[0];
  } else {
    sizes_idx_0 = 0;
  }

  i = c->size[0] * c->size[1];
  c->size[0] = md2 + sizes_idx_0;
  c->size[1] = nx;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  for (i = 0; i < nx; i++) {
    for (i1 = 0; i1 < md2; i1++) {
      c_data[i1 + c->size[0] * i] = seqmat_data[i1 + md2 * i];
    }
  }

  for (i = 0; i < nx; i++) {
    for (i1 = 0; i1 < sizes_idx_0; i1++) {
      c_data[(i1 + md2) + c->size[0] * i] = kmat_data[i1 + sizes_idx_0 * i];
    }
  }

  /* 'mapTF:148' C = numel(c)/k_svm; */
  GC = (double)(c->size[0] * c->size[1]) / varargin_8;

  /* 'mapTF:149' seqmat = zeros(C,l_svm); */
  i = (int)GC;
  i1 = seqmat->size[0] * seqmat->size[1];
  seqmat->size[0] = (int)GC;
  varargin_7_idx_0 = (int)varargin_7;
  seqmat->size[1] = (int)varargin_7;
  emxEnsureCapacity_real_T(seqmat, i1);
  seqmat_data = seqmat->data;
  loop_ub = (int)GC * (int)varargin_7;
  for (i1 = 0; i1 < loop_ub; i1++) {
    seqmat_data[i1] = 0.0;
  }

  /* 'mapTF:150' for i = 1:C */
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:151' seqmat(i,c(i,:)) = 1; */
    loop_ub = c->size[1];
    i1 = r2->size[0];
    r2->size[0] = c->size[1];
    emxEnsureCapacity_int32_T(r2, i1);
    r3 = r2->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      r3[i1] = (int)c_data[b_i + c->size[0] * i1];
    }

    loop_ub = r2->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      seqmat_data[b_i + seqmat->size[0] * (r3[i1] - 1)] = 1.0;
    }
  }

  emxInit_cell_wrap_6(&Smat, 1);

  /* 'mapTF:153' Smat = cell(l_svm,1); */
  i = Smat->size[0];
  Smat->size[0] = (int)varargin_7;
  emxEnsureCapacity_cell_wrap_6(Smat, i);
  Smat_data = Smat->data;
  for (i = 0; i < varargin_7_idx_0; i++) {
    Smat_data[i].f1->size[0] = 0;
    Smat_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:154' Smat = coder.nullcopy(Smat); */
  /* 'mapTF:155' for i = 1:l_svm */
  for (b_i = 0; b_i < varargin_7_idx_0; b_i++) {
    /* 'mapTF:156' f = find(prod(c-i,2)==0); */
    i = kmat->size[0] * kmat->size[1];
    kmat->size[0] = c->size[0];
    kmat->size[1] = c->size[1];
    emxEnsureCapacity_real_T(kmat, i);
    kmat_data = kmat->data;
    loop_ub = c->size[0] * c->size[1];
    for (i = 0; i < loop_ub; i++) {
      kmat_data[i] = c_data[i] - ((double)b_i + 1.0);
    }

    prod(kmat, minnorm);
    LEN_data = minnorm->data;
    i = b_lab->size[0];
    b_lab->size[0] = minnorm->size[0];
    emxEnsureCapacity_boolean_T(b_lab, i);
    b_lab_data = b_lab->data;
    loop_ub = minnorm->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_lab_data[i] = (LEN_data[i] == 0.0);
    }

    eml_find(b_lab, r2);
    r3 = r2->data;
    i = f->size[0];
    f->size[0] = r2->size[0];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = r2->size[0];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = r3[i];
    }

    /* 'mapTF:157' Smat{i} = zeros(length(f), l_svm); */
    i = Smat_data[b_i].f1->size[0] * Smat_data[b_i].f1->size[1];
    Smat_data[b_i].f1->size[0] = f->size[0];
    Smat_data[b_i].f1->size[1] = (int)varargin_7;
    emxEnsureCapacity_real_T(Smat_data[b_i].f1, i);
    loop_ub = f->size[0] * (int)varargin_7;
    for (i = 0; i < loop_ub; i++) {
      Smat_data[b_i].f1->data[i] = 0.0;
    }

    /* 'mapTF:158' for j = 1:length(f) */
    i = f->size[0];
    for (j = 0; j < i; j++) {
      /* 'mapTF:159' Smat{i}(j,c(f(j),:)) = 1; */
      loop_ub = c->size[1];
      i1 = r2->size[0];
      r2->size[0] = c->size[1];
      emxEnsureCapacity_int32_T(r2, i1);
      r3 = r2->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        r3[i1] = (int)c_data[((int)f_data[j] + c->size[0] * i1) - 1];
      }

      loop_ub = r2->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        Smat_data[b_i].f1->data[j + Smat_data[b_i].f1->size[0] * (r3[i1] - 1)] =
          1.0;
      }
    }
  }

  emxFree_boolean_T(&b_lab);
  emxFree_int32_T(&r2);

  /* 'mapTF:162' L = length(ss); */
  /* 'mapTF:163' B = b-1; */
  /* 'mapTF:164' maxnorm = zeros(B,1); */
  /* 'mapTF:165' minnorm = zeros(B,1); */
  /* 'mapTF:166' vec = zeros(l_svm,1); */
  /* 'mapTF:167' IND = zeros(4^l_svm,1); */
  GC = rt_powd_snf(4.0, varargin_7);
  nx = (int)rt_powd_snf(4.0, varargin_7);
  i = lab->size[0];
  lab->size[0] = (int)GC;
  emxEnsureCapacity_real_T(lab, i);
  lab_data = lab->data;
  for (i = 0; i < nx; i++) {
    lab_data[i] = 0.0;
  }

  /* 'mapTF:168' kmat = zeros(B,4^l_svm); */
  i = kmat->size[0] * kmat->size[1];
  kmat->size[0] = (int)((unsigned int)b - 1U);
  kmat->size[1] = (int)GC;
  emxEnsureCapacity_real_T(kmat, i);
  kmat_data = kmat->data;
  loop_ub = (int)((unsigned int)b - 1U) * (int)GC;
  for (i = 0; i < loop_ub; i++) {
    kmat_data[i] = 0.0;
  }

  /* 'mapTF:169' if LS */
  if (varargin_10 != 0.0) {
    /* 'mapTF:170' kmat = zeros(B/2, 4^l_svm); */
    i = (int)((b - 1.0) / 2.0);
    i1 = kmat->size[0] * kmat->size[1];
    kmat->size[0] = i;
    kmat->size[1] = (int)GC;
    emxEnsureCapacity_real_T(kmat, i1);
    kmat_data = kmat->data;
    loop_ub = i * (int)GC;
    for (i = 0; i < loop_ub; i++) {
      kmat_data[i] = 0.0;
    }
  }

  /* 'mapTF:173' for j = 1:B */
  i = (int)(unsigned int)b;
  i1 = f->size[0];
  f->size[0] = (int)((unsigned int)b - 1U);
  emxEnsureCapacity_real_T(f, i1);
  f_data = f->data;
  i1 = minnorm->size[0];
  minnorm->size[0] = (int)((unsigned int)b - 1U);
  emxEnsureCapacity_real_T(minnorm, i1);
  LEN_data = minnorm->data;
  emxInit_real_T(&x, 2);
  emxInit_real_T(&b_lpwm, 2);
  for (j = 0; j <= i - 2; j++) {
    /* 'mapTF:174' vec = max(lpwm{j}'); */
    /* 'mapTF:175' vec2 = min(lpwm{j}'); */
    /* 'mapTF:176' maxnorm(j) = sum(exp(seqmat*vec')); */
    i1 = b_lpwm->size[0] * b_lpwm->size[1];
    b_lpwm->size[0] = 4;
    loop_ub = PWM2_data[j].f1->size[0];
    b_lpwm->size[1] = PWM2_data[j].f1->size[0];
    emxEnsureCapacity_real_T(b_lpwm, i1);
    c_data = b_lpwm->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      for (i2 = 0; i2 < 4; i2++) {
        c_data[i2 + 4 * i1] = PWM2_data[j].f1->data[i1 + PWM2_data[j].f1->size[0]
          * i2];
      }
    }

    c_maximum(b_lpwm, y);
    mtimes(seqmat, y, c);
    c_data = c->data;
    i1 = x->size[0] * x->size[1];
    x->size[0] = c->size[0];
    x->size[1] = c->size[1];
    emxEnsureCapacity_real_T(x, i1);
    shift_data = x->data;
    loop_ub = c->size[0] * c->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      shift_data[i1] = c_data[i1];
    }

    nx = c->size[0];
    for (loop_ub = 0; loop_ub < nx; loop_ub++) {
      shift_data[loop_ub] = exp(shift_data[loop_ub]);
    }

    b_sum(x, y);
    c_data = y->data;
    f_data[j] = c_data[0];

    /* 'mapTF:177' minnorm(j) = sum(exp(seqmat*vec2')); */
    i1 = b_lpwm->size[0] * b_lpwm->size[1];
    b_lpwm->size[0] = 4;
    loop_ub = PWM2_data[j].f1->size[0];
    b_lpwm->size[1] = PWM2_data[j].f1->size[0];
    emxEnsureCapacity_real_T(b_lpwm, i1);
    c_data = b_lpwm->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      for (i2 = 0; i2 < 4; i2++) {
        c_data[i2 + 4 * i1] = PWM2_data[j].f1->data[i1 + PWM2_data[j].f1->size[0]
          * i2];
      }
    }

    b_minimum(b_lpwm, y);
    c_data = y->data;
    if ((seqmat->size[0] == 0) || (seqmat->size[1] == 0) || (y->size[1] == 0)) {
      i1 = ff->size[0];
      ff->size[0] = seqmat->size[0];
      emxEnsureCapacity_real_T(ff, i1);
      ff_data = ff->data;
      loop_ub = seqmat->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        ff_data[i1] = 0.0;
      }
    } else {
      i1 = ff->size[0];
      ff->size[0] = seqmat->size[0];
      emxEnsureCapacity_real_T(ff, i1);
      ff_data = ff->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)seqmat->
                  size[0], (blasint)1, (blasint)seqmat->size[1], 1.0,
                  &seqmat_data[0], (blasint)seqmat->size[0], &c_data[0],
                  (blasint)1, 0.0, &ff_data[0], (blasint)seqmat->size[0]);
    }

    nx = ff->size[0];
    for (loop_ub = 0; loop_ub < nx; loop_ub++) {
      ff_data[loop_ub] = exp(ff_data[loop_ub]);
    }

    LEN_data[j] = blockedSummation(ff, ff->size[0]);
  }

  emxFree_real_T(&b_lpwm);
  emxFree_real_T(&x);
  emxFree_real_T(&y);
  emxFree_real_T(&c);

  /* 'mapTF:179' dnorm = maxnorm-minnorm; */
  if (f->size[0] == minnorm->size[0]) {
    loop_ub = f->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[i1] -= LEN_data[i1];
    }
  } else {
    minus(f, minnorm);
    f_data = f->data;
  }

  emxInit_real_T(&vec, 1);

  /* 'mapTF:180' vec = zeros(l_svm,1); */
  i1 = vec->size[0];
  vec->size[0] = (int)varargin_7;
  emxEnsureCapacity_real_T(vec, i1);
  shift_data = vec->data;
  for (i1 = 0; i1 < varargin_7_idx_0; i1++) {
    shift_data[i1] = 0.0;
  }

  emxInit_cell_wrap_5(&LL);

  /* 'mapTF:181' all_pwm = zeros(4, l_svm, B); */
  /* 'mapTF:182' for cur_idx=1:B */
  /* 'mapTF:186' LL = cell(length(V),1); */
  nx = V->size[0];
  i1 = LL->size[0];
  LL->size[0] = V->size[0];
  emxEnsureCapacity_cell_wrap_5(LL, i1);
  pwm_data = LL->data;
  for (i1 = 0; i1 < nx; i1++) {
    pwm_data[i1].f1->size[0] = 0;
    pwm_data[i1].f1->size[1] = 4;
  }

  emxInit_cell_wrap_0(&VV);

  /* 'mapTF:187' LL = coder.nullcopy(LL); */
  /* 'mapTF:188' VV = cell(length(V),1); */
  nx = V->size[0];
  i1 = VV->size[0];
  VV->size[0] = V->size[0];
  emxEnsureCapacity_cell_wrap_0(VV, i1);
  VV_data = VV->data;
  for (i1 = 0; i1 < nx; i1++) {
    VV_data[i1].f1->size[0] = 0;
  }

  emxInit_cell_wrap_7(&NN);

  /* 'mapTF:189' VV = coder.nullcopy(VV); */
  /* 'mapTF:190' NN = cell(length(V),1); */
  nx = V->size[0];
  i1 = NN->size[0];
  NN->size[0] = V->size[0];
  emxEnsureCapacity_cell_wrap_7(NN, i1);
  NN_data = NN->data;
  for (i1 = 0; i1 < nx; i1++) {
    NN_data[i1].f1->size[0] = 0;
  }

  /* 'mapTF:191' NN = coder.nullcopy(NN); */
  /* 'mapTF:192' fprintf('Mapping motifs\n'); */
  printf("Mapping motifs\n");
  fflush(stdout);

  /* 'mapTF:193' tic */
  tic();

  /* 'mapTF:194' for I = 1:length(ss) */
  i1 = ss->size[0];
  emxInit_real_T(&pwm_prob, 2);
  for (sizes_idx_0 = 0; sizes_idx_0 < i1; sizes_idx_0++) {
    /* 'mapTF:195' seq2 = ss{I}; */
    /*      ss_onehot = zeros(4, length(seq2)); */
    /*      for idx=1:length(seq2) */
    /*          ss_onehot(seq2(idx), idx) = 1; */
    /*      end */
    /*      pwm_prob = zeros(b-1,length(seqindmat{I})); */
    /*   */
    /*      for i = 1:length(seqindmat{I}) */
    /*          ind = seqindmat{I}(i); */
    /*          if IND(ind) == 0 */
    /*              IND(ind) = 1; */
    /*              % SEQ = seq2(i:i+l_svm-1); */
    /*              % for j = 1:B */
    /*              %     for jj = 1:l_svm */
    /*              %         vec(jj) = lpwm{j}(jj,SEQ(jj)); */
    /*              %     end */
    /*              %     kmat(j,ind) = sum(exp(seqmat*vec)); */
    /*              % end */
    /*              seq3 = ss_onehot(:,i:i+l_svm-1); */
    /*              vec3 = reshape(nonzeros(all_pwm .* seq3), l_svm, b-1); */
    /*              for j = 1:B */
    /*                  kmat(j,ind) = sum(exp(seqmat*vec3(:,j))); */
    /*              end */
    /*              kmat(:,ind) = log((kmat(:,ind)-minnorm)./dnorm); */
    /*              pwm_prob(:,i) = kmat(:,ind); */
    /*          else */
    /*              pwm_prob(:,i) = kmat(:,ind); */
    /*          end */
    /*      end */
    /* 'mapTF:225' pwm_prob = zeros(B,numel(seqindmat{I})/2); */
    i2 = pwm_prob->size[0] * pwm_prob->size[1];
    pwm_prob->size[0] = (int)(b - 1.0);
    pwm_prob->size[1] = (int)((double)(seqindmat_data[sizes_idx_0].f1->size[0] <<
      1) / 2.0);
    emxEnsureCapacity_real_T(pwm_prob, i2);
    LEN_2_data = pwm_prob->data;
    loop_ub = (int)(b - 1.0) * (int)((double)(seqindmat_data[sizes_idx_0]
      .f1->size[0] << 1) / 2.0);
    for (i2 = 0; i2 < loop_ub; i2++) {
      LEN_2_data[i2] = 0.0;
    }

    /* 'mapTF:226' for i = 1:numel(seqindmat{I})/2 */
    i2 = (int)((double)(seqindmat_data[sizes_idx_0].f1->size[0] << 1) / 2.0);
    for (b_i = 0; b_i < i2; b_i++) {
      /* 'mapTF:227' ind = seqindmat{I}(i,:); */
      /* 'mapTF:228' if LS */
      if (varargin_10 != 0.0) {
        /* 'mapTF:229' if IND(ind(1)) == 0 */
        md2 = (int)seqindmat_data[sizes_idx_0].f1->data[b_i];
        if (lab_data[md2 - 1] == 0.0) {
          /* 'mapTF:230' IND(ind(1)) = 1; */
          lab_data[md2 - 1] = 1.0;

          /* 'mapTF:231' SEQ = seq2(i:i+l_svm-1); */
          if ((double)b_i + 1.0 > (((double)b_i + 1.0) + varargin_7) - 1.0) {
            i3 = 0;
          } else {
            i3 = b_i;
          }

          /* 'mapTF:232' for j = 1:pl */
          for (j = 0; j <= pl; j++) {
            /* 'mapTF:233' for jj = 1:l_svm */
            for (b_loop_ub = 0; b_loop_ub < varargin_7_idx_0; b_loop_ub++) {
              /* 'mapTF:234' vec(jj) = lpwm{p1(j)}(jj,SEQ(jj)); */
              shift_data[b_loop_ub] = PWM2_data[p1_data[j] - 1].f1->
                data[b_loop_ub + PWM2_data[p1_data[j] - 1].f1->size[0] * ((int)
                ss_data[sizes_idx_0].f1->data[i3 + b_loop_ub] - 1)];
            }

            /* 'mapTF:236' kmat(j,ind(1)) = sum(exp(seqmat*vec)); */
            loop_ub = seqmat->size[0];
            if ((seqmat->size[0] == 0) || (seqmat->size[1] == 0) || (vec->size[0]
                 == 0)) {
              b_loop_ub = ff->size[0];
              ff->size[0] = seqmat->size[0];
              emxEnsureCapacity_real_T(ff, b_loop_ub);
              ff_data = ff->data;
              for (b_loop_ub = 0; b_loop_ub < loop_ub; b_loop_ub++) {
                ff_data[b_loop_ub] = 0.0;
              }
            } else {
              b_loop_ub = ff->size[0];
              ff->size[0] = seqmat->size[0];
              emxEnsureCapacity_real_T(ff, b_loop_ub);
              ff_data = ff->data;
              cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                          seqmat->size[0], (blasint)1, (blasint)seqmat->size[1],
                          1.0, &seqmat_data[0], (blasint)seqmat->size[0],
                          &shift_data[0], (blasint)vec->size[0], 0.0, &ff_data[0],
                          (blasint)seqmat->size[0]);
            }

            nx = ff->size[0];
            for (loop_ub = 0; loop_ub < nx; loop_ub++) {
              ff_data[loop_ub] = exp(ff_data[loop_ub]);
            }

            kmat_data[j + kmat->size[0] * (md2 - 1)] = blockedSummation(ff,
              ff->size[0]);
          }

          /* 'mapTF:238' kmat(:,ind(1)) = log((kmat(:,ind(1))-minnorm(p1))./dnorm(p1)); */
          loop_ub = kmat->size[0];
          if (kmat->size[0] == 1) {
            nx = p1->size[0];
          } else {
            nx = kmat->size[0];
          }

          if ((kmat->size[0] == p1->size[0]) && (nx == p1->size[0])) {
            i3 = ff->size[0];
            ff->size[0] = kmat->size[0];
            emxEnsureCapacity_real_T(ff, i3);
            ff_data = ff->data;
            for (i3 = 0; i3 < loop_ub; i3++) {
              ff_data[i3] = (kmat_data[i3 + kmat->size[0] * (md2 - 1)] -
                             LEN_data[p1_data[i3] - 1]) / f_data[p1_data[i3] - 1];
            }
          } else {
            b_binary_expand_op(ff, kmat, seqindmat, sizes_idx_0, b_i, minnorm,
                               p1, f);
            ff_data = ff->data;
          }

          nx = ff->size[0];
          for (loop_ub = 0; loop_ub < nx; loop_ub++) {
            ff_data[loop_ub] = log(ff_data[loop_ub]);
          }

          loop_ub = ff->size[0];
          for (i3 = 0; i3 < loop_ub; i3++) {
            kmat_data[i3 + kmat->size[0] * (md2 - 1)] = ff_data[i3];
          }
        }

        /* 'mapTF:240' pwm_prob(p1,i) = kmat(:,ind(1)); */
        loop_ub = kmat->size[0];
        for (i3 = 0; i3 < loop_ub; i3++) {
          LEN_2_data[(p1_data[i3] + pwm_prob->size[0] * b_i) - 1] = kmat_data[i3
            + kmat->size[0] * (md2 - 1)];
        }

        /* 'mapTF:241' if IND(ind(2)) == 0 */
        md2 = (int)seqindmat_data[sizes_idx_0].f1->data[b_i +
          seqindmat_data[sizes_idx_0].f1->size[0]];
        if (lab_data[md2 - 1] == 0.0) {
          /* 'mapTF:242' IND(ind(2)) = 1; */
          lab_data[md2 - 1] = 1.0;

          /* 'mapTF:243' SEQ = seq2(i:i+l_svm-1); */
          if ((double)b_i + 1.0 > (((double)b_i + 1.0) + varargin_7) - 1.0) {
            i3 = 0;
          } else {
            i3 = b_i;
          }

          /* 'mapTF:244' for j = 1:pl */
          for (j = 0; j <= pl; j++) {
            /* 'mapTF:245' for jj = 1:l_svm */
            for (b_loop_ub = 0; b_loop_ub < varargin_7_idx_0; b_loop_ub++) {
              /* 'mapTF:246' vec(jj) = lpwm{p2(j)}(jj,SEQ(jj)); */
              shift_data[b_loop_ub] = PWM2_data[(int)p2_data[j] - 1].f1->
                data[b_loop_ub + PWM2_data[(int)p2_data[j] - 1].f1->size[0] *
                ((int)ss_data[sizes_idx_0].f1->data[i3 + b_loop_ub] - 1)];
            }

            /* 'mapTF:248' kmat(j,ind(2)) = sum(exp(seqmat*vec)); */
            loop_ub = seqmat->size[0];
            if ((seqmat->size[0] == 0) || (seqmat->size[1] == 0) || (vec->size[0]
                 == 0)) {
              b_loop_ub = ff->size[0];
              ff->size[0] = seqmat->size[0];
              emxEnsureCapacity_real_T(ff, b_loop_ub);
              ff_data = ff->data;
              for (b_loop_ub = 0; b_loop_ub < loop_ub; b_loop_ub++) {
                ff_data[b_loop_ub] = 0.0;
              }
            } else {
              b_loop_ub = ff->size[0];
              ff->size[0] = seqmat->size[0];
              emxEnsureCapacity_real_T(ff, b_loop_ub);
              ff_data = ff->data;
              cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                          seqmat->size[0], (blasint)1, (blasint)seqmat->size[1],
                          1.0, &seqmat_data[0], (blasint)seqmat->size[0],
                          &shift_data[0], (blasint)vec->size[0], 0.0, &ff_data[0],
                          (blasint)seqmat->size[0]);
            }

            nx = ff->size[0];
            for (loop_ub = 0; loop_ub < nx; loop_ub++) {
              ff_data[loop_ub] = exp(ff_data[loop_ub]);
            }

            kmat_data[j + kmat->size[0] * (md2 - 1)] = blockedSummation(ff,
              ff->size[0]);
          }

          /* 'mapTF:250' kmat(:,ind(2)) = log((kmat(:,ind(2))-minnorm(p2))./dnorm(p2)); */
          loop_ub = kmat->size[0];
          if (kmat->size[0] == 1) {
            nx = p2->size[0];
          } else {
            nx = kmat->size[0];
          }

          if ((kmat->size[0] == p2->size[0]) && (nx == p2->size[0])) {
            i3 = ff->size[0];
            ff->size[0] = kmat->size[0];
            emxEnsureCapacity_real_T(ff, i3);
            ff_data = ff->data;
            for (i3 = 0; i3 < loop_ub; i3++) {
              nx = (int)p2_data[i3] - 1;
              ff_data[i3] = (kmat_data[i3 + kmat->size[0] * (md2 - 1)] -
                             LEN_data[nx]) / f_data[nx];
            }
          } else {
            binary_expand_op(ff, kmat, seqindmat, sizes_idx_0, b_i, minnorm, p2,
                             f);
            ff_data = ff->data;
          }

          nx = ff->size[0];
          for (loop_ub = 0; loop_ub < nx; loop_ub++) {
            ff_data[loop_ub] = log(ff_data[loop_ub]);
          }

          loop_ub = ff->size[0];
          for (i3 = 0; i3 < loop_ub; i3++) {
            kmat_data[i3 + kmat->size[0] * (md2 - 1)] = ff_data[i3];
          }
        }

        /* 'mapTF:252' pwm_prob(p2,i) = kmat(:,ind(2)); */
        loop_ub = kmat->size[0];
        for (i3 = 0; i3 < loop_ub; i3++) {
          LEN_2_data[((int)p2_data[i3] + pwm_prob->size[0] * b_i) - 1] =
            kmat_data[i3 + kmat->size[0] * (md2 - 1)];
        }
      } else {
        /* 'mapTF:253' else */
        /* 'mapTF:254' SEQ = seq2(i:i+l_svm-1); */
        if ((double)b_i + 1.0 > (((double)b_i + 1.0) + varargin_7) - 1.0) {
          i3 = 0;
        } else {
          i3 = b_i;
        }

        /* 'mapTF:255' for j = 1:B */
        for (j = 0; j <= i - 2; j++) {
          /* 'mapTF:256' for jj = 1:l_svm */
          for (b_loop_ub = 0; b_loop_ub < varargin_7_idx_0; b_loop_ub++) {
            /* 'mapTF:257' vec(jj) = lpwm{j}(jj,SEQ(jj)); */
            shift_data[b_loop_ub] = PWM2_data[j].f1->data[b_loop_ub +
              PWM2_data[j].f1->size[0] * ((int)ss_data[sizes_idx_0].f1->data[i3
              + b_loop_ub] - 1)];
          }

          /* 'mapTF:259' pwm_prob(j,i) = log((sum(exp(seqmat*vec))-minnorm(j))/dnorm(j)); */
          loop_ub = seqmat->size[0];
          if ((seqmat->size[0] == 0) || (seqmat->size[1] == 0) || (vec->size[0] ==
               0)) {
            b_loop_ub = ff->size[0];
            ff->size[0] = seqmat->size[0];
            emxEnsureCapacity_real_T(ff, b_loop_ub);
            ff_data = ff->data;
            for (b_loop_ub = 0; b_loop_ub < loop_ub; b_loop_ub++) {
              ff_data[b_loop_ub] = 0.0;
            }
          } else {
            b_loop_ub = ff->size[0];
            ff->size[0] = seqmat->size[0];
            emxEnsureCapacity_real_T(ff, b_loop_ub);
            ff_data = ff->data;
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                        seqmat->size[0], (blasint)1, (blasint)seqmat->size[1],
                        1.0, &seqmat_data[0], (blasint)seqmat->size[0],
                        &shift_data[0], (blasint)vec->size[0], 0.0, &ff_data[0],
                        (blasint)seqmat->size[0]);
          }

          nx = ff->size[0];
          for (loop_ub = 0; loop_ub < nx; loop_ub++) {
            ff_data[loop_ub] = exp(ff_data[loop_ub]);
          }

          LEN_2_data[j + pwm_prob->size[0] * b_i] = log((blockedSummation(ff,
            ff->size[0]) - LEN_data[j]) / f_data[j]);
        }
      }
    }

    /* 'mapTF:263' [LL{I}, NN{I}] = MAPTF(fn, ss{I}, pwm_prob, l_svm, k_svm, LEN, LEN_2, shift, P{I}, names, a, b); */
    MAPTF(ss_data[sizes_idx_0].f1, pwm_prob, varargin_7, LEN, LEN_2, shift,
          P_data[sizes_idx_0].f1, names, len->size[0], pwm_data[sizes_idx_0].f1,
          NN_data[sizes_idx_0].f1);

    /* 'mapTF:264' if numel(LL{I}) > 0 */
    if ((pwm_data[sizes_idx_0].f1->size[0] << 2) > 0) {
      /* 'mapTF:265' VV{I} = scoreseqkmer(PWM2, lPWM2, LL{I}, ss{I}, Smat, l_svm, k_svm, ofn, V{I}); */
      scoreseqkmer(PWM2, lPWM2, pwm_data[sizes_idx_0].f1, ss_data[sizes_idx_0].
                   f1, Smat, varargin_7, V_data[sizes_idx_0].f1,
                   VV_data[sizes_idx_0].f1);
    }

    /* 'mapTF:267' if mod(I,100)==0 */
    if (b_mod((double)sizes_idx_0 + 1.0, 100.0) == 0.0) {
      /* 'mapTF:268' fprintf('%d out of %d sequences done...\n', int32(I), int32(length(ss))); */
      printf("%d out of %d sequences done...\n", sizes_idx_0 + 1, ss->size[0]);
      fflush(stdout);

      /* 'mapTF:269' toc */
      toc();
    }
  }

  emxFree_real_T(&len);
  emxFree_cell_wrap_4(&names);
  emxFree_cell_wrap_2(&seqindmat);
  emxFree_cell_wrap_1(&V);
  emxFree_cell_wrap_0(&P);
  emxFree_real_T(&pwm_prob);
  emxFree_real_T(&kmat);
  emxFree_real_T(&vec);
  emxFree_real_T(&minnorm);
  emxFree_cell_wrap_6(&Smat);
  emxFree_real_T(&seqmat);
  emxFree_real_T(&p2);
  emxFree_int32_T(&p1);
  emxFree_real_T(&ff);
  emxFree_real_T(&f);
  emxFree_real_T(&lab);
  emxFree_cell_wrap_5(&lpwm);
  emxFree_real_T(&shift);
  emxFree_real_T(&LEN_2);
  emxFree_real_T(&LEN);
  emxFree_cell_wrap_5(&lPWM2);
  emxFree_cell_wrap_5(&PWM2);

  /* 'mapTF:272' fprintf('%d out of %d sequences done...\n', int32(length(ss)), int32(length(ss))); */
  printf("%d out of %d sequences done...\n", ss->size[0], ss->size[0]);
  fflush(stdout);

  /*  clear kmat */
  /* 'mapTF:274' PWM_corr(ofn, VV, NN, LL, seq, dsvmcut); */
  PWM_corr(ofn, VV, NN, LL, seq, varargin_12);
  emxFree_cell_wrap_4(&seq);
  emxFree_cell_wrap_3(&ss);
  emxFree_cell_wrap_7(&NN);
  emxFree_cell_wrap_0(&VV);
  emxFree_cell_wrap_5(&LL);
  emxFree_char_T(&ofn);
}

/* End of code generation (mapTF.c) */
