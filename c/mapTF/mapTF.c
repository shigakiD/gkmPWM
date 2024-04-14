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
#include "nonzeros.h"
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
#include <emmintrin.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Variable Definitions */
static const char cv[11] = { '_', 'm', 'o', 't', 'i', 'f', 's', '.', 'o', 'u',
  't' };

/* Function Declarations */
static void MAPTF(const emxArray_real_T *ss, const emxArray_real_T *pwm_prob,
                  double l_svm, const emxArray_cell_wrap_4 *PWM2,
                  emxArray_real_T *LEN, const emxArray_real_T *LEN_2, const
                  emxArray_real_T *shift, const emxArray_real_T *gkmprob, const
                  emxArray_cell_wrap_3 *names, double a, emxArray_real_T *Lmat,
                  emxArray_cell_wrap_3 *NAME);
static void PWM_corr(const emxArray_char_T *fn, const emxArray_cell_wrap_0 *VV,
                     const emxArray_cell_wrap_6 *NN, const emxArray_cell_wrap_4 *
                     LL, const emxArray_cell_wrap_3 *seq);
static void b_binary_expand_op(emxArray_real_T *f, const emxArray_real_T
  *all_pwm, const emxArray_int8_T *ss_onehot, int i2, int i3, int i4);
static double countGC(const emxArray_cell_wrap_2 *s);
static void d_binary_expand_op(emxArray_real_T *mat, int i1, int i, int i2, int
  i3, const emxArray_real_T *pwm_prob, int i4, int i5, const emxArray_real_T
  *pos);
static void getMOTIF(const emxArray_char_T *fn, emxArray_cell_wrap_4 *mat,
                     emxArray_cell_wrap_3 *names, emxArray_real_T *len);
static void getdenovomotif(const emxArray_char_T *filename, emxArray_cell_wrap_5
  *mat, emxArray_real_T *w);
static void i_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *mat,
  const emxArray_real_T *r1);
static void letterconvert(const emxArray_char_T *s, emxArray_real_T *en);
static void ppmsim(emxArray_cell_wrap_5 *mot, const emxArray_real_T *lenvec,
                   double *ind, double *M);
static void process_motifs(const emxArray_char_T *dfn, const emxArray_char_T
  *lfn, const emxArray_char_T *memefn, const emxArray_char_T *ofn);
static void scoreseqkmer(const emxArray_cell_wrap_4 *PWM3, const
  emxArray_cell_wrap_4 *lPWM3, const emxArray_real_T *Lmat, const
  emxArray_real_T *ss, const emxArray_cell_wrap_5 *Smat, double l_svm, const
  emxArray_real_T *dsvm, emxArray_real_T *varscore);
static void seq2pv(const emxArray_char_T *sfn, const emxArray_char_T *wfn,
                   double l_svm, emxArray_cell_wrap_0 *P, emxArray_cell_wrap_1
                   *V, emxArray_cell_wrap_0 *seqindmat, emxArray_cell_wrap_2
                   *seqout, emxArray_cell_wrap_3 *seq);
static void trim_pwm(emxArray_cell_wrap_5 *p, emxArray_real_T *info,
                     emxArray_real_T *len);

/* Function Definitions */
/*
 * function [Lmat, NAME] = MAPTF(fn, ss, GC, pwm_prob, l_svm, k_svm, ofn, PWM, PWM2, pwm, lpwm, lab, LEN, LEN_2, shift, gkmprob, names, len, a, b, rnum)
 */
static void MAPTF(const emxArray_real_T *ss, const emxArray_real_T *pwm_prob,
                  double l_svm, const emxArray_cell_wrap_4 *PWM2,
                  emxArray_real_T *LEN, const emxArray_real_T *LEN_2, const
                  emxArray_real_T *shift, const emxArray_real_T *gkmprob, const
                  emxArray_cell_wrap_3 *names, double a, emxArray_real_T *Lmat,
                  emxArray_cell_wrap_3 *NAME)
{
  __m128d r;
  const cell_wrap_3 *names_data;
  cell_wrap_3 *NAME_data;
  emxArray_boolean_T *b_pos;
  emxArray_int32_T *iidx;
  emxArray_real_T *C;
  emxArray_real_T *L2;
  emxArray_real_T *b_L2;
  emxArray_real_T *b_mat;
  emxArray_real_T *ind;
  emxArray_real_T *mat;
  emxArray_real_T *neg;
  emxArray_real_T *pos;
  emxArray_real_T *vec;
  const double *LEN_2_data;
  const double *gkmprob_data;
  const double *pwm_prob_data;
  const double *shift_data;
  double L;
  double d;
  double name_len;
  double pos_tmp;
  double y;
  double *C_data;
  double *LEN_data;
  double *b_mat_data;
  double *ind_data;
  double *mat_data;
  double *neg_data;
  double *pos_data;
  double *vec_data;
  int b_i;
  int b_loop_ub;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int j;
  int k;
  int loop_ub;
  int nx;
  int pwm_prob_tmp;
  int *iidx_data;
  bool *b_pos_data;
  names_data = names->data;
  gkmprob_data = gkmprob->data;
  shift_data = shift->data;
  LEN_2_data = LEN_2->data;
  LEN_data = LEN->data;
  pwm_prob_data = pwm_prob->data;
  emxInit_real_T(&mat, 2);

  /* 'mapTF:202' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  /* 'mapTF:203' L = length(ss)-l_svm+1; */
  L = ((double)ss->size[1] - l_svm) + 1.0;

  /* 'mapTF:204' ACGT = 'ACGT'; */
  /* 'mapTF:205' lab = [lab; 0]; */
  /* 'mapTF:206' n = sum(LEN)+1; */
  y = blockedSummation(LEN, LEN->size[0]);

  /* 'mapTF:207' mat = zeros(n,L)-Inf; */
  i = mat->size[0] * mat->size[1];
  mat->size[0] = (int)(y + 1.0);
  loop_ub = (int)L;
  mat->size[1] = (int)L;
  emxEnsureCapacity_real_T(mat, i);
  mat_data = mat->data;
  nx = (int)(y + 1.0) * (int)L;
  for (i = 0; i < nx; i++) {
    mat_data[i] = rtMinusInf;
  }

  emxInit_real_T(&ind, 2);

  /* 'mapTF:208' ind = zeros(n,L); */
  i = ind->size[0] * ind->size[1];
  ind->size[0] = (int)(y + 1.0);
  ind->size[1] = (int)L;
  emxEnsureCapacity_real_T(ind, i);
  ind_data = ind->data;
  for (i = 0; i < nx; i++) {
    ind_data[i] = 0.0;
  }

  emxInit_real_T(&b_mat, 1);

  /* 'mapTF:209' LEN = [0;LEN]; */
  i = b_mat->size[0];
  b_mat->size[0] = LEN->size[0] + 1;
  emxEnsureCapacity_real_T(b_mat, i);
  b_mat_data = b_mat->data;
  b_mat_data[0] = 0.0;
  b_loop_ub = LEN->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    b_mat_data[i + 1] = LEN_data[i];
  }

  i = LEN->size[0];
  LEN->size[0] = b_mat->size[0];
  emxEnsureCapacity_real_T(LEN, i);
  LEN_data = LEN->data;
  b_loop_ub = b_mat->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    LEN_data[i] = b_mat_data[i];
  }

  emxInit_real_T(&C, 1);

  /* 'mapTF:210' C = cumsum(LEN); */
  i = C->size[0];
  C->size[0] = LEN->size[0];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  b_loop_ub = LEN->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    C_data[i] = LEN_data[i];
  }

  if (LEN->size[0] != 1) {
    i = LEN->size[0];
    for (k = 0; k <= i - 2; k++) {
      C_data[k + 1] += C_data[k];
    }
  }

  emxInit_real_T(&pos, 1);

  /* 'mapTF:211' pos = log(gkmprob); */
  i = pos->size[0];
  pos->size[0] = gkmprob->size[0];
  emxEnsureCapacity_real_T(pos, i);
  pos_data = pos->data;
  b_loop_ub = gkmprob->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    pos_data[i] = gkmprob_data[i];
  }

  nx = gkmprob->size[0];
  for (k = 0; k < nx; k++) {
    pos_data[k] = log(pos_data[k]);
  }

  emxInit_real_T(&neg, 1);

  /* 'mapTF:212' neg = log(1-gkmprob); */
  i = neg->size[0];
  neg->size[0] = gkmprob->size[0];
  emxEnsureCapacity_real_T(neg, i);
  neg_data = neg->data;
  b_loop_ub = gkmprob->size[0];
  nx = (gkmprob->size[0] / 2) << 1;
  k = nx - 2;
  for (i = 0; i <= k; i += 2) {
    _mm_storeu_pd(&neg_data[i], _mm_sub_pd(_mm_set1_pd(1.0), _mm_loadu_pd
      (&gkmprob_data[i])));
  }

  for (i = nx; i < b_loop_ub; i++) {
    neg_data[i] = 1.0 - gkmprob_data[i];
  }

  nx = neg->size[0];
  for (k = 0; k < nx; k++) {
    neg_data[k] = log(neg_data[k]);
  }

  /* 'mapTF:213' for i = 1:a */
  i = (int)a;
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:214' mat(C(i)+1,1) = pwm_prob(C(i)+1,1); */
    nx = (int)(C_data[b_i] + 1.0) - 1;
    mat_data[nx] = pwm_prob_data[nx];
  }

  /* 'mapTF:216' mat(n,1) = neg(1); */
  mat_data[(int)(y + 1.0) - 1] = neg_data[0];

  /* 'mapTF:217' for i = 2:L */
  i1 = (int)(L + -1.0);
  emxInit_real_T(&vec, 2);
  vec_data = vec->data;
  for (b_i = 0; b_i < i1; b_i++) {
    /* 'mapTF:218' for j = 1:a */
    for (j = 0; j < i; j++) {
      /* 'mapTF:219' [mat(C(j)+1,i),ind(C(j)+1,i)] = max(mat(1:end,i-1)+pwm_prob(C(j)+1,i)+pos(i)); */
      if (1 > mat->size[0]) {
        b_loop_ub = 0;
      } else {
        b_loop_ub = mat->size[0];
      }

      pwm_prob_tmp = (int)(C_data[j] + 1.0) - 1;
      name_len = pwm_prob_data[pwm_prob_tmp + pwm_prob->size[0] * (b_i + 1)];
      pos_tmp = pos_data[b_i + 1];
      i2 = b_mat->size[0];
      b_mat->size[0] = b_loop_ub;
      emxEnsureCapacity_real_T(b_mat, i2);
      b_mat_data = b_mat->data;
      nx = (b_loop_ub / 2) << 1;
      k = nx - 2;
      for (i2 = 0; i2 <= k; i2 += 2) {
        r = _mm_loadu_pd(&mat_data[i2 + mat->size[0] * b_i]);
        _mm_storeu_pd(&b_mat_data[i2], _mm_add_pd(_mm_add_pd(r, _mm_set1_pd
          (name_len)), _mm_set1_pd(pos_tmp)));
      }

      for (i2 = nx; i2 < b_loop_ub; i2++) {
        b_mat_data[i2] = (mat_data[i2 + mat->size[0] * b_i] + name_len) +
          pos_tmp;
      }

      d_maximum(b_mat, &d, &nx);
      mat_data[((int)(C_data[j] + 1.0) + mat->size[0] * (b_i + 1)) - 1] = d;
      ind_data[pwm_prob_tmp + ind->size[0] * (b_i + 1)] = nx;

      /* 'mapTF:220' mat(C(j)+2:LEN(j+1),i) = mat(C(j)+1:(LEN(j+1)-1),i-1)+pwm_prob(C(j)+2:LEN(j+1),i)+pos(i); */
      d = LEN_data[j + 1];
      if (C_data[j] + 1.0 > d - 1.0) {
        pwm_prob_tmp = 0;
        i2 = 0;
      } else {
        i2 = (int)(d - 1.0);
      }

      if (C_data[j] + 2.0 > d) {
        i3 = 0;
        nx = 0;
        i4 = 0;
      } else {
        i3 = (int)(C_data[j] + 2.0) - 1;
        nx = (int)d;
        i4 = (int)(C_data[j] + 2.0) - 1;
      }

      b_loop_ub = i2 - pwm_prob_tmp;
      if (b_loop_ub == nx - i3) {
        i2 = b_mat->size[0];
        b_mat->size[0] = b_loop_ub;
        emxEnsureCapacity_real_T(b_mat, i2);
        b_mat_data = b_mat->data;
        nx = (b_loop_ub / 2) << 1;
        k = nx - 2;
        for (i2 = 0; i2 <= k; i2 += 2) {
          r = _mm_loadu_pd(&mat_data[(pwm_prob_tmp + i2) + mat->size[0] * b_i]);
          _mm_storeu_pd(&b_mat_data[i2], _mm_add_pd(_mm_add_pd(r, _mm_loadu_pd
            (&pwm_prob_data[(i3 + i2) + pwm_prob->size[0] * (b_i + 1)])),
            _mm_set1_pd(pos_tmp)));
        }

        for (i2 = nx; i2 < b_loop_ub; i2++) {
          b_mat_data[i2] = (mat_data[(pwm_prob_tmp + i2) + mat->size[0] * b_i] +
                            pwm_prob_data[(i3 + i2) + pwm_prob->size[0] * (b_i +
            1)]) + pos_tmp;
        }

        b_loop_ub = b_mat->size[0];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          mat_data[(i4 + i2) + mat->size[0] * (b_i + 1)] = b_mat_data[i2];
        }
      } else {
        d_binary_expand_op(mat, i4, b_i, pwm_prob_tmp, i2 - 1, pwm_prob, i3, nx
                           - 1, pos);
        mat_data = mat->data;
      }

      /* 'mapTF:221' ind(C(j)+2:LEN(j+1),i) = C(j)+1:(LEN(j+1)-1); */
      if (C_data[j] + 2.0 > d) {
        i2 = 0;
        i3 = 0;
      } else {
        i2 = (int)(C_data[j] + 2.0) - 1;
        i3 = (int)LEN_data[j + 1];
      }

      if (rtIsNaN(C_data[j] + 1.0) || rtIsNaN(d - 1.0)) {
        nx = vec->size[0] * vec->size[1];
        vec->size[0] = 1;
        vec->size[1] = 1;
        emxEnsureCapacity_real_T(vec, nx);
        vec_data = vec->data;
        vec_data[0] = rtNaN;
      } else if (d - 1.0 < C_data[j] + 1.0) {
        vec->size[0] = 1;
        vec->size[1] = 0;
      } else if ((rtIsInf(C_data[j] + 1.0) || rtIsInf(d - 1.0)) && (C_data[j] +
                  1.0 == d - 1.0)) {
        nx = vec->size[0] * vec->size[1];
        vec->size[0] = 1;
        vec->size[1] = 1;
        emxEnsureCapacity_real_T(vec, nx);
        vec_data = vec->data;
        vec_data[0] = rtNaN;
      } else if (floor(C_data[j] + 1.0) == C_data[j] + 1.0) {
        nx = vec->size[0] * vec->size[1];
        vec->size[0] = 1;
        b_loop_ub = (int)floor((d - 1.0) - (C_data[j] + 1.0));
        vec->size[1] = b_loop_ub + 1;
        emxEnsureCapacity_real_T(vec, nx);
        vec_data = vec->data;
        for (nx = 0; nx <= b_loop_ub; nx++) {
          vec_data[nx] = (C_data[j] + 1.0) + (double)nx;
        }
      } else {
        eml_float_colon(C_data[j] + 1.0, d - 1.0, vec);
        vec_data = vec->data;
      }

      b_loop_ub = i3 - i2;
      for (i3 = 0; i3 < b_loop_ub; i3++) {
        ind_data[(i2 + i3) + ind->size[0] * (b_i + 1)] = vec_data[i3];
      }
    }

    /* 'mapTF:223' [mat(n,i),ind(n,i)] = max(mat(1:end,i-1)+neg(i)); */
    if (1 > mat->size[0]) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = mat->size[0];
    }

    name_len = neg_data[b_i + 1];
    i2 = b_mat->size[0];
    b_mat->size[0] = b_loop_ub;
    emxEnsureCapacity_real_T(b_mat, i2);
    b_mat_data = b_mat->data;
    nx = (b_loop_ub / 2) << 1;
    k = nx - 2;
    for (i2 = 0; i2 <= k; i2 += 2) {
      r = _mm_loadu_pd(&mat_data[i2 + mat->size[0] * b_i]);
      _mm_storeu_pd(&b_mat_data[i2], _mm_add_pd(r, _mm_set1_pd(name_len)));
    }

    for (i2 = nx; i2 < b_loop_ub; i2++) {
      b_mat_data[i2] = mat_data[i2 + mat->size[0] * b_i] + name_len;
    }

    d_maximum(b_mat, &d, &nx);
    mat_data[((int)(y + 1.0) + mat->size[0] * (b_i + 1)) - 1] = d;
    ind_data[((int)(y + 1.0) + ind->size[0] * (b_i + 1)) - 1] = nx;
  }

  emxFree_real_T(&mat);

  /* 'mapTF:225' path = zeros(L,1); */
  i = pos->size[0];
  pos->size[0] = (int)L;
  emxEnsureCapacity_real_T(pos, i);
  pos_data = pos->data;
  for (i = 0; i < loop_ub; i++) {
    pos_data[i] = 0.0;
  }

  /* 'mapTF:226' path(end) = n; */
  pos_data[(int)L - 1] = y + 1.0;

  /* 'mapTF:228' for i = L-1:-1:1 */
  i = (int)-((-1.0 - (L - 1.0)) + 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    name_len = (L - 1.0) + -(double)b_i;

    /* 'mapTF:229' path(i) = ind(path(i+1),i+1); */
    pos_data[(int)name_len - 1] = ind_data[((int)pos_data[(int)(unsigned int)
      name_len] + ind->size[0] * (int)(unsigned int)name_len) - 1];
  }

  emxFree_real_T(&ind);
  emxInit_real_T(&L2, 2);
  mat_data = L2->data;

  /* 'mapTF:232' total_len = length(ss); */
  /* 'mapTF:233' for i = 1:a */
  /* 'mapTF:246' PATH = zeros(total_len,1); */
  /* 'mapTF:247' C_2 = cumsum([0;LEN_2]); */
  /* 'mapTF:248' for i = 1:a */
  /* 'mapTF:257' N =sum(LEN_2)+1; */
  /* 'mapTF:258' PATH(PATH==0) = N; */
  /* 'mapTF:259' L2 = []; */
  L2->size[0] = 0;
  L2->size[1] = 4;

  /* 'mapTF:260' for i = 1:length(PWM2) */
  i = PWM2->size[0];
  emxInit_int32_T(&iidx, 1);
  emxInit_boolean_T(&b_pos, 1);
  emxInit_real_T(&b_L2, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:261' f = find(path==C(i)+1); */
    loop_ub = pos->size[0];
    i1 = b_pos->size[0];
    b_pos->size[0] = pos->size[0];
    emxEnsureCapacity_boolean_T(b_pos, i1);
    b_pos_data = b_pos->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_pos_data[i1] = (pos_data[i1] == C_data[b_i] + 1.0);
    }

    eml_find(b_pos, iidx);
    iidx_data = iidx->data;
    i1 = neg->size[0];
    neg->size[0] = iidx->size[0];
    emxEnsureCapacity_real_T(neg, i1);
    neg_data = neg->data;
    loop_ub = iidx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      neg_data[i1] = iidx_data[i1];
    }

    /* 'mapTF:262' if ~isempty(f) */
    if (neg->size[0] != 0) {
      /* 'mapTF:263' for j = 1:length(f) */
      i1 = neg->size[0];
      for (j = 0; j < i1; j++) {
        /* 'mapTF:264' if f(j)+shift(i)+LEN_2(i)< length(ss)-l_svm+1 && f(j)+shift(i) > l_svm-1 */
        d = neg_data[j] + shift_data[b_i];
        pos_tmp = d + LEN_2_data[b_i];
        if ((pos_tmp < ((double)ss->size[1] - l_svm) + 1.0) && (d > l_svm - 1.0))
        {
          /* 'mapTF:265' [~,ind] = max(gkmprob(f(j):f(j)+2*shift(i)+LEN_2(i)-l_svm)); */
          name_len = ((neg_data[j] + 2.0 * shift_data[b_i]) + LEN_2_data[b_i]) -
            l_svm;
          if (neg_data[j] > name_len) {
            i2 = 0;
            i3 = 0;
          } else {
            i2 = (int)neg_data[j] - 1;
            i3 = (int)name_len;
          }

          loop_ub = i3 - i2;
          i3 = b_mat->size[0];
          b_mat->size[0] = loop_ub;
          emxEnsureCapacity_real_T(b_mat, i3);
          b_mat_data = b_mat->data;
          for (i3 = 0; i3 < loop_ub; i3++) {
            b_mat_data[i3] = gkmprob_data[i2 + i3];
          }

          d_maximum(b_mat, &name_len, &nx);

          /* 'mapTF:266' if f(j)-1+ind < 3 */
          name_len = (neg_data[j] - 1.0) + (double)nx;
          if (name_len < 3.0) {
            /* 'mapTF:267' R = 1:5; */
            i2 = vec->size[0] * vec->size[1];
            vec->size[0] = 1;
            vec->size[1] = 5;
            emxEnsureCapacity_real_T(vec, i2);
            vec_data = vec->data;
            for (i2 = 0; i2 < 5; i2++) {
              vec_data[i2] = (double)i2 + 1.0;
            }
          } else if (name_len > L - 2.0) {
            /* 'mapTF:268' elseif f(j)-1+ind > L-2 */
            /* 'mapTF:269' R = L-4:L; */
            if (L < L - 4.0) {
              vec->size[0] = 1;
              vec->size[1] = 0;
            } else if ((rtIsInf(L - 4.0) || rtIsInf(L)) && (L - 4.0 == L)) {
              i2 = vec->size[0] * vec->size[1];
              vec->size[0] = 1;
              vec->size[1] = 1;
              emxEnsureCapacity_real_T(vec, i2);
              vec_data = vec->data;
              vec_data[0] = rtNaN;
            } else if (L - 4.0 == L - 4.0) {
              i2 = vec->size[0] * vec->size[1];
              vec->size[0] = 1;
              loop_ub = (int)floor(L - (L - 4.0));
              vec->size[1] = loop_ub + 1;
              emxEnsureCapacity_real_T(vec, i2);
              vec_data = vec->data;
              for (i2 = 0; i2 <= loop_ub; i2++) {
                vec_data[i2] = (L - 4.0) + (double)i2;
              }
            } else {
              eml_float_colon(L - 4.0, L, vec);
              vec_data = vec->data;
            }
          } else {
            /* 'mapTF:270' else */
            /* 'mapTF:271' R = f(j)+ind-3:f(j)+ind+1; */
            name_len = neg_data[j] + (double)nx;
            if (name_len + 1.0 < name_len - 3.0) {
              vec->size[0] = 1;
              vec->size[1] = 0;
            } else if (name_len - 3.0 == name_len - 3.0) {
              i2 = vec->size[0] * vec->size[1];
              vec->size[0] = 1;
              loop_ub = (int)((name_len + 1.0) - (name_len - 3.0));
              vec->size[1] = loop_ub + 1;
              emxEnsureCapacity_real_T(vec, i2);
              vec_data = vec->data;
              for (i2 = 0; i2 <= loop_ub; i2++) {
                vec_data[i2] = (name_len - 3.0) + (double)i2;
              }
            } else {
              eml_float_colon(name_len - 3.0, name_len + 1.0, vec);
              vec_data = vec->data;
            }
          }

          /* 'mapTF:273' L2 = [L2; f(j)+shift(i) f(j)+shift(i)+LEN_2(i)-1 mean(gkmprob(R)) i]; */
          i2 = b_mat->size[0];
          b_mat->size[0] = vec->size[1];
          emxEnsureCapacity_real_T(b_mat, i2);
          b_mat_data = b_mat->data;
          loop_ub = vec->size[1];
          for (i2 = 0; i2 < loop_ub; i2++) {
            b_mat_data[i2] = gkmprob_data[(int)vec_data[i2] - 1];
          }

          loop_ub = L2->size[0];
          i2 = b_L2->size[0] * b_L2->size[1];
          b_L2->size[0] = L2->size[0] + 1;
          b_L2->size[1] = 4;
          emxEnsureCapacity_real_T(b_L2, i2);
          LEN_data = b_L2->data;
          for (i2 = 0; i2 < 4; i2++) {
            for (i3 = 0; i3 < loop_ub; i3++) {
              LEN_data[i3 + b_L2->size[0] * i2] = mat_data[i3 + L2->size[0] * i2];
            }
          }

          LEN_data[L2->size[0]] = d;
          LEN_data[L2->size[0] + b_L2->size[0]] = pos_tmp - 1.0;
          LEN_data[L2->size[0] + b_L2->size[0] * 2] = blockedSummation(b_mat,
            vec->size[1]) / (double)vec->size[1];
          LEN_data[L2->size[0] + b_L2->size[0] * 3] = (double)b_i + 1.0;
          i2 = L2->size[0] * L2->size[1];
          L2->size[0] = b_L2->size[0];
          L2->size[1] = 4;
          emxEnsureCapacity_real_T(L2, i2);
          mat_data = L2->data;
          loop_ub = b_L2->size[0] * 4;
          for (i2 = 0; i2 < loop_ub; i2++) {
            mat_data[i2] = LEN_data[i2];
          }
        }
      }
    }
  }

  emxFree_boolean_T(&b_pos);
  emxFree_real_T(&b_mat);
  emxFree_real_T(&vec);
  emxFree_real_T(&neg);
  emxFree_real_T(&C);

  /* 'mapTF:278' name_len = numel(L2)/4; */
  name_len = (double)(L2->size[0] << 2) / 4.0;

  /* 'mapTF:279' NAME = cell(name_len,1); */
  /* 'mapTF:280' Lmat = zeros(name_len,4); */
  i = (int)name_len;
  i1 = Lmat->size[0] * Lmat->size[1];
  Lmat->size[0] = (int)name_len;
  Lmat->size[1] = 4;
  emxEnsureCapacity_real_T(Lmat, i1);
  LEN_data = Lmat->data;
  loop_ub = (int)name_len << 2;
  for (i1 = 0; i1 < loop_ub; i1++) {
    LEN_data[i1] = 0.0;
  }

  /* 'mapTF:281' for i=1:name_len */
  i1 = NAME->size[0];
  NAME->size[0] = (int)name_len;
  emxEnsureCapacity_cell_wrap_3(NAME, i1);
  NAME_data = NAME->data;
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:282' NAME{i} = ''; */
    NAME_data[b_i].f1->size[0] = 1;
    NAME_data[b_i].f1->size[1] = 0;
  }

  /* 'mapTF:284' if ~isempty(L2) */
  if (L2->size[0] != 0) {
    /* 'mapTF:285' [~,b] = sort(L2(:,1)); */
    loop_ub = L2->size[0];
    i1 = pos->size[0];
    pos->size[0] = L2->size[0];
    emxEnsureCapacity_real_T(pos, i1);
    pos_data = pos->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      pos_data[i1] = mat_data[i1];
    }

    sort(pos, iidx);
    iidx_data = iidx->data;

    /* 'mapTF:286' L2 = L2(b,:); */
    i1 = b_L2->size[0] * b_L2->size[1];
    b_L2->size[0] = iidx->size[0];
    b_L2->size[1] = 4;
    emxEnsureCapacity_real_T(b_L2, i1);
    LEN_data = b_L2->data;
    loop_ub = iidx->size[0];
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        LEN_data[i2 + b_L2->size[0] * i1] = mat_data[(iidx_data[i2] + L2->size[0]
          * i1) - 1];
      }
    }

    i1 = L2->size[0] * L2->size[1];
    L2->size[0] = b_L2->size[0];
    L2->size[1] = 4;
    emxEnsureCapacity_real_T(L2, i1);
    mat_data = L2->data;
    loop_ub = b_L2->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      mat_data[i1] = LEN_data[i1];
    }

    /* 'mapTF:287' for i = 1:name_len */
    i1 = NAME->size[0];
    NAME->size[0] = (int)name_len;
    emxEnsureCapacity_cell_wrap_3(NAME, i1);
    NAME_data = NAME->data;
    i1 = Lmat->size[0] * Lmat->size[1];
    Lmat->size[0] = (int)name_len;
    Lmat->size[1] = 4;
    emxEnsureCapacity_real_T(Lmat, i1);
    LEN_data = Lmat->data;
    for (b_i = 0; b_i < i; b_i++) {
      /* 'mapTF:288' NAME{i} = names{L2(i,4)}; */
      i1 = NAME_data[b_i].f1->size[0] * NAME_data[b_i].f1->size[1];
      NAME_data[b_i].f1->size[0] = 1;
      NAME_data[b_i].f1->size[1] = names_data[(int)mat_data[b_i + L2->size[0] *
        3] - 1].f1->size[1];
      emxEnsureCapacity_char_T(NAME_data[b_i].f1, i1);
      loop_ub = names_data[(int)mat_data[b_i + L2->size[0] * 3] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        NAME_data[b_i].f1->data[i1] = names_data[(int)mat_data[b_i + L2->size[0]
          * 3] - 1].f1->data[i1];
      }

      /* 'mapTF:289' Lmat(i,:) = [L2(i,4) L2(i,1) L2(i,2) L2(i,3)]; */
      LEN_data[b_i] = mat_data[b_i + L2->size[0] * 3];
      LEN_data[b_i + Lmat->size[0]] = mat_data[b_i];
      LEN_data[b_i + Lmat->size[0] * 2] = mat_data[b_i + L2->size[0]];
      LEN_data[b_i + Lmat->size[0] * 3] = mat_data[b_i + L2->size[0] * 2];
    }
  }

  emxFree_real_T(&b_L2);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&L2);
  emxFree_real_T(&pos);
}

/*
 * function PWM_corr(fn, VV, NN, LL, seq)
 */
static void PWM_corr(const emxArray_char_T *fn, const emxArray_cell_wrap_0 *VV,
                     const emxArray_cell_wrap_6 *NN, const emxArray_cell_wrap_4 *
                     LL, const emxArray_cell_wrap_3 *seq)
{
  static const char b_cv[19] = { '_', 'T', 'F', 'B', 'S', '_', 'l', 'o', 'c',
    'a', 't', 'i', 'o', 'n', 's', '.', 'o', 'u', 't' };

  FILE* b_NULL;
  FILE* filestar;
  const cell_wrap_0 *VV_data;
  const cell_wrap_3 *seq_data;
  const cell_wrap_4 *LL_data;
  const cell_wrap_6 *NN_data;
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

  /* 'mapTF:392' n = length(VV); */
  /* 'mapTF:393' fid1 = fopen([fn '_TFBS_locations.out'],'w'); */
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

  /* 'mapTF:394' for j = 1:n */
  i = VV->size[0];
  emxInit_char_T(&varargin_2, 2);
  emxInit_char_T(&varargin_8, 2);
  for (j = 0; j < i; j++) {
    /* 'mapTF:395' NAME = NN{j}; */
    /* 'mapTF:396' Lmat = LL{j}; */
    /* 'mapTF:397' varscore = VV{j}; */
    /* 'mapTF:398' s = seq{j*2}; */
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

    /* 'mapTF:399' L = length(NAME); */
    /* 'mapTF:400' for i = 1:L */
    i1 = NN_data[j].f1->size[0];
    for (b_i = 0; b_i < i1; b_i++) {
      /* 'mapTF:401' if varscore(i) > 0.6 */
      if (VV_data[j].f1->data[b_i] > 0.6) {
        /* 'mapTF:402' fprintf(fid1, '%d\t%s\t%d\t%d\t%d\t%f\t%f\t%s\n', int32(j), NAME{i}, ... */
        /* 'mapTF:403'                 int32(Lmat(i,1)), int32(Lmat(i,2)), ... */
        /* 'mapTF:404'                 int32(Lmat(i,3)), Lmat(i,4), varscore(i), s(Lmat(i,2):Lmat(i,3))); */
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

  /* 'mapTF:408' fclose(fid1); */
  cfclose(fileid);
}

static void b_binary_expand_op(emxArray_real_T *f, const emxArray_real_T
  *all_pwm, const emxArray_int8_T *ss_onehot, int i2, int i3, int i4)
{
  emxArray_real_T *b_all_pwm;
  const double *all_pwm_data;
  double *b_all_pwm_data;
  int all_pwm_tmp;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  int unnamed_idx_1;
  const signed char *ss_onehot_data;
  ss_onehot_data = ss_onehot->data;
  all_pwm_data = all_pwm->data;
  emxInit_real_T(&b_all_pwm, 3);
  unnamed_idx_1 = i3 - i4;
  i = b_all_pwm->size[0] * b_all_pwm->size[1] * b_all_pwm->size[2];
  b_all_pwm->size[0] = 4;
  if (unnamed_idx_1 == 1) {
    b_all_pwm->size[1] = all_pwm->size[1];
  } else {
    b_all_pwm->size[1] = unnamed_idx_1;
  }

  b_all_pwm->size[2] = all_pwm->size[2];
  emxEnsureCapacity_real_T(b_all_pwm, i);
  b_all_pwm_data = b_all_pwm->data;
  stride_0_1 = (all_pwm->size[1] != 1);
  stride_1_1 = (unnamed_idx_1 != 1);
  loop_ub = all_pwm->size[2];
  for (i = 0; i < loop_ub; i++) {
    aux_0_1 = 0;
    aux_1_1 = 0;
    if (unnamed_idx_1 == 1) {
      b_loop_ub = all_pwm->size[1];
    } else {
      b_loop_ub = unnamed_idx_1;
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      all_pwm_tmp = i2 + aux_1_1;
      b_all_pwm_data[4 * i1 + 4 * b_all_pwm->size[1] * i] = all_pwm_data[4 *
        aux_0_1 + 4 * all_pwm->size[1] * i] * (double)ss_onehot_data[4 *
        all_pwm_tmp];
      b_all_pwm_data[(4 * i1 + 4 * b_all_pwm->size[1] * i) + 1] = all_pwm_data
        [(4 * aux_0_1 + 4 * all_pwm->size[1] * i) + 1] * (double)ss_onehot_data
        [4 * all_pwm_tmp + 1];
      b_all_pwm_data[(4 * i1 + 4 * b_all_pwm->size[1] * i) + 2] = all_pwm_data
        [(4 * aux_0_1 + 4 * all_pwm->size[1] * i) + 2] * (double)ss_onehot_data
        [4 * all_pwm_tmp + 2];
      b_all_pwm_data[(4 * i1 + 4 * b_all_pwm->size[1] * i) + 3] = all_pwm_data
        [(4 * aux_0_1 + 4 * all_pwm->size[1] * i) + 3] * (double)ss_onehot_data
        [4 * all_pwm_tmp + 3];
      aux_1_1 += stride_1_1;
      aux_0_1 += stride_0_1;
    }
  }

  nonzeros(b_all_pwm, f);
  emxFree_real_T(&b_all_pwm);
}

/*
 * function GC = countGC(s)
 */
static double countGC(const emxArray_cell_wrap_2 *s)
{
  const cell_wrap_2 *s_data;
  double GC;
  double L;
  double d;
  int b_i;
  int i;
  int i1;
  int j;
  s_data = s->data;

  /* 'mapTF:560' GC = 0; */
  GC = 0.0;

  /* 'mapTF:561' L = 0; */
  L = 0.0;

  /* 'mapTF:562' for i = 1:length(s) */
  i = s->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:563' S = s{i}; */
    /* 'mapTF:564' N = length(S); */
    /* 'mapTF:565' L = N + L; */
    i1 = s_data[b_i].f1->size[1];
    L += (double)s_data[b_i].f1->size[1];

    /* 'mapTF:566' for j = 1:N */
    for (j = 0; j < i1; j++) {
      /* 'mapTF:567' if S(j) == 2 || S(j) == 3 */
      d = s_data[b_i].f1->data[j];
      if ((d == 2.0) || (d == 3.0)) {
        /* 'mapTF:568' GC = GC + 1; */
        GC++;
      }
    }
  }

  /* 'mapTF:572' GC = GC/L; */
  GC /= L;
  return GC;
}

static void d_binary_expand_op(emxArray_real_T *mat, int i1, int i, int i2, int
  i3, const emxArray_real_T *pwm_prob, int i4, int i5, const emxArray_real_T
  *pos)
{
  emxArray_real_T *b_mat;
  const double *pos_data;
  const double *pwm_prob_data;
  double b_pos;
  double *b_mat_data;
  double *mat_data;
  int b_i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  pos_data = pos->data;
  pwm_prob_data = pwm_prob->data;
  mat_data = mat->data;
  emxInit_real_T(&b_mat, 1);
  b_pos = pos_data[i + 1];
  b_i = b_mat->size[0];
  if ((i5 - i4) + 1 == 1) {
    b_mat->size[0] = (i3 - i2) + 1;
  } else {
    b_mat->size[0] = (i5 - i4) + 1;
  }

  emxEnsureCapacity_real_T(b_mat, b_i);
  b_mat_data = b_mat->data;
  stride_0_0 = ((i3 - i2) + 1 != 1);
  stride_1_0 = ((i5 - i4) + 1 != 1);
  if ((i5 - i4) + 1 == 1) {
    loop_ub = (i3 - i2) + 1;
  } else {
    loop_ub = (i5 - i4) + 1;
  }

  for (b_i = 0; b_i < loop_ub; b_i++) {
    b_mat_data[b_i] = (mat_data[(i2 + b_i * stride_0_0) + mat->size[0] * i] +
                       pwm_prob_data[(i4 + b_i * stride_1_0) + pwm_prob->size[0]
                       * (i + 1)]) + b_pos;
  }

  loop_ub = b_mat->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    mat_data[(i1 + b_i) + mat->size[0] * (i + 1)] = b_mat_data[b_i];
  }

  emxFree_real_T(&b_mat);
}

/*
 * function [mat, names, len] = getMOTIF(fn)
 */
static void getMOTIF(const emxArray_char_T *fn, emxArray_cell_wrap_4 *mat,
                     emxArray_cell_wrap_3 *names, emxArray_real_T *len)
{
  static const char b[5] = { 'M', 'O', 'T', 'I', 'F' };

  static const char b_cv[3] = { 'a', 'l', 'l' };

  FILE* b_NULL;
  FILE* filestar;
  int st;
  int wherefrom;
  long position_t;
  cell_wrap_3 *names_data;
  cell_wrap_4 *mat_data;
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

  /* 'mapTF:343' fid = fopen(fn); */
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

  /* 'mapTF:344' if fid < 0 */
  if (fid < 0) {
    /* 'mapTF:345' fprintf('ERROR: Cannot open gkmPWM motif files\n'); */
    printf("ERROR: Cannot open gkmPWM motif files\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:347' curr_pos = ftell(fid); */
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

  /* 'mapTF:348' idx=0; */
  idx = 0.0;

  /* 'mapTF:349' while ~feof(fid) */
  b_NULL = NULL;
  emxInit_char_T(&line, 2);
  do {
    exitg1 = 0;
    getfilestar(fid, &filestar, &b_bool);
    if (filestar == b_NULL) {
      i = 0;
    } else {
      st = feof(filestar);
      i = ((int)st != 0);
    }

    if (i == 0) {
      /* 'mapTF:350' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF:351' if length(line) >= 5 && strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF:352' idx=idx+1; */
          idx++;
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /* 'mapTF:355' fseek(fid, curr_pos, 'bof'); */
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

  /* 'mapTF:357' mat = cell(idx, 1); */
  ret = (int)idx;
  i = mat->size[0];
  mat->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_4(mat, i);
  mat_data = mat->data;

  /* 'mapTF:358' mat = coder.nullcopy(mat); */
  /* 'mapTF:359' names = cell(idx, 1); */
  i = names->size[0];
  names->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_3(names, i);
  names_data = names->data;
  for (i = 0; i < ret; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 4;
    names_data[i].f1->size[0] = 1;
    names_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:360' names = coder.nullcopy(names); */
  /* 'mapTF:361' len = zeros(idx, 1); */
  i = len->size[0];
  len->size[0] = (int)idx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;
  for (i = 0; i < ret; i++) {
    len_data[i] = 0.0;
  }

  /* 'mapTF:363' i=0; */
  b_i = 0U;

  /* 'mapTF:364' while ~feof(fid) */
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
      i = 0;
    } else {
      st = feof(filestar);
      i = ((int)st != 0);
    }

    if (i == 0) {
      /* 'mapTF:365' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF:366' if length(line) < 5 || ~strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF:369' i = i + 1; */
          b_i++;

          /* 'mapTF:370' names{i} = fgetl(fid); */
          fgetl(fid, names_data[(int)b_i - 1].f1);

          /* 'mapTF:371' len(i) = real(str2double(fgetl(fid))); */
          fgetl(fid, line);
          dc = str2double(line);
          len_data[(int)b_i - 1] = dc.re;

          /* 'mapTF:372' tmp = zeros(len(i), 4); */
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

          /* 'mapTF:373' line_count = 1; */
          line_count = 1U;

          /* 'mapTF:374' line = fgetl(fid); */
          fgetl(fid, line);

          /* 'mapTF:375' while ~isempty(line) */
          while (line->size[1] != 0) {
            /* 'mapTF:376' [v1, remain] = strtok(line); */
            b_strtok(line, v1, remain);

            /* 'mapTF:377' [v2, remain] = strtok(remain); */
            b_strtok(remain, v2, b_remain);

            /* 'mapTF:378' [v3, v4] = strtok(remain); */
            b_strtok(b_remain, line, remain);

            /* 'mapTF:379' vals = [real(str2double(v1)),  */
            /* 'mapTF:380'                 real(str2double(v2)),  */
            /* 'mapTF:381'                 real(str2double(v3)),  */
            /* 'mapTF:382'                 real(str2double(v4))]'; */
            dc = str2double(v1);
            dc1 = str2double(v2);
            dc2 = str2double(line);
            dc3 = str2double(remain);
            tmp_data[(int)line_count - 1] = dc.re;
            tmp_data[((int)line_count + tmp->size[0]) - 1] = dc1.re;
            tmp_data[((int)line_count + tmp->size[0] * 2) - 1] = dc2.re;
            tmp_data[((int)line_count + tmp->size[0] * 3) - 1] = dc3.re;

            /* 'mapTF:383' tmp(line_count,:) = vals; */
            /* 'mapTF:384' line = fgetl(fid); */
            fgetl(fid, line);

            /* 'mapTF:385' line_count = line_count + 1; */
            line_count++;
          }

          /* 'mapTF:387' mat{i} = tmp; */
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

  /* 'mapTF:389' fclose(fid); */
  cfclose(fid);
}

/*
 * function [mat,w] = getdenovomotif(filename)
 */
static void getdenovomotif(const emxArray_char_T *filename, emxArray_cell_wrap_5
  *mat, emxArray_real_T *w)
{
  static const char cv1[6] = { 'w', 'e', 'i', 'g', 'h', 't' };

  static const char b[5] = { 'M', 'O', 'T', 'I', 'F' };

  static const char b_cv[3] = { 'a', 'l', 'l' };

  FILE* b_NULL;
  FILE* filestar;
  long position_t;
  cell_wrap_5 *mat_data;
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
  unsigned int b_i;
  int c_i;
  int exitg1;
  int fid;
  int i;
  int j;
  unsigned int line_count;
  int ret;
  int st;
  int wherefrom;
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
  /* 'mapTF:853' fid = fopen(filename); */
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

  /* 'mapTF:854' if fid < 0 */
  if (fid < 0) {
    /* 'mapTF:855' fprintf('ERROR: Cannot open gkmPWM motif files\n'); */
    printf("ERROR: Cannot open gkmPWM motif files\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:857' curr_pos = ftell(fid); */
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

  /* 'mapTF:858' idx=0; */
  idx = 0.0;

  /* 'mapTF:859' while ~feof(fid) */
  b_NULL = NULL;
  emxInit_char_T(&line, 2);
  do {
    exitg1 = 0;
    getfilestar(fid, &filestar, &b_bool);
    if (filestar == b_NULL) {
      i = 0;
    } else {
      st = feof(filestar);
      i = ((int)st != 0);
    }

    if (i == 0) {
      /* 'mapTF:860' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF:861' if length(line) >= 5 && strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF:862' idx=idx+1; */
          idx++;
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /* 'mapTF:865' fseek(fid, curr_pos, 'bof'); */
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

  /* 'mapTF:867' mat = cell(idx, 1); */
  ret = (int)idx;
  i = mat->size[0];
  mat->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_5(mat, i);
  mat_data = mat->data;
  for (i = 0; i < ret; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:868' mat = coder.nullcopy(mat); */
  /* 'mapTF:869' w = zeros(idx, 1); */
  i = w->size[0];
  w->size[0] = (int)idx;
  emxEnsureCapacity_real_T(w, i);
  w_data = w->data;
  for (i = 0; i < ret; i++) {
    w_data[i] = 0.0;
  }

  /* 'mapTF:871' i = 0; */
  b_i = 0U;

  /* 'mapTF:872' while ~feof(fid) */
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
      i = 0;
    } else {
      st = feof(filestar);
      i = ((int)st != 0);
    }

    if (i == 0) {
      /* 'mapTF:873' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF:874' if length(line) < 5 || ~strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF:877' i = i + 1; */
          b_i++;

          /* 'mapTF:878' line = fgetl(fid); */
          fgetl(fid, line);
          line_data = line->data;

          /* 'mapTF:880' [~, remain] = strtok(line); */
          b_strtok(line, a__9, remain);

          /* 'mapTF:881' [curr_w, remain] = strtok(remain); */
          b_strtok(remain, curr_w, b_remain);

          /* 'mapTF:882' [~, remain] = strtok(remain); */
          b_strtok(b_remain, a__9, remain);

          /* 'mapTF:883' [curr_alphabet, remain] = strtok(remain); */
          b_strtok(remain, curr_alphabet, b_remain);

          /* 'mapTF:884' [~, remain] = strtok(remain); */
          b_strtok(b_remain, a__9, remain);

          /* 'mapTF:885' [curr_length, ~] = strtok(remain); */
          b_strtok(remain, a__9, a__12);

          /* 'mapTF:886' w(i) = real(str2double(curr_w)); */
          dc = str2double(curr_w);
          w_data[(int)b_i - 1] = dc.re;

          /* 'mapTF:887' curr_alphabet = real(str2double(curr_alphabet)); */
          dc = str2double(curr_alphabet);

          /* 'mapTF:888' curr_length = real(str2double(curr_length)); */
          dc1 = str2double(a__9);

          /* 'mapTF:889' tmp = zeros(curr_length, curr_alphabet); */
          i = tmp->size[0] * tmp->size[1];
          tmp->size[0] = (int)dc1.re;
          tmp->size[1] = (int)dc.re;
          emxEnsureCapacity_real_T(tmp, i);
          tmp_data = tmp->data;
          ret = (int)dc1.re * (int)dc.re;
          for (i = 0; i < ret; i++) {
            tmp_data[i] = 0.0;
          }

          /* 'mapTF:890' line_count = 1; */
          line_count = 1U;

          /* 'mapTF:891' while ~isempty(line) */
          while (line->size[1] != 0) {
            /* 'mapTF:892' if strfind(line, 'weight') */
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

            if (match_out->size[1] != 0) {
              /* 'mapTF:893' line = fgetl(fid); */
              fgetl(fid, line);
            }

            /* 'mapTF:895' [v1, remain] = strtok(line); */
            b_strtok(line, curr_w, remain);

            /* 'mapTF:896' [v2, remain] = strtok(remain); */
            b_strtok(remain, curr_alphabet, b_remain);

            /* 'mapTF:897' [v3, v4] = strtok(remain); */
            b_strtok(b_remain, a__9, a__12);

            /* 'mapTF:898' vals = [real(str2double(v1)),  */
            /* 'mapTF:899'                 real(str2double(v2)),  */
            /* 'mapTF:900'                 real(str2double(v3)),  */
            /* 'mapTF:901'                 real(str2double(v4))]'; */
            dc = str2double(curr_w);
            dc1 = str2double(curr_alphabet);
            dc2 = str2double(a__9);
            dc3 = str2double(a__12);
            vals[0] = dc.re;
            vals[1] = dc1.re;
            vals[2] = dc2.re;
            vals[3] = dc3.re;

            /* 'mapTF:902' tmp(line_count,:) = vals; */
            ret = tmp->size[1];
            for (i = 0; i < ret; i++) {
              tmp_data[((int)line_count + tmp->size[0] * i) - 1] = vals[i];
            }

            /* 'mapTF:903' line = fgetl(fid); */
            fgetl(fid, line);
            line_data = line->data;

            /* 'mapTF:904' line_count = line_count + 1; */
            line_count++;
          }

          /* 'mapTF:906' mat{i} = tmp; */
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

  /* 'mapTF:908' fclose(fid); */
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

  /* 'mapTF:325' l = length(s); */
  /* 'mapTF:326' en = zeros(1,l); */
  i = en->size[0] * en->size[1];
  en->size[0] = 1;
  en->size[1] = s->size[1];
  emxEnsureCapacity_real_T(en, i);
  en_data = en->data;

  /* 'mapTF:327' for i = 1:l */
  i = s->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:328' if strcmp(s(i),'A') || strcmp(s(i), 'a') */
    c = s_data[b_i];
    if ((!(c != 'A')) || (!(c != 'a'))) {
      /* 'mapTF:329' en(i) = 0; */
      en_data[b_i] = 0.0;
    } else if ((!(c != 'C')) || (!(c != 'c'))) {
      /* 'mapTF:330' elseif strcmp(s(i),'C') || strcmp(s(i),'c') */
      /* 'mapTF:331' en(i) = 1; */
      en_data[b_i] = 1.0;
    } else if ((!(c != 'G')) || (!(c != 'g'))) {
      /* 'mapTF:332' elseif strcmp(s(i),'G') || strcmp(s(i),'g') */
      /* 'mapTF:333' en(i) = 2; */
      en_data[b_i] = 2.0;
    } else {
      /* 'mapTF:334' else */
      /* 'mapTF:335' en(i) = 3; */
      en_data[b_i] = 3.0;
    }
  }
}

/*
 * function [ind, M] = ppmsim(mot,lenvec)
 */
static void ppmsim(emxArray_cell_wrap_5 *mot, const emxArray_real_T *lenvec,
                   double *ind, double *M)
{
  __m128d r;
  cell_wrap_5 *mot_data;
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
  int loop_ub;
  int scalarLB;
  int vectorUB;
  mot_data = mot->data;

  /* 'mapTF:826' n = length(lenvec)-1; */
  /* 'mapTF:827' simmat = ones(n-1,1); */
  /* 'mapTF:828' for i = 1:n+1 */
  i = lenvec->size[0];
  emxInit_real_T(&diag_a, 2);
  emxInit_real_T(&A, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:829' mot{i} = mot{i}-1/4; */
    loop_ub = mot_data[b_i].f1->size[0] * mot_data[b_i].f1->size[1];
    scalarLB = (loop_ub / 2) << 1;
    vectorUB = scalarLB - 2;
    for (i1 = 0; i1 <= vectorUB; i1 += 2) {
      r = _mm_loadu_pd(&mot_data[b_i].f1->data[i1]);
      _mm_storeu_pd(&mot_data[b_i].f1->data[i1], _mm_sub_pd(r, _mm_set1_pd(0.25)));
    }

    for (i1 = scalarLB; i1 < loop_ub; i1++) {
      mot_data[b_i].f1->data[i1] -= 0.25;
    }

    /* 'mapTF:830' mot{i} = mot{i}/sqrt(sum(sum(mot{i}.^2))); */
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
    scalarLB = (loop_ub / 2) << 1;
    vectorUB = scalarLB - 2;
    for (i1 = 0; i1 <= vectorUB; i1 += 2) {
      r = _mm_loadu_pd(&mot_data[b_i].f1->data[i1]);
      _mm_storeu_pd(&mot_data[b_i].f1->data[i1], _mm_div_pd(r, _mm_set1_pd(MM)));
    }

    for (i1 = scalarLB; i1 < loop_ub; i1++) {
      mot_data[b_i].f1->data[i1] /= MM;
    }
  }

  /* 'mapTF:832' M = 0; */
  *M = 0.0;

  /* 'mapTF:833' ind = 1; */
  *ind = 1.0;

  /* 'mapTF:834' for j = 2:n+1 */
  i = lenvec->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&rmat, 2);
  emxInit_real_T(&diag_b, 2);
  emxInit_real_T(&diag_overall, 1);
  for (scalarLB = 0; scalarLB <= i - 2; scalarLB++) {
    /* 'mapTF:835' mat = mot{1}*mot{j}'; */
    if ((mot_data[0].f1->size[0] == 0) || (mot_data[0].f1->size[1] == 0) ||
        (mot_data[scalarLB + 1].f1->size[0] == 0) || (mot_data[scalarLB + 1].
         f1->size[1] == 0)) {
      i1 = mat->size[0] * mat->size[1];
      mat->size[0] = mot_data[0].f1->size[0];
      mat->size[1] = mot_data[scalarLB + 1].f1->size[0];
      emxEnsureCapacity_real_T(mat, i1);
      mat_data = mat->data;
      loop_ub = mot_data[0].f1->size[0] * mot_data[scalarLB + 1].f1->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        mat_data[i1] = 0.0;
      }
    } else {
      i1 = mat->size[0] * mat->size[1];
      mat->size[0] = mot_data[0].f1->size[0];
      mat->size[1] = mot_data[scalarLB + 1].f1->size[0];
      emxEnsureCapacity_real_T(mat, i1);
      mat_data = mat->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)mot_data[0].
                  f1->size[0], (blasint)mot_data[scalarLB + 1].f1->size[0],
                  (blasint)mot_data[0].f1->size[1], 1.0, &mot_data[0].f1->data[0],
                  (blasint)mot_data[0].f1->size[0], &mot_data[scalarLB + 1]
                  .f1->data[0], (blasint)mot_data[scalarLB + 1].f1->size[0], 0.0,
                  &mat_data[0], (blasint)mot_data[0].f1->size[0]);
    }

    /* 'mapTF:836' rmat = rot90(mot{1},2)*mot{j}'; */
    rot90(mot_data[0].f1, A);
    mat_data = A->data;
    if ((A->size[0] == 0) || (A->size[1] == 0) || (mot_data[scalarLB + 1]
         .f1->size[0] == 0) || (mot_data[scalarLB + 1].f1->size[1] == 0)) {
      i1 = rmat->size[0] * rmat->size[1];
      rmat->size[0] = A->size[0];
      rmat->size[1] = mot_data[scalarLB + 1].f1->size[0];
      emxEnsureCapacity_real_T(rmat, i1);
      rmat_data = rmat->data;
      loop_ub = A->size[0] * mot_data[scalarLB + 1].f1->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        rmat_data[i1] = 0.0;
      }
    } else {
      i1 = rmat->size[0] * rmat->size[1];
      rmat->size[0] = A->size[0];
      rmat->size[1] = mot_data[scalarLB + 1].f1->size[0];
      emxEnsureCapacity_real_T(rmat, i1);
      rmat_data = rmat->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)A->size[0],
                  (blasint)mot_data[scalarLB + 1].f1->size[0], (blasint)A->size
                  [1], 1.0, &mat_data[0], (blasint)A->size[0],
                  &mot_data[scalarLB + 1].f1->data[0], (blasint)
                  mot_data[scalarLB + 1].f1->size[0], 0.0, &rmat_data[0],
                  (blasint)A->size[0]);
    }

    /*  MM = max([sum(spdiags(mat)) sum(spdiags(rmat))]); */
    /* 'mapTF:838' diag_a = sum(spdiags(mat)); */
    spdiags(mat, A);
    b_sum(A, diag_a);
    mat_data = diag_a->data;

    /* 'mapTF:839' diag_b = sum(spdiags(rmat)); */
    spdiags(rmat, A);
    b_sum(A, diag_b);
    rmat_data = diag_b->data;

    /* 'mapTF:840' diag_overall = zeros(length(diag_a)+length(diag_b),1); */
    i1 = diag_overall->size[0];
    diag_overall->size[0] = (int)((unsigned int)diag_a->size[1] + diag_b->size[1]);
    emxEnsureCapacity_real_T(diag_overall, i1);
    diag_overall_data = diag_overall->data;
    loop_ub = (int)((unsigned int)diag_a->size[1] + diag_b->size[1]);
    for (i1 = 0; i1 < loop_ub; i1++) {
      diag_overall_data[i1] = 0.0;
    }

    /* 'mapTF:841' diag_overall(1:length(diag_a)) = diag_a; */
    loop_ub = diag_a->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      diag_overall_data[i1] = mat_data[i1];
    }

    /* 'mapTF:842' diag_overall(length(diag_a)+1:end) = diag_b; */
    if (diag_a->size[1] + 1U > (unsigned int)diag_overall->size[0]) {
      i1 = 0;
    } else {
      i1 = diag_a->size[1];
    }

    loop_ub = diag_b->size[1];
    for (vectorUB = 0; vectorUB < loop_ub; vectorUB++) {
      diag_overall_data[i1 + vectorUB] = rmat_data[vectorUB];
    }

    /* 'mapTF:843' MM = max(diag_overall); */
    MM = maximum(diag_overall);

    /* 'mapTF:844' if MM > M */
    if (MM > *M) {
      /* 'mapTF:845' M = MM; */
      *M = MM;

      /* 'mapTF:846' ind = j-1; */
      *ind = ((double)scalarLB + 2.0) - 1.0;
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
 * function process_motifs(dfn, lfn, memefn, ofn)
 */
static void process_motifs(const emxArray_char_T *dfn, const emxArray_char_T
  *lfn, const emxArray_char_T *memefn, const emxArray_char_T *ofn)
{
  static const char b_cv[5] = { 'M', 'O', 'T', 'I', 'F' };

  FILE* b_NULL;
  FILE* c_NULL;
  FILE* d_NULL;
  FILE* filestar;
  cell_wrap_3 *names_data;
  cell_wrap_5 r1;
  cell_wrap_5 *PWM2_data;
  cell_wrap_5 *P_data;
  cell_wrap_5 *cur_PWM_data;
  cell_wrap_5 *cur_PWM_tmp_data;
  emxArray_boolean_T *b_clus;
  emxArray_cell_wrap_3 *names;
  emxArray_cell_wrap_5 *P;
  emxArray_cell_wrap_5 *PWM2;
  emxArray_cell_wrap_5 *b_cur_PWM_tmp;
  emxArray_cell_wrap_5 *b_new_PWM2;
  emxArray_cell_wrap_5 *c_cur_PWM_tmp;
  emxArray_cell_wrap_5 *cur_PWM;
  emxArray_cell_wrap_5 *cur_PWM_tmp;
  emxArray_cell_wrap_5 *d_cur_PWM_tmp;
  emxArray_cell_wrap_5 *new_PWM2;
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
  double a;
  double curr_pos;
  double d;
  double idx;
  double *I_2_data;
  double *clus_data;
  double *f_data;
  double *lasso_weight_data;
  double *uid_data;
  double *vec2_data;
  double *vec_data;
  double *w_data;
  double *zscore_data;
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
  emxInit_cell_wrap_5(&cur_PWM_tmp, 1);

  /* dfn: file name for denovo motifs */
  /* lfn: file name for lasso motifs */
  /* memefn: file name for the meme input for gkmPWMlasso */
  /* ofn: output filename */
  /* 'mapTF:579' a = 1; */
  a = 1.0;

  /* 'mapTF:580' LEN = zeros(1,1); */
  /* 'mapTF:581' shift = zeros(1,1); */
  /* 'mapTF:582' [p,w] = getdenovomotif(dfn); */
  getdenovomotif(dfn, cur_PWM_tmp, f);
  f_data = f->data;
  cur_PWM_tmp_data = cur_PWM_tmp->data;

  /* 'mapTF:583' N = numel(w); */
  N = f->size[0];

  /*  fid = fopen(lfn,'r'); */
  /*  X = textscan(fid,'%f\t%f\t%s\t%f\t%f\t%f\n','delimiter', '\t', 'headerlines', 4); */
  /*  fclose(fid); */
  /* 'mapTF:588' fid = fopen(lfn, 'r'); */
  fileid = cfopen(lfn, "rb");

  /* 'mapTF:589' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF:590' fprintf("ERROR: gkmPWMlasso output file cannot be opened.\n") */
    printf("ERROR: gkmPWMlasso output file cannot be opened.\n");
    fflush(stdout);
    exit(1);
  }

  emxInit_char_T(&b_fileid, 2);
  emxInit_char_T(&c_fileid, 2);
  emxInit_char_T(&d_fileid, 2);
  emxInit_char_T(&e_fileid, 2);

  /* 'mapTF:593' fgetl(fid); */
  b_fgets(fileid, b_fileid);

  /* 'mapTF:594' fgetl(fid); */
  b_fgets(fileid, c_fileid);

  /* 'mapTF:595' fgetl(fid); */
  b_fgets(fileid, d_fileid);

  /* 'mapTF:596' fgetl(fid); */
  b_fgets(fileid, e_fileid);

  /* 'mapTF:597' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF:598' idx=0; */
  idx = 0.0;

  /* 'mapTF:599' while ~feof(fid) */
  emxFree_char_T(&e_fileid);
  emxFree_char_T(&d_fileid);
  emxFree_char_T(&c_fileid);
  emxFree_char_T(&b_fileid);
  emxInit_char_T(&f_fileid, 2);
  do {
    exitg1 = 0;
    d = b_feof(fileid);
    if (!(d != 0.0)) {
      /* 'mapTF:600' idx=idx+1; */
      idx++;

      /* 'mapTF:601' fgetl(fid); */
      b_fgets(fileid, f_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&f_fileid);
  emxInit_real_T(&clus, 1);

  /* 'mapTF:603' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF:604' clus = zeros(idx, 1); */
  match_idx = (int)idx;
  i = clus->size[0];
  clus->size[0] = (int)idx;
  emxEnsureCapacity_real_T(clus, i);
  clus_data = clus->data;
  for (i = 0; i < match_idx; i++) {
    clus_data[i] = 0.0;
  }

  emxInit_real_T(&uid, 1);

  /* 'mapTF:605' uid = zeros(idx, 1); */
  i = uid->size[0];
  uid->size[0] = (int)idx;
  emxEnsureCapacity_real_T(uid, i);
  uid_data = uid->data;
  for (i = 0; i < match_idx; i++) {
    uid_data[i] = 0.0;
  }

  emxInit_real_T(&lasso_weight, 1);

  /* 'mapTF:606' lasso_weight = zeros(idx, 1); */
  i = lasso_weight->size[0];
  lasso_weight->size[0] = (int)idx;
  emxEnsureCapacity_real_T(lasso_weight, i);
  lasso_weight_data = lasso_weight->data;
  for (i = 0; i < match_idx; i++) {
    lasso_weight_data[i] = 0.0;
  }

  emxInit_real_T(&zscore, 1);

  /* 'mapTF:607' zscore = zeros(idx, 1); */
  i = zscore->size[0];
  zscore->size[0] = (int)idx;
  emxEnsureCapacity_real_T(zscore, i);
  zscore_data = zscore->data;
  for (i = 0; i < match_idx; i++) {
    zscore_data[i] = 0.0;
  }

  /* 'mapTF:609' for cur_idx=1:idx */
  cur_idx = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&p1, 2);
  emxInit_char_T(&remain, 2);
  emxInit_char_T(&p2, 2);
  emxInit_char_T(&b_remain, 2);
  emxInit_char_T(&p4, 2);
  exitg2 = false;
  while ((!exitg2) && (cur_idx <= (int)idx - 1)) {
    /* 'mapTF:610' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF:611' if cur_line == -1 */
    hal = (cur_line->size[1] != 0);
    if (hal) {
      hal = (0 > cur_line->size[1] - 1);
    }

    if (hal) {
      exitg2 = true;
    } else {
      /* 'mapTF:614' [p1, remain] = strtok(cur_line, char(9)); */
      c_strtok(cur_line, p1, remain);

      /* 'mapTF:615' [p2, remain] = strtok(remain, char(9)); */
      c_strtok(remain, p2, b_remain);

      /* 'mapTF:616' [p3, remain] = strtok(remain, char(9)); */
      c_strtok(b_remain, cur_line, remain);

      /* 'mapTF:617' [p4, remain] = strtok(remain, char(9)); */
      c_strtok(remain, p4, b_remain);

      /* 'mapTF:618' [p5, p6] = strtok(remain, char(9)); */
      c_strtok(b_remain, cur_line, remain);

      /* 'mapTF:619' clus(cur_idx, 1) = real(str2double(p1)); */
      dc = str2double(p1);
      clus_data[cur_idx] = dc.re;

      /* 'mapTF:620' uid(cur_idx, 1) = real(str2double(p2)); */
      dc = str2double(p2);
      uid_data[cur_idx] = dc.re;

      /* 'mapTF:621' lasso_weight(cur_idx, 1) = real(str2double(p4)); */
      dc = str2double(p4);
      lasso_weight_data[cur_idx] = dc.re;

      /* 'mapTF:622' zscore(cur_idx, 1) = real(str2double(p5)); */
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

  /* 'mapTF:624' fclose(fid); */
  cfclose(fileid);

  /*  n = X{1}(end); */
  /* 'mapTF:628' n = clus(end); */
  /* 'mapTF:629' vec = zeros(n,1); */
  i = vec->size[0];
  vec->size[0] = (int)clus_data[clus->size[0] - 1];
  emxEnsureCapacity_real_T(vec, i);
  vec_data = vec->data;

  /* 'mapTF:630' vec2 = zeros(n,1); */
  i = vec2->size[0];
  vec2->size[0] = (int)clus_data[clus->size[0] - 1];
  emxEnsureCapacity_real_T(vec2, i);
  vec2_data = vec2->data;

  /* 'mapTF:631' w = [w; zeros(n,1)]; */
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

  /* 'mapTF:632' for i = 1:n */
  i = (int)clus_data[clus->size[0] - 1];
  emxInit_int32_T(&r, 1);
  emxInit_boolean_T(&b_clus, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /*      f = find(X{1}==i); */
    /*      vec(i) = X{2}(f(1)); */
    /*      vec2(i) = X{5}(f(1)); */
    /*      w(i+N) = X{4}(f(1)); */
    /* 'mapTF:637' f = find(clus == i); */
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

    /* 'mapTF:638' vec(i) = uid(f(1)); */
    vec_data[b_i] = uid_data[(int)f_data[0] - 1];

    /* 'mapTF:639' vec2(i) = zscore(f(1)); */
    vec2_data[b_i] = zscore_data[(int)f_data[0] - 1];

    /* 'mapTF:640' w(i+N) = lasso_weight(f(1)); */
    w_data[(int)((unsigned int)b_i + N)] = lasso_weight_data[(int)f_data[0] - 1];
  }

  emxFree_boolean_T(&b_clus);
  emxFree_int32_T(&r);
  emxFree_real_T(&zscore);
  emxInit_cell_wrap_5(&P, 1);
  emxInit_cell_wrap_5(&cur_PWM, 1);

  /* 'mapTF:642' P = getmotif(memefn,vec); */
  getmotif(memefn, vec, P);
  P_data = P->data;

  /*  [pp, info, len] = trim_pwm([p;P],0.25); */
  /* 'mapTF:645' denovo_len = length(p); */
  /* 'mapTF:646' database_len = length(P); */
  /* 'mapTF:647' cur_PWM = cell(denovo_len+database_len,1); */
  unnamed_idx_0 = (int)((unsigned int)cur_PWM_tmp->size[0] + P->size[0]);
  i = cur_PWM->size[0];
  cur_PWM->size[0] = (int)((unsigned int)cur_PWM_tmp->size[0] + P->size[0]);
  emxEnsureCapacity_cell_wrap_5(cur_PWM, i);
  cur_PWM_data = cur_PWM->data;
  emxFree_real_T(&vec);
  for (i = 0; i < unnamed_idx_0; i++) {
    cur_PWM_data[i].f1->size[0] = 0;
    cur_PWM_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:648' cur_PWM = coder.nullcopy(cur_PWM); */
  /* 'mapTF:649' for cur_idx=1:denovo_len */
  i = cur_PWM_tmp->size[0];
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'mapTF:650' cur_PWM{cur_idx} = p{cur_idx}; */
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

  /* 'mapTF:652' for cur_idx=denovo_len+1:denovo_len+database_len */
  i = P->size[0];
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    b_cur_idx = ((unsigned int)cur_PWM_tmp->size[0] + cur_idx) + 1U;

    /* 'mapTF:653' cur_PWM{cur_idx} = P{cur_idx-denovo_len}; */
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

  emxFree_cell_wrap_5(&P);
  emxInitStruct_cell_wrap_5(&r1);

  /* 'mapTF:655' [pp, info, len] = trim_pwm(cur_PWM,0.25); */
  trim_pwm(cur_PWM, lasso_weight, uid);
  uid_data = uid->data;
  lasso_weight_data = lasso_weight->data;
  cur_PWM_data = cur_PWM->data;

  /* 'mapTF:657' coder.varsize("PWM2"); */
  /* 'mapTF:658' PWM2 = {pp{1}}; */
  i = r1.f1->size[0] * r1.f1->size[1];
  r1.f1->size[0] = cur_PWM_data[0].f1->size[0];
  r1.f1->size[1] = cur_PWM_data[0].f1->size[1];
  emxEnsureCapacity_real_T(r1.f1, i);
  loop_ub = cur_PWM_data[0].f1->size[0] * cur_PWM_data[0].f1->size[1];
  for (i = 0; i < loop_ub; i++) {
    r1.f1->data[i] = cur_PWM_data[0].f1->data[i];
  }

  emxInit_cell_wrap_5(&PWM2, 2);
  emxInit_real_T(&LEN_2, 2);
  emxInit_real_T(&I_2, 2);
  i = PWM2->size[0] * PWM2->size[1];
  PWM2->size[0] = 1;
  PWM2->size[1] = 1;
  emxEnsureCapacity_cell_wrap_5(PWM2, i);
  PWM2_data = PWM2->data;
  emxCopyStruct_cell_wrap_5(&PWM2_data[0], &r1);

  /* 'mapTF:659' coder.varsize("LEN_2"); */
  /* 'mapTF:660' LEN_2 = zeros(1,1); */
  i = LEN_2->size[0] * LEN_2->size[1];
  LEN_2->size[0] = 1;
  LEN_2->size[1] = 1;
  emxEnsureCapacity_real_T(LEN_2, i);
  vec_data = LEN_2->data;
  vec_data[0] = 0.0;

  /* 'mapTF:661' coder.varsize("I_2"); */
  /* 'mapTF:662' I_2 = zeros(1,1); */
  i = I_2->size[0] * I_2->size[1];
  I_2->size[0] = 1;
  I_2->size[1] = 1;
  emxEnsureCapacity_real_T(I_2, i);
  I_2_data = I_2->data;
  I_2_data[0] = 0.0;

  /* 'mapTF:663' hal = true; */
  hal = true;

  /* 'mapTF:664' for ii = 1:length(w) */
  i = w->size[0];
  emxFreeStruct_cell_wrap_5(&r1);
  emxInit_cell_wrap_5(&b_cur_PWM_tmp, 1);
  emxInit_cell_wrap_5(&new_PWM2, 1);
  emxInit_cell_wrap_5(&c_cur_PWM_tmp, 1);
  emxInit_cell_wrap_5(&b_new_PWM2, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:665' if ii > N */
    if (b_i + 1 > N) {
      /* 'mapTF:666' hal = true; */
      hal = true;

      /*  [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]); */
      /* 'mapTF:668' PWM2_len = length(PWM2); */
      u0 = PWM2->size[0];
      if (u0 < 1) {
        u0 = 1;
      }

      /* 'mapTF:669' cur_PWM_tmp = cell(PWM2_len+1,1); */
      unnamed_idx_0 = u0 + 1;
      i1 = c_cur_PWM_tmp->size[0];
      c_cur_PWM_tmp->size[0] = u0 + 1;
      emxEnsureCapacity_cell_wrap_5(c_cur_PWM_tmp, i1);
      cur_PWM_tmp_data = c_cur_PWM_tmp->data;
      for (i1 = 0; i1 < unnamed_idx_0; i1++) {
        cur_PWM_tmp_data[i1].f1->size[0] = 0;
        cur_PWM_tmp_data[i1].f1->size[1] = 0;
      }

      /* 'mapTF:670' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
      i1 = cur_PWM_tmp->size[0];
      cur_PWM_tmp->size[0] = c_cur_PWM_tmp->size[0];
      emxEnsureCapacity_cell_wrap_5(cur_PWM_tmp, i1);
      cur_PWM_tmp_data = cur_PWM_tmp->data;

      /* 'mapTF:671' cur_PWM_tmp{1} = pp{ii}; */
      i1 = cur_PWM_tmp_data[0].f1->size[0] * cur_PWM_tmp_data[0].f1->size[1];
      cur_PWM_tmp_data[0].f1->size[0] = cur_PWM_data[b_i].f1->size[0];
      cur_PWM_tmp_data[0].f1->size[1] = cur_PWM_data[b_i].f1->size[1];
      emxEnsureCapacity_real_T(cur_PWM_tmp_data[0].f1, i1);
      loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cur_PWM_tmp_data[0].f1->data[i1] = cur_PWM_data[b_i].f1->data[i1];
      }

      /* 'mapTF:672' for cur_idx=1:PWM2_len */
      for (cur_idx = 0; cur_idx < u0; cur_idx++) {
        /* 'mapTF:673' cur_PWM_tmp{cur_idx+1} = PWM2{cur_idx}; */
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

      /* 'mapTF:675' cur_length = zeros(length(LEN_2)+1,1); */
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

      /* 'mapTF:676' cur_length(1) = len(ii); */
      f_data[0] = uid_data[b_i];

      /* 'mapTF:677' cur_length(2:end) = LEN_2; */
      for (i1 = 0; i1 < u0; i1++) {
        f_data[i1 + 1] = vec_data[i1];
      }

      /* 'mapTF:678' [~,cor] = ppmsim(cur_PWM_tmp, cur_length); */
      ppmsim(cur_PWM_tmp, f, &curr_pos, &idx);

      /* 'mapTF:680' if cor > 0.8 || vec2(ii-N) < 1.5 */
      if ((idx > 0.8) || (vec2_data[b_i - N] < 1.5)) {
        /* 'mapTF:681' hal = false; */
        hal = false;
      }
    } else if ((b_i + 1 > 1) && (a > 2.0)) {
      /* 'mapTF:683' elseif ii > 1 && a > 2 */
      /* 'mapTF:684' hal = true; */
      hal = true;

      /*  [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]); */
      /* 'mapTF:686' PWM2_len = length(PWM2); */
      u0 = PWM2->size[0];
      if (u0 < 1) {
        u0 = 1;
      }

      /* 'mapTF:687' cur_PWM_tmp = cell(PWM2_len+1,1); */
      unnamed_idx_0 = u0 + 1;
      i1 = b_cur_PWM_tmp->size[0];
      b_cur_PWM_tmp->size[0] = u0 + 1;
      emxEnsureCapacity_cell_wrap_5(b_cur_PWM_tmp, i1);
      P_data = b_cur_PWM_tmp->data;
      for (i1 = 0; i1 < unnamed_idx_0; i1++) {
        P_data[i1].f1->size[0] = 0;
        P_data[i1].f1->size[1] = 0;
      }

      /* 'mapTF:688' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
      i1 = cur_PWM_tmp->size[0];
      cur_PWM_tmp->size[0] = b_cur_PWM_tmp->size[0];
      emxEnsureCapacity_cell_wrap_5(cur_PWM_tmp, i1);
      cur_PWM_tmp_data = cur_PWM_tmp->data;

      /* 'mapTF:689' cur_PWM_tmp{1} = pp{ii}; */
      i1 = cur_PWM_tmp_data[0].f1->size[0] * cur_PWM_tmp_data[0].f1->size[1];
      cur_PWM_tmp_data[0].f1->size[0] = cur_PWM_data[b_i].f1->size[0];
      cur_PWM_tmp_data[0].f1->size[1] = cur_PWM_data[b_i].f1->size[1];
      emxEnsureCapacity_real_T(cur_PWM_tmp_data[0].f1, i1);
      loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cur_PWM_tmp_data[0].f1->data[i1] = cur_PWM_data[b_i].f1->data[i1];
      }

      /* 'mapTF:690' for cur_idx=1:PWM2_len */
      for (cur_idx = 0; cur_idx < u0; cur_idx++) {
        /* 'mapTF:691' cur_PWM_tmp{cur_idx+1} = PWM2{cur_idx}; */
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

      /* 'mapTF:693' cur_length = zeros(length(LEN_2)+1,1); */
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

      /* 'mapTF:694' cur_length(1) = len(ii); */
      f_data[0] = uid_data[b_i];

      /* 'mapTF:695' cur_length(2:end) = LEN_2; */
      for (i1 = 0; i1 < u0; i1++) {
        f_data[i1 + 1] = vec_data[i1];
      }

      /* 'mapTF:696' [~,cor] = ppmsim(cur_PWM_tmp, cur_length); */
      ppmsim(cur_PWM_tmp, f, &curr_pos, &idx);

      /* 'mapTF:698' if cor > 0.8 */
      if (idx > 0.8) {
        /* 'mapTF:699' hal = false; */
        hal = false;
      }
    }

    /* 'mapTF:702' if hal && w(ii) > 0 && len(ii) >= 6 */
    if (hal && (w_data[b_i] > 0.0) && (uid_data[b_i] >= 6.0)) {
      /* 'mapTF:703' if len(ii) > 10 */
      if (uid_data[b_i] > 10.0) {
        /* 'mapTF:704' if info(ii)/len(ii) > 0.7 */
        if (lasso_weight_data[b_i] / uid_data[b_i] > 0.7) {
          /*  PWM2{a} = pp{ii}; */
          /*  LEN_2(a) = len(ii); */
          /*  I_2(a) = info(ii); */
          /*  a = a+1; */
          /*  PWM2{a} = rot90(PWM2{a-1},2); */
          /*  LEN_2(a) = LEN_2(a-1); */
          /*  I_2(a) = I_2(a-1); */
          /*  a = a+1; */
          /* 'mapTF:714' new_len = a+1; */
          /* 'mapTF:715' new_PWM2 = cell(new_len,1); */
          unnamed_idx_0 = (int)(a + 1.0);
          i1 = b_new_PWM2->size[0];
          b_new_PWM2->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_cell_wrap_5(b_new_PWM2, i1);
          P_data = b_new_PWM2->data;
          for (i1 = 0; i1 < unnamed_idx_0; i1++) {
            P_data[i1].f1->size[0] = 0;
            P_data[i1].f1->size[1] = 0;
          }

          /* 'mapTF:716' new_PWM2 = coder.nullcopy(new_PWM2); */
          i1 = cur_PWM_tmp->size[0];
          cur_PWM_tmp->size[0] = b_new_PWM2->size[0];
          emxEnsureCapacity_cell_wrap_5(cur_PWM_tmp, i1);
          cur_PWM_tmp_data = cur_PWM_tmp->data;

          /* 'mapTF:717' for idx = 1:new_len */
          for (match_idx = 0; match_idx < unnamed_idx_0; match_idx++) {
            /* 'mapTF:718' if idx <= a-1 */
            if ((double)match_idx + 1.0 <= a - 1.0) {
              /* 'mapTF:719' new_PWM2{idx} = PWM2{idx}; */
              i1 = cur_PWM_tmp_data[match_idx].f1->size[0] *
                cur_PWM_tmp_data[match_idx].f1->size[1];
              cur_PWM_tmp_data[match_idx].f1->size[0] = PWM2_data[match_idx].
                f1->size[0];
              cur_PWM_tmp_data[match_idx].f1->size[1] = PWM2_data[match_idx].
                f1->size[1];
              emxEnsureCapacity_real_T(cur_PWM_tmp_data[match_idx].f1, i1);
              loop_ub = PWM2_data[match_idx].f1->size[0] * PWM2_data[match_idx].
                f1->size[1];
              for (i1 = 0; i1 < loop_ub; i1++) {
                cur_PWM_tmp_data[match_idx].f1->data[i1] = PWM2_data[match_idx].
                  f1->data[i1];
              }
            } else if ((double)match_idx + 1.0 == a) {
              /* 'mapTF:720' elseif idx == a */
              /* 'mapTF:721' new_PWM2{idx} = pp{ii}; */
              i1 = cur_PWM_tmp_data[match_idx].f1->size[0] *
                cur_PWM_tmp_data[match_idx].f1->size[1];
              cur_PWM_tmp_data[match_idx].f1->size[0] = cur_PWM_data[b_i]
                .f1->size[0];
              cur_PWM_tmp_data[match_idx].f1->size[1] = cur_PWM_data[b_i]
                .f1->size[1];
              emxEnsureCapacity_real_T(cur_PWM_tmp_data[match_idx].f1, i1);
              loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i]
                .f1->size[1];
              for (i1 = 0; i1 < loop_ub; i1++) {
                cur_PWM_tmp_data[match_idx].f1->data[i1] = cur_PWM_data[b_i].
                  f1->data[i1];
              }
            } else if ((double)match_idx + 1.0 == a + 1.0) {
              /* 'mapTF:722' elseif idx == a+1 */
              /* 'mapTF:723' new_PWM2{a+1} = rot90(pp{ii},2); */
              rot90(cur_PWM_data[b_i].f1, cur_PWM_tmp_data[(int)(a + 1.0) - 1].
                    f1);
            }
          }

          /* 'mapTF:726' PWM2_len = length(new_PWM2); */
          /* 'mapTF:727' PWM2 = cell(PWM2_len, 1); */
          match_idx = cur_PWM_tmp->size[0];
          i1 = PWM2->size[0] * PWM2->size[1];
          PWM2->size[0] = cur_PWM_tmp->size[0];
          PWM2->size[1] = 1;
          emxEnsureCapacity_cell_wrap_5(PWM2, i1);
          PWM2_data = PWM2->data;
          for (i1 = 0; i1 < match_idx; i1++) {
            PWM2_data[i1].f1->size[0] = 0;
            PWM2_data[i1].f1->size[1] = 0;
          }

          /* 'mapTF:728' PWM2 = coder.nullcopy(PWM2); */
          /* 'mapTF:729' for cur_idx = 1:PWM2_len */
          i1 = cur_PWM_tmp->size[0];
          for (cur_idx = 0; cur_idx < i1; cur_idx++) {
            /* 'mapTF:730' PWM2{cur_idx} = new_PWM2{cur_idx}; */
            match_idx = PWM2_data[cur_idx].f1->size[0] * PWM2_data[cur_idx]
              .f1->size[1];
            PWM2_data[cur_idx].f1->size[0] = cur_PWM_tmp_data[cur_idx].f1->size
              [0];
            PWM2_data[cur_idx].f1->size[1] = cur_PWM_tmp_data[cur_idx].f1->size
              [1];
            emxEnsureCapacity_real_T(PWM2_data[cur_idx].f1, match_idx);
            loop_ub = cur_PWM_tmp_data[cur_idx].f1->size[0] *
              cur_PWM_tmp_data[cur_idx].f1->size[1];
            for (match_idx = 0; match_idx < loop_ub; match_idx++) {
              PWM2_data[cur_idx].f1->data[match_idx] = cur_PWM_tmp_data[cur_idx]
                .f1->data[match_idx];
            }
          }

          /* 'mapTF:732' LEN_2_new = zeros(a+1,1); */
          i1 = f->size[0];
          f->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_real_T(f, i1);
          f_data = f->data;
          for (i1 = 0; i1 < unnamed_idx_0; i1++) {
            f_data[i1] = 0.0;
          }

          /* 'mapTF:733' LEN_2_new(1:a-1) = LEN_2(1:a-1); */
          if (1.0 > a - 1.0) {
            loop_ub = 0;
          } else {
            loop_ub = (int)(a - 1.0);
          }

          for (i1 = 0; i1 < loop_ub; i1++) {
            f_data[i1] = vec_data[i1];
          }

          /* 'mapTF:734' LEN_2_new(a:a+1) = [len(ii), len(ii)]; */
          f_data[(int)a - 1] = uid_data[b_i];
          f_data[(int)(a + 1.0) - 1] = uid_data[b_i];

          /* 'mapTF:736' I_2_new = zeros(a+1,1); */
          i1 = clus->size[0];
          clus->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_real_T(clus, i1);
          clus_data = clus->data;
          for (i1 = 0; i1 < unnamed_idx_0; i1++) {
            clus_data[i1] = 0.0;
          }

          /* 'mapTF:737' I_2_new(1:a-1) = I_2(1:a-1); */
          if (1.0 > a - 1.0) {
            loop_ub = 0;
          } else {
            loop_ub = (int)(a - 1.0);
          }

          for (i1 = 0; i1 < loop_ub; i1++) {
            clus_data[i1] = I_2_data[i1];
          }

          /* 'mapTF:738' I_2_new(a:a+1) = [info(ii), info(ii)]; */
          clus_data[(int)a - 1] = lasso_weight_data[b_i];
          clus_data[(int)(a + 1.0) - 1] = lasso_weight_data[b_i];

          /* 'mapTF:740' LEN_2 = LEN_2_new; */
          i1 = LEN_2->size[0] * LEN_2->size[1];
          LEN_2->size[0] = f->size[0];
          LEN_2->size[1] = 1;
          emxEnsureCapacity_real_T(LEN_2, i1);
          vec_data = LEN_2->data;
          loop_ub = f->size[0];
          for (i1 = 0; i1 < loop_ub; i1++) {
            vec_data[i1] = f_data[i1];
          }

          /* 'mapTF:741' I_2 = I_2_new; */
          i1 = I_2->size[0] * I_2->size[1];
          I_2->size[0] = clus->size[0];
          I_2->size[1] = 1;
          emxEnsureCapacity_real_T(I_2, i1);
          I_2_data = I_2->data;
          loop_ub = clus->size[0];
          for (i1 = 0; i1 < loop_ub; i1++) {
            I_2_data[i1] = clus_data[i1];
          }

          /* 'mapTF:742' a = a+1; */
          a++;

          /* 'mapTF:743' a = a+1; */
          a++;
        }
      } else if ((lasso_weight_data[b_i] > 6.0) || (lasso_weight_data[b_i] /
                  uid_data[b_i] > 1.0)) {
        /* 'mapTF:745' elseif info(ii) > 6 || info(ii)/len(ii) > 1 */
        /*  PWM2{a} = pp{ii}; */
        /*  LEN_2(a) = len(ii); */
        /*  I_2(a) = info(ii); */
        /*  a = a+1; */
        /*  PWM2{a} = rot90(PWM2{a-1},2); */
        /*  LEN_2(a) = LEN_2(a-1); */
        /*  I_2(a) = I_2(a-1); */
        /*  a = a+1; */
        /* 'mapTF:754' new_len = a+1; */
        /* 'mapTF:755' new_PWM2 = cell(new_len,1); */
        unnamed_idx_0 = (int)(a + 1.0);
        i1 = new_PWM2->size[0];
        new_PWM2->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_cell_wrap_5(new_PWM2, i1);
        P_data = new_PWM2->data;
        for (i1 = 0; i1 < unnamed_idx_0; i1++) {
          P_data[i1].f1->size[0] = 0;
          P_data[i1].f1->size[1] = 0;
        }

        /* 'mapTF:756' new_PWM2 = coder.nullcopy(new_PWM2); */
        i1 = cur_PWM_tmp->size[0];
        cur_PWM_tmp->size[0] = new_PWM2->size[0];
        emxEnsureCapacity_cell_wrap_5(cur_PWM_tmp, i1);
        cur_PWM_tmp_data = cur_PWM_tmp->data;

        /* 'mapTF:757' for idx = 1:new_len */
        for (match_idx = 0; match_idx < unnamed_idx_0; match_idx++) {
          /* 'mapTF:758' if idx <= a-1 */
          if ((double)match_idx + 1.0 <= a - 1.0) {
            /* 'mapTF:759' new_PWM2{idx} = PWM2{idx}; */
            i1 = cur_PWM_tmp_data[match_idx].f1->size[0] *
              cur_PWM_tmp_data[match_idx].f1->size[1];
            cur_PWM_tmp_data[match_idx].f1->size[0] = PWM2_data[match_idx]
              .f1->size[0];
            cur_PWM_tmp_data[match_idx].f1->size[1] = PWM2_data[match_idx]
              .f1->size[1];
            emxEnsureCapacity_real_T(cur_PWM_tmp_data[match_idx].f1, i1);
            loop_ub = PWM2_data[match_idx].f1->size[0] * PWM2_data[match_idx].
              f1->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              cur_PWM_tmp_data[match_idx].f1->data[i1] = PWM2_data[match_idx].
                f1->data[i1];
            }
          } else if ((double)match_idx + 1.0 == a) {
            /* 'mapTF:760' elseif idx == a */
            /* 'mapTF:761' new_PWM2{idx} = pp{ii}; */
            i1 = cur_PWM_tmp_data[match_idx].f1->size[0] *
              cur_PWM_tmp_data[match_idx].f1->size[1];
            cur_PWM_tmp_data[match_idx].f1->size[0] = cur_PWM_data[b_i].f1->
              size[0];
            cur_PWM_tmp_data[match_idx].f1->size[1] = cur_PWM_data[b_i].f1->
              size[1];
            emxEnsureCapacity_real_T(cur_PWM_tmp_data[match_idx].f1, i1);
            loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i].f1->
              size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              cur_PWM_tmp_data[match_idx].f1->data[i1] = cur_PWM_data[b_i]
                .f1->data[i1];
            }
          } else if ((double)match_idx + 1.0 == a + 1.0) {
            /* 'mapTF:762' elseif idx == a+1 */
            /* 'mapTF:763' new_PWM2{a+1} = rot90(pp{ii},2); */
            rot90(cur_PWM_data[b_i].f1, cur_PWM_tmp_data[(int)(a + 1.0) - 1].f1);
          }
        }

        /* 'mapTF:766' PWM2_len = length(new_PWM2); */
        /* 'mapTF:767' PWM2 = cell(PWM2_len, 1); */
        match_idx = cur_PWM_tmp->size[0];
        i1 = PWM2->size[0] * PWM2->size[1];
        PWM2->size[0] = cur_PWM_tmp->size[0];
        PWM2->size[1] = 1;
        emxEnsureCapacity_cell_wrap_5(PWM2, i1);
        PWM2_data = PWM2->data;
        for (i1 = 0; i1 < match_idx; i1++) {
          PWM2_data[i1].f1->size[0] = 0;
          PWM2_data[i1].f1->size[1] = 0;
        }

        /* 'mapTF:768' PWM2 = coder.nullcopy(PWM2); */
        /* 'mapTF:769' for cur_idx = 1:PWM2_len */
        i1 = cur_PWM_tmp->size[0];
        for (cur_idx = 0; cur_idx < i1; cur_idx++) {
          /* 'mapTF:770' PWM2{cur_idx} = new_PWM2{cur_idx}; */
          match_idx = PWM2_data[cur_idx].f1->size[0] * PWM2_data[cur_idx]
            .f1->size[1];
          PWM2_data[cur_idx].f1->size[0] = cur_PWM_tmp_data[cur_idx].f1->size[0];
          PWM2_data[cur_idx].f1->size[1] = cur_PWM_tmp_data[cur_idx].f1->size[1];
          emxEnsureCapacity_real_T(PWM2_data[cur_idx].f1, match_idx);
          loop_ub = cur_PWM_tmp_data[cur_idx].f1->size[0] *
            cur_PWM_tmp_data[cur_idx].f1->size[1];
          for (match_idx = 0; match_idx < loop_ub; match_idx++) {
            PWM2_data[cur_idx].f1->data[match_idx] = cur_PWM_tmp_data[cur_idx].
              f1->data[match_idx];
          }
        }

        /* 'mapTF:772' LEN_2_new = zeros(a+1,1); */
        i1 = f->size[0];
        f->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        for (i1 = 0; i1 < unnamed_idx_0; i1++) {
          f_data[i1] = 0.0;
        }

        /* 'mapTF:773' LEN_2_new(1:a-1) = LEN_2(1:a-1); */
        if (1.0 > a - 1.0) {
          loop_ub = 0;
        } else {
          loop_ub = (int)(a - 1.0);
        }

        for (i1 = 0; i1 < loop_ub; i1++) {
          f_data[i1] = vec_data[i1];
        }

        /* 'mapTF:774' LEN_2_new(a:a+1) = [len(ii), len(ii)]; */
        f_data[(int)a - 1] = uid_data[b_i];
        f_data[(int)(a + 1.0) - 1] = uid_data[b_i];

        /* 'mapTF:776' I_2_new = zeros(a+1,1); */
        i1 = clus->size[0];
        clus->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_real_T(clus, i1);
        clus_data = clus->data;
        for (i1 = 0; i1 < unnamed_idx_0; i1++) {
          clus_data[i1] = 0.0;
        }

        /* 'mapTF:777' I_2_new(1:a-1) = I_2(1:a-1); */
        if (1.0 > a - 1.0) {
          loop_ub = 0;
        } else {
          loop_ub = (int)(a - 1.0);
        }

        for (i1 = 0; i1 < loop_ub; i1++) {
          clus_data[i1] = I_2_data[i1];
        }

        /* 'mapTF:778' I_2_new(a:a+1) = [info(ii), info(ii)]; */
        clus_data[(int)a - 1] = lasso_weight_data[b_i];
        clus_data[(int)(a + 1.0) - 1] = lasso_weight_data[b_i];

        /* 'mapTF:780' LEN_2 = LEN_2_new; */
        i1 = LEN_2->size[0] * LEN_2->size[1];
        LEN_2->size[0] = f->size[0];
        LEN_2->size[1] = 1;
        emxEnsureCapacity_real_T(LEN_2, i1);
        vec_data = LEN_2->data;
        loop_ub = f->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          vec_data[i1] = f_data[i1];
        }

        /* 'mapTF:781' I_2 = I_2_new; */
        i1 = I_2->size[0] * I_2->size[1];
        I_2->size[0] = clus->size[0];
        I_2->size[1] = 1;
        emxEnsureCapacity_real_T(I_2, i1);
        I_2_data = I_2->data;
        loop_ub = clus->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          I_2_data[i1] = clus_data[i1];
        }

        /* 'mapTF:782' a = a+1; */
        a++;

        /* 'mapTF:783' a = a+1; */
        a++;
      }
    }
  }

  emxFree_cell_wrap_5(&b_new_PWM2);
  emxFree_cell_wrap_5(&c_cur_PWM_tmp);
  emxFree_cell_wrap_5(&new_PWM2);
  emxFree_cell_wrap_5(&cur_PWM);
  emxFree_real_T(&vec2);
  emxFree_real_T(&uid);
  emxFree_real_T(&w);

  /* 'mapTF:788' num2 = length(strfind(fileread(memefn),'MOTIF')); */
  fileread(memefn, cur_line);
  cur_line_data = cur_line->data;
  if (cur_line->size[1] == 0) {
    match_idx = 0;
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
    match_idx = match_out->size[1];
    emxFree_int32_T(&match_out);
  }

  /* 'mapTF:789' [p,names] = getmotif(memefn,1:num2); */
  emxInit_real_T(&y, 2);
  if (match_idx < 1) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = match_idx;
    emxEnsureCapacity_real_T(y, i);
    zscore_data = y->data;
    loop_ub = match_idx - 1;
    for (i = 0; i <= loop_ub; i++) {
      zscore_data[i] = (double)i + 1.0;
    }
  }

  emxInit_cell_wrap_3(&names);
  b_getmotif(memefn, y, cur_PWM_tmp, names);
  names_data = names->data;

  /* 'mapTF:790' [p,info,lenvec] = trim_pwm(p,0.25); */
  trim_pwm(cur_PWM_tmp, lasso_weight, f);
  f_data = f->data;
  cur_PWM_tmp_data = cur_PWM_tmp->data;

  /* 'mapTF:791' fid = fopen([ofn '_motifs.out'], 'w'); */
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

  /* 'mapTF:792' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF:793' fprintf("ERROR: Cannot create the combined motif file\n"); */
    printf("ERROR: Cannot create the combined motif file\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:795' a = 1; */
  a = 1.0;

  /* 'mapTF:796' for i = 1:length(PWM2) */
  u0 = PWM2->size[0];
  if (u0 < 1) {
    u0 = 1;
  }

  unnamed_idx_0 = cur_PWM_tmp->size[0] + 1;
  i = cur_PWM_tmp->size[0];
  loop_ub = f->size[0];
  emxInit_cell_wrap_5(&d_cur_PWM_tmp, 1);
  for (b_i = 0; b_i < u0; b_i++) {
    /*  [ind, r] = ppmsim([PWM2{i};p], [LEN_2(i);lenvec]); */
    /* 'mapTF:798' p_len = length(p); */
    /* 'mapTF:799' cur_PWM_tmp = cell(p_len+1,1); */
    i1 = d_cur_PWM_tmp->size[0];
    d_cur_PWM_tmp->size[0] = unnamed_idx_0;
    emxEnsureCapacity_cell_wrap_5(d_cur_PWM_tmp, i1);
    P_data = d_cur_PWM_tmp->data;
    for (i1 = 0; i1 < unnamed_idx_0; i1++) {
      P_data[i1].f1->size[0] = 0;
      P_data[i1].f1->size[1] = 0;
    }

    /* 'mapTF:800' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
    i1 = b_cur_PWM_tmp->size[0];
    b_cur_PWM_tmp->size[0] = d_cur_PWM_tmp->size[0];
    emxEnsureCapacity_cell_wrap_5(b_cur_PWM_tmp, i1);
    P_data = b_cur_PWM_tmp->data;

    /* 'mapTF:801' cur_PWM_tmp{1} = PWM2{i}; */
    i1 = P_data[0].f1->size[0] * P_data[0].f1->size[1];
    P_data[0].f1->size[0] = PWM2_data[b_i].f1->size[0];
    P_data[0].f1->size[1] = PWM2_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(P_data[0].f1, i1);
    match_idx = PWM2_data[b_i].f1->size[0] * PWM2_data[b_i].f1->size[1];
    for (i1 = 0; i1 < match_idx; i1++) {
      P_data[0].f1->data[i1] = PWM2_data[b_i].f1->data[i1];
    }

    /* 'mapTF:802' for cur_idx=1:p_len */
    for (cur_idx = 0; cur_idx < i; cur_idx++) {
      /* 'mapTF:803' cur_PWM_tmp{cur_idx+1} = p{cur_idx}; */
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

    /* 'mapTF:805' [ind, r] = ppmsim(cur_PWM_tmp, [LEN_2(i);lenvec]); */
    i1 = clus->size[0];
    clus->size[0] = f->size[0] + 1;
    emxEnsureCapacity_real_T(clus, i1);
    clus_data = clus->data;
    clus_data[0] = vec_data[b_i];
    for (i1 = 0; i1 < loop_ub; i1++) {
      clus_data[i1 + 1] = f_data[i1];
    }

    ppmsim(b_cur_PWM_tmp, clus, &curr_pos, &idx);

    /* 'mapTF:807' if r > 0.80 */
    if (idx > 0.8) {
      /* 'mapTF:808' fprintf(fid,'MOTIF %d\n%s\n%d\n', int32(a), names{ind}, int32(LEN_2(i))); */
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

      /* 'mapTF:809' a = a+1; */
      a++;

      /* 'mapTF:810' for j = 1:LEN_2(i) */
      i1 = (int)vec_data[b_i];
      if (0 <= i1 - 1) {
        d_NULL = NULL;
      }

      for (N = 0; N < i1; N++) {
        /* 'mapTF:811' fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4)); */
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

      /* 'mapTF:813' fprintf(fid, '\n'); */
      b_NULL = NULL;
      getfilestar(fileid, &filestar, &hal);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "\n");
        if (hal) {
          fflush(filestar);
        }
      }
    } else {
      d = vec_data[b_i];
      if (I_2_data[b_i] / d > 1.0) {
        /* 'mapTF:814' elseif I_2(i)/LEN_2(i) > 1 */
        /* 'mapTF:815' fprintf(fid,'MOTIF %d\n%s\n%d\n', int32(a), consen(PWM2{i}, LEN_2(i)), int32(LEN_2(i))); */
        b_NULL = NULL;
        getfilestar(fileid, &filestar, &hal);
        if (!(filestar == b_NULL)) {
          fprintf(filestar, "MOTIF %d\n%s\n%d\n", (int)a, "", (int)rt_roundd_snf
                  (d));
          if (hal) {
            fflush(filestar);
          }
        }

        /* 'mapTF:816' a = a+1; */
        a++;

        /* 'mapTF:817' for j = 1:LEN_2(i) */
        i1 = (int)d;
        if (0 <= (int)d - 1) {
          c_NULL = NULL;
        }

        for (N = 0; N < i1; N++) {
          /* 'mapTF:818' fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4)); */
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

        /* 'mapTF:820' fprintf(fid, '\n'); */
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

  emxFree_cell_wrap_5(&d_cur_PWM_tmp);
  emxFree_cell_wrap_5(&cur_PWM_tmp);
  emxFree_cell_wrap_3(&names);
  emxFree_cell_wrap_5(&b_cur_PWM_tmp);
  emxFree_real_T(&I_2);
  emxFree_real_T(&LEN_2);
  emxFree_cell_wrap_5(&PWM2);
  emxFree_real_T(&f);
  emxFree_char_T(&cur_line);
  emxFree_real_T(&clus);

  /* 'mapTF:823' fclose(fid); */
  cfclose(fileid);
}

/*
 * function varscore = scoreseqkmer(PWM3, lPWM3, Lmat, ss, Smat, l_svm, k_svm, ofn, dsvm)
 */
static void scoreseqkmer(const emxArray_cell_wrap_4 *PWM3, const
  emxArray_cell_wrap_4 *lPWM3, const emxArray_real_T *Lmat, const
  emxArray_real_T *ss, const emxArray_cell_wrap_5 *Smat, double l_svm, const
  emxArray_real_T *dsvm, emxArray_real_T *varscore)
{
  static const signed char varc[12] = { 2, 1, 1, 1, 3, 3, 2, 2, 4, 4, 4, 3 };

  const cell_wrap_4 *PWM3_data;
  const cell_wrap_4 *lPWM3_data;
  const cell_wrap_5 *Smat_data;
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
  double Lind;
  double M;
  double d;
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
  lPWM3_data = lPWM3->data;
  PWM3_data = PWM3->data;

  /* 'mapTF:294' varc = [2 3 4; 1 3 4;1 2 4;1 2 3]; */
  /* 'mapTF:295' O = ones(1,l_svm-1); */
  /* 'mapTF:296' n = numel(Lmat)/4; */
  n = (double)(Lmat->size[0] << 2) / 4.0;

  /* 'mapTF:297' varscore = zeros(n,1); */
  i = varscore->size[0];
  i1 = (int)n;
  varscore->size[0] = (int)n;
  emxEnsureCapacity_real_T(varscore, i);
  varscore_data = varscore->data;

  /* 'mapTF:298' for i = 1:n */
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
    /* 'mapTF:299' M = Lmat(i,1); */
    M = Lmat_data[b_i];

    /* 'mapTF:300' ind = [O ss(Lmat(i,2):Lmat(i,3)) O]; */
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

    /* 'mapTF:301' DSVM = dsvm(Lmat(i,2):Lmat(i,3),:); */
    d = Lmat_data[b_i + Lmat->size[0]];
    n = Lmat_data[b_i + Lmat->size[0] * 2];
    if (d > n) {
      i = -1;
      i2 = -1;
    } else {
      i = (int)d - 2;
      i2 = (int)n - 1;
    }

    /* 'mapTF:302' L = Lmat(i,3)-Lmat(i,2)+1; */
    L = (n - d) + 1.0;

    /* 'mapTF:303' matscore = zeros(L,3); */
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

    /* 'mapTF:304' for ii = 1:L */
    if (0 <= (int)L - 1) {
      loop_ub_tmp = (int)(2.0 * l_svm - 1.0);
      i5 = (int)l_svm;
    }

    for (ii = 0; ii < i3; ii++) {
      /* 'mapTF:305' Lind = l_svm-1+ii; */
      Lind = (l_svm - 1.0) + ((double)ii + 1.0);

      /* 'mapTF:306' scores = zeros(2*l_svm-1,1); */
      i4 = scores->size[0];
      scores->size[0] = (int)(2.0 * l_svm - 1.0);
      emxEnsureCapacity_real_T(scores, i4);
      scores_data = scores->data;
      for (i4 = 0; i4 < loop_ub_tmp; i4++) {
        scores_data[i4] = 0.0;
      }

      /* 'mapTF:307' for j = 1:2*l_svm-1 */
      for (j = 0; j < loop_ub_tmp; j++) {
        /* 'mapTF:308' scores(j) = lPWM3{M}(Lind-l_svm+j,ind(Lind-l_svm+j)); */
        nx = (int)((Lind - l_svm) + ((double)j + 1.0));
        scores_data[j] = lPWM3_data[(int)M - 1].f1->data[(nx + lPWM3_data[(int)M
          - 1].f1->size[0] * ((int)ind_data[nx - 1] - 1)) - 1];
      }

      /* 'mapTF:310' evec = 0; */
      L = 0.0;

      /* 'mapTF:311' scores(l_svm) = 0; */
      scores_data[(int)l_svm - 1] = 0.0;

      /* 'mapTF:312' for j = 1:l_svm */
      for (j = 0; j < i5; j++) {
        /* 'mapTF:313' evec = evec+sum(exp(Smat{l_svm-j+1}*scores(j:l_svm+j-1))); */
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

        L += blockedSummation(y, y->size[0]);
      }

      /* 'mapTF:315' for iii = 1:3 */
      /* 'mapTF:316' V = PWM3{M}(Lind,varc(ind(Lind),iii))-PWM3{M}(Lind,ind(Lind)); */
      /* 'mapTF:317' matscore(ii,iii) = evec*V; */
      nx = (int)ind_data[(int)Lind - 1];
      n = PWM3_data[(int)M - 1].f1->data[((int)Lind + PWM3_data[(int)M - 1]
        .f1->size[0] * (nx - 1)) - 1];
      matscore_data[ii] = L * (PWM3_data[(int)M - 1].f1->data[((int)Lind +
        PWM3_data[(int)M - 1].f1->size[0] * (varc[nx - 1] - 1)) - 1] - n);

      /* 'mapTF:316' V = PWM3{M}(Lind,varc(ind(Lind),iii))-PWM3{M}(Lind,ind(Lind)); */
      /* 'mapTF:317' matscore(ii,iii) = evec*V; */
      matscore_data[ii + matscore->size[0]] = L * (PWM3_data[(int)M - 1]
        .f1->data[((int)Lind + PWM3_data[(int)M - 1].f1->size[0] * (varc[nx + 3]
        - 1)) - 1] - n);

      /* 'mapTF:316' V = PWM3{M}(Lind,varc(ind(Lind),iii))-PWM3{M}(Lind,ind(Lind)); */
      /* 'mapTF:317' matscore(ii,iii) = evec*V; */
      matscore_data[ii + matscore->size[0] * 2] = L * (PWM3_data[(int)M - 1].
        f1->data[((int)Lind + PWM3_data[(int)M - 1].f1->size[0] * (varc[nx + 7]
        - 1)) - 1] - n);
    }

    /* 'mapTF:320' varscore(i) = ip(matscore(:), DSVM(:)); */
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

    /* 'mapTF:340' c = x'*y/sqrt(x'*x)/sqrt(y'*y); */
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
      Lind = 0.0;
    } else {
      Lind = cblas_ddot((blasint)((i2 - i) * 3), &y_data[0], (blasint)1,
                        &y_data[0], (blasint)1);
    }

    L = sqrt(L);
    Lind = sqrt(Lind);
    varscore_data[b_i] = n / L / Lind;
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
 * function [P,V,seqindmat,seqout,seq] = seq2pv(sfn, wfn, l_svm)
 */
static void seq2pv(const emxArray_char_T *sfn, const emxArray_char_T *wfn,
                   double l_svm, emxArray_cell_wrap_0 *P, emxArray_cell_wrap_1
                   *V, emxArray_cell_wrap_0 *seqindmat, emxArray_cell_wrap_2
                   *seqout, emxArray_cell_wrap_3 *seq)
{
  static const signed char mat[12] = { 1, 0, 0, 0, 2, 2, 1, 1, 3, 3, 3, 2 };

  __m128d r;
  cell_wrap_0 *P_data;
  cell_wrap_0 *seqindmat_data;
  cell_wrap_1 *V_data;
  cell_wrap_2 *seqout_data;
  cell_wrap_3 *seq_data;
  cell_wrap_3 *sequences_data;
  emxArray_cell_wrap_3 *sequences;
  emxArray_char_T *b_fileid;
  emxArray_char_T *c_fileid;
  emxArray_char_T *cur_alpha;
  emxArray_char_T *cur_line;
  emxArray_char_T *cur_seq;
  emxArray_real_T *R;
  emxArray_real_T *RR;
  emxArray_real_T *S;
  emxArray_real_T *alpha;
  emxArray_real_T *b_pow;
  emxArray_real_T *b_w;
  emxArray_real_T *p;
  emxArray_real_T *ss;
  emxArray_real_T *v;
  emxArray_real_T *w;
  creal_T dc;
  double b_j[2];
  double b_varargin_1_tmp[2];
  double varargin_1_tmp[2];
  double b_I;
  double curr_pos;
  double idx;
  double *RR_data;
  double *R_data;
  double *S_data;
  double *alpha_data;
  double *p_data;
  double *pow_data;
  double *ss_data;
  double *v_data;
  double *w_data;
  int b_i;
  int b_loop_ub;
  int c_loop_ub;
  int cen;
  int cen2;
  int d_loop_ub;
  int exitg1;
  int i;
  int i1;
  int ii;
  int j;
  int l;
  int loop_ub;
  int unnamed_idx_0_tmp_tmp;
  int varargin_2;
  int vectorUB;
  signed char fileid;
  bool exitg2;
  bool y;

  /* sfn: fasta file */
  /* wfn: kmer weight file */
  /* ofn: output prefix */
  /* 'mapTF:414' fid = fopen(wfn, 'r'); */
  fileid = cfopen(wfn, "rb");

  /* 'mapTF:415' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF:416' fprintf("ERROR: Weight file cannot be opened.\n") */
    printf("ERROR: Weight file cannot be opened.\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:419' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF:420' idx=0; */
  idx = 0.0;

  /* 'mapTF:421' while ~feof(fid) */
  emxInit_char_T(&b_fileid, 2);
  do {
    exitg1 = 0;
    b_I = b_feof(fileid);
    if (!(b_I != 0.0)) {
      /* 'mapTF:422' idx=idx+1; */
      idx++;

      /* 'mapTF:423' fgetl(fid); */
      b_fgets(fileid, b_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_fileid);
  emxInit_cell_wrap_3(&sequences);

  /* 'mapTF:425' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF:426' sequences = cell(idx, 1); */
  cen = (int)idx;
  i = sequences->size[0];
  sequences->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_3(sequences, i);
  sequences_data = sequences->data;
  for (i = 0; i < cen; i++) {
    sequences_data[i].f1->size[0] = 1;
    sequences_data[i].f1->size[1] = 0;
  }

  emxInit_real_T(&alpha, 1);

  /* 'mapTF:427' sequences = coder.nullcopy(sequences); */
  /* 'mapTF:428' alpha = zeros(idx, 1); */
  i = alpha->size[0];
  alpha->size[0] = (int)idx;
  emxEnsureCapacity_real_T(alpha, i);
  alpha_data = alpha->data;
  for (i = 0; i < cen; i++) {
    alpha_data[i] = 0.0;
  }

  /* 'mapTF:429' for cur_idx=1:idx */
  cen = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&cur_seq, 2);
  emxInit_char_T(&cur_alpha, 2);
  exitg2 = false;
  while ((!exitg2) && (cen <= (int)idx - 1)) {
    /* 'mapTF:430' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF:431' if cur_line == -1 */
    y = (cur_line->size[1] != 0);
    if (y) {
      y = (0 > cur_line->size[1] - 1);
    }

    if (y) {
      exitg2 = true;
    } else {
      /* 'mapTF:434' [cur_seq, cur_alpha] = strtok(cur_line, char(9)); */
      c_strtok(cur_line, cur_seq, cur_alpha);

      /* 'mapTF:435' alpha(cur_idx,1) = real(str2double(cur_alpha)); */
      dc = str2double(cur_alpha);
      alpha_data[cen] = dc.re;

      /* 'mapTF:436' sequences{cur_idx} = (strip(cur_seq)); */
      strip(cur_seq, sequences_data[cen].f1);
      cen++;
    }
  }

  emxFree_char_T(&cur_alpha);
  emxFree_char_T(&cur_seq);

  /* 'mapTF:438' fclose(fid); */
  cfclose(fileid);

  /* 'mapTF:440' l = length(sequences{1}); */
  varargin_2 = sequences_data[0].f1->size[1];
  l = sequences_data[0].f1->size[1];

  /* 'mapTF:441' if l ~= l_svm */
  if (sequences_data[0].f1->size[1] != l_svm) {
    /* 'mapTF:442' fprintf("ERROR: L must be the same as the length of k-mer in the weight file\n"); */
    printf("ERROR: L must be the same as the length of k-mer in the weight file\n");
    fflush(stdout);
    exit(1);
  }

  emxInit_real_T(&w, 1);

  /* 'mapTF:445' w = zeros(4^l,1); */
  loop_ub = (int)rt_powd_snf(4.0, sequences_data[0].f1->size[1]);
  i = w->size[0];
  w->size[0] = (int)rt_powd_snf(4.0, sequences_data[0].f1->size[1]);
  emxEnsureCapacity_real_T(w, i);
  w_data = w->data;
  for (i = 0; i < loop_ub; i++) {
    w_data[i] = 0.0;
  }

  /* 'mapTF:446' pow = (4.^(0:(l-1)))'; */
  emxInit_real_T(&R, 2);
  if (sequences_data[0].f1->size[1] - 1 < 0) {
    R->size[1] = 0;
  } else {
    i = R->size[0] * R->size[1];
    R->size[0] = 1;
    R->size[1] = sequences_data[0].f1->size[1];
    emxEnsureCapacity_real_T(R, i);
    R_data = R->data;
    loop_ub = sequences_data[0].f1->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      R_data[i] = i;
    }
  }

  i = R->size[0] * R->size[1];
  R->size[0] = 1;
  emxEnsureCapacity_real_T(R, i);
  R_data = R->data;
  loop_ub = R->size[1] - 1;
  for (i = 0; i <= loop_ub; i++) {
    curr_pos = R_data[i];
    R_data[i] = rt_powd_snf(4.0, curr_pos);
  }

  emxInit_real_T(&b_pow, 1);
  i = b_pow->size[0];
  b_pow->size[0] = R->size[1];
  emxEnsureCapacity_real_T(b_pow, i);
  pow_data = b_pow->data;
  loop_ub = R->size[1];
  for (i = 0; i < loop_ub; i++) {
    pow_data[i] = R_data[i];
  }

  /* 'mapTF:447' fprintf('calculating indices\n'); */
  printf("calculating indices\n");
  fflush(stdout);

  /* 'mapTF:448' for i = 1:numel(alpha) */
  i = alpha->size[0];
  emxInit_real_T(&ss, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:449' ss = letterconvert(sequences{i}); */
    letterconvert(sequences_data[b_i].f1, ss);
    ss_data = ss->data;

    /* 'mapTF:450' rs = 3-fliplr(ss); */
    /* 'mapTF:451' w(ss*pow+1) = alpha(i); */
    if (ss->size[1] < 1) {
      curr_pos = 0.0;
    } else {
      curr_pos = cblas_ddot((blasint)ss->size[1], &ss_data[0], (blasint)1,
                            &pow_data[0], (blasint)1);
    }

    w_data[(int)(curr_pos + 1.0) - 1] = alpha_data[b_i];

    /* 'mapTF:452' w(rs*pow+1) = alpha(i); */
    fliplr(ss);
    i1 = ss->size[0] * ss->size[1];
    ss->size[0] = 1;
    emxEnsureCapacity_real_T(ss, i1);
    ss_data = ss->data;
    loop_ub = ss->size[1] - 1;
    cen = (ss->size[1] / 2) << 1;
    vectorUB = cen - 2;
    for (i1 = 0; i1 <= vectorUB; i1 += 2) {
      r = _mm_loadu_pd(&ss_data[i1]);
      _mm_storeu_pd(&ss_data[i1], _mm_sub_pd(_mm_set1_pd(3.0), r));
    }

    for (i1 = cen; i1 <= loop_ub; i1++) {
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

  emxInit_real_T(&b_w, 1);

  /* 'mapTF:454' m = mean(alpha); */
  curr_pos = blockedSummation(alpha, alpha->size[0]) / (double)alpha->size[0];

  /* 'mapTF:455' s = std(alpha); */
  idx = b_std(alpha);

  /* 'mapTF:456' W = (1/2)*(1+erf((w-m)/s/sqrt(2))); */
  i = b_w->size[0];
  b_w->size[0] = w->size[0];
  emxEnsureCapacity_real_T(b_w, i);
  alpha_data = b_w->data;
  loop_ub = w->size[0];
  cen = (w->size[0] / 2) << 1;
  vectorUB = cen - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&w_data[i]);
    _mm_storeu_pd(&alpha_data[i], _mm_div_pd(_mm_div_pd(_mm_sub_pd(r,
      _mm_set1_pd(curr_pos)), _mm_set1_pd(idx)), _mm_set1_pd(1.4142135623730951)));
  }

  for (i = cen; i < loop_ub; i++) {
    alpha_data[i] = (w_data[i] - curr_pos) / idx / 1.4142135623730951;
  }

  applyScalarFunction(b_w, alpha);
  alpha_data = alpha->data;
  loop_ub = alpha->size[0];
  cen = (alpha->size[0] / 2) << 1;
  vectorUB = cen - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r = _mm_loadu_pd(&alpha_data[i]);
    _mm_storeu_pd(&alpha_data[i], _mm_mul_pd(_mm_set1_pd(0.5), _mm_add_pd(r,
      _mm_set1_pd(1.0))));
  }

  for (i = cen; i < loop_ub; i++) {
    alpha_data[i] = 0.5 * (alpha_data[i] + 1.0);
  }

  /* 'mapTF:457' W(W>0.99) = 0.99; */
  cen = alpha->size[0];
  for (b_i = 0; b_i < cen; b_i++) {
    if (alpha_data[b_i] > 0.99) {
      alpha_data[b_i] = 0.99;
    }
  }

  /* W(W<0.01) = 0.01; */
  /*  seq = importdata(sfn); */
  /* 'mapTF:461' fid = fopen(sfn, 'r'); */
  fileid = cfopen(sfn, "rb");

  /* 'mapTF:462' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF:463' fprintf("ERROR: Sequence file (.fa or .fasta) cannot be opened.\n") */
    printf("ERROR: Sequence file (.fa or .fasta) cannot be opened.\n");
    fflush(stdout);
    exit(1);
  }

  /* 'mapTF:466' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF:467' idx=0; */
  idx = 0.0;

  /* 'mapTF:468' while ~feof(fid) */
  emxInit_char_T(&c_fileid, 2);
  do {
    exitg1 = 0;
    b_I = b_feof(fileid);
    if (!(b_I != 0.0)) {
      /* 'mapTF:469' idx=idx+1; */
      idx++;

      /* 'mapTF:470' fgetl(fid); */
      b_fgets(fileid, c_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&c_fileid);

  /* 'mapTF:472' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF:473' seq = cell(idx, 1); */
  cen = (int)idx;
  i = seq->size[0];
  seq->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_3(seq, i);
  seq_data = seq->data;
  for (i = 0; i < cen; i++) {
    seq_data[i].f1->size[0] = 1;
    seq_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:474' seq = coder.nullcopy(seq); */
  /* 'mapTF:475' for cur_idx=1:idx */
  cen = 0;
  exitg2 = false;
  while ((!exitg2) && (cen <= (int)idx - 1)) {
    /* 'mapTF:476' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF:477' if cur_line == -1 */
    y = (cur_line->size[1] != 0);
    if (y) {
      y = (0 > cur_line->size[1] - 1);
    }

    if (y) {
      exitg2 = true;
    } else {
      /* 'mapTF:480' seq{cur_idx} = (strip(cur_line)); */
      strip(cur_line, seq_data[cen].f1);
      cen++;
    }
  }

  emxFree_char_T(&cur_line);

  /* 'mapTF:482' fclose(fid); */
  cfclose(fileid);

  /* 'mapTF:484' n = length(seq)/2; */
  curr_pos = (double)seq->size[0] / 2.0;

  /* 'mapTF:485' seqout = cell(n,1); */
  /* 'mapTF:486' seqindmat = cell(n,1); */
  unnamed_idx_0_tmp_tmp = (int)curr_pos;
  i = seqindmat->size[0];
  seqindmat->size[0] = (int)curr_pos;
  emxEnsureCapacity_cell_wrap_0(seqindmat, i);
  seqindmat_data = seqindmat->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    seqindmat_data[i].f1->size[0] = 0;
  }

  /* 'mapTF:487' seqindmat = coder.nullcopy(seqindmat); */
  /* 'mapTF:488' fprintf('converting kmers to probabilities\n'); */
  printf("converting kmers to probabilities\n");
  fflush(stdout);

  /* 'mapTF:489' P = cell(n,1); */
  /* 'mapTF:490' for i = 1:n */
  i = seqout->size[0];
  seqout->size[0] = (int)curr_pos;
  emxEnsureCapacity_cell_wrap_2(seqout, i);
  seqout_data = seqout->data;
  i = P->size[0];
  P->size[0] = (int)curr_pos;
  emxEnsureCapacity_cell_wrap_0(P, i);
  P_data = P->data;
  if (0 <= (int)curr_pos - 1) {
    if (1 > sequences_data[0].f1->size[1]) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = sequences_data[0].f1->size[1];
    }
  }

  emxInit_real_T(&p, 1);
  for (b_i = 0; b_i < unnamed_idx_0_tmp_tmp; b_i++) {
    /* 'mapTF:491' if mod(i,1000)==0 */
    if (fmod((double)b_i + 1.0, 1000.0) == 0.0) {
      /* 'mapTF:492' fprintf('%d sequences converted\n', int32(i)); */
      printf("%d sequences converted\n", b_i + 1);
      fflush(stdout);
    }

    /* 'mapTF:494' L = length(seq{2*i})-l+1; */
    cen = ((b_i + 1) << 1) - 1;
    idx = (double)(seq_data[cen].f1->size[1] - varargin_2) + 1.0;

    /* 'mapTF:495' seqindmat{i} = zeros(L,1); */
    cen2 = (int)idx;
    i = seqindmat_data[b_i].f1->size[0];
    seqindmat_data[b_i].f1->size[0] = (int)idx;
    emxEnsureCapacity_real_T(seqindmat_data[b_i].f1, i);
    for (i = 0; i < cen2; i++) {
      seqindmat_data[b_i].f1->data[i] = 0.0;
    }

    /* 'mapTF:496' ss = letterconvert(seq{2*i}); */
    letterconvert(seq_data[cen].f1, ss);
    ss_data = ss->data;

    /* 'mapTF:497' seqout{i} = ss+1; */
    i = seqout_data[b_i].f1->size[0] * seqout_data[b_i].f1->size[1];
    seqout_data[b_i].f1->size[0] = 1;
    seqout_data[b_i].f1->size[1] = ss->size[1];
    emxEnsureCapacity_real_T(seqout_data[b_i].f1, i);
    loop_ub = ss->size[1];
    cen = (ss->size[1] / 2) << 1;
    vectorUB = cen - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      r = _mm_loadu_pd(&ss_data[i]);
      _mm_storeu_pd(&seqout_data[b_i].f1->data[i], _mm_add_pd(r, _mm_set1_pd(1.0)));
    }

    for (i = cen; i < loop_ub; i++) {
      seqout_data[b_i].f1->data[i] = ss_data[i] + 1.0;
    }

    /* 'mapTF:498' p = zeros(L,1); */
    i = p->size[0];
    p->size[0] = (int)idx;
    emxEnsureCapacity_real_T(p, i);
    p_data = p->data;
    for (i = 0; i < cen2; i++) {
      p_data[i] = 0.0;
    }

    /* 'mapTF:499' I = ss(1:l)*pow; */
    i = R->size[0] * R->size[1];
    R->size[0] = 1;
    R->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(R, i);
    R_data = R->data;
    for (i = 0; i < b_loop_ub; i++) {
      R_data[i] = ss_data[i];
    }

    if (b_loop_ub < 1) {
      b_I = 0.0;
    } else {
      b_I = cblas_ddot((blasint)b_loop_ub, &R_data[0], (blasint)1, &pow_data[0],
                       (blasint)1);
    }

    /* 'mapTF:500' seqindmat{i}(1) = I+1; */
    seqindmat_data[b_i].f1->data[0] = b_I + 1.0;

    /* 'mapTF:501' p(1) = W(I+1); */
    p_data[0] = alpha_data[(int)(b_I + 1.0) - 1];

    /* 'mapTF:502' for j = 2:L */
    i = (int)(idx + -1.0);
    for (j = 0; j < i; j++) {
      /* 'mapTF:503' I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1); */
      b_I = (b_I - ss_data[j]) / 4.0 + rt_powd_snf(4.0, (double)varargin_2 - 1.0)
        * ss_data[(int)((unsigned int)j + varargin_2)];

      /* 'mapTF:504' seqindmat{i}(j) = I+1; */
      seqindmat_data[b_i].f1->data[j + 1] = b_I + 1.0;

      /* 'mapTF:505' p(j) = W(I+1); */
      p_data[j + 1] = alpha_data[(int)(b_I + 1.0) - 1];
    }

    /* 'mapTF:507' P{i} = p; */
    i = P_data[b_i].f1->size[0];
    P_data[b_i].f1->size[0] = p->size[0];
    emxEnsureCapacity_real_T(P_data[b_i].f1, i);
    loop_ub = p->size[0];
    for (i = 0; i < loop_ub; i++) {
      P_data[b_i].f1->data[i] = p_data[i];
    }
  }

  emxFree_real_T(&alpha);

  /* 'mapTF:510' fprintf('Running dsvm\n'); */
  printf("Running dsvm\n");
  fflush(stdout);

  /* 'mapTF:511' mat = [1 2 3;0 2 3;0 1 3;0 1 2]; */
  /* 'mapTF:512' O = ones(1,l); */
  /* 'mapTF:513' V = cell(n,1); */
  /* 'mapTF:514' for i = 1:n */
  i = V->size[0];
  V->size[0] = (int)curr_pos;
  emxEnsureCapacity_cell_wrap_1(V, i);
  V_data = V->data;
  if (0 <= (int)curr_pos - 1) {
    if (1 > sequences_data[0].f1->size[1]) {
      c_loop_ub = 0;
    } else {
      c_loop_ub = sequences_data[0].f1->size[1];
    }
  }

  emxFree_cell_wrap_3(&sequences);
  emxInit_real_T(&v, 2);
  emxInit_real_T(&RR, 2);
  RR_data = RR->data;
  emxInit_real_T(&S, 2);
  for (b_i = 0; b_i < unnamed_idx_0_tmp_tmp; b_i++) {
    /* 'mapTF:515' if mod(i,1000)==0 */
    if (fmod((double)b_i + 1.0, 1000.0) == 0.0) {
      /* 'mapTF:516' fprintf('%d sequences converted\n', int32(i)); */
      printf("%d sequences converted\n", b_i + 1);
      fflush(stdout);
    }

    /* 'mapTF:518' L = length(seq{2*i}); */
    cen2 = ((b_i + 1) << 1) - 1;
    b_loop_ub = seq_data[cen2].f1->size[1];

    /* 'mapTF:519' ss = seqout{i}-1; */
    i = ss->size[0] * ss->size[1];
    ss->size[0] = 1;
    ss->size[1] = seqout_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(ss, i);
    ss_data = ss->data;
    loop_ub = seqout_data[b_i].f1->size[1];
    cen = (loop_ub / 2) << 1;
    vectorUB = cen - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      r = _mm_loadu_pd(&seqout_data[b_i].f1->data[i]);
      _mm_storeu_pd(&ss_data[i], _mm_sub_pd(r, _mm_set1_pd(1.0)));
    }

    for (i = cen; i < loop_ub; i++) {
      ss_data[i] = seqout_data[b_i].f1->data[i] - 1.0;
    }

    /* 'mapTF:520' p = zeros(L+l-1,1); */
    loop_ub = (int)((double)((unsigned int)seq_data[cen2].f1->size[1] +
      varargin_2) - 1.0);
    i = p->size[0];
    p->size[0] = loop_ub;
    emxEnsureCapacity_real_T(p, i);
    p_data = p->data;
    for (i = 0; i < loop_ub; i++) {
      p_data[i] = 0.0;
    }

    /* 'mapTF:521' I = ss(1:l)*pow; */
    i = R->size[0] * R->size[1];
    R->size[0] = 1;
    R->size[1] = c_loop_ub;
    emxEnsureCapacity_real_T(R, i);
    R_data = R->data;
    for (i = 0; i < c_loop_ub; i++) {
      R_data[i] = ss_data[i];
    }

    vectorUB = c_loop_ub;
    if (c_loop_ub < 1) {
      b_I = 0.0;
    } else {
      b_I = cblas_ddot((blasint)c_loop_ub, &R_data[0], (blasint)1, &pow_data[0],
                       (blasint)1);
    }

    /* 'mapTF:522' p(1) = w(I+1); */
    p_data[0] = w_data[(int)(b_I + 1.0) - 1];

    /* 'mapTF:523' for j = 2:L-l+1 */
    i = seq_data[cen2].f1->size[1] - l;
    for (j = 0; j < i; j++) {
      /* 'mapTF:524' I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1); */
      b_I = (b_I - ss_data[j]) / 4.0 + rt_powd_snf(4.0, (double)varargin_2 - 1.0)
        * ss_data[(int)((unsigned int)j + varargin_2)];

      /* 'mapTF:525' p(j) = w(I+1); */
      p_data[j + 1] = w_data[(int)(b_I + 1.0) - 1];
    }

    /* 'mapTF:527' v = zeros(L,3); */
    i = seq_data[cen2].f1->size[1];
    i1 = v->size[0] * v->size[1];
    v->size[0] = i;
    v->size[1] = 3;
    emxEnsureCapacity_real_T(v, i1);
    v_data = v->data;
    loop_ub = i * 3;
    for (i = 0; i < loop_ub; i++) {
      v_data[i] = 0.0;
    }

    /* 'mapTF:528' for j = 1:L */
    i = seq_data[cen2].f1->size[1];
    if (0 <= i - 1) {
      varargin_1_tmp[0] = 1.0;
      b_varargin_1_tmp[1] = (unsigned int)b_loop_ub;
      b_j[1] = ((double)b_loop_ub + (double)l) - 1.0;
      if (1 > l) {
        d_loop_ub = 0;
      } else {
        d_loop_ub = varargin_2;
      }

      vectorUB = d_loop_ub;
    }

    for (j = 0; j < i; j++) {
      /* 'mapTF:529' R = max([1 j-l+1]):min([j+l-1 L]); */
      varargin_1_tmp[1] = (double)((j - varargin_2) + 1) + 1.0;
      idx = b_maximum(varargin_1_tmp);
      b_varargin_1_tmp[0] = (unsigned int)j + varargin_2;
      b_I = minimum(b_varargin_1_tmp);
      if (rtIsNaN(idx) || rtIsNaN(b_I)) {
        i1 = R->size[0] * R->size[1];
        R->size[0] = 1;
        R->size[1] = 1;
        emxEnsureCapacity_real_T(R, i1);
        R_data = R->data;
        R_data[0] = rtNaN;
      } else if (b_I < idx) {
        R->size[0] = 1;
        R->size[1] = 0;
      } else if ((rtIsInf(idx) || rtIsInf(b_I)) && (idx == b_I)) {
        i1 = R->size[0] * R->size[1];
        R->size[0] = 1;
        R->size[1] = 1;
        emxEnsureCapacity_real_T(R, i1);
        R_data = R->data;
        R_data[0] = rtNaN;
      } else if (floor(idx) == idx) {
        i1 = R->size[0] * R->size[1];
        R->size[0] = 1;
        loop_ub = (int)floor(b_I - idx);
        R->size[1] = loop_ub + 1;
        emxEnsureCapacity_real_T(R, i1);
        R_data = R->data;
        for (i1 = 0; i1 <= loop_ub; i1++) {
          R_data[i1] = idx + (double)i1;
        }
      } else {
        eml_float_colon(idx, b_I, R);
        R_data = R->data;
      }

      /* 'mapTF:530' RR = max([1 j-l+1]):min([j L+l-1]); */
      b_j[0] = (double)j + 1.0;
      curr_pos = minimum(b_j);
      if (rtIsNaN(idx) || rtIsNaN(curr_pos)) {
        i1 = RR->size[0] * RR->size[1];
        RR->size[0] = 1;
        RR->size[1] = 1;
        emxEnsureCapacity_real_T(RR, i1);
        RR_data = RR->data;
        RR_data[0] = rtNaN;
      } else if (curr_pos < idx) {
        RR->size[0] = 1;
        RR->size[1] = 0;
      } else if ((rtIsInf(idx) || rtIsInf(curr_pos)) && (idx == curr_pos)) {
        i1 = RR->size[0] * RR->size[1];
        RR->size[0] = 1;
        RR->size[1] = 1;
        emxEnsureCapacity_real_T(RR, i1);
        RR_data = RR->data;
        RR_data[0] = rtNaN;
      } else if (floor(idx) == idx) {
        i1 = RR->size[0] * RR->size[1];
        RR->size[0] = 1;
        loop_ub = (int)floor(curr_pos - idx);
        RR->size[1] = loop_ub + 1;
        emxEnsureCapacity_real_T(RR, i1);
        RR_data = RR->data;
        for (i1 = 0; i1 <= loop_ub; i1++) {
          RR_data[i1] = idx + (double)i1;
        }
      } else {
        eml_float_colon(idx, curr_pos, RR);
        RR_data = RR->data;
      }

      /* 'mapTF:531' S = ss(R); */
      i1 = S->size[0] * S->size[1];
      S->size[0] = 1;
      S->size[1] = R->size[1];
      emxEnsureCapacity_real_T(S, i1);
      S_data = S->data;
      loop_ub = R->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        S_data[i1] = ss_data[(int)R_data[i1] - 1];
      }

      /* 'mapTF:532' ref = sum(p(RR)); */
      i1 = b_w->size[0];
      b_w->size[0] = RR->size[1];
      emxEnsureCapacity_real_T(b_w, i1);
      alpha_data = b_w->data;
      loop_ub = RR->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        alpha_data[i1] = p_data[(int)RR_data[i1] - 1];
      }

      curr_pos = blockedSummation(b_w, RR->size[1]);

      /* 'mapTF:533' if max([1 j-l+1]) == 1 */
      if (idx == 1.0) {
        /* 'mapTF:534' cen = j; */
        cen = j + 1;

        /* 'mapTF:535' cen2 = j; */
        cen2 = j + 1;
      } else if (b_I == b_loop_ub) {
        /* 'mapTF:536' elseif min([j+l-1 L]) == L */
        /* 'mapTF:537' cen = l; */
        cen = varargin_2;

        /* 'mapTF:538' cen2 = L-j+1; */
        cen2 = b_loop_ub - j;
      } else {
        /* 'mapTF:539' else */
        /* 'mapTF:540' cen = l; */
        cen = varargin_2;

        /* 'mapTF:541' cen2 = l; */
        cen2 = varargin_2;
      }

      /* 'mapTF:543' for ii = 1:3 */
      for (ii = 0; ii < 3; ii++) {
        /* 'mapTF:544' S(cen) = mat(ss(j)+1,ii); */
        S_data[cen - 1] = mat[((int)(ss_data[j] + 1.0) + (ii << 2)) - 1];

        /* 'mapTF:545' I = S(1:l)*pow; */
        i1 = R->size[0] * R->size[1];
        R->size[0] = 1;
        R->size[1] = d_loop_ub;
        emxEnsureCapacity_real_T(R, i1);
        R_data = R->data;
        for (i1 = 0; i1 < d_loop_ub; i1++) {
          R_data[i1] = S_data[i1];
        }

        if (vectorUB < 1) {
          b_I = 0.0;
        } else {
          b_I = cblas_ddot((blasint)d_loop_ub, &R_data[0], (blasint)1,
                           &pow_data[0], (blasint)1);
        }

        /* 'mapTF:546' v(j,ii) = w(I+1); */
        v_data[j + v->size[0] * ii] = w_data[(int)(b_I + 1.0) - 1];

        /* 'mapTF:547' if length(RR) > 1 */
        if (RR->size[1] > 1) {
          /* 'mapTF:548' for jj = 2:cen2 */
          for (loop_ub = 0; loop_ub <= cen2 - 2; loop_ub++) {
            /* 'mapTF:549' I = (I-S(jj-1))/4+4^(l-1)*S(jj+l-1); */
            b_I = (b_I - S_data[loop_ub]) / 4.0 + rt_powd_snf(4.0, (double)
              varargin_2 - 1.0) * S_data[(int)((unsigned int)loop_ub +
              varargin_2)];

            /* 'mapTF:550' v(j,ii)=v(j,ii)+w(I+1); */
            v_data[j + v->size[0] * ii] += w_data[(int)(b_I + 1.0) - 1];
          }
        }

        /* 'mapTF:553' v(j,ii) = v(j,ii)-ref; */
        v_data[j + v->size[0] * ii] -= curr_pos;
      }
    }

    /* 'mapTF:556' V{i} = v; */
    i = V_data[b_i].f1->size[0] * V_data[b_i].f1->size[1];
    V_data[b_i].f1->size[0] = v->size[0];
    V_data[b_i].f1->size[1] = 3;
    emxEnsureCapacity_real_T(V_data[b_i].f1, i);
    loop_ub = v->size[0] * 3;
    for (i = 0; i < loop_ub; i++) {
      V_data[b_i].f1->data[i] = v_data[i];
    }
  }

  emxFree_real_T(&b_w);
  emxFree_real_T(&S);
  emxFree_real_T(&RR);
  emxFree_real_T(&R);
  emxFree_real_T(&v);
  emxFree_real_T(&p);
  emxFree_real_T(&ss);
  emxFree_real_T(&b_pow);
  emxFree_real_T(&w);
}

/*
 * function [pp, info, len] = trim_pwm(p,cut)
 */
static void trim_pwm(emxArray_cell_wrap_5 *p, emxArray_real_T *info,
                     emxArray_real_T *len)
{
  __m128d r2;
  __m128d r3;
  cell_wrap_5 *p_data;
  emxArray_real_T *b_mat;
  emxArray_real_T *mat;
  emxArray_real_T *r;
  emxArray_real_T *vec;
  double *b_mat_data;
  double *info_data;
  double *len_data;
  double *mat_data;
  double *r1;
  int b_i;
  int c_i;
  int i;
  int idx;
  int j;
  int nrows;
  int nx;
  int nxin;
  p_data = p->data;

  /* 'mapTF:911' l = length(p); */
  /* 'mapTF:912' info = zeros(l, 1); */
  nx = p->size[0];
  i = info->size[0];
  info->size[0] = nx;
  emxEnsureCapacity_real_T(info, i);
  info_data = info->data;

  /* 'mapTF:913' len = zeros(l,1); */
  i = len->size[0];
  len->size[0] = nx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;

  /* 'mapTF:914' for i = 1:l */
  i = p->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&vec, 1);
  emxInit_real_T(&r, 2);
  emxInit_real_T(&b_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:915' mat = p{i}+(p{i}==0); */
    idx = mat->size[0] * mat->size[1];
    mat->size[0] = p_data[b_i].f1->size[0];
    mat->size[1] = p_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(mat, idx);
    mat_data = mat->data;
    nrows = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
    for (idx = 0; idx < nrows; idx++) {
      mat_data[idx] = p_data[b_i].f1->data[idx] + (double)(p_data[b_i].f1->
        data[idx] == 0.0);
    }

    /* 'mapTF:916' vec = 2+sum(mat.*log(mat)/log(2),2); */
    idx = r->size[0] * r->size[1];
    r->size[0] = mat->size[0];
    r->size[1] = mat->size[1];
    emxEnsureCapacity_real_T(r, idx);
    r1 = r->data;
    nrows = mat->size[0] * mat->size[1];
    for (idx = 0; idx < nrows; idx++) {
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
      nrows = mat->size[0] * mat->size[1];
      nx = (nrows / 2) << 1;
      nxin = nx - 2;
      for (idx = 0; idx <= nxin; idx += 2) {
        r2 = _mm_loadu_pd(&mat_data[idx]);
        r3 = _mm_loadu_pd(&r1[idx]);
        _mm_storeu_pd(&b_mat_data[idx], _mm_div_pd(_mm_mul_pd(r2, r3),
          _mm_set1_pd(0.69314718055994529)));
      }

      for (idx = nx; idx < nrows; idx++) {
        b_mat_data[idx] = mat_data[idx] * r1[idx] / 0.69314718055994529;
      }

      sum(b_mat, vec);
      mat_data = vec->data;
    } else {
      i_binary_expand_op(vec, mat, r);
      mat_data = vec->data;
    }

    nrows = vec->size[0];
    nx = (vec->size[0] / 2) << 1;
    nxin = nx - 2;
    for (idx = 0; idx <= nxin; idx += 2) {
      r2 = _mm_loadu_pd(&mat_data[idx]);
      _mm_storeu_pd(&mat_data[idx], _mm_add_pd(r2, _mm_set1_pd(2.0)));
    }

    for (idx = nx; idx < nrows; idx++) {
      mat_data[idx] += 2.0;
    }

    /* 'mapTF:917' while (vec(1) < cut || mean(vec(1:3)) < cut || mean(vec(2:4)) < cut) && length(vec) > 4 */
    while (((mat_data[0] < 0.25) || (((mat_data[0] + mat_data[1]) + mat_data[2])
             / 3.0 < 0.25) || (((mat_data[1] + mat_data[2]) + mat_data[3]) / 3.0
             < 0.25)) && (vec->size[0] > 4)) {
      /* 'mapTF:918' p{i}(1,:) = []; */
      nx = p_data[b_i].f1->size[0] - 2;
      nxin = p_data[b_i].f1->size[1];
      nrows = p_data[b_i].f1->size[0] - 1;
      for (j = 0; j < nxin; j++) {
        for (c_i = 0; c_i < nrows; c_i++) {
          p_data[b_i].f1->data[c_i + p_data[b_i].f1->size[0] * j] = p_data[b_i].
            f1->data[(c_i + p_data[b_i].f1->size[0] * j) + 1];
        }
      }

      if (1 > nrows) {
        nrows = 0;
      } else {
        nrows = nx + 1;
      }

      nx = p_data[b_i].f1->size[1] - 1;
      for (idx = 0; idx <= nx; idx++) {
        for (nxin = 0; nxin < nrows; nxin++) {
          p_data[b_i].f1->data[nxin + nrows * idx] = p_data[b_i].f1->data[nxin +
            p_data[b_i].f1->size[0] * idx];
        }
      }

      idx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
      p_data[b_i].f1->size[0] = nrows;
      p_data[b_i].f1->size[1] = nx + 1;
      emxEnsureCapacity_real_T(p_data[b_i].f1, idx);

      /* 'mapTF:919' vec(1) = []; */
      nxin = vec->size[0];
      nx = vec->size[0];
      for (nrows = 0; nrows <= nx - 2; nrows++) {
        mat_data[nrows] = mat_data[nrows + 1];
      }

      idx = vec->size[0];
      vec->size[0] = nxin - 1;
      emxEnsureCapacity_real_T(vec, idx);
      mat_data = vec->data;
    }

    /* 'mapTF:921' while (vec(end) < cut || mean(vec(end-2:end)) < cut || mean(vec(end-3:end-1)) < cut) && length(vec) > 4 */
    while (((mat_data[vec->size[0] - 1] < 0.25) || (((mat_data[vec->size[0] - 3]
               + mat_data[vec->size[0] - 2]) + mat_data[vec->size[0] - 1]) / 3.0
             < 0.25) || (((mat_data[vec->size[0] - 4] + mat_data[vec->size[0] -
                           3]) + mat_data[vec->size[0] - 2]) / 3.0 < 0.25)) &&
           (vec->size[0] > 4)) {
      /* 'mapTF:922' vec(end) = []; */
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

      /* 'mapTF:923' p{i}(end,:) = []; */
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
        nrows = 0;
      } else {
        nrows = nx + 1;
      }

      nx = p_data[b_i].f1->size[1] - 1;
      for (idx = 0; idx <= nx; idx++) {
        for (nxin = 0; nxin < nrows; nxin++) {
          p_data[b_i].f1->data[nxin + nrows * idx] = p_data[b_i].f1->data[nxin +
            p_data[b_i].f1->size[0] * idx];
        }
      }

      idx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
      p_data[b_i].f1->size[0] = nrows;
      p_data[b_i].f1->size[1] = nx + 1;
      emxEnsureCapacity_real_T(p_data[b_i].f1, idx);
    }

    /* 'mapTF:925' info(i) = sum(vec); */
    info_data[b_i] = blockedSummation(vec, vec->size[0]);

    /* 'mapTF:926' [len(i), ~] = size(p{i}); */
    len_data[b_i] = p_data[b_i].f1->size[0];
  }

  emxFree_real_T(&b_mat);
  emxFree_real_T(&r);
  emxFree_real_T(&vec);
  emxFree_real_T(&mat);

  /* 'mapTF:928' pp = p; */
}

/*
 * function mapTF(varargin)
 */
void mapTF(const emxArray_char_T *varargin_1, const emxArray_char_T *varargin_2,
           const emxArray_char_T *varargin_3, const emxArray_char_T *varargin_4,
           const emxArray_char_T *varargin_5, const emxArray_char_T *varargin_6,
           double varargin_7, double varargin_8, double varargin_9)
{
  __m128d r2;
  __m128d r5;
  __m128d r6;
  cell_wrap_0 *P_data;
  cell_wrap_0 *VV_data;
  cell_wrap_0 *seqindmat_data;
  cell_wrap_1 *V_data;
  cell_wrap_2 *ss_data;
  cell_wrap_4 *PWM2_data;
  cell_wrap_4 *PWM3_data;
  cell_wrap_4 *PWM_data;
  cell_wrap_4 *lPWM3_data;
  cell_wrap_4 *pwm_data;
  cell_wrap_5 *Smat_data;
  cell_wrap_6 *NN_data;
  emxArray_boolean_T *b_f;
  emxArray_cell_wrap_0 *P;
  emxArray_cell_wrap_0 *VV;
  emxArray_cell_wrap_0 *seqindmat;
  emxArray_cell_wrap_1 *V;
  emxArray_cell_wrap_2 *ss;
  emxArray_cell_wrap_3 *names;
  emxArray_cell_wrap_3 *seq;
  emxArray_cell_wrap_4 *LL;
  emxArray_cell_wrap_4 *PWM;
  emxArray_cell_wrap_4 *PWM2;
  emxArray_cell_wrap_4 *PWM3;
  emxArray_cell_wrap_4 *lPWM3;
  emxArray_cell_wrap_4 *lpwm;
  emxArray_cell_wrap_4 *p;
  emxArray_cell_wrap_4 *pwm;
  emxArray_cell_wrap_5 *Smat;
  emxArray_cell_wrap_6 *NN;
  emxArray_char_T *b_varargin_6;
  emxArray_int32_T *r3;
  emxArray_int8_T *IND;
  emxArray_int8_T *ss_onehot;
  emxArray_real_T *B;
  emxArray_real_T *LEN;
  emxArray_real_T *LEN_2;
  emxArray_real_T *a__3;
  emxArray_real_T *all_pwm;
  emxArray_real_T *b_B;
  emxArray_real_T *b_all_pwm;
  emxArray_real_T *b_lpwm;
  emxArray_real_T *c;
  emxArray_real_T *f;
  emxArray_real_T *len;
  emxArray_real_T *maxnorm;
  emxArray_real_T *minnorm;
  emxArray_real_T *pwm_prob;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  emxArray_real_T *seqmat;
  emxArray_real_T *shift;
  emxArray_real_T *x;
  double GCmat[4];
  double b_varargin_7[2];
  double GC;
  double b;
  double *LEN_2_data;
  double *LEN_data;
  double *all_pwm_data;
  double *c_data;
  double *f_data;
  double *maxnorm_data;
  double *minnorm_data;
  double *pwm_prob_data;
  double *seqmat_data;
  double *shift_data;
  int b_I;
  int b_i;
  int b_loop_ub;
  int c_loop_ub;
  int d_loop_ub;
  int i;
  int i1;
  int i2;
  int i3;
  int j;
  int loop_ub;
  int nx;
  int sizes_idx_0;
  int unnamed_idx_0;
  int varargin_7_idx_0;
  int vectorUB;
  int *r4;
  const char *varargin_6_data;
  signed char *IND_data;
  char *b_varargin_6_data;
  signed char *ss_onehot_data;
  bool empty_non_axis_sizes;
  bool *b_f_data;
  if (!isInitialized_mapTF) {
    mapTF_initialize();
  }

  varargin_6_data = varargin_6->data;
  emxInit_cell_wrap_0(&P);
  emxInit_cell_wrap_1(&V);
  emxInit_cell_wrap_0(&seqindmat);
  emxInit_cell_wrap_2(&ss);
  emxInit_cell_wrap_3(&seq);
  emxInit_char_T(&b_varargin_6, 2);

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
  /* 'mapTF:48' if nargin < 6 */
  /* 'mapTF:52' fn = varargin{1}; */
  /* 'mapTF:53' wfn = varargin{2}; */
  /* 'mapTF:54' mfn1 = varargin{3}; */
  /* 'mapTF:55' mfn2 = varargin{4}; */
  /* 'mapTF:56' memefn = varargin{5}; */
  /* 'mapTF:57' ofn = varargin{6}; */
  /* 'mapTF:58' l_svm = varargin{7}; */
  /* 'mapTF:59' k_svm = varargin{8}; */
  /* 'mapTF:60' nfrac = varargin{9}; */
  /* 'mapTF:62' fprintf('Processing Motifs\n'); */
  printf("Processing Motifs\n");
  fflush(stdout);

  /* 'mapTF:63' process_motifs(mfn1, mfn2, memefn, ofn) */
  process_motifs(varargin_3, varargin_4, varargin_5, varargin_6);

  /* 'mapTF:64' mfn = [ofn '_motifs.out']; */
  /* 'mapTF:65' [P,V, seqindmat, ss, seq] = seq2pv(fn, wfn,l_svm); */
  seq2pv(varargin_1, varargin_2, varargin_7, P, V, seqindmat, ss, seq);
  ss_data = ss->data;
  seqindmat_data = seqindmat->data;
  V_data = V->data;
  P_data = P->data;

  /* 'mapTF:66' GC = countGC(ss); */
  GC = countGC(ss);

  /* 'mapTF:67' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  b = 0.5 - GC / 2.0;
  GCmat[0] = b;
  GCmat[1] = GC / 2.0;
  GCmat[2] = GC / 2.0;
  GCmat[3] = b;

  /* 'mapTF:68' b = 1; */
  b = 1.0;

  /* 'mapTF:69' [p,names,len] = getMOTIF(mfn); */
  i = b_varargin_6->size[0] * b_varargin_6->size[1];
  b_varargin_6->size[0] = 1;
  b_varargin_6->size[1] = varargin_6->size[1] + 11;
  emxEnsureCapacity_char_T(b_varargin_6, i);
  b_varargin_6_data = b_varargin_6->data;
  loop_ub = varargin_6->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_varargin_6_data[i] = varargin_6_data[i];
  }

  for (i = 0; i < 11; i++) {
    b_varargin_6_data[i + varargin_6->size[1]] = cv[i];
  }

  emxInit_cell_wrap_4(&PWM);
  emxInit_cell_wrap_4(&PWM2);
  emxInit_cell_wrap_4(&PWM3);
  emxInit_cell_wrap_4(&lPWM3);
  emxInit_cell_wrap_4(&p);
  emxInit_cell_wrap_3(&names);
  emxInit_real_T(&len, 1);
  getMOTIF(b_varargin_6, p, names, len);
  f_data = len->data;
  pwm_data = p->data;

  /* 'mapTF:70' a = numel(len); */
  /* 'mapTF:71' PWM = cell(a,1); */
  unnamed_idx_0 = len->size[0];
  i = PWM->size[0];
  PWM->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_4(PWM, i);
  PWM_data = PWM->data;

  /* 'mapTF:72' PWM = coder.nullcopy(PWM); */
  /* 'mapTF:73' PWM2 = cell(a,1); */
  i = PWM2->size[0];
  PWM2->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_4(PWM2, i);
  PWM2_data = PWM2->data;

  /* 'mapTF:74' PWM2 = coder.nullcopy(PWM2); */
  /* 'mapTF:75' PWM3 = cell(a,1); */
  i = PWM3->size[0];
  PWM3->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_4(PWM3, i);
  PWM3_data = PWM3->data;

  /* 'mapTF:76' PWM3 = coder.nullcopy(PWM3); */
  /* 'mapTF:77' lPWM3 = cell(a,1); */
  i = lPWM3->size[0];
  lPWM3->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_4(lPWM3, i);
  lPWM3_data = lPWM3->data;
  emxFree_char_T(&b_varargin_6);
  for (i = 0; i < unnamed_idx_0; i++) {
    PWM_data[i].f1->size[0] = 0;
    PWM_data[i].f1->size[1] = 4;
    PWM2_data[i].f1->size[0] = 0;
    PWM2_data[i].f1->size[1] = 4;
    PWM3_data[i].f1->size[0] = 0;
    PWM3_data[i].f1->size[1] = 4;
    lPWM3_data[i].f1->size[0] = 0;
    lPWM3_data[i].f1->size[1] = 4;
  }

  emxInit_real_T(&LEN, 1);

  /* 'mapTF:78' lPWM3 = coder.nullcopy(lPWM3); */
  /* 'mapTF:80' LEN = zeros(a,1); */
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

  /* 'mapTF:81' LEN_2 = zeros(a,1); */
  /* 'mapTF:82' shift = zeros(a,1); */
  /* 'mapTF:84' for i = 1:a */
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
    c_loop_ub = (int)(varargin_7 - 1.0);
  }

  emxInit_real_T(&r, 2);
  emxInit_real_T(&r1, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:85' shift(i) = max([l_svm-len(i) 4]); */
    b_varargin_7[0] = varargin_7 - f_data[b_i];
    GC = b_maximum(b_varargin_7);
    shift_data[b_i] = GC;

    /* 'mapTF:86' PWM{i} = [repmat(GCmat,shift(i), 1); p{i} ;repmat(GCmat,shift(i), 1)]; */
    repmat(GCmat, GC, r);
    all_pwm_data = r->data;
    repmat(GCmat, GC, r1);
    pwm_prob_data = r1->data;
    loop_ub = pwm_data[b_i].f1->size[0];
    i1 = PWM_data[b_i].f1->size[0] * PWM_data[b_i].f1->size[1];
    PWM_data[b_i].f1->size[0] = (r->size[0] + pwm_data[b_i].f1->size[0]) +
      r1->size[0];
    PWM_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(PWM_data[b_i].f1, i1);
    nx = r->size[0];
    unnamed_idx_0 = r1->size[0];
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < nx; i2++) {
        PWM_data[b_i].f1->data[i2 + PWM_data[b_i].f1->size[0] * i1] =
          all_pwm_data[i2 + r->size[0] * i1];
      }

      for (i2 = 0; i2 < loop_ub; i2++) {
        PWM_data[b_i].f1->data[(i2 + r->size[0]) + PWM_data[b_i].f1->size[0] *
          i1] = pwm_data[b_i].f1->data[i2 + pwm_data[b_i].f1->size[0] * i1];
      }

      for (i2 = 0; i2 < unnamed_idx_0; i2++) {
        PWM_data[b_i].f1->data[((i2 + r->size[0]) + pwm_data[b_i].f1->size[0]) +
          PWM_data[b_i].f1->size[0] * i1] = pwm_prob_data[i2 + r1->size[0] * i1];
      }
    }

    /* 'mapTF:87' PWM2{i} = p{i}; */
    i1 = PWM2_data[b_i].f1->size[0] * PWM2_data[b_i].f1->size[1];
    PWM2_data[b_i].f1->size[0] = pwm_data[b_i].f1->size[0];
    PWM2_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(PWM2_data[b_i].f1, i1);
    loop_ub = pwm_data[b_i].f1->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      PWM2_data[b_i].f1->data[i1] = pwm_data[b_i].f1->data[i1];
    }

    /* 'mapTF:88' PWM3{i} = [ones(l_svm-1,4)/4; p{i} ;ones(l_svm-1,4)/4]; */
    i1 = PWM3_data[b_i].f1->size[0] * PWM3_data[b_i].f1->size[1];
    PWM3_data[b_i].f1->size[0] = (varargin_7_idx_0 + pwm_data[b_i].f1->size[0])
      + sizes_idx_0;
    PWM3_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(PWM3_data[b_i].f1, i1);
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        PWM3_data[b_i].f1->data[i2 + PWM3_data[b_i].f1->size[0] * i1] = 0.25;
      }
    }

    loop_ub = pwm_data[b_i].f1->size[0];
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        PWM3_data[b_i].f1->data[(i2 + varargin_7_idx_0) + PWM3_data[b_i]
          .f1->size[0] * i1] = pwm_data[b_i].f1->data[i2 + pwm_data[b_i]
          .f1->size[0] * i1];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < c_loop_ub; i2++) {
        PWM3_data[b_i].f1->data[((i2 + varargin_7_idx_0) + pwm_data[b_i]
          .f1->size[0]) + PWM3_data[b_i].f1->size[0] * i1] = 0.25;
      }
    }

    /* 'mapTF:89' lPWM3{i} = log((PWM3{i}+10^-10)/(1+4*10^-10)); */
    i1 = lPWM3_data[b_i].f1->size[0] * lPWM3_data[b_i].f1->size[1];
    lPWM3_data[b_i].f1->size[0] = PWM3_data[b_i].f1->size[0];
    lPWM3_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(lPWM3_data[b_i].f1, i1);
    loop_ub = PWM3_data[b_i].f1->size[0] * 4;
    unnamed_idx_0 = (loop_ub / 2) << 1;
    vectorUB = unnamed_idx_0 - 2;
    for (i1 = 0; i1 <= vectorUB; i1 += 2) {
      r2 = _mm_loadu_pd(&PWM3_data[b_i].f1->data[i1]);
      _mm_storeu_pd(&lPWM3_data[b_i].f1->data[i1], _mm_div_pd(_mm_add_pd(r2,
        _mm_set1_pd(1.0E-10)), _mm_set1_pd(1.0000000004)));
    }

    for (i1 = unnamed_idx_0; i1 < loop_ub; i1++) {
      lPWM3_data[b_i].f1->data[i1] = (PWM3_data[b_i].f1->data[i1] + 1.0E-10) /
        1.0000000004;
    }

    nx = lPWM3_data[b_i].f1->size[0] << 2;
    for (unnamed_idx_0 = 0; unnamed_idx_0 < nx; unnamed_idx_0++) {
      lPWM3_data[b_i].f1->data[unnamed_idx_0] = log(lPWM3_data[b_i].f1->
        data[unnamed_idx_0]);
    }

    /* 'mapTF:90' LEN_2(i) = len(i); */
    LEN_2_data[b_i] = f_data[b_i];

    /* 'mapTF:91' LEN(i) = len(i)+2*shift(i)-l_svm+1; */
    GC = ((f_data[b_i] + 2.0 * shift_data[b_i]) - varargin_7) + 1.0;
    LEN_data[b_i] = GC;

    /* 'mapTF:92' for j = 1:LEN(i) */
    i1 = (int)GC;
    for (j = 0; j < i1; j++) {
      /* 'mapTF:93' b = b+1; */
      b++;
    }
  }

  emxFree_real_T(&r1);
  emxFree_real_T(&r);
  emxFree_cell_wrap_4(&p);
  emxInit_cell_wrap_4(&pwm);
  emxInit_cell_wrap_4(&lpwm);

  /* 'mapTF:97' pwm = cell(b,1); */
  unnamed_idx_0 = (int)b;
  i = pwm->size[0];
  pwm->size[0] = (int)b;
  emxEnsureCapacity_cell_wrap_4(pwm, i);
  pwm_data = pwm->data;

  /* 'mapTF:98' pwm = coder.nullcopy(pwm); */
  /* 'mapTF:99' lpwm = cell(b,1); */
  i = lpwm->size[0];
  lpwm->size[0] = (int)b;
  emxEnsureCapacity_cell_wrap_4(lpwm, i);
  PWM2_data = lpwm->data;
  for (i = 0; i < unnamed_idx_0; i++) {
    pwm_data[i].f1->size[0] = 0;
    pwm_data[i].f1->size[1] = 4;
    PWM2_data[i].f1->size[0] = 0;
    PWM2_data[i].f1->size[1] = 4;
  }

  /* 'mapTF:100' lpwm = coder.nullcopy(lpwm); */
  /* 'mapTF:101' lab = zeros(b,1); */
  /* 'mapTF:102' b = 1; */
  b = 1.0;

  /* 'mapTF:103' for i = 1:a */
  i = len->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:104' for j = 1:LEN(i) */
    i1 = (int)LEN_data[b_i];
    for (j = 0; j < i1; j++) {
      /* 'mapTF:105' pwm{b} = PWM{i}(j:j+l_svm-1,:); */
      GC = (((double)j + 1.0) + varargin_7) - 1.0;
      if ((double)j + 1.0 > GC) {
        i2 = 0;
        sizes_idx_0 = 0;
      } else {
        i2 = j;
        sizes_idx_0 = (int)GC;
      }

      loop_ub = sizes_idx_0 - i2;
      sizes_idx_0 = (int)(b + (double)j) - 1;
      unnamed_idx_0 = pwm_data[sizes_idx_0].f1->size[0] * pwm_data[sizes_idx_0].
        f1->size[1];
      pwm_data[(int)(b + (double)j) - 1].f1->size[0] = loop_ub;
      pwm_data[(int)(b + (double)j) - 1].f1->size[1] = 4;
      emxEnsureCapacity_real_T(pwm_data[(int)(b + (double)j) - 1].f1,
        unnamed_idx_0);
      for (unnamed_idx_0 = 0; unnamed_idx_0 < 4; unnamed_idx_0++) {
        for (vectorUB = 0; vectorUB < loop_ub; vectorUB++) {
          pwm_data[(int)(b + (double)j) - 1].f1->data[vectorUB +
            pwm_data[sizes_idx_0].f1->size[0] * unnamed_idx_0] = PWM_data[b_i].
            f1->data[(i2 + vectorUB) + PWM_data[b_i].f1->size[0] * unnamed_idx_0];
        }
      }

      /* 'mapTF:106' lpwm{b} = log((pwm{b}+10^-10)/(1+4*10^-10)); */
      i2 = PWM2_data[sizes_idx_0].f1->size[0] * PWM2_data[sizes_idx_0].f1->size
        [1];
      PWM2_data[(int)(b + (double)j) - 1].f1->size[0] = pwm_data[sizes_idx_0].
        f1->size[0];
      PWM2_data[(int)(b + (double)j) - 1].f1->size[1] = 4;
      emxEnsureCapacity_real_T(PWM2_data[(int)(b + (double)j) - 1].f1, i2);
      loop_ub = pwm_data[sizes_idx_0].f1->size[0] * 4;
      unnamed_idx_0 = (loop_ub / 2) << 1;
      vectorUB = unnamed_idx_0 - 2;
      for (i2 = 0; i2 <= vectorUB; i2 += 2) {
        r2 = _mm_loadu_pd(&pwm_data[sizes_idx_0].f1->data[i2]);
        _mm_storeu_pd(&PWM2_data[sizes_idx_0].f1->data[i2], _mm_div_pd
                      (_mm_add_pd(r2, _mm_set1_pd(1.0E-10)), _mm_set1_pd
                       (1.0000000004)));
      }

      for (i2 = unnamed_idx_0; i2 < loop_ub; i2++) {
        PWM2_data[(int)(b + (double)j) - 1].f1->data[i2] = (pwm_data[sizes_idx_0]
          .f1->data[i2] + 1.0E-10) / 1.0000000004;
      }

      nx = PWM2_data[sizes_idx_0].f1->size[0] << 2;
      for (unnamed_idx_0 = 0; unnamed_idx_0 < nx; unnamed_idx_0++) {
        PWM2_data[(int)(b + (double)j) - 1].f1->data[unnamed_idx_0] = log
          (PWM2_data[sizes_idx_0].f1->data[unnamed_idx_0]);
      }

      /* 'mapTF:107' lab(b) = i; */
      /* 'mapTF:108' b = b+1; */
    }

    b += (double)(i1 - 1) + 1.0;
  }

  emxFree_cell_wrap_4(&pwm);
  emxFree_cell_wrap_4(&PWM);
  emxInit_real_T(&c, 2);
  emxInit_real_T(&a__3, 1);
  emxInit_real_T(&seqmat, 2);
  emxInit_real_T(&f, 1);
  emxInit_real_T(&pwm_prob, 2);

  /* 'mapTF:112' [c,~,~,~,~,rcnum] = genIndex(l_svm,k_svm,nfrac); */
  genIndex(varargin_7, varargin_8, varargin_9, seqmat, pwm_prob, f, a__3, c, &GC);
  seqmat_data = seqmat->data;

  /* 'mapTF:113' c2 = c(1:numel(c)/k_svm-rcnum,:); */
  GC = (double)(seqmat->size[0] * seqmat->size[1]) / varargin_8 - GC;
  if (1.0 > GC) {
    loop_ub = 0;
  } else {
    loop_ub = (int)GC;
  }

  /* 'mapTF:114' c = [c;l_svm+1-fliplr(c2)]; */
  b_loop_ub = seqmat->size[1];
  i = pwm_prob->size[0] * pwm_prob->size[1];
  pwm_prob->size[0] = loop_ub;
  pwm_prob->size[1] = seqmat->size[1];
  emxEnsureCapacity_real_T(pwm_prob, i);
  pwm_prob_data = pwm_prob->data;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      pwm_prob_data[i1 + pwm_prob->size[0] * i] = seqmat_data[i1 + seqmat->size
        [0] * i];
    }
  }

  b_fliplr(pwm_prob);
  pwm_prob_data = pwm_prob->data;
  loop_ub = pwm_prob->size[0] * pwm_prob->size[1];
  unnamed_idx_0 = (loop_ub / 2) << 1;
  vectorUB = unnamed_idx_0 - 2;
  for (i = 0; i <= vectorUB; i += 2) {
    r2 = _mm_loadu_pd(&pwm_prob_data[i]);
    _mm_storeu_pd(&pwm_prob_data[i], _mm_sub_pd(_mm_set1_pd(varargin_7 + 1.0),
      r2));
  }

  for (i = unnamed_idx_0; i < loop_ub; i++) {
    pwm_prob_data[i] = (varargin_7 + 1.0) - pwm_prob_data[i];
  }

  if ((seqmat->size[0] != 0) && (seqmat->size[1] != 0)) {
    unnamed_idx_0 = seqmat->size[1];
  } else if ((pwm_prob->size[0] != 0) && (pwm_prob->size[1] != 0)) {
    unnamed_idx_0 = pwm_prob->size[1];
  } else {
    unnamed_idx_0 = seqmat->size[1];
    if (pwm_prob->size[1] > seqmat->size[1]) {
      unnamed_idx_0 = pwm_prob->size[1];
    }
  }

  empty_non_axis_sizes = (unnamed_idx_0 == 0);
  if (empty_non_axis_sizes || ((seqmat->size[0] != 0) && (seqmat->size[1] != 0)))
  {
    nx = seqmat->size[0];
  } else {
    nx = 0;
  }

  if (empty_non_axis_sizes || ((pwm_prob->size[0] != 0) && (pwm_prob->size[1] !=
        0))) {
    sizes_idx_0 = pwm_prob->size[0];
  } else {
    sizes_idx_0 = 0;
  }

  varargin_7_idx_0 = sizes_idx_0;
  i = c->size[0] * c->size[1];
  c->size[0] = nx + sizes_idx_0;
  c->size[1] = unnamed_idx_0;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  for (i = 0; i < unnamed_idx_0; i++) {
    for (i1 = 0; i1 < nx; i1++) {
      c_data[i1 + c->size[0] * i] = seqmat_data[i1 + nx * i];
    }
  }

  for (i = 0; i < unnamed_idx_0; i++) {
    for (i1 = 0; i1 < sizes_idx_0; i1++) {
      c_data[(i1 + nx) + c->size[0] * i] = pwm_prob_data[i1 + sizes_idx_0 * i];
    }
  }

  /* 'mapTF:115' C = numel(c)/k_svm; */
  GC = (double)(c->size[0] * c->size[1]) / varargin_8;

  /* 'mapTF:116' seqmat = zeros(C,l_svm); */
  i = (int)GC;
  i1 = seqmat->size[0] * seqmat->size[1];
  seqmat->size[0] = (int)GC;
  nx = (int)varargin_7;
  seqmat->size[1] = (int)varargin_7;
  emxEnsureCapacity_real_T(seqmat, i1);
  seqmat_data = seqmat->data;
  loop_ub = (int)GC * (int)varargin_7;
  for (i1 = 0; i1 < loop_ub; i1++) {
    seqmat_data[i1] = 0.0;
  }

  /* 'mapTF:117' for i = 1:C */
  emxInit_int32_T(&r3, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF:118' seqmat(i,c(i,:)) = 1; */
    loop_ub = c->size[1];
    i1 = r3->size[0];
    r3->size[0] = c->size[1];
    emxEnsureCapacity_int32_T(r3, i1);
    r4 = r3->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      r4[i1] = (int)c_data[b_i + c->size[0] * i1];
    }

    loop_ub = r3->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      seqmat_data[b_i + seqmat->size[0] * (r4[i1] - 1)] = 1.0;
    }
  }

  emxInit_cell_wrap_5(&Smat, 1);

  /* 'mapTF:120' Smat = cell(l_svm,1); */
  i = Smat->size[0];
  Smat->size[0] = (int)varargin_7;
  emxEnsureCapacity_cell_wrap_5(Smat, i);
  Smat_data = Smat->data;
  for (i = 0; i < nx; i++) {
    Smat_data[i].f1->size[0] = 0;
    Smat_data[i].f1->size[1] = 0;
  }

  /* 'mapTF:121' Smat = coder.nullcopy(Smat); */
  /* 'mapTF:122' for i = 1:l_svm */
  emxInit_boolean_T(&b_f, 1);
  for (b_i = 0; b_i < nx; b_i++) {
    /* 'mapTF:123' f = find(prod(c-i,2)==0); */
    i = pwm_prob->size[0] * pwm_prob->size[1];
    pwm_prob->size[0] = c->size[0];
    pwm_prob->size[1] = c->size[1];
    emxEnsureCapacity_real_T(pwm_prob, i);
    pwm_prob_data = pwm_prob->data;
    loop_ub = c->size[0] * c->size[1];
    unnamed_idx_0 = (loop_ub / 2) << 1;
    vectorUB = unnamed_idx_0 - 2;
    for (i = 0; i <= vectorUB; i += 2) {
      r2 = _mm_loadu_pd(&c_data[i]);
      _mm_storeu_pd(&pwm_prob_data[i], _mm_sub_pd(r2, _mm_set1_pd((double)b_i +
        1.0)));
    }

    for (i = unnamed_idx_0; i < loop_ub; i++) {
      pwm_prob_data[i] = c_data[i] - ((double)b_i + 1.0);
    }

    prod(pwm_prob, f);
    f_data = f->data;
    i = b_f->size[0];
    b_f->size[0] = f->size[0];
    emxEnsureCapacity_boolean_T(b_f, i);
    b_f_data = b_f->data;
    loop_ub = f->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_f_data[i] = (f_data[i] == 0.0);
    }

    eml_find(b_f, r3);
    r4 = r3->data;
    i = f->size[0];
    f->size[0] = r3->size[0];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = r3->size[0];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = r4[i];
    }

    /* 'mapTF:124' Smat{i} = zeros(length(f), l_svm); */
    i = Smat_data[b_i].f1->size[0] * Smat_data[b_i].f1->size[1];
    Smat_data[b_i].f1->size[0] = f->size[0];
    Smat_data[b_i].f1->size[1] = (int)varargin_7;
    emxEnsureCapacity_real_T(Smat_data[b_i].f1, i);
    loop_ub = f->size[0] * (int)varargin_7;
    for (i = 0; i < loop_ub; i++) {
      Smat_data[b_i].f1->data[i] = 0.0;
    }

    /* 'mapTF:125' for j = 1:length(f) */
    i = f->size[0];
    for (j = 0; j < i; j++) {
      /* 'mapTF:126' Smat{i}(j,c(f(j),:)) = 1; */
      loop_ub = c->size[1];
      i1 = r3->size[0];
      r3->size[0] = c->size[1];
      emxEnsureCapacity_int32_T(r3, i1);
      r4 = r3->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        r4[i1] = (int)c_data[((int)f_data[j] + c->size[0] * i1) - 1];
      }

      loop_ub = r3->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        Smat_data[b_i].f1->data[j + Smat_data[b_i].f1->size[0] * (r4[i1] - 1)] =
          1.0;
      }
    }
  }

  emxFree_boolean_T(&b_f);
  emxFree_int32_T(&r3);
  emxInit_int8_T(&IND, 1);

  /* 'mapTF:129' L = length(ss); */
  /* 'mapTF:130' B = b-1; */
  /* 'mapTF:131' maxnorm = zeros(B,1); */
  /* 'mapTF:132' minnorm = zeros(B,1); */
  /* 'mapTF:133' vec = zeros(l_svm,1); */
  /* 'mapTF:134' IND = zeros(4^l_svm,1); */
  GC = rt_powd_snf(4.0, varargin_7);
  unnamed_idx_0 = (int)rt_powd_snf(4.0, varargin_7);
  i = IND->size[0];
  IND->size[0] = (int)GC;
  emxEnsureCapacity_int8_T(IND, i);
  IND_data = IND->data;
  for (i = 0; i < unnamed_idx_0; i++) {
    IND_data[i] = 0;
  }

  /* 'mapTF:135' kmat = zeros(B,4^l_svm); */
  i = c->size[0] * c->size[1];
  c->size[0] = (int)((unsigned int)b - 1U);
  c->size[1] = (int)GC;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  loop_ub = (int)((unsigned int)b - 1U) * (int)GC;
  for (i = 0; i < loop_ub; i++) {
    c_data[i] = 0.0;
  }

  emxInit_real_T(&maxnorm, 1);
  emxInit_real_T(&minnorm, 1);

  /* 'mapTF:136' for j = 1:B */
  i = (int)(unsigned int)b;
  i1 = maxnorm->size[0];
  maxnorm->size[0] = (int)((unsigned int)b - 1U);
  emxEnsureCapacity_real_T(maxnorm, i1);
  maxnorm_data = maxnorm->data;
  i1 = minnorm->size[0];
  minnorm->size[0] = (int)((unsigned int)b - 1U);
  emxEnsureCapacity_real_T(minnorm, i1);
  minnorm_data = minnorm->data;
  emxInit_real_T(&x, 2);
  emxInit_real_T(&B, 2);
  emxInit_real_T(&b_lpwm, 2);
  for (j = 0; j <= i - 2; j++) {
    /* 'mapTF:137' vec = max(lpwm{j}'); */
    /* 'mapTF:138' vec2 = min(lpwm{j}'); */
    /* 'mapTF:139' maxnorm(j) = sum(exp(seqmat*vec')); */
    i1 = b_lpwm->size[0] * b_lpwm->size[1];
    b_lpwm->size[0] = 4;
    loop_ub = PWM2_data[j].f1->size[0];
    b_lpwm->size[1] = PWM2_data[j].f1->size[0];
    emxEnsureCapacity_real_T(b_lpwm, i1);
    f_data = b_lpwm->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[4 * i1] = PWM2_data[j].f1->data[i1];
      f_data[4 * i1 + 1] = PWM2_data[j].f1->data[i1 + PWM2_data[j].f1->size[0]];
      f_data[4 * i1 + 2] = PWM2_data[j].f1->data[i1 + PWM2_data[j].f1->size[0] *
        2];
      f_data[4 * i1 + 3] = PWM2_data[j].f1->data[i1 + PWM2_data[j].f1->size[0] *
        3];
    }

    c_maximum(b_lpwm, B);
    mtimes(seqmat, B, pwm_prob);
    pwm_prob_data = pwm_prob->data;
    i1 = x->size[0] * x->size[1];
    x->size[0] = pwm_prob->size[0];
    x->size[1] = pwm_prob->size[1];
    emxEnsureCapacity_real_T(x, i1);
    f_data = x->data;
    loop_ub = pwm_prob->size[0] * pwm_prob->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[i1] = pwm_prob_data[i1];
    }

    nx = pwm_prob->size[0];
    for (unnamed_idx_0 = 0; unnamed_idx_0 < nx; unnamed_idx_0++) {
      f_data[unnamed_idx_0] = exp(f_data[unnamed_idx_0]);
    }

    b_sum(x, B);
    shift_data = B->data;
    maxnorm_data[j] = shift_data[0];

    /* 'mapTF:140' minnorm(j) = sum(exp(seqmat*vec2')); */
    i1 = b_lpwm->size[0] * b_lpwm->size[1];
    b_lpwm->size[0] = 4;
    loop_ub = PWM2_data[j].f1->size[0];
    b_lpwm->size[1] = PWM2_data[j].f1->size[0];
    emxEnsureCapacity_real_T(b_lpwm, i1);
    f_data = b_lpwm->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[4 * i1] = PWM2_data[j].f1->data[i1];
      f_data[4 * i1 + 1] = PWM2_data[j].f1->data[i1 + PWM2_data[j].f1->size[0]];
      f_data[4 * i1 + 2] = PWM2_data[j].f1->data[i1 + PWM2_data[j].f1->size[0] *
        2];
      f_data[4 * i1 + 3] = PWM2_data[j].f1->data[i1 + PWM2_data[j].f1->size[0] *
        3];
    }

    b_minimum(b_lpwm, B);
    shift_data = B->data;
    if ((seqmat->size[0] == 0) || (seqmat->size[1] == 0) || (B->size[1] == 0)) {
      i1 = a__3->size[0];
      a__3->size[0] = seqmat->size[0];
      emxEnsureCapacity_real_T(a__3, i1);
      LEN_2_data = a__3->data;
      loop_ub = seqmat->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        LEN_2_data[i1] = 0.0;
      }
    } else {
      i1 = a__3->size[0];
      a__3->size[0] = seqmat->size[0];
      emxEnsureCapacity_real_T(a__3, i1);
      LEN_2_data = a__3->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)seqmat->
                  size[0], (blasint)1, (blasint)seqmat->size[1], 1.0,
                  &seqmat_data[0], (blasint)seqmat->size[0], &shift_data[0],
                  (blasint)1, 0.0, &LEN_2_data[0], (blasint)seqmat->size[0]);
    }

    nx = a__3->size[0];
    for (unnamed_idx_0 = 0; unnamed_idx_0 < nx; unnamed_idx_0++) {
      LEN_2_data[unnamed_idx_0] = exp(LEN_2_data[unnamed_idx_0]);
    }

    minnorm_data[j] = blockedSummation(a__3, a__3->size[0]);
  }

  emxFree_real_T(&b_lpwm);
  emxFree_real_T(&B);
  emxFree_real_T(&x);
  emxInit_real_T(&all_pwm, 3);

  /* 'mapTF:142' dnorm = maxnorm; */
  /* 'mapTF:143' vec = zeros(l_svm,1); */
  /* 'mapTF:144' all_pwm = zeros(4, l_svm, B); */
  i1 = all_pwm->size[0] * all_pwm->size[1] * all_pwm->size[2];
  all_pwm->size[0] = 4;
  all_pwm->size[1] = (int)varargin_7;
  all_pwm->size[2] = (int)((unsigned int)b - 1U);
  emxEnsureCapacity_real_T(all_pwm, i1);
  all_pwm_data = all_pwm->data;
  loop_ub = ((int)varargin_7 << 2) * (int)((unsigned int)b - 1U);
  for (i1 = 0; i1 < loop_ub; i1++) {
    all_pwm_data[i1] = 0.0;
  }

  /* 'mapTF:145' for cur_idx=1:B */
  for (unnamed_idx_0 = 0; unnamed_idx_0 <= i - 2; unnamed_idx_0++) {
    /* 'mapTF:146' all_pwm(:,:,cur_idx) = lpwm{cur_idx}'; */
    loop_ub = PWM2_data[unnamed_idx_0].f1->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      all_pwm_data[4 * i1 + 4 * all_pwm->size[1] * unnamed_idx_0] =
        PWM2_data[unnamed_idx_0].f1->data[i1];
      all_pwm_data[(4 * i1 + 4 * all_pwm->size[1] * unnamed_idx_0) + 1] =
        PWM2_data[unnamed_idx_0].f1->data[i1 + PWM2_data[unnamed_idx_0].f1->
        size[0]];
      all_pwm_data[(4 * i1 + 4 * all_pwm->size[1] * unnamed_idx_0) + 2] =
        PWM2_data[unnamed_idx_0].f1->data[i1 + PWM2_data[unnamed_idx_0].f1->
        size[0] * 2];
      all_pwm_data[(4 * i1 + 4 * all_pwm->size[1] * unnamed_idx_0) + 3] =
        PWM2_data[unnamed_idx_0].f1->data[i1 + PWM2_data[unnamed_idx_0].f1->
        size[0] * 3];
    }
  }

  emxFree_cell_wrap_4(&lpwm);
  emxInit_cell_wrap_4(&LL);

  /* 'mapTF:149' LL = cell(length(V),1); */
  unnamed_idx_0 = V->size[0];
  i1 = LL->size[0];
  LL->size[0] = V->size[0];
  emxEnsureCapacity_cell_wrap_4(LL, i1);
  pwm_data = LL->data;
  for (i1 = 0; i1 < unnamed_idx_0; i1++) {
    pwm_data[i1].f1->size[0] = 0;
    pwm_data[i1].f1->size[1] = 4;
  }

  emxInit_cell_wrap_0(&VV);

  /* 'mapTF:150' LL = coder.nullcopy(LL); */
  /* 'mapTF:151' VV = cell(length(V),1); */
  unnamed_idx_0 = V->size[0];
  i1 = VV->size[0];
  VV->size[0] = V->size[0];
  emxEnsureCapacity_cell_wrap_0(VV, i1);
  VV_data = VV->data;
  for (i1 = 0; i1 < unnamed_idx_0; i1++) {
    VV_data[i1].f1->size[0] = 0;
  }

  emxInit_cell_wrap_6(&NN);

  /* 'mapTF:152' VV = coder.nullcopy(VV); */
  /* 'mapTF:153' NN = cell(length(V),1); */
  unnamed_idx_0 = V->size[0];
  i1 = NN->size[0];
  NN->size[0] = V->size[0];
  emxEnsureCapacity_cell_wrap_6(NN, i1);
  NN_data = NN->data;
  for (i1 = 0; i1 < unnamed_idx_0; i1++) {
    NN_data[i1].f1->size[0] = 0;
  }

  /* 'mapTF:154' NN = coder.nullcopy(NN); */
  /* 'mapTF:155' fprintf('Mapping motifs\n'); */
  printf("Mapping motifs\n");
  fflush(stdout);

  /* 'mapTF:156' tic */
  tic();

  /* 'mapTF:157' for I = 1:length(ss) */
  i1 = ss->size[0];
  emxInit_int8_T(&ss_onehot, 2);
  emxInit_real_T(&b_B, 1);
  emxInit_real_T(&b_all_pwm, 3);
  for (b_I = 0; b_I < i1; b_I++) {
    /* 'mapTF:158' seq2 = ss{I}; */
    /* 'mapTF:160' ss_onehot = zeros(4, length(seq2)); */
    i2 = ss_onehot->size[0] * ss_onehot->size[1];
    ss_onehot->size[0] = 4;
    ss_onehot->size[1] = ss_data[b_I].f1->size[1];
    emxEnsureCapacity_int8_T(ss_onehot, i2);
    ss_onehot_data = ss_onehot->data;
    loop_ub = ss_data[b_I].f1->size[1] << 2;
    for (i2 = 0; i2 < loop_ub; i2++) {
      ss_onehot_data[i2] = 0;
    }

    /* 'mapTF:161' for idx=1:length(seq2) */
    i2 = ss_data[b_I].f1->size[1];
    for (unnamed_idx_0 = 0; unnamed_idx_0 < i2; unnamed_idx_0++) {
      /* 'mapTF:162' ss_onehot(seq2(idx), idx) = 1; */
      ss_onehot_data[((int)ss_data[b_I].f1->data[unnamed_idx_0] + 4 *
                      unnamed_idx_0) - 1] = 1;
    }

    /* 'mapTF:164' pwm_prob = zeros(b-1,length(seqindmat{I})); */
    i2 = pwm_prob->size[0] * pwm_prob->size[1];
    pwm_prob->size[0] = (int)((unsigned int)b - 1U);
    pwm_prob->size[1] = seqindmat_data[b_I].f1->size[0];
    emxEnsureCapacity_real_T(pwm_prob, i2);
    pwm_prob_data = pwm_prob->data;
    loop_ub = (int)((unsigned int)b - 1U) * seqindmat_data[b_I].f1->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      pwm_prob_data[i2] = 0.0;
    }

    /* 'mapTF:166' for i = 1:length(seqindmat{I}) */
    i2 = seqindmat_data[b_I].f1->size[0];
    for (b_i = 0; b_i < i2; b_i++) {
      /* 'mapTF:167' ind = seqindmat{I}(i); */
      /* 'mapTF:168' if IND(ind) == 0 */
      c_loop_ub = (int)seqindmat_data[b_I].f1->data[b_i];
      if (IND_data[c_loop_ub - 1] == 0) {
        /* 'mapTF:169' IND(ind) = 1; */
        IND_data[c_loop_ub - 1] = 1;

        /*  SEQ = seq2(i:i+l_svm-1); */
        /*  for j = 1:B */
        /*      for jj = 1:l_svm */
        /*          vec(jj) = lpwm{j}(jj,SEQ(jj)); */
        /*      end */
        /*      kmat(j,ind) = sum(exp(seqmat*vec)); */
        /*  end */
        /* 'mapTF:177' seq3 = ss_onehot(:,i:i+l_svm-1); */
        GC = (((double)b_i + 1.0) + varargin_7) - 1.0;
        if ((double)b_i + 1.0 > GC) {
          sizes_idx_0 = 0;
          unnamed_idx_0 = -1;
        } else {
          sizes_idx_0 = b_i;
          unnamed_idx_0 = (int)GC - 1;
        }

        /* 'mapTF:178' vec3 = reshape(nonzeros(all_pwm .* seq3), l_svm, b-1); */
        if (all_pwm->size[1] == (unnamed_idx_0 - sizes_idx_0) + 1) {
          unnamed_idx_0 = b_all_pwm->size[0] * b_all_pwm->size[1] *
            b_all_pwm->size[2];
          b_all_pwm->size[0] = 4;
          b_all_pwm->size[1] = all_pwm->size[1];
          loop_ub = all_pwm->size[2];
          b_all_pwm->size[2] = all_pwm->size[2];
          emxEnsureCapacity_real_T(b_all_pwm, unnamed_idx_0);
          f_data = b_all_pwm->data;
          for (unnamed_idx_0 = 0; unnamed_idx_0 < loop_ub; unnamed_idx_0++) {
            b_loop_ub = all_pwm->size[1];
            for (vectorUB = 0; vectorUB < b_loop_ub; vectorUB++) {
              nx = sizes_idx_0 + vectorUB;
              f_data[4 * vectorUB + 4 * b_all_pwm->size[1] * unnamed_idx_0] =
                all_pwm_data[4 * vectorUB + 4 * all_pwm->size[1] * unnamed_idx_0]
                * (double)ss_onehot_data[4 * nx];
              f_data[(4 * vectorUB + 4 * b_all_pwm->size[1] * unnamed_idx_0) + 1]
                = all_pwm_data[(4 * vectorUB + 4 * all_pwm->size[1] *
                                unnamed_idx_0) + 1] * (double)ss_onehot_data[4 *
                nx + 1];
              f_data[(4 * vectorUB + 4 * b_all_pwm->size[1] * unnamed_idx_0) + 2]
                = all_pwm_data[(4 * vectorUB + 4 * all_pwm->size[1] *
                                unnamed_idx_0) + 2] * (double)ss_onehot_data[4 *
                nx + 2];
              f_data[(4 * vectorUB + 4 * b_all_pwm->size[1] * unnamed_idx_0) + 3]
                = all_pwm_data[(4 * vectorUB + 4 * all_pwm->size[1] *
                                unnamed_idx_0) + 3] * (double)ss_onehot_data[4 *
                nx + 3];
            }
          }

          nonzeros(b_all_pwm, f);
          f_data = f->data;
        } else {
          b_binary_expand_op(f, all_pwm, ss_onehot, sizes_idx_0, unnamed_idx_0,
                             sizes_idx_0 - 1);
          f_data = f->data;
        }

        /* 'mapTF:179' for j = 1:B */
        if (0 <= (int)(unsigned int)b - 2) {
          varargin_7_idx_0 = (int)varargin_7;
          i3 = (int)varargin_7;
          d_loop_ub = (int)varargin_7;
        }

        for (j = 0; j <= i - 2; j++) {
          /* 'mapTF:180' kmat(j,ind) = sum(exp(seqmat*vec3(:,j))); */
          sizes_idx_0 = b_B->size[0];
          b_B->size[0] = i3;
          emxEnsureCapacity_real_T(b_B, sizes_idx_0);
          shift_data = b_B->data;
          for (sizes_idx_0 = 0; sizes_idx_0 < d_loop_ub; sizes_idx_0++) {
            shift_data[sizes_idx_0] = f_data[sizes_idx_0 + varargin_7_idx_0 * j];
          }

          loop_ub = seqmat->size[0];
          if ((seqmat->size[0] == 0) || (seqmat->size[1] == 0) || ((int)
               varargin_7 == 0)) {
            sizes_idx_0 = a__3->size[0];
            a__3->size[0] = seqmat->size[0];
            emxEnsureCapacity_real_T(a__3, sizes_idx_0);
            LEN_2_data = a__3->data;
            for (sizes_idx_0 = 0; sizes_idx_0 < loop_ub; sizes_idx_0++) {
              LEN_2_data[sizes_idx_0] = 0.0;
            }
          } else {
            sizes_idx_0 = a__3->size[0];
            a__3->size[0] = seqmat->size[0];
            emxEnsureCapacity_real_T(a__3, sizes_idx_0);
            LEN_2_data = a__3->data;
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                        seqmat->size[0], (blasint)1, (blasint)seqmat->size[1],
                        1.0, &seqmat_data[0], (blasint)seqmat->size[0],
                        &shift_data[0], (blasint)(int)varargin_7, 0.0,
                        &LEN_2_data[0], (blasint)seqmat->size[0]);
          }

          nx = a__3->size[0];
          for (unnamed_idx_0 = 0; unnamed_idx_0 < nx; unnamed_idx_0++) {
            LEN_2_data[unnamed_idx_0] = exp(LEN_2_data[unnamed_idx_0]);
          }

          c_data[j + c->size[0] * (c_loop_ub - 1)] = blockedSummation(a__3,
            a__3->size[0]);
        }

        /* 'mapTF:182' kmat(:,ind) = log((kmat(:,ind)-minnorm)./dnorm); */
        loop_ub = c->size[0];
        if (c->size[0] == 1) {
          b_loop_ub = minnorm->size[0];
        } else {
          b_loop_ub = c->size[0];
        }

        if ((c->size[0] == minnorm->size[0]) && (b_loop_ub == maxnorm->size[0]))
        {
          sizes_idx_0 = f->size[0];
          f->size[0] = c->size[0];
          emxEnsureCapacity_real_T(f, sizes_idx_0);
          f_data = f->data;
          unnamed_idx_0 = (c->size[0] / 2) << 1;
          vectorUB = unnamed_idx_0 - 2;
          for (sizes_idx_0 = 0; sizes_idx_0 <= vectorUB; sizes_idx_0 += 2) {
            r2 = _mm_loadu_pd(&c_data[sizes_idx_0 + c->size[0] * (c_loop_ub - 1)]);
            r5 = _mm_loadu_pd(&minnorm_data[sizes_idx_0]);
            r6 = _mm_loadu_pd(&maxnorm_data[sizes_idx_0]);
            _mm_storeu_pd(&f_data[sizes_idx_0], _mm_div_pd(_mm_sub_pd(r2, r5),
              r6));
          }

          for (sizes_idx_0 = unnamed_idx_0; sizes_idx_0 < loop_ub; sizes_idx_0++)
          {
            f_data[sizes_idx_0] = (c_data[sizes_idx_0 + c->size[0] * (c_loop_ub
              - 1)] - minnorm_data[sizes_idx_0]) / maxnorm_data[sizes_idx_0];
          }
        } else {
          binary_expand_op(f, c, seqindmat, b_I, b_i, minnorm, maxnorm);
          f_data = f->data;
        }

        nx = f->size[0];
        for (unnamed_idx_0 = 0; unnamed_idx_0 < nx; unnamed_idx_0++) {
          f_data[unnamed_idx_0] = log(f_data[unnamed_idx_0]);
        }

        loop_ub = f->size[0];
        for (sizes_idx_0 = 0; sizes_idx_0 < loop_ub; sizes_idx_0++) {
          c_data[sizes_idx_0 + c->size[0] * (c_loop_ub - 1)] =
            f_data[sizes_idx_0];
        }

        /* 'mapTF:183' pwm_prob(:,i) = kmat(:,ind); */
        loop_ub = c->size[0];
        for (sizes_idx_0 = 0; sizes_idx_0 < loop_ub; sizes_idx_0++) {
          pwm_prob_data[sizes_idx_0 + pwm_prob->size[0] * b_i] =
            c_data[sizes_idx_0 + c->size[0] * (c_loop_ub - 1)];
        }
      } else {
        /* 'mapTF:184' else */
        /* 'mapTF:185' pwm_prob(:,i) = kmat(:,ind); */
        loop_ub = c->size[0];
        for (sizes_idx_0 = 0; sizes_idx_0 < loop_ub; sizes_idx_0++) {
          pwm_prob_data[sizes_idx_0 + pwm_prob->size[0] * b_i] =
            c_data[sizes_idx_0 + c->size[0] * (c_loop_ub - 1)];
        }
      }
    }

    /* 'mapTF:188' [LL{I}, NN{I}] = MAPTF(fn, ss{I}, GC, pwm_prob, l_svm, k_svm, ofn, PWM, PWM2, pwm, lpwm, lab, LEN, LEN_2, shift, P{I}, names, len, a,b, I); */
    i2 = f->size[0];
    f->size[0] = LEN->size[0];
    emxEnsureCapacity_real_T(f, i2);
    f_data = f->data;
    loop_ub = LEN->size[0] - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
      f_data[i2] = LEN_data[i2];
    }

    MAPTF(ss_data[b_I].f1, pwm_prob, varargin_7, PWM2, f, LEN_2, shift,
          P_data[b_I].f1, names, len->size[0], pwm_data[b_I].f1, NN_data[b_I].f1);

    /* 'mapTF:189' if numel(LL{I}) > 0 */
    if ((pwm_data[b_I].f1->size[0] << 2) > 0) {
      /* 'mapTF:190' VV{I} = scoreseqkmer(PWM3, lPWM3, LL{I}, ss{I}, Smat, l_svm, k_svm, ofn, V{I}); */
      scoreseqkmer(PWM3, lPWM3, pwm_data[b_I].f1, ss_data[b_I].f1, Smat,
                   varargin_7, V_data[b_I].f1, VV_data[b_I].f1);
    }

    /* 'mapTF:192' if mod(I,100)==0 */
    if (b_mod((double)b_I + 1.0, 100.0) == 0.0) {
      /* 'mapTF:193' fprintf('%d out of %d sequences done...\n', int32(I), int32(length(ss))); */
      printf("%d out of %d sequences done...\n", b_I + 1, ss->size[0]);
      fflush(stdout);

      /* 'mapTF:194' toc */
      toc();
    }
  }

  emxFree_real_T(&b_all_pwm);
  emxFree_real_T(&b_B);
  emxFree_real_T(&len);
  emxFree_cell_wrap_3(&names);
  emxFree_cell_wrap_0(&seqindmat);
  emxFree_cell_wrap_1(&V);
  emxFree_cell_wrap_0(&P);
  emxFree_real_T(&pwm_prob);
  emxFree_int8_T(&ss_onehot);
  emxFree_real_T(&all_pwm);
  emxFree_int8_T(&IND);
  emxFree_real_T(&minnorm);
  emxFree_real_T(&maxnorm);
  emxFree_real_T(&f);
  emxFree_cell_wrap_5(&Smat);
  emxFree_real_T(&seqmat);
  emxFree_real_T(&a__3);
  emxFree_real_T(&c);
  emxFree_real_T(&shift);
  emxFree_real_T(&LEN_2);
  emxFree_real_T(&LEN);
  emxFree_cell_wrap_4(&lPWM3);
  emxFree_cell_wrap_4(&PWM3);
  emxFree_cell_wrap_4(&PWM2);

  /* 'mapTF:197' fprintf('%d out of %d sequences done...\n', int32(length(ss)), int32(length(ss))); */
  printf("%d out of %d sequences done...\n", ss->size[0], ss->size[0]);
  fflush(stdout);

  /*  clear kmat */
  /* 'mapTF:199' PWM_corr(ofn, VV, NN, LL, seq); */
  PWM_corr(varargin_6, VV, NN, LL, seq);
  emxFree_cell_wrap_3(&seq);
  emxFree_cell_wrap_2(&ss);
  emxFree_cell_wrap_6(&NN);
  emxFree_cell_wrap_0(&VV);
  emxFree_cell_wrap_4(&LL);
}

/* End of code generation (mapTF.c) */
