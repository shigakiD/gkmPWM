/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mapTF2_ls.c
 *
 * Code generation for function 'mapTF2_ls'
 *
 */

/* Include files */
#include "mapTF2_ls.h"
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
#include "mapTF2_ls_data.h"
#include "mapTF2_ls_emxutil.h"
#include "mapTF2_ls_initialize.h"
#include "mapTF2_ls_rtwutil.h"
#include "mapTF2_ls_types.h"
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
#include "unique.h"
#include "cblas.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_cell_wrap_15
#define typedef_cell_wrap_15

typedef struct {
  double f1[3];
} cell_wrap_15;

#endif                                 /* typedef_cell_wrap_15 */

/* Variable Definitions */
static const char cv[11] = { '_', 'm', 'o', 't', 'i', 'f', 's', '.', 'o', 'u',
  't' };

/* Function Declarations */
static void MAPTF(const emxArray_char_T *s, const emxArray_real_T *ss, const
                  emxArray_real_T *pwm_prob, const emxArray_cell_wrap_4 *Smat,
                  double l_svm, const emxArray_char_T *ofn, const
                  emxArray_cell_wrap_3 *PWM2, emxArray_real_T *LEN, const
                  emxArray_real_T *LEN_2, const emxArray_real_T *shift,
                  emxArray_real_T *gkmprob, const emxArray_real_T *dsvm, const
                  emxArray_cell_wrap_2 *names, double a, double rnum);
static void PWM_corr(emxArray_real_T *kmer, emxArray_real_T *omat,
                     emxArray_real_T *dsvm, const emxArray_cell_wrap_2 *NAME,
                     const emxArray_real_T *Lmat, emxArray_cell_wrap_2 *NAME2,
                     emxArray_real_T *Lmat2);
static void PWM_corr2(const emxArray_char_T *fn, double NUM, emxArray_real_T
                      *kmer, emxArray_real_T *dsvm, const emxArray_cell_wrap_2
                      *NAME, const emxArray_real_T *Lmat, const emxArray_char_T *
                      s);
static void b_binary_expand_op(emxArray_real_T *f, const emxArray_real_T
  *all_pwm, const emxArray_int8_T *ss_onehot, int i2, int i3, int i4);
static void getMOTIF(const emxArray_char_T *fn, emxArray_cell_wrap_3 *mat,
                     emxArray_cell_wrap_2 *names, emxArray_real_T *len);
static void getdenovomotif(const emxArray_char_T *filename, emxArray_cell_wrap_4
  *mat, emxArray_real_T *w);
static void j_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *mat,
  const emxArray_real_T *x);
static void letterconvert(const emxArray_char_T *s, emxArray_real_T *en);
static void minus(emxArray_real_T *maxnorm, const emxArray_real_T *minnorm);
static void ppmsim(emxArray_cell_wrap_4 *mot, const emxArray_real_T *lenvec,
                   double *ind, double *M);
static void process_motifs(const emxArray_char_T *dfn, const emxArray_char_T
  *lfn, const emxArray_char_T *memefn, const emxArray_char_T *ofn);
static void scoreseqkmer(const emxArray_real_T *ss, const emxArray_cell_wrap_4
  *Smat, double l, const emxArray_real_T *mat, const emxArray_real_T *dsvm,
  emxArray_real_T *omat);
static void seq2pv(const emxArray_char_T *sfn, const emxArray_char_T *wfn,
                   double input_l, emxArray_cell_wrap_0 *P, emxArray_cell_wrap_0
                   *V, emxArray_cell_wrap_0 *seqindmat, emxArray_cell_wrap_1
                   *seqout, emxArray_cell_wrap_2 *seq);
static void trim_pwm(emxArray_cell_wrap_4 *p, emxArray_real_T *info,
                     emxArray_real_T *len);

/* Function Definitions */
/*
 * function MAPTF(fn, mfn, s,ss, GC, pwm_prob, Smat, l_svm, k_svm, ofn, PWM, PWM2, pwm, lpwm, lab, LEN, LEN_2, shift, gkmprob, dsvm, names, len, a, b,rnum)
 */
static void MAPTF(const emxArray_char_T *s, const emxArray_real_T *ss, const
                  emxArray_real_T *pwm_prob, const emxArray_cell_wrap_4 *Smat,
                  double l_svm, const emxArray_char_T *ofn, const
                  emxArray_cell_wrap_3 *PWM2, emxArray_real_T *LEN, const
                  emxArray_real_T *LEN_2, const emxArray_real_T *shift,
                  emxArray_real_T *gkmprob, const emxArray_real_T *dsvm, const
                  emxArray_cell_wrap_2 *names, double a, double rnum)
{
  const cell_wrap_2 *names_data;
  cell_wrap_2 *NAME_data;
  const cell_wrap_3 *PWM2_data;
  emxArray_boolean_T *b_path;
  emxArray_cell_wrap_2 *NAME;
  emxArray_cell_wrap_2 *b_NAME;
  emxArray_int32_T *f;
  emxArray_int32_T *r;
  emxArray_real_T *C;
  emxArray_real_T *Cmat;
  emxArray_real_T *F;
  emxArray_real_T *MAT;
  emxArray_real_T *PATH;
  emxArray_real_T *TRANS;
  emxArray_real_T *b_F;
  emxArray_real_T *b_L;
  emxArray_real_T *b_MAT;
  emxArray_real_T *b_neg;
  emxArray_real_T *c_L;
  emxArray_real_T *ind;
  emxArray_real_T *mat;
  emxArray_real_T *neg;
  emxArray_real_T *path;
  emxArray_real_T *vec;
  const double *LEN_2_data;
  const double *dsvm_data;
  const double *pwm_prob_data;
  const double *shift_data;
  double L;
  double TRANS_tmp;
  double y;
  double *C_data;
  double *LEN_data;
  double *L_data;
  double *PATH_data;
  double *TRANS_data;
  double *gkmprob_data;
  double *ind_data;
  double *mat_data;
  double *neg_data;
  double *path_data;
  int b_i;
  int b_loop_ub;
  int i;
  int i1;
  int i2;
  int k;
  int loop_ub;
  int nx;
  int result;
  int *f_data;
  int *r1;
  signed char input_sizes_idx_1;
  signed char sizes_idx_1;
  bool empty_non_axis_sizes;
  bool *b_path_data;
  names_data = names->data;
  dsvm_data = dsvm->data;
  gkmprob_data = gkmprob->data;
  shift_data = shift->data;
  LEN_2_data = LEN_2->data;
  LEN_data = LEN->data;
  PWM2_data = PWM2->data;
  pwm_prob_data = pwm_prob->data;
  emxInit_real_T(&mat, 2);

  /* 'mapTF2_ls:153' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  /* 'mapTF2_ls:154' L = length(ss)-l_svm+1; */
  L = ((double)ss->size[1] - l_svm) + 1.0;

  /* 'mapTF2_ls:155' ACGT = 'ACGT'; */
  /* 'mapTF2_ls:156' lab = [lab; 0]; */
  /* 'mapTF2_ls:157' n = sum(LEN)+1; */
  y = blockedSummation(LEN, LEN->size[0]);

  /* 'mapTF2_ls:158' mat = zeros(n,L)-Inf; */
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

  /* 'mapTF2_ls:159' ind = zeros(n,L); */
  i = ind->size[0] * ind->size[1];
  ind->size[0] = (int)(y + 1.0);
  ind->size[1] = (int)L;
  emxEnsureCapacity_real_T(ind, i);
  ind_data = ind->data;
  for (i = 0; i < nx; i++) {
    ind_data[i] = 0.0;
  }

  emxInit_real_T(&TRANS, 2);

  /* 'mapTF2_ls:160' TRANS = zeros(n); */
  i = TRANS->size[0] * TRANS->size[1];
  TRANS->size[0] = (int)(y + 1.0);
  TRANS->size[1] = (int)(y + 1.0);
  emxEnsureCapacity_real_T(TRANS, i);
  TRANS_data = TRANS->data;
  b_loop_ub = (int)(y + 1.0) * (int)(y + 1.0);
  for (i = 0; i < b_loop_ub; i++) {
    TRANS_data[i] = 0.0;
  }

  emxInit_real_T(&PATH, 1);

  /* 'mapTF2_ls:161' LEN = [0;LEN]; */
  i = PATH->size[0];
  PATH->size[0] = LEN->size[0] + 1;
  emxEnsureCapacity_real_T(PATH, i);
  PATH_data = PATH->data;
  PATH_data[0] = 0.0;
  b_loop_ub = LEN->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    PATH_data[i + 1] = LEN_data[i];
  }

  i = LEN->size[0];
  LEN->size[0] = PATH->size[0];
  emxEnsureCapacity_real_T(LEN, i);
  LEN_data = LEN->data;
  b_loop_ub = PATH->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    LEN_data[i] = PATH_data[i];
  }

  emxInit_real_T(&C, 1);

  /* 'mapTF2_ls:162' C = cumsum(LEN); */
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

  /* 'mapTF2_ls:164' gkmprob(gkmprob>0.99) = 0.99; */
  nx = gkmprob->size[0];
  for (b_i = 0; b_i < nx; b_i++) {
    if (gkmprob_data[b_i] > 0.99) {
      gkmprob_data[b_i] = 0.99;
    }
  }

  emxInit_real_T(&neg, 1);

  /* 'mapTF2_ls:165' pos = gkmprob; */
  /* 'mapTF2_ls:166' neg = log(1-gkmprob); */
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

  /* 'mapTF2_ls:167' for i = 1:a */
  i = (int)a;
  emxInit_int32_T(&f, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:168' for j = 1:LEN(i+1)-1 */
    i1 = (int)(LEN_data[b_i + 1] - 1.0);
    for (k = 0; k < i1; k++) {
      /* 'mapTF2_ls:169' TRANS(C(i)+j,C(i)+j+1)=1; */
      TRANS_tmp = C_data[b_i] + ((double)k + 1.0);
      TRANS_data[((int)TRANS_tmp + TRANS->size[0] * ((int)(TRANS_tmp + 1.0) - 1))
        - 1] = 1.0;
    }

    /* 'mapTF2_ls:171' TRANS(C(i+1),C(1:end-1)+1) = 1; */
    if (1 > C->size[0] - 1) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = C->size[0] - 1;
    }

    i1 = f->size[0];
    f->size[0] = b_loop_ub;
    emxEnsureCapacity_int32_T(f, i1);
    f_data = f->data;
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      f_data[i1] = (int)(C_data[i1] + 1.0);
    }

    b_loop_ub = f->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      TRANS_data[((int)C_data[b_i + 1] + TRANS->size[0] * (f_data[i1] - 1)) - 1]
        = 1.0;
    }

    /* 'mapTF2_ls:172' TRANS(C(i+1),n)=1; */
    TRANS_data[((int)C_data[b_i + 1] + TRANS->size[0] * ((int)(y + 1.0) - 1)) -
      1] = 1.0;
  }

  /* 'mapTF2_ls:174' TRANS(n,C(1:end-1)+1) = 1; */
  if (1 > C->size[0] - 1) {
    b_loop_ub = 0;
  } else {
    b_loop_ub = C->size[0] - 1;
  }

  i1 = f->size[0];
  f->size[0] = b_loop_ub;
  emxEnsureCapacity_int32_T(f, i1);
  f_data = f->data;
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    f_data[i1] = (int)(C_data[i1] + 1.0);
  }

  b_loop_ub = f->size[0];
  for (i1 = 0; i1 < b_loop_ub; i1++) {
    TRANS_data[((int)(y + 1.0) + TRANS->size[0] * (f_data[i1] - 1)) - 1] = 1.0;
  }

  /* 'mapTF2_ls:175' TRANS(n,n) = 1; */
  TRANS_data[((int)(y + 1.0) + TRANS->size[0] * ((int)(y + 1.0) - 1)) - 1] = 1.0;

  /* 'mapTF2_ls:176' TRANS = log(TRANS); */
  nx = TRANS->size[0] * TRANS->size[1];
  for (k = 0; k < nx; k++) {
    TRANS_data[k] = log(TRANS_data[k]);
  }

  /* 'mapTF2_ls:177' for i = 1:a */
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:178' mat(C(i)+1,1) = pwm_prob(C(i)+1,1)+pos(1); */
    nx = (int)(C_data[b_i] + 1.0) - 1;
    mat_data[nx] = pwm_prob_data[nx] + gkmprob_data[0];
  }

  /* 'mapTF2_ls:180' mat(n,1) = neg(1); */
  mat_data[(int)(y + 1.0) - 1] = neg_data[0];

  /* 'mapTF2_ls:181' for i = 2:L */
  i1 = (int)(L + -1.0);
  for (b_i = 0; b_i < i1; b_i++) {
    /* 'mapTF2_ls:182' for j = 1:a */
    for (k = 0; k < i; k++) {
      /* 'mapTF2_ls:183' [mat(C(j)+1,i),ind(C(j)+1,i)] = max(mat(1:end,i-1)+TRANS(1:end,C(j)+1)+pwm_prob(C(j)+1,i)); */
      if (1 > mat->size[0]) {
        b_loop_ub = 0;
      } else {
        b_loop_ub = mat->size[0];
      }

      if (1 > TRANS->size[0]) {
        i2 = 0;
      } else {
        i2 = TRANS->size[0];
      }

      if (b_loop_ub == i2) {
        nx = (int)(C_data[k] + 1.0) - 1;
        TRANS_tmp = pwm_prob_data[nx + pwm_prob->size[0] * (b_i + 1)];
        i2 = PATH->size[0];
        PATH->size[0] = b_loop_ub;
        emxEnsureCapacity_real_T(PATH, i2);
        PATH_data = PATH->data;
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          PATH_data[i2] = (mat_data[i2 + mat->size[0] * b_i] + TRANS_data[i2 +
                           TRANS->size[0] * nx]) + TRANS_tmp;
        }

        c_maximum(PATH, &TRANS_tmp, &nx);
        mat_data[((int)(C_data[k] + 1.0) + mat->size[0] * (b_i + 1)) - 1] =
          TRANS_tmp;
      } else {
        nx = d_binary_expand_op(mat, b_loop_ub - 1, b_i, TRANS, i2 - 1, C, k,
          pwm_prob);
        mat_data = mat->data;
      }

      ind_data[((int)(C_data[k] + 1.0) + ind->size[0] * (b_i + 1)) - 1] = nx;

      /* 'mapTF2_ls:184' for jj = 1:LEN(j+1)-1 */
      i2 = (int)(LEN_data[k + 1] - 1.0);
      for (nx = 0; nx < i2; nx++) {
        /* 'mapTF2_ls:185' mat(C(j)+jj+1,i) = mat(C(j)+jj,i-1)+TRANS(C(j)+jj,C(j)+jj+1)+pwm_prob(C(j)+jj+1,i); */
        TRANS_tmp = C_data[k] + ((double)nx + 1.0);
        mat_data[((int)(TRANS_tmp + 1.0) + mat->size[0] * (b_i + 1)) - 1] =
          (mat_data[((int)TRANS_tmp + mat->size[0] * b_i) - 1] + TRANS_data
           [((int)TRANS_tmp + TRANS->size[0] * ((int)(TRANS_tmp + 1.0) - 1)) - 1])
          + pwm_prob_data[((int)(TRANS_tmp + 1.0) + pwm_prob->size[0] * (b_i + 1))
          - 1];

        /* 'mapTF2_ls:186' ind(C(j)+jj+1,i) = C(j)+jj; */
        ind_data[((int)(TRANS_tmp + 1.0) + ind->size[0] * (b_i + 1)) - 1] =
          TRANS_tmp;
      }
    }

    /* 'mapTF2_ls:189' [mat(n,i),ind(n,i)] = max(mat(1:end,i-1)+TRANS(1:end,end)+neg(i-1)); */
    if (1 > mat->size[0]) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = mat->size[0];
    }

    if (1 > TRANS->size[0]) {
      i2 = 0;
    } else {
      i2 = TRANS->size[0];
    }

    if (b_loop_ub == i2) {
      i2 = PATH->size[0];
      PATH->size[0] = b_loop_ub;
      emxEnsureCapacity_real_T(PATH, i2);
      PATH_data = PATH->data;
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        PATH_data[i2] = (mat_data[i2 + mat->size[0] * b_i] + TRANS_data[i2 +
                         TRANS->size[0] * (TRANS->size[1] - 1)]) + neg_data[b_i];
      }

      c_maximum(PATH, &TRANS_tmp, &nx);
      mat_data[((int)(y + 1.0) + mat->size[0] * (b_i + 1)) - 1] = TRANS_tmp;
    } else {
      nx = e_binary_expand_op(mat, b_loop_ub - 1, b_i, TRANS, i2 - 1, neg, y);
      mat_data = mat->data;
    }

    ind_data[((int)(y + 1.0) + ind->size[0] * (b_i + 1)) - 1] = nx;
  }

  emxFree_real_T(&TRANS);
  emxInit_real_T(&path, 1);

  /* 'mapTF2_ls:191' path = zeros(L,1); */
  i1 = path->size[0];
  path->size[0] = (int)L;
  emxEnsureCapacity_real_T(path, i1);
  path_data = path->data;
  for (i1 = 0; i1 < loop_ub; i1++) {
    path_data[i1] = 0.0;
  }

  /* 'mapTF2_ls:192' path(end) = n; */
  path_data[(int)L - 1] = y + 1.0;

  /*  for i = fliplr(1:L-1) */
  /* 'mapTF2_ls:194' for i = L-1:-1:1 */
  i1 = (int)-((-1.0 - (L - 1.0)) + 1.0);
  for (b_i = 0; b_i < i1; b_i++) {
    TRANS_tmp = (L - 1.0) + -(double)b_i;

    /* 'mapTF2_ls:195' path(i) = ind(path(i+1),i+1); */
    path_data[(int)TRANS_tmp - 1] = ind_data[((int)path_data[(int)(unsigned int)
      TRANS_tmp] + ind->size[0] * (int)(unsigned int)TRANS_tmp) - 1];
  }

  emxFree_real_T(&ind);

  /* 'mapTF2_ls:197' PATH = zeros(length(s),1); */
  i1 = PATH->size[0];
  PATH->size[0] = s->size[1];
  emxEnsureCapacity_real_T(PATH, i1);
  PATH_data = PATH->data;
  loop_ub = s->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    PATH_data[i1] = 0.0;
  }

  /* 'mapTF2_ls:198' C_2 = cumsum([0;LEN_2]); */
  i1 = neg->size[0];
  neg->size[0] = LEN_2->size[0] + 1;
  emxEnsureCapacity_real_T(neg, i1);
  neg_data = neg->data;
  neg_data[0] = 0.0;
  loop_ub = LEN_2->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    neg_data[i1 + 1] = LEN_2_data[i1];
  }

  if (neg->size[0] != 1) {
    i1 = neg->size[0];
    for (k = 0; k <= i1 - 2; k++) {
      neg_data[k + 1] += neg_data[k];
    }
  }

  /* 'mapTF2_ls:199' for i = 1:a */
  emxInit_real_T(&vec, 2);
  LEN_data = vec->data;
  emxInit_int32_T(&r, 2);
  emxInit_boolean_T(&b_path, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:200' f = find(path==C(i)+1); */
    loop_ub = path->size[0];
    i1 = b_path->size[0];
    b_path->size[0] = path->size[0];
    emxEnsureCapacity_boolean_T(b_path, i1);
    b_path_data = b_path->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_path_data[i1] = (path_data[i1] == C_data[b_i] + 1.0);
    }

    eml_find(b_path, f);
    f_data = f->data;

    /* 'mapTF2_ls:201' if ~isempty(f) */
    if (f->size[0] != 0) {
      /* 'mapTF2_ls:202' for j = 1:length(f) */
      i1 = f->size[0];
      for (k = 0; k < i1; k++) {
        /* 'mapTF2_ls:203' vec = f(j)+shift(i):f(j)+shift(i)+LEN_2(i)-1; */
        TRANS_tmp = (double)f_data[k] + shift_data[b_i];
        y = (TRANS_tmp + LEN_2_data[b_i]) - 1.0;
        if (rtIsNaN(TRANS_tmp) || rtIsNaN(y)) {
          i2 = vec->size[0] * vec->size[1];
          vec->size[0] = 1;
          vec->size[1] = 1;
          emxEnsureCapacity_real_T(vec, i2);
          LEN_data = vec->data;
          LEN_data[0] = rtNaN;
        } else if (y < TRANS_tmp) {
          vec->size[0] = 1;
          vec->size[1] = 0;
        } else if ((rtIsInf(TRANS_tmp) || rtIsInf(y)) && (TRANS_tmp == y)) {
          i2 = vec->size[0] * vec->size[1];
          vec->size[0] = 1;
          vec->size[1] = 1;
          emxEnsureCapacity_real_T(vec, i2);
          LEN_data = vec->data;
          LEN_data[0] = rtNaN;
        } else if (floor(TRANS_tmp) == TRANS_tmp) {
          i2 = vec->size[0] * vec->size[1];
          vec->size[0] = 1;
          loop_ub = (int)floor(y - TRANS_tmp);
          vec->size[1] = loop_ub + 1;
          emxEnsureCapacity_real_T(vec, i2);
          LEN_data = vec->data;
          for (i2 = 0; i2 <= loop_ub; i2++) {
            LEN_data[i2] = TRANS_tmp + (double)i2;
          }
        } else {
          eml_float_colon(TRANS_tmp, y, vec);
          LEN_data = vec->data;
        }

        /* 'mapTF2_ls:204' PATH(vec) = C_2(i)+1:C_2(i)+LEN_2(i); */
        i2 = r->size[0] * r->size[1];
        r->size[0] = 1;
        r->size[1] = vec->size[1];
        emxEnsureCapacity_int32_T(r, i2);
        r1 = r->data;
        loop_ub = vec->size[1];
        for (i2 = 0; i2 < loop_ub; i2++) {
          r1[i2] = (int)LEN_data[i2];
        }

        y = neg_data[b_i] + LEN_2_data[b_i];
        if (rtIsNaN(neg_data[b_i] + 1.0) || rtIsNaN(y)) {
          PATH_data[r1[0] - 1] = rtNaN;
        } else if (!(y < neg_data[b_i] + 1.0)) {
          if ((rtIsInf(neg_data[b_i] + 1.0) || rtIsInf(y)) && (neg_data[b_i] +
               1.0 == y)) {
            PATH_data[r1[0] - 1] = rtNaN;
          } else if (floor(neg_data[b_i] + 1.0) == neg_data[b_i] + 1.0) {
            loop_ub = (int)floor(y - (neg_data[b_i] + 1.0));
            for (i2 = 0; i2 <= loop_ub; i2++) {
              PATH_data[r1[i2] - 1] = (neg_data[b_i] + 1.0) + (double)i2;
            }
          } else {
            eml_float_colon(neg_data[b_i] + 1.0, y, vec);
            LEN_data = vec->data;
            loop_ub = vec->size[1];
            for (i2 = 0; i2 < loop_ub; i2++) {
              PATH_data[r1[i2] - 1] = LEN_data[i2];
            }
          }
        }
      }
    }
  }

  emxFree_int32_T(&r);
  emxFree_real_T(&vec);

  /* 'mapTF2_ls:208' N =sum(LEN_2)+1; */
  y = blockedSummation(LEN_2, LEN_2->size[0]);

  /* 'mapTF2_ls:209' PATH(PATH==0) = N; */
  nx = PATH->size[0];
  for (b_i = 0; b_i < nx; b_i++) {
    if (PATH_data[b_i] == 0.0) {
      PATH_data[b_i] = y + 1.0;
    }
  }

  emxInit_real_T(&MAT, 2);
  ind_data = MAT->data;

  /* 'mapTF2_ls:210' MAT = []; */
  MAT->size[0] = 0;
  MAT->size[1] = 4;

  /* 'mapTF2_ls:211' lab = []; */
  neg->size[0] = 0;

  /* 'mapTF2_ls:212' for i = 1:a */
  emxInit_real_T(&b_MAT, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:213' MAT = [MAT;PWM2{i}]; */
    loop_ub = MAT->size[0];
    i1 = b_MAT->size[0] * b_MAT->size[1];
    b_MAT->size[0] = MAT->size[0] + PWM2_data[b_i].f1->size[0];
    b_MAT->size[1] = 4;
    emxEnsureCapacity_real_T(b_MAT, i1);
    LEN_data = b_MAT->data;
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        LEN_data[i2 + b_MAT->size[0] * i1] = ind_data[i2 + MAT->size[0] * i1];
      }
    }

    loop_ub = PWM2_data[b_i].f1->size[0];
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        LEN_data[(i2 + MAT->size[0]) + b_MAT->size[0] * i1] = PWM2_data[b_i].
          f1->data[i2 + PWM2_data[b_i].f1->size[0] * i1];
      }
    }

    i1 = MAT->size[0] * MAT->size[1];
    MAT->size[0] = b_MAT->size[0];
    MAT->size[1] = 4;
    emxEnsureCapacity_real_T(MAT, i1);
    ind_data = MAT->data;
    loop_ub = b_MAT->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      ind_data[i1] = LEN_data[i1];
    }

    /* 'mapTF2_ls:214' lab = [lab;ones(LEN_2(i),1)*i]; */
    i1 = neg->size[0];
    loop_ub = (int)LEN_2_data[b_i];
    i2 = neg->size[0];
    neg->size[0] += loop_ub;
    emxEnsureCapacity_real_T(neg, i2);
    neg_data = neg->data;
    for (i2 = 0; i2 < loop_ub; i2++) {
      neg_data[i1 + i2] = (unsigned int)(b_i + 1);
    }
  }

  /* 'mapTF2_ls:216' lab = [lab;0]; */
  i = neg->size[0];
  i1 = neg->size[0];
  neg->size[0]++;
  emxEnsureCapacity_real_T(neg, i1);
  neg_data = neg->data;
  neg_data[i] = 0.0;

  /* 'mapTF2_ls:217' MAT = [MAT;ones(1,4)/4]; */
  i = b_MAT->size[0] * b_MAT->size[1];
  b_MAT->size[0] = MAT->size[0] + 1;
  b_MAT->size[1] = 4;
  emxEnsureCapacity_real_T(b_MAT, i);
  LEN_data = b_MAT->data;
  loop_ub = MAT->size[0];
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      LEN_data[i1 + b_MAT->size[0] * i] = ind_data[i1 + MAT->size[0] * i];
    }

    LEN_data[MAT->size[0] + b_MAT->size[0] * i] = 0.25;
  }

  i = MAT->size[0] * MAT->size[1];
  MAT->size[0] = b_MAT->size[0];
  MAT->size[1] = 4;
  emxEnsureCapacity_real_T(MAT, i);
  ind_data = MAT->data;
  loop_ub = b_MAT->size[0] * 4;
  for (i = 0; i < loop_ub; i++) {
    ind_data[i] = LEN_data[i];
  }

  /* 'mapTF2_ls:218' n = sum(LEN_2)+1; */
  y = blockedSummation(LEN_2, LEN_2->size[0]);

  /* 'mapTF2_ls:219' PWM = ones(L,4)/4; */
  /* 'mapTF2_ls:220' PWM_alt = ones(L,4)/4; */
  /* 'mapTF2_ls:222' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  /* 'mapTF2_ls:224' I = lab(PATH); */
  /* 'mapTF2_ls:225' Cmat = [I PATH ones(length(I),4)/4]; */
  if (PATH->size[0] != 0) {
    result = PATH->size[0];
  } else {
    result = 0;
  }

  empty_non_axis_sizes = (result == 0);
  if (empty_non_axis_sizes || (PATH->size[0] != 0)) {
    nx = 1;
  } else {
    nx = 0;
  }

  if (empty_non_axis_sizes || (PATH->size[0] != 0)) {
    input_sizes_idx_1 = 1;
  } else {
    input_sizes_idx_1 = 0;
  }

  if (empty_non_axis_sizes || (PATH->size[0] != 0)) {
    sizes_idx_1 = 4;
  } else {
    sizes_idx_1 = 0;
  }

  emxInit_real_T(&b_neg, 1);
  i = b_neg->size[0];
  b_neg->size[0] = PATH->size[0];
  emxEnsureCapacity_real_T(b_neg, i);
  LEN_data = b_neg->data;
  loop_ub = PATH->size[0];
  for (i = 0; i < loop_ub; i++) {
    LEN_data[i] = neg_data[(int)PATH_data[i] - 1];
  }

  emxInit_real_T(&Cmat, 2);
  i = Cmat->size[0] * Cmat->size[1];
  Cmat->size[0] = result;
  Cmat->size[1] = (nx + input_sizes_idx_1) + sizes_idx_1;
  emxEnsureCapacity_real_T(Cmat, i);
  TRANS_data = Cmat->data;
  for (i = 0; i < nx; i++) {
    for (i1 = 0; i1 < result; i1++) {
      TRANS_data[i1] = LEN_data[i1];
    }
  }

  emxFree_real_T(&b_neg);
  loop_ub = input_sizes_idx_1;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < result; i1++) {
      TRANS_data[i1 + Cmat->size[0] * nx] = PATH_data[i1];
    }
  }

  loop_ub = sizes_idx_1;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < result; i1++) {
      TRANS_data[i1 + Cmat->size[0] * ((i + nx) + input_sizes_idx_1)] = 0.25;
    }
  }

  /* 'mapTF2_ls:226' for i = 1:L+l_svm-1 */
  i = (int)((L + l_svm) - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:227' if PATH(i)~=n */
    if (PATH_data[b_i] != y + 1.0) {
      /* 'mapTF2_ls:228' Cmat(i,3:6) = MAT(PATH(i),:); */
      nx = (int)PATH_data[b_i] - 1;
      TRANS_data[b_i + Cmat->size[0] * 2] = ind_data[nx];
      TRANS_data[b_i + Cmat->size[0] * 3] = ind_data[nx + MAT->size[0]];
      TRANS_data[b_i + Cmat->size[0] * 4] = ind_data[nx + MAT->size[0] * 2];
      TRANS_data[b_i + Cmat->size[0] * 5] = ind_data[nx + MAT->size[0] * 3];
    } else {
      /* 'mapTF2_ls:229' else */
      /* 'mapTF2_ls:230' Cmat(i,2) = 0; */
      TRANS_data[b_i + Cmat->size[0]] = 0.0;
    }
  }

  /* 'mapTF2_ls:233' Cmat(:,3:6) = (Cmat(:,3:6)+0.0001)/(1.0001); */
  nx = Cmat->size[0] - 1;
  i = b_MAT->size[0] * b_MAT->size[1];
  b_MAT->size[0] = Cmat->size[0];
  b_MAT->size[1] = 4;
  emxEnsureCapacity_real_T(b_MAT, i);
  LEN_data = b_MAT->data;
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 <= nx; i1++) {
      LEN_data[i1 + b_MAT->size[0] * i] = (TRANS_data[i1 + Cmat->size[0] * (i +
        2)] + 0.0001) / 1.0001;
    }
  }

  loop_ub = b_MAT->size[0];
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      TRANS_data[i1 + Cmat->size[0] * (i + 2)] = LEN_data[i1 + b_MAT->size[0] *
        i];
    }
  }

  emxInit_real_T(&b_L, 2);
  L_data = b_L->data;

  /* 'mapTF2_ls:234' LEN(1) = []; */
  /* 'mapTF2_ls:235' L = []; */
  b_L->size[0] = 0;
  b_L->size[1] = 0;

  /* 'mapTF2_ls:236' for i = 1:length(PWM2) */
  i = PWM2->size[0];
  emxInit_real_T(&F, 2);
  TRANS_data = F->data;
  emxInit_real_T(&b_F, 2);
  emxInit_real_T(&c_L, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:237' f = find(path==C(i)+1); */
    loop_ub = path->size[0];
    i1 = b_path->size[0];
    b_path->size[0] = path->size[0];
    emxEnsureCapacity_boolean_T(b_path, i1);
    b_path_data = b_path->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_path_data[i1] = (path_data[i1] == C_data[b_i] + 1.0);
    }

    eml_find(b_path, f);
    f_data = f->data;

    /* 'mapTF2_ls:238' F = []; */
    F->size[0] = 0;
    F->size[1] = 3;

    /* 'mapTF2_ls:239' a = 0; */
    a = 0.0;

    /* 'mapTF2_ls:240' if ~isempty(f) */
    if (f->size[0] != 0) {
      /* 'mapTF2_ls:241' for j = 1:length(f) */
      i1 = f->size[0];
      for (k = 0; k < i1; k++) {
        /* 'mapTF2_ls:242' if f(j)+shift(i)+LEN_2(i)< length(ss)-l_svm+1 && f(j)+shift(i) > l_svm-1 */
        TRANS_tmp = (double)f_data[k] + shift_data[b_i];
        if ((TRANS_tmp + LEN_2_data[b_i] < ((double)ss->size[1] - l_svm) + 1.0) &&
            (TRANS_tmp > l_svm - 1.0)) {
          /* 'mapTF2_ls:243' a = a+1; */
          a++;

          /* 'mapTF2_ls:244' F = [F; f(j)+shift(i) f(j)+shift(i)+LEN_2(i)-1 mean(gkmprob(f(j):f(j)+2*shift(i)+LEN_2(i)-l_svm))]; */
          i2 = f_data[k];
          TRANS_tmp = (((double)i2 + 2.0 * shift_data[b_i]) + LEN_2_data[b_i]) -
            l_svm;
          if (i2 > TRANS_tmp) {
            nx = -1;
            result = -1;
          } else {
            nx = f_data[k] - 2;
            result = (int)TRANS_tmp - 1;
          }

          loop_ub = result - nx;
          result = PATH->size[0];
          PATH->size[0] = loop_ub;
          emxEnsureCapacity_real_T(PATH, result);
          PATH_data = PATH->data;
          for (result = 0; result < loop_ub; result++) {
            PATH_data[result] = gkmprob_data[(nx + result) + 1];
          }

          b_loop_ub = F->size[0];
          nx = b_F->size[0] * b_F->size[1];
          b_F->size[0] = F->size[0] + 1;
          b_F->size[1] = 3;
          emxEnsureCapacity_real_T(b_F, nx);
          LEN_data = b_F->data;
          for (nx = 0; nx < 3; nx++) {
            for (result = 0; result < b_loop_ub; result++) {
              LEN_data[result + b_F->size[0] * nx] = TRANS_data[result + F->
                size[0] * nx];
            }
          }

          TRANS_tmp = (double)i2 + shift_data[b_i];
          LEN_data[F->size[0]] = TRANS_tmp;
          LEN_data[F->size[0] + b_F->size[0]] = (TRANS_tmp + LEN_2_data[b_i]) -
            1.0;
          LEN_data[F->size[0] + b_F->size[0] * 2] = blockedSummation(PATH,
            loop_ub) / (double)loop_ub;
          i2 = F->size[0] * F->size[1];
          F->size[0] = b_F->size[0];
          F->size[1] = 3;
          emxEnsureCapacity_real_T(F, i2);
          TRANS_data = F->data;
          loop_ub = b_F->size[0] * 3;
          for (i2 = 0; i2 < loop_ub; i2++) {
            TRANS_data[i2] = LEN_data[i2];
          }
        }
      }

      /* 'mapTF2_ls:247' if a > 0 */
      if (a > 0.0) {
        /* 'mapTF2_ls:248' L = [L;F ones(a,1)*i]; */
        if (F->size[0] != 0) {
          result = F->size[0];
        } else if ((int)a != 0) {
          result = (int)a;
        } else {
          result = 0;
          if ((int)a > 0) {
            result = (int)a;
          }
        }

        empty_non_axis_sizes = (result == 0);
        if (empty_non_axis_sizes || (F->size[0] != 0)) {
          input_sizes_idx_1 = 3;
        } else {
          input_sizes_idx_1 = 0;
        }

        if (empty_non_axis_sizes || ((int)a != 0)) {
          sizes_idx_1 = 1;
        } else {
          sizes_idx_1 = 0;
        }

        nx = result;
        i1 = input_sizes_idx_1 + sizes_idx_1;
        if ((b_L->size[0] != 0) && (b_L->size[1] != 0)) {
          k = b_L->size[1];
        } else if ((result != 0) && (i1 != 0)) {
          k = i1;
        } else {
          k = b_L->size[1];
          if (i1 > b_L->size[1]) {
            k = i1;
          }
        }

        empty_non_axis_sizes = (k == 0);
        if (empty_non_axis_sizes || ((b_L->size[0] != 0) && (b_L->size[1] != 0)))
        {
          loop_ub = b_L->size[0];
        } else {
          loop_ub = 0;
        }

        if ((!empty_non_axis_sizes) && ((result == 0) || (i1 == 0))) {
          nx = 0;
        }

        i2 = mat->size[0] * mat->size[1];
        mat->size[0] = result;
        mat->size[1] = i1;
        emxEnsureCapacity_real_T(mat, i2);
        mat_data = mat->data;
        b_loop_ub = input_sizes_idx_1;
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          for (i2 = 0; i2 < result; i2++) {
            mat_data[i2 + mat->size[0] * i1] = TRANS_data[i2 + result * i1];
          }
        }

        b_loop_ub = sizes_idx_1;
        for (i1 = 0; i1 < b_loop_ub; i1++) {
          for (i2 = 0; i2 < result; i2++) {
            mat_data[i2 + mat->size[0] * input_sizes_idx_1] = (unsigned int)(b_i
              + 1);
          }
        }

        i1 = c_L->size[0] * c_L->size[1];
        c_L->size[0] = loop_ub + nx;
        c_L->size[1] = k;
        emxEnsureCapacity_real_T(c_L, i1);
        LEN_data = c_L->data;
        for (i1 = 0; i1 < k; i1++) {
          for (i2 = 0; i2 < loop_ub; i2++) {
            LEN_data[i2 + c_L->size[0] * i1] = L_data[i2 + loop_ub * i1];
          }
        }

        for (i1 = 0; i1 < k; i1++) {
          for (i2 = 0; i2 < nx; i2++) {
            LEN_data[(i2 + loop_ub) + c_L->size[0] * i1] = mat_data[i2 + nx * i1];
          }
        }

        i1 = b_L->size[0] * b_L->size[1];
        b_L->size[0] = c_L->size[0];
        b_L->size[1] = c_L->size[1];
        emxEnsureCapacity_real_T(b_L, i1);
        L_data = b_L->data;
        loop_ub = c_L->size[0] * c_L->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          L_data[i1] = LEN_data[i1];
        }
      }
    }
  }

  emxFree_real_T(&b_F);
  emxFree_boolean_T(&b_path);
  emxFree_real_T(&path);
  emxFree_real_T(&C);
  emxFree_real_T(&mat);
  emxInit_cell_wrap_2(&NAME);

  /* 'mapTF2_ls:252' NAME = cell(numel(L)/4,1); */
  nx = (int)((double)(b_L->size[0] * b_L->size[1]) / 4.0);
  i = NAME->size[0];
  NAME->size[0] = (int)((double)(b_L->size[0] * b_L->size[1]) / 4.0);
  emxEnsureCapacity_cell_wrap_2(NAME, i);
  NAME_data = NAME->data;
  for (i = 0; i < nx; i++) {
    NAME_data[i].f1->size[0] = 1;
    NAME_data[i].f1->size[1] = 0;
  }

  /* 'mapTF2_ls:253' NAME = coder.nullcopy(NAME); */
  /* 'mapTF2_ls:254' Lmat = zeros(numel(L)/4,4); */
  i = MAT->size[0] * MAT->size[1];
  MAT->size[0] = (int)((double)(b_L->size[0] * b_L->size[1]) / 4.0);
  MAT->size[1] = 4;
  emxEnsureCapacity_real_T(MAT, i);
  ind_data = MAT->data;
  loop_ub = (int)((double)(b_L->size[0] * b_L->size[1]) / 4.0) << 2;
  for (i = 0; i < loop_ub; i++) {
    ind_data[i] = 0.0;
  }

  /* 'mapTF2_ls:255' if ~isempty(L) */
  if ((b_L->size[0] != 0) && (b_L->size[1] != 0)) {
    /* 'mapTF2_ls:256' [~,b] = sort(L(:,1)); */
    loop_ub = b_L->size[0];
    i = neg->size[0];
    neg->size[0] = b_L->size[0];
    emxEnsureCapacity_real_T(neg, i);
    neg_data = neg->data;
    for (i = 0; i < loop_ub; i++) {
      neg_data[i] = L_data[i];
    }

    sort(neg, f);
    f_data = f->data;

    /* 'mapTF2_ls:257' L = L(b,:); */
    nx = b_L->size[1] - 1;
    i = c_L->size[0] * c_L->size[1];
    c_L->size[0] = f->size[0];
    c_L->size[1] = b_L->size[1];
    emxEnsureCapacity_real_T(c_L, i);
    LEN_data = c_L->data;
    for (i = 0; i <= nx; i++) {
      loop_ub = f->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        LEN_data[i1 + c_L->size[0] * i] = L_data[(f_data[i1] + b_L->size[0] * i)
          - 1];
      }
    }

    i = b_L->size[0] * b_L->size[1];
    b_L->size[0] = c_L->size[0];
    b_L->size[1] = c_L->size[1];
    emxEnsureCapacity_real_T(b_L, i);
    L_data = b_L->data;
    loop_ub = c_L->size[0] * c_L->size[1];
    for (i = 0; i < loop_ub; i++) {
      L_data[i] = LEN_data[i];
    }

    /* 'mapTF2_ls:258' for i = 1:numel(L)/4 */
    i = (int)((double)(b_L->size[0] * b_L->size[1]) / 4.0);
    for (b_i = 0; b_i < i; b_i++) {
      /* 'mapTF2_ls:259' NAME{i} = names{L(i,4)}; */
      i1 = NAME_data[b_i].f1->size[0] * NAME_data[b_i].f1->size[1];
      NAME_data[b_i].f1->size[0] = 1;
      NAME_data[b_i].f1->size[1] = names_data[(int)L_data[b_i + b_L->size[0] * 3]
        - 1].f1->size[1];
      emxEnsureCapacity_char_T(NAME_data[b_i].f1, i1);
      loop_ub = names_data[(int)L_data[b_i + b_L->size[0] * 3] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        NAME_data[b_i].f1->data[i1] = names_data[(int)L_data[b_i + b_L->size[0] *
          3] - 1].f1->data[i1];
      }

      /* 'mapTF2_ls:260' Lmat(i,:) = [L(i,4) L(i,1) L(i,2) L(i,3)]; */
      ind_data[b_i] = L_data[b_i + b_L->size[0] * 3];
      ind_data[b_i + MAT->size[0]] = L_data[b_i];
      ind_data[b_i + MAT->size[0] * 2] = L_data[b_i + b_L->size[0]];
      ind_data[b_i + MAT->size[0] * 3] = L_data[b_i + b_L->size[0] * 2];
    }
  }

  emxFree_real_T(&c_L);
  emxFree_real_T(&b_L);
  emxFree_int32_T(&f);
  emxFree_real_T(&neg);

  /* 'mapTF2_ls:263' kmer = scoreseqkmer(fn, rnum, ss, Smat, l_svm,k_svm, ofn,Cmat, dsvm); */
  scoreseqkmer(ss, Smat, l_svm, Cmat, dsvm, F);

  /* 'mapTF2_ls:264' [omat, NAME, Lmat] = PWM_corr(ofn, rnum, mfn,kmer,Cmat, dsvm, NAME, Lmat,PWM2); */
  i = PATH->size[0];
  PATH->size[0] = dsvm->size[0];
  emxEnsureCapacity_real_T(PATH, i);
  PATH_data = PATH->data;
  loop_ub = dsvm->size[0] - 1;
  for (i = 0; i <= loop_ub; i++) {
    PATH_data[i] = dsvm_data[i];
  }

  emxInit_cell_wrap_2(&b_NAME);
  PWM_corr(F, Cmat, PATH, NAME, MAT, b_NAME, b_MAT);

  /* 'mapTF2_ls:265' kmer = scoreseqkmer(fn, rnum, ss, Smat, l_svm, k_svm,ofn,omat, dsvm); */
  scoreseqkmer(ss, Smat, l_svm, Cmat, dsvm, F);

  /* 'mapTF2_ls:266' PWM_corr2(ofn, rnum, kmer, dsvm, NAME, Lmat, s) */
  i = PATH->size[0];
  PATH->size[0] = dsvm->size[0];
  emxEnsureCapacity_real_T(PATH, i);
  PATH_data = PATH->data;
  loop_ub = dsvm->size[0] - 1;
  emxFree_cell_wrap_2(&NAME);
  emxFree_real_T(&Cmat);
  emxFree_real_T(&MAT);
  for (i = 0; i <= loop_ub; i++) {
    PATH_data[i] = dsvm_data[i];
  }

  PWM_corr2(ofn, rnum, F, PATH, b_NAME, b_MAT, s);
  emxFree_real_T(&b_MAT);
  emxFree_cell_wrap_2(&b_NAME);
  emxFree_real_T(&F);
  emxFree_real_T(&PATH);
}

/*
 * function [omat, NAME2, Lmat2] = PWM_corr(fn,NUM,mfn,kmer,omat, dsvm, NAME, Lmat, p)
 */
static void PWM_corr(emxArray_real_T *kmer, emxArray_real_T *omat,
                     emxArray_real_T *dsvm, const emxArray_cell_wrap_2 *NAME,
                     const emxArray_real_T *Lmat, emxArray_cell_wrap_2 *NAME2,
                     emxArray_real_T *Lmat2)
{
  const cell_wrap_2 *NAME_data;
  cell_wrap_2 *NAME2_data;
  emxArray_boolean_T *b_LEN;
  emxArray_int32_T *ib;
  emxArray_int32_T *iidx;
  emxArray_real_T *LEN;
  emxArray_real_T *b_dsvm;
  emxArray_real_T *b_kmer;
  emxArray_real_T *f;
  emxArray_real_T *y;
  emxArray_uint32_T *n;
  const double *Lmat_data;
  double b_Lmat;
  double b_x;
  double x;
  double *LEN_data;
  double *Lmat2_data;
  double *dsvm_data;
  double *f_data;
  double *kmer_data;
  double *omat_data;
  double *y_data;
  int b_i;
  int dsvm_tmp;
  int i;
  int i1;
  int i2;
  int loop_ub;
  int *ib_data;
  int *iidx_data;
  unsigned int *n_data;
  bool *b_LEN_data;
  Lmat_data = Lmat->data;
  NAME_data = NAME->data;
  dsvm_data = dsvm->data;
  omat_data = omat->data;
  kmer_data = kmer->data;
  emxInit_real_T(&LEN, 1);

  /* 'mapTF2_ls:375' [coords, ind] = sort(kmer(:,1)); */
  loop_ub = kmer->size[0];
  i = LEN->size[0];
  LEN->size[0] = kmer->size[0];
  emxEnsureCapacity_real_T(LEN, i);
  LEN_data = LEN->data;
  for (i = 0; i < loop_ub; i++) {
    LEN_data[i] = kmer_data[i];
  }

  emxInit_int32_T(&iidx, 1);
  emxInit_real_T(&b_kmer, 2);
  sort(LEN, iidx);
  iidx_data = iidx->data;
  LEN_data = LEN->data;

  /* 'mapTF2_ls:376' kmer = kmer(ind, :); */
  i = b_kmer->size[0] * b_kmer->size[1];
  b_kmer->size[0] = iidx->size[0];
  b_kmer->size[1] = 3;
  emxEnsureCapacity_real_T(b_kmer, i);
  Lmat2_data = b_kmer->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      Lmat2_data[i1 + b_kmer->size[0] * i] = kmer_data[(iidx_data[i1] +
        kmer->size[0] * i) - 1];
    }
  }

  i = kmer->size[0] * kmer->size[1];
  kmer->size[0] = b_kmer->size[0];
  kmer->size[1] = 3;
  emxEnsureCapacity_real_T(kmer, i);
  kmer_data = kmer->data;
  loop_ub = b_kmer->size[0] * 3;
  for (i = 0; i < loop_ub; i++) {
    kmer_data[i] = Lmat2_data[i];
  }

  emxFree_real_T(&b_kmer);
  emxInit_real_T(&b_dsvm, 1);

  /* 'mapTF2_ls:377' dsvm = dsvm(ind, :); */
  i = b_dsvm->size[0];
  b_dsvm->size[0] = iidx->size[0];
  emxEnsureCapacity_real_T(b_dsvm, i);
  Lmat2_data = b_dsvm->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    Lmat2_data[i] = dsvm_data[iidx_data[i] - 1];
  }

  i = dsvm->size[0];
  dsvm->size[0] = b_dsvm->size[0];
  emxEnsureCapacity_real_T(dsvm, i);
  dsvm_data = dsvm->data;
  loop_ub = b_dsvm->size[0];
  for (i = 0; i < loop_ub; i++) {
    dsvm_data[i] = Lmat2_data[i];
  }

  emxInit_uint32_T(&n);
  n_data = n->data;

  /* 'mapTF2_ls:378' L = length(NAME); */
  /* 'mapTF2_ls:379' ind = ones(L,1); */
  /* 'mapTF2_ls:380' n = []; */
  n->size[0] = 0;

  /* 'mapTF2_ls:381' for i = 1:L */
  i = NAME->size[0];
  emxInit_real_T(&f, 1);
  emxInit_int32_T(&ib, 1);
  emxInit_real_T(&y, 1);
  emxInit_boolean_T(&b_LEN, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:382' F = find(coords > Lmat(i,2)); */
    /* 'mapTF2_ls:383' FF = find(coords <= Lmat(i,3)); */
    /* 'mapTF2_ls:384' f = intersect(F,FF); */
    b_Lmat = Lmat_data[b_i + Lmat->size[0]];
    loop_ub = LEN->size[0];
    i1 = b_LEN->size[0];
    b_LEN->size[0] = LEN->size[0];
    emxEnsureCapacity_boolean_T(b_LEN, i1);
    b_LEN_data = b_LEN->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_LEN_data[i1] = (LEN_data[i1] > b_Lmat);
    }

    eml_find(b_LEN, iidx);
    iidx_data = iidx->data;
    b_Lmat = Lmat_data[b_i + Lmat->size[0] * 2];
    loop_ub = LEN->size[0];
    i1 = b_LEN->size[0];
    b_LEN->size[0] = LEN->size[0];
    emxEnsureCapacity_boolean_T(b_LEN, i1);
    b_LEN_data = b_LEN->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_LEN_data[i1] = (LEN_data[i1] <= b_Lmat);
    }

    eml_find(b_LEN, ib);
    ib_data = ib->data;
    i1 = b_dsvm->size[0];
    b_dsvm->size[0] = iidx->size[0];
    emxEnsureCapacity_real_T(b_dsvm, i1);
    Lmat2_data = b_dsvm->data;
    loop_ub = iidx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      Lmat2_data[i1] = iidx_data[i1];
    }

    i1 = y->size[0];
    y->size[0] = ib->size[0];
    emxEnsureCapacity_real_T(y, i1);
    y_data = y->data;
    loop_ub = ib->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      y_data[i1] = ib_data[i1];
    }

    do_vectors(b_dsvm, y, f, iidx, ib);
    f_data = f->data;

    /* 'mapTF2_ls:385' a1 = ip(kmer(f,end),dsvm(f,end)); */
    loop_ub = f->size[0];
    i1 = b_dsvm->size[0];
    b_dsvm->size[0] = f->size[0];
    emxEnsureCapacity_real_T(b_dsvm, i1);
    Lmat2_data = b_dsvm->data;
    i1 = y->size[0];
    y->size[0] = f->size[0];
    emxEnsureCapacity_real_T(y, i1);
    y_data = y->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      dsvm_tmp = (int)f_data[i1] - 1;
      Lmat2_data[i1] = kmer_data[dsvm_tmp + kmer->size[0] * 2];
      y_data[i1] = dsvm_data[dsvm_tmp];
    }

    /* 'mapTF2_ls:426' c = x'*y/sqrt(x'*x)/sqrt(y'*y); */
    if (f->size[0] < 1) {
      b_Lmat = 0.0;
      x = 0.0;
      b_x = 0.0;
    } else {
      b_Lmat = cblas_ddot((blasint)f->size[0], &Lmat2_data[0], (blasint)1,
                          &y_data[0], (blasint)1);
      x = cblas_ddot((blasint)f->size[0], &Lmat2_data[0], (blasint)1,
                     &Lmat2_data[0], (blasint)1);
      b_x = cblas_ddot((blasint)f->size[0], &y_data[0], (blasint)1, &y_data[0],
                       (blasint)1);
    }

    x = sqrt(x);
    b_x = sqrt(b_x);

    /* 'mapTF2_ls:386' if a1 > 0.6 */
    if (b_Lmat / x / b_x > 0.6) {
      /* 'mapTF2_ls:387' n = [n;i]; */
      i1 = n->size[0];
      dsvm_tmp = n->size[0];
      n->size[0]++;
      emxEnsureCapacity_uint32_T(n, dsvm_tmp);
      n_data = n->data;
      n_data[i1] = (unsigned int)(b_i + 1);
    } else {
      /* 'mapTF2_ls:388' else */
      /* 'mapTF2_ls:389' ind(i) = 0; */
      /* 'mapTF2_ls:390' omat(Lmat(i,2):Lmat(i,3),1:2) = 0; */
      b_Lmat = Lmat_data[b_i + Lmat->size[0]];
      x = Lmat_data[b_i + Lmat->size[0] * 2];
      if (b_Lmat > x) {
        i1 = 0;
        dsvm_tmp = 0;
      } else {
        i1 = (int)b_Lmat - 1;
        dsvm_tmp = (int)x;
      }

      loop_ub = dsvm_tmp - i1;
      for (dsvm_tmp = 0; dsvm_tmp < 2; dsvm_tmp++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          omat_data[(i1 + i2) + omat->size[0] * dsvm_tmp] = 0.0;
        }
      }

      /* 'mapTF2_ls:391' omat(Lmat(i,2):Lmat(i,3),3:6) = 0.25; */
      b_Lmat = Lmat_data[b_i + Lmat->size[0]];
      x = Lmat_data[b_i + Lmat->size[0] * 2];
      if (b_Lmat > x) {
        i1 = 0;
        dsvm_tmp = 0;
      } else {
        i1 = (int)b_Lmat - 1;
        dsvm_tmp = (int)x;
      }

      loop_ub = dsvm_tmp - i1;
      for (dsvm_tmp = 0; dsvm_tmp < 4; dsvm_tmp++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          omat_data[(i1 + i2) + omat->size[0] * (dsvm_tmp + 2)] = 0.25;
        }
      }
    }
  }

  emxFree_boolean_T(&b_LEN);
  emxFree_real_T(&b_dsvm);
  emxFree_real_T(&y);
  emxFree_int32_T(&ib);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&LEN);
  emxFree_real_T(&f);

  /* 'mapTF2_ls:394' Lmat2 = Lmat(n,:); */
  i = Lmat2->size[0] * Lmat2->size[1];
  Lmat2->size[0] = n->size[0];
  Lmat2->size[1] = 4;
  emxEnsureCapacity_real_T(Lmat2, i);
  Lmat2_data = Lmat2->data;
  loop_ub = n->size[0];
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      Lmat2_data[i1 + Lmat2->size[0] * i] = Lmat_data[((int)n_data[i1] +
        Lmat->size[0] * i) - 1];
    }
  }

  /*  NAME2 = NAME(n); */
  /* 'mapTF2_ls:396' NAME2_new_length = length(n); */
  /* 'mapTF2_ls:397' NAME2_new = cell(NAME2_new_length,1); */
  /* 'mapTF2_ls:398' for idx=1:NAME2_new_length */
  i = n->size[0];
  i1 = NAME2->size[0];
  NAME2->size[0] = n->size[0];
  emxEnsureCapacity_cell_wrap_2(NAME2, i1);
  NAME2_data = NAME2->data;
  for (dsvm_tmp = 0; dsvm_tmp < i; dsvm_tmp++) {
    /* 'mapTF2_ls:399' NAME2_new{idx} = NAME{n(idx)}; */
    i1 = NAME2_data[dsvm_tmp].f1->size[0] * NAME2_data[dsvm_tmp].f1->size[1];
    NAME2_data[dsvm_tmp].f1->size[0] = 1;
    NAME2_data[dsvm_tmp].f1->size[1] = NAME_data[(int)n_data[dsvm_tmp] - 1]
      .f1->size[1];
    emxEnsureCapacity_char_T(NAME2_data[dsvm_tmp].f1, i1);
    loop_ub = NAME_data[(int)n_data[dsvm_tmp] - 1].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      NAME2_data[dsvm_tmp].f1->data[i1] = NAME_data[(int)n_data[dsvm_tmp] - 1].
        f1->data[i1];
    }
  }

  emxFree_uint32_T(&n);

  /* 'mapTF2_ls:401' NAME2 = NAME2_new; */
  /* 'mapTF2_ls:402' LEN = zeros(length(p),1); */
  /* 'mapTF2_ls:403' for i = 1:length(p) */
  /* 'mapTF2_ls:406' C = cumsum([0;LEN]); */
  /* 'mapTF2_ls:407' mat = [zeros(length(omat),2) ones(length(omat),4)/4]; */
  /* 'mapTF2_ls:408' for i = 1:L-length(n) */
}

/*
 * function PWM_corr2(fn,NUM,kmer, dsvm, NAME, Lmat, s)
 */
static void PWM_corr2(const emxArray_char_T *fn, double NUM, emxArray_real_T
                      *kmer, emxArray_real_T *dsvm, const emxArray_cell_wrap_2
                      *NAME, const emxArray_real_T *Lmat, const emxArray_char_T *
                      s)
{
  static const char b_cv[18] = { '_', 'k', 'm', 'e', 'r', '_', 'P', 'W', 'M',
    '_', 'l', 'o', 'c', 's', '.', 'o', 'u', 't' };

  FILE* b_NULL;
  FILE* filestar;
  const cell_wrap_2 *NAME_data;
  emxArray_boolean_T *b_coords;
  emxArray_char_T *b_fn;
  emxArray_char_T *varargin_8;
  emxArray_int32_T *ib;
  emxArray_int32_T *iidx;
  emxArray_real_T *b_dsvm;
  emxArray_real_T *b_kmer;
  emxArray_real_T *coords;
  emxArray_real_T *f;
  emxArray_real_T *y;
  const double *Lmat_data;
  double b_Lmat;
  double b_x;
  double d;
  double d1;
  double x;
  double *b_kmer_data;
  double *coords_data;
  double *dsvm_data;
  double *f_data;
  double *kmer_data;
  double *y_data;
  int b_i;
  int b_loop_ub;
  int c_loop_ub;
  int dsvm_tmp;
  int fid1;
  int i;
  int i1;
  int i2;
  int loop_ub;
  unsigned int unnamed_idx_1;
  int *ib_data;
  int *iidx_data;
  const char *fn_data;
  const char *s_data;
  signed char fileid;
  char *b_fn_data;
  char *varargin_8_data;
  bool autoflush;
  bool *b_coords_data;
  s_data = s->data;
  Lmat_data = Lmat->data;
  NAME_data = NAME->data;
  dsvm_data = dsvm->data;
  kmer_data = kmer->data;
  fn_data = fn->data;
  emxInit_real_T(&coords, 1);

  /* 'mapTF2_ls:511' [coords, ind] = sort(kmer(:,1)); */
  loop_ub = kmer->size[0];
  i = coords->size[0];
  coords->size[0] = kmer->size[0];
  emxEnsureCapacity_real_T(coords, i);
  coords_data = coords->data;
  for (i = 0; i < loop_ub; i++) {
    coords_data[i] = kmer_data[i];
  }

  emxInit_real_T(&f, 1);
  emxInit_int32_T(&iidx, 1);
  sort(coords, iidx);
  iidx_data = iidx->data;
  coords_data = coords->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_real_T(f, i);
  f_data = f->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = iidx_data[i];
  }

  emxInit_real_T(&b_kmer, 2);

  /* 'mapTF2_ls:512' kmer = kmer(ind, :); */
  i = b_kmer->size[0] * b_kmer->size[1];
  b_kmer->size[0] = f->size[0];
  b_kmer->size[1] = 3;
  emxEnsureCapacity_real_T(b_kmer, i);
  b_kmer_data = b_kmer->data;
  loop_ub = f->size[0];
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_kmer_data[i1 + b_kmer->size[0] * i] = kmer_data[((int)f_data[i1] +
        kmer->size[0] * i) - 1];
    }
  }

  i = kmer->size[0] * kmer->size[1];
  kmer->size[0] = b_kmer->size[0];
  kmer->size[1] = 3;
  emxEnsureCapacity_real_T(kmer, i);
  kmer_data = kmer->data;
  loop_ub = b_kmer->size[0] * 3;
  for (i = 0; i < loop_ub; i++) {
    kmer_data[i] = b_kmer_data[i];
  }

  emxFree_real_T(&b_kmer);
  emxInit_real_T(&b_dsvm, 1);

  /* 'mapTF2_ls:513' dsvm = dsvm(ind, :); */
  i = b_dsvm->size[0];
  b_dsvm->size[0] = f->size[0];
  emxEnsureCapacity_real_T(b_dsvm, i);
  b_kmer_data = b_dsvm->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_kmer_data[i] = dsvm_data[(int)f_data[i] - 1];
  }

  i = dsvm->size[0];
  dsvm->size[0] = b_dsvm->size[0];
  emxEnsureCapacity_real_T(dsvm, i);
  dsvm_data = dsvm->data;
  loop_ub = b_dsvm->size[0];
  for (i = 0; i < loop_ub; i++) {
    dsvm_data[i] = b_kmer_data[i];
  }

  /* 'mapTF2_ls:514' L = length(NAME); */
  /* 'mapTF2_ls:515' if NUM ==1 */
  emxInit_char_T(&b_fn, 2);
  if (NUM == 1.0) {
    /* 'mapTF2_ls:516' fid1 = fopen([fn '_kmer_PWM_locs.out'],'w'); */
    i = b_fn->size[0] * b_fn->size[1];
    b_fn->size[0] = 1;
    b_fn->size[1] = fn->size[1] + 18;
    emxEnsureCapacity_char_T(b_fn, i);
    b_fn_data = b_fn->data;
    loop_ub = fn->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_fn_data[i] = fn_data[i];
    }

    for (i = 0; i < 18; i++) {
      b_fn_data[i + fn->size[1]] = b_cv[i];
    }

    fileid = cfopen(b_fn, "wb");
    fid1 = fileid;
  } else {
    /* 'mapTF2_ls:517' else */
    /* 'mapTF2_ls:518' fid1 = fopen([fn '_kmer_PWM_locs.out'],'a+'); */
    i = b_fn->size[0] * b_fn->size[1];
    b_fn->size[0] = 1;
    b_fn->size[1] = fn->size[1] + 18;
    emxEnsureCapacity_char_T(b_fn, i);
    b_fn_data = b_fn->data;
    loop_ub = fn->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_fn_data[i] = fn_data[i];
    }

    for (i = 0; i < 18; i++) {
      b_fn_data[i + fn->size[1]] = b_cv[i];
    }

    fileid = cfopen(b_fn, "ab+");
    fid1 = fileid;
  }

  /* 'mapTF2_ls:520' ind = ones(L,1); */
  /* 'mapTF2_ls:521' for i = 1:L */
  i = NAME->size[0];
  if (0 <= NAME->size[0] - 1) {
    b_loop_ub = coords->size[0];
    c_loop_ub = coords->size[0];
    b_NULL = NULL;
  }

  emxInit_int32_T(&ib, 1);
  emxInit_char_T(&varargin_8, 2);
  emxInit_real_T(&y, 1);
  emxInit_boolean_T(&b_coords, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:522' F = find(coords > Lmat(i,2)); */
    /* 'mapTF2_ls:523' FF = find(coords <= Lmat(i,3)); */
    /* 'mapTF2_ls:524' f = intersect(F,FF); */
    b_Lmat = Lmat_data[b_i + Lmat->size[0]];
    i1 = b_coords->size[0];
    b_coords->size[0] = coords->size[0];
    emxEnsureCapacity_boolean_T(b_coords, i1);
    b_coords_data = b_coords->data;
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_coords_data[i1] = (coords_data[i1] > b_Lmat);
    }

    eml_find(b_coords, iidx);
    iidx_data = iidx->data;
    b_Lmat = Lmat_data[b_i + Lmat->size[0] * 2];
    i1 = b_coords->size[0];
    b_coords->size[0] = coords->size[0];
    emxEnsureCapacity_boolean_T(b_coords, i1);
    b_coords_data = b_coords->data;
    for (i1 = 0; i1 < c_loop_ub; i1++) {
      b_coords_data[i1] = (coords_data[i1] <= b_Lmat);
    }

    eml_find(b_coords, ib);
    ib_data = ib->data;
    i1 = b_dsvm->size[0];
    b_dsvm->size[0] = iidx->size[0];
    emxEnsureCapacity_real_T(b_dsvm, i1);
    b_kmer_data = b_dsvm->data;
    loop_ub = iidx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_kmer_data[i1] = iidx_data[i1];
    }

    i1 = y->size[0];
    y->size[0] = ib->size[0];
    emxEnsureCapacity_real_T(y, i1);
    y_data = y->data;
    loop_ub = ib->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      y_data[i1] = ib_data[i1];
    }

    do_vectors(b_dsvm, y, f, iidx, ib);
    f_data = f->data;

    /* 'mapTF2_ls:525' a1 = ip(kmer(f,end),dsvm(f,end)); */
    loop_ub = f->size[0];
    i1 = b_dsvm->size[0];
    b_dsvm->size[0] = f->size[0];
    emxEnsureCapacity_real_T(b_dsvm, i1);
    b_kmer_data = b_dsvm->data;
    i1 = y->size[0];
    y->size[0] = f->size[0];
    emxEnsureCapacity_real_T(y, i1);
    y_data = y->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      dsvm_tmp = (int)f_data[i1] - 1;
      b_kmer_data[i1] = kmer_data[dsvm_tmp + kmer->size[0] * 2];
      y_data[i1] = dsvm_data[dsvm_tmp];
    }

    /* 'mapTF2_ls:426' c = x'*y/sqrt(x'*x)/sqrt(y'*y); */
    if (f->size[0] < 1) {
      b_Lmat = 0.0;
      x = 0.0;
      b_x = 0.0;
    } else {
      b_Lmat = cblas_ddot((blasint)f->size[0], &b_kmer_data[0], (blasint)1,
                          &y_data[0], (blasint)1);
      x = cblas_ddot((blasint)f->size[0], &b_kmer_data[0], (blasint)1,
                     &b_kmer_data[0], (blasint)1);
      b_x = cblas_ddot((blasint)f->size[0], &y_data[0], (blasint)1, &y_data[0],
                       (blasint)1);
    }

    x = sqrt(x);
    b_x = sqrt(b_x);

    /* 'mapTF2_ls:526' fprintf(fid1, '%d\t%s\t%d\t%d\t%d\t%f\t%f\t%s\n', int32(NUM), NAME{i}, ... */
    /* 'mapTF2_ls:527'         int32(Lmat(i,1)), int32(Lmat(i,2)), int32(Lmat(i,3)), Lmat(i,4), a1, s(Lmat(i,2):Lmat(i,3))); */
    d = Lmat_data[b_i + Lmat->size[0]];
    d1 = Lmat_data[b_i + Lmat->size[0] * 2];
    if (d > d1) {
      i1 = -1;
      dsvm_tmp = 0;
    } else {
      i1 = (int)d - 2;
      dsvm_tmp = (int)d1;
    }

    i2 = b_fn->size[0] * b_fn->size[1];
    b_fn->size[0] = 1;
    b_fn->size[1] = NAME_data[b_i].f1->size[1] + 1;
    emxEnsureCapacity_char_T(b_fn, i2);
    b_fn_data = b_fn->data;
    loop_ub = NAME_data[b_i].f1->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      b_fn_data[i2] = NAME_data[b_i].f1->data[i2];
    }

    b_fn_data[NAME_data[b_i].f1->size[1]] = '\x00';
    unnamed_idx_1 = (unsigned int)((dsvm_tmp - i1) - 1);
    dsvm_tmp = varargin_8->size[0] * varargin_8->size[1];
    varargin_8->size[0] = 1;
    varargin_8->size[1] = (int)unnamed_idx_1 + 1;
    emxEnsureCapacity_char_T(varargin_8, dsvm_tmp);
    varargin_8_data = varargin_8->data;
    loop_ub = (int)unnamed_idx_1;
    for (dsvm_tmp = 0; dsvm_tmp < loop_ub; dsvm_tmp++) {
      varargin_8_data[dsvm_tmp] = s_data[(i1 + dsvm_tmp) + 1];
    }

    varargin_8_data[(int)unnamed_idx_1] = '\x00';
    getfilestar(fid1, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "%d\t%s\t%d\t%d\t%d\t%f\t%f\t%s\n", (int)rt_roundd_snf
              (NUM), &b_fn_data[0], (int)rt_roundd_snf(Lmat_data[b_i]), (int)
              rt_roundd_snf(Lmat_data[b_i + Lmat->size[0]]), (int)rt_roundd_snf
              (Lmat_data[b_i + Lmat->size[0] * 2]), Lmat_data[b_i + Lmat->size[0]
              * 3], b_Lmat / x / b_x, &varargin_8_data[0]);
      if (autoflush) {
        fflush(filestar);
      }
    }
  }

  emxFree_boolean_T(&b_coords);
  emxFree_char_T(&b_fn);
  emxFree_real_T(&b_dsvm);
  emxFree_real_T(&y);
  emxFree_char_T(&varargin_8);
  emxFree_int32_T(&ib);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&coords);
  emxFree_real_T(&f);

  /* 'mapTF2_ls:529' fclose(fid1); */
  cfclose(fid1);
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
 * function [mat, names, len] = getMOTIF(fn)
 */
static void getMOTIF(const emxArray_char_T *fn, emxArray_cell_wrap_3 *mat,
                     emxArray_cell_wrap_2 *names, emxArray_real_T *len)
{
  static const char b[5] = { 'M', 'O', 'T', 'I', 'F' };

  static const char b_cv[3] = { 'a', 'l', 'l' };

  FILE* b_NULL;
  FILE* filestar;
  int st;
  int wherefrom;
  long position_t;
  cell_wrap_2 *names_data;
  cell_wrap_3 *mat_data;
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

  /* 'mapTF2_ls:429' fid = fopen(fn); */
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

  /* 'mapTF2_ls:430' if fid < 0 */
  if (fid < 0) {
    /* 'mapTF2_ls:431' fprintf('ERROR: Cannot open gkmPWM motif files\n'); */
    printf("ERROR: Cannot open gkmPWM motif files\n");
    fflush(stdout);
  }

  /* 'mapTF2_ls:433' curr_pos = ftell(fid); */
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

  /* 'mapTF2_ls:434' idx=0; */
  idx = 0.0;

  /* 'mapTF2_ls:435' while ~feof(fid) */
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
      /* 'mapTF2_ls:436' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF2_ls:437' if length(line) >= 5 && strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF2_ls:438' idx=idx+1; */
          idx++;
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /* 'mapTF2_ls:441' fseek(fid, curr_pos, 'bof'); */
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

  /* 'mapTF2_ls:443' mat = cell(idx, 1); */
  ret = (int)idx;
  i = mat->size[0];
  mat->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_3(mat, i);
  mat_data = mat->data;

  /* 'mapTF2_ls:444' mat = coder.nullcopy(mat); */
  /* 'mapTF2_ls:445' names = cell(idx, 1); */
  i = names->size[0];
  names->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_2(names, i);
  names_data = names->data;
  for (i = 0; i < ret; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 4;
    names_data[i].f1->size[0] = 1;
    names_data[i].f1->size[1] = 0;
  }

  /* 'mapTF2_ls:446' names = coder.nullcopy(names); */
  /* 'mapTF2_ls:447' len = zeros(idx, 1); */
  i = len->size[0];
  len->size[0] = (int)idx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;
  for (i = 0; i < ret; i++) {
    len_data[i] = 0.0;
  }

  /* 'mapTF2_ls:449' i=0; */
  b_i = 0U;

  /* 'mapTF2_ls:450' while ~feof(fid) */
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
      /* 'mapTF2_ls:451' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF2_ls:452' if length(line) < 5 || ~strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF2_ls:455' i = i + 1; */
          b_i++;

          /* 'mapTF2_ls:456' names{i} = fgetl(fid); */
          fgetl(fid, names_data[(int)b_i - 1].f1);

          /* 'mapTF2_ls:457' len(i) = real(str2double(fgetl(fid))); */
          fgetl(fid, line);
          dc = str2double(line);
          len_data[(int)b_i - 1] = dc.re;

          /* 'mapTF2_ls:458' tmp = zeros(len(i), 4); */
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

          /* 'mapTF2_ls:459' line_count = 1; */
          line_count = 1U;

          /* 'mapTF2_ls:460' line = fgetl(fid); */
          fgetl(fid, line);

          /* 'mapTF2_ls:461' while ~isempty(line) */
          while (line->size[1] != 0) {
            /* 'mapTF2_ls:462' [v1, remain] = strtok(line); */
            b_strtok(line, v1, remain);

            /* 'mapTF2_ls:463' [v2, remain] = strtok(remain); */
            b_strtok(remain, v2, b_remain);

            /* 'mapTF2_ls:464' [v3, v4] = strtok(remain); */
            b_strtok(b_remain, line, remain);

            /* 'mapTF2_ls:465' vals = [real(str2double(v1)),  */
            /* 'mapTF2_ls:466'                 real(str2double(v2)),  */
            /* 'mapTF2_ls:467'                 real(str2double(v3)),  */
            /* 'mapTF2_ls:468'                 real(str2double(v4))]'; */
            dc = str2double(v1);
            dc1 = str2double(v2);
            dc2 = str2double(line);
            dc3 = str2double(remain);
            tmp_data[(int)line_count - 1] = dc.re;
            tmp_data[((int)line_count + tmp->size[0]) - 1] = dc1.re;
            tmp_data[((int)line_count + tmp->size[0] * 2) - 1] = dc2.re;
            tmp_data[((int)line_count + tmp->size[0] * 3) - 1] = dc3.re;

            /* 'mapTF2_ls:469' tmp(line_count,:) = vals; */
            /* 'mapTF2_ls:470' line = fgetl(fid); */
            fgetl(fid, line);

            /* 'mapTF2_ls:471' line_count = line_count + 1; */
            line_count++;
          }

          /* 'mapTF2_ls:473' mat{i} = tmp; */
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

  /*  mat = {}; */
  /*  names = {}; */
  /*  len = []; */
  /*  fid = fopen(fn); */
  /*  if fid < 0 */
  /*      fprintf('ERROR: Cannot open the combined motif file\n'); */
  /*  else */
  /*     i=0; */
  /*     while ~feof(fid) */
  /*          line = fgetl(fid); */
  /*          if length(line) >= 5 */
  /*              if strcmp(line(1:5), 'MOTIF') */
  /*                  i = i+1; */
  /*                  names{i} = fgetl(fid); */
  /*                  len = [len;str2double(fgetl(fid))]; */
  /*                  mat{i} = zeros(len(i),4); */
  /*                  for j = 1:len(i) */
  /*                      line = fgetl(fid); */
  /*                      [v1, remain] = strtok(line); */
  /*                      [v2, remain] = strtok(remain); */
  /*                      [v3, v4] = strtok(remain); */
  /*                      vals = [str2double(v1),  */
  /*                          str2double(v2),  */
  /*                          str2double(v3),  */
  /*                          str2double(v4)]; */
  /*                      mat{i}(j,:) = vals; */
  /*                  end */
  /*              end */
  /*          end */
  /*      end */
  /*  end */
  /* 'mapTF2_ls:508' fclose(fid); */
  cfclose(fid);
}

/*
 * function [mat,w] = getdenovomotif(filename)
 */
static void getdenovomotif(const emxArray_char_T *filename, emxArray_cell_wrap_4
  *mat, emxArray_real_T *w)
{
  static const char cv1[6] = { 'w', 'e', 'i', 'g', 'h', 't' };

  static const char b[5] = { 'M', 'O', 'T', 'I', 'F' };

  static const char b_cv[3] = { 'a', 'l', 'l' };

  FILE* b_NULL;
  FILE* filestar;
  long position_t;
  cell_wrap_4 *mat_data;
  emxArray_char_T *a__10;
  emxArray_char_T *a__13;
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
  /* 'mapTF2_ls:964' fid = fopen(filename); */
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

  /* 'mapTF2_ls:965' if fid < 0 */
  if (fid < 0) {
    /* 'mapTF2_ls:966' fprintf('ERROR: Cannot open gkmPWM motif files\n'); */
    printf("ERROR: Cannot open gkmPWM motif files\n");
    fflush(stdout);
  }

  /* 'mapTF2_ls:968' curr_pos = ftell(fid); */
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

  /* 'mapTF2_ls:969' idx=0; */
  idx = 0.0;

  /* 'mapTF2_ls:970' while ~feof(fid) */
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
      /* 'mapTF2_ls:971' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF2_ls:972' if length(line) >= 5 && strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF2_ls:973' idx=idx+1; */
          idx++;
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  /* 'mapTF2_ls:976' fseek(fid, curr_pos, 'bof'); */
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

  /* 'mapTF2_ls:978' mat = cell(idx, 1); */
  ret = (int)idx;
  i = mat->size[0];
  mat->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_4(mat, i);
  mat_data = mat->data;
  for (i = 0; i < ret; i++) {
    mat_data[i].f1->size[0] = 0;
    mat_data[i].f1->size[1] = 0;
  }

  /* 'mapTF2_ls:979' mat = coder.nullcopy(mat); */
  /* 'mapTF2_ls:980' w = zeros(idx, 1); */
  i = w->size[0];
  w->size[0] = (int)idx;
  emxEnsureCapacity_real_T(w, i);
  w_data = w->data;
  for (i = 0; i < ret; i++) {
    w_data[i] = 0.0;
  }

  /* 'mapTF2_ls:982' i = 0; */
  b_i = 0U;

  /* 'mapTF2_ls:983' while ~feof(fid) */
  b_NULL = NULL;
  emxInit_char_T(&a__10, 2);
  emxInit_char_T(&a__13, 2);
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
      /* 'mapTF2_ls:984' line = fgetl(fid); */
      fgetl(fid, line);
      line_data = line->data;

      /* 'mapTF2_ls:985' if length(line) < 5 || ~strcmp(line(1:5), 'MOTIF') */
      if (line->size[1] >= 5) {
        for (i = 0; i < 5; i++) {
          a[i] = line_data[i];
        }

        ret = memcmp(&a[0], &b[0], 5);
        if (ret == 0) {
          /* 'mapTF2_ls:988' i = i + 1; */
          b_i++;

          /* 'mapTF2_ls:989' line = fgetl(fid); */
          fgetl(fid, line);
          line_data = line->data;

          /* 'mapTF2_ls:991' [~, remain] = strtok(line); */
          b_strtok(line, a__10, remain);

          /* 'mapTF2_ls:992' [curr_w, remain] = strtok(remain); */
          b_strtok(remain, curr_w, b_remain);

          /* 'mapTF2_ls:993' [~, remain] = strtok(remain); */
          b_strtok(b_remain, a__10, remain);

          /* 'mapTF2_ls:994' [curr_alphabet, remain] = strtok(remain); */
          b_strtok(remain, curr_alphabet, b_remain);

          /* 'mapTF2_ls:995' [~, remain] = strtok(remain); */
          b_strtok(b_remain, a__10, remain);

          /* 'mapTF2_ls:996' [curr_length, ~] = strtok(remain); */
          b_strtok(remain, a__10, a__13);

          /* 'mapTF2_ls:997' w(i) = real(str2double(curr_w)); */
          dc = str2double(curr_w);
          w_data[(int)b_i - 1] = dc.re;

          /* 'mapTF2_ls:998' curr_alphabet = real(str2double(curr_alphabet)); */
          dc = str2double(curr_alphabet);

          /* 'mapTF2_ls:999' curr_length = real(str2double(curr_length)); */
          dc1 = str2double(a__10);

          /* 'mapTF2_ls:1000' tmp = zeros(curr_length, curr_alphabet); */
          i = tmp->size[0] * tmp->size[1];
          tmp->size[0] = (int)dc1.re;
          tmp->size[1] = (int)dc.re;
          emxEnsureCapacity_real_T(tmp, i);
          tmp_data = tmp->data;
          ret = (int)dc1.re * (int)dc.re;
          for (i = 0; i < ret; i++) {
            tmp_data[i] = 0.0;
          }

          /* 'mapTF2_ls:1001' line_count = 1; */
          line_count = 1U;

          /* 'mapTF2_ls:1002' while ~isempty(line) */
          while (line->size[1] != 0) {
            /* 'mapTF2_ls:1003' if strfind(line, 'weight') */
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
              /* 'mapTF2_ls:1004' line = fgetl(fid); */
              fgetl(fid, line);
            }

            /* 'mapTF2_ls:1006' [v1, remain] = strtok(line); */
            b_strtok(line, curr_w, remain);

            /* 'mapTF2_ls:1007' [v2, remain] = strtok(remain); */
            b_strtok(remain, curr_alphabet, b_remain);

            /* 'mapTF2_ls:1008' [v3, v4] = strtok(remain); */
            b_strtok(b_remain, a__10, a__13);

            /* 'mapTF2_ls:1009' vals = [real(str2double(v1)),  */
            /* 'mapTF2_ls:1010'                 real(str2double(v2)),  */
            /* 'mapTF2_ls:1011'                 real(str2double(v3)),  */
            /* 'mapTF2_ls:1012'                 real(str2double(v4))]'; */
            dc = str2double(curr_w);
            dc1 = str2double(curr_alphabet);
            dc2 = str2double(a__10);
            dc3 = str2double(a__13);
            vals[0] = dc.re;
            vals[1] = dc1.re;
            vals[2] = dc2.re;
            vals[3] = dc3.re;

            /* 'mapTF2_ls:1013' tmp(line_count,:) = vals; */
            ret = tmp->size[1];
            for (i = 0; i < ret; i++) {
              tmp_data[((int)line_count + tmp->size[0] * i) - 1] = vals[i];
            }

            /* 'mapTF2_ls:1014' line = fgetl(fid); */
            fgetl(fid, line);
            line_data = line->data;

            /* 'mapTF2_ls:1015' line_count = line_count + 1; */
            line_count++;
          }

          /* 'mapTF2_ls:1017' mat{i} = tmp; */
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
  emxFree_char_T(&a__13);
  emxFree_char_T(&a__10);
  emxFree_char_T(&line);

  /*   */
  /*  else */
  /*     i=0; */
  /*     while ~feof(fid) */
  /*          line = fgetl(fid); */
  /*          if length(line) >= 5 */
  /*              if strcmp(line(1:5), 'MOTIF') */
  /*                  i = i+1; */
  /*                  line = fgetl(fid); */
  /*                  mat{i} = []; */
  /*                  [~, tmp] = strtok(line); */
  /*                  [tmp] = strtok(tmp); */
  /*                  w = [w;str2double(tmp)]; */
  /*                  while ~isempty(line) */
  /*                      if strfind(line, 'weight') */
  /*                          line = fgetl(fid); */
  /*                      end */
  /*                      [v1, remain] = strtok(line); */
  /*                      [v2, remain] = strtok(remain); */
  /*                      [v3, v4] = strtok(remain); */
  /*                      vals = [str2double(v1),  */
  /*                          str2double(v2),  */
  /*                          str2double(v3),  */
  /*                          str2double(v4)]'; */
  /*                      mat{i} = [mat{i}; vals]; */
  /*                      line = fgetl(fid); */
  /*                  end */
  /*              end */
  /*          end */
  /*      end */
  /*  end */
  /*  mat = mat'; */
  /* 'mapTF2_ls:1051' fclose(fid); */
  cfclose(fid);
}

static void j_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *mat,
  const emxArray_real_T *x)
{
  emxArray_real_T *b_mat;
  const double *mat_data;
  const double *x_data;
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
  x_data = x->data;
  mat_data = mat->data;
  emxInit_real_T(&b_mat, 2);
  i = b_mat->size[0] * b_mat->size[1];
  if (x->size[0] == 1) {
    b_mat->size[0] = mat->size[0];
  } else {
    b_mat->size[0] = x->size[0];
  }

  if (x->size[1] == 1) {
    b_mat->size[1] = mat->size[1];
  } else {
    b_mat->size[1] = x->size[1];
  }

  emxEnsureCapacity_real_T(b_mat, i);
  b_mat_data = b_mat->data;
  stride_0_0 = (mat->size[0] != 1);
  stride_0_1 = (mat->size[1] != 1);
  stride_1_0 = (x->size[0] != 1);
  stride_1_1 = (x->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (x->size[1] == 1) {
    loop_ub = mat->size[1];
  } else {
    loop_ub = x->size[1];
  }

  for (i = 0; i < loop_ub; i++) {
    if (x->size[0] == 1) {
      b_loop_ub = mat->size[0];
    } else {
      b_loop_ub = x->size[0];
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_mat_data[i1 + b_mat->size[0] * i] = mat_data[i1 * stride_0_0 + mat->
        size[0] * aux_0_1] * x_data[i1 * stride_1_0 + x->size[0] * aux_1_1] /
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

  /* 'mapTF2_ls:360' l = length(s); */
  /* 'mapTF2_ls:361' en = zeros(1,l); */
  i = en->size[0] * en->size[1];
  en->size[0] = 1;
  en->size[1] = s->size[1];
  emxEnsureCapacity_real_T(en, i);
  en_data = en->data;

  /* 'mapTF2_ls:362' for i = 1:l */
  i = s->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:363' if strcmp(s(i),'A') || strcmp(s(i), 'a') */
    c = s_data[b_i];
    if ((!(c != 'A')) || (!(c != 'a'))) {
      /* 'mapTF2_ls:364' en(i) = 0; */
      en_data[b_i] = 0.0;
    } else if ((!(c != 'C')) || (!(c != 'c'))) {
      /* 'mapTF2_ls:365' elseif strcmp(s(i),'C') || strcmp(s(i),'c') */
      /* 'mapTF2_ls:366' en(i) = 1; */
      en_data[b_i] = 1.0;
    } else if ((!(c != 'G')) || (!(c != 'g'))) {
      /* 'mapTF2_ls:367' elseif strcmp(s(i),'G') || strcmp(s(i),'g') */
      /* 'mapTF2_ls:368' en(i) = 2; */
      en_data[b_i] = 2.0;
    } else {
      /* 'mapTF2_ls:369' else */
      /* 'mapTF2_ls:370' en(i) = 3; */
      en_data[b_i] = 3.0;
    }
  }
}

static void minus(emxArray_real_T *maxnorm, const emxArray_real_T *minnorm)
{
  emxArray_real_T *b_maxnorm;
  const double *minnorm_data;
  double *b_maxnorm_data;
  double *maxnorm_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  minnorm_data = minnorm->data;
  maxnorm_data = maxnorm->data;
  emxInit_real_T(&b_maxnorm, 1);
  i = b_maxnorm->size[0];
  if (minnorm->size[0] == 1) {
    b_maxnorm->size[0] = maxnorm->size[0];
  } else {
    b_maxnorm->size[0] = minnorm->size[0];
  }

  emxEnsureCapacity_real_T(b_maxnorm, i);
  b_maxnorm_data = b_maxnorm->data;
  stride_0_0 = (maxnorm->size[0] != 1);
  stride_1_0 = (minnorm->size[0] != 1);
  if (minnorm->size[0] == 1) {
    loop_ub = maxnorm->size[0];
  } else {
    loop_ub = minnorm->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_maxnorm_data[i] = maxnorm_data[i * stride_0_0] - minnorm_data[i *
      stride_1_0];
  }

  i = maxnorm->size[0];
  maxnorm->size[0] = b_maxnorm->size[0];
  emxEnsureCapacity_real_T(maxnorm, i);
  maxnorm_data = maxnorm->data;
  loop_ub = b_maxnorm->size[0];
  for (i = 0; i < loop_ub; i++) {
    maxnorm_data[i] = b_maxnorm_data[i];
  }

  emxFree_real_T(&b_maxnorm);
}

/*
 * function [ind, M] = ppmsim(mot,lenvec)
 */
static void ppmsim(emxArray_cell_wrap_4 *mot, const emxArray_real_T *lenvec,
                   double *ind, double *M)
{
  cell_wrap_4 *mot_data;
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

  /* 'mapTF2_ls:936' n = length(lenvec)-1; */
  /* 'mapTF2_ls:937' simmat = ones(n-1,1); */
  /* 'mapTF2_ls:938' for i = 1:n+1 */
  i = lenvec->size[0];
  emxInit_real_T(&diag_a, 2);
  emxInit_real_T(&A, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:939' mot{i} = mot{i}-1/4; */
    loop_ub = mot_data[b_i].f1->size[0] * mot_data[b_i].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      mot_data[b_i].f1->data[i1] -= 0.25;
    }

    /* 'mapTF2_ls:940' mot{i} = mot{i}/sqrt(sum(sum(mot{i}.^2))); */
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

  /* 'mapTF2_ls:942' M = 0; */
  *M = 0.0;

  /* 'mapTF2_ls:943' ind = 1; */
  *ind = 1.0;

  /* 'mapTF2_ls:944' for j = 2:n+1 */
  i = lenvec->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&rmat, 2);
  emxInit_real_T(&diag_b, 2);
  emxInit_real_T(&diag_overall, 1);
  for (b_i = 0; b_i <= i - 2; b_i++) {
    /* 'mapTF2_ls:945' mat = mot{1}*mot{j}'; */
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

    /* 'mapTF2_ls:946' rmat = rot90(mot{1},2)*mot{j}'; */
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

    /* 'mapTF2_ls:947' diag_a = sum(spdiags(mat)); */
    spdiags(mat, A);
    b_sum(A, diag_a);
    mat_data = diag_a->data;

    /* 'mapTF2_ls:948' diag_b = sum(spdiags(rmat)); */
    spdiags(rmat, A);
    b_sum(A, diag_b);
    rmat_data = diag_b->data;

    /* 'mapTF2_ls:949' diag_overall = zeros(length(diag_a)+length(diag_b),1); */
    i1 = diag_overall->size[0];
    diag_overall->size[0] = (int)((unsigned int)diag_a->size[1] + diag_b->size[1]);
    emxEnsureCapacity_real_T(diag_overall, i1);
    diag_overall_data = diag_overall->data;
    loop_ub = (int)((unsigned int)diag_a->size[1] + diag_b->size[1]);
    for (i1 = 0; i1 < loop_ub; i1++) {
      diag_overall_data[i1] = 0.0;
    }

    /* 'mapTF2_ls:950' diag_overall(1:length(diag_a)) = diag_a; */
    loop_ub = diag_a->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      diag_overall_data[i1] = mat_data[i1];
    }

    /* 'mapTF2_ls:951' diag_overall(length(diag_a)+1:end) = diag_b; */
    if (diag_a->size[1] + 1U > (unsigned int)diag_overall->size[0]) {
      i1 = 0;
    } else {
      i1 = diag_a->size[1];
    }

    loop_ub = diag_b->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      diag_overall_data[i1 + i2] = rmat_data[i2];
    }

    /* 'mapTF2_ls:952' MM = max(diag_overall); */
    MM = maximum(diag_overall);

    /*  MM = max([sum(spdiags(mat)) sum(spdiags(rmat))]); */
    /* 'mapTF2_ls:954' if MM > M */
    if (MM > *M) {
      /* 'mapTF2_ls:955' M = MM; */
      *M = MM;

      /* 'mapTF2_ls:956' ind = j-1; */
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
  cell_wrap_2 *names_data;
  cell_wrap_4 r1;
  cell_wrap_4 *PWM2_data;
  cell_wrap_4 *P_data;
  cell_wrap_4 *cur_PWM_data;
  cell_wrap_4 *cur_PWM_tmp_data;
  emxArray_boolean_T *b_clus;
  emxArray_cell_wrap_2 *names;
  emxArray_cell_wrap_4 *P;
  emxArray_cell_wrap_4 *PWM2;
  emxArray_cell_wrap_4 *b_cur_PWM_tmp;
  emxArray_cell_wrap_4 *b_new_PWM2;
  emxArray_cell_wrap_4 *c_cur_PWM_tmp;
  emxArray_cell_wrap_4 *cur_PWM;
  emxArray_cell_wrap_4 *cur_PWM_tmp;
  emxArray_cell_wrap_4 *d_cur_PWM_tmp;
  emxArray_cell_wrap_4 *new_PWM2;
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
  int b_loop_ub;
  int cur_idx;
  int exitg1;
  int i;
  int i1;
  int i2;
  int j;
  int loop_ub;
  int match_idx;
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
  emxInit_cell_wrap_4(&cur_PWM_tmp, 1);

  /* dfn: file name for denovo motifs */
  /* lfn: file name for lasso motifs */
  /* memefn: file name for the meme input for gkmPWMlasso */
  /* ofn: output filename */
  /* 'mapTF2_ls:689' a = 1; */
  a = 1.0;

  /* 'mapTF2_ls:690' LEN = zeros(1,1); */
  /* 'mapTF2_ls:691' shift = zeros(1,1); */
  /* 'mapTF2_ls:692' [p,w] = getdenovomotif(dfn); */
  getdenovomotif(dfn, cur_PWM_tmp, f);
  f_data = f->data;
  cur_PWM_tmp_data = cur_PWM_tmp->data;

  /* 'mapTF2_ls:693' N = numel(w); */
  N = f->size[0];

  /*  fid = fopen(lfn,'r'); */
  /*  X = textscan(fid,'%f\t%f\t%s\t%f\t%f\t%f\n','delimiter', '\t', 'headerlines', 4); */
  /*  fclose(fid); */
  /* 'mapTF2_ls:698' fid = fopen(lfn, 'r'); */
  fileid = cfopen(lfn, "rb");

  /* 'mapTF2_ls:699' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF2_ls:700' fprintf("ERROR: gkmPWMlasso output file cannot be opened.\n") */
    printf("ERROR: gkmPWMlasso output file cannot be opened.\n");
    fflush(stdout);
  }

  emxInit_char_T(&b_fileid, 2);
  emxInit_char_T(&c_fileid, 2);
  emxInit_char_T(&d_fileid, 2);
  emxInit_char_T(&e_fileid, 2);

  /* 'mapTF2_ls:703' fgetl(fid); */
  b_fgets(fileid, b_fileid);

  /* 'mapTF2_ls:704' fgetl(fid); */
  b_fgets(fileid, c_fileid);

  /* 'mapTF2_ls:705' fgetl(fid); */
  b_fgets(fileid, d_fileid);

  /* 'mapTF2_ls:706' fgetl(fid); */
  b_fgets(fileid, e_fileid);

  /* 'mapTF2_ls:707' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF2_ls:708' idx=0; */
  idx = 0.0;

  /* 'mapTF2_ls:709' while ~feof(fid) */
  emxFree_char_T(&e_fileid);
  emxFree_char_T(&d_fileid);
  emxFree_char_T(&c_fileid);
  emxFree_char_T(&b_fileid);
  emxInit_char_T(&f_fileid, 2);
  do {
    exitg1 = 0;
    d = b_feof(fileid);
    if (!(d != 0.0)) {
      /* 'mapTF2_ls:710' idx=idx+1; */
      idx++;

      /* 'mapTF2_ls:711' fgetl(fid); */
      b_fgets(fileid, f_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&f_fileid);
  emxInit_real_T(&clus, 1);

  /* 'mapTF2_ls:713' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF2_ls:714' clus = zeros(idx, 1); */
  match_idx = (int)idx;
  i = clus->size[0];
  clus->size[0] = (int)idx;
  emxEnsureCapacity_real_T(clus, i);
  clus_data = clus->data;
  for (i = 0; i < match_idx; i++) {
    clus_data[i] = 0.0;
  }

  emxInit_real_T(&uid, 1);

  /* 'mapTF2_ls:715' uid = zeros(idx, 1); */
  i = uid->size[0];
  uid->size[0] = (int)idx;
  emxEnsureCapacity_real_T(uid, i);
  uid_data = uid->data;
  for (i = 0; i < match_idx; i++) {
    uid_data[i] = 0.0;
  }

  emxInit_real_T(&lasso_weight, 1);

  /* 'mapTF2_ls:716' lasso_weight = zeros(idx, 1); */
  i = lasso_weight->size[0];
  lasso_weight->size[0] = (int)idx;
  emxEnsureCapacity_real_T(lasso_weight, i);
  lasso_weight_data = lasso_weight->data;
  for (i = 0; i < match_idx; i++) {
    lasso_weight_data[i] = 0.0;
  }

  emxInit_real_T(&zscore, 1);

  /* 'mapTF2_ls:717' zscore = zeros(idx, 1); */
  i = zscore->size[0];
  zscore->size[0] = (int)idx;
  emxEnsureCapacity_real_T(zscore, i);
  zscore_data = zscore->data;
  for (i = 0; i < match_idx; i++) {
    zscore_data[i] = 0.0;
  }

  /* 'mapTF2_ls:719' for cur_idx=1:idx */
  cur_idx = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&p1, 2);
  emxInit_char_T(&remain, 2);
  emxInit_char_T(&p2, 2);
  emxInit_char_T(&b_remain, 2);
  emxInit_char_T(&p4, 2);
  exitg2 = false;
  while ((!exitg2) && (cur_idx <= (int)idx - 1)) {
    /* 'mapTF2_ls:720' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF2_ls:721' if cur_line == -1 */
    hal = (cur_line->size[1] != 0);
    if (hal) {
      hal = (0 > cur_line->size[1] - 1);
    }

    if (hal) {
      exitg2 = true;
    } else {
      /* 'mapTF2_ls:724' [p1, remain] = strtok(cur_line, char(9)); */
      c_strtok(cur_line, p1, remain);

      /* 'mapTF2_ls:725' [p2, remain] = strtok(remain, char(9)); */
      c_strtok(remain, p2, b_remain);

      /* 'mapTF2_ls:726' [p3, remain] = strtok(remain, char(9)); */
      c_strtok(b_remain, cur_line, remain);

      /* 'mapTF2_ls:727' [p4, remain] = strtok(remain, char(9)); */
      c_strtok(remain, p4, b_remain);

      /* 'mapTF2_ls:728' [p5, p6] = strtok(remain, char(9)); */
      c_strtok(b_remain, cur_line, remain);

      /* 'mapTF2_ls:729' clus(cur_idx, 1) = real(str2double(p1)); */
      dc = str2double(p1);
      clus_data[cur_idx] = dc.re;

      /* 'mapTF2_ls:730' uid(cur_idx, 1) = real(str2double(p2)); */
      dc = str2double(p2);
      uid_data[cur_idx] = dc.re;

      /* 'mapTF2_ls:731' lasso_weight(cur_idx, 1) = real(str2double(p4)); */
      dc = str2double(p4);
      lasso_weight_data[cur_idx] = dc.re;

      /* 'mapTF2_ls:732' zscore(cur_idx, 1) = real(str2double(p5)); */
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

  /* 'mapTF2_ls:734' fclose(fid); */
  cfclose(fileid);

  /*  n = X{1}(end); */
  /* 'mapTF2_ls:737' n = clus(end); */
  /* 'mapTF2_ls:738' vec = zeros(n,1); */
  i = vec->size[0];
  vec->size[0] = (int)clus_data[clus->size[0] - 1];
  emxEnsureCapacity_real_T(vec, i);
  vec_data = vec->data;

  /* 'mapTF2_ls:739' vec2 = zeros(n,1); */
  i = vec2->size[0];
  vec2->size[0] = (int)clus_data[clus->size[0] - 1];
  emxEnsureCapacity_real_T(vec2, i);
  vec2_data = vec2->data;

  /* 'mapTF2_ls:740' w = [w; zeros(n,1)]; */
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

  /* 'mapTF2_ls:741' for i = 1:n */
  i = (int)clus_data[clus->size[0] - 1];
  emxInit_int32_T(&r, 1);
  emxInit_boolean_T(&b_clus, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /*      f = find(X{1}==i); */
    /*      vec(i) = X{2}(f(1)); */
    /*      vec2(i) = X{5}(f(1)); */
    /*      w(i+N) = X{4}(f(1)); */
    /* 'mapTF2_ls:746' f = find(clus == i); */
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

    /* 'mapTF2_ls:747' vec(i) = uid(f(1)); */
    vec_data[b_i] = uid_data[(int)f_data[0] - 1];

    /* 'mapTF2_ls:748' vec2(i) = zscore(f(1)); */
    vec2_data[b_i] = zscore_data[(int)f_data[0] - 1];

    /* 'mapTF2_ls:749' w(i+N) = lasso_weight(f(1)); */
    w_data[(int)((unsigned int)b_i + N)] = lasso_weight_data[(int)f_data[0] - 1];
  }

  emxFree_boolean_T(&b_clus);
  emxFree_int32_T(&r);
  emxFree_real_T(&zscore);
  emxInit_cell_wrap_4(&P, 1);
  emxInit_cell_wrap_4(&cur_PWM, 1);

  /* 'mapTF2_ls:751' P = getmotif(memefn,vec); */
  getmotif(memefn, vec, P);
  P_data = P->data;

  /* 'mapTF2_ls:753' denovo_len = length(p); */
  /* 'mapTF2_ls:754' database_len = length(P); */
  /* 'mapTF2_ls:755' cur_PWM = cell(denovo_len+database_len,1); */
  unnamed_idx_0 = (int)((unsigned int)cur_PWM_tmp->size[0] + P->size[0]);
  i = cur_PWM->size[0];
  cur_PWM->size[0] = (int)((unsigned int)cur_PWM_tmp->size[0] + P->size[0]);
  emxEnsureCapacity_cell_wrap_4(cur_PWM, i);
  cur_PWM_data = cur_PWM->data;
  emxFree_real_T(&vec);
  for (i = 0; i < unnamed_idx_0; i++) {
    cur_PWM_data[i].f1->size[0] = 0;
    cur_PWM_data[i].f1->size[1] = 0;
  }

  /* 'mapTF2_ls:756' cur_PWM = coder.nullcopy(cur_PWM); */
  /* 'mapTF2_ls:757' for cur_idx=1:denovo_len */
  i = cur_PWM_tmp->size[0];
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'mapTF2_ls:758' cur_PWM{cur_idx} = p{cur_idx}; */
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

  /* 'mapTF2_ls:760' for cur_idx=denovo_len+1:denovo_len+database_len */
  i = P->size[0];
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    b_cur_idx = ((unsigned int)cur_PWM_tmp->size[0] + cur_idx) + 1U;

    /* 'mapTF2_ls:761' cur_PWM{cur_idx} = P{cur_idx-denovo_len}; */
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

  emxFree_cell_wrap_4(&P);
  emxInitStruct_cell_wrap_4(&r1);

  /*  [pp, info, len] = trim_pwm([p;P],0.25); */
  /* 'mapTF2_ls:765' [pp, info, len] = trim_pwm(cur_PWM,0.25); */
  trim_pwm(cur_PWM, lasso_weight, uid);
  uid_data = uid->data;
  lasso_weight_data = lasso_weight->data;
  cur_PWM_data = cur_PWM->data;

  /* 'mapTF2_ls:767' coder.varsize("PWM2"); */
  /* 'mapTF2_ls:768' PWM2 = {pp{1}}; */
  i = r1.f1->size[0] * r1.f1->size[1];
  r1.f1->size[0] = cur_PWM_data[0].f1->size[0];
  r1.f1->size[1] = cur_PWM_data[0].f1->size[1];
  emxEnsureCapacity_real_T(r1.f1, i);
  loop_ub = cur_PWM_data[0].f1->size[0] * cur_PWM_data[0].f1->size[1];
  for (i = 0; i < loop_ub; i++) {
    r1.f1->data[i] = cur_PWM_data[0].f1->data[i];
  }

  emxInit_cell_wrap_4(&PWM2, 2);
  emxInit_real_T(&LEN_2, 2);
  emxInit_real_T(&I_2, 2);
  i = PWM2->size[0] * PWM2->size[1];
  PWM2->size[0] = 1;
  PWM2->size[1] = 1;
  emxEnsureCapacity_cell_wrap_4(PWM2, i);
  PWM2_data = PWM2->data;
  emxCopyStruct_cell_wrap_4(&PWM2_data[0], &r1);

  /*  PWM2 = coder.nullcopy(PWM2); */
  /*  PWM2{1} = pp{1}; */
  /* 'mapTF2_ls:772' coder.varsize("LEN_2"); */
  /* 'mapTF2_ls:773' coder.varsize("I_2"); */
  /* 'mapTF2_ls:774' LEN_2 = zeros(1,1); */
  i = LEN_2->size[0] * LEN_2->size[1];
  LEN_2->size[0] = 1;
  LEN_2->size[1] = 1;
  emxEnsureCapacity_real_T(LEN_2, i);
  vec_data = LEN_2->data;
  vec_data[0] = 0.0;

  /* 'mapTF2_ls:775' I_2 = zeros(1,1); */
  i = I_2->size[0] * I_2->size[1];
  I_2->size[0] = 1;
  I_2->size[1] = 1;
  emxEnsureCapacity_real_T(I_2, i);
  I_2_data = I_2->data;
  I_2_data[0] = 0.0;

  /* 'mapTF2_ls:776' hal = true; */
  hal = true;

  /* 'mapTF2_ls:777' for ii = 1:length(w) */
  i = w->size[0];
  emxFreeStruct_cell_wrap_4(&r1);
  emxInit_cell_wrap_4(&b_cur_PWM_tmp, 1);
  emxInit_cell_wrap_4(&new_PWM2, 1);
  emxInit_cell_wrap_4(&c_cur_PWM_tmp, 1);
  emxInit_cell_wrap_4(&b_new_PWM2, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:778' if ii > N */
    if (b_i + 1 > N) {
      /* 'mapTF2_ls:779' hal = true; */
      hal = true;

      /* 'mapTF2_ls:780' PWM2_len = length(PWM2); */
      match_idx = PWM2->size[0];
      if (match_idx < 1) {
        match_idx = 1;
      }

      if (PWM2->size[0] == 0) {
        match_idx = 0;
      }

      /* 'mapTF2_ls:781' cur_PWM_tmp = cell(PWM2_len+1,1); */
      unnamed_idx_0 = match_idx + 1;
      i1 = c_cur_PWM_tmp->size[0];
      c_cur_PWM_tmp->size[0] = match_idx + 1;
      emxEnsureCapacity_cell_wrap_4(c_cur_PWM_tmp, i1);
      cur_PWM_tmp_data = c_cur_PWM_tmp->data;
      for (i1 = 0; i1 < unnamed_idx_0; i1++) {
        cur_PWM_tmp_data[i1].f1->size[0] = 0;
        cur_PWM_tmp_data[i1].f1->size[1] = 0;
      }

      /* 'mapTF2_ls:782' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
      i1 = cur_PWM_tmp->size[0];
      cur_PWM_tmp->size[0] = c_cur_PWM_tmp->size[0];
      emxEnsureCapacity_cell_wrap_4(cur_PWM_tmp, i1);
      cur_PWM_tmp_data = cur_PWM_tmp->data;

      /* 'mapTF2_ls:783' cur_PWM_tmp{1} = pp{ii}; */
      i1 = cur_PWM_tmp_data[0].f1->size[0] * cur_PWM_tmp_data[0].f1->size[1];
      cur_PWM_tmp_data[0].f1->size[0] = cur_PWM_data[b_i].f1->size[0];
      cur_PWM_tmp_data[0].f1->size[1] = cur_PWM_data[b_i].f1->size[1];
      emxEnsureCapacity_real_T(cur_PWM_tmp_data[0].f1, i1);
      loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cur_PWM_tmp_data[0].f1->data[i1] = cur_PWM_data[b_i].f1->data[i1];
      }

      /* 'mapTF2_ls:784' for cur_idx=1:PWM2_len */
      for (cur_idx = 0; cur_idx < match_idx; cur_idx++) {
        /* 'mapTF2_ls:785' cur_PWM_tmp{cur_idx+1} = PWM2{cur_idx}; */
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

      /* 'mapTF2_ls:787' cur_length = zeros(length(LEN_2)+1,1); */
      match_idx = LEN_2->size[0];
      if (match_idx < 1) {
        match_idx = 1;
      }

      if (LEN_2->size[0] == 0) {
        match_idx = 0;
      }

      i1 = f->size[0];
      f->size[0] = match_idx + 1;
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      for (i1 = 0; i1 <= match_idx; i1++) {
        f_data[i1] = 0.0;
      }

      /* 'mapTF2_ls:788' cur_length(1) = len(ii); */
      f_data[0] = uid_data[b_i];

      /* 'mapTF2_ls:789' cur_length(2:end) = LEN_2; */
      for (i1 = 0; i1 < match_idx; i1++) {
        f_data[i1 + 1] = vec_data[i1];
      }

      /* 'mapTF2_ls:790' [~,cor] = ppmsim(cur_PWM_tmp, cur_length); */
      ppmsim(cur_PWM_tmp, f, &curr_pos, &idx);

      /*  [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]); */
      /* 'mapTF2_ls:792' if cor > 0.7 || vec2(ii-N) < 1.5 */
      if ((idx > 0.7) || (vec2_data[b_i - N] < 1.5)) {
        /* 'mapTF2_ls:793' hal = false; */
        hal = false;
      }
    } else if ((b_i + 1 > 1) && (a > 2.0)) {
      /* 'mapTF2_ls:795' elseif ii > 1 && a > 2 */
      /* 'mapTF2_ls:796' hal = true; */
      hal = true;

      /* 'mapTF2_ls:797' PWM2_len = length(PWM2); */
      match_idx = PWM2->size[0];
      if (match_idx < 1) {
        match_idx = 1;
      }

      if (PWM2->size[0] == 0) {
        match_idx = 0;
      }

      /* 'mapTF2_ls:798' cur_PWM_tmp = cell(PWM2_len+1,1); */
      unnamed_idx_0 = match_idx + 1;
      i1 = b_cur_PWM_tmp->size[0];
      b_cur_PWM_tmp->size[0] = match_idx + 1;
      emxEnsureCapacity_cell_wrap_4(b_cur_PWM_tmp, i1);
      P_data = b_cur_PWM_tmp->data;
      for (i1 = 0; i1 < unnamed_idx_0; i1++) {
        P_data[i1].f1->size[0] = 0;
        P_data[i1].f1->size[1] = 0;
      }

      /* 'mapTF2_ls:799' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
      i1 = cur_PWM_tmp->size[0];
      cur_PWM_tmp->size[0] = b_cur_PWM_tmp->size[0];
      emxEnsureCapacity_cell_wrap_4(cur_PWM_tmp, i1);
      cur_PWM_tmp_data = cur_PWM_tmp->data;

      /* 'mapTF2_ls:800' cur_PWM_tmp{1} = pp{ii}; */
      i1 = cur_PWM_tmp_data[0].f1->size[0] * cur_PWM_tmp_data[0].f1->size[1];
      cur_PWM_tmp_data[0].f1->size[0] = cur_PWM_data[b_i].f1->size[0];
      cur_PWM_tmp_data[0].f1->size[1] = cur_PWM_data[b_i].f1->size[1];
      emxEnsureCapacity_real_T(cur_PWM_tmp_data[0].f1, i1);
      loop_ub = cur_PWM_data[b_i].f1->size[0] * cur_PWM_data[b_i].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        cur_PWM_tmp_data[0].f1->data[i1] = cur_PWM_data[b_i].f1->data[i1];
      }

      /* 'mapTF2_ls:801' for cur_idx=1:PWM2_len */
      for (cur_idx = 0; cur_idx < match_idx; cur_idx++) {
        /* 'mapTF2_ls:802' cur_PWM_tmp{cur_idx+1} = PWM2{cur_idx}; */
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

      /* 'mapTF2_ls:804' cur_length = zeros(length(LEN_2)+1,1); */
      match_idx = LEN_2->size[0];
      if (match_idx < 1) {
        match_idx = 1;
      }

      if (LEN_2->size[0] == 0) {
        match_idx = 0;
      }

      i1 = f->size[0];
      f->size[0] = match_idx + 1;
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      for (i1 = 0; i1 <= match_idx; i1++) {
        f_data[i1] = 0.0;
      }

      /* 'mapTF2_ls:805' cur_length(1) = len(ii); */
      f_data[0] = uid_data[b_i];

      /* 'mapTF2_ls:806' cur_length(2:end) = LEN_2; */
      for (i1 = 0; i1 < match_idx; i1++) {
        f_data[i1 + 1] = vec_data[i1];
      }

      /* 'mapTF2_ls:807' [~,cor] = ppmsim(cur_PWM_tmp, cur_length); */
      ppmsim(cur_PWM_tmp, f, &curr_pos, &idx);

      /*  [~,cor] = ppmsim([pp{ii} PWM2], [len(ii) LEN_2]); */
      /* 'mapTF2_ls:809' if cor > 0.7 */
      if (idx > 0.7) {
        /* 'mapTF2_ls:810' hal = false; */
        hal = false;
      }
    }

    /* 'mapTF2_ls:813' if hal && w(ii) > 0 && len(ii) >= 6 */
    if (hal && (w_data[b_i] > 0.0) && (uid_data[b_i] >= 6.0)) {
      /* 'mapTF2_ls:814' if len(ii) > 10 */
      if (uid_data[b_i] > 10.0) {
        /* 'mapTF2_ls:815' if info(ii)/len(ii) > 0.7 */
        if (lasso_weight_data[b_i] / uid_data[b_i] > 0.7) {
          /* 'mapTF2_ls:816' new_len = a+1; */
          /* 'mapTF2_ls:817' new_PWM2 = cell(new_len,1); */
          j = (int)(a + 1.0);
          unnamed_idx_0 = (int)(a + 1.0);
          i1 = b_new_PWM2->size[0];
          b_new_PWM2->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_cell_wrap_4(b_new_PWM2, i1);
          P_data = b_new_PWM2->data;
          for (i1 = 0; i1 < j; i1++) {
            P_data[i1].f1->size[0] = 0;
            P_data[i1].f1->size[1] = 0;
          }

          /* 'mapTF2_ls:818' new_PWM2 = coder.nullcopy(new_PWM2); */
          i1 = cur_PWM_tmp->size[0];
          cur_PWM_tmp->size[0] = b_new_PWM2->size[0];
          emxEnsureCapacity_cell_wrap_4(cur_PWM_tmp, i1);
          cur_PWM_tmp_data = cur_PWM_tmp->data;

          /* 'mapTF2_ls:819' for idx = 1:new_len */
          for (match_idx = 0; match_idx < j; match_idx++) {
            /* 'mapTF2_ls:820' if idx <= a-1 */
            if ((double)match_idx + 1.0 <= a - 1.0) {
              /* 'mapTF2_ls:821' new_PWM2{idx} = PWM2{idx}; */
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
              /* 'mapTF2_ls:822' elseif idx == a */
              /* 'mapTF2_ls:823' new_PWM2{idx} = pp{ii}; */
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
              /* 'mapTF2_ls:824' elseif idx == a+1 */
              /* 'mapTF2_ls:825' new_PWM2{a+1} = rot90(pp{ii},2); */
              rot90(cur_PWM_data[b_i].f1, cur_PWM_tmp_data[(int)(a + 1.0) - 1].
                    f1);
            }
          }

          /* 'mapTF2_ls:828' PWM2_len = length(new_PWM2); */
          /* 'mapTF2_ls:829' PWM2 = cell(PWM2_len, 1); */
          match_idx = cur_PWM_tmp->size[0];
          i1 = PWM2->size[0] * PWM2->size[1];
          PWM2->size[0] = cur_PWM_tmp->size[0];
          PWM2->size[1] = 1;
          emxEnsureCapacity_cell_wrap_4(PWM2, i1);
          PWM2_data = PWM2->data;
          for (i1 = 0; i1 < match_idx; i1++) {
            PWM2_data[i1].f1->size[0] = 0;
            PWM2_data[i1].f1->size[1] = 0;
          }

          /* 'mapTF2_ls:830' PWM2 = coder.nullcopy(PWM2); */
          /* 'mapTF2_ls:831' for cur_idx = 1:PWM2_len */
          i1 = cur_PWM_tmp->size[0];
          for (cur_idx = 0; cur_idx < i1; cur_idx++) {
            /* 'mapTF2_ls:832' PWM2{cur_idx} = new_PWM2{cur_idx}; */
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

          /*  PWM2{end+1} = pp{ii}; */
          /*  LEN_2(a) = len(ii); */
          /*  I_2(a) = info(ii); */
          /*  PWM2{end+1} = rot90(PWM2{a-1},2); */
          /*  LEN_2(a) = LEN_2(a-1); */
          /*  I_2(a) = I_2(a-1); */
          /* 'mapTF2_ls:843' LEN_2_new = zeros(a+1,1); */
          i1 = f->size[0];
          f->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_real_T(f, i1);
          f_data = f->data;
          for (i1 = 0; i1 < j; i1++) {
            f_data[i1] = 0.0;
          }

          /* 'mapTF2_ls:844' LEN_2_new(1:a-1) = LEN_2(1:a-1); */
          if (1.0 > a - 1.0) {
            loop_ub = 0;
          } else {
            loop_ub = (int)(a - 1.0);
          }

          for (i1 = 0; i1 < loop_ub; i1++) {
            f_data[i1] = vec_data[i1];
          }

          /* 'mapTF2_ls:845' LEN_2_new(a:a+1) = [len(ii), len(ii)]; */
          f_data[(int)a - 1] = uid_data[b_i];
          f_data[(int)(a + 1.0) - 1] = uid_data[b_i];

          /* 'mapTF2_ls:847' I_2_new = zeros(a+1,1); */
          i1 = clus->size[0];
          clus->size[0] = (int)(a + 1.0);
          emxEnsureCapacity_real_T(clus, i1);
          clus_data = clus->data;
          for (i1 = 0; i1 < j; i1++) {
            clus_data[i1] = 0.0;
          }

          /* 'mapTF2_ls:848' I_2_new(1:a-1) = I_2(1:a-1); */
          if (1.0 > a - 1.0) {
            loop_ub = 0;
          } else {
            loop_ub = (int)(a - 1.0);
          }

          for (i1 = 0; i1 < loop_ub; i1++) {
            clus_data[i1] = I_2_data[i1];
          }

          /* 'mapTF2_ls:849' I_2_new(a:a+1) = [info(ii), info(ii)]; */
          clus_data[(int)a - 1] = lasso_weight_data[b_i];
          clus_data[(int)(a + 1.0) - 1] = lasso_weight_data[b_i];

          /* 'mapTF2_ls:851' LEN_2 = LEN_2_new; */
          i1 = LEN_2->size[0] * LEN_2->size[1];
          LEN_2->size[0] = f->size[0];
          LEN_2->size[1] = 1;
          emxEnsureCapacity_real_T(LEN_2, i1);
          vec_data = LEN_2->data;
          loop_ub = f->size[0];
          for (i1 = 0; i1 < loop_ub; i1++) {
            vec_data[i1] = f_data[i1];
          }

          /* 'mapTF2_ls:852' I_2 = I_2_new; */
          i1 = I_2->size[0] * I_2->size[1];
          I_2->size[0] = clus->size[0];
          I_2->size[1] = 1;
          emxEnsureCapacity_real_T(I_2, i1);
          I_2_data = I_2->data;
          loop_ub = clus->size[0];
          for (i1 = 0; i1 < loop_ub; i1++) {
            I_2_data[i1] = clus_data[i1];
          }

          /* 'mapTF2_ls:853' a = a+1; */
          a++;

          /* 'mapTF2_ls:854' a = a+1; */
          a++;
        }
      } else if ((lasso_weight_data[b_i] > 6.0) || (lasso_weight_data[b_i] /
                  uid_data[b_i] > 1.0)) {
        /* 'mapTF2_ls:856' elseif info(ii) > 6 || info(ii)/len(ii) > 1 */
        /* 'mapTF2_ls:857' new_len = a+1; */
        /* 'mapTF2_ls:858' new_PWM2 = cell(new_len,1); */
        j = (int)(a + 1.0);
        unnamed_idx_0 = (int)(a + 1.0);
        i1 = new_PWM2->size[0];
        new_PWM2->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_cell_wrap_4(new_PWM2, i1);
        P_data = new_PWM2->data;
        for (i1 = 0; i1 < j; i1++) {
          P_data[i1].f1->size[0] = 0;
          P_data[i1].f1->size[1] = 0;
        }

        /* 'mapTF2_ls:859' new_PWM2 = coder.nullcopy(new_PWM2); */
        i1 = cur_PWM_tmp->size[0];
        cur_PWM_tmp->size[0] = new_PWM2->size[0];
        emxEnsureCapacity_cell_wrap_4(cur_PWM_tmp, i1);
        cur_PWM_tmp_data = cur_PWM_tmp->data;

        /* 'mapTF2_ls:860' for idx = 1:new_len */
        for (match_idx = 0; match_idx < j; match_idx++) {
          /* 'mapTF2_ls:861' if idx <= a-1 */
          if ((double)match_idx + 1.0 <= a - 1.0) {
            /* 'mapTF2_ls:862' new_PWM2{idx} = PWM2{idx}; */
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
            /* 'mapTF2_ls:863' elseif idx == a */
            /* 'mapTF2_ls:864' new_PWM2{idx} = pp{ii}; */
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
            /* 'mapTF2_ls:865' elseif idx == a+1 */
            /* 'mapTF2_ls:866' new_PWM2{a+1} = rot90(pp{ii},2); */
            rot90(cur_PWM_data[b_i].f1, cur_PWM_tmp_data[(int)(a + 1.0) - 1].f1);
          }
        }

        /* 'mapTF2_ls:869' PWM2_len = length(new_PWM2); */
        /* 'mapTF2_ls:870' PWM2 = cell(PWM2_len, 1); */
        match_idx = cur_PWM_tmp->size[0];
        i1 = PWM2->size[0] * PWM2->size[1];
        PWM2->size[0] = cur_PWM_tmp->size[0];
        PWM2->size[1] = 1;
        emxEnsureCapacity_cell_wrap_4(PWM2, i1);
        PWM2_data = PWM2->data;
        for (i1 = 0; i1 < match_idx; i1++) {
          PWM2_data[i1].f1->size[0] = 0;
          PWM2_data[i1].f1->size[1] = 0;
        }

        /* 'mapTF2_ls:871' PWM2 = coder.nullcopy(PWM2); */
        /* 'mapTF2_ls:872' for cur_idx = 1:PWM2_len */
        i1 = cur_PWM_tmp->size[0];
        for (cur_idx = 0; cur_idx < i1; cur_idx++) {
          /* 'mapTF2_ls:873' PWM2{cur_idx} = new_PWM2{cur_idx}; */
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

        /*  PWM2{end+1} = pp{ii}; */
        /*  LEN_2(a) = len(ii); */
        /*  I_2(a) = info(ii); */
        /*  PWM2{end+1} = rot90(PWM2{a-1},2); */
        /*  LEN_2(a) = LEN_2(a-1); */
        /*  I_2(a) = I_2(a-1); */
        /* 'mapTF2_ls:884' LEN_2_new = zeros(a+1,1); */
        i1 = f->size[0];
        f->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        for (i1 = 0; i1 < j; i1++) {
          f_data[i1] = 0.0;
        }

        /* 'mapTF2_ls:885' LEN_2_new(1:a-1) = LEN_2(1:a-1); */
        if (1.0 > a - 1.0) {
          loop_ub = 0;
        } else {
          loop_ub = (int)(a - 1.0);
        }

        for (i1 = 0; i1 < loop_ub; i1++) {
          f_data[i1] = vec_data[i1];
        }

        /* 'mapTF2_ls:886' LEN_2_new(a:a+1) = [len(ii), len(ii)]; */
        f_data[(int)a - 1] = uid_data[b_i];
        f_data[(int)(a + 1.0) - 1] = uid_data[b_i];

        /* 'mapTF2_ls:888' I_2_new = zeros(a+1,1); */
        i1 = clus->size[0];
        clus->size[0] = (int)(a + 1.0);
        emxEnsureCapacity_real_T(clus, i1);
        clus_data = clus->data;
        for (i1 = 0; i1 < j; i1++) {
          clus_data[i1] = 0.0;
        }

        /* 'mapTF2_ls:889' I_2_new(1:a-1) = I_2(1:a-1); */
        if (1.0 > a - 1.0) {
          loop_ub = 0;
        } else {
          loop_ub = (int)(a - 1.0);
        }

        for (i1 = 0; i1 < loop_ub; i1++) {
          clus_data[i1] = I_2_data[i1];
        }

        /* 'mapTF2_ls:890' I_2_new(a:a+1) = [info(ii), info(ii)]; */
        clus_data[(int)a - 1] = lasso_weight_data[b_i];
        clus_data[(int)(a + 1.0) - 1] = lasso_weight_data[b_i];

        /* 'mapTF2_ls:892' LEN_2 = LEN_2_new; */
        i1 = LEN_2->size[0] * LEN_2->size[1];
        LEN_2->size[0] = f->size[0];
        LEN_2->size[1] = 1;
        emxEnsureCapacity_real_T(LEN_2, i1);
        vec_data = LEN_2->data;
        loop_ub = f->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          vec_data[i1] = f_data[i1];
        }

        /* 'mapTF2_ls:893' I_2 = I_2_new; */
        i1 = I_2->size[0] * I_2->size[1];
        I_2->size[0] = clus->size[0];
        I_2->size[1] = 1;
        emxEnsureCapacity_real_T(I_2, i1);
        I_2_data = I_2->data;
        loop_ub = clus->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          I_2_data[i1] = clus_data[i1];
        }

        /* 'mapTF2_ls:894' a = a+1; */
        a++;

        /* 'mapTF2_ls:895' a = a+1; */
        a++;
      }
    }
  }

  emxFree_cell_wrap_4(&b_new_PWM2);
  emxFree_cell_wrap_4(&c_cur_PWM_tmp);
  emxFree_cell_wrap_4(&new_PWM2);
  emxFree_cell_wrap_4(&cur_PWM);
  emxFree_real_T(&vec2);
  emxFree_real_T(&uid);
  emxFree_real_T(&w);

  /* 'mapTF2_ls:899' num2 = length(strfind(fileread(memefn),'MOTIF')); */
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
      j = 1;
      while ((j <= 5) && (cur_line_data[(b_i + j) - 1] == b_cv[j - 1])) {
        j++;
      }

      if (j > 5) {
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

  /* 'mapTF2_ls:900' [p,names] = getmotif(memefn,1:num2); */
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

  emxInit_cell_wrap_2(&names);
  b_getmotif(memefn, y, cur_PWM_tmp, names);
  names_data = names->data;

  /* 'mapTF2_ls:901' [p,info,lenvec] = trim_pwm(p,0.25); */
  trim_pwm(cur_PWM_tmp, lasso_weight, f);
  f_data = f->data;
  cur_PWM_tmp_data = cur_PWM_tmp->data;

  /* 'mapTF2_ls:902' fid = fopen([ofn '_motifs.out'], 'w'); */
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

  /* 'mapTF2_ls:903' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF2_ls:904' fprintf("ERROR: Cannot create the combined motif file\n"); */
    printf("ERROR: Cannot create the combined motif file\n");
    fflush(stdout);
  }

  /* 'mapTF2_ls:906' a = 1; */
  a = 1.0;

  /* 'mapTF2_ls:907' for i = 1:length(PWM2) */
  match_idx = PWM2->size[0];
  if (match_idx < 1) {
    match_idx = 1;
  }

  if (PWM2->size[0] == 0) {
    match_idx = 0;
  }

  if (0 <= match_idx - 1) {
    unnamed_idx_0 = cur_PWM_tmp->size[0] + 1;
    i2 = cur_PWM_tmp->size[0];
    b_loop_ub = f->size[0];
  }

  emxInit_cell_wrap_4(&d_cur_PWM_tmp, 1);
  for (b_i = 0; b_i < match_idx; b_i++) {
    /* 'mapTF2_ls:908' p_len = length(p); */
    /* 'mapTF2_ls:909' cur_PWM_tmp = cell(p_len+1,1); */
    i = d_cur_PWM_tmp->size[0];
    d_cur_PWM_tmp->size[0] = unnamed_idx_0;
    emxEnsureCapacity_cell_wrap_4(d_cur_PWM_tmp, i);
    P_data = d_cur_PWM_tmp->data;
    for (i = 0; i < unnamed_idx_0; i++) {
      P_data[i].f1->size[0] = 0;
      P_data[i].f1->size[1] = 0;
    }

    /* 'mapTF2_ls:910' cur_PWM_tmp = coder.nullcopy(cur_PWM_tmp); */
    i = b_cur_PWM_tmp->size[0];
    b_cur_PWM_tmp->size[0] = d_cur_PWM_tmp->size[0];
    emxEnsureCapacity_cell_wrap_4(b_cur_PWM_tmp, i);
    P_data = b_cur_PWM_tmp->data;

    /* 'mapTF2_ls:911' cur_PWM_tmp{1} = PWM2{i}; */
    i = P_data[0].f1->size[0] * P_data[0].f1->size[1];
    P_data[0].f1->size[0] = PWM2_data[b_i].f1->size[0];
    P_data[0].f1->size[1] = PWM2_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(P_data[0].f1, i);
    loop_ub = PWM2_data[b_i].f1->size[0] * PWM2_data[b_i].f1->size[1];
    for (i = 0; i < loop_ub; i++) {
      P_data[0].f1->data[i] = PWM2_data[b_i].f1->data[i];
    }

    /* 'mapTF2_ls:912' for cur_idx=1:p_len */
    for (cur_idx = 0; cur_idx < i2; cur_idx++) {
      /* 'mapTF2_ls:913' cur_PWM_tmp{cur_idx+1} = p{cur_idx}; */
      i = P_data[cur_idx + 1].f1->size[0] * P_data[cur_idx + 1].f1->size[1];
      P_data[cur_idx + 1].f1->size[0] = cur_PWM_tmp_data[cur_idx].f1->size[0];
      P_data[cur_idx + 1].f1->size[1] = cur_PWM_tmp_data[cur_idx].f1->size[1];
      emxEnsureCapacity_real_T(P_data[cur_idx + 1].f1, i);
      loop_ub = cur_PWM_tmp_data[cur_idx].f1->size[0] * cur_PWM_tmp_data[cur_idx]
        .f1->size[1];
      for (i = 0; i < loop_ub; i++) {
        P_data[cur_idx + 1].f1->data[i] = cur_PWM_tmp_data[cur_idx].f1->data[i];
      }
    }

    /* 'mapTF2_ls:915' [ind, r] = ppmsim(cur_PWM_tmp, [LEN_2(i);lenvec]); */
    i = clus->size[0];
    clus->size[0] = f->size[0] + 1;
    emxEnsureCapacity_real_T(clus, i);
    clus_data = clus->data;
    clus_data[0] = vec_data[b_i];
    for (i = 0; i < b_loop_ub; i++) {
      clus_data[i + 1] = f_data[i];
    }

    ppmsim(b_cur_PWM_tmp, clus, &curr_pos, &idx);

    /*  [ind, r] = ppmsim([PWM2{i};p], [LEN_2(i);lenvec]); */
    /* 'mapTF2_ls:917' if r > 0.80 */
    if (idx > 0.8) {
      /* 'mapTF2_ls:918' fprintf(fid,'MOTIF %d\n%s\n%d\n', int32(a), names{ind}, int32(LEN_2(i))); */
      i = cur_line->size[0] * cur_line->size[1];
      cur_line->size[0] = 1;
      cur_line->size[1] = names_data[(int)curr_pos - 1].f1->size[1] + 1;
      emxEnsureCapacity_char_T(cur_line, i);
      cur_line_data = cur_line->data;
      loop_ub = names_data[(int)curr_pos - 1].f1->size[1];
      for (i = 0; i < loop_ub; i++) {
        cur_line_data[i] = names_data[(int)curr_pos - 1].f1->data[i];
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

      /* 'mapTF2_ls:919' a = a+1; */
      a++;

      /* 'mapTF2_ls:920' for j = 1:LEN_2(i) */
      i = (int)vec_data[b_i];
      if (0 <= i - 1) {
        d_NULL = NULL;
      }

      for (j = 0; j < i; j++) {
        /* 'mapTF2_ls:921' fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4)); */
        print_processing(PWM2_data[b_i].f1->data[j], PWM2_data[b_i].f1->data[j +
                         PWM2_data[b_i].f1->size[0]], PWM2_data[b_i].f1->data[j
                         + PWM2_data[b_i].f1->size[0] * 2], PWM2_data[b_i]
                         .f1->data[j + PWM2_data[b_i].f1->size[0] * 3],
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

      /* 'mapTF2_ls:923' fprintf(fid, '\n'); */
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
        /* 'mapTF2_ls:924' elseif I_2(i)/LEN_2(i) > 1 */
        /* 'mapTF2_ls:925' fprintf(fid,'MOTIF %d\n%s\n%d\n', int32(a), consen(PWM2{i}, LEN_2(i)), int32(LEN_2(i))); */
        b_NULL = NULL;
        getfilestar(fileid, &filestar, &hal);
        if (!(filestar == b_NULL)) {
          fprintf(filestar, "MOTIF %d\n%s\n%d\n", (int)a, "", (int)rt_roundd_snf
                  (d));
          if (hal) {
            fflush(filestar);
          }
        }

        /* 'mapTF2_ls:926' a = a+1; */
        a++;

        /* 'mapTF2_ls:927' for j = 1:LEN_2(i) */
        i = (int)d;
        if (0 <= (int)d - 1) {
          c_NULL = NULL;
        }

        for (j = 0; j < i; j++) {
          /* 'mapTF2_ls:928' fprintf(fid,'%0.3f %0.3f %0.3f %0.3f\n',PWM2{i}(j,1),PWM2{i}(j,2),PWM2{i}(j,3),PWM2{i}(j,4)); */
          print_processing(PWM2_data[b_i].f1->data[j], PWM2_data[b_i].f1->data[j
                           + PWM2_data[b_i].f1->size[0]], PWM2_data[b_i]
                           .f1->data[j + PWM2_data[b_i].f1->size[0] * 2],
                           PWM2_data[b_i].f1->data[j + PWM2_data[b_i].f1->size[0]
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

        /* 'mapTF2_ls:930' fprintf(fid, '\n'); */
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

  emxFree_cell_wrap_4(&d_cur_PWM_tmp);
  emxFree_cell_wrap_4(&cur_PWM_tmp);
  emxFree_cell_wrap_2(&names);
  emxFree_cell_wrap_4(&b_cur_PWM_tmp);
  emxFree_real_T(&I_2);
  emxFree_real_T(&LEN_2);
  emxFree_cell_wrap_4(&PWM2);
  emxFree_real_T(&f);
  emxFree_char_T(&cur_line);
  emxFree_real_T(&clus);

  /* 'mapTF2_ls:933' fclose(fid); */
  cfclose(fileid);
}

/*
 * function omat = scoreseqkmer(fn, NUM, ss, Smat, l,k, ofn,mat, dsvm)
 */
static void scoreseqkmer(const emxArray_real_T *ss, const emxArray_cell_wrap_4
  *Smat, double l, const emxArray_real_T *mat, const emxArray_real_T *dsvm,
  emxArray_real_T *omat)
{
  cell_wrap_15 varc[4];
  const cell_wrap_4 *Smat_data;
  emxArray_boolean_T *b_lab;
  emxArray_int32_T *f;
  emxArray_int32_T *r;
  emxArray_int8_T *varvec;
  emxArray_real_T *PWM;
  emxArray_real_T *a;
  emxArray_real_T *b;
  emxArray_real_T *evec;
  emxArray_real_T *lab;
  emxArray_real_T *loc;
  emxArray_real_T *path_ref;
  emxArray_real_T *r2;
  emxArray_real_T *scores;
  emxArray_real_T *u;
  emxArray_real_T *vec;
  emxArray_real_T *y;
  const double *dsvm_data;
  const double *mat_data;
  const double *ss_data;
  double b_a;
  double b_b;
  double b_x_idx_1_tmp;
  double d;
  double x_idx_1_tmp;
  double *PWM_data;
  double *a_data;
  double *evec_data;
  double *lab_data;
  double *loc_data;
  double *path_ref_data;
  double *scores_data;
  double *vec_data;
  double *y_data;
  int L;
  int b_i;
  int dimSize;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int j;
  int nx;
  int tmp2;
  int *f_data;
  int *r1;
  signed char *varvec_data;
  bool *b_lab_data;
  dsvm_data = dsvm->data;
  mat_data = mat->data;
  Smat_data = Smat->data;
  ss_data = ss->data;
  emxInit_real_T(&loc, 1);

  /* fn: fasta file */
  /* NUM: sequence number */
  /* l,k: gapped kmer parameters */
  /* ofn: output header */
  /* 'mapTF2_ls:273' L = length(ss); */
  L = ss->size[1];

  /* 'mapTF2_ls:274' nvar = numel(dsvm); */
  /* 'mapTF2_ls:275' loc = zeros(nvar, 1); */
  i = loc->size[0];
  loc->size[0] = dsvm->size[0];
  emxEnsureCapacity_real_T(loc, i);
  loc_data = loc->data;
  dimSize = dsvm->size[0];
  for (i = 0; i < dimSize; i++) {
    loc_data[i] = 0.0;
  }

  emxInit_int8_T(&varvec, 1);

  /* 'mapTF2_ls:276' varvec = zeros(nvar, 1); */
  i = varvec->size[0];
  varvec->size[0] = dsvm->size[0];
  emxEnsureCapacity_int8_T(varvec, i);
  varvec_data = varvec->data;

  /* 'mapTF2_ls:277' varc = cell(4,1); */
  /* 'mapTF2_ls:278' varc{1} = [2 3 4]; */
  /* 'mapTF2_ls:279' varc{2} = [1 3 4]; */
  /* 'mapTF2_ls:280' varc{3} = [1 2 4]; */
  /* 'mapTF2_ls:281' varc{4} = [1 2 3]; */
  varc[0].f1[0] = 2.0;
  varc[1].f1[0] = 1.0;
  varc[2].f1[0] = 1.0;
  varc[3].f1[0] = 1.0;
  varc[0].f1[1] = 3.0;
  varc[1].f1[1] = 3.0;
  varc[2].f1[1] = 2.0;
  varc[3].f1[1] = 2.0;
  varc[0].f1[2] = 4.0;
  varc[1].f1[2] = 4.0;
  varc[2].f1[2] = 4.0;
  varc[3].f1[2] = 3.0;

  /* 'mapTF2_ls:282' for i = 1:nvar */
  i = dsvm->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:283' loc(i) = floor((i-1)/3)+l; */
    loc_data[b_i] = floor((((double)b_i + 1.0) - 1.0) / 3.0) + l;

    /* 'mapTF2_ls:284' varvec(i) = varc{ss(loc(i))}(mod(i-1,3)+1); */
    varvec_data[b_i] = (signed char)varc[(int)ss_data[(int)loc_data[b_i] - 1] -
      1].f1[(int)(b_mod(((double)b_i + 1.0) - 1.0, 3.0) + 1.0) - 1];
  }

  emxInit_real_T(&lab, 1);

  /* 'mapTF2_ls:286' lab = mat(:,1); */
  /* 'mapTF2_ls:287' lab = [lab; zeros(l,1)]; */
  dimSize = mat->size[0];
  i = lab->size[0];
  lab->size[0] = mat->size[0] + (int)l;
  emxEnsureCapacity_real_T(lab, i);
  lab_data = lab->data;
  for (i = 0; i < dimSize; i++) {
    lab_data[i] = mat_data[i];
  }

  nx = (int)l;
  for (i = 0; i < nx; i++) {
    lab_data[i + dimSize] = 0.0;
  }

  emxInit_real_T(&path_ref, 1);

  /* 'mapTF2_ls:288' path_ref = mat(:,2); */
  /* 'mapTF2_ls:289' path_ref = [path_ref; zeros(l,1)]; */
  dimSize = mat->size[0];
  i = path_ref->size[0];
  path_ref->size[0] = mat->size[0] + (int)l;
  emxEnsureCapacity_real_T(path_ref, i);
  path_ref_data = path_ref->data;
  for (i = 0; i < dimSize; i++) {
    path_ref_data[i] = mat_data[i + mat->size[0]];
  }

  nx = (int)l;
  for (i = 0; i < nx; i++) {
    path_ref_data[i + dimSize] = 0.0;
  }

  emxInit_real_T(&PWM, 2);

  /* 'mapTF2_ls:290' PWM = log(mat(:,3:6)); */
  dimSize = mat->size[0];
  i = PWM->size[0] * PWM->size[1];
  PWM->size[0] = mat->size[0];
  PWM->size[1] = 4;
  emxEnsureCapacity_real_T(PWM, i);
  PWM_data = PWM->data;
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < dimSize; i1++) {
      PWM_data[i1 + PWM->size[0] * i] = mat_data[i1 + mat->size[0] * (i + 2)];
    }
  }

  nx = mat->size[0] << 2;
  for (tmp2 = 0; tmp2 < nx; tmp2++) {
    PWM_data[tmp2] = log(PWM_data[tmp2]);
  }

  emxInit_real_T(&evec, 1);

  /*  c = combnk(1:l,k); */
  /* 'mapTF2_ls:292' c = flip(nchoosek(1:l,k)); */
  /* 'mapTF2_ls:293' [c1,~] = size(c); */
  /* 'mapTF2_ls:294' O = ones(1,k); */
  /* 'mapTF2_ls:295' evec = zeros(nvar, 1); */
  i = evec->size[0];
  evec->size[0] = dsvm->size[0];
  emxEnsureCapacity_real_T(evec, i);
  evec_data = evec->data;
  dimSize = dsvm->size[0];
  for (i = 0; i < dimSize; i++) {
    evec_data[i] = 0.0;
  }

  emxInit_boolean_T(&b_lab, 1);

  /* 'mapTF2_ls:296' f = find(lab==0); */
  i = b_lab->size[0];
  b_lab->size[0] = lab->size[0];
  emxEnsureCapacity_boolean_T(b_lab, i);
  b_lab_data = b_lab->data;
  dimSize = lab->size[0];
  for (i = 0; i < dimSize; i++) {
    b_lab_data[i] = (lab_data[i] == 0.0);
  }

  emxInit_int32_T(&f, 1);
  emxInit_int32_T(&r, 1);
  eml_find(b_lab, r);
  r1 = r->data;
  i = f->size[0];
  f->size[0] = r->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  dimSize = r->size[0];
  for (i = 0; i < dimSize; i++) {
    f_data[i] = r1[i];
  }

  /* 'mapTF2_ls:297' GCvec = PWM(f(1),:); */
  /* 'mapTF2_ls:298' for i = 1:nvar */
  i = dsvm->size[0];
  emxInit_real_T(&vec, 2);
  vec_data = vec->data;
  emxInit_real_T(&scores, 1);
  scores_data = scores->data;
  emxInit_real_T(&y, 1);
  emxInit_real_T(&a, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:299' if lab(loc(i)) > 0 */
    i1 = (int)loc_data[b_i] - 1;
    d = lab_data[i1];
    if (d > 0.0) {
      /* 'mapTF2_ls:300' vec = (max([1 loc(i)-l+1]):min([L loc(i)+l-1])); */
      x_idx_1_tmp = (loc_data[b_i] - l) + 1.0;
      if (1.0 < x_idx_1_tmp) {
        b_a = x_idx_1_tmp;
      } else {
        b_a = 1.0;
      }

      b_x_idx_1_tmp = (loc_data[b_i] + l) - 1.0;
      if (L > b_x_idx_1_tmp) {
        b_b = b_x_idx_1_tmp;
      } else {
        b_b = L;
      }

      if (b_b < b_a) {
        vec->size[0] = 1;
        vec->size[1] = 0;
      } else if (floor(b_a) == b_a) {
        i2 = vec->size[0] * vec->size[1];
        vec->size[0] = 1;
        dimSize = (int)floor(b_b - b_a);
        vec->size[1] = dimSize + 1;
        emxEnsureCapacity_real_T(vec, i2);
        vec_data = vec->data;
        for (i2 = 0; i2 <= dimSize; i2++) {
          vec_data[i2] = b_a + (double)i2;
        }
      } else {
        eml_float_colon(b_a, b_b, vec);
        vec_data = vec->data;
      }

      /* 'mapTF2_ls:301' LAB = lab(loc(i)); */
      /* 'mapTF2_ls:302' path_vec = path_ref(vec); */
      /* 'mapTF2_ls:303' path_loc = path_ref(loc(i)); */
      /* 'mapTF2_ls:304' scores = zeros(length(vec),1); */
      i2 = scores->size[0];
      scores->size[0] = vec->size[1];
      emxEnsureCapacity_real_T(scores, i2);
      scores_data = scores->data;
      dimSize = vec->size[1];
      for (i2 = 0; i2 < dimSize; i2++) {
        scores_data[i2] = 0.0;
      }

      /* 'mapTF2_ls:305' V = exp(PWM(loc(i),varvec(i)))-exp(PWM(loc(i), ss(loc(i)))); */
      /* 'mapTF2_ls:306' for j = 1:length(vec) */
      i2 = vec->size[1];
      for (j = 0; j < i2; j++) {
        /* 'mapTF2_ls:307' if lab(vec(j)) == LAB && path_loc-path_vec(j)==loc(i)-vec(j) */
        i3 = (int)(unsigned int)vec_data[j] - 1;
        if ((lab_data[i3] == d) && (path_ref_data[i1] - path_ref_data[(int)
             vec_data[j] - 1] == loc_data[b_i] - (double)(unsigned int)
             vec_data[j])) {
          /* 'mapTF2_ls:308' scores(j) = PWM(vec(j),ss(vec(j))); */
          scores_data[j] = PWM_data[i3 + PWM->size[0] * ((int)ss_data[i3] - 1)];
        } else {
          /* 'mapTF2_ls:309' else */
          /* 'mapTF2_ls:310' scores(j) = GCvec(ss(vec(j))); */
          scores_data[j] = PWM_data[(f_data[0] + PWM->size[0] * ((int)ss_data[i3]
            - 1)) - 1];
        }
      }

      /* 'mapTF2_ls:313' if loc(i)-l+1 < 1 */
      if (x_idx_1_tmp < 1.0) {
        /* 'mapTF2_ls:314' a = l-loc(i); */
        b_a = l - loc_data[b_i];

        /* 'mapTF2_ls:315' scores(loc(i)) = 0; */
        scores_data[i1] = 0.0;
      } else {
        /* 'mapTF2_ls:316' else */
        /* 'mapTF2_ls:317' a = 0; */
        b_a = 0.0;

        /* 'mapTF2_ls:318' scores(l) = 0; */
        scores_data[(int)l - 1] = 0.0;
      }

      /* 'mapTF2_ls:321' if loc(i)+l-1 > L */
      if (b_x_idx_1_tmp > L) {
        /* 'mapTF2_ls:322' b = l-1+L-loc(i); */
        b_b = ((l - 1.0) + (double)L) - loc_data[b_i];
      } else {
        /* 'mapTF2_ls:323' else */
        /* 'mapTF2_ls:324' b = l-1; */
        b_b = l - 1.0;
      }

      /* 'mapTF2_ls:326' for j = a:b */
      i2 = (int)(b_b + (1.0 - b_a));
      for (j = 0; j < i2; j++) {
        x_idx_1_tmp = b_a + (double)j;

        /* 'mapTF2_ls:327' evec(i) = evec(i)+sum(exp(Smat{l-j}*scores(1+j-a:l+j-a))); */
        d = (x_idx_1_tmp + 1.0) - b_a;
        b_x_idx_1_tmp = (l + x_idx_1_tmp) - b_a;
        if (d > b_x_idx_1_tmp) {
          i3 = -1;
          i4 = -1;
        } else {
          i3 = (int)d - 2;
          i4 = (int)b_x_idx_1_tmp - 1;
        }

        nx = a->size[0] * a->size[1];
        tmp2 = (int)(l - x_idx_1_tmp) - 1;
        a->size[0] = Smat_data[tmp2].f1->size[0];
        a->size[1] = Smat_data[tmp2].f1->size[1];
        emxEnsureCapacity_real_T(a, nx);
        a_data = a->data;
        dimSize = Smat_data[tmp2].f1->size[0] * Smat_data[tmp2].f1->size[1];
        for (nx = 0; nx < dimSize; nx++) {
          a_data[nx] = Smat_data[tmp2].f1->data[nx];
        }

        nx = vec->size[0] * vec->size[1];
        vec->size[0] = 1;
        dimSize = i4 - i3;
        vec->size[1] = dimSize;
        emxEnsureCapacity_real_T(vec, nx);
        vec_data = vec->data;
        for (i4 = 0; i4 < dimSize; i4++) {
          vec_data[i4] = scores_data[(i3 + i4) + 1];
        }

        if ((Smat_data[tmp2].f1->size[0] == 0) || (Smat_data[tmp2].f1->size[1] ==
             0) || (dimSize == 0)) {
          i3 = y->size[0];
          y->size[0] = Smat_data[tmp2].f1->size[0];
          emxEnsureCapacity_real_T(y, i3);
          y_data = y->data;
          dimSize = Smat_data[tmp2].f1->size[0];
          for (i3 = 0; i3 < dimSize; i3++) {
            y_data[i3] = 0.0;
          }
        } else {
          i3 = y->size[0];
          y->size[0] = Smat_data[tmp2].f1->size[0];
          emxEnsureCapacity_real_T(y, i3);
          y_data = y->data;
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)
                      Smat_data[(int)(l - x_idx_1_tmp) - 1].f1->size[0],
                      (blasint)1, (blasint)Smat_data[(int)(l - x_idx_1_tmp) - 1]
                      .f1->size[1], 1.0, &a_data[0], (blasint)Smat_data[(int)(l
            - x_idx_1_tmp) - 1].f1->size[0], &vec_data[0], (blasint)1, 0.0,
                      &y_data[0], (blasint)Smat_data[(int)(l - x_idx_1_tmp) - 1]
                      .f1->size[0]);
        }

        nx = y->size[0];
        for (tmp2 = 0; tmp2 < nx; tmp2++) {
          y_data[tmp2] = exp(y_data[tmp2]);
        }

        evec_data[b_i] += blockedSummation(y, y->size[0]);

        /* for jj = 1:c1 */
        /*     if ismember(l-j , c(jj,:)); */
        /*         evec(i) = evec(i) + prod(scores(c(jj,:)+j-a)); */
        /*     end */
        /* end */
      }

      /* 'mapTF2_ls:334' evec(i) = evec(i)*V; */
      evec_data[b_i] *= exp(PWM_data[i1 + PWM->size[0] * (varvec_data[b_i] - 1)])
        - exp(PWM_data[i1 + PWM->size[0] * ((int)ss_data[i1] - 1)]);
    }
  }

  emxFree_real_T(&a);
  emxFree_real_T(&vec);
  emxFree_real_T(&PWM);
  emxInit_real_T(&u, 1);

  /* 'mapTF2_ls:337' u = unique(lab); */
  unique_vector(lab, u);
  vec_data = u->data;

  /* 'mapTF2_ls:338' for i = 1:numel(u) */
  i = u->size[0];
  emxInit_real_T(&r2, 1);
  emxInit_real_T(&b, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:339' if u(i) ~= 0 */
    if (vec_data[b_i] != 0.0) {
      /* 'mapTF2_ls:340' d = diff(lab==u(i)); */
      dimSize = lab->size[0];
      i1 = b_lab->size[0];
      b_lab->size[0] = lab->size[0];
      emxEnsureCapacity_boolean_T(b_lab, i1);
      b_lab_data = b_lab->data;
      for (i1 = 0; i1 < dimSize; i1++) {
        b_lab_data[i1] = (lab_data[i1] == vec_data[b_i]);
      }

      dimSize = b_lab->size[0];
      if (b_lab->size[0] == 0) {
        path_ref->size[0] = 0;
      } else {
        nx = b_lab->size[0] - 1;
        if (nx > 1) {
          nx = 1;
        }

        if (nx < 1) {
          path_ref->size[0] = 0;
        } else {
          i1 = path_ref->size[0];
          path_ref->size[0] = b_lab->size[0] - 1;
          emxEnsureCapacity_real_T(path_ref, i1);
          path_ref_data = path_ref->data;
          if (b_lab->size[0] - 1 != 0) {
            nx = b_lab_data[0];
            for (L = 2; L <= dimSize; L++) {
              tmp2 = nx;
              nx = b_lab_data[L - 1];
              path_ref_data[L - 2] = nx - tmp2;
            }
          }
        }
      }

      /* 'mapTF2_ls:341' f = find(d==1); */
      i1 = b_lab->size[0];
      b_lab->size[0] = path_ref->size[0];
      emxEnsureCapacity_boolean_T(b_lab, i1);
      b_lab_data = b_lab->data;
      dimSize = path_ref->size[0];
      for (i1 = 0; i1 < dimSize; i1++) {
        b_lab_data[i1] = (path_ref_data[i1] == 1.0);
      }

      eml_find(b_lab, r);
      r1 = r->data;
      i1 = f->size[0];
      f->size[0] = r->size[0];
      emxEnsureCapacity_int32_T(f, i1);
      f_data = f->data;
      dimSize = r->size[0];
      for (i1 = 0; i1 < dimSize; i1++) {
        f_data[i1] = r1[i1];
      }

      /* 'mapTF2_ls:342' F = find(d==-1); */
      i1 = b_lab->size[0];
      b_lab->size[0] = path_ref->size[0];
      emxEnsureCapacity_boolean_T(b_lab, i1);
      b_lab_data = b_lab->data;
      dimSize = path_ref->size[0];
      for (i1 = 0; i1 < dimSize; i1++) {
        b_lab_data[i1] = (path_ref_data[i1] == -1.0);
      }

      eml_find(b_lab, r);
      r1 = r->data;
      i1 = path_ref->size[0];
      path_ref->size[0] = r->size[0];
      emxEnsureCapacity_real_T(path_ref, i1);
      path_ref_data = path_ref->data;
      dimSize = r->size[0];
      for (i1 = 0; i1 < dimSize; i1++) {
        path_ref_data[i1] = r1[i1];
      }

      /* 'mapTF2_ls:343' if numel(f) ~= numel(F) */
      if (f->size[0] != path_ref->size[0]) {
        /* 'mapTF2_ls:344' F = [F;numel(dsvm)]; */
        i1 = path_ref->size[0];
        i2 = path_ref->size[0];
        path_ref->size[0]++;
        emxEnsureCapacity_real_T(path_ref, i2);
        path_ref_data = path_ref->data;
        path_ref_data[i1] = dsvm->size[0];
      }

      /* 'mapTF2_ls:346' for j = 1:numel(f) */
      i1 = f->size[0];
      for (j = 0; j < i1; j++) {
        /* 'mapTF2_ls:347' a = []; */
        scores->size[0] = 0;

        /* 'mapTF2_ls:348' for jj = f(j)+1:F(j) */
        i2 = (int)path_ref_data[j] - f_data[j];
        for (tmp2 = 0; tmp2 < i2; tmp2++) {
          x_idx_1_tmp = ((double)f_data[j] + 1.0) + (double)tmp2;

          /* 'mapTF2_ls:349' a = [a;find(loc==jj)]; */
          dimSize = loc->size[0];
          i3 = b_lab->size[0];
          b_lab->size[0] = loc->size[0];
          emxEnsureCapacity_boolean_T(b_lab, i3);
          b_lab_data = b_lab->data;
          for (i3 = 0; i3 < dimSize; i3++) {
            b_lab_data[i3] = (loc_data[i3] == x_idx_1_tmp);
          }

          eml_find(b_lab, r);
          r1 = r->data;
          i3 = scores->size[0];
          i4 = scores->size[0];
          scores->size[0] += r->size[0];
          emxEnsureCapacity_real_T(scores, i4);
          scores_data = scores->data;
          dimSize = r->size[0];
          for (i4 = 0; i4 < dimSize; i4++) {
            scores_data[i3 + i4] = r1[i4];
          }
        }

        /* 'mapTF2_ls:351' evec(a) = evec(a)*(evec(a)'*dsvm(a))/(evec(a)'*evec(a)); */
        i2 = y->size[0];
        y->size[0] = scores->size[0];
        emxEnsureCapacity_real_T(y, i2);
        y_data = y->data;
        dimSize = scores->size[0];
        for (i2 = 0; i2 < dimSize; i2++) {
          y_data[i2] = evec_data[(int)scores_data[i2] - 1];
        }

        i2 = r2->size[0];
        r2->size[0] = scores->size[0];
        emxEnsureCapacity_real_T(r2, i2);
        a_data = r2->data;
        dimSize = scores->size[0];
        for (i2 = 0; i2 < dimSize; i2++) {
          a_data[i2] = evec_data[(int)scores_data[i2] - 1];
        }

        i2 = b->size[0];
        b->size[0] = scores->size[0];
        emxEnsureCapacity_real_T(b, i2);
        PWM_data = b->data;
        dimSize = scores->size[0];
        for (i2 = 0; i2 < dimSize; i2++) {
          PWM_data[i2] = dsvm_data[(int)scores_data[i2] - 1];
        }

        if (scores->size[0] < 1) {
          b_b = 0.0;
        } else {
          b_b = cblas_ddot((blasint)scores->size[0], &y_data[0], (blasint)1,
                           &PWM_data[0], (blasint)1);
        }

        i2 = b->size[0];
        b->size[0] = scores->size[0];
        emxEnsureCapacity_real_T(b, i2);
        PWM_data = b->data;
        dimSize = scores->size[0];
        for (i2 = 0; i2 < dimSize; i2++) {
          PWM_data[i2] = evec_data[(int)scores_data[i2] - 1];
        }

        if (scores->size[0] < 1) {
          x_idx_1_tmp = 0.0;
        } else {
          x_idx_1_tmp = cblas_ddot((blasint)scores->size[0], &a_data[0],
            (blasint)1, &PWM_data[0], (blasint)1);
        }

        i2 = y->size[0];
        y->size[0] = scores->size[0];
        emxEnsureCapacity_real_T(y, i2);
        y_data = y->data;
        dimSize = scores->size[0];
        for (i2 = 0; i2 < dimSize; i2++) {
          y_data[i2] = evec_data[(int)scores_data[i2] - 1] * b_b / x_idx_1_tmp;
        }

        dimSize = y->size[0];
        for (i2 = 0; i2 < dimSize; i2++) {
          evec_data[(int)scores_data[i2] - 1] = y_data[i2];
        }
      }
    }
  }

  emxFree_boolean_T(&b_lab);
  emxFree_int32_T(&r);
  emxFree_real_T(&y);
  emxFree_real_T(&b);
  emxFree_real_T(&r2);
  emxFree_real_T(&u);
  emxFree_real_T(&scores);
  emxFree_int32_T(&f);
  emxFree_real_T(&path_ref);
  emxFree_real_T(&lab);

  /* 'mapTF2_ls:355' omat = [loc varvec evec]; */
  i = omat->size[0] * omat->size[1];
  omat->size[0] = loc->size[0];
  omat->size[1] = 3;
  emxEnsureCapacity_real_T(omat, i);
  vec_data = omat->data;
  dimSize = loc->size[0];
  for (i = 0; i < dimSize; i++) {
    vec_data[i] = loc_data[i];
  }

  emxFree_real_T(&loc);
  dimSize = varvec->size[0];
  for (i = 0; i < dimSize; i++) {
    vec_data[i + omat->size[0]] = varvec_data[i];
  }

  emxFree_int8_T(&varvec);
  dimSize = evec->size[0];
  for (i = 0; i < dimSize; i++) {
    vec_data[i + omat->size[0] * 2] = evec_data[i];
  }

  emxFree_real_T(&evec);
}

/*
 * function [P,V,seqindmat,seqout,seq] = seq2pv(sfn, wfn, input_l)
 */
static void seq2pv(const emxArray_char_T *sfn, const emxArray_char_T *wfn,
                   double input_l, emxArray_cell_wrap_0 *P, emxArray_cell_wrap_0
                   *V, emxArray_cell_wrap_0 *seqindmat, emxArray_cell_wrap_1
                   *seqout, emxArray_cell_wrap_2 *seq)
{
  static const signed char mat[12] = { 1, 0, 0, 0, 2, 2, 1, 1, 3, 3, 3, 2 };

  cell_wrap_0 *P_data;
  cell_wrap_0 *seqindmat_data;
  cell_wrap_1 *seqout_data;
  cell_wrap_2 *seq_data;
  cell_wrap_2 *sequences_data;
  emxArray_cell_wrap_2 *sequences;
  emxArray_char_T *b_fileid;
  emxArray_char_T *c_fileid;
  emxArray_char_T *cur_alpha;
  emxArray_char_T *cur_line;
  emxArray_char_T *cur_seq;
  emxArray_real_T *O;
  emxArray_real_T *S;
  emxArray_real_T *alpha;
  emxArray_real_T *b_pow;
  emxArray_real_T *b_y;
  emxArray_real_T *p;
  emxArray_real_T *ss;
  emxArray_real_T *w;
  creal_T dc;
  double b_I;
  double curr_pos;
  double idx;
  double *O_data;
  double *S_data;
  double *alpha_data;
  double *p_data;
  double *pow_data;
  double *ss_data;
  double *w_data;
  double *y_data;
  unsigned int a;
  int b_i;
  unsigned int b_j;
  int b_j1;
  int b_loop_ub;
  int c_loop_ub;
  int exitg1;
  int i;
  int i1;
  int j;
  int j2;
  int loop_ub;
  int n;
  int nd2;
  int varargin_2;
  signed char fileid;
  bool exitg2;
  bool y;

  /* sfn: fasta file */
  /* wfn: kmer weight file */
  /* ofn: output header */
  /*  fid = fopen(wfn, 'r'); */
  /*  X = textscan(fid,'%s\t%f\n'); */
  /*  fclose(fid); */
  /* 'mapTF2_ls:538' fid = fopen(wfn, 'r'); */
  fileid = cfopen(wfn, "rb");

  /* 'mapTF2_ls:539' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF2_ls:540' fprintf("ERROR: Weight file cannot be opened.\n") */
    printf("ERROR: Weight file cannot be opened.\n");
    fflush(stdout);
  }

  /* 'mapTF2_ls:543' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF2_ls:544' idx=0; */
  idx = 0.0;

  /* 'mapTF2_ls:545' while ~feof(fid) */
  emxInit_char_T(&b_fileid, 2);
  do {
    exitg1 = 0;
    b_I = b_feof(fileid);
    if (!(b_I != 0.0)) {
      /* 'mapTF2_ls:546' idx=idx+1; */
      idx++;

      /* 'mapTF2_ls:547' fgetl(fid); */
      b_fgets(fileid, b_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_fileid);
  emxInit_cell_wrap_2(&sequences);

  /* 'mapTF2_ls:549' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF2_ls:550' sequences = cell(idx, 1); */
  n = (int)idx;
  i = sequences->size[0];
  sequences->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_2(sequences, i);
  sequences_data = sequences->data;
  for (i = 0; i < n; i++) {
    sequences_data[i].f1->size[0] = 1;
    sequences_data[i].f1->size[1] = 0;
  }

  emxInit_real_T(&alpha, 1);

  /* 'mapTF2_ls:551' sequences = coder.nullcopy(sequences); */
  /* 'mapTF2_ls:552' alpha = zeros(idx, 1); */
  i = alpha->size[0];
  alpha->size[0] = (int)idx;
  emxEnsureCapacity_real_T(alpha, i);
  alpha_data = alpha->data;
  for (i = 0; i < n; i++) {
    alpha_data[i] = 0.0;
  }

  /* 'mapTF2_ls:553' for cur_idx=1:idx */
  n = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&cur_seq, 2);
  emxInit_char_T(&cur_alpha, 2);
  exitg2 = false;
  while ((!exitg2) && (n <= (int)idx - 1)) {
    /* 'mapTF2_ls:554' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF2_ls:555' if cur_line == -1 */
    y = (cur_line->size[1] != 0);
    if (y) {
      y = (0 > cur_line->size[1] - 1);
    }

    if (y) {
      exitg2 = true;
    } else {
      /* 'mapTF2_ls:558' [cur_seq, cur_alpha] = strtok(cur_line, char(9)); */
      c_strtok(cur_line, cur_seq, cur_alpha);

      /* 'mapTF2_ls:559' alpha(cur_idx,1) = real(str2double(cur_alpha)); */
      dc = str2double(cur_alpha);
      alpha_data[n] = dc.re;

      /* 'mapTF2_ls:560' sequences{cur_idx} = (strip(cur_seq)); */
      strip(cur_seq, sequences_data[n].f1);
      n++;
    }
  }

  emxFree_char_T(&cur_alpha);
  emxFree_char_T(&cur_seq);

  /* 'mapTF2_ls:562' fclose(fid); */
  cfclose(fileid);

  /* 'mapTF2_ls:564' l = length(sequences{1}); */
  varargin_2 = sequences_data[0].f1->size[1];

  /* 'mapTF2_ls:565' if l ~= input_l */
  if (sequences_data[0].f1->size[1] != input_l) {
    /* 'mapTF2_ls:566' fprintf("ERROR: L must be the same as the length of k-mer in the weight file\n"); */
    printf("ERROR: L must be the same as the length of k-mer in the weight file\n");
    fflush(stdout);
  }

  emxInit_real_T(&w, 1);

  /* 'mapTF2_ls:570' w = zeros(4^l,1); */
  j2 = (int)rt_powd_snf(4.0, sequences_data[0].f1->size[1]);
  i = w->size[0];
  w->size[0] = (int)rt_powd_snf(4.0, sequences_data[0].f1->size[1]);
  emxEnsureCapacity_real_T(w, i);
  w_data = w->data;
  for (i = 0; i < j2; i++) {
    w_data[i] = 0.0;
  }

  /* 'mapTF2_ls:571' pow = (4.^(0:(l-1)))'; */
  emxInit_real_T(&b_y, 2);
  if (sequences_data[0].f1->size[1] - 1 < 0) {
    b_y->size[1] = 0;
  } else {
    i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = sequences_data[0].f1->size[1];
    emxEnsureCapacity_real_T(b_y, i);
    y_data = b_y->data;
    j2 = sequences_data[0].f1->size[1] - 1;
    for (i = 0; i <= j2; i++) {
      y_data[i] = i;
    }
  }

  i = b_y->size[0] * b_y->size[1];
  b_y->size[0] = 1;
  emxEnsureCapacity_real_T(b_y, i);
  y_data = b_y->data;
  j2 = b_y->size[1] - 1;
  for (i = 0; i <= j2; i++) {
    curr_pos = y_data[i];
    y_data[i] = rt_powd_snf(4.0, curr_pos);
  }

  emxInit_real_T(&b_pow, 1);
  i = b_pow->size[0];
  b_pow->size[0] = b_y->size[1];
  emxEnsureCapacity_real_T(b_pow, i);
  pow_data = b_pow->data;
  j2 = b_y->size[1];
  for (i = 0; i < j2; i++) {
    pow_data[i] = y_data[i];
  }

  /* 'mapTF2_ls:572' fprintf('calculating indices\n'); */
  printf("calculating indices\n");
  fflush(stdout);

  /* 'mapTF2_ls:573' for i = 1:numel(alpha) */
  i = alpha->size[0];
  emxInit_real_T(&ss, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:574' ss = letterconvert(sequences{i}); */
    letterconvert(sequences_data[b_i].f1, ss);
    ss_data = ss->data;

    /* 'mapTF2_ls:575' rs = 3-fliplr(ss); */
    /* 'mapTF2_ls:576' w(ss*pow+1) = alpha(i); */
    if (ss->size[1] < 1) {
      curr_pos = 0.0;
    } else {
      curr_pos = cblas_ddot((blasint)ss->size[1], &ss_data[0], (blasint)1,
                            &pow_data[0], (blasint)1);
    }

    w_data[(int)(curr_pos + 1.0) - 1] = alpha_data[b_i];

    /* 'mapTF2_ls:577' w(rs*pow+1) = alpha(i); */
    n = ss->size[1] - 1;
    nd2 = ss->size[1] >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = n - b_j1;
      curr_pos = ss_data[b_j1];
      ss_data[b_j1] = ss_data[j2];
      ss_data[j2] = curr_pos;
    }

    i1 = ss->size[0] * ss->size[1];
    ss->size[0] = 1;
    emxEnsureCapacity_real_T(ss, i1);
    ss_data = ss->data;
    j2 = ss->size[1] - 1;
    for (i1 = 0; i1 <= j2; i1++) {
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

  emxInit_real_T(&p, 1);

  /* 'mapTF2_ls:579' m = mean(alpha); */
  curr_pos = blockedSummation(alpha, alpha->size[0]) / (double)alpha->size[0];

  /* 'mapTF2_ls:580' s = std(alpha); */
  idx = b_std(alpha);

  /* 'mapTF2_ls:581' W = (1/2)*(1+erf((w-m)/s/sqrt(2))); */
  i = p->size[0];
  p->size[0] = w->size[0];
  emxEnsureCapacity_real_T(p, i);
  p_data = p->data;
  j2 = w->size[0];
  for (i = 0; i < j2; i++) {
    p_data[i] = (w_data[i] - curr_pos) / idx / 1.4142135623730951;
  }

  applyScalarFunction(p, alpha);
  alpha_data = alpha->data;
  j2 = alpha->size[0];
  for (i = 0; i < j2; i++) {
    alpha_data[i] = 0.5 * (alpha_data[i] + 1.0);
  }

  /*  seq = importdata(sfn); */
  /* 'mapTF2_ls:584' fid = fopen(sfn, 'r'); */
  fileid = cfopen(sfn, "rb");

  /* 'mapTF2_ls:585' if fid == -1 */
  if (fileid == -1) {
    /* 'mapTF2_ls:586' fprintf("ERROR: Sequence file (.fa or .fasta) cannot be opened.\n") */
    printf("ERROR: Sequence file (.fa or .fasta) cannot be opened.\n");
    fflush(stdout);
  }

  /* 'mapTF2_ls:589' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'mapTF2_ls:590' idx=0; */
  idx = 0.0;

  /* 'mapTF2_ls:591' while ~feof(fid) */
  emxInit_char_T(&c_fileid, 2);
  do {
    exitg1 = 0;
    b_I = b_feof(fileid);
    if (!(b_I != 0.0)) {
      /* 'mapTF2_ls:592' idx=idx+1; */
      idx++;

      /* 'mapTF2_ls:593' fgetl(fid); */
      b_fgets(fileid, c_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&c_fileid);

  /* 'mapTF2_ls:595' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'mapTF2_ls:596' seq = cell(idx, 1); */
  n = (int)idx;
  i = seq->size[0];
  seq->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_2(seq, i);
  seq_data = seq->data;
  for (i = 0; i < n; i++) {
    seq_data[i].f1->size[0] = 1;
    seq_data[i].f1->size[1] = 0;
  }

  /* 'mapTF2_ls:597' seq = coder.nullcopy(seq); */
  /* 'mapTF2_ls:598' for cur_idx=1:idx */
  n = 0;
  exitg2 = false;
  while ((!exitg2) && (n <= (int)idx - 1)) {
    /* 'mapTF2_ls:599' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'mapTF2_ls:600' if cur_line == -1 */
    y = (cur_line->size[1] != 0);
    if (y) {
      y = (0 > cur_line->size[1] - 1);
    }

    if (y) {
      exitg2 = true;
    } else {
      /* 'mapTF2_ls:603' seq{cur_idx} = (strip(cur_line)); */
      strip(cur_line, seq_data[n].f1);
      n++;
    }
  }

  emxFree_char_T(&cur_line);

  /* 'mapTF2_ls:605' fclose(fid); */
  cfclose(fileid);

  /* 'mapTF2_ls:607' n = length(seq)/2; */
  curr_pos = (double)seq->size[0] / 2.0;

  /* 'mapTF2_ls:608' seqout = cell(n,1); */
  /* 'mapTF2_ls:609' seqindmat = cell(n,1); */
  b_j1 = (int)curr_pos;
  i = seqindmat->size[0];
  seqindmat->size[0] = (int)curr_pos;
  emxEnsureCapacity_cell_wrap_0(seqindmat, i);
  seqindmat_data = seqindmat->data;
  for (i = 0; i < b_j1; i++) {
    seqindmat_data[i].f1->size[0] = 0;
  }

  /* 'mapTF2_ls:610' seqindmat = coder.nullcopy(seqindmat); */
  /* 'mapTF2_ls:611' fprintf('converting kmers to probabilities\n'); */
  printf("converting kmers to probabilities\n");
  fflush(stdout);

  /* 'mapTF2_ls:612' P = cell(n,1); */
  /* 'mapTF2_ls:613' for i = 1:n */
  i = seqout->size[0];
  seqout->size[0] = (int)curr_pos;
  emxEnsureCapacity_cell_wrap_1(seqout, i);
  seqout_data = seqout->data;
  i = P->size[0];
  P->size[0] = (int)curr_pos;
  emxEnsureCapacity_cell_wrap_0(P, i);
  P_data = P->data;
  if (0 <= (int)curr_pos - 1) {
    if (1 > sequences_data[0].f1->size[1]) {
      loop_ub = 0;
    } else {
      loop_ub = sequences_data[0].f1->size[1];
    }
  }

  for (b_i = 0; b_i < b_j1; b_i++) {
    /* 'mapTF2_ls:614' if mod(i,1000)==0 */
    if (fmod((double)b_i + 1.0, 1000.0) == 0.0) {
      /* 'mapTF2_ls:615' fprintf('%d sequences converted\n', int32(i)); */
      printf("%d sequences converted\n", b_i + 1);
      fflush(stdout);
    }

    /* 'mapTF2_ls:617' L = length(seq{2*i})-l+1; */
    n = ((b_i + 1) << 1) - 1;
    idx = (double)(seq_data[n].f1->size[1] - varargin_2) + 1.0;

    /* 'mapTF2_ls:618' seqindmat{i} = zeros(L,1); */
    nd2 = (int)idx;
    i = seqindmat_data[b_i].f1->size[0];
    seqindmat_data[b_i].f1->size[0] = (int)idx;
    emxEnsureCapacity_real_T(seqindmat_data[b_i].f1, i);
    for (i = 0; i < nd2; i++) {
      seqindmat_data[b_i].f1->data[i] = 0.0;
    }

    /* 'mapTF2_ls:619' ss = letterconvert(seq{2*i}); */
    letterconvert(seq_data[n].f1, ss);
    ss_data = ss->data;

    /* 'mapTF2_ls:620' seqout{i} = ss+1; */
    i = seqout_data[b_i].f1->size[0] * seqout_data[b_i].f1->size[1];
    seqout_data[b_i].f1->size[0] = 1;
    seqout_data[b_i].f1->size[1] = ss->size[1];
    emxEnsureCapacity_real_T(seqout_data[b_i].f1, i);
    j2 = ss->size[1];
    for (i = 0; i < j2; i++) {
      seqout_data[b_i].f1->data[i] = ss_data[i] + 1.0;
    }

    /* 'mapTF2_ls:621' p = zeros(L,1); */
    i = p->size[0];
    p->size[0] = (int)idx;
    emxEnsureCapacity_real_T(p, i);
    p_data = p->data;
    for (i = 0; i < nd2; i++) {
      p_data[i] = 0.0;
    }

    /* 'mapTF2_ls:622' I = ss(1:l)*pow; */
    i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = loop_ub;
    emxEnsureCapacity_real_T(b_y, i);
    y_data = b_y->data;
    for (i = 0; i < loop_ub; i++) {
      y_data[i] = ss_data[i];
    }

    if (loop_ub < 1) {
      b_I = 0.0;
    } else {
      b_I = cblas_ddot((blasint)loop_ub, &y_data[0], (blasint)1, &pow_data[0],
                       (blasint)1);
    }

    /* 'mapTF2_ls:623' seqindmat{i}(1) = I+1; */
    seqindmat_data[b_i].f1->data[0] = b_I + 1.0;

    /* 'mapTF2_ls:624' p(1) = W(I+1); */
    p_data[0] = alpha_data[(int)(b_I + 1.0) - 1];

    /* 'mapTF2_ls:625' for j = 2:L */
    i = (int)(idx + -1.0);
    for (j = 0; j < i; j++) {
      /* 'mapTF2_ls:626' I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1); */
      b_I = (b_I - ss_data[j]) / 4.0 + rt_powd_snf(4.0, (double)varargin_2 - 1.0)
        * ss_data[(int)((unsigned int)j + varargin_2)];

      /* 'mapTF2_ls:627' seqindmat{i}(j) = I+1; */
      seqindmat_data[b_i].f1->data[j + 1] = b_I + 1.0;

      /* 'mapTF2_ls:628' p(j) = W(I+1); */
      p_data[j + 1] = alpha_data[(int)(b_I + 1.0) - 1];
    }

    /* 'mapTF2_ls:630' P{i} = p; */
    i = P_data[b_i].f1->size[0];
    P_data[b_i].f1->size[0] = p->size[0];
    emxEnsureCapacity_real_T(P_data[b_i].f1, i);
    j2 = p->size[0];
    for (i = 0; i < j2; i++) {
      P_data[b_i].f1->data[i] = p_data[i];
    }
  }

  emxInit_real_T(&O, 2);

  /* 'mapTF2_ls:632' fprintf('Running dsvm\n') */
  printf("Running dsvm\n");
  fflush(stdout);

  /* 'mapTF2_ls:633' mat = [1 2 3;0 2 3;0 1 3;0 1 2]; */
  /* 'mapTF2_ls:634' O = ones(1,l); */
  i = O->size[0] * O->size[1];
  O->size[0] = 1;
  O->size[1] = sequences_data[0].f1->size[1];
  emxEnsureCapacity_real_T(O, i);
  O_data = O->data;
  j2 = sequences_data[0].f1->size[1];
  for (i = 0; i < j2; i++) {
    O_data[i] = 1.0;
  }

  /* 'mapTF2_ls:635' V = cell(n,1); */
  /* 'mapTF2_ls:636' for i = 1:n */
  i = V->size[0];
  V->size[0] = (int)curr_pos;
  emxEnsureCapacity_cell_wrap_0(V, i);
  seqindmat_data = V->data;
  if (0 <= (int)curr_pos - 1) {
    if (1 > sequences_data[0].f1->size[1]) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = sequences_data[0].f1->size[1];
    }
  }

  emxInit_real_T(&S, 2);
  for (b_i = 0; b_i < b_j1; b_i++) {
    /* 'mapTF2_ls:637' if mod(i,1000)==0 */
    if (fmod((double)b_i + 1.0, 1000.0) == 0.0) {
      /* 'mapTF2_ls:638' fprintf('%d sequences converted\n', int32(i)); */
      printf("%d sequences converted\n", b_i + 1);
      fflush(stdout);
    }

    /* 'mapTF2_ls:640' L = length(seq{2*i})-2*l+2; */
    n = ((b_i + 1) << 1) - 1;
    idx = ((double)seq_data[n].f1->size[1] - 2.0 * (double)varargin_2) + 2.0;

    /* 'mapTF2_ls:641' ss = letterconvert(seq{2*i}); */
    letterconvert(seq_data[n].f1, ss);
    ss_data = ss->data;

    /* 'mapTF2_ls:642' p = zeros(L+l-1,1); */
    j2 = (int)((idx + (double)varargin_2) - 1.0);
    i = p->size[0];
    p->size[0] = j2;
    emxEnsureCapacity_real_T(p, i);
    p_data = p->data;
    for (i = 0; i < j2; i++) {
      p_data[i] = 0.0;
    }

    /* 'mapTF2_ls:643' I = ss(1:l)*pow; */
    i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(b_y, i);
    y_data = b_y->data;
    for (i = 0; i < b_loop_ub; i++) {
      y_data[i] = ss_data[i];
    }

    nd2 = b_loop_ub;
    if (b_loop_ub < 1) {
      b_I = 0.0;
    } else {
      b_I = cblas_ddot((blasint)b_loop_ub, &y_data[0], (blasint)1, &pow_data[0],
                       (blasint)1);
    }

    /* 'mapTF2_ls:644' p(1) = w(I+1); */
    p_data[0] = w_data[(int)(b_I + 1.0) - 1];

    /* 'mapTF2_ls:645' for j = 2:L+l-1 */
    i = (int)(((idx + (double)sequences_data[0].f1->size[1]) - 1.0) + -1.0);
    for (j = 0; j < i; j++) {
      /* 'mapTF2_ls:646' I = (I-ss(j-1))/4+4^(l-1)*ss(j+l-1); */
      b_I = (b_I - ss_data[j]) / 4.0 + rt_powd_snf(4.0, (double)varargin_2 - 1.0)
        * ss_data[(int)((((double)j + 2.0) + (double)varargin_2) - 1.0) - 1];

      /* 'mapTF2_ls:647' p(j) = w(I+1); */
      p_data[j + 1] = w_data[(int)(b_I + 1.0) - 1];
    }

    /* 'mapTF2_ls:649' v = zeros(3*L,1); */
    j2 = (int)(3.0 * idx);
    i = alpha->size[0];
    alpha->size[0] = j2;
    emxEnsureCapacity_real_T(alpha, i);
    alpha_data = alpha->data;
    for (i = 0; i < j2; i++) {
      alpha_data[i] = 0.0;
    }

    /* 'mapTF2_ls:650' a = 1; */
    a = 1U;

    /* 'mapTF2_ls:651' for j = l:L */
    b_I = fmod(idx, 4.294967296E+9);
    if (b_I < 0.0) {
      n = -(int)(unsigned int)-b_I;
    } else {
      n = (int)(unsigned int)b_I;
    }

    i = n - varargin_2;
    if (0 <= i) {
      c_loop_ub = varargin_2;
      nd2 = varargin_2;
    }

    for (j = 0; j <= i; j++) {
      b_j = (unsigned int)varargin_2 + j;

      /* 'mapTF2_ls:652' S = ss(j-l+1:j+l-1); */
      b_I = ((double)b_j - (double)varargin_2) + 1.0;
      curr_pos = ((double)b_j + (double)varargin_2) - 1.0;
      if (b_I > curr_pos) {
        i1 = 0;
        loop_ub = 0;
      } else {
        i1 = (int)b_I - 1;
        loop_ub = (int)curr_pos;
      }

      n = S->size[0] * S->size[1];
      S->size[0] = 1;
      j2 = loop_ub - i1;
      S->size[1] = j2;
      emxEnsureCapacity_real_T(S, n);
      S_data = S->data;
      for (loop_ub = 0; loop_ub < j2; loop_ub++) {
        S_data[loop_ub] = ss_data[i1 + loop_ub];
      }

      /* 'mapTF2_ls:653' ref = O*p(j-l+1:j); */
      if (b_I > b_j) {
        i1 = 0;
        loop_ub = 0;
      } else {
        i1 = (int)b_I - 1;
        loop_ub = (int)b_j;
      }

      n = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      j2 = loop_ub - i1;
      b_y->size[1] = j2;
      emxEnsureCapacity_real_T(b_y, n);
      y_data = b_y->data;
      for (loop_ub = 0; loop_ub < j2; loop_ub++) {
        y_data[loop_ub] = p_data[i1 + loop_ub];
      }

      if (sequences_data[0].f1->size[1] < 1) {
        curr_pos = 0.0;
      } else {
        curr_pos = cblas_ddot((blasint)sequences_data[0].f1->size[1], &O_data[0],
                              (blasint)1, &y_data[0], (blasint)1);
      }

      /* 'mapTF2_ls:654' for ii = 1:3 */
      for (loop_ub = 0; loop_ub < 3; loop_ub++) {
        /* 'mapTF2_ls:655' S(l) = mat(ss(j)+1,ii); */
        S_data[varargin_2 - 1] = mat[((int)(ss_data[(int)b_j - 1] + 1.0) +
          (loop_ub << 2)) - 1];

        /* 'mapTF2_ls:656' I = S(1:l)*pow; */
        i1 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        b_y->size[1] = c_loop_ub;
        emxEnsureCapacity_real_T(b_y, i1);
        y_data = b_y->data;
        for (i1 = 0; i1 < c_loop_ub; i1++) {
          y_data[i1] = S_data[i1];
        }

        if (nd2 < 1) {
          b_I = 0.0;
        } else {
          b_I = cblas_ddot((blasint)varargin_2, &y_data[0], (blasint)1,
                           &pow_data[0], (blasint)1);
        }

        /* 'mapTF2_ls:657' v(a) = w(I+1); */
        n = (int)(a + loop_ub) - 1;
        alpha_data[n] = w_data[(int)(b_I + 1.0) - 1];

        /* 'mapTF2_ls:658' for jj = 2:l */
        for (j2 = 0; j2 <= c_loop_ub - 2; j2++) {
          /* 'mapTF2_ls:659' I = (I-S(jj-1))/4+4^(l-1)*S(jj+l-1); */
          b_I = (b_I - S_data[j2]) / 4.0 + rt_powd_snf(4.0, (double)varargin_2 -
            1.0) * S_data[(int)((unsigned int)j2 + varargin_2)];

          /* 'mapTF2_ls:660' v(a)=v(a)+w(I+1); */
          alpha_data[n] += w_data[(int)(b_I + 1.0) - 1];
        }

        /* 'mapTF2_ls:662' v(a) = v(a)-ref; */
        alpha_data[n] -= curr_pos;

        /* 'mapTF2_ls:663' a = a+1; */
      }

      a += 3U;
    }

    /* 'mapTF2_ls:666' V{i} = v; */
    i = seqindmat_data[b_i].f1->size[0];
    seqindmat_data[b_i].f1->size[0] = alpha->size[0];
    emxEnsureCapacity_real_T(seqindmat_data[b_i].f1, i);
    j2 = alpha->size[0];
    for (i = 0; i < j2; i++) {
      seqindmat_data[b_i].f1->data[i] = alpha_data[i];
    }
  }

  emxFree_real_T(&b_y);
  emxFree_real_T(&S);
  emxFree_real_T(&O);
  emxFree_real_T(&p);
  emxFree_real_T(&ss);
  emxFree_real_T(&b_pow);
  emxFree_real_T(&w);
  emxFree_real_T(&alpha);
  emxFree_cell_wrap_2(&sequences);
}

/*
 * function [pp, info, len] = trim_pwm(p,cut)
 */
static void trim_pwm(emxArray_cell_wrap_4 *p, emxArray_real_T *info,
                     emxArray_real_T *len)
{
  cell_wrap_4 *p_data;
  emxArray_real_T *b_mat;
  emxArray_real_T *mat;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  double *b_mat_data;
  double *info_data;
  double *len_data;
  double *mat_data;
  double *x_data;
  int b_i;
  int c_i;
  int i;
  int idx;
  int j;
  int nrows;
  int nx;
  int nxin;
  p_data = p->data;

  /* 'mapTF2_ls:1054' l = length(p); */
  /* 'mapTF2_ls:1055' info = zeros(l, 1); */
  nx = p->size[0];
  i = info->size[0];
  info->size[0] = nx;
  emxEnsureCapacity_real_T(info, i);
  info_data = info->data;

  /* 'mapTF2_ls:1056' len = zeros(l,1); */
  i = len->size[0];
  len->size[0] = nx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;

  /* 'mapTF2_ls:1057' for i = 1:l */
  i = p->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&vec, 1);
  emxInit_real_T(&x, 2);
  emxInit_real_T(&b_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:1058' mat = p{i}+(p{i}==0); */
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

    /* 'mapTF2_ls:1059' vec = 2+sum(mat.*log(mat)/log(2),2); */
    idx = x->size[0] * x->size[1];
    x->size[0] = mat->size[0];
    x->size[1] = mat->size[1];
    emxEnsureCapacity_real_T(x, idx);
    x_data = x->data;
    nx = mat->size[0] * mat->size[1];
    for (idx = 0; idx < nx; idx++) {
      x_data[idx] = mat_data[idx];
    }

    nx = mat->size[0] * mat->size[1];
    for (nrows = 0; nrows < nx; nrows++) {
      x_data[nrows] = log(x_data[nrows]);
    }

    if ((mat->size[0] == x->size[0]) && (mat->size[1] == x->size[1])) {
      idx = b_mat->size[0] * b_mat->size[1];
      b_mat->size[0] = mat->size[0];
      b_mat->size[1] = mat->size[1];
      emxEnsureCapacity_real_T(b_mat, idx);
      b_mat_data = b_mat->data;
      nx = mat->size[0] * mat->size[1];
      for (idx = 0; idx < nx; idx++) {
        b_mat_data[idx] = mat_data[idx] * x_data[idx] / 0.69314718055994529;
      }

      sum(b_mat, vec);
      mat_data = vec->data;
    } else {
      j_binary_expand_op(vec, mat, x);
      mat_data = vec->data;
    }

    nx = vec->size[0];
    for (idx = 0; idx < nx; idx++) {
      mat_data[idx] += 2.0;
    }

    /* 'mapTF2_ls:1060' while (vec(1) < cut || mean(vec(1:3)) < cut || mean(vec(2:4)) < cut) && length(vec) > 4 */
    while (((mat_data[0] < 0.25) || (((mat_data[0] + mat_data[1]) + mat_data[2])
             / 3.0 < 0.25) || (((mat_data[1] + mat_data[2]) + mat_data[3]) / 3.0
             < 0.25)) && (vec->size[0] > 4)) {
      /* 'mapTF2_ls:1061' p{i}(1,:) = []; */
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
        nx = 0;
      } else {
        nx++;
      }

      nxin = p_data[b_i].f1->size[1] - 1;
      for (idx = 0; idx <= nxin; idx++) {
        for (nrows = 0; nrows < nx; nrows++) {
          p_data[b_i].f1->data[nrows + nx * idx] = p_data[b_i].f1->data[nrows +
            p_data[b_i].f1->size[0] * idx];
        }
      }

      idx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
      p_data[b_i].f1->size[0] = nx;
      p_data[b_i].f1->size[1] = nxin + 1;
      emxEnsureCapacity_real_T(p_data[b_i].f1, idx);

      /* 'mapTF2_ls:1062' vec(1) = []; */
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

    /* 'mapTF2_ls:1064' while (vec(end) < cut || mean(vec(end-2:end)) < cut || mean(vec(end-3:end-1)) < cut) && length(vec) > 4 */
    while (((mat_data[vec->size[0] - 1] < 0.25) || (((mat_data[vec->size[0] - 3]
               + mat_data[vec->size[0] - 2]) + mat_data[vec->size[0] - 1]) / 3.0
             < 0.25) || (((mat_data[vec->size[0] - 4] + mat_data[vec->size[0] -
                           3]) + mat_data[vec->size[0] - 2]) / 3.0 < 0.25)) &&
           (vec->size[0] > 4)) {
      /* 'mapTF2_ls:1065' vec(end) = []; */
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

      /* 'mapTF2_ls:1066' p{i}(end,:) = []; */
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
          p_data[b_i].f1->data[nrows + nx * idx] = p_data[b_i].f1->data[nrows +
            p_data[b_i].f1->size[0] * idx];
        }
      }

      idx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
      p_data[b_i].f1->size[0] = nx;
      p_data[b_i].f1->size[1] = nxin + 1;
      emxEnsureCapacity_real_T(p_data[b_i].f1, idx);
    }

    /* 'mapTF2_ls:1068' info(i) = sum(vec); */
    info_data[b_i] = blockedSummation(vec, vec->size[0]);

    /* 'mapTF2_ls:1069' [len(i), ~] = size(p{i}); */
    len_data[b_i] = p_data[b_i].f1->size[0];
  }

  emxFree_real_T(&b_mat);
  emxFree_real_T(&x);
  emxFree_real_T(&vec);
  emxFree_real_T(&mat);

  /* 'mapTF2_ls:1071' pp = p; */
}

/*
 * function mapTF2(varargin)
 */
void mapTF2_ls(const emxArray_char_T *varargin_1, const emxArray_char_T
               *varargin_2, const emxArray_char_T *varargin_3, const
               emxArray_char_T *varargin_4, const emxArray_char_T *varargin_5,
               const emxArray_char_T *varargin_6, double varargin_7, double
               varargin_8, double varargin_9)
{
  cell_wrap_0 *P_data;
  cell_wrap_0 *V_data;
  cell_wrap_0 *seqindmat_data;
  cell_wrap_1 *ss_data;
  cell_wrap_2 *seq_data;
  cell_wrap_3 *PWM2_data;
  cell_wrap_3 *PWM_data;
  cell_wrap_3 *pwm_data;
  cell_wrap_4 *Smat_data;
  emxArray_boolean_T *b_f;
  emxArray_cell_wrap_0 *P;
  emxArray_cell_wrap_0 *V;
  emxArray_cell_wrap_0 *seqindmat;
  emxArray_cell_wrap_1 *ss;
  emxArray_cell_wrap_2 *names;
  emxArray_cell_wrap_2 *seq;
  emxArray_cell_wrap_3 *PWM;
  emxArray_cell_wrap_3 *PWM2;
  emxArray_cell_wrap_3 *lpwm;
  emxArray_cell_wrap_3 *p;
  emxArray_cell_wrap_3 *pwm;
  emxArray_cell_wrap_4 *Smat;
  emxArray_char_T *b_varargin_6;
  emxArray_int32_T *r2;
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
  double GC;
  double L;
  double d;
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
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int input_sizes_idx_0;
  int j;
  int loop_ub;
  int nx;
  int sizes_idx_0;
  int *r3;
  const char *varargin_6_data;
  signed char *IND_data;
  char *b_varargin_6_data;
  signed char *ss_onehot_data;
  bool empty_non_axis_sizes;
  bool *b_f_data;
  if (!isInitialized_mapTF2_ls) {
    mapTF2_ls_initialize();
  }

  varargin_6_data = varargin_6->data;

  /*  if nargin < 6 */
  /*      error('Need at least 6 inputs') */
  /*  end */
  /* 'mapTF2_ls:6' fn = varargin{1}; */
  /* 'mapTF2_ls:7' wfn = varargin{2}; */
  /* 'mapTF2_ls:8' mfn1 = varargin{3}; */
  /* 'mapTF2_ls:9' mfn2 = varargin{4}; */
  /* 'mapTF2_ls:10' memefn = varargin{5}; */
  /* 'mapTF2_ls:11' ofn = varargin{6}; */
  /* 'mapTF2_ls:12' l_svm = varargin{7}; */
  /* 'mapTF2_ls:13' k_svm = varargin{8}; */
  /* 'mapTF2_ls:14' nfrac = varargin{9}; */
  /* 'mapTF2_ls:16' if l_svm < k_svm */
  if (varargin_7 < varargin_8) {
    /* 'mapTF2_ls:17' fprintf('ERROR: L must be greater or equal to the K\n'); */
    printf("ERROR: L must be greater or equal to the K\n");
    fflush(stdout);
  }

  emxInit_cell_wrap_0(&P);
  emxInit_cell_wrap_0(&V);
  emxInit_cell_wrap_0(&seqindmat);
  emxInit_cell_wrap_1(&ss);
  emxInit_cell_wrap_2(&seq);

  /*  if nargin == 7 */
  /*      l_svm = varargin{7}; */
  /*      k_svm = varargin{8}; */
  /*      if l_svm < k_svm */
  /*          error('7th argument must be greater or equal to the 8th argument') */
  /*      end */
  /*  end */
  /* 'mapTF2_ls:28' fprintf('Processing Motifs\n'); */
  printf("Processing Motifs\n");
  fflush(stdout);

  /* 'mapTF2_ls:29' process_motifs(mfn1, mfn2, memefn, ofn) */
  process_motifs(varargin_3, varargin_4, varargin_5, varargin_6);

  /* 'mapTF2_ls:30' mfn = [ofn '_motifs.out']; */
  /* 'mapTF2_ls:31' [P,V, seqindmat, ss, seq] = seq2pv(fn, wfn, l_svm); */
  seq2pv(varargin_1, varargin_2, varargin_7, P, V, seqindmat, ss, seq);
  seq_data = seq->data;
  ss_data = ss->data;
  seqindmat_data = seqindmat->data;
  V_data = V->data;
  P_data = P->data;

  /* 'mapTF2_ls:32' GC = countGC(ss); */
  /* 'mapTF2_ls:670' GC = 0; */
  GC = 0.0;

  /* 'mapTF2_ls:671' L = 0; */
  L = 0.0;

  /* 'mapTF2_ls:672' for i = 1:length(s) */
  i = ss->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:673' S = s{i}; */
    /* 'mapTF2_ls:674' N = length(S); */
    /* 'mapTF2_ls:675' L = N + L; */
    i1 = ss_data[b_i].f1->size[1];
    L += (double)ss_data[b_i].f1->size[1];

    /* 'mapTF2_ls:676' for j = 1:N */
    for (j = 0; j < i1; j++) {
      /* 'mapTF2_ls:677' if S(j) == 2 || S(j) == 3 */
      d = ss_data[b_i].f1->data[j];
      if ((d == 2.0) || (d == 3.0)) {
        /* 'mapTF2_ls:678' GC = GC + 1; */
        GC++;
      }
    }
  }

  emxInit_char_T(&b_varargin_6, 2);

  /* 'mapTF2_ls:682' GC = GC/L; */
  GC /= L;

  /* 'mapTF2_ls:33' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  L = 0.5 - GC / 2.0;
  GCmat[0] = L;
  GCmat[1] = GC / 2.0;
  GCmat[2] = GC / 2.0;
  GCmat[3] = L;

  /* 'mapTF2_ls:34' b = 1; */
  GC = 1.0;

  /* 'mapTF2_ls:35' [p,names,len] = getMOTIF(mfn); */
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

  emxInit_cell_wrap_3(&PWM);
  emxInit_cell_wrap_3(&PWM2);
  emxInit_cell_wrap_3(&p);
  emxInit_cell_wrap_2(&names);
  emxInit_real_T(&len, 1);
  getMOTIF(b_varargin_6, p, names, len);
  f_data = len->data;
  pwm_data = p->data;

  /* 'mapTF2_ls:36' a = numel(len); */
  /* 'mapTF2_ls:37' PWM = cell(a,1); */
  nx = len->size[0];
  i = PWM->size[0];
  PWM->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_3(PWM, i);
  PWM_data = PWM->data;

  /* 'mapTF2_ls:38' PWM = coder.nullcopy(PWM); */
  /* 'mapTF2_ls:39' PWM2 = cell(a,1); */
  i = PWM2->size[0];
  PWM2->size[0] = len->size[0];
  emxEnsureCapacity_cell_wrap_3(PWM2, i);
  PWM2_data = PWM2->data;
  emxFree_char_T(&b_varargin_6);
  for (i = 0; i < nx; i++) {
    PWM_data[i].f1->size[0] = 0;
    PWM_data[i].f1->size[1] = 4;
    PWM2_data[i].f1->size[0] = 0;
    PWM2_data[i].f1->size[1] = 4;
  }

  emxInit_real_T(&LEN, 1);

  /* 'mapTF2_ls:40' PWM2 = coder.nullcopy(PWM2); */
  /* 'mapTF2_ls:42' LEN = zeros(a,1); */
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

  /* 'mapTF2_ls:43' LEN_2 = zeros(a,1); */
  /* 'mapTF2_ls:44' shift = zeros(a,1); */
  /* 'mapTF2_ls:45' for i = 1:a */
  i = len->size[0];
  i1 = shift->size[0];
  shift->size[0] = len->size[0];
  emxEnsureCapacity_real_T(shift, i1);
  shift_data = shift->data;
  i1 = LEN_2->size[0];
  LEN_2->size[0] = len->size[0];
  emxEnsureCapacity_real_T(LEN_2, i1);
  LEN_2_data = LEN_2->data;
  emxInit_real_T(&r, 2);
  emxInit_real_T(&r1, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:46' shift(i) = max([l_svm-len(i) 4]); */
    d = varargin_7 - f_data[b_i];
    if ((d < 4.0) || rtIsNaN(d)) {
      d = 4.0;
      shift_data[b_i] = 4.0;
    } else {
      shift_data[b_i] = d;
    }

    /* 'mapTF2_ls:47' PWM{i} = [repmat(GCmat,shift(i), 1); p{i} ;repmat(GCmat,shift(i), 1)]; */
    repmat(GCmat, d, r);
    all_pwm_data = r->data;
    repmat(GCmat, d, r1);
    pwm_prob_data = r1->data;
    loop_ub = pwm_data[b_i].f1->size[0];
    i1 = PWM_data[b_i].f1->size[0] * PWM_data[b_i].f1->size[1];
    PWM_data[b_i].f1->size[0] = (r->size[0] + pwm_data[b_i].f1->size[0]) +
      r1->size[0];
    PWM_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(PWM_data[b_i].f1, i1);
    j = r->size[0];
    nx = r1->size[0];
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < j; i2++) {
        PWM_data[b_i].f1->data[i2 + PWM_data[b_i].f1->size[0] * i1] =
          all_pwm_data[i2 + r->size[0] * i1];
      }

      for (i2 = 0; i2 < loop_ub; i2++) {
        PWM_data[b_i].f1->data[(i2 + r->size[0]) + PWM_data[b_i].f1->size[0] *
          i1] = pwm_data[b_i].f1->data[i2 + pwm_data[b_i].f1->size[0] * i1];
      }

      for (i2 = 0; i2 < nx; i2++) {
        PWM_data[b_i].f1->data[((i2 + r->size[0]) + pwm_data[b_i].f1->size[0]) +
          PWM_data[b_i].f1->size[0] * i1] = pwm_prob_data[i2 + r1->size[0] * i1];
      }
    }

    /* 'mapTF2_ls:48' PWM2{i} = p{i}; */
    i1 = PWM2_data[b_i].f1->size[0] * PWM2_data[b_i].f1->size[1];
    PWM2_data[b_i].f1->size[0] = pwm_data[b_i].f1->size[0];
    PWM2_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(PWM2_data[b_i].f1, i1);
    loop_ub = pwm_data[b_i].f1->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      PWM2_data[b_i].f1->data[i1] = pwm_data[b_i].f1->data[i1];
    }

    /* 'mapTF2_ls:49' LEN_2(i) = len(i); */
    LEN_2_data[b_i] = f_data[b_i];

    /* 'mapTF2_ls:50' LEN(i) = len(i)+2*shift(i)-l_svm+1; */
    d = ((f_data[b_i] + 2.0 * shift_data[b_i]) - varargin_7) + 1.0;
    LEN_data[b_i] = d;

    /* 'mapTF2_ls:51' for j = 1:LEN(i) */
    i1 = (int)d;
    for (j = 0; j < i1; j++) {
      /* 'mapTF2_ls:52' b = b+1; */
      GC++;
    }
  }

  emxFree_real_T(&r1);
  emxFree_real_T(&r);
  emxFree_cell_wrap_3(&p);
  emxInit_cell_wrap_3(&pwm);
  emxInit_cell_wrap_3(&lpwm);

  /* 'mapTF2_ls:56' pwm = cell(b,1); */
  nx = (int)GC;
  i = pwm->size[0];
  pwm->size[0] = (int)GC;
  emxEnsureCapacity_cell_wrap_3(pwm, i);
  pwm_data = pwm->data;

  /* 'mapTF2_ls:57' pwm = coder.nullcopy(pwm); */
  /* 'mapTF2_ls:58' lpwm = cell(b,1); */
  i = lpwm->size[0];
  lpwm->size[0] = (int)GC;
  emxEnsureCapacity_cell_wrap_3(lpwm, i);
  PWM2_data = lpwm->data;
  for (i = 0; i < nx; i++) {
    pwm_data[i].f1->size[0] = 0;
    pwm_data[i].f1->size[1] = 4;
    PWM2_data[i].f1->size[0] = 0;
    PWM2_data[i].f1->size[1] = 4;
  }

  /* 'mapTF2_ls:59' lpwm = coder.nullcopy(lpwm); */
  /* 'mapTF2_ls:60' lab = zeros(b,1); */
  /* 'mapTF2_ls:61' b = 1; */
  GC = 1.0;

  /* 'mapTF2_ls:62' for i = 1:a */
  i = len->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:63' for j = 1:LEN(i) */
    i1 = (int)LEN_data[b_i];
    for (j = 0; j < i1; j++) {
      /* 'mapTF2_ls:64' pwm{b} = PWM{i}(j:j+l_svm-1,:); */
      d = (((double)j + 1.0) + varargin_7) - 1.0;
      if ((double)j + 1.0 > d) {
        i2 = 0;
        i3 = 0;
      } else {
        i2 = j;
        i3 = (int)d;
      }

      loop_ub = i3 - i2;
      i3 = (int)(GC + (double)j) - 1;
      i4 = pwm_data[i3].f1->size[0] * pwm_data[i3].f1->size[1];
      pwm_data[(int)(GC + (double)j) - 1].f1->size[0] = loop_ub;
      pwm_data[(int)(GC + (double)j) - 1].f1->size[1] = 4;
      emxEnsureCapacity_real_T(pwm_data[(int)(GC + (double)j) - 1].f1, i4);
      for (i4 = 0; i4 < 4; i4++) {
        for (i5 = 0; i5 < loop_ub; i5++) {
          pwm_data[(int)(GC + (double)j) - 1].f1->data[i5 + pwm_data[i3]
            .f1->size[0] * i4] = PWM_data[b_i].f1->data[(i2 + i5) + PWM_data[b_i]
            .f1->size[0] * i4];
        }
      }

      /* 'mapTF2_ls:65' lpwm{b} = log((pwm{b}+10^-10)/(1+4*10^-10)); */
      i2 = PWM2_data[i3].f1->size[0] * PWM2_data[i3].f1->size[1];
      PWM2_data[(int)(GC + (double)j) - 1].f1->size[0] = pwm_data[i3].f1->size[0];
      PWM2_data[(int)(GC + (double)j) - 1].f1->size[1] = 4;
      emxEnsureCapacity_real_T(PWM2_data[(int)(GC + (double)j) - 1].f1, i2);
      loop_ub = pwm_data[i3].f1->size[0] * 4;
      for (i2 = 0; i2 < loop_ub; i2++) {
        PWM2_data[(int)(GC + (double)j) - 1].f1->data[i2] = (pwm_data[i3]
          .f1->data[i2] + 1.0E-10) / 1.0000000004;
      }

      nx = PWM2_data[i3].f1->size[0] << 2;
      for (sizes_idx_0 = 0; sizes_idx_0 < nx; sizes_idx_0++) {
        PWM2_data[(int)(GC + (double)j) - 1].f1->data[sizes_idx_0] = log
          (PWM2_data[i3].f1->data[sizes_idx_0]);
      }

      /* 'mapTF2_ls:66' lab(b) = i; */
      /* 'mapTF2_ls:67' b = b+1; */
    }

    GC += (double)(i1 - 1) + 1.0;
  }

  emxFree_cell_wrap_3(&pwm);
  emxFree_cell_wrap_3(&PWM);
  emxInit_real_T(&c, 2);
  emxInit_real_T(&a__3, 1);
  emxInit_real_T(&seqmat, 2);
  emxInit_real_T(&f, 1);
  emxInit_real_T(&pwm_prob, 2);

  /* 'mapTF2_ls:72' [c,~,~,~,~,rcnum] = genIndex(l_svm,k_svm,nfrac); */
  genIndex(varargin_7, varargin_8, varargin_9, seqmat, pwm_prob, f, a__3, c, &L);
  seqmat_data = seqmat->data;

  /* 'mapTF2_ls:73' c2 = c(1:numel(c)/k_svm-rcnum,:); */
  d = (double)(seqmat->size[0] * seqmat->size[1]) / varargin_8 - L;
  if (1.0 > d) {
    loop_ub = 0;
  } else {
    loop_ub = (int)d;
  }

  /* 'mapTF2_ls:74' c = [c;l_svm+1-fliplr(c2)]; */
  j = seqmat->size[1];
  i = pwm_prob->size[0] * pwm_prob->size[1];
  pwm_prob->size[0] = loop_ub;
  pwm_prob->size[1] = seqmat->size[1];
  emxEnsureCapacity_real_T(pwm_prob, i);
  pwm_prob_data = pwm_prob->data;
  for (i = 0; i < j; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      pwm_prob_data[i1 + pwm_prob->size[0] * i] = seqmat_data[i1 + seqmat->size
        [0] * i];
    }
  }

  fliplr(pwm_prob);
  pwm_prob_data = pwm_prob->data;
  loop_ub = pwm_prob->size[0] * pwm_prob->size[1];
  for (i = 0; i < loop_ub; i++) {
    pwm_prob_data[i] = (varargin_7 + 1.0) - pwm_prob_data[i];
  }

  if ((seqmat->size[0] != 0) && (seqmat->size[1] != 0)) {
    nx = seqmat->size[1];
  } else if ((pwm_prob->size[0] != 0) && (pwm_prob->size[1] != 0)) {
    nx = pwm_prob->size[1];
  } else {
    nx = seqmat->size[1];
    if (pwm_prob->size[1] > seqmat->size[1]) {
      nx = pwm_prob->size[1];
    }
  }

  empty_non_axis_sizes = (nx == 0);
  if (empty_non_axis_sizes || ((seqmat->size[0] != 0) && (seqmat->size[1] != 0)))
  {
    input_sizes_idx_0 = seqmat->size[0];
  } else {
    input_sizes_idx_0 = 0;
  }

  if (empty_non_axis_sizes || ((pwm_prob->size[0] != 0) && (pwm_prob->size[1] !=
        0))) {
    sizes_idx_0 = pwm_prob->size[0];
  } else {
    sizes_idx_0 = 0;
  }

  j = input_sizes_idx_0;
  input_sizes_idx_0 = sizes_idx_0;
  i = c->size[0] * c->size[1];
  c->size[0] = j + sizes_idx_0;
  c->size[1] = nx;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  for (i = 0; i < nx; i++) {
    for (i1 = 0; i1 < j; i1++) {
      c_data[i1 + c->size[0] * i] = seqmat_data[i1 + j * i];
    }
  }

  for (i = 0; i < nx; i++) {
    for (i1 = 0; i1 < sizes_idx_0; i1++) {
      c_data[(i1 + j) + c->size[0] * i] = pwm_prob_data[i1 + sizes_idx_0 * i];
    }
  }

  /* 'mapTF2_ls:75' C = numel(c)/k_svm; */
  L = (double)(c->size[0] * c->size[1]) / varargin_8;

  /* 'mapTF2_ls:76' seqmat = zeros(C,l_svm); */
  i = (int)L;
  i1 = seqmat->size[0] * seqmat->size[1];
  seqmat->size[0] = (int)L;
  nx = (int)varargin_7;
  seqmat->size[1] = (int)varargin_7;
  emxEnsureCapacity_real_T(seqmat, i1);
  seqmat_data = seqmat->data;
  loop_ub = (int)L * (int)varargin_7;
  for (i1 = 0; i1 < loop_ub; i1++) {
    seqmat_data[i1] = 0.0;
  }

  /* 'mapTF2_ls:77' for i = 1:C */
  emxInit_int32_T(&r2, 1);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'mapTF2_ls:78' seqmat(i,c(i,:)) = 1; */
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

  emxInit_cell_wrap_4(&Smat, 1);

  /* 'mapTF2_ls:80' Smat = cell(l_svm,1); */
  i = Smat->size[0];
  Smat->size[0] = (int)varargin_7;
  emxEnsureCapacity_cell_wrap_4(Smat, i);
  Smat_data = Smat->data;
  for (i = 0; i < nx; i++) {
    Smat_data[i].f1->size[0] = 0;
    Smat_data[i].f1->size[1] = 0;
  }

  /* 'mapTF2_ls:81' Smat = coder.nullcopy(Smat); */
  /* 'mapTF2_ls:82' for i = 1:l_svm */
  emxInit_boolean_T(&b_f, 1);
  for (b_i = 0; b_i < nx; b_i++) {
    /* 'mapTF2_ls:83' f = find(prod(c-i,2)==0); */
    i = pwm_prob->size[0] * pwm_prob->size[1];
    pwm_prob->size[0] = c->size[0];
    pwm_prob->size[1] = c->size[1];
    emxEnsureCapacity_real_T(pwm_prob, i);
    pwm_prob_data = pwm_prob->data;
    loop_ub = c->size[0] * c->size[1];
    for (i = 0; i < loop_ub; i++) {
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

    eml_find(b_f, r2);
    r3 = r2->data;
    i = f->size[0];
    f->size[0] = r2->size[0];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = r2->size[0];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = r3[i];
    }

    /* 'mapTF2_ls:84' Smat{i} = zeros(length(f), l_svm); */
    i = Smat_data[b_i].f1->size[0] * Smat_data[b_i].f1->size[1];
    Smat_data[b_i].f1->size[0] = f->size[0];
    Smat_data[b_i].f1->size[1] = (int)varargin_7;
    emxEnsureCapacity_real_T(Smat_data[b_i].f1, i);
    loop_ub = f->size[0] * (int)varargin_7;
    for (i = 0; i < loop_ub; i++) {
      Smat_data[b_i].f1->data[i] = 0.0;
    }

    /* 'mapTF2_ls:85' for j = 1:length(f) */
    i = f->size[0];
    for (j = 0; j < i; j++) {
      /* 'mapTF2_ls:86' Smat{i}(j,c(f(j),:)) = 1; */
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

  emxFree_boolean_T(&b_f);
  emxFree_int32_T(&r2);
  emxInit_int8_T(&IND, 1);

  /* 'mapTF2_ls:89' L = length(ss); */
  /* 'mapTF2_ls:90' B = b-1; */
  /* 'mapTF2_ls:91' maxnorm = zeros(B,1); */
  /* 'mapTF2_ls:92' minnorm = zeros(B,1); */
  /* 'mapTF2_ls:93' vec = zeros(l_svm,1); */
  /* 'mapTF2_ls:94' IND = zeros(4^l_svm,1); */
  L = rt_powd_snf(4.0, varargin_7);
  nx = (int)rt_powd_snf(4.0, varargin_7);
  i = IND->size[0];
  IND->size[0] = (int)L;
  emxEnsureCapacity_int8_T(IND, i);
  IND_data = IND->data;
  for (i = 0; i < nx; i++) {
    IND_data[i] = 0;
  }

  /* 'mapTF2_ls:95' kmat = zeros(B,4^l_svm); */
  i = c->size[0] * c->size[1];
  c->size[0] = (int)((unsigned int)GC - 1U);
  c->size[1] = (int)L;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  loop_ub = (int)((unsigned int)GC - 1U) * (int)L;
  for (i = 0; i < loop_ub; i++) {
    c_data[i] = 0.0;
  }

  emxInit_real_T(&maxnorm, 1);
  emxInit_real_T(&minnorm, 1);

  /* 'mapTF2_ls:96' for j = 1:B */
  i = (int)(unsigned int)GC;
  i1 = maxnorm->size[0];
  maxnorm->size[0] = (int)((unsigned int)GC - 1U);
  emxEnsureCapacity_real_T(maxnorm, i1);
  maxnorm_data = maxnorm->data;
  i1 = minnorm->size[0];
  minnorm->size[0] = (int)((unsigned int)GC - 1U);
  emxEnsureCapacity_real_T(minnorm, i1);
  minnorm_data = minnorm->data;
  emxInit_real_T(&x, 2);
  emxInit_real_T(&B, 2);
  emxInit_real_T(&b_lpwm, 2);
  for (j = 0; j <= i - 2; j++) {
    /* 'mapTF2_ls:97' vec = max(lpwm{j}'); */
    /* 'mapTF2_ls:98' vec2 = min(lpwm{j}'); */
    /* 'mapTF2_ls:99' maxnorm(j) = sum(exp(seqmat*vec')); */
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

    b_maximum(b_lpwm, B);
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
    for (sizes_idx_0 = 0; sizes_idx_0 < nx; sizes_idx_0++) {
      f_data[sizes_idx_0] = exp(f_data[sizes_idx_0]);
    }

    b_sum(x, B);
    shift_data = B->data;
    maxnorm_data[j] = shift_data[0];

    /* 'mapTF2_ls:100' minnorm(j) = sum(exp(seqmat*vec2')); */
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

    minimum(b_lpwm, B);
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
    for (sizes_idx_0 = 0; sizes_idx_0 < nx; sizes_idx_0++) {
      LEN_2_data[sizes_idx_0] = exp(LEN_2_data[sizes_idx_0]);
    }

    minnorm_data[j] = blockedSummation(a__3, a__3->size[0]);
  }

  emxFree_real_T(&b_lpwm);
  emxFree_real_T(&B);
  emxFree_real_T(&x);

  /* 'mapTF2_ls:102' dnorm = maxnorm-minnorm; */
  if (maxnorm->size[0] == minnorm->size[0]) {
    loop_ub = maxnorm->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      maxnorm_data[i1] -= minnorm_data[i1];
    }
  } else {
    minus(maxnorm, minnorm);
    maxnorm_data = maxnorm->data;
  }

  emxInit_real_T(&all_pwm, 3);

  /* 'mapTF2_ls:103' vec = zeros(l_svm,1); */
  /* 'mapTF2_ls:105' all_pwm = zeros(4, l_svm, b-1); */
  i1 = all_pwm->size[0] * all_pwm->size[1] * all_pwm->size[2];
  all_pwm->size[0] = 4;
  all_pwm->size[1] = (int)varargin_7;
  all_pwm->size[2] = (int)((unsigned int)GC - 1U);
  emxEnsureCapacity_real_T(all_pwm, i1);
  all_pwm_data = all_pwm->data;
  loop_ub = ((int)varargin_7 << 2) * (int)((unsigned int)GC - 1U);
  for (i1 = 0; i1 < loop_ub; i1++) {
    all_pwm_data[i1] = 0.0;
  }

  /* 'mapTF2_ls:106' for cur_idx=1:b-1 */
  for (sizes_idx_0 = 0; sizes_idx_0 <= i - 2; sizes_idx_0++) {
    /* 'mapTF2_ls:107' all_pwm(:,:,cur_idx) = lpwm{cur_idx}'; */
    loop_ub = PWM2_data[sizes_idx_0].f1->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      all_pwm_data[4 * i1 + 4 * all_pwm->size[1] * sizes_idx_0] =
        PWM2_data[sizes_idx_0].f1->data[i1];
      all_pwm_data[(4 * i1 + 4 * all_pwm->size[1] * sizes_idx_0) + 1] =
        PWM2_data[sizes_idx_0].f1->data[i1 + PWM2_data[sizes_idx_0].f1->size[0]];
      all_pwm_data[(4 * i1 + 4 * all_pwm->size[1] * sizes_idx_0) + 2] =
        PWM2_data[sizes_idx_0].f1->data[i1 + PWM2_data[sizes_idx_0].f1->size[0] *
        2];
      all_pwm_data[(4 * i1 + 4 * all_pwm->size[1] * sizes_idx_0) + 3] =
        PWM2_data[sizes_idx_0].f1->data[i1 + PWM2_data[sizes_idx_0].f1->size[0] *
        3];
    }
  }

  emxFree_cell_wrap_3(&lpwm);

  /* 'mapTF2_ls:110' fprintf('Mapping motifs\n'); */
  printf("Mapping motifs\n");
  fflush(stdout);

  /* 'mapTF2_ls:111' tic */
  tic();

  /* 'mapTF2_ls:112' for I = 1:length(ss) */
  i1 = ss->size[0];
  emxInit_int8_T(&ss_onehot, 2);
  emxInit_real_T(&b_B, 1);
  emxInit_real_T(&b_all_pwm, 3);
  for (b_I = 0; b_I < i1; b_I++) {
    /* 'mapTF2_ls:113' seq2 = ss{I}; */
    /* 'mapTF2_ls:115' ss_onehot = zeros(4, length(seq2)); */
    i2 = ss_onehot->size[0] * ss_onehot->size[1];
    ss_onehot->size[0] = 4;
    ss_onehot->size[1] = ss_data[b_I].f1->size[1];
    emxEnsureCapacity_int8_T(ss_onehot, i2);
    ss_onehot_data = ss_onehot->data;
    loop_ub = ss_data[b_I].f1->size[1] << 2;
    for (i2 = 0; i2 < loop_ub; i2++) {
      ss_onehot_data[i2] = 0;
    }

    /* 'mapTF2_ls:116' for idx=1:length(seq2) */
    i2 = ss_data[b_I].f1->size[1];
    for (sizes_idx_0 = 0; sizes_idx_0 < i2; sizes_idx_0++) {
      /* 'mapTF2_ls:117' ss_onehot(seq2(idx), idx) = 1; */
      ss_onehot_data[((int)ss_data[b_I].f1->data[sizes_idx_0] + 4 * sizes_idx_0)
        - 1] = 1;
    }

    /* 'mapTF2_ls:119' pwm_prob = zeros(b-1,length(seqindmat{I})); */
    i2 = pwm_prob->size[0] * pwm_prob->size[1];
    pwm_prob->size[0] = (int)((unsigned int)GC - 1U);
    pwm_prob->size[1] = seqindmat_data[b_I].f1->size[0];
    emxEnsureCapacity_real_T(pwm_prob, i2);
    pwm_prob_data = pwm_prob->data;
    loop_ub = (int)((unsigned int)GC - 1U) * seqindmat_data[b_I].f1->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      pwm_prob_data[i2] = 0.0;
    }

    /* 'mapTF2_ls:121' for i = 1:length(seqindmat{I}) */
    i2 = seqindmat_data[b_I].f1->size[0];
    for (b_i = 0; b_i < i2; b_i++) {
      /* 'mapTF2_ls:122' ind = seqindmat{I}(i); */
      /* 'mapTF2_ls:123' if IND(ind) == 0 */
      i3 = (int)seqindmat_data[b_I].f1->data[b_i] - 1;
      if (IND_data[i3] == 0) {
        /* 'mapTF2_ls:124' IND(ind) = 1; */
        IND_data[i3] = 1;

        /*              SEQ = seq2(i:i+l_svm-1); */
        /*              for j = 1:b-1 */
        /*                  for jj = 1:l_svm */
        /*                      vec(jj) = lpwm{j}(jj,SEQ(jj)); */
        /*                  end */
        /*  %                 kmat(j,ind) = sum(exp(seqmat*vec)); */
        /*              end */
        /* 'mapTF2_ls:132' seq3 = ss_onehot(:,i:i+l_svm-1); */
        d = (((double)b_i + 1.0) + varargin_7) - 1.0;
        if ((double)b_i + 1.0 > d) {
          i4 = 0;
          i5 = -1;
        } else {
          i4 = b_i;
          i5 = (int)d - 1;
        }

        /* 'mapTF2_ls:133' vec3 = reshape(nonzeros(all_pwm .* seq3), l_svm, b-1); */
        if (all_pwm->size[1] == (i5 - i4) + 1) {
          i5 = b_all_pwm->size[0] * b_all_pwm->size[1] * b_all_pwm->size[2];
          b_all_pwm->size[0] = 4;
          b_all_pwm->size[1] = all_pwm->size[1];
          loop_ub = all_pwm->size[2];
          b_all_pwm->size[2] = all_pwm->size[2];
          emxEnsureCapacity_real_T(b_all_pwm, i5);
          f_data = b_all_pwm->data;
          for (i5 = 0; i5 < loop_ub; i5++) {
            j = all_pwm->size[1];
            for (nx = 0; nx < j; nx++) {
              sizes_idx_0 = i4 + nx;
              f_data[4 * nx + 4 * b_all_pwm->size[1] * i5] = all_pwm_data[4 * nx
                + 4 * all_pwm->size[1] * i5] * (double)ss_onehot_data[4 *
                sizes_idx_0];
              f_data[(4 * nx + 4 * b_all_pwm->size[1] * i5) + 1] = all_pwm_data
                [(4 * nx + 4 * all_pwm->size[1] * i5) + 1] * (double)
                ss_onehot_data[4 * sizes_idx_0 + 1];
              f_data[(4 * nx + 4 * b_all_pwm->size[1] * i5) + 2] = all_pwm_data
                [(4 * nx + 4 * all_pwm->size[1] * i5) + 2] * (double)
                ss_onehot_data[4 * sizes_idx_0 + 2];
              f_data[(4 * nx + 4 * b_all_pwm->size[1] * i5) + 3] = all_pwm_data
                [(4 * nx + 4 * all_pwm->size[1] * i5) + 3] * (double)
                ss_onehot_data[4 * sizes_idx_0 + 3];
            }
          }

          nonzeros(b_all_pwm, f);
          f_data = f->data;
        } else {
          b_binary_expand_op(f, all_pwm, ss_onehot, i4, i5, i4 - 1);
          f_data = f->data;
        }

        /* 'mapTF2_ls:134' for j = 1:b-1 */
        if (0 <= (int)(unsigned int)GC - 2) {
          input_sizes_idx_0 = (int)varargin_7;
          i6 = (int)varargin_7;
          b_loop_ub = (int)varargin_7;
        }

        for (j = 0; j <= i - 2; j++) {
          /* 'mapTF2_ls:135' kmat(j,ind) = sum(exp(seqmat*vec3(:,j))); */
          i4 = b_B->size[0];
          b_B->size[0] = i6;
          emxEnsureCapacity_real_T(b_B, i4);
          shift_data = b_B->data;
          for (i4 = 0; i4 < b_loop_ub; i4++) {
            shift_data[i4] = f_data[i4 + input_sizes_idx_0 * j];
          }

          loop_ub = seqmat->size[0];
          if ((seqmat->size[0] == 0) || (seqmat->size[1] == 0) || ((int)
               varargin_7 == 0)) {
            i4 = a__3->size[0];
            a__3->size[0] = seqmat->size[0];
            emxEnsureCapacity_real_T(a__3, i4);
            LEN_2_data = a__3->data;
            for (i4 = 0; i4 < loop_ub; i4++) {
              LEN_2_data[i4] = 0.0;
            }
          } else {
            i4 = a__3->size[0];
            a__3->size[0] = seqmat->size[0];
            emxEnsureCapacity_real_T(a__3, i4);
            LEN_2_data = a__3->data;
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                        seqmat->size[0], (blasint)1, (blasint)seqmat->size[1],
                        1.0, &seqmat_data[0], (blasint)seqmat->size[0],
                        &shift_data[0], (blasint)(int)varargin_7, 0.0,
                        &LEN_2_data[0], (blasint)seqmat->size[0]);
          }

          nx = a__3->size[0];
          for (sizes_idx_0 = 0; sizes_idx_0 < nx; sizes_idx_0++) {
            LEN_2_data[sizes_idx_0] = exp(LEN_2_data[sizes_idx_0]);
          }

          c_data[j + c->size[0] * i3] = blockedSummation(a__3, a__3->size[0]);
        }

        /*              kmat(:,ind) = sum(exp(seqmat * vec3)); */
        /* 'mapTF2_ls:139' kmat(:,ind) = log((kmat(:,ind)-minnorm)./dnorm); */
        loop_ub = c->size[0];
        if (c->size[0] == 1) {
          j = minnorm->size[0];
        } else {
          j = c->size[0];
        }

        if ((c->size[0] == minnorm->size[0]) && (j == maxnorm->size[0])) {
          i4 = f->size[0];
          f->size[0] = c->size[0];
          emxEnsureCapacity_real_T(f, i4);
          f_data = f->data;
          for (i4 = 0; i4 < loop_ub; i4++) {
            f_data[i4] = (c_data[i4 + c->size[0] * i3] - minnorm_data[i4]) /
              maxnorm_data[i4];
          }
        } else {
          binary_expand_op(f, c, seqindmat, b_I, b_i, minnorm, maxnorm);
          f_data = f->data;
        }

        nx = f->size[0];
        for (sizes_idx_0 = 0; sizes_idx_0 < nx; sizes_idx_0++) {
          f_data[sizes_idx_0] = log(f_data[sizes_idx_0]);
        }

        loop_ub = f->size[0];
        for (i4 = 0; i4 < loop_ub; i4++) {
          c_data[i4 + c->size[0] * i3] = f_data[i4];
        }

        /* 'mapTF2_ls:140' pwm_prob(:,i) = kmat(:,ind); */
        loop_ub = c->size[0];
        for (i4 = 0; i4 < loop_ub; i4++) {
          pwm_prob_data[i4 + pwm_prob->size[0] * b_i] = c_data[i4 + c->size[0] *
            i3];
        }
      } else {
        /* 'mapTF2_ls:141' else */
        /* 'mapTF2_ls:142' pwm_prob(:,i) = kmat(:,ind); */
        loop_ub = c->size[0];
        for (i4 = 0; i4 < loop_ub; i4++) {
          pwm_prob_data[i4 + pwm_prob->size[0] * b_i] = c_data[i4 + c->size[0] *
            i3];
        }
      }
    }

    /* 'mapTF2_ls:145' MAPTF(fn, mfn, seq{I*2}, ss{I}, GC, pwm_prob, Smat, l_svm, k_svm, ofn, PWM, PWM2, pwm, lpwm, lab, LEN, LEN_2, shift, P{I}, V{I}, names, len, a,b, I) */
    i2 = f->size[0];
    f->size[0] = LEN->size[0];
    emxEnsureCapacity_real_T(f, i2);
    f_data = f->data;
    loop_ub = LEN->size[0] - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
      f_data[i2] = LEN_data[i2];
    }

    i2 = a__3->size[0];
    a__3->size[0] = P_data[b_I].f1->size[0];
    emxEnsureCapacity_real_T(a__3, i2);
    LEN_2_data = a__3->data;
    loop_ub = P_data[b_I].f1->size[0] - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
      LEN_2_data[i2] = P_data[b_I].f1->data[i2];
    }

    MAPTF(seq_data[((b_I + 1) << 1) - 1].f1, ss_data[b_I].f1, pwm_prob, Smat,
          varargin_7, varargin_6, PWM2, f, LEN_2, shift, a__3, V_data[b_I].f1,
          names, len->size[0], (double)b_I + 1.0);

    /* 'mapTF2_ls:146' if mod(I,100)==0 */
    if (b_mod((double)b_I + 1.0, 100.0) == 0.0) {
      /* 'mapTF2_ls:147' fprintf('%d out of %d done...\n', int32(I), int32(length(ss))); */
      printf("%d out of %d done...\n", b_I + 1, ss->size[0]);
      fflush(stdout);

      /* 'mapTF2_ls:148' toc */
      toc();
    }
  }

  emxFree_real_T(&b_all_pwm);
  emxFree_real_T(&b_B);
  emxFree_real_T(&len);
  emxFree_cell_wrap_2(&names);
  emxFree_cell_wrap_2(&seq);
  emxFree_cell_wrap_1(&ss);
  emxFree_cell_wrap_0(&seqindmat);
  emxFree_cell_wrap_0(&V);
  emxFree_cell_wrap_0(&P);
  emxFree_real_T(&pwm_prob);
  emxFree_int8_T(&ss_onehot);
  emxFree_real_T(&all_pwm);
  emxFree_int8_T(&IND);
  emxFree_real_T(&minnorm);
  emxFree_real_T(&maxnorm);
  emxFree_real_T(&f);
  emxFree_cell_wrap_4(&Smat);
  emxFree_real_T(&seqmat);
  emxFree_real_T(&a__3);
  emxFree_real_T(&c);
  emxFree_real_T(&shift);
  emxFree_real_T(&LEN_2);
  emxFree_real_T(&LEN);
  emxFree_cell_wrap_3(&PWM2);
}

/* End of code generation (mapTF2_ls.c) */
