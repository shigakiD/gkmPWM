/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWM.c
 *
 * Code generation for function 'gkmPWM'
 *
 */

/* Include files */
#include "gkmPWM.h"
#include "PWM2kmers.h"
#include "PWM2kmers_norc.h"
#include "blockedSummation.h"
#include "colon.h"
#include "combineVectorElements.h"
#include "corrcoef.h"
#include "diff.h"
#include "eig.h"
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
#include "getgkmcounts.h"
#include "getmotif.h"
#include "gkmPWM_data.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_initialize.h"
#include "gkmPWM_rtwutil.h"
#include "gkmPWM_types.h"
#include "minOrMax.h"
#include "mod.h"
#include "mpower.h"
#include "mtimes.h"
#include "nullAssignment.h"
#include "randperm.h"
#include "repmat.h"
#include "rot90.h"
#include "sort.h"
#include "spdiags.h"
#include "sprintf.h"
#include "std.h"
#include "str2double.h"
#include "strip.h"
#include "strtok.h"
#include "sum.h"
#include "tic.h"
#include "toc.h"
#include "cblas.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Function Declarations */
static void ab_binary_expand_op(emxArray_real_T *kweig, const emxArray_real_T
  *indvec, int iii, const emxArray_cell_wrap_14 *ktree2, double k);
static double adjust_PWM(emxArray_real_T *p, const double GC[4]);
static void avg_info(const emxArray_cell_wrap_2 *p, double l_svm,
                     emxArray_real_T *info);
static void b_plus(emxArray_real_T *b, const emxArray_real_T *res);
static void b_seed_kmers(const emxArray_char_T *fn, double num, const
  emxArray_cell_wrap_1 *ik, emxArray_cell_wrap_0 *p, emxArray_cell_wrap_1 *mat,
  emxArray_cell_wrap_2 *pwms, double *c);
static void bb_binary_expand_op(emxArray_real_T *kweig, const emxArray_real_T
  *indvec2, int iii, const emxArray_cell_wrap_14 *ktree2, double k);
static void binary_expand_op(emxArray_real_T *cfile, const emxArray_real_T
  *negvec, double y, double b);
static void cb_binary_expand_op(emxArray_real_T *mat, const emxArray_real_T *x);
static void createMEME(const emxArray_char_T *fileh_Value, const
  emxArray_cell_wrap_2 *PWM, const emxArray_char_T *memefile, double GC, const
  emxArray_real_T *C, double r, const emxArray_real_T *R, double rcorr, const
  emxArray_real_T *E, const emxArray_real_T *Rd);
static void db_binary_expand_op(emxArray_real_T *mat, const emxArray_cell_wrap_2
  *p, int i, int i2, int i3, int i4, int i5);
static void f_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3);
static void g_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3);
static void getEMprob_v3(const emxArray_real_T *PWM, const emxArray_real_T *res,
  const double negmat[16], const emxArray_cell_wrap_14 *poscell, const
  emxArray_real_T *rc, const emxArray_real_T *diffc, const emxArray_real_T *indc,
  const emxArray_real_T *indloc, const emxArray_real_T *xc, double reg, double
  l_svm, double k_svm, double rcnum, double RC, emxArray_real_T *kweig, double
  P_data[], int *P_size);
static void gkmPWM_lagrange(const emxArray_real_T *kweig, const double negmat[16],
  const emxArray_cell_wrap_2 *PWM, const emxArray_real_T *negvec, double n,
  double rcorr, double reg, double l_svm, double k_svm, double RC, const
  emxArray_real_T *rc, const emxArray_real_T *diffc, const emxArray_real_T *indc,
  const emxArray_real_T *xc, double rcnum, emxArray_cell_wrap_2 *b_PWM,
  emxArray_real_T *scorevec, emxArray_real_T *C, double *r, emxArray_real_T *R,
  emxArray_real_T *E, emxArray_real_T *Rd);
static void h_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3);
static void i_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3);
static void j_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3);
static void k_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3);
static void l_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, const emxArray_real_T *rs);
static void ls_kweigtree(const emxArray_real_T *mat, const double negmat[16],
  const emxArray_cell_wrap_14 *poscell, const emxArray_real_T *c, const
  emxArray_real_T *s, const emxArray_real_T *ind, const emxArray_real_T *indloc,
  const emxArray_real_T *x, double l, double k, double rcnum, emxArray_real_T
  *kweig);
static void ls_kweigtree_norc(const emxArray_real_T *mat, const double negmat[16],
  const emxArray_cell_wrap_14 *poscell, const emxArray_real_T *c, const
  emxArray_real_T *s, const emxArray_real_T *ind, const emxArray_real_T *indloc,
  const emxArray_real_T *x, double l, double k, emxArray_real_T *kweig);
static void m_binary_expand_op(emxArray_real_T *kmat, int i, const
  emxArray_real_T *y, const emxArray_real_T *negvec, double b);
static void minus(emxArray_real_T *res, const emxArray_real_T *kweig);
static void n_binary_expand_op(emxArray_real_T *res, const emxArray_real_T *f,
  const emxArray_real_T *C, const emxArray_real_T *ord, int ii);
static void o_binary_expand_op(emxArray_real_T *kmat, const emxArray_real_T *ord,
  int ii, const emxArray_real_T *f);
static void p_binary_expand_op(emxArray_real_T *kmat, const emxArray_real_T *f,
  int jj, const emxArray_real_T *y, const emxArray_real_T *negvec, double b);
static void seed_kmers(const emxArray_char_T *fn, double num,
  emxArray_cell_wrap_0 *p, emxArray_cell_wrap_1 *mat, emxArray_cell_wrap_2 *pwms,
  double *c);
static void w_binary_expand_op(double p_data[], int *p_size, const creal_T
  ps_data[], const int *ps_size, const creal_T t, const creal_T B_data[], const
  int B_size[2], const creal_T E);
static void x_binary_expand_op(creal_T x[4], const double MAT_data[], const int
  MAT_size[2], const creal_T E, const signed char b[4]);
static void y_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *b,
  const emxArray_real_T *A, int iindx);

/* Function Definitions */
static void ab_binary_expand_op(emxArray_real_T *kweig, const emxArray_real_T
  *indvec, int iii, const emxArray_cell_wrap_14 *ktree2, double k)
{
  const cell_wrap_14 *ktree2_data;
  emxArray_real_T *b_kweig;
  const double *indvec_data;
  double *b_kweig_data;
  double *kweig_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  ktree2_data = ktree2->data;
  indvec_data = indvec->data;
  kweig_data = kweig->data;
  emxInit_real_T(&b_kweig, 1);
  i = b_kweig->size[0];
  if (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] == 1) {
    b_kweig->size[0] = indvec->size[0];
  } else {
    b_kweig->size[0] = ktree2_data[(int)(k - 1.0) - 1].f1->size[0];
  }

  emxEnsureCapacity_real_T(b_kweig, i);
  b_kweig_data = b_kweig->data;
  stride_0_0 = (indvec->size[0] != 1);
  stride_1_0 = (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] != 1);
  if (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] == 1) {
    loop_ub = indvec->size[0];
  } else {
    loop_ub = ktree2_data[(int)(k - 1.0) - 1].f1->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_kweig_data[i] = kweig_data[((int)indvec_data[i * stride_0_0] + kweig->
      size[0] * (3 - iii)) - 1] + ktree2_data[(int)(k - 1.0) - 1].f1->data[i *
      stride_1_0];
  }

  loop_ub = b_kweig->size[0];
  for (i = 0; i < loop_ub; i++) {
    kweig_data[((int)indvec_data[i] + kweig->size[0] * (3 - iii)) - 1] =
      b_kweig_data[i];
  }

  emxFree_real_T(&b_kweig);
}

/*
 * function [pp, len] = adjust_PWM(p,GC)
 */
static double adjust_PWM(emxArray_real_T *p, const double GC[4])
{
  emxArray_real_T *mat;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  double b_vec_data;
  double len;
  double *mat_data;
  double *p_data;
  double *vec_data;
  int exitg1;
  int i;
  int k;
  int nx;
  int vstride;
  int xj;
  bool b_x[3];
  bool b;
  bool b2;
  bool guard1 = false;
  bool guard2 = false;
  bool guard3 = false;
  bool guard4 = false;
  p_data = p->data;
  emxInit_real_T(&mat, 2);

  /* extends or truncates the PWM based on information.  I try to do it intelligently, so you may disagree on the condition required for adjustment. */
  /* 'gkmPWM:924' [len,~] = size(p); */
  len = p->size[0];

  /* 'gkmPWM:925' info = zeros(len, 1); */
  /* 'gkmPWM:926' cut = 0.2; */
  /* maximum information of a column to truncate */
  /* 'gkmPWM:927' ext = 0.7; */
  /* minimum information needed to extend */
  /* The rest of the code is easy enough to read through quickly */
  /* 'gkmPWM:929' mat = p+(p==0); */
  xj = mat->size[0] * mat->size[1];
  mat->size[0] = p->size[0];
  mat->size[1] = 4;
  emxEnsureCapacity_real_T(mat, xj);
  mat_data = mat->data;
  nx = p->size[0] * 4;
  for (xj = 0; xj < nx; xj++) {
    mat_data[xj] = p_data[xj] + (double)(p_data[xj] == 0.0);
  }

  emxInit_real_T(&x, 2);

  /* 'gkmPWM:930' vec = 2+sum(mat.*log(mat)/log(2),2); */
  xj = x->size[0] * x->size[1];
  x->size[0] = mat->size[0];
  x->size[1] = 4;
  emxEnsureCapacity_real_T(x, xj);
  vec_data = x->data;
  nx = mat->size[0] * 4;
  for (xj = 0; xj < nx; xj++) {
    vec_data[xj] = mat_data[xj];
  }

  nx = mat->size[0] << 2;
  for (k = 0; k < nx; k++) {
    vec_data[k] = log(vec_data[k]);
  }

  if (mat->size[0] == x->size[0]) {
    nx = mat->size[0] * 4;
    xj = mat->size[0] * mat->size[1];
    mat->size[1] = 4;
    emxEnsureCapacity_real_T(mat, xj);
    mat_data = mat->data;
    for (xj = 0; xj < nx; xj++) {
      mat_data[xj] = mat_data[xj] * vec_data[xj] / 0.69314718055994529;
    }
  } else {
    cb_binary_expand_op(mat, x);
    mat_data = mat->data;
  }

  emxFree_real_T(&x);
  emxInit_real_T(&vec, 1);
  vec_data = vec->data;
  if (mat->size[0] == 0) {
    vec->size[0] = 0;
  } else {
    vstride = mat->size[0];
    xj = vec->size[0];
    vec->size[0] = mat->size[0];
    emxEnsureCapacity_real_T(vec, xj);
    vec_data = vec->data;
    for (xj = 0; xj < vstride; xj++) {
      vec_data[xj] = mat_data[xj];
    }

    for (k = 0; k < 3; k++) {
      nx = (k + 1) * vstride;
      for (xj = 0; xj < vstride; xj++) {
        vec_data[xj] += mat_data[nx + xj];
      }
    }
  }

  nx = vec->size[0];
  for (xj = 0; xj < nx; xj++) {
    vec_data[xj] += 2.0;
  }

  /* 'gkmPWM:931' b = true; */
  b = true;

  /* 'gkmPWM:932' b2 = true; */
  b2 = true;

  /* 'gkmPWM:933' while ((vec(1) < cut && max(vec(2:3)) <= ext) || mean(vec(1:3) < cut)) && len > 12 */
  do {
    exitg1 = 0;
    guard1 = false;
    guard2 = false;
    guard3 = false;
    guard4 = false;
    if (vec_data[0] < 0.2) {
      if (vec_data[1] < vec_data[2]) {
        b_vec_data = vec_data[2];
      } else {
        b_vec_data = vec_data[1];
      }

      if (b_vec_data <= 0.7) {
        guard3 = true;
      } else {
        guard4 = true;
      }
    } else {
      guard4 = true;
    }

    if (guard4) {
      for (xj = 0; xj < 3; xj++) {
        b_x[xj] = (vec_data[xj] < 0.2);
      }

      b_vec_data = b_x[0];
      for (k = 0; k < 2; k++) {
        b_vec_data += (double)b_x[k + 1];
      }

      if (b_vec_data / 3.0 != 0.0) {
        guard3 = true;
      } else {
        exitg1 = 1;
      }
    }

    if (guard3) {
      if (len > 12.0) {
        /* 'gkmPWM:934' p(1,:) = []; */
        nx = p->size[0] - 2;
        vstride = p->size[0] - 1;
        for (k = 0; k < 4; k++) {
          for (i = 0; i < vstride; i++) {
            p_data[i + p->size[0] * k] = p_data[(i + p->size[0] * k) + 1];
          }
        }

        if (1 > vstride) {
          nx = 0;
        } else {
          nx++;
        }

        for (xj = 0; xj < 4; xj++) {
          for (vstride = 0; vstride < nx; vstride++) {
            p_data[vstride + nx * xj] = p_data[vstride + p->size[0] * xj];
          }
        }

        xj = p->size[0] * p->size[1];
        p->size[0] = nx;
        p->size[1] = 4;
        emxEnsureCapacity_real_T(p, xj);
        p_data = p->data;

        /* 'gkmPWM:935' vec(1) = []; */
        nx = vec->size[0];
        vstride = vec->size[0] - 1;
        for (k = 0; k < vstride; k++) {
          vec_data[k] = vec_data[k + 1];
        }

        xj = vec->size[0];
        if (1 > vstride) {
          vec->size[0] = 0;
        } else {
          vec->size[0] = nx - 1;
        }

        emxEnsureCapacity_real_T(vec, xj);
        vec_data = vec->data;

        /* 'gkmPWM:936' len = len-1; */
        len--;

        /* 'gkmPWM:937' b = false; */
        b = false;

        /* 'gkmPWM:938' if ((vec(end) < cut && max(vec(end-2:end-1)) <= ext) || mean(vec(end-2:end) < cut)) && len > 12 */
        if (vec_data[vec->size[0] - 1] < 0.2) {
          if (vec_data[vec->size[0] - 3] < vec_data[vec->size[0] - 2]) {
            b_vec_data = vec_data[vec->size[0] - 2];
          } else {
            b_vec_data = vec_data[vec->size[0] - 3];
          }

          if (b_vec_data <= 0.7) {
            guard1 = true;
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      } else {
        exitg1 = 1;
      }
    }

    if (guard2) {
      for (xj = 0; xj < 3; xj++) {
        b_x[xj] = (vec_data[(xj + vec->size[0]) - 3] < 0.2);
      }

      b_vec_data = b_x[0];
      for (k = 0; k < 2; k++) {
        b_vec_data += (double)b_x[k + 1];
      }

      if (b_vec_data / 3.0 != 0.0) {
        guard1 = true;
      }
    }

    if (guard1 && (len > 12.0)) {
      /* 'gkmPWM:939' vec(end) = []; */
      xj = vec->size[0];
      nx = vec->size[0];
      vstride = vec->size[0] - 1;
      for (k = xj; k <= vstride; k++) {
        vec_data[k - 1] = vec_data[k];
      }

      xj = vec->size[0];
      if (1 > vstride) {
        vec->size[0] = 0;
      } else {
        vec->size[0] = nx - 1;
      }

      emxEnsureCapacity_real_T(vec, xj);
      vec_data = vec->data;

      /* 'gkmPWM:940' p(end,:) = []; */
      xj = p->size[0];
      nx = p->size[0] - 2;
      vstride = p->size[0] - 1;
      for (k = 0; k < 4; k++) {
        for (i = xj; i <= vstride; i++) {
          p_data[(i + p->size[0] * k) - 1] = p_data[i + p->size[0] * k];
        }
      }

      if (1 > vstride) {
        nx = 0;
      } else {
        nx++;
      }

      for (xj = 0; xj < 4; xj++) {
        for (vstride = 0; vstride < nx; vstride++) {
          p_data[vstride + nx * xj] = p_data[vstride + p->size[0] * xj];
        }
      }

      xj = p->size[0] * p->size[1];
      p->size[0] = nx;
      p->size[1] = 4;
      emxEnsureCapacity_real_T(p, xj);
      p_data = p->data;

      /* 'gkmPWM:941' len = len-1; */
      len--;

      /* 'gkmPWM:942' b2 = false; */
      b2 = false;
    }
  } while (exitg1 == 0);

  /* 'gkmPWM:945' if b && min(vec(1:2)) > ext && len < 20 */
  if (b) {
    if (vec_data[0] > vec_data[1]) {
      b_vec_data = vec_data[1];
    } else {
      b_vec_data = vec_data[0];
    }

    if ((b_vec_data > 0.7) && (len < 20.0)) {
      /* 'gkmPWM:946' mat = [GC ; mat]; */
      /* 'gkmPWM:947' p = [GC;p]; */
      xj = mat->size[0] * mat->size[1];
      mat->size[0] = p->size[0] + 1;
      mat->size[1] = 4;
      emxEnsureCapacity_real_T(mat, xj);
      mat_data = mat->data;
      nx = p->size[0];
      for (xj = 0; xj < 4; xj++) {
        mat_data[mat->size[0] * xj] = GC[xj];
        for (vstride = 0; vstride < nx; vstride++) {
          mat_data[(vstride + mat->size[0] * xj) + 1] = p_data[vstride + p->
            size[0] * xj];
        }
      }

      xj = p->size[0] * p->size[1];
      p->size[0] = mat->size[0];
      p->size[1] = 4;
      emxEnsureCapacity_real_T(p, xj);
      p_data = p->data;
      nx = mat->size[0];
      for (xj = 0; xj < 4; xj++) {
        for (vstride = 0; vstride < nx; vstride++) {
          p_data[vstride + p->size[0] * xj] = mat_data[vstride + mat->size[0] *
            xj];
        }
      }

      /* 'gkmPWM:948' len = len+1; */
      len++;
    }
  }

  /* 'gkmPWM:950' while ((vec(end) < cut && max(vec(end-2:end-1)) <= ext) || mean(vec(end-2:end) < cut)) && len > 12 */
  do {
    exitg1 = 0;
    guard1 = false;
    guard2 = false;
    if (vec_data[vec->size[0] - 1] < 0.2) {
      if (vec_data[vec->size[0] - 3] < vec_data[vec->size[0] - 2]) {
        b_vec_data = vec_data[vec->size[0] - 2];
      } else {
        b_vec_data = vec_data[vec->size[0] - 3];
      }

      if (b_vec_data <= 0.7) {
        guard1 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2) {
      for (xj = 0; xj < 3; xj++) {
        b_x[xj] = (vec_data[(xj + vec->size[0]) - 3] < 0.2);
      }

      b_vec_data = b_x[0];
      for (k = 0; k < 2; k++) {
        b_vec_data += (double)b_x[k + 1];
      }

      if (b_vec_data / 3.0 != 0.0) {
        guard1 = true;
      } else {
        exitg1 = 1;
      }
    }

    if (guard1) {
      if (len > 12.0) {
        /* 'gkmPWM:951' vec(end) = []; */
        xj = vec->size[0];
        nx = vec->size[0];
        vstride = vec->size[0] - 1;
        for (k = xj; k <= vstride; k++) {
          vec_data[k - 1] = vec_data[k];
        }

        xj = vec->size[0];
        if (1 > vstride) {
          vec->size[0] = 0;
        } else {
          vec->size[0] = nx - 1;
        }

        emxEnsureCapacity_real_T(vec, xj);
        vec_data = vec->data;

        /* 'gkmPWM:952' p(end,:) = []; */
        xj = p->size[0];
        nx = p->size[0] - 2;
        vstride = p->size[0] - 1;
        for (k = 0; k < 4; k++) {
          for (i = xj; i <= vstride; i++) {
            p_data[(i + p->size[0] * k) - 1] = p_data[i + p->size[0] * k];
          }
        }

        if (1 > vstride) {
          nx = 0;
        } else {
          nx++;
        }

        for (xj = 0; xj < 4; xj++) {
          for (vstride = 0; vstride < nx; vstride++) {
            p_data[vstride + nx * xj] = p_data[vstride + p->size[0] * xj];
          }
        }

        xj = p->size[0] * p->size[1];
        p->size[0] = nx;
        p->size[1] = 4;
        emxEnsureCapacity_real_T(p, xj);
        p_data = p->data;

        /* 'gkmPWM:953' len = len-1; */
        len--;

        /* 'gkmPWM:954' b2 = false; */
        b2 = false;
      } else {
        exitg1 = 1;
      }
    }
  } while (exitg1 == 0);

  /* 'gkmPWM:956' if b2 && min(vec(end-1:end)) > ext && len < 20 */
  if (b2) {
    if (vec_data[vec->size[0] - 2] > vec_data[vec->size[0] - 1]) {
      b_vec_data = vec_data[vec->size[0] - 1];
    } else {
      b_vec_data = vec_data[vec->size[0] - 2];
    }

    if ((b_vec_data > 0.7) && (len < 20.0)) {
      /* 'gkmPWM:957' p = [p;GC]; */
      xj = mat->size[0] * mat->size[1];
      mat->size[0] = p->size[0] + 1;
      mat->size[1] = 4;
      emxEnsureCapacity_real_T(mat, xj);
      mat_data = mat->data;
      nx = p->size[0];
      for (xj = 0; xj < 4; xj++) {
        for (vstride = 0; vstride < nx; vstride++) {
          mat_data[vstride + mat->size[0] * xj] = p_data[vstride + p->size[0] *
            xj];
        }
      }

      for (xj = 0; xj < 4; xj++) {
        mat_data[p->size[0] + mat->size[0] * xj] = GC[xj];
      }

      xj = p->size[0] * p->size[1];
      p->size[0] = mat->size[0];
      p->size[1] = 4;
      emxEnsureCapacity_real_T(p, xj);
      p_data = p->data;
      nx = mat->size[0];
      for (xj = 0; xj < 4; xj++) {
        for (vstride = 0; vstride < nx; vstride++) {
          p_data[vstride + p->size[0] * xj] = mat_data[vstride + mat->size[0] *
            xj];
        }
      }

      /* 'gkmPWM:958' len = len+1; */
      len++;
    }
  }

  emxFree_real_T(&vec);
  emxFree_real_T(&mat);

  /* 'gkmPWM:960' pp = p; */
  return len;
}

/*
 * function info = avg_info(p,l_svm)
 */
static void avg_info(const emxArray_cell_wrap_2 *p, double l_svm,
                     emxArray_real_T *info)
{
  const cell_wrap_2 *p_data;
  emxArray_real_T *mat;
  emxArray_real_T *x;
  emxArray_real_T *y;
  double varargin_1;
  double *info_data;
  double *mat_data;
  double *x_data;
  double *y_data;
  int b_i;
  int i;
  int i1;
  int k;
  int nx;
  int vstride;
  int xj;
  p_data = p->data;

  /* 'gkmPWM:963' info = zeros(length(p),1); */
  i = info->size[0];
  info->size[0] = p->size[0];
  emxEnsureCapacity_real_T(info, i);
  info_data = info->data;
  nx = p->size[0];
  for (i = 0; i < nx; i++) {
    info_data[i] = 0.0;
  }

  /* 'gkmPWM:964' for i = 1:length(p) */
  i = p->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&x, 2);
  emxInit_real_T(&y, 1);
  y_data = y->data;
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:965' [l,~] = size(p{i}); */
    /* 'gkmPWM:966' mat = p{i}(l_svm:end-l_svm+1,:)+(p{i}(l_svm:end-l_svm+1,:)==0.0); */
    varargin_1 = ((double)p_data[b_i].f1->size[0] - l_svm) + 1.0;
    if (l_svm > varargin_1) {
      i1 = 0;
      xj = 0;
    } else {
      i1 = (int)l_svm - 1;
      xj = (int)varargin_1;
    }

    varargin_1 = ((double)p_data[b_i].f1->size[0] - l_svm) + 1.0;
    if (l_svm > varargin_1) {
      vstride = 0;
      k = 0;
    } else {
      vstride = (int)l_svm - 1;
      k = (int)varargin_1;
    }

    nx = xj - i1;
    if (nx == k - vstride) {
      xj = mat->size[0] * mat->size[1];
      mat->size[0] = nx;
      mat->size[1] = 4;
      emxEnsureCapacity_real_T(mat, xj);
      mat_data = mat->data;
      for (xj = 0; xj < 4; xj++) {
        for (k = 0; k < nx; k++) {
          mat_data[k + mat->size[0] * xj] = p_data[b_i].f1->data[(i1 + k) +
            p_data[b_i].f1->size[0] * xj] + (double)(p_data[b_i].f1->data
            [(vstride + k) + p_data[b_i].f1->size[0] * xj] == 0.0);
        }
      }
    } else {
      db_binary_expand_op(mat, p, b_i, i1, xj - 1, vstride, k - 1);
      mat_data = mat->data;
    }

    /* 'gkmPWM:967' info(i) = sum(2+sum(mat.*log(max(mat,0))/log(2),2)); */
    i1 = x->size[0] * x->size[1];
    x->size[0] = mat->size[0];
    x->size[1] = 4;
    emxEnsureCapacity_real_T(x, i1);
    x_data = x->data;
    nx = mat->size[0] * 4;
    for (i1 = 0; i1 < nx; i1++) {
      varargin_1 = mat_data[i1];
      x_data[i1] = fmax(varargin_1, 0.0);
    }

    nx = x->size[0] << 2;
    for (k = 0; k < nx; k++) {
      x_data[k] = log(x_data[k]);
    }

    if (mat->size[0] == x->size[0]) {
      nx = mat->size[0] * 4;
      i1 = mat->size[0] * mat->size[1];
      mat->size[1] = 4;
      emxEnsureCapacity_real_T(mat, i1);
      mat_data = mat->data;
      for (i1 = 0; i1 < nx; i1++) {
        mat_data[i1] = mat_data[i1] * x_data[i1] / 0.69314718055994529;
      }
    } else {
      cb_binary_expand_op(mat, x);
      mat_data = mat->data;
    }

    if (mat->size[0] == 0) {
      y->size[0] = 0;
    } else {
      vstride = mat->size[0];
      i1 = y->size[0];
      y->size[0] = mat->size[0];
      emxEnsureCapacity_real_T(y, i1);
      y_data = y->data;
      for (xj = 0; xj < vstride; xj++) {
        y_data[xj] = mat_data[xj];
      }

      for (k = 0; k < 3; k++) {
        nx = (k + 1) * vstride;
        for (xj = 0; xj < vstride; xj++) {
          y_data[xj] += mat_data[nx + xj];
        }
      }
    }

    nx = y->size[0];
    for (i1 = 0; i1 < nx; i1++) {
      y_data[i1] += 2.0;
    }

    info_data[b_i] = blockedSummation(y, y->size[0]);
  }

  emxFree_real_T(&y);
  emxFree_real_T(&x);
  emxFree_real_T(&mat);
}

static void b_plus(emxArray_real_T *b, const emxArray_real_T *res)
{
  emxArray_real_T *b_res;
  const double *res_data;
  double *b_data;
  double *b_res_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  res_data = res->data;
  b_data = b->data;
  emxInit_real_T(&b_res, 1);
  i = b_res->size[0];
  if (b->size[0] == 1) {
    b_res->size[0] = res->size[0];
  } else {
    b_res->size[0] = b->size[0];
  }

  emxEnsureCapacity_real_T(b_res, i);
  b_res_data = b_res->data;
  stride_0_0 = (res->size[0] != 1);
  stride_1_0 = (b->size[0] != 1);
  if (b->size[0] == 1) {
    loop_ub = res->size[0];
  } else {
    loop_ub = b->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_res_data[i] = res_data[i * stride_0_0] + b_data[i * stride_1_0];
  }

  i = b->size[0];
  b->size[0] = b_res->size[0];
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  loop_ub = b_res->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_data[i] = b_res_data[i];
  }

  emxFree_real_T(&b_res);
}

/*
 * function [p, mat,pwms, c] = seed_kmers(fn, num, pn, ik)
 */
static void b_seed_kmers(const emxArray_char_T *fn, double num, const
  emxArray_cell_wrap_1 *ik, emxArray_cell_wrap_0 *p, emxArray_cell_wrap_1 *mat,
  emxArray_cell_wrap_2 *pwms, double *c)
{
  cell_wrap_0 *p_data;
  cell_wrap_0 *s_data;
  cell_wrap_0 *sequences_data;
  const cell_wrap_1 *ik_data;
  cell_wrap_1 *b_mat_data;
  cell_wrap_1 *mat_data;
  cell_wrap_2 *b_pwms_data;
  cell_wrap_2 *pwms_data;
  emxArray_boolean_T *b_x;
  emxArray_cell_wrap_0 *s;
  emxArray_cell_wrap_0 *sequences;
  emxArray_cell_wrap_1 *b_mat;
  emxArray_cell_wrap_2 *b_pwms;
  emxArray_char_T *b_fileid;
  emxArray_char_T *cur_alpha;
  emxArray_char_T *cur_line;
  emxArray_char_T *cur_seq;
  emxArray_int32_T *D;
  emxArray_real_T *DD;
  emxArray_real_T *M;
  emxArray_real_T *alpha;
  emxArray_real_T *rs;
  emxArray_real_T *ss;
  creal_T dc;
  double h_y[9];
  double m[2];
  double curr_pos;
  double d;
  double idx;
  double *DD_data;
  double *M_data;
  double *alpha_data;
  double *rs_data;
  double *ss_data;
  int CC[9];
  int x_size[2];
  int b_d;
  int b_i;
  int b_loop_ub;
  int b_varargin_2;
  int b_y;
  int c_loop_ub;
  int c_y;
  int d_loop_ub;
  int d_y;
  int e_loop_ub;
  int e_y;
  int exitg1;
  int f_y;
  int g_y;
  int i;
  int i1;
  int i10;
  int i11;
  int i12;
  int i13;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int i8;
  int i9;
  int iindx;
  int j;
  int k;
  int l;
  int loop_ub;
  int loop_ub_tmp;
  int unnamed_idx_0_tmp_tmp;
  int varargin_2;
  int vlen;
  int x;
  int *D_data;
  signed char B[9];
  signed char BB[9];
  signed char fileid;
  bool exitg2;
  bool y;
  bool *x_data;
  ik_data = ik->data;

  /* 'gkmPWM:172' fid = fopen(fn, 'r'); */
  fileid = cfopen(fn, "rb");

  /* 'gkmPWM:173' if fid == -1 */
  if (fileid == -1) {
    /* 'gkmPWM:174' fprintf("ERROR: Weight file cannot be opened.\n") */
    printf("ERROR: Weight file cannot be opened.\n");
    fflush(stdout);
    exit(1);
  }

  /*  a = textscan(fid, '%s\t%f\n'); */
  /* 'gkmPWM:177' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'gkmPWM:178' idx=0; */
  idx = 0.0;

  /* 'gkmPWM:179' while ~feof(fid) */
  emxInit_char_T(&b_fileid, 2);
  do {
    exitg1 = 0;
    d = b_feof(fileid);
    if (d == 0.0) {
      /* 'gkmPWM:180' idx=idx+1; */
      idx++;

      /* 'gkmPWM:181' fgetl(fid); */
      b_fgets(fileid, b_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_fileid);
  emxInit_cell_wrap_0(&sequences, 1);

  /* 'gkmPWM:183' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'gkmPWM:184' sequences = cell(idx, 1); */
  unnamed_idx_0_tmp_tmp = (int)idx;
  i = sequences->size[0];
  sequences->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_0(sequences, i);
  sequences_data = sequences->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    sequences_data[i].f1->size[0] = 1;
    sequences_data[i].f1->size[1] = 0;
  }

  emxInit_real_T(&alpha, 1);

  /* 'gkmPWM:185' sequences = coder.nullcopy(sequences); */
  /* 'gkmPWM:186' alpha = zeros(idx, 1); */
  i = alpha->size[0];
  alpha->size[0] = (int)idx;
  emxEnsureCapacity_real_T(alpha, i);
  alpha_data = alpha->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    alpha_data[i] = 0.0;
  }

  /* 'gkmPWM:187' for cur_idx=1:idx */
  b_d = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&cur_seq, 2);
  emxInit_char_T(&cur_alpha, 2);
  exitg2 = false;
  while ((!exitg2) && (b_d <= (int)idx - 1)) {
    /* 'gkmPWM:188' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'gkmPWM:189' if cur_line == -1 */
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
      /* 'gkmPWM:192' [cur_seq, cur_alpha] = strtok(cur_line, char(9)); */
      b_strtok(cur_line, cur_seq, cur_alpha);

      /* 'gkmPWM:193' alpha(cur_idx,1) = real(str2double(cur_alpha)); */
      dc = str2double(cur_alpha);
      alpha_data[b_d] = dc.re;

      /* 'gkmPWM:194' sequences{cur_idx} = (strip(cur_seq)); */
      strip(cur_seq, sequences_data[b_d].f1);
      b_d++;
    }
  }

  emxFree_char_T(&cur_alpha);
  emxFree_char_T(&cur_seq);
  emxFree_char_T(&cur_line);
  emxInit_real_T(&M, 1);
  emxInit_int32_T(&D, 1);

  /* 'gkmPWM:196' fclose(fid); */
  cfclose(fileid);

  /*  [w, ind] = sort(a{2}, pn); */
  /*  s = a{1}(ind(1:min([100000 length(a{1})]))); */
  /* 'gkmPWM:201' [w, ind] = sort(alpha, pn); */
  d_sort(alpha, D);
  D_data = D->data;
  alpha_data = alpha->data;
  i = M->size[0];
  M->size[0] = D->size[0];
  emxEnsureCapacity_real_T(M, i);
  M_data = M->data;
  loop_ub = D->size[0];
  for (i = 0; i < loop_ub; i++) {
    M_data[i] = D_data[i];
  }

  emxInit_cell_wrap_0(&s, 1);

  /* 'gkmPWM:202' s_len = min([100000 length(sequences)]); */
  m[0] = 100000.0;
  m[1] = sequences->size[0];

  /* 'gkmPWM:203' s = cell(s_len, 1); */
  unnamed_idx_0_tmp_tmp = (int)b_minimum(m);
  i = s->size[0];
  s->size[0] = unnamed_idx_0_tmp_tmp;
  emxEnsureCapacity_cell_wrap_0(s, i);
  s_data = s->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    s_data[i].f1->size[0] = 1;
    s_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:204' s = coder.nullcopy(s); */
  /* 'gkmPWM:205' for cur_idx=1:s_len */
  for (b_d = 0; b_d < unnamed_idx_0_tmp_tmp; b_d++) {
    /* 'gkmPWM:206' s{cur_idx} = sequences{ind(cur_idx)}; */
    i = s_data[b_d].f1->size[0] * s_data[b_d].f1->size[1];
    s_data[b_d].f1->size[0] = 1;
    s_data[b_d].f1->size[1] = sequences_data[(int)M_data[b_d] - 1].f1->size[1];
    emxEnsureCapacity_char_T(s_data[b_d].f1, i);
    loop_ub = sequences_data[(int)M_data[b_d] - 1].f1->size[1];
    for (i = 0; i < loop_ub; i++) {
      s_data[b_d].f1->data[i] = sequences_data[(int)M_data[b_d] - 1].f1->data[i];
    }
  }

  /* 'gkmPWM:210' l = length(s{1}); */
  varargin_2 = s_data[0].f1->size[1] - 1;
  l = s_data[0].f1->size[1] - 4;

  /* 'gkmPWM:211' k = round(l/2)+1; */
  x = (int)rt_roundd((double)s_data[0].f1->size[1] / 2.0);

  /* 'gkmPWM:212' ikl = length(ik); */
  b_varargin_2 = ik->size[0];

  /* 'gkmPWM:213' p = cell(num,1); */
  unnamed_idx_0_tmp_tmp = (int)num;
  i = p->size[0];
  p->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_0(p, i);
  p_data = p->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    p_data[i].f1->size[0] = 1;
    p_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:214' p = coder.nullcopy(p); */
  i = sequences->size[0];
  sequences->size[0] = p->size[0];
  emxEnsureCapacity_cell_wrap_0(sequences, i);
  sequences_data = sequences->data;

  /* 'gkmPWM:215' c = ikl+1; */
  *c = (double)ik->size[0] + 1.0;

  /* 'gkmPWM:216' p{1} = s{1}; */
  i = sequences_data[0].f1->size[0] * sequences_data[0].f1->size[1];
  sequences_data[0].f1->size[0] = 1;
  sequences_data[0].f1->size[1] = s_data[0].f1->size[1];
  emxEnsureCapacity_char_T(sequences_data[0].f1, i);
  loop_ub = s_data[0].f1->size[1];
  for (i = 0; i < loop_ub; i++) {
    sequences_data[0].f1->data[i] = s_data[0].f1->data[i];
  }

  emxInit_cell_wrap_1(&b_mat);

  /* 'gkmPWM:217' mat = cell(ikl+num,1); */
  b_d = (int)((double)ik->size[0] + num);
  i = b_mat->size[0];
  b_mat->size[0] = (int)((double)ik->size[0] + num);
  emxEnsureCapacity_cell_wrap_1(b_mat, i);
  mat_data = b_mat->data;
  for (i = 0; i < b_d; i++) {
    mat_data[i].f1->size[0] = 1;
    mat_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:218' mat = coder.nullcopy(mat); */
  /*  mat(1:ikl) = ik; */
  /* 'gkmPWM:220' for cur_idx=1:length(ik) */
  i = ik->size[0];
  for (b_d = 0; b_d < i; b_d++) {
    /* 'gkmPWM:221' mat{cur_idx} = ik{cur_idx}; */
    b_y = mat_data[b_d].f1->size[0] * mat_data[b_d].f1->size[1];
    mat_data[b_d].f1->size[0] = 1;
    mat_data[b_d].f1->size[1] = ik_data[b_d].f1->size[1];
    emxEnsureCapacity_real_T(mat_data[b_d].f1, b_y);
    loop_ub = ik_data[b_d].f1->size[1];
    for (b_y = 0; b_y < loop_ub; b_y++) {
      mat_data[b_d].f1->data[b_y] = ik_data[b_d].f1->data[b_y];
    }
  }

  emxInit_cell_wrap_2(&b_pwms);

  /* 'gkmPWM:223' mat{c} = letterconvert(s{1}); */
  letterconvert(s_data[0].f1, mat_data[ik->size[0]].f1);

  /* 'gkmPWM:224' pwms = cell(num,1); */
  i = b_pwms->size[0];
  b_pwms->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_2(b_pwms, i);
  pwms_data = b_pwms->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    pwms_data[i].f1->size[0] = 0;
    pwms_data[i].f1->size[1] = 4;
  }

  /* 'gkmPWM:225' pwms = coder.nullcopy(pwms); */
  /* 'gkmPWM:226' for i = 1:num */
  for (b_i = 0; b_i < unnamed_idx_0_tmp_tmp; b_i++) {
    /* 'gkmPWM:227' pwms{i} = zeros(l,4); */
    i = pwms_data[b_i].f1->size[0] * pwms_data[b_i].f1->size[1];
    pwms_data[b_i].f1->size[0] = varargin_2 + 1;
    pwms_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(pwms_data[b_i].f1, i);
    loop_ub = (varargin_2 + 1) << 2;
    for (i = 0; i < loop_ub; i++) {
      pwms_data[b_i].f1->data[i] = 0.0;
    }
  }

  /* 'gkmPWM:229' for i = 1:l */
  i = s_data[0].f1->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:230' pwms{1}(i,mat{c}(i)+1) = pwms{1}(i,mat{c}(i)+1)+w(i); */
    d = mat_data[b_varargin_2].f1->data[b_i];
    pwms_data[0].f1->data[b_i + pwms_data[0].f1->size[0] * ((int)(d + 1.0) - 1)]
      += alpha_data[b_i];
  }

  /* 'gkmPWM:232' B = zeros(9,1); */
  /* 'gkmPWM:233' BB = zeros(9,1); */
  for (b_i = 0; b_i < 9; b_i++) {
    B[b_i] = 0;
    BB[b_i] = 0;
  }

  /* 'gkmPWM:234' B(1:5) = (0:4)'; */
  for (i = 0; i < 5; i++) {
    B[i] = (signed char)i;
  }

  /* 'gkmPWM:235' B(6:9) = 0; */
  for (b_i = 0; b_i < 4; b_i++) {
    B[b_i + 5] = 0;
  }

  /* 'gkmPWM:236' BB(1:5) = 0; */
  for (b_i = 0; b_i < 5; b_i++) {
    BB[b_i] = 0;
  }

  /* 'gkmPWM:237' BB(6:9) = (1:4)'; */
  for (i = 0; i < 4; i++) {
    BB[i + 5] = (signed char)(i + 1);
  }

  /* 'gkmPWM:238' CC = [l l-1 l-2 l-3 l-4 l-1 l-2 l-3 l-4]; */
  CC[0] = s_data[0].f1->size[1];
  CC[1] = s_data[0].f1->size[1] - 1;
  CC[2] = s_data[0].f1->size[1] - 2;
  CC[3] = s_data[0].f1->size[1] - 3;
  CC[4] = s_data[0].f1->size[1] - 4;
  CC[5] = s_data[0].f1->size[1] - 1;
  CC[6] = s_data[0].f1->size[1] - 2;
  CC[7] = s_data[0].f1->size[1] - 3;
  CC[8] = s_data[0].f1->size[1] - 4;

  /* this process picks kmers to seed the PWMs.  kmers that match one of the seeds by round(l/2)+1 or more are added that particular seed.  Otherwise, it becomes another seed. */
  /* 'gkmPWM:240' for i = 2:100000 */
  b_i = 1;
  emxInit_real_T(&ss, 2);
  emxInit_real_T(&rs, 2);
  emxInit_real_T(&DD, 1);
  emxInit_boolean_T(&b_x, 2);
  exitg2 = false;
  while ((!exitg2) && (b_i - 1 < 99999)) {
    /* 'gkmPWM:241' ss = letterconvert(s{i}); */
    letterconvert(s_data[b_i].f1, ss);
    ss_data = ss->data;

    /* 'gkmPWM:242' rs = 3-fliplr(ss); */
    i = rs->size[0] * rs->size[1];
    rs->size[0] = 1;
    rs->size[1] = ss->size[1];
    emxEnsureCapacity_real_T(rs, i);
    rs_data = rs->data;
    loop_ub = ss->size[1];
    for (i = 0; i < loop_ub; i++) {
      rs_data[i] = ss_data[i];
    }

    fliplr(rs);
    i = rs->size[0] * rs->size[1];
    rs->size[0] = 1;
    emxEnsureCapacity_real_T(rs, i);
    rs_data = rs->data;
    loop_ub = rs->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      rs_data[i] = 3.0 - rs_data[i];
    }

    /* 'gkmPWM:243' M = zeros(c,1); */
    loop_ub_tmp = (int)*c;
    i = M->size[0];
    M->size[0] = (int)*c;
    emxEnsureCapacity_real_T(M, i);
    M_data = M->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      M_data[i] = 0.0;
    }

    /* 'gkmPWM:244' D = zeros(c,1); */
    i = D->size[0];
    D->size[0] = (int)*c;
    emxEnsureCapacity_int32_T(D, i);
    D_data = D->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      D_data[i] = 0;
    }

    /* 'gkmPWM:245' DD = zeros(c,1); */
    i = DD->size[0];
    DD->size[0] = (int)*c;
    emxEnsureCapacity_real_T(DD, i);
    DD_data = DD->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      DD_data[i] = 0.0;
    }

    /* 'gkmPWM:246' for j = 1:c */
    for (j = 0; j < loop_ub_tmp; j++) {
      /* 'gkmPWM:247' [m,d] = max([sum(mat{j}==ss) sum(mat{j}(2:end)==ss(1:l-1)) sum(mat{j}(3:end)==ss(1:l-2)) sum(mat{j}(4:end)==ss(1:l-3)) sum(mat{j}(5:end)==ss(1:l-4)) sum(mat{j}(1:l-1)==ss(2:end)) sum(mat{j}(1:l-2)==ss(3:end)) sum(mat{j}(1:l-3)==ss(4:end)) sum(mat{j}(1:l-4)==ss(5:end))]); */
      if (2 > mat_data[j].f1->size[1]) {
        i = 0;
        b_y = 0;
      } else {
        i = 1;
        b_y = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2) {
        unnamed_idx_0_tmp_tmp = -1;
      } else {
        unnamed_idx_0_tmp_tmp = l + 2;
      }

      if (3 > mat_data[j].f1->size[1]) {
        b_d = 0;
        c_y = 0;
      } else {
        b_d = 2;
        c_y = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 1) {
        d_y = -1;
      } else {
        d_y = l + 1;
      }

      if (4 > mat_data[j].f1->size[1]) {
        e_y = 0;
        i1 = 0;
      } else {
        e_y = 3;
        i1 = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 2) {
        i2 = -1;
      } else {
        i2 = l;
      }

      if (5 > mat_data[j].f1->size[1]) {
        i3 = 0;
        i4 = 0;
      } else {
        i3 = 4;
        i4 = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 3) {
        i5 = 0;
      } else {
        i5 = l;
      }

      if (1 > varargin_2) {
        loop_ub = -1;
      } else {
        loop_ub = l + 2;
      }

      if (2 > ss->size[1]) {
        i6 = 0;
        i7 = 0;
      } else {
        i6 = 1;
        i7 = ss->size[1];
      }

      if (1 > varargin_2 - 1) {
        b_loop_ub = -1;
      } else {
        b_loop_ub = l + 1;
      }

      if (3 > ss->size[1]) {
        i8 = 0;
        i9 = 0;
      } else {
        i8 = 2;
        i9 = ss->size[1];
      }

      if (1 > varargin_2 - 2) {
        c_loop_ub = -1;
      } else {
        c_loop_ub = l;
      }

      if (4 > ss->size[1]) {
        i10 = 0;
        i11 = 0;
      } else {
        i10 = 3;
        i11 = ss->size[1];
      }

      if (1 > varargin_2 - 3) {
        d_loop_ub = 0;
      } else {
        d_loop_ub = l;
      }

      if (5 > ss->size[1]) {
        i12 = 0;
        i13 = 0;
      } else {
        i12 = 4;
        i13 = ss->size[1];
      }

      if (mat_data[j].f1->size[1] == ss->size[1]) {
        f_y = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = mat_data[j].f1->size[1];
        emxEnsureCapacity_boolean_T(b_x, f_y);
        x_data = b_x->data;
        e_loop_ub = mat_data[j].f1->size[1];
        for (f_y = 0; f_y < e_loop_ub; f_y++) {
          x_data[f_y] = (mat_data[j].f1->data[f_y] == ss_data[f_y]);
        }
      } else {
        l_binary_expand_op(b_x, b_mat, j, ss);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        g_y = 0;
      } else {
        g_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          g_y += x_data[k - 1];
        }
      }

      e_loop_ub = b_y - i;
      if (e_loop_ub == unnamed_idx_0_tmp_tmp + 1) {
        b_y = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, b_y);
        x_data = b_x->data;
        for (b_y = 0; b_y < e_loop_ub; b_y++) {
          x_data[b_y] = (mat_data[j].f1->data[i + b_y] == ss_data[b_y]);
        }
      } else {
        k_binary_expand_op(b_x, b_mat, j, i, b_y - 1, ss, unnamed_idx_0_tmp_tmp
                           - 2);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        b_y = 0;
      } else {
        b_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          b_y += x_data[k - 1];
        }
      }

      e_loop_ub = c_y - b_d;
      if (e_loop_ub == d_y + 1) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[b_d + i] == ss_data[i]);
        }
      } else {
        j_binary_expand_op(b_x, b_mat, j, b_d, c_y - 1, ss, d_y - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        f_y = 0;
      } else {
        f_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          f_y += x_data[k - 1];
        }
      }

      e_loop_ub = i1 - e_y;
      if (e_loop_ub == i2 + 1) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[e_y + i] == ss_data[i]);
        }
      } else {
        i_binary_expand_op(b_x, b_mat, j, e_y, i1 - 1, ss, i2);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        e_y = 0;
      } else {
        e_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          e_y += x_data[k - 1];
        }
      }

      e_loop_ub = i4 - i3;
      if (e_loop_ub == i5) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i3 + i] == ss_data[i]);
        }
      } else {
        i_binary_expand_op(b_x, b_mat, j, i3, i4 - 1, ss, i5 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        e_loop_ub = 0;
      } else {
        e_loop_ub = x_data[0];
        for (k = 2; k <= vlen; k++) {
          e_loop_ub += x_data[k - 1];
        }
      }

      if (loop_ub + 1 == i7 - i6) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == ss_data[i6 + i]);
        }
      } else {
        h_binary_expand_op(b_x, b_mat, j, loop_ub - 2, ss, i6, i7 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        d_y = 0;
      } else {
        d_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          d_y += x_data[k - 1];
        }
      }

      if (b_loop_ub + 1 == i9 - i8) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = b_loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= b_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == ss_data[i8 + i]);
        }
      } else {
        g_binary_expand_op(b_x, b_mat, j, b_loop_ub - 1, ss, i8, i9 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        c_y = 0;
      } else {
        c_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          c_y += x_data[k - 1];
        }
      }

      if (c_loop_ub + 1 == i11 - i10) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = c_loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= c_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == ss_data[i10 + i]);
        }
      } else {
        f_binary_expand_op(b_x, b_mat, j, c_loop_ub, ss, i10, i11 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        unnamed_idx_0_tmp_tmp = 0;
      } else {
        unnamed_idx_0_tmp_tmp = x_data[0];
        for (k = 2; k <= vlen; k++) {
          unnamed_idx_0_tmp_tmp += x_data[k - 1];
        }
      }

      if (d_loop_ub == i13 - i12) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = d_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < d_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == ss_data[i12 + i]);
        }
      } else {
        f_binary_expand_op(b_x, b_mat, j, d_loop_ub - 1, ss, i12, i13 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        b_d = 0;
      } else {
        b_d = x_data[0];
        for (k = 2; k <= vlen; k++) {
          b_d += x_data[k - 1];
        }
      }

      h_y[0] = g_y;
      h_y[1] = b_y;
      h_y[2] = f_y;
      h_y[3] = e_y;
      h_y[4] = e_loop_ub;
      h_y[5] = d_y;
      h_y[6] = c_y;
      h_y[7] = unnamed_idx_0_tmp_tmp;
      h_y[8] = b_d;
      b_maximum(h_y, &curr_pos, &iindx);

      /* 'gkmPWM:248' [mm,dd] = max([sum(mat{j}==rs) sum(mat{j}(2:end)==rs(1:l-1)) sum(mat{j}(3:end)==rs(1:l-2)) sum(mat{j}(4:end)==rs(1:l-3)) sum(mat{j}(5:end)==rs(1:l-4)) sum(mat{j}(1:l-1)==rs(2:end)) sum(mat{j}(1:l-2)==rs(3:end)) sum(mat{j}(1:l-3)==rs(4:end)) sum(mat{j}(1:l-4)==rs(5:end))]); */
      if (2 > mat_data[j].f1->size[1]) {
        i = 0;
        b_y = 0;
      } else {
        i = 1;
        b_y = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2) {
        unnamed_idx_0_tmp_tmp = -1;
      } else {
        unnamed_idx_0_tmp_tmp = l + 2;
      }

      if (3 > mat_data[j].f1->size[1]) {
        b_d = 0;
        c_y = 0;
      } else {
        b_d = 2;
        c_y = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 1) {
        d_y = -1;
      } else {
        d_y = l + 1;
      }

      if (4 > mat_data[j].f1->size[1]) {
        e_y = 0;
        i1 = 0;
      } else {
        e_y = 3;
        i1 = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 2) {
        i2 = -1;
      } else {
        i2 = l;
      }

      if (5 > mat_data[j].f1->size[1]) {
        i3 = 0;
        i4 = 0;
      } else {
        i3 = 4;
        i4 = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 3) {
        i5 = 0;
      } else {
        i5 = l;
      }

      if (1 > varargin_2) {
        loop_ub = -1;
      } else {
        loop_ub = l + 2;
      }

      if (2 > rs->size[1]) {
        i6 = 0;
        i7 = 0;
      } else {
        i6 = 1;
        i7 = rs->size[1];
      }

      if (1 > varargin_2 - 1) {
        b_loop_ub = -1;
      } else {
        b_loop_ub = l + 1;
      }

      if (3 > rs->size[1]) {
        i8 = 0;
        i9 = 0;
      } else {
        i8 = 2;
        i9 = rs->size[1];
      }

      if (1 > varargin_2 - 2) {
        c_loop_ub = -1;
      } else {
        c_loop_ub = l;
      }

      if (4 > rs->size[1]) {
        i10 = 0;
        i11 = 0;
      } else {
        i10 = 3;
        i11 = rs->size[1];
      }

      if (1 > varargin_2 - 3) {
        d_loop_ub = 0;
      } else {
        d_loop_ub = l;
      }

      if (5 > rs->size[1]) {
        i12 = 0;
        i13 = 0;
      } else {
        i12 = 4;
        i13 = rs->size[1];
      }

      if (mat_data[j].f1->size[1] == rs->size[1]) {
        f_y = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = mat_data[j].f1->size[1];
        emxEnsureCapacity_boolean_T(b_x, f_y);
        x_data = b_x->data;
        e_loop_ub = mat_data[j].f1->size[1];
        for (f_y = 0; f_y < e_loop_ub; f_y++) {
          x_data[f_y] = (mat_data[j].f1->data[f_y] == rs_data[f_y]);
        }
      } else {
        l_binary_expand_op(b_x, b_mat, j, rs);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        g_y = 0;
      } else {
        g_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          g_y += x_data[k - 1];
        }
      }

      e_loop_ub = b_y - i;
      if (e_loop_ub == unnamed_idx_0_tmp_tmp + 1) {
        b_y = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, b_y);
        x_data = b_x->data;
        for (b_y = 0; b_y < e_loop_ub; b_y++) {
          x_data[b_y] = (mat_data[j].f1->data[i + b_y] == rs_data[b_y]);
        }
      } else {
        k_binary_expand_op(b_x, b_mat, j, i, b_y - 1, rs, unnamed_idx_0_tmp_tmp
                           - 2);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        b_y = 0;
      } else {
        b_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          b_y += x_data[k - 1];
        }
      }

      e_loop_ub = c_y - b_d;
      if (e_loop_ub == d_y + 1) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[b_d + i] == rs_data[i]);
        }
      } else {
        j_binary_expand_op(b_x, b_mat, j, b_d, c_y - 1, rs, d_y - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        f_y = 0;
      } else {
        f_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          f_y += x_data[k - 1];
        }
      }

      e_loop_ub = i1 - e_y;
      if (e_loop_ub == i2 + 1) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[e_y + i] == rs_data[i]);
        }
      } else {
        i_binary_expand_op(b_x, b_mat, j, e_y, i1 - 1, rs, i2);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        e_y = 0;
      } else {
        e_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          e_y += x_data[k - 1];
        }
      }

      e_loop_ub = i4 - i3;
      if (e_loop_ub == i5) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i3 + i] == rs_data[i]);
        }
      } else {
        i_binary_expand_op(b_x, b_mat, j, i3, i4 - 1, rs, i5 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        e_loop_ub = 0;
      } else {
        e_loop_ub = x_data[0];
        for (k = 2; k <= vlen; k++) {
          e_loop_ub += x_data[k - 1];
        }
      }

      if (loop_ub + 1 == i7 - i6) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == rs_data[i6 + i]);
        }
      } else {
        h_binary_expand_op(b_x, b_mat, j, loop_ub - 2, rs, i6, i7 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        d_y = 0;
      } else {
        d_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          d_y += x_data[k - 1];
        }
      }

      if (b_loop_ub + 1 == i9 - i8) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = b_loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= b_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == rs_data[i8 + i]);
        }
      } else {
        g_binary_expand_op(b_x, b_mat, j, b_loop_ub - 1, rs, i8, i9 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        c_y = 0;
      } else {
        c_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          c_y += x_data[k - 1];
        }
      }

      if (c_loop_ub + 1 == i11 - i10) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = c_loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= c_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == rs_data[i10 + i]);
        }
      } else {
        f_binary_expand_op(b_x, b_mat, j, c_loop_ub, rs, i10, i11 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        unnamed_idx_0_tmp_tmp = 0;
      } else {
        unnamed_idx_0_tmp_tmp = x_data[0];
        for (k = 2; k <= vlen; k++) {
          unnamed_idx_0_tmp_tmp += x_data[k - 1];
        }
      }

      if (d_loop_ub == i13 - i12) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = d_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < d_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == rs_data[i12 + i]);
        }
      } else {
        f_binary_expand_op(b_x, b_mat, j, d_loop_ub - 1, rs, i12, i13 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        b_d = 0;
      } else {
        b_d = x_data[0];
        for (k = 2; k <= vlen; k++) {
          b_d += x_data[k - 1];
        }
      }

      h_y[0] = g_y;
      h_y[1] = b_y;
      h_y[2] = f_y;
      h_y[3] = e_y;
      h_y[4] = e_loop_ub;
      h_y[5] = d_y;
      h_y[6] = c_y;
      h_y[7] = unnamed_idx_0_tmp_tmp;
      h_y[8] = b_d;
      b_maximum(h_y, &idx, &b_d);

      /* 'gkmPWM:249' [M(j),ddd] = max([m mm]); */
      m[0] = curr_pos;
      m[1] = idx;
      c_maximum(m, &M_data[j], &unnamed_idx_0_tmp_tmp);

      /* 'gkmPWM:250' if ddd == 1 */
      if (unnamed_idx_0_tmp_tmp == 1) {
        /* 'gkmPWM:251' D(j) = d; */
        D_data[j] = iindx;

        /* 'gkmPWM:252' DD(j) = 1; */
        DD_data[j] = 1.0;
      } else {
        /* 'gkmPWM:253' else */
        /* 'gkmPWM:254' D(j) = dd; */
        D_data[j] = b_d;

        /* 'gkmPWM:255' DD(j) = 2; */
        DD_data[j] = 2.0;
      }
    }

    /* 'gkmPWM:258' if max(M) < k */
    if (maximum(M) < (double)x + 1.0) {
      /* 'gkmPWM:259' c = c+1; */
      (*c)++;

      /* 'gkmPWM:260' p{c-ikl} = s{i}; */
      i = (int)(*c - (double)b_varargin_2) - 1;
      i = sequences_data[i].f1->size[0] * sequences_data[i].f1->size[1];
      sequences_data[(int)(*c - (double)b_varargin_2) - 1].f1->size[0] = 1;
      sequences_data[(int)(*c - (double)b_varargin_2) - 1].f1->size[1] =
        s_data[b_i].f1->size[1];
      emxEnsureCapacity_char_T(sequences_data[(int)(*c - (double)b_varargin_2) -
        1].f1, i);
      loop_ub = s_data[b_i].f1->size[1];
      for (i = 0; i < loop_ub; i++) {
        sequences_data[(int)(*c - (double)b_varargin_2) - 1].f1->data[i] =
          s_data[b_i].f1->data[i];
      }

      /* 'gkmPWM:261' mat{c} = ss; */
      i = mat_data[(int)*c - 1].f1->size[0] * mat_data[(int)*c - 1].f1->size[1];
      mat_data[(int)*c - 1].f1->size[0] = 1;
      mat_data[(int)*c - 1].f1->size[1] = ss->size[1];
      emxEnsureCapacity_real_T(mat_data[(int)*c - 1].f1, i);
      loop_ub = ss->size[1];
      for (i = 0; i < loop_ub; i++) {
        mat_data[(int)*c - 1].f1->data[i] = ss_data[i];
      }

      /* 'gkmPWM:262' ss = ss+1; */
      i = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity_real_T(ss, i);
      ss_data = ss->data;
      loop_ub = ss->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        ss_data[i]++;
      }

      /* 'gkmPWM:263' for j = 1:l */
      for (j = 0; j <= varargin_2; j++) {
        /* 'gkmPWM:264' pwms{c-ikl}(j,ss(j)) = pwms{c-ikl}(j,ss(j))+w(i); */
        i = (int)ss_data[j] - 1;
        pwms_data[(int)(*c - (double)b_varargin_2) - 1].f1->data[j + pwms_data
          [(int)(*c - (double)b_varargin_2) - 1].f1->size[0] * i] = pwms_data
          [(int)(*c - (double)ik->size[0]) - 1].f1->data[j + pwms_data[(int)(*c
          - (double)ik->size[0]) - 1].f1->size[0] * i] + alpha_data[b_i];
      }
    } else {
      /* 'gkmPWM:266' else */
      /* 'gkmPWM:267' [~,d] = max(M); */
      d_maximum(M, &curr_pos, &iindx);

      /* 'gkmPWM:268' if DD(d) == 1 && d > ikl */
      d = DD_data[iindx - 1];
      if ((d == 1.0) && (iindx > b_varargin_2)) {
        /* 'gkmPWM:269' ss = ss+1; */
        i = ss->size[0] * ss->size[1];
        ss->size[0] = 1;
        emxEnsureCapacity_real_T(ss, i);
        ss_data = ss->data;
        loop_ub = ss->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          ss_data[i]++;
        }

        /* 'gkmPWM:270' d = d-ikl; */
        b_d = (iindx - b_varargin_2) - 1;

        /* 'gkmPWM:271' for j = 1:CC(D(d)) */
        i = CC[D_data[b_d] - 1];
        for (j = 0; j < i; j++) {
          /* 'gkmPWM:272' pwms{d}(j+B(D(d)),ss(j+BB(D(d)))) = pwms{d}(j+B(D(d)),ss(j+BB(D(d))))+w(i); */
          b_y = (int)((unsigned int)j + B[D_data[b_d] - 1]);
          unnamed_idx_0_tmp_tmp = (int)ss_data[(int)((unsigned int)j +
            BB[D_data[b_d] - 1])] - 1;
          pwms_data[b_d].f1->data[b_y + pwms_data[b_d].f1->size[0] *
            unnamed_idx_0_tmp_tmp] += alpha_data[b_i];
        }
      } else if ((d == 2.0) && (iindx > b_varargin_2)) {
        /* 'gkmPWM:274' elseif DD(d) == 2 && d > ikl */
        /* 'gkmPWM:275' rs = rs+1; */
        i = rs->size[0] * rs->size[1];
        rs->size[0] = 1;
        emxEnsureCapacity_real_T(rs, i);
        rs_data = rs->data;
        loop_ub = rs->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          rs_data[i]++;
        }

        /* 'gkmPWM:276' d = d-ikl; */
        b_d = (iindx - b_varargin_2) - 1;

        /* 'gkmPWM:277' for j = 1:CC(D(d)) */
        i = CC[D_data[b_d] - 1];
        for (j = 0; j < i; j++) {
          /* 'gkmPWM:278' pwms{d}(j+B(D(d)),rs(j+BB(D(d)))) = pwms{d}(j+B(D(d)),rs(j+BB(D(d))))+w(i); */
          b_y = (int)((unsigned int)j + B[D_data[b_d] - 1]);
          unnamed_idx_0_tmp_tmp = (int)rs_data[(int)((unsigned int)j +
            BB[D_data[b_d] - 1])] - 1;
          pwms_data[b_d].f1->data[b_y + pwms_data[b_d].f1->size[0] *
            unnamed_idx_0_tmp_tmp] += alpha_data[b_i];
        }
      }
    }

    /* 'gkmPWM:282' if c == num+ikl */
    if (*c == num + (double)b_varargin_2) {
      exitg2 = true;
    } else {
      b_i++;
    }
  }

  emxFree_boolean_T(&b_x);
  emxFree_real_T(&DD);
  emxFree_int32_T(&D);
  emxFree_real_T(&M);
  emxFree_real_T(&rs);
  emxFree_real_T(&ss);
  emxFree_cell_wrap_0(&s);
  emxFree_real_T(&alpha);

  /*  mat = mat(1:c); */
  /* 'gkmPWM:287' new_mat = cell(c, 1); */
  /* 'gkmPWM:288' for cur_idx=1:c */
  i = (int)*c;
  b_y = mat->size[0];
  mat->size[0] = (int)*c;
  emxEnsureCapacity_cell_wrap_1(mat, b_y);
  b_mat_data = mat->data;
  for (b_d = 0; b_d < i; b_d++) {
    /* 'gkmPWM:289' new_mat{cur_idx} = mat{cur_idx}; */
    b_y = b_mat_data[b_d].f1->size[0] * b_mat_data[b_d].f1->size[1];
    b_mat_data[b_d].f1->size[0] = 1;
    b_mat_data[b_d].f1->size[1] = mat_data[b_d].f1->size[1];
    emxEnsureCapacity_real_T(b_mat_data[b_d].f1, b_y);
    loop_ub = mat_data[b_d].f1->size[1];
    for (b_y = 0; b_y < loop_ub; b_y++) {
      b_mat_data[b_d].f1->data[b_y] = mat_data[b_d].f1->data[b_y];
    }
  }

  emxFree_cell_wrap_1(&b_mat);

  /* 'gkmPWM:291' mat = new_mat; */
  /*  p = p(1:c-ikl); */
  /* 'gkmPWM:294' p_len = c-ikl; */
  /* 'gkmPWM:295' new_p = cell(p_len, 1); */
  /* 'gkmPWM:296' for cur_idx=1:p_len */
  i = (int)(*c - (double)ik->size[0]);
  b_y = p->size[0];
  p->size[0] = (int)(*c - (double)ik->size[0]);
  emxEnsureCapacity_cell_wrap_0(p, b_y);
  p_data = p->data;
  for (b_d = 0; b_d < i; b_d++) {
    /* 'gkmPWM:297' new_p{cur_idx} = p{cur_idx}; */
    b_y = p_data[b_d].f1->size[0] * p_data[b_d].f1->size[1];
    p_data[b_d].f1->size[0] = 1;
    p_data[b_d].f1->size[1] = sequences_data[b_d].f1->size[1];
    emxEnsureCapacity_char_T(p_data[b_d].f1, b_y);
    loop_ub = sequences_data[b_d].f1->size[1];
    for (b_y = 0; b_y < loop_ub; b_y++) {
      p_data[b_d].f1->data[b_y] = sequences_data[b_d].f1->data[b_y];
    }
  }

  emxFree_cell_wrap_0(&sequences);

  /* 'gkmPWM:299' p = new_p; */
  /*  pwms = pwms(1:c-ikl); */
  /* 'gkmPWM:302' pwms_len = c-ikl; */
  /* 'gkmPWM:303' new_pwms = cell(pwms_len, 1); */
  /* 'gkmPWM:304' for cur_idx=1:pwms_len */
  i = (int)(*c - (double)ik->size[0]);
  b_y = pwms->size[0];
  pwms->size[0] = (int)(*c - (double)ik->size[0]);
  emxEnsureCapacity_cell_wrap_2(pwms, b_y);
  b_pwms_data = pwms->data;
  for (b_d = 0; b_d < i; b_d++) {
    /* 'gkmPWM:305' new_pwms{cur_idx} = pwms{cur_idx}; */
    b_y = b_pwms_data[b_d].f1->size[0] * b_pwms_data[b_d].f1->size[1];
    b_pwms_data[b_d].f1->size[0] = pwms_data[b_d].f1->size[0];
    b_pwms_data[b_d].f1->size[1] = 4;
    emxEnsureCapacity_real_T(b_pwms_data[b_d].f1, b_y);
    loop_ub = pwms_data[b_d].f1->size[0] * 4;
    for (b_y = 0; b_y < loop_ub; b_y++) {
      b_pwms_data[b_d].f1->data[b_y] = pwms_data[b_d].f1->data[b_y];
    }
  }

  emxFree_cell_wrap_2(&b_pwms);

  /* 'gkmPWM:307' pwms = new_pwms; */
  /* 'gkmPWM:309' for i = 1:c-ikl */
  i = (int)(*c - (double)ik->size[0]);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:310' for j = 1:l */
    for (j = 0; j <= varargin_2; j++) {
      /* 'gkmPWM:311' pwms{i}(j,:) = pwms{i}(j,:)/sum(pwms{i}(j,:)); */
      curr_pos = b_pwms_data[b_i].f1->data[j];
      for (k = 0; k < 3; k++) {
        curr_pos += b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0] *
          (k + 1)];
      }

      for (b_y = 0; b_y < 4; b_y++) {
        b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0] * b_y] /=
          curr_pos;
      }
    }
  }
}

static void bb_binary_expand_op(emxArray_real_T *kweig, const emxArray_real_T
  *indvec2, int iii, const emxArray_cell_wrap_14 *ktree2, double k)
{
  const cell_wrap_14 *ktree2_data;
  emxArray_real_T *b_kweig;
  const double *indvec2_data;
  double *b_kweig_data;
  double *kweig_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  ktree2_data = ktree2->data;
  indvec2_data = indvec2->data;
  kweig_data = kweig->data;
  emxInit_real_T(&b_kweig, 1);
  i = b_kweig->size[0];
  if (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] == 1) {
    b_kweig->size[0] = indvec2->size[0];
  } else {
    b_kweig->size[0] = ktree2_data[(int)(k - 1.0) - 1].f1->size[0];
  }

  emxEnsureCapacity_real_T(b_kweig, i);
  b_kweig_data = b_kweig->data;
  stride_0_0 = (indvec2->size[0] != 1);
  stride_1_0 = (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] != 1);
  if (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] == 1) {
    loop_ub = indvec2->size[0];
  } else {
    loop_ub = ktree2_data[(int)(k - 1.0) - 1].f1->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_kweig_data[i] = kweig_data[((int)indvec2_data[i * stride_0_0] +
      kweig->size[0] * iii) - 1] + ktree2_data[(int)(k - 1.0) - 1].f1->data[i *
      stride_1_0];
  }

  loop_ub = b_kweig->size[0];
  for (i = 0; i < loop_ub; i++) {
    kweig_data[((int)indvec2_data[i] + kweig->size[0] * iii) - 1] =
      b_kweig_data[i];
  }

  emxFree_real_T(&b_kweig);
}

static void binary_expand_op(emxArray_real_T *cfile, const emxArray_real_T
  *negvec, double y, double b)
{
  emxArray_real_T *b_cfile;
  const double *negvec_data;
  double *b_cfile_data;
  double *cfile_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  negvec_data = negvec->data;
  cfile_data = cfile->data;
  emxInit_real_T(&b_cfile, 1);
  i = b_cfile->size[0];
  if (negvec->size[0] == 1) {
    b_cfile->size[0] = cfile->size[0];
  } else {
    b_cfile->size[0] = negvec->size[0];
  }

  emxEnsureCapacity_real_T(b_cfile, i);
  b_cfile_data = b_cfile->data;
  stride_0_0 = (cfile->size[0] != 1);
  stride_1_0 = (negvec->size[0] != 1);
  if (negvec->size[0] == 1) {
    loop_ub = cfile->size[0];
  } else {
    loop_ub = negvec->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_cfile_data[i] = cfile_data[i * stride_0_0] - negvec_data[i * stride_1_0] /
      y * b;
  }

  i = cfile->size[0];
  cfile->size[0] = b_cfile->size[0];
  emxEnsureCapacity_real_T(cfile, i);
  cfile_data = cfile->data;
  loop_ub = b_cfile->size[0];
  for (i = 0; i < loop_ub; i++) {
    cfile_data[i] = b_cfile_data[i];
  }

  emxFree_real_T(&b_cfile);
}

static void cb_binary_expand_op(emxArray_real_T *mat, const emxArray_real_T *x)
{
  emxArray_real_T *b_mat;
  const double *x_data;
  double *b_mat_data;
  double *mat_data;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  x_data = x->data;
  mat_data = mat->data;
  emxInit_real_T(&b_mat, 2);
  i = b_mat->size[0] * b_mat->size[1];
  if (x->size[0] == 1) {
    b_mat->size[0] = mat->size[0];
  } else {
    b_mat->size[0] = x->size[0];
  }

  b_mat->size[1] = 4;
  emxEnsureCapacity_real_T(b_mat, i);
  b_mat_data = b_mat->data;
  stride_0_0 = (mat->size[0] != 1);
  stride_1_0 = (x->size[0] != 1);
  if (x->size[0] == 1) {
    loop_ub = mat->size[0];
  } else {
    loop_ub = x->size[0];
  }

  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_mat_data[i1 + b_mat->size[0] * i] = mat_data[i1 * stride_0_0 + mat->
        size[0] * i] * x_data[i1 * stride_1_0 + x->size[0] * i] /
        0.69314718055994529;
    }
  }

  i = mat->size[0] * mat->size[1];
  mat->size[0] = b_mat->size[0];
  mat->size[1] = 4;
  emxEnsureCapacity_real_T(mat, i);
  mat_data = mat->data;
  loop_ub = b_mat->size[0];
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      mat_data[i1 + mat->size[0] * i] = b_mat_data[i1 + b_mat->size[0] * i];
    }
  }

  emxFree_real_T(&b_mat);
}

/*
 * function createMEME(fileh,PWM, memefile, GC, C, r, R, rcorr, E, Rd)
 */
static void createMEME(const emxArray_char_T *fileh_Value, const
  emxArray_cell_wrap_2 *PWM, const emxArray_char_T *memefile, double GC, const
  emxArray_real_T *C, double r, const emxArray_real_T *R, double rcorr, const
  emxArray_real_T *E, const emxArray_real_T *Rd)
{
  static const char cv[5] = { 'M', 'O', 'T', 'I', 'F' };

  FILE* b_NULL;
  FILE* filestar;
  cell_wrap_0 *names_data;
  const cell_wrap_2 *PWM_data;
  cell_wrap_4 *cur_PWM_data;
  cell_wrap_4 *p_data;
  emxArray_cell_wrap_0 *names;
  emxArray_cell_wrap_4 *b_cur_PWM;
  emxArray_cell_wrap_4 *cur_PWM;
  emxArray_cell_wrap_4 *p;
  emxArray_char_T *b_memefile;
  emxArray_char_T *text;
  emxArray_int32_T *match_out;
  emxArray_int32_T *matches;
  emxArray_real_T *A;
  emxArray_real_T *b_r;
  emxArray_real_T *mat;
  emxArray_real_T *rmat;
  emxArray_real_T *varargin_1;
  emxArray_real_T *y;
  emxArray_uint32_T *lenvec;
  const double *C_data;
  const double *E_data;
  const double *R_data;
  const double *Rd_data;
  double MM;
  double a;
  double b_varargin_1;
  double *mat_data;
  double *rmat_data;
  double *varargin_1_data;
  int size_tmp[2];
  int b_i;
  int i;
  int i1;
  int i2;
  int j;
  int loop_ub;
  int match_idx;
  unsigned int *lenvec_data;
  int *match_out_data;
  int *matches_data;
  signed char b_fileid;
  signed char fileid;
  char *text_data;
  bool autoflush;
  Rd_data = Rd->data;
  E_data = E->data;
  R_data = R->data;
  C_data = C->data;
  PWM_data = PWM->data;
  emxInit_char_T(&text, 2);
  emxInit_char_T(&b_memefile, 2);

  /* 'gkmPWM:118' num = numel(C); */
  /* 'gkmPWM:119' GC = round(GC*100)/100; */
  GC = rt_roundd(GC * 100.0) / 100.0;

  /* 'gkmPWM:120' GCvec = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  /* 'gkmPWM:121' f = strfind(fileread(memefile),'MOTIF'); */
  fileread(memefile, b_memefile);

  /* 'gkmPWM:122' num2 = length(strfind(fileread(memefile),'MOTIF')); */
  fileread(memefile, text);
  text_data = text->data;
  emxFree_char_T(&b_memefile);
  if (text->size[1] == 0) {
    for (i = 0; i < 2; i++) {
      size_tmp[i] = 1 - i;
    }
  } else {
    emxInit_int32_T(&match_out, 2);
    emxInit_int32_T(&matches, 2);
    b_i = matches->size[0] * matches->size[1];
    matches->size[0] = 1;
    matches->size[1] = text->size[1];
    emxEnsureCapacity_int32_T(matches, b_i);
    matches_data = matches->data;
    match_idx = 0;
    b_i = text->size[1];
    for (i = 0; i <= b_i - 5; i++) {
      j = 1;
      while ((j <= 5) && (text_data[(i + j) - 1] == cv[j - 1])) {
        j++;
      }

      if (j > 5) {
        matches_data[match_idx] = i + 1;
        match_idx++;
      }
    }

    b_i = match_out->size[0] * match_out->size[1];
    match_out->size[0] = 1;
    match_out->size[1] = match_idx;
    emxEnsureCapacity_int32_T(match_out, b_i);
    match_out_data = match_out->data;
    for (i = 0; i < match_idx; i++) {
      match_out_data[i] = matches_data[i];
    }

    emxFree_int32_T(&matches);
    for (b_i = 0; b_i < 2; b_i++) {
      size_tmp[b_i] = match_out->size[b_i];
    }

    emxFree_int32_T(&match_out);
  }

  /* 'gkmPWM:123' [p,names] = getmotif(memefile,1:num2); */
  emxInit_real_T(&y, 2);
  if (size_tmp[1] < 1) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    b_i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = size_tmp[1];
    emxEnsureCapacity_real_T(y, b_i);
    mat_data = y->data;
    loop_ub = size_tmp[1] - 1;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      mat_data[b_i] = (double)b_i + 1.0;
    }
  }

  emxInit_cell_wrap_4(&p);
  emxInit_cell_wrap_0(&names, 1);
  getmotif(memefile, y, p, names);
  names_data = names->data;
  p_data = p->data;

  /* 'gkmPWM:124' fid = fopen(sprintf('%s_denovo.meme', fileh), 'w'); */
  b_sprintf(fileh_Value, text);
  fileid = cfopen(text, "wb");

  /* 'gkmPWM:125' if fid == -1 */
  if (fileid == -1) {
    /* 'gkmPWM:126' fprintf("ERROR: Cannot create gkmPWM denovo motif file.\n"); */
    printf("ERROR: Cannot create gkmPWM denovo motif file.\n");
    fflush(stdout);
    exit(1);
  }

  /* 'gkmPWM:128' fid2 = fopen(sprintf('%s_gkmPWM.out', fileh), 'w'); */
  c_sprintf(fileh_Value, text);
  b_fileid = cfopen(text, "wb");

  /* 'gkmPWM:129' if fid == -1 */
  if (fileid == -1) {
    /* 'gkmPWM:130' fprintf("ERROR: Cannot create gkmPWM output file.\n"); */
    printf("ERROR: Cannot create gkmPWM output file.\n");
    fflush(stdout);
    exit(1);
  }

  emxInit_uint32_T(&lenvec);

  /* 'gkmPWM:133' lenvec = zeros(num2,1); */
  b_i = lenvec->size[0];
  lenvec->size[0] = size_tmp[1];
  emxEnsureCapacity_uint32_T(lenvec, b_i);
  lenvec_data = lenvec->data;

  /* 'gkmPWM:134' for i = 1:num2 */
  b_i = size_tmp[1];
  for (i = 0; i < b_i; i++) {
    /* 'gkmPWM:135' [lenvec(i),~] = size(p{i}); */
    lenvec_data[i] = (unsigned int)p_data[i].f1->size[0];
  }

  /* 'gkmPWM:137' fprintf(fid, 'MEME\n\n'); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "MEME\n\n");
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:138' fprintf(fid, 'ALPHABET= ACGT\n\n'); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "ALPHABET= ACGT\n\n");
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:139' fprintf(fid, 'Correlation with SVM weight vector: %0.3f\n\n', r); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Correlation with SVM weight vector: %0.3f\n\n", r);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:140' fprintf(fid, 'Max PWM Correlation: %0.3f\n\n', rcorr); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Max PWM Correlation: %0.3f\n\n", rcorr);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:141' fprintf(fid, 'Background letter frequencies (from negative set)\n'); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Background letter frequencies (from negative set)\n");
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:142' fprintf(fid, 'A %0.2f C %0.2f G %0.2f T %0.2f\n\n', GCvec(1), GCvec(2), GCvec(3), GCvec(4)); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "A %0.2f C %0.2f G %0.2f T %0.2f\n\n", 0.5 - GC / 2.0, GC /
            2.0, GC / 2.0, 0.5 - GC / 2.0);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:143' fprintf(fid2, 'Correlation with SVM weight vector:\t%0.3f\n', r); */
  b_NULL = NULL;
  getfilestar(b_fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Correlation with SVM weight vector:\t%0.3f\n", r);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:144' fprintf(fid2, 'Max PWM Correlation:\t%0.3f\n', rcorr); */
  b_NULL = NULL;
  getfilestar(b_fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Max PWM Correlation:\t%0.3f\n", rcorr);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:145' fprintf(fid2, 'MOTIF\tID\tSimilarity\tRedundancy\tWeight\tZ\tError\n'); */
  b_NULL = NULL;
  getfilestar(b_fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "MOTIF\tID\tSimilarity\tRedundancy\tWeight\tZ\tError\n");
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:146' a = 1; */
  a = 1.0;

  /* 'gkmPWM:147' b = 1; */
  /* 'gkmPWM:148' for i = 1:num */
  b_i = C->size[0];
  emxInit_cell_wrap_4(&cur_PWM);
  emxInit_cell_wrap_4(&b_cur_PWM);
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&rmat, 2);
  emxInit_real_T(&varargin_1, 2);
  emxInit_real_T(&A, 2);
  emxInit_real_T(&b_r, 2);
  for (i = 0; i < b_i; i++) {
    /* 'gkmPWM:149' [len,~] = size(PWM{i}); */
    /* 'gkmPWM:150' p_len = length(p); */
    /* 'gkmPWM:151' cur_PWM = cell(p_len+1,1); */
    match_idx = p->size[0] + 1;
    i1 = cur_PWM->size[0];
    cur_PWM->size[0] = p->size[0] + 1;
    emxEnsureCapacity_cell_wrap_4(cur_PWM, i1);
    cur_PWM_data = cur_PWM->data;
    for (i1 = 0; i1 < match_idx; i1++) {
      cur_PWM_data[i1].f1->size[0] = 0;
      cur_PWM_data[i1].f1->size[1] = 0;
    }

    /* 'gkmPWM:152' cur_PWM = coder.nullcopy(cur_PWM); */
    i1 = b_cur_PWM->size[0];
    b_cur_PWM->size[0] = cur_PWM->size[0];
    emxEnsureCapacity_cell_wrap_4(b_cur_PWM, i1);
    cur_PWM_data = b_cur_PWM->data;

    /* 'gkmPWM:153' cur_PWM{1} = PWM{i}; */
    i1 = cur_PWM_data[0].f1->size[0] * cur_PWM_data[0].f1->size[1];
    cur_PWM_data[0].f1->size[0] = PWM_data[i].f1->size[0];
    cur_PWM_data[0].f1->size[1] = 4;
    emxEnsureCapacity_real_T(cur_PWM_data[0].f1, i1);
    loop_ub = PWM_data[i].f1->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      cur_PWM_data[0].f1->data[i1] = PWM_data[i].f1->data[i1];
    }

    /* 'gkmPWM:154' for cur_idx=1:p_len */
    i1 = p->size[0];
    for (match_idx = 0; match_idx < i1; match_idx++) {
      /* 'gkmPWM:155' cur_PWM{cur_idx+1} = p{cur_idx}; */
      i2 = cur_PWM_data[match_idx + 1].f1->size[0] * cur_PWM_data[match_idx + 1]
        .f1->size[1];
      cur_PWM_data[match_idx + 1].f1->size[0] = p_data[match_idx].f1->size[0];
      cur_PWM_data[match_idx + 1].f1->size[1] = p_data[match_idx].f1->size[1];
      emxEnsureCapacity_real_T(cur_PWM_data[match_idx + 1].f1, i2);
      loop_ub = p_data[match_idx].f1->size[0] * p_data[match_idx].f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        cur_PWM_data[match_idx + 1].f1->data[i2] = p_data[match_idx].f1->data[i2];
      }
    }

    /* 'gkmPWM:157' [M, MM] = matchMotif(cur_PWM, [len;lenvec]); */
    /* 'gkmPWM:333' n = length(lenvec)-1; */
    /* 'gkmPWM:334' simmat = ones(n-1,1); */
    /* 'gkmPWM:335' for i = 1:n+1 */
    i1 = lenvec->size[0] + 1;
    for (match_idx = 0; match_idx < i1; match_idx++) {
      /* 'gkmPWM:336' mot{i} = mot{i}-1/4; */
      loop_ub = cur_PWM_data[match_idx].f1->size[0] * cur_PWM_data[match_idx].
        f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        cur_PWM_data[match_idx].f1->data[i2] -= 0.25;
      }

      /* 'gkmPWM:337' mot{i} = mot{i}/sqrt(sum(sum(mot{i}.^2))); */
      i2 = A->size[0] * A->size[1];
      A->size[0] = cur_PWM_data[match_idx].f1->size[0];
      A->size[1] = cur_PWM_data[match_idx].f1->size[1];
      emxEnsureCapacity_real_T(A, i2);
      mat_data = A->data;
      loop_ub = cur_PWM_data[match_idx].f1->size[0] * cur_PWM_data[match_idx].
        f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_varargin_1 = cur_PWM_data[match_idx].f1->data[i2];
        mat_data[i2] = pow(b_varargin_1, 2.0);
      }

      c_sum(A, b_r);
      b_varargin_1 = sqrt(sum(b_r));
      loop_ub = cur_PWM_data[match_idx].f1->size[0] * cur_PWM_data[match_idx].
        f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        cur_PWM_data[match_idx].f1->data[i2] /= b_varargin_1;
      }
    }

    /* 'gkmPWM:339' M = 0; */
    b_varargin_1 = 0.0;

    /* 'gkmPWM:340' ind = 1; */
    match_idx = 0;

    /* 'gkmPWM:341' for j = 2:n+1 */
    i1 = lenvec->size[0] + 1;
    for (j = 0; j <= i1 - 2; j++) {
      /* 'gkmPWM:342' mat = mot{1}*mot{j}'; */
      if ((cur_PWM_data[0].f1->size[0] == 0) || (cur_PWM_data[0].f1->size[1] ==
           0) || (cur_PWM_data[j + 1].f1->size[0] == 0) || (cur_PWM_data[j + 1].
           f1->size[1] == 0)) {
        i2 = mat->size[0] * mat->size[1];
        mat->size[0] = cur_PWM_data[0].f1->size[0];
        mat->size[1] = cur_PWM_data[j + 1].f1->size[0];
        emxEnsureCapacity_real_T(mat, i2);
        mat_data = mat->data;
        loop_ub = cur_PWM_data[0].f1->size[0] * cur_PWM_data[j + 1].f1->size[0];
        for (i2 = 0; i2 < loop_ub; i2++) {
          mat_data[i2] = 0.0;
        }
      } else {
        i2 = mat->size[0] * mat->size[1];
        mat->size[0] = cur_PWM_data[0].f1->size[0];
        mat->size[1] = cur_PWM_data[j + 1].f1->size[0];
        emxEnsureCapacity_real_T(mat, i2);
        mat_data = mat->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)
                    cur_PWM_data[0].f1->size[0], (blasint)cur_PWM_data[j + 1].
                    f1->size[0], (blasint)cur_PWM_data[0].f1->size[1], 1.0,
                    &cur_PWM_data[0].f1->data[0], (blasint)cur_PWM_data[0]
                    .f1->size[0], &cur_PWM_data[j + 1].f1->data[0], (blasint)
                    cur_PWM_data[j + 1].f1->size[0], 0.0, &mat_data[0], (blasint)
                    cur_PWM_data[0].f1->size[0]);
      }

      /* 'gkmPWM:343' rmat = rot90(mot{1},2)*mot{j}'; */
      rot90(cur_PWM_data[0].f1, A);
      mat_data = A->data;
      if ((A->size[0] == 0) || (A->size[1] == 0) || (cur_PWM_data[j + 1]
           .f1->size[0] == 0) || (cur_PWM_data[j + 1].f1->size[1] == 0)) {
        i2 = rmat->size[0] * rmat->size[1];
        rmat->size[0] = A->size[0];
        rmat->size[1] = cur_PWM_data[j + 1].f1->size[0];
        emxEnsureCapacity_real_T(rmat, i2);
        rmat_data = rmat->data;
        loop_ub = A->size[0] * cur_PWM_data[j + 1].f1->size[0];
        for (i2 = 0; i2 < loop_ub; i2++) {
          rmat_data[i2] = 0.0;
        }
      } else {
        i2 = rmat->size[0] * rmat->size[1];
        rmat->size[0] = A->size[0];
        rmat->size[1] = cur_PWM_data[j + 1].f1->size[0];
        emxEnsureCapacity_real_T(rmat, i2);
        rmat_data = rmat->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)A->size[0],
                    (blasint)cur_PWM_data[j + 1].f1->size[0], (blasint)A->size[1],
                    1.0, &mat_data[0], (blasint)A->size[0], &cur_PWM_data[j + 1]
                    .f1->data[0], (blasint)cur_PWM_data[j + 1].f1->size[0], 0.0,
                    &rmat_data[0], (blasint)A->size[0]);
      }

      /* 'gkmPWM:344' MM = max([sum(spdiags(mat)) sum(spdiags(rmat))]); */
      spdiags(mat, A);
      c_sum(A, y);
      mat_data = y->data;
      spdiags(rmat, A);
      c_sum(A, b_r);
      rmat_data = b_r->data;
      i2 = varargin_1->size[0] * varargin_1->size[1];
      varargin_1->size[0] = 1;
      varargin_1->size[1] = y->size[1] + b_r->size[1];
      emxEnsureCapacity_real_T(varargin_1, i2);
      varargin_1_data = varargin_1->data;
      loop_ub = y->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        varargin_1_data[i2] = mat_data[i2];
      }

      loop_ub = b_r->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        varargin_1_data[i2 + y->size[1]] = rmat_data[i2];
      }

      MM = h_maximum(varargin_1);

      /* 'gkmPWM:345' if MM > M */
      if (MM > b_varargin_1) {
        /* 'gkmPWM:346' M = MM; */
        b_varargin_1 = MM;

        /* 'gkmPWM:347' ind = j-1; */
        match_idx = j;
      }
    }

    /*  matches motifs to the best motif in our database */
    /*  [M, MM] = matchMotif([PWM{i}; p], [len;lenvec]);% matches motifs to the best motif */
    /* 'gkmPWM:159' fprintf(fid, 'MOTIF %d %s\n', int32(a), names{MM}); */
    i1 = text->size[0] * text->size[1];
    text->size[0] = 1;
    text->size[1] = names_data[match_idx].f1->size[1] + 1;
    emxEnsureCapacity_char_T(text, i1);
    text_data = text->data;
    loop_ub = names_data[match_idx].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      text_data[i1] = names_data[match_idx].f1->data[i1];
    }

    text_data[names_data[match_idx].f1->size[1]] = '\x00';
    b_NULL = NULL;
    getfilestar(fileid, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "MOTIF %d %s\n", (int)a, &text_data[0]);
      if (autoflush) {
        fflush(filestar);
      }
    }

    /* 'gkmPWM:160' fprintf(fid2,'%s\t%d\t%0.3f\t%0.3f\t%0.2f\t%0.3f\t%0.3f\n',names{MM},int32(MM),M,Rd(i),C(i),R(i),E(i)); */
    i1 = text->size[0] * text->size[1];
    text->size[0] = 1;
    text->size[1] = names_data[match_idx].f1->size[1] + 1;
    emxEnsureCapacity_char_T(text, i1);
    text_data = text->data;
    loop_ub = names_data[match_idx].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      text_data[i1] = names_data[match_idx].f1->data[i1];
    }

    text_data[names_data[match_idx].f1->size[1]] = '\x00';
    b_NULL = NULL;
    getfilestar(b_fileid, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "%s\t%d\t%0.3f\t%0.3f\t%0.2f\t%0.3f\t%0.3f\n",
              &text_data[0], match_idx + 1, b_varargin_1, Rd_data[i], C_data[i],
              R_data[i], E_data[i]);
      if (autoflush) {
        fflush(filestar);
      }
    }

    /* 'gkmPWM:161' fprintf(fid, 'weight= %0.3f l= 4 w= %d z-score= %0.2f motifsim= %0.3f\n', C(i), int32(len), R(i), M); */
    b_NULL = NULL;
    getfilestar(fileid, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar,
              "weight= %0.3f l= 4 w= %d z-score= %0.2f motifsim= %0.3f\n",
              C_data[i], PWM_data[i].f1->size[0], R_data[i], b_varargin_1);
      if (autoflush) {
        fflush(filestar);
      }
    }

    /* 'gkmPWM:162' for j = 1:len */
    i1 = PWM_data[i].f1->size[0];
    for (j = 0; j < i1; j++) {
      /* 'gkmPWM:163' fprintf(fid, '%0.3f %0.3f %0.3f %0.3f\n',PWM{i}(j,1),PWM{i}(j,2),PWM{i}(j,3),PWM{i}(j,4)); */
      b_NULL = NULL;
      getfilestar(fileid, &filestar, &autoflush);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "%0.3f %0.3f %0.3f %0.3f\n", PWM_data[i].f1->data[j],
                PWM_data[i].f1->data[j + PWM_data[i].f1->size[0]], PWM_data[i].
                f1->data[j + PWM_data[i].f1->size[0] * 2], PWM_data[i].f1->
                data[j + PWM_data[i].f1->size[0] * 3]);
        if (autoflush) {
          fflush(filestar);
        }
      }
    }

    /* 'gkmPWM:165' fprintf(fid, '\n'); */
    b_NULL = NULL;
    getfilestar(fileid, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "\n");
      if (autoflush) {
        fflush(filestar);
      }
    }

    /* 'gkmPWM:166' a = a+1; */
    a++;
  }

  emxFree_real_T(&b_r);
  emxFree_real_T(&A);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&rmat);
  emxFree_real_T(&mat);
  emxFree_real_T(&y);
  emxFree_char_T(&text);
  emxFree_cell_wrap_4(&b_cur_PWM);
  emxFree_cell_wrap_0(&names);
  emxFree_cell_wrap_4(&p);
  emxFree_cell_wrap_4(&cur_PWM);
  emxFree_uint32_T(&lenvec);

  /* 'gkmPWM:168' fclose(fid); */
  cfclose(fileid);

  /* 'gkmPWM:169' fclose(fid2); */
  cfclose(b_fileid);
}

static void db_binary_expand_op(emxArray_real_T *mat, const emxArray_cell_wrap_2
  *p, int i, int i2, int i3, int i4, int i5)
{
  const cell_wrap_2 *p_data;
  double *mat_data;
  int b_i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  p_data = p->data;
  b_i = mat->size[0] * mat->size[1];
  if ((i5 - i4) + 1 == 1) {
    mat->size[0] = (i3 - i2) + 1;
  } else {
    mat->size[0] = (i5 - i4) + 1;
  }

  mat->size[1] = 4;
  emxEnsureCapacity_real_T(mat, b_i);
  mat_data = mat->data;
  stride_0_0 = ((i3 - i2) + 1 != 1);
  stride_1_0 = ((i5 - i4) + 1 != 1);
  if ((i5 - i4) + 1 == 1) {
    loop_ub = (i3 - i2) + 1;
  } else {
    loop_ub = (i5 - i4) + 1;
  }

  for (b_i = 0; b_i < 4; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      mat_data[i1 + mat->size[0] * b_i] = p_data[i].f1->data[(i2 + i1 *
        stride_0_0) + p_data[i].f1->size[0] * b_i] + (double)(p_data[i].f1->
        data[(i4 + i1 * stride_1_0) + p_data[i].f1->size[0] * b_i] == 0.0);
    }
  }
}

static void f_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3)
{
  const cell_wrap_1 *mat_data;
  const double *rs_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  bool *x_data;
  rs_data = rs->data;
  mat_data = mat->data;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  if ((i3 - i2) + 1 == 1) {
    x->size[1] = i1 + 1;
  } else {
    x->size[1] = (i3 - i2) + 1;
  }

  emxEnsureCapacity_boolean_T(x, i);
  x_data = x->data;
  stride_0_1 = (i1 + 1 != 1);
  stride_1_1 = ((i3 - i2) + 1 != 1);
  if ((i3 - i2) + 1 == 1) {
    loop_ub = i1 + 1;
  } else {
    loop_ub = (i3 - i2) + 1;
  }

  for (i = 0; i < loop_ub; i++) {
    x_data[i] = (mat_data[j].f1->data[i * stride_0_1] == rs_data[i2 + i *
                 stride_1_1]);
  }
}

static void g_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3)
{
  const cell_wrap_1 *mat_data;
  const double *rs_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  bool *x_data;
  rs_data = rs->data;
  mat_data = mat->data;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  if ((i3 - i2) + 1 == 1) {
    x->size[1] = i1 + 2;
  } else {
    x->size[1] = (i3 - i2) + 1;
  }

  emxEnsureCapacity_boolean_T(x, i);
  x_data = x->data;
  stride_0_1 = (i1 + 2 != 1);
  stride_1_1 = ((i3 - i2) + 1 != 1);
  if ((i3 - i2) + 1 == 1) {
    loop_ub = i1 + 2;
  } else {
    loop_ub = (i3 - i2) + 1;
  }

  for (i = 0; i < loop_ub; i++) {
    x_data[i] = (mat_data[j].f1->data[i * stride_0_1] == rs_data[i2 + i *
                 stride_1_1]);
  }
}

/*
 * function [kweig,P] = getEMprob_v3(PWM,res,negmat,poscell,rc,diffc,indc,indloc,xc,reg,l_svm,k_svm,rcnum,RC)
 */
static void getEMprob_v3(const emxArray_real_T *PWM, const emxArray_real_T *res,
  const double negmat[16], const emxArray_cell_wrap_14 *poscell, const
  emxArray_real_T *rc, const emxArray_real_T *diffc, const emxArray_real_T *indc,
  const emxArray_real_T *indloc, const emxArray_real_T *xc, double reg, double
  l_svm, double k_svm, double rcnum, double RC, emxArray_real_T *kweig, double
  P_data[], int *P_size)
{
  static const signed char b_b[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 1 };

  static const signed char d_b[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  static const signed char c_b[4] = { 1, 0, 0, 1 };

  emxArray_int8_T *b_x;
  emxArray_int8_T *x;
  emxArray_real_T b_MAT_data;
  emxArray_real_T b_p_data;
  emxArray_real_T e_x_data;
  emxArray_real_T *A;
  emxArray_real_T *b;
  emxArray_real_T *c_x;
  emxArray_real_T *vec;
  creal_T B[16];
  creal_T b_mat[16];
  creal_T d_x[9];
  creal_T dcv[9];
  creal_T t_data[8];
  creal_T MAT[4];
  creal_T b_y[4];
  creal_T ps[4];
  creal_T dcv2[3];
  creal_T dcv1[2];
  creal_T E;
  creal_T M;
  creal_T t;
  double MAT_data[16];
  double mat[16];
  double tmp_data[16];
  double d_x_data[12];
  double Posvec_data[4];
  double p_data[4];
  double y[4];
  const double *PWM_data;
  const double *res_data;
  double B_re_tmp;
  double b_B_re_tmp;
  double bim;
  double brm;
  double c_B_re_tmp;
  double d_B_re_tmp;
  double e_re;
  double im;
  double re;
  double *A_data;
  double *b_data;
  double *c_x_data;
  double *vec_data;
  int b_tmp_data[4];
  int B_size[2];
  int MAT_size[2];
  int Posvec_size[2];
  int b_x_size[2];
  int tmp_size[2];
  int x_size[2];
  int b_I;
  int b_i;
  int i;
  int i1;
  int iindx;
  int k;
  int last;
  int loop_ub;
  int p_size;
  signed char ind[6];
  signed char c_tmp_data[3];
  signed char *b_x_data;
  signed char *x_data;
  res_data = res->data;
  PWM_data = PWM->data;

  /* Lagrange optimization (see paper for derivation) */
  /* 'gkmPWM:595' a = true; */
  /* 'gkmPWM:596' posvec = 1:4; */
  /* 'gkmPWM:597' if RC */
  emxInit_real_T(&A, 2);
  if (RC != 0.0) {
    /* 'gkmPWM:598' A =  ls_kweigtree(PWM,negmat,poscell,rc,diffc,indc,indloc,xc,l_svm,k_svm,rcnum); */
    ls_kweigtree(PWM, negmat, poscell, rc, diffc, indc, indloc, xc, l_svm, k_svm,
                 rcnum, A);
    A_data = A->data;
  } else {
    /* 'gkmPWM:599' else */
    /* 'gkmPWM:600' A =  ls_kweigtree_norc(PWM,negmat,poscell,rc,diffc,indc,indloc,xc,l_svm,k_svm,rcnum); */
    ls_kweigtree_norc(PWM, negmat, poscell, rc, diffc, indc, indloc, xc, l_svm,
                      k_svm, A);
    A_data = A->data;
  }

  /* 'gkmPWM:602' b = res+A*PWM(l_svm,:)'; */
  for (i = 0; i < 4; i++) {
    y[i] = PWM_data[((int)l_svm + PWM->size[0] * i) - 1];
  }

  emxInit_real_T(&b, 1);
  if (A->size[0] == 0) {
    b->size[0] = 0;
  } else {
    i = b->size[0];
    b->size[0] = A->size[0];
    emxEnsureCapacity_real_T(b, i);
    b_data = b->data;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, (blasint)A->size[0],
                (blasint)1, (blasint)4, 1.0, &A_data[0], (blasint)A->size[0],
                &y[0], (blasint)1, 0.0, &b_data[0], (blasint)A->size[0]);
  }

  if (res->size[0] == b->size[0]) {
    i = b->size[0];
    b->size[0] = res->size[0];
    emxEnsureCapacity_real_T(b, i);
    b_data = b->data;
    loop_ub = res->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_data[i] += res_data[i];
    }
  } else {
    b_plus(b, res);
    b_data = b->data;
  }

  /* 'gkmPWM:603' mat = A'*A; */
  if (A->size[0] == 0) {
    memset(&mat[0], 0, 16U * sizeof(double));
  } else {
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)4, (blasint)4,
                (blasint)A->size[0], 1.0, &A_data[0], (blasint)A->size[0],
                &A_data[0], (blasint)A->size[0], 0.0, &mat[0], (blasint)4);
  }

  /* 'gkmPWM:604' y = A'*b; */
  c_mtimes(A, b, y);

  /* 'gkmPWM:605' if reg > 0 */
  if (reg > 0.0) {
    /* 'gkmPWM:606' M = min(eig(mat)); */
    eig(mat, b_y);
    M = c_minimum(b_y);

    /* 'gkmPWM:607' B=(mat-reg*M*eye(4))^-1; */
    E.re = reg * M.re;
    E.im = reg * M.im;
    for (i = 0; i < 16; i++) {
      i1 = b_b[i];
      b_mat[i].re = mat[i] - E.re * (double)i1;
      b_mat[i].im = 0.0 - E.im * (double)i1;
    }

    b_mpower(b_mat, B);
  } else {
    /* 'gkmPWM:608' else */
    /* 'gkmPWM:609' M = 0; */
    M.re = 0.0;
    M.im = 0.0;

    /* 'gkmPWM:610' B = mat^-1; */
    c_mpower(mat, tmp_data);
    for (i = 0; i < 16; i++) {
      B[i].re = tmp_data[i];
      B[i].im = 0.0;
    }
  }

  /* 'gkmPWM:612' ps = B*y; */
  for (i = 0; i < 4; i++) {
    b_y[i].re = y[i];
    b_y[i].im = 0.0;
  }

  for (i = 0; i < 4; i++) {
    re = 0.0;
    im = 0.0;
    for (i1 = 0; i1 < 4; i1++) {
      last = i + (i1 << 2);
      B_re_tmp = B[last].re;
      b_B_re_tmp = b_y[i1].im;
      c_B_re_tmp = B[last].im;
      d_B_re_tmp = b_y[i1].re;
      re += B_re_tmp * d_B_re_tmp - c_B_re_tmp * b_B_re_tmp;
      im += B_re_tmp * b_B_re_tmp + c_B_re_tmp * d_B_re_tmp;
    }

    ps[i].re = re;
    ps[i].im = im;
  }

  /* 'gkmPWM:613' coder.varsize('p'); */
  /* 'gkmPWM:614' p = real(ps+(1-sum(ps))*B/sum(sum(B))*ones(4,1)); */
  E = ps[0];
  for (k = 0; k < 3; k++) {
    E.re += ps[k + 1].re;
    E.im += ps[k + 1].im;
  }

  t.re = 1.0 - E.re;
  t.im = 0.0 - E.im;
  e_sum(B, MAT);
  E = MAT[0];
  for (k = 0; k < 3; k++) {
    E.re += MAT[k + 1].re;
    E.im += MAT[k + 1].im;
  }

  for (i = 0; i < 16; i++) {
    im = B[i].im;
    re = B[i].re;
    b_B_re_tmp = t.re * re - t.im * im;
    im = t.re * im + t.im * re;
    if (E.im == 0.0) {
      if (im == 0.0) {
        B[i].re = b_B_re_tmp / E.re;
        B[i].im = 0.0;
      } else if (b_B_re_tmp == 0.0) {
        B[i].re = 0.0;
        B[i].im = im / E.re;
      } else {
        B[i].re = b_B_re_tmp / E.re;
        B[i].im = im / E.re;
      }
    } else if (E.re == 0.0) {
      if (b_B_re_tmp == 0.0) {
        B[i].re = im / E.im;
        B[i].im = 0.0;
      } else if (im == 0.0) {
        B[i].re = 0.0;
        B[i].im = -(b_B_re_tmp / E.im);
      } else {
        B[i].re = im / E.im;
        B[i].im = -(b_B_re_tmp / E.im);
      }
    } else {
      brm = fabs(E.re);
      bim = fabs(E.im);
      if (brm > bim) {
        B_re_tmp = E.im / E.re;
        re = E.re + B_re_tmp * E.im;
        B[i].re = (b_B_re_tmp + B_re_tmp * im) / re;
        B[i].im = (im - B_re_tmp * b_B_re_tmp) / re;
      } else if (bim == brm) {
        if (E.re > 0.0) {
          re = 0.5;
        } else {
          re = -0.5;
        }

        if (E.im > 0.0) {
          B_re_tmp = 0.5;
        } else {
          B_re_tmp = -0.5;
        }

        B[i].re = (b_B_re_tmp * re + im * B_re_tmp) / brm;
        B[i].im = (im * re - b_B_re_tmp * B_re_tmp) / brm;
      } else {
        B_re_tmp = E.re / E.im;
        re = E.im + B_re_tmp * E.re;
        B[i].re = (B_re_tmp * b_B_re_tmp + im) / re;
        B[i].im = (B_re_tmp * im - b_B_re_tmp) / re;
      }
    }
  }

  for (i = 0; i < 4; i++) {
    b_y[i].re = 1.0;
    b_y[i].im = 0.0;
  }

  p_size = 4;

  /* 'gkmPWM:615' coder.varsize('P'); */
  /* 'gkmPWM:616' P = zeros(4,1); */
  *P_size = 4;
  for (i = 0; i < 4; i++) {
    re = 0.0;
    for (i1 = 0; i1 < 4; i1++) {
      last = i + (i1 << 2);
      re += B[last].re - B[last].im * b_y[i1].im;
    }

    p_data[i] = ps[i].re + re;
    P_data[i] = 0.0;
  }

  /* 'gkmPWM:617' e = 1e5; */
  e_re = 100000.0;

  /* The solution to the lagrange optimization problem */
  /* The following deals with the case if the optimal solution has a negative probability.  I cheat a little by only considering the cases where the base with the maximum solution is non-zero.  This works just fine in practice since the sum(p) = 1 constraint forces one of the bases to be positive.  It speeds up computation almost two fold. */
  /* 'gkmPWM:620' if min(p) < 0 */
  im = p_data[0];
  for (k = 0; k < 3; k++) {
    re = p_data[k + 1];
    if (im > re) {
      im = re;
    }
  }

  if (im < 0.0) {
    emxInit_int8_T(&x, 2);

    /* 'gkmPWM:621' I = 0; */
    b_I = 0;

    /* 'gkmPWM:622' [~,a] = max(p); */
    b_p_data.data = &p_data[0];
    b_p_data.size = &p_size;
    b_p_data.allocatedSize = 4;
    b_p_data.numDimensions = 1;
    b_p_data.canFreeData = false;
    d_maximum(&b_p_data, &im, &iindx);

    /* 'gkmPWM:623' nvec = posvec; */
    /* 'gkmPWM:624' nvec(a) = []; */
    i = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = 4;
    emxEnsureCapacity_int8_T(x, i);
    x_data = x->data;
    for (i = 0; i < 4; i++) {
      x_data[i] = (signed char)(i + 1);
    }

    for (k = iindx; k < 4; k++) {
      x_data[k - 1] = x_data[k];
    }

    i = x->size[0] * x->size[1];
    x->size[1] = 3;
    emxEnsureCapacity_int8_T(x, i);
    x_data = x->data;

    /* Check cases where one of the probabilities is zero */
    /* 'gkmPWM:626' for i = 1:3 */
    emxInit_real_T(&vec, 1);
    vec_data = vec->data;
    emxInit_int8_T(&b_x, 2);
    emxInit_real_T(&c_x, 2);
    Posvec_size[1] = 3;
    for (b_i = 0; b_i < 3; b_i++) {
      /* 'gkmPWM:627' Posvec = posvec; */
      /* 'gkmPWM:628' Posvec(nvec(i)) = []; */
      i = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 1;
      b_x->size[1] = 4;
      emxEnsureCapacity_int8_T(b_x, i);
      b_x_data = b_x->data;
      for (i = 0; i < 4; i++) {
        b_x_data[i] = (signed char)(i + 1);
      }

      last = x_data[b_i];
      for (k = last; k < 4; k++) {
        b_x_data[k - 1] = b_x_data[k];
      }

      i = b_x->size[0] * b_x->size[1];
      b_x->size[1] = 3;
      emxEnsureCapacity_int8_T(b_x, i);
      b_x_data = b_x->data;

      /* 'gkmPWM:629' MAT = mat; */
      /* 'gkmPWM:630' MAT(:,nvec(i)) =[]; */
      /* 'gkmPWM:631' MAT(nvec(i),:) =[]; */
      i = c_x->size[0] * c_x->size[1];
      c_x->size[0] = 4;
      c_x->size[1] = 4;
      emxEnsureCapacity_real_T(c_x, i);
      c_x_data = c_x->data;
      for (i = 0; i < 16; i++) {
        c_x_data[i] = mat[i];
      }

      for (k = last; k < 4; k++) {
        for (loop_ub = 0; loop_ub < 4; loop_ub++) {
          c_x_data[loop_ub + c_x->size[0] * (k - 1)] = c_x_data[loop_ub +
            c_x->size[0] * k];
        }
      }

      for (i = 0; i < 3; i++) {
        for (i1 = 0; i1 < 4; i1++) {
          c_x_data[i1 + (i << 2)] = c_x_data[i1 + c_x->size[0] * i];
        }
      }

      i = c_x->size[0] * c_x->size[1];
      c_x->size[0] = 4;
      c_x->size[1] = 3;
      emxEnsureCapacity_real_T(c_x, i);
      c_x_data = c_x->data;
      b_x_size[0] = 4;
      b_x_size[1] = 3;
      for (i = 0; i < 12; i++) {
        d_x_data[i] = c_x_data[i];
      }

      for (k = 0; k < 3; k++) {
        for (loop_ub = last; loop_ub < 4; loop_ub++) {
          d_x_data[(loop_ub + b_x_size[0] * k) - 1] = d_x_data[loop_ub +
            b_x_size[0] * k];
        }
      }

      for (i = 0; i < 3; i++) {
        for (i1 = 0; i1 < 3; i1++) {
          d_x_data[i1 + 3 * i] = d_x_data[i1 + b_x_size[0] * i];
        }
      }

      b_x_size[0] = 3;
      b_x_size[1] = 3;

      /* 'gkmPWM:632' Y = y(Posvec); */
      /* 'gkmPWM:633' if reg > 0 */
      if (reg > 0.0) {
        /* 'gkmPWM:634' B=(MAT-reg*M*eye(3))^-1; */
        E.re = reg * M.re;
        E.im = reg * M.im;
        for (i = 0; i < 3; i++) {
          for (i1 = 0; i1 < 3; i1++) {
            last = i1 + 3 * i;
            k = d_b[last];
            d_x[last].re = d_x_data[i1 + b_x_size[0] * i] - E.re * (double)k;
            d_x[last].im = 0.0 - E.im * (double)k;
          }
        }

        d_mpower(d_x, dcv);
        B_size[0] = 3;
        B_size[1] = 3;
        memcpy(&B[0], &dcv[0], 9U * sizeof(creal_T));
      } else {
        /* 'gkmPWM:635' else */
        /* 'gkmPWM:636' B=(MAT)^-1; */
        e_x_data.data = &d_x_data[0];
        e_x_data.size = &b_x_size[0];
        e_x_data.allocatedSize = 12;
        e_x_data.numDimensions = 2;
        e_x_data.canFreeData = false;
        mpower(&e_x_data, c_x);
        c_x_data = c_x->data;
        B_size[0] = c_x->size[0];
        B_size[1] = c_x->size[1];
        loop_ub = c_x->size[0] * c_x->size[1];
        for (i = 0; i < loop_ub; i++) {
          B[i].re = c_x_data[i];
          B[i].im = 0.0;
        }
      }

      /* 'gkmPWM:638' ps = B*Y; */
      for (i = 0; i < 3; i++) {
        b_y[i].re = y[b_x_data[i] - 1];
        b_y[i].im = 0.0;
      }

      loop_ub = B_size[0];
      for (i = 0; i < loop_ub; i++) {
        ps[i].re = 0.0;
        ps[i].im = 0.0;
        last = B_size[1];
        for (i1 = 0; i1 < last; i1++) {
          k = i + B_size[0] * i1;
          B_re_tmp = B[k].re;
          b_B_re_tmp = b_y[i1].im;
          c_B_re_tmp = B[k].im;
          d_B_re_tmp = b_y[i1].re;
          ps[i].re += B_re_tmp * d_B_re_tmp - c_B_re_tmp * b_B_re_tmp;
          ps[i].im += B_re_tmp * b_B_re_tmp + c_B_re_tmp * d_B_re_tmp;
        }
      }

      /* 'gkmPWM:639' p = real(ps+(1-sum(ps))*B/sum(sum(B))*ones(3,1)); */
      E = ps[0];
      for (k = 0; k < 2; k++) {
        E.re += ps[k + 1].re;
        E.im += ps[k + 1].im;
      }

      t.re = 1.0 - E.re;
      t.im = 0.0 - E.im;
      f_sum(B, B_size, b_y, x_size);
      last = x_size[1];
      if (x_size[1] == 0) {
        E.re = 0.0;
        E.im = 0.0;
      } else {
        E = b_y[0];
        for (k = 2; k <= last; k++) {
          E.re += b_y[k - 1].re;
          E.im += b_y[k - 1].im;
        }
      }

      if (B_size[0] == 3) {
        loop_ub = B_size[0] * B_size[1];
        for (i = 0; i < loop_ub; i++) {
          re = B[i].re;
          im = B[i].im;
          b_B_re_tmp = t.re * re - t.im * im;
          im = t.re * im + t.im * re;
          if (E.im == 0.0) {
            if (im == 0.0) {
              d_x[i].re = b_B_re_tmp / E.re;
              d_x[i].im = 0.0;
            } else if (b_B_re_tmp == 0.0) {
              d_x[i].re = 0.0;
              d_x[i].im = im / E.re;
            } else {
              d_x[i].re = b_B_re_tmp / E.re;
              d_x[i].im = im / E.re;
            }
          } else if (E.re == 0.0) {
            if (b_B_re_tmp == 0.0) {
              d_x[i].re = im / E.im;
              d_x[i].im = 0.0;
            } else if (im == 0.0) {
              d_x[i].re = 0.0;
              d_x[i].im = -(b_B_re_tmp / E.im);
            } else {
              d_x[i].re = im / E.im;
              d_x[i].im = -(b_B_re_tmp / E.im);
            }
          } else {
            brm = fabs(E.re);
            bim = fabs(E.im);
            if (brm > bim) {
              B_re_tmp = E.im / E.re;
              re = E.re + B_re_tmp * E.im;
              d_x[i].re = (b_B_re_tmp + B_re_tmp * im) / re;
              d_x[i].im = (im - B_re_tmp * b_B_re_tmp) / re;
            } else if (bim == brm) {
              if (E.re > 0.0) {
                re = 0.5;
              } else {
                re = -0.5;
              }

              if (E.im > 0.0) {
                B_re_tmp = 0.5;
              } else {
                B_re_tmp = -0.5;
              }

              d_x[i].re = (b_B_re_tmp * re + im * B_re_tmp) / brm;
              d_x[i].im = (im * re - b_B_re_tmp * B_re_tmp) / brm;
            } else {
              B_re_tmp = E.re / E.im;
              re = E.im + B_re_tmp * E.re;
              d_x[i].re = (B_re_tmp * b_B_re_tmp + im) / re;
              d_x[i].im = (B_re_tmp * im - b_B_re_tmp) / re;
            }
          }
        }

        for (i = 0; i < 3; i++) {
          dcv2[i].re = 1.0;
          dcv2[i].im = 0.0;
        }

        p_size = 3;
        for (i = 0; i < 3; i++) {
          re = 0.0;
          for (i1 = 0; i1 < 3; i1++) {
            last = i + 3 * i1;
            re += d_x[last].re - d_x[last].im * dcv2[i1].im;
          }

          p_data[i] = ps[i].re + re;
        }
      } else {
        w_binary_expand_op(p_data, &p_size, ps, &B_size[0], t, B, B_size, E);
      }

      /* solution */
      /* Checks if solution is permitted.  If so, makes sure that it creates a smaller error than other cases */
      /* 'gkmPWM:641' if min(p) >= 0 */
      im = p_data[0];
      for (k = 0; k < 2; k++) {
        re = p_data[k + 1];
        if (im > re) {
          im = re;
        }
      }

      if (im >= 0.0) {
        /* 'gkmPWM:642' vec = b-A(:,Posvec)*p; */
        loop_ub = A->size[0];
        i = c_x->size[0] * c_x->size[1];
        c_x->size[0] = A->size[0];
        c_x->size[1] = 3;
        emxEnsureCapacity_real_T(c_x, i);
        c_x_data = c_x->data;
        for (i = 0; i < 3; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            c_x_data[i1 + c_x->size[0] * i] = A_data[i1 + A->size[0] *
              (b_x_data[i] - 1)];
          }
        }

        loop_ub = A->size[0];
        if (A->size[0] == 0) {
          vec->size[0] = 0;
          for (i = 0; i < loop_ub; i++) {
            vec_data[i] = 0.0;
          }
        } else {
          i = vec->size[0];
          vec->size[0] = A->size[0];
          emxEnsureCapacity_real_T(vec, i);
          vec_data = vec->data;
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                      A->size[0], (blasint)1, (blasint)3, 1.0, &c_x_data[0],
                      (blasint)A->size[0], &p_data[0], (blasint)3, 0.0,
                      &vec_data[0], (blasint)A->size[0]);
        }

        loop_ub = b->size[0];
        if (b->size[0] == vec->size[0]) {
          i = vec->size[0];
          vec->size[0] = b->size[0];
          emxEnsureCapacity_real_T(vec, i);
          vec_data = vec->data;
          for (i = 0; i < loop_ub; i++) {
            vec_data[i] = b_data[i] - vec_data[i];
          }
        } else {
          minus(vec, b);
          vec_data = vec->data;
        }

        /* 'gkmPWM:643' if I == 0 */
        if (b_I == 0) {
          /* 'gkmPWM:644' I = 1; */
          b_I = 1;

          /* 'gkmPWM:645' P = zeros(4,1); */
          *P_size = 4;
          for (i = 0; i < 4; i++) {
            P_data[i] = 0.0;
          }

          /* 'gkmPWM:646' P(Posvec) = real(p); */
          for (i = 0; i < 3; i++) {
            c_tmp_data[i] = b_x_data[i];
          }

          for (i = 0; i < 3; i++) {
            P_data[c_tmp_data[i] - 1] = p_data[i];
          }

          /* 'gkmPWM:647' e = vec'*vec-reg*M*p'*p; */
          E.re = reg * M.re;
          if (vec->size[0] < 1) {
            im = 0.0;
          } else {
            im = cblas_ddot((blasint)vec->size[0], &vec_data[0], (blasint)1,
                            &vec_data[0], (blasint)1);
          }

          re = 0.0;
          loop_ub = p_size;
          for (i = 0; i < loop_ub; i++) {
            re += p_data[i] * E.re * p_data[i];
          }

          e_re = im - re;
        } else {
          /* 'gkmPWM:648' else */
          /* 'gkmPWM:649' E = vec'*vec-reg*M*p'*p; */
          E.re = reg * M.re;
          if (vec->size[0] < 1) {
            im = 0.0;
          } else {
            im = cblas_ddot((blasint)vec->size[0], &vec_data[0], (blasint)1,
                            &vec_data[0], (blasint)1);
          }

          re = 0.0;
          loop_ub = p_size;
          for (i = 0; i < loop_ub; i++) {
            re += p_data[i] * E.re * p_data[i];
          }

          E.re = im - re;

          /* 'gkmPWM:650' if E < e */
          if (E.re < e_re) {
            /* 'gkmPWM:651' e = E; */
            e_re = E.re;

            /* 'gkmPWM:652' P = zeros(4,1); */
            *P_size = 4;
            for (i = 0; i < 4; i++) {
              P_data[i] = 0.0;
            }

            /* 'gkmPWM:653' P(Posvec) = real(p); */
            for (i = 0; i < 3; i++) {
              c_tmp_data[i] = b_x_data[i];
            }

            for (i = 0; i < 3; i++) {
              P_data[c_tmp_data[i] - 1] = p_data[i];
            }
          }
        }
      }
    }

    emxFree_int8_T(&b_x);

    /* Check cases where two of the probabilities are zero */
    /* 'gkmPWM:659' ind = [nvec(1) nvec(2);nvec(1) nvec(3);nvec(2) nvec(3)]; */
    ind[0] = x_data[0];
    ind[3] = x_data[1];
    ind[1] = x_data[0];
    ind[4] = x_data[2];
    ind[2] = x_data[1];
    ind[5] = x_data[2];

    /* 'gkmPWM:660' for i = 1:3 */
    emxFree_int8_T(&x);
    for (b_i = 0; b_i < 3; b_i++) {
      /* 'gkmPWM:661' Posvec = posvec; */
      /* 'gkmPWM:662' Posvec(ind(i,:)) = []; */
      for (i = 0; i < 2; i++) {
        x_size[i] = ind[b_i + 3 * i];
      }

      Posvec_size[0] = 1;
      Posvec_size[1] = 4;
      for (i = 0; i < 4; i++) {
        Posvec_data[i] = (double)i + 1.0;
      }

      c_nullAssignment(Posvec_data, Posvec_size, x_size);

      /* 'gkmPWM:663' MAT = mat; */
      /* 'gkmPWM:664' MAT(:,ind(i,:)) =[]; */
      /* 'gkmPWM:665' MAT(ind(i,:),:) =[]; */
      tmp_size[0] = 4;
      tmp_size[1] = 4;
      memcpy(&tmp_data[0], &mat[0], 16U * sizeof(double));
      d_nullAssignment(tmp_data, tmp_size, x_size);
      MAT_size[0] = 4;
      MAT_size[1] = tmp_size[1];
      loop_ub = 4 * tmp_size[1];
      if (0 <= loop_ub - 1) {
        memcpy(&MAT_data[0], &tmp_data[0], loop_ub * sizeof(double));
      }

      e_nullAssignment(MAT_data, MAT_size, x_size);

      /* 'gkmPWM:666' Y = y(Posvec); */
      /* 'gkmPWM:667' if reg > 0 */
      if (reg > 0.0) {
        /* 'gkmPWM:668' B=(MAT-reg*M*eye(2))^-1; */
        E.re = reg * M.re;
        E.im = reg * M.im;
        if ((MAT_size[0] == 2) && (MAT_size[1] == 2)) {
          for (i = 0; i < 2; i++) {
            for (i1 = 0; i1 < 2; i1++) {
              last = i1 + (i << 1);
              k = c_b[last];
              MAT[last].re = MAT_data[i1 + MAT_size[0] * i] - E.re * (double)k;
              MAT[last].im = 0.0 - E.im * (double)k;
            }
          }
        } else {
          x_binary_expand_op(MAT, MAT_data, MAT_size, E, c_b);
        }

        brm = fabs(MAT[0].re);
        c_B_re_tmp = fabs(MAT[1].re);
        bim = fabs(MAT[0].im);
        d_B_re_tmp = fabs(MAT[1].im);
        if (c_B_re_tmp + d_B_re_tmp > brm + bim) {
          if (MAT[1].im == 0.0) {
            if (MAT[0].im == 0.0) {
              E.re = MAT[0].re / MAT[1].re;
              E.im = 0.0;
            } else if (MAT[0].re == 0.0) {
              E.re = 0.0;
              E.im = MAT[0].im / MAT[1].re;
            } else {
              E.re = MAT[0].re / MAT[1].re;
              E.im = MAT[0].im / MAT[1].re;
            }
          } else if (MAT[1].re == 0.0) {
            if (MAT[0].re == 0.0) {
              E.re = MAT[0].im / MAT[1].im;
              E.im = 0.0;
            } else if (MAT[0].im == 0.0) {
              E.re = 0.0;
              E.im = -(MAT[0].re / MAT[1].im);
            } else {
              E.re = MAT[0].im / MAT[1].im;
              E.im = -(MAT[0].re / MAT[1].im);
            }
          } else if (c_B_re_tmp > d_B_re_tmp) {
            B_re_tmp = MAT[1].im / MAT[1].re;
            re = MAT[1].re + B_re_tmp * MAT[1].im;
            E.re = (MAT[0].re + B_re_tmp * MAT[0].im) / re;
            E.im = (MAT[0].im - B_re_tmp * MAT[0].re) / re;
          } else if (d_B_re_tmp == c_B_re_tmp) {
            if (MAT[1].re > 0.0) {
              re = 0.5;
            } else {
              re = -0.5;
            }

            if (MAT[1].im > 0.0) {
              B_re_tmp = 0.5;
            } else {
              B_re_tmp = -0.5;
            }

            E.re = (MAT[0].re * re + MAT[0].im * B_re_tmp) / c_B_re_tmp;
            E.im = (MAT[0].im * re - MAT[0].re * B_re_tmp) / c_B_re_tmp;
          } else {
            B_re_tmp = MAT[1].re / MAT[1].im;
            re = MAT[1].im + B_re_tmp * MAT[1].re;
            E.re = (B_re_tmp * MAT[0].re + MAT[0].im) / re;
            E.im = (B_re_tmp * MAT[0].im - MAT[0].re) / re;
          }

          im = (E.re * MAT[3].re - E.im * MAT[3].im) - MAT[2].re;
          re = (E.re * MAT[3].im + E.im * MAT[3].re) - MAT[2].im;
          if (re == 0.0) {
            t.re = 1.0 / im;
            t.im = 0.0;
          } else if (im == 0.0) {
            t.re = 0.0;
            t.im = -(1.0 / re);
          } else {
            brm = fabs(im);
            bim = fabs(re);
            if (brm > bim) {
              B_re_tmp = re / im;
              re = im + B_re_tmp * re;
              t.re = 1.0 / re;
              t.im = (0.0 - B_re_tmp) / re;
            } else if (bim == brm) {
              if (im > 0.0) {
                im = 0.5;
              } else {
                im = -0.5;
              }

              t.re = im / brm;
              if (re > 0.0) {
                re = 0.5;
              } else {
                re = -0.5;
              }

              t.im = (0.0 - re) / brm;
            } else {
              B_re_tmp = im / re;
              re += B_re_tmp * im;
              t.re = B_re_tmp / re;
              t.im = -1.0 / re;
            }
          }

          if (MAT[1].im == 0.0) {
            if (MAT[3].im == 0.0) {
              b_B_re_tmp = MAT[3].re / MAT[1].re;
              im = 0.0;
            } else if (MAT[3].re == 0.0) {
              b_B_re_tmp = 0.0;
              im = MAT[3].im / MAT[1].re;
            } else {
              b_B_re_tmp = MAT[3].re / MAT[1].re;
              im = MAT[3].im / MAT[1].re;
            }
          } else if (MAT[1].re == 0.0) {
            if (MAT[3].re == 0.0) {
              b_B_re_tmp = MAT[3].im / MAT[1].im;
              im = 0.0;
            } else if (MAT[3].im == 0.0) {
              b_B_re_tmp = 0.0;
              im = -(MAT[3].re / MAT[1].im);
            } else {
              b_B_re_tmp = MAT[3].im / MAT[1].im;
              im = -(MAT[3].re / MAT[1].im);
            }
          } else if (c_B_re_tmp > d_B_re_tmp) {
            B_re_tmp = MAT[1].im / MAT[1].re;
            re = MAT[1].re + B_re_tmp * MAT[1].im;
            b_B_re_tmp = (MAT[3].re + B_re_tmp * MAT[3].im) / re;
            im = (MAT[3].im - B_re_tmp * MAT[3].re) / re;
          } else if (d_B_re_tmp == c_B_re_tmp) {
            if (MAT[1].re > 0.0) {
              re = 0.5;
            } else {
              re = -0.5;
            }

            if (MAT[1].im > 0.0) {
              B_re_tmp = 0.5;
            } else {
              B_re_tmp = -0.5;
            }

            b_B_re_tmp = (MAT[3].re * re + MAT[3].im * B_re_tmp) / c_B_re_tmp;
            im = (MAT[3].im * re - MAT[3].re * B_re_tmp) / c_B_re_tmp;
          } else {
            B_re_tmp = MAT[1].re / MAT[1].im;
            re = MAT[1].im + B_re_tmp * MAT[1].re;
            b_B_re_tmp = (B_re_tmp * MAT[3].re + MAT[3].im) / re;
            im = (B_re_tmp * MAT[3].im - MAT[3].re) / re;
          }

          ps[0].re = b_B_re_tmp * t.re - im * t.im;
          ps[0].im = b_B_re_tmp * t.im + im * t.re;
          ps[1].re = -t.re;
          ps[1].im = -t.im;
          if (MAT[1].im == 0.0) {
            if (-MAT[2].im == 0.0) {
              b_B_re_tmp = -MAT[2].re / MAT[1].re;
              im = 0.0;
            } else if (-MAT[2].re == 0.0) {
              b_B_re_tmp = 0.0;
              im = -MAT[2].im / MAT[1].re;
            } else {
              b_B_re_tmp = -MAT[2].re / MAT[1].re;
              im = -MAT[2].im / MAT[1].re;
            }
          } else if (MAT[1].re == 0.0) {
            if (-MAT[2].re == 0.0) {
              b_B_re_tmp = -MAT[2].im / MAT[1].im;
              im = 0.0;
            } else if (-MAT[2].im == 0.0) {
              b_B_re_tmp = 0.0;
              im = -(-MAT[2].re / MAT[1].im);
            } else {
              b_B_re_tmp = -MAT[2].im / MAT[1].im;
              im = -(-MAT[2].re / MAT[1].im);
            }
          } else if (c_B_re_tmp > d_B_re_tmp) {
            B_re_tmp = MAT[1].im / MAT[1].re;
            re = MAT[1].re + B_re_tmp * MAT[1].im;
            b_B_re_tmp = (-MAT[2].re + B_re_tmp * -MAT[2].im) / re;
            im = (-MAT[2].im - B_re_tmp * -MAT[2].re) / re;
          } else if (d_B_re_tmp == c_B_re_tmp) {
            if (MAT[1].re > 0.0) {
              re = 0.5;
            } else {
              re = -0.5;
            }

            if (MAT[1].im > 0.0) {
              B_re_tmp = 0.5;
            } else {
              B_re_tmp = -0.5;
            }

            b_B_re_tmp = (-MAT[2].re * re + -MAT[2].im * B_re_tmp) / c_B_re_tmp;
            im = (-MAT[2].im * re - -MAT[2].re * B_re_tmp) / c_B_re_tmp;
          } else {
            B_re_tmp = MAT[1].re / MAT[1].im;
            re = MAT[1].im + B_re_tmp * MAT[1].re;
            b_B_re_tmp = (B_re_tmp * -MAT[2].re + -MAT[2].im) / re;
            im = (B_re_tmp * -MAT[2].im - (-MAT[2].re)) / re;
          }

          ps[2].re = b_B_re_tmp * t.re - im * t.im;
          ps[2].im = b_B_re_tmp * t.im + im * t.re;
          ps[3].re = E.re * t.re - E.im * t.im;
          ps[3].im = E.re * t.im + E.im * t.re;
        } else {
          if (MAT[0].im == 0.0) {
            if (MAT[1].im == 0.0) {
              E.re = MAT[1].re / MAT[0].re;
              E.im = 0.0;
            } else if (MAT[1].re == 0.0) {
              E.re = 0.0;
              E.im = MAT[1].im / MAT[0].re;
            } else {
              E.re = MAT[1].re / MAT[0].re;
              E.im = MAT[1].im / MAT[0].re;
            }
          } else if (MAT[0].re == 0.0) {
            if (MAT[1].re == 0.0) {
              E.re = MAT[1].im / MAT[0].im;
              E.im = 0.0;
            } else if (MAT[1].im == 0.0) {
              E.re = 0.0;
              E.im = -(MAT[1].re / MAT[0].im);
            } else {
              E.re = MAT[1].im / MAT[0].im;
              E.im = -(MAT[1].re / MAT[0].im);
            }
          } else if (brm > bim) {
            B_re_tmp = MAT[0].im / MAT[0].re;
            re = MAT[0].re + B_re_tmp * MAT[0].im;
            E.re = (MAT[1].re + B_re_tmp * MAT[1].im) / re;
            E.im = (MAT[1].im - B_re_tmp * MAT[1].re) / re;
          } else if (bim == brm) {
            if (MAT[0].re > 0.0) {
              re = 0.5;
            } else {
              re = -0.5;
            }

            if (MAT[0].im > 0.0) {
              B_re_tmp = 0.5;
            } else {
              B_re_tmp = -0.5;
            }

            E.re = (MAT[1].re * re + MAT[1].im * B_re_tmp) / brm;
            E.im = (MAT[1].im * re - MAT[1].re * B_re_tmp) / brm;
          } else {
            B_re_tmp = MAT[0].re / MAT[0].im;
            re = MAT[0].im + B_re_tmp * MAT[0].re;
            E.re = (B_re_tmp * MAT[1].re + MAT[1].im) / re;
            E.im = (B_re_tmp * MAT[1].im - MAT[1].re) / re;
          }

          im = MAT[3].re - (E.re * MAT[2].re - E.im * MAT[2].im);
          re = MAT[3].im - (E.re * MAT[2].im + E.im * MAT[2].re);
          if (re == 0.0) {
            t.re = 1.0 / im;
            t.im = 0.0;
          } else if (im == 0.0) {
            t.re = 0.0;
            t.im = -(1.0 / re);
          } else {
            c_B_re_tmp = fabs(im);
            d_B_re_tmp = fabs(re);
            if (c_B_re_tmp > d_B_re_tmp) {
              B_re_tmp = re / im;
              re = im + B_re_tmp * re;
              t.re = 1.0 / re;
              t.im = (0.0 - B_re_tmp) / re;
            } else if (d_B_re_tmp == c_B_re_tmp) {
              if (im > 0.0) {
                im = 0.5;
              } else {
                im = -0.5;
              }

              t.re = im / c_B_re_tmp;
              if (re > 0.0) {
                re = 0.5;
              } else {
                re = -0.5;
              }

              t.im = (0.0 - re) / c_B_re_tmp;
            } else {
              B_re_tmp = im / re;
              re += B_re_tmp * im;
              t.re = B_re_tmp / re;
              t.im = -1.0 / re;
            }
          }

          if (MAT[0].im == 0.0) {
            if (MAT[3].im == 0.0) {
              b_B_re_tmp = MAT[3].re / MAT[0].re;
              im = 0.0;
            } else if (MAT[3].re == 0.0) {
              b_B_re_tmp = 0.0;
              im = MAT[3].im / MAT[0].re;
            } else {
              b_B_re_tmp = MAT[3].re / MAT[0].re;
              im = MAT[3].im / MAT[0].re;
            }
          } else if (MAT[0].re == 0.0) {
            if (MAT[3].re == 0.0) {
              b_B_re_tmp = MAT[3].im / MAT[0].im;
              im = 0.0;
            } else if (MAT[3].im == 0.0) {
              b_B_re_tmp = 0.0;
              im = -(MAT[3].re / MAT[0].im);
            } else {
              b_B_re_tmp = MAT[3].im / MAT[0].im;
              im = -(MAT[3].re / MAT[0].im);
            }
          } else if (brm > bim) {
            B_re_tmp = MAT[0].im / MAT[0].re;
            re = MAT[0].re + B_re_tmp * MAT[0].im;
            b_B_re_tmp = (MAT[3].re + B_re_tmp * MAT[3].im) / re;
            im = (MAT[3].im - B_re_tmp * MAT[3].re) / re;
          } else if (bim == brm) {
            if (MAT[0].re > 0.0) {
              re = 0.5;
            } else {
              re = -0.5;
            }

            if (MAT[0].im > 0.0) {
              B_re_tmp = 0.5;
            } else {
              B_re_tmp = -0.5;
            }

            b_B_re_tmp = (MAT[3].re * re + MAT[3].im * B_re_tmp) / brm;
            im = (MAT[3].im * re - MAT[3].re * B_re_tmp) / brm;
          } else {
            B_re_tmp = MAT[0].re / MAT[0].im;
            re = MAT[0].im + B_re_tmp * MAT[0].re;
            b_B_re_tmp = (B_re_tmp * MAT[3].re + MAT[3].im) / re;
            im = (B_re_tmp * MAT[3].im - MAT[3].re) / re;
          }

          ps[0].re = b_B_re_tmp * t.re - im * t.im;
          ps[0].im = b_B_re_tmp * t.im + im * t.re;
          ps[1].re = -E.re * t.re - -E.im * t.im;
          ps[1].im = -E.re * t.im + -E.im * t.re;
          if (MAT[0].im == 0.0) {
            if (-MAT[2].im == 0.0) {
              b_B_re_tmp = -MAT[2].re / MAT[0].re;
              im = 0.0;
            } else if (-MAT[2].re == 0.0) {
              b_B_re_tmp = 0.0;
              im = -MAT[2].im / MAT[0].re;
            } else {
              b_B_re_tmp = -MAT[2].re / MAT[0].re;
              im = -MAT[2].im / MAT[0].re;
            }
          } else if (MAT[0].re == 0.0) {
            if (-MAT[2].re == 0.0) {
              b_B_re_tmp = -MAT[2].im / MAT[0].im;
              im = 0.0;
            } else if (-MAT[2].im == 0.0) {
              b_B_re_tmp = 0.0;
              im = -(-MAT[2].re / MAT[0].im);
            } else {
              b_B_re_tmp = -MAT[2].im / MAT[0].im;
              im = -(-MAT[2].re / MAT[0].im);
            }
          } else if (brm > bim) {
            B_re_tmp = MAT[0].im / MAT[0].re;
            re = MAT[0].re + B_re_tmp * MAT[0].im;
            b_B_re_tmp = (-MAT[2].re + B_re_tmp * -MAT[2].im) / re;
            im = (-MAT[2].im - B_re_tmp * -MAT[2].re) / re;
          } else if (bim == brm) {
            if (MAT[0].re > 0.0) {
              re = 0.5;
            } else {
              re = -0.5;
            }

            if (MAT[0].im > 0.0) {
              B_re_tmp = 0.5;
            } else {
              B_re_tmp = -0.5;
            }

            b_B_re_tmp = (-MAT[2].re * re + -MAT[2].im * B_re_tmp) / brm;
            im = (-MAT[2].im * re - -MAT[2].re * B_re_tmp) / brm;
          } else {
            B_re_tmp = MAT[0].re / MAT[0].im;
            re = MAT[0].im + B_re_tmp * MAT[0].re;
            b_B_re_tmp = (B_re_tmp * -MAT[2].re + -MAT[2].im) / re;
            im = (B_re_tmp * -MAT[2].im - (-MAT[2].re)) / re;
          }

          ps[2].re = b_B_re_tmp * t.re - im * t.im;
          ps[2].im = b_B_re_tmp * t.im + im * t.re;
          ps[3] = t;
        }

        B_size[0] = 2;
        B_size[1] = 2;
        memcpy(&B[0], &ps[0], 4U * sizeof(creal_T));
      } else {
        /* 'gkmPWM:669' else */
        /* 'gkmPWM:670' B=MAT^-1; */
        b_MAT_data.data = &MAT_data[0];
        b_MAT_data.size = &MAT_size[0];
        b_MAT_data.allocatedSize = 16;
        b_MAT_data.numDimensions = 2;
        b_MAT_data.canFreeData = false;
        mpower(&b_MAT_data, c_x);
        c_x_data = c_x->data;
        B_size[0] = c_x->size[0];
        B_size[1] = c_x->size[1];
        loop_ub = c_x->size[0] * c_x->size[1];
        for (i = 0; i < loop_ub; i++) {
          B[i].re = c_x_data[i];
          B[i].im = 0.0;
        }
      }

      /* 'gkmPWM:672' ps = B*Y; */
      loop_ub = Posvec_size[1];
      for (i = 0; i < loop_ub; i++) {
        b_y[i].re = y[(int)Posvec_data[i] - 1];
        b_y[i].im = 0.0;
      }

      loop_ub = B_size[0];
      for (i = 0; i < loop_ub; i++) {
        ps[i].re = 0.0;
        ps[i].im = 0.0;
        last = B_size[1];
        for (i1 = 0; i1 < last; i1++) {
          k = i + B_size[0] * i1;
          B_re_tmp = B[k].re;
          b_B_re_tmp = b_y[i1].im;
          c_B_re_tmp = B[k].im;
          d_B_re_tmp = b_y[i1].re;
          ps[i].re += B_re_tmp * d_B_re_tmp - c_B_re_tmp * b_B_re_tmp;
          ps[i].im += B_re_tmp * b_B_re_tmp + c_B_re_tmp * d_B_re_tmp;
        }
      }

      /* 'gkmPWM:673' p = real(ps+(1-sum(ps))*B/sum(sum(B))*ones(2,1)); */
      if (B_size[0] == 0) {
        E.re = 0.0;
        E.im = 0.0;
      } else {
        E = ps[0];
        for (k = 2; k <= loop_ub; k++) {
          E.re += ps[k - 1].re;
          E.im += ps[k - 1].im;
        }
      }

      t.re = 1.0 - E.re;
      t.im = 0.0 - E.im;
      f_sum(B, B_size, b_y, x_size);
      last = x_size[1];
      if (x_size[1] == 0) {
        E.re = 0.0;
        E.im = 0.0;
      } else {
        E = b_y[0];
        for (k = 2; k <= last; k++) {
          E.re += b_y[k - 1].re;
          E.im += b_y[k - 1].im;
        }
      }

      k = B_size[0];
      loop_ub = B_size[0] * B_size[1];
      for (i = 0; i < loop_ub; i++) {
        re = B[i].re;
        im = B[i].im;
        b_B_re_tmp = t.re * re - t.im * im;
        im = t.re * im + t.im * re;
        if (E.im == 0.0) {
          if (im == 0.0) {
            t_data[i].re = b_B_re_tmp / E.re;
            t_data[i].im = 0.0;
          } else if (b_B_re_tmp == 0.0) {
            t_data[i].re = 0.0;
            t_data[i].im = im / E.re;
          } else {
            t_data[i].re = b_B_re_tmp / E.re;
            t_data[i].im = im / E.re;
          }
        } else if (E.re == 0.0) {
          if (b_B_re_tmp == 0.0) {
            t_data[i].re = im / E.im;
            t_data[i].im = 0.0;
          } else if (im == 0.0) {
            t_data[i].re = 0.0;
            t_data[i].im = -(b_B_re_tmp / E.im);
          } else {
            t_data[i].re = im / E.im;
            t_data[i].im = -(b_B_re_tmp / E.im);
          }
        } else {
          brm = fabs(E.re);
          bim = fabs(E.im);
          if (brm > bim) {
            B_re_tmp = E.im / E.re;
            re = E.re + B_re_tmp * E.im;
            t_data[i].re = (b_B_re_tmp + B_re_tmp * im) / re;
            t_data[i].im = (im - B_re_tmp * b_B_re_tmp) / re;
          } else if (bim == brm) {
            if (E.re > 0.0) {
              re = 0.5;
            } else {
              re = -0.5;
            }

            if (E.im > 0.0) {
              B_re_tmp = 0.5;
            } else {
              B_re_tmp = -0.5;
            }

            t_data[i].re = (b_B_re_tmp * re + im * B_re_tmp) / brm;
            t_data[i].im = (im * re - b_B_re_tmp * B_re_tmp) / brm;
          } else {
            B_re_tmp = E.re / E.im;
            re = E.im + B_re_tmp * E.re;
            t_data[i].re = (B_re_tmp * b_B_re_tmp + im) / re;
            t_data[i].im = (B_re_tmp * im - b_B_re_tmp) / re;
          }
        }
      }

      for (i = 0; i < 2; i++) {
        dcv1[i].re = 1.0;
        dcv1[i].im = 0.0;
      }

      p_size = B_size[0];
      for (i = 0; i < k; i++) {
        re = 0.0;
        for (i1 = 0; i1 < 2; i1++) {
          last = i + k * i1;
          re += t_data[last].re - t_data[last].im * dcv1[i1].im;
        }

        p_data[i] = ps[i].re + re;
      }

      /* solution */
      /* Checks if solution is permitted.  If so, makes sure that it creates a smaller error than other cases */
      /* 'gkmPWM:675' if min(p) >= 0 */
      last = p_size;
      if (p_size <= 2) {
        if (p_size == 1) {
          im = p_data[0];
        } else if (p_data[0] > p_data[p_size - 1]) {
          im = p_data[p_size - 1];
        } else {
          im = p_data[0];
        }
      } else {
        im = p_data[0];
        for (k = 2; k <= last; k++) {
          re = p_data[k - 1];
          if (im > re) {
            im = re;
          }
        }
      }

      if (im >= 0.0) {
        /* 'gkmPWM:676' vec = b-A(:,Posvec)*p; */
        loop_ub = A->size[0];
        i = c_x->size[0] * c_x->size[1];
        c_x->size[0] = A->size[0];
        last = Posvec_size[1];
        c_x->size[1] = Posvec_size[1];
        emxEnsureCapacity_real_T(c_x, i);
        c_x_data = c_x->data;
        for (i = 0; i < last; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            c_x_data[i1 + c_x->size[0] * i] = A_data[i1 + A->size[0] * ((int)
              Posvec_data[i] - 1)];
          }
        }

        loop_ub = A->size[0];
        if ((A->size[0] == 0) || (Posvec_size[1] == 0) || (p_size == 0)) {
          i = vec->size[0];
          vec->size[0] = A->size[0];
          emxEnsureCapacity_real_T(vec, i);
          vec_data = vec->data;
          for (i = 0; i < loop_ub; i++) {
            vec_data[i] = 0.0;
          }
        } else {
          i = vec->size[0];
          vec->size[0] = A->size[0];
          emxEnsureCapacity_real_T(vec, i);
          vec_data = vec->data;
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                      A->size[0], (blasint)1, (blasint)Posvec_size[1], 1.0,
                      &c_x_data[0], (blasint)A->size[0], &p_data[0], (blasint)
                      p_size, 0.0, &vec_data[0], (blasint)A->size[0]);
        }

        loop_ub = b->size[0];
        if (b->size[0] == vec->size[0]) {
          i = vec->size[0];
          vec->size[0] = b->size[0];
          emxEnsureCapacity_real_T(vec, i);
          vec_data = vec->data;
          for (i = 0; i < loop_ub; i++) {
            vec_data[i] = b_data[i] - vec_data[i];
          }
        } else {
          minus(vec, b);
          vec_data = vec->data;
        }

        /* 'gkmPWM:677' if I == 0 */
        if (b_I == 0) {
          /* 'gkmPWM:678' I = 1; */
          b_I = 1;

          /* 'gkmPWM:679' P = zeros(4,1); */
          *P_size = 4;
          for (i = 0; i < 4; i++) {
            P_data[i] = 0.0;
          }

          /* 'gkmPWM:680' P(Posvec) = real(p); */
          loop_ub = Posvec_size[1];
          for (i = 0; i < loop_ub; i++) {
            b_tmp_data[i] = (int)Posvec_data[i];
          }

          for (i = 0; i < loop_ub; i++) {
            P_data[b_tmp_data[i] - 1] = p_data[i];
          }

          /* 'gkmPWM:681' e = vec'*vec-reg*M*p'*p; */
          E.re = reg * M.re;
          if (vec->size[0] < 1) {
            im = 0.0;
          } else {
            im = cblas_ddot((blasint)vec->size[0], &vec_data[0], (blasint)1,
                            &vec_data[0], (blasint)1);
          }

          re = 0.0;
          loop_ub = p_size;
          for (i = 0; i < loop_ub; i++) {
            re += p_data[i] * E.re * p_data[i];
          }

          e_re = im - re;
        } else {
          /* 'gkmPWM:682' else */
          /* 'gkmPWM:683' E = vec'*vec-reg*M*p'*p; */
          E.re = reg * M.re;
          if (vec->size[0] < 1) {
            im = 0.0;
          } else {
            im = cblas_ddot((blasint)vec->size[0], &vec_data[0], (blasint)1,
                            &vec_data[0], (blasint)1);
          }

          re = 0.0;
          loop_ub = p_size;
          for (i = 0; i < loop_ub; i++) {
            re += p_data[i] * E.re * p_data[i];
          }

          E.re = im - re;

          /* 'gkmPWM:684' if E < e */
          if (E.re < e_re) {
            /* 'gkmPWM:685' e = E; */
            e_re = E.re;

            /* 'gkmPWM:686' P = zeros(4,1); */
            *P_size = 4;
            for (i = 0; i < 4; i++) {
              P_data[i] = 0.0;
            }

            /* 'gkmPWM:687' P(Posvec) = real(p); */
            loop_ub = Posvec_size[1];
            for (i = 0; i < loop_ub; i++) {
              b_tmp_data[i] = (int)Posvec_data[i];
            }

            for (i = 0; i < loop_ub; i++) {
              P_data[b_tmp_data[i] - 1] = p_data[i];
            }
          }
        }
      }
    }

    emxFree_real_T(&c_x);

    /* Checks to see if one non-zero case is better than the other cases */
    /* 'gkmPWM:693' if I == 0 */
    if (b_I == 0) {
      /* 'gkmPWM:694' P = zeros(4,1); */
      *P_size = 4;
      for (i = 0; i < 4; i++) {
        P_data[i] = 0.0;
      }

      /* 'gkmPWM:695' P(a) = 1; */
      P_data[iindx - 1] = 1.0;
    } else {
      /* 'gkmPWM:696' else */
      /* 'gkmPWM:697' vec = b-A(:,a); */
      if (b->size[0] == A->size[0]) {
        i = vec->size[0];
        vec->size[0] = b->size[0];
        emxEnsureCapacity_real_T(vec, i);
        vec_data = vec->data;
        loop_ub = b->size[0];
        for (i = 0; i < loop_ub; i++) {
          vec_data[i] = b_data[i] - A_data[i + A->size[0] * (iindx - 1)];
        }
      } else {
        y_binary_expand_op(vec, b, A, iindx);
        vec_data = vec->data;
      }

      /* 'gkmPWM:698' E = vec'*vec-reg*M; */
      /* 'gkmPWM:699' if E < e */
      if (vec->size[0] < 1) {
        im = 0.0;
      } else {
        im = cblas_ddot((blasint)vec->size[0], &vec_data[0], (blasint)1,
                        &vec_data[0], (blasint)1);
      }

      if (im - reg * M.re < e_re) {
        /* 'gkmPWM:700' e = E; */
        /* 'gkmPWM:701' P = zeros(4,1); */
        *P_size = 4;
        for (i = 0; i < 4; i++) {
          P_data[i] = 0.0;
        }

        /* 'gkmPWM:702' P(Posvec) = real(p); */
        last = Posvec_size[1];
        loop_ub = Posvec_size[1];
        for (i = 0; i < loop_ub; i++) {
          b_tmp_data[i] = (int)Posvec_data[i];
        }

        for (i = 0; i < last; i++) {
          P_data[b_tmp_data[i] - 1] = p_data[i];
        }
      }
    }

    emxFree_real_T(&vec);
  } else {
    /* 'gkmPWM:705' else */
    /* 'gkmPWM:706' P = real(p); */
    *P_size = 4;
    for (i = 0; i < 4; i++) {
      P_data[i] = p_data[i];
    }
  }

  emxFree_real_T(&b);

  /* 'gkmPWM:708' kweig = A*(P-PWM(l_svm,:)'); */
  for (i = 0; i < 4; i++) {
    y[i] = P_data[i] - PWM_data[((int)l_svm + PWM->size[0] * i) - 1];
  }

  d_mtimes(A, y, kweig);

  /* 'gkmPWM:709' P = real(P); */
  /* 'gkmPWM:710' kweig = real(kweig); */
  /* 'gkmPWM:711' p = zeros(4,1); */
  emxFree_real_T(&A);
}

/*
 * function [PWM, scorevec, C, r, R, E, Rd] = gkmPWM_lagrange(kweig,negmat,PWM,negvec,n,rcorr,reg,l_svm,k_svm,RC,rc,diffc,indc,xc,rcnum)
 */
static void gkmPWM_lagrange(const emxArray_real_T *kweig, const double negmat[16],
  const emxArray_cell_wrap_2 *PWM, const emxArray_real_T *negvec, double n,
  double rcorr, double reg, double l_svm, double k_svm, double RC, const
  emxArray_real_T *rc, const emxArray_real_T *diffc, const emxArray_real_T *indc,
  const emxArray_real_T *xc, double rcnum, emxArray_cell_wrap_2 *b_PWM,
  emxArray_real_T *scorevec, emxArray_real_T *C, double *r, emxArray_real_T *R,
  emxArray_real_T *E, emxArray_real_T *Rd)
{
  cell_wrap_14 *loc_data;
  cell_wrap_14 *poscell_data;
  const cell_wrap_2 *PWM_data;
  cell_wrap_2 *b_PWM_data;
  cell_wrap_2 *new_PWM_data;
  emxArray_boolean_T *b_ct;
  emxArray_boolean_T *d_ct;
  emxArray_cell_wrap_14 *b_new_loc;
  emxArray_cell_wrap_14 *loc;
  emxArray_cell_wrap_14 *new_loc;
  emxArray_cell_wrap_14 *poscell;
  emxArray_cell_wrap_2 *b_new_PWM;
  emxArray_cell_wrap_2 *c_PWM;
  emxArray_cell_wrap_2 *new_PWM;
  emxArray_int32_T *c_iidx;
  emxArray_int32_T *iidx;
  emxArray_real_T *CT;
  emxArray_real_T *b_a;
  emxArray_real_T *b_iidx;
  emxArray_real_T *b_loc;
  emxArray_real_T *b_mat;
  emxArray_real_T *c_ct;
  emxArray_real_T *ct;
  emxArray_real_T *d_PWM;
  emxArray_real_T *diffC;
  emxArray_real_T *f;
  emxArray_real_T *kmat;
  emxArray_real_T *lenvec;
  emxArray_real_T *mat;
  emxArray_real_T *ord;
  emxArray_real_T *res;
  emxArray_real_T *temp;
  emxArray_real_T *vec;
  double b_varargin_1[9];
  double b[4];
  double tmp_data[4];
  double varargin_1[4];
  double b_lenvec[2];
  const double *kweig_data;
  const double *negvec_data;
  const double *rc_data;
  double M;
  double S;
  double c;
  double scount;
  double *CT_data;
  double *C_data;
  double *b_iidx_data;
  double *ct_data;
  double *diffC_data;
  double *f_data;
  double *kmat_data;
  double *lenvec_data;
  double *mat_data;
  double *ord_data;
  double *res_data;
  double *scorevec_data;
  double *temp_data;
  unsigned int a;
  int acount;
  int b_loop_ub;
  int c_loop_ub;
  int d_loop_ub;
  int e_loop_ub;
  int exitg1;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int j;
  int jj;
  int lcomb;
  int loop_ub;
  int nx;
  int u0;
  int varargin_2;
  int *iidx_data;
  bool x[9];
  bool guard1 = false;
  bool *b_ct_data;
  rc_data = rc->data;
  negvec_data = negvec->data;
  PWM_data = PWM->data;
  kweig_data = kweig->data;
  emxInit_cell_wrap_2(&c_PWM);
  i = c_PWM->size[0];
  c_PWM->size[0] = PWM->size[0];
  emxEnsureCapacity_cell_wrap_2(c_PWM, i);
  b_PWM_data = c_PWM->data;
  loop_ub = PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    emxCopyStruct_cell_wrap_2(&b_PWM_data[i], &PWM_data[i]);
  }

  emxInit_real_T(&diffC, 2);

  /* Note: This code is rather messy.  I block commmented to the best of my ability, so hopefully this makes enough sense.  If something seems non-trivial, then I probably found a mathematical trick to speed up computation time (in particular dynamic programming). */
  /* 'gkmPWM:354' GC = PWM{1}(1,:); */
  /* 'gkmPWM:355' lcomb = length(diffc); */
  lcomb = diffc->size[0];

  /* 'gkmPWM:356' diffC = zeros(lcomb,l_svm); */
  i = diffC->size[0] * diffC->size[1];
  diffC->size[0] = diffc->size[0];
  i1 = (int)l_svm;
  diffC->size[1] = (int)l_svm;
  emxEnsureCapacity_real_T(diffC, i);
  diffC_data = diffC->data;
  loop_ub = diffc->size[0] * (int)l_svm;
  for (i = 0; i < loop_ub; i++) {
    diffC_data[i] = 0.0;
  }

  /* 'gkmPWM:357' for i = 1:l_svm */
  emxInit_real_T(&ct, 2);
  emxInit_real_T(&f, 1);
  emxInit_real_T(&CT, 2);
  emxInit_int32_T(&iidx, 1);
  emxInit_boolean_T(&b_ct, 2);
  emxInit_real_T(&b_iidx, 1);
  emxInit_real_T(&c_ct, 2);
  for (acount = 0; acount < i1; acount++) {
    /* 'gkmPWM:358' ct = rc+i-1; */
    i = ct->size[0] * ct->size[1];
    ct->size[0] = rc->size[0];
    ct->size[1] = rc->size[1];
    emxEnsureCapacity_real_T(ct, i);
    ct_data = ct->data;
    loop_ub = rc->size[0] * rc->size[1];
    for (i = 0; i < loop_ub; i++) {
      ct_data[i] = (rc_data[i] + ((double)acount + 1.0)) - 1.0;
    }

    /* 'gkmPWM:359' f = find(sum(ct==l_svm,2)); */
    i = b_ct->size[0] * b_ct->size[1];
    b_ct->size[0] = ct->size[0];
    b_ct->size[1] = ct->size[1];
    emxEnsureCapacity_boolean_T(b_ct, i);
    b_ct_data = b_ct->data;
    loop_ub = ct->size[0] * ct->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_ct_data[i] = (ct_data[i] == l_svm);
    }

    combineVectorElements(b_ct, iidx);
    iidx_data = iidx->data;
    i = b_iidx->size[0];
    b_iidx->size[0] = iidx->size[0];
    emxEnsureCapacity_real_T(b_iidx, i);
    b_iidx_data = b_iidx->data;
    loop_ub = iidx->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_iidx_data[i] = iidx_data[i];
    }

    c_eml_find(b_iidx, iidx);
    iidx_data = iidx->data;
    i = f->size[0];
    f->size[0] = iidx->size[0];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = iidx->size[0];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = iidx_data[i];
    }

    /* 'gkmPWM:360' ct = ct(f,:); */
    nx = ct->size[1] - 1;
    i = c_ct->size[0] * c_ct->size[1];
    c_ct->size[0] = f->size[0];
    c_ct->size[1] = ct->size[1];
    emxEnsureCapacity_real_T(c_ct, i);
    mat_data = c_ct->data;
    for (i = 0; i <= nx; i++) {
      loop_ub = f->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        mat_data[i2 + c_ct->size[0] * i] = ct_data[((int)f_data[i2] + ct->size[0]
          * i) - 1];
      }
    }

    i = ct->size[0] * ct->size[1];
    ct->size[0] = c_ct->size[0];
    ct->size[1] = c_ct->size[1];
    emxEnsureCapacity_real_T(ct, i);
    ct_data = ct->data;
    loop_ub = c_ct->size[0] * c_ct->size[1];
    for (i = 0; i < loop_ub; i++) {
      ct_data[i] = mat_data[i];
    }

    /* 'gkmPWM:361' CT = zeros(length(ct),k_svm-1); */
    if ((ct->size[0] == 0) || (ct->size[1] == 0)) {
      nx = 0;
    } else {
      u0 = ct->size[0];
      nx = ct->size[1];
      if (u0 >= nx) {
        nx = u0;
      }
    }

    i = CT->size[0] * CT->size[1];
    CT->size[0] = nx;
    CT->size[1] = (int)(k_svm - 1.0);
    emxEnsureCapacity_real_T(CT, i);
    CT_data = CT->data;
    loop_ub = nx * (int)(k_svm - 1.0);
    for (i = 0; i < loop_ub; i++) {
      CT_data[i] = 0.0;
    }

    /* 'gkmPWM:362' for j = 1:numel(ct)/k_svm */
    i = (int)((double)(ct->size[0] * ct->size[1]) / k_svm);
    for (j = 0; j < i; j++) {
      /* 'gkmPWM:363' a = 1; */
      a = 1U;

      /* 'gkmPWM:364' for jj = 1:k_svm */
      i2 = (int)k_svm;
      for (jj = 0; jj < i2; jj++) {
        /* 'gkmPWM:365' if ct(j,jj) ~= l_svm */
        S = ct_data[j + ct->size[0] * jj];
        if (S != l_svm) {
          /* 'gkmPWM:366' CT(j,a) = ct(j,jj); */
          CT_data[j + CT->size[0] * ((int)a - 1)] = S;

          /* 'gkmPWM:367' a = a+1; */
          a++;
        }
      }
    }

    /* 'gkmPWM:371' for j = 2:length(f) */
    i = f->size[0];
    for (j = 0; j <= i - 2; j++) {
      /* 'gkmPWM:372' a = 1; */
      /* 'gkmPWM:373' while CT(j,a)==CT(j-1,a) */
      for (a = 1U; CT_data[(j + CT->size[0] * ((int)a - 1)) + 1] == CT_data[j +
           CT->size[0] * ((int)a - 1)]; a++) {
        /* 'gkmPWM:374' a = a+1; */
      }

      /* 'gkmPWM:376' if a < 2 */
      if ((int)a < 2) {
        /* 'gkmPWM:377' a = 2; */
        a = 2U;
      }

      /* 'gkmPWM:379' diffC(f(j),i)=a; */
      diffC_data[((int)f_data[j + 1] + diffC->size[0] * acount) - 1] = a;
    }

    /* 'gkmPWM:381' diffC(f(1),i)=2; */
    diffC_data[((int)f_data[0] + diffC->size[0] * acount) - 1] = 2.0;
  }

  emxFree_boolean_T(&b_ct);

  /* 'gkmPWM:383' m = length(PWM); */
  varargin_2 = PWM->size[0] - 1;

  /* 'gkmPWM:384' scorevec = zeros(1,n); */
  i = scorevec->size[0] * scorevec->size[1];
  scorevec->size[0] = 1;
  nx = (int)n;
  scorevec->size[1] = (int)n;
  emxEnsureCapacity_real_T(scorevec, i);
  scorevec_data = scorevec->data;
  for (i = 0; i < nx; i++) {
    scorevec_data[i] = 0.0;
  }

  emxInit_real_T(&lenvec, 1);
  emxInit_cell_wrap_14(&loc);

  /* 'gkmPWM:385' lenvec = zeros(m,1); */
  /* 'gkmPWM:386' loc = cell(m, 1); */
  /* 'gkmPWM:387' for i = 1:m */
  i = PWM->size[0];
  i1 = lenvec->size[0];
  lenvec->size[0] = PWM->size[0];
  emxEnsureCapacity_real_T(lenvec, i1);
  lenvec_data = lenvec->data;
  i1 = loc->size[0];
  loc->size[0] = PWM->size[0];
  emxEnsureCapacity_cell_wrap_14(loc, i1);
  loc_data = loc->data;
  for (acount = 0; acount < i; acount++) {
    /* 'gkmPWM:388' lenvec(i) = length(PWM{i})-l_svm*2+2; */
    u0 = PWM_data[acount].f1->size[0];
    if (u0 < 4) {
      u0 = 4;
    }

    if (PWM_data[acount].f1->size[0] == 0) {
      u0 = 0;
    }

    lenvec_data[acount] = ((double)u0 - l_svm * 2.0) + 2.0;

    /* 'gkmPWM:389' loc{i} = zeros(length(PWM{i}), 1); */
    u0 = PWM_data[acount].f1->size[0];
    if (u0 < 4) {
      u0 = 4;
    }

    if (PWM_data[acount].f1->size[0] == 0) {
      nx = 0;
    } else {
      nx = u0;
    }

    i1 = loc_data[acount].f1->size[0];
    loc_data[acount].f1->size[0] = nx;
    emxEnsureCapacity_real_T(loc_data[acount].f1, i1);
    for (i1 = 0; i1 < nx; i1++) {
      loc_data[acount].f1->data[i1] = 0.0;
    }

    /* 'gkmPWM:390' loc{i}(l_svm:lenvec(i)+l_svm-1) = 1; */
    S = (lenvec_data[acount] + l_svm) - 1.0;
    if (l_svm > S) {
      i1 = -1;
      i2 = 0;
    } else {
      i1 = (int)l_svm - 2;
      i2 = (int)S;
    }

    loop_ub = (i2 - i1) - 1;
    for (i2 = 0; i2 < loop_ub; i2++) {
      loc_data[acount].f1->data[(i1 + i2) + 1] = 1.0;
    }
  }

  emxInit_real_T(&kmat, 2);

  /* 'gkmPWM:393' kmat = zeros(lcomb*4^k_svm, m); */
  M = pow(4.0, k_svm);
  i = kmat->size[0] * kmat->size[1];
  kmat->size[0] = (int)((double)diffc->size[0] * M);
  kmat->size[1] = PWM->size[0];
  emxEnsureCapacity_real_T(kmat, i);
  kmat_data = kmat->data;
  loop_ub = (int)((double)diffc->size[0] * M) * PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    kmat_data[i] = 0.0;
  }

  /* 'gkmPWM:394' KMAT = zeros(lcomb*4^k_svm, m); */
  i = CT->size[0] * CT->size[1];
  CT->size[0] = (int)((double)diffc->size[0] * M);
  CT->size[1] = PWM->size[0];
  emxEnsureCapacity_real_T(CT, i);
  CT_data = CT->data;
  loop_ub = (int)((double)diffc->size[0] * M) * PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    CT_data[i] = 0.0;
  }

  /* 'gkmPWM:396' fprintf('Mapping PWMs to gkm space\n'); */
  printf("Mapping PWMs to gkm space\n");
  fflush(stdout);

  /* 'gkmPWM:397' if RC */
  if (RC != 0.0) {
    /* 'gkmPWM:398' for i = 1:m */
    i = PWM->size[0];
    for (acount = 0; acount < i; acount++) {
      /* 'gkmPWM:399' kmat(:,i) = PWM2kmers(PWM{i},negmat,rc,diffc,indc,loc{i},xc,l_svm,k_svm,rcnum)-negvec*(lenvec(i)+l_svm-1); */
      PWM2kmers(PWM_data[acount].f1, negmat, rc, diffc, indc, loc_data[acount].
                f1, xc, l_svm, k_svm, rcnum, b_iidx);
      b_iidx_data = b_iidx->data;
      S = (lenvec_data[acount] + l_svm) - 1.0;
      if (b_iidx->size[0] == negvec->size[0]) {
        loop_ub = b_iidx->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          kmat_data[i1 + kmat->size[0] * acount] = b_iidx_data[i1] -
            negvec_data[i1] * S;
        }
      } else {
        m_binary_expand_op(kmat, acount, b_iidx, negvec, S);
        kmat_data = kmat->data;
      }

      /* map PWMs to gapped kmers */
    }
  } else {
    /* 'gkmPWM:401' else */
    /* 'gkmPWM:402' for i = 1:m */
    i = PWM->size[0];
    for (acount = 0; acount < i; acount++) {
      /* 'gkmPWM:403' kmat(:,i) = PWM2kmers_norc(PWM{i},negmat,rc,diffc,indc,loc{i},xc,l_svm,k_svm,rcnum)-negvec*(lenvec(i)+l_svm-1); */
      PWM2kmers_norc(PWM_data[acount].f1, negmat, rc, diffc, indc,
                     loc_data[acount].f1, xc, l_svm, k_svm, b_iidx);
      b_iidx_data = b_iidx->data;
      S = (lenvec_data[acount] + l_svm) - 1.0;
      if (b_iidx->size[0] == negvec->size[0]) {
        loop_ub = b_iidx->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          kmat_data[i1 + kmat->size[0] * acount] = b_iidx_data[i1] -
            negvec_data[i1] * S;
        }
      } else {
        m_binary_expand_op(kmat, acount, b_iidx, negvec, S);
        kmat_data = kmat->data;
      }

      /* map PWMs to gapped kmers */
    }
  }

  emxInit_cell_wrap_14(&poscell);

  /* the following loop creates indices for the PWM column optimize to utilize dynamic programming. */
  /* 'gkmPWM:408' poscell = cell(k_svm,1); */
  /* 'gkmPWM:409' for i = 1:k_svm */
  i = (int)k_svm;
  i1 = poscell->size[0];
  poscell->size[0] = (int)k_svm;
  emxEnsureCapacity_cell_wrap_14(poscell, i1);
  poscell_data = poscell->data;
  if (0 <= (int)k_svm - 1) {
    c = pow(4.0, k_svm - 1.0);
    scount = M;
  }

  emxInit_real_T(&temp, 2);
  emxInit_real_T(&vec, 2);
  for (acount = 0; acount < i; acount++) {
    /* 'gkmPWM:410' temp = zeros(4^(k_svm-1),1)'; */
    i1 = temp->size[0] * temp->size[1];
    temp->size[0] = 1;
    nx = (int)c;
    temp->size[1] = (int)c;
    emxEnsureCapacity_real_T(temp, i1);
    temp_data = temp->data;
    for (i1 = 0; i1 < nx; i1++) {
      temp_data[i1] = 0.0;
    }

    /* 'gkmPWM:411' vec = 1:4^(i):4^(k_svm); */
    M = pow(4.0, (double)acount + 1.0);
    if (floor(M) == M) {
      i1 = vec->size[0] * vec->size[1];
      vec->size[0] = 1;
      loop_ub = (int)floor((scount - 1.0) / M);
      vec->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(vec, i1);
      diffC_data = vec->data;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        diffC_data[i1] = M * (double)i1 + 1.0;
      }
    } else {
      eml_float_colon(M, scount, vec);
      diffC_data = vec->data;
    }

    /* 'gkmPWM:412' for ii = 1:4^(i-1) */
    i1 = (int)pow(4.0, ((double)acount + 1.0) - 1.0);
    if (0 <= i1 - 1) {
      b_loop_ub = vec->size[1];
    }

    for (j = 0; j < i1; j++) {
      /* 'gkmPWM:413' t = length(vec); */
      /* 'gkmPWM:414' temp(1+(ii-1)*t:t+(ii-1)*t) = vec+ii-1; */
      S = (((double)j + 1.0) - 1.0) * (double)vec->size[1] + 1.0;
      if (S > (double)vec->size[1] + (((double)j + 1.0) - 1.0) * (double)
          vec->size[1]) {
        i2 = 1;
      } else {
        i2 = (int)S;
      }

      for (i3 = 0; i3 < b_loop_ub; i3++) {
        temp_data[(i2 + i3) - 1] = (diffC_data[i3] + ((double)j + 1.0)) - 1.0;
      }
    }

    /* 'gkmPWM:416' poscell{i} = sort(temp)'; */
    b_sort(temp);
    temp_data = temp->data;
    i1 = poscell_data[acount].f1->size[0];
    poscell_data[acount].f1->size[0] = temp->size[1];
    emxEnsureCapacity_real_T(poscell_data[acount].f1, i1);
    loop_ub = temp->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      poscell_data[acount].f1->data[i1] = temp_data[i1];
    }
  }

  emxFree_real_T(&vec);

  /* 'gkmPWM:418' fprintf('Running Recursion\n'); */
  printf("Running Recursion\n");
  fflush(stdout);

  /* 'gkmPWM:419' acount = 0; */
  acount = 0;

  /* 'gkmPWM:420' i = 0; */
  a = 0U;

  /* 'gkmPWM:421' scount = 0; */
  scount = 0.0;

  /* 'gkmPWM:422' C = (kmat'*kmat)^(-1)*(kmat'*kweig); */
  b_mtimes(kmat, kmat, c_ct);
  if ((kmat->size[0] == 0) || (kmat->size[1] == 0) || (kweig->size[0] == 0)) {
    i = b_iidx->size[0];
    b_iidx->size[0] = kmat->size[1];
    emxEnsureCapacity_real_T(b_iidx, i);
    b_iidx_data = b_iidx->data;
    loop_ub = kmat->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_iidx_data[i] = 0.0;
    }
  } else {
    i = b_iidx->size[0];
    b_iidx->size[0] = kmat->size[1];
    emxEnsureCapacity_real_T(b_iidx, i);
    b_iidx_data = b_iidx->data;
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)kmat->size[1],
                (blasint)1, (blasint)kmat->size[0], 1.0, &kmat_data[0], (blasint)
                kmat->size[0], &kweig_data[0], (blasint)kweig->size[0], 0.0,
                &b_iidx_data[0], (blasint)kmat->size[1]);
  }

  emxInit_real_T(&b_a, 2);
  mpower(c_ct, b_a);
  temp_data = b_a->data;
  if ((b_a->size[0] == 0) || (b_a->size[1] == 0) || (b_iidx->size[0] == 0)) {
    i = C->size[0];
    C->size[0] = b_a->size[0];
    emxEnsureCapacity_real_T(C, i);
    C_data = C->data;
    loop_ub = b_a->size[0];
    for (i = 0; i < loop_ub; i++) {
      C_data[i] = 0.0;
    }
  } else {
    i = C->size[0];
    C->size[0] = b_a->size[0];
    emxEnsureCapacity_real_T(C, i);
    C_data = C->data;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)b_a->size[0],
                (blasint)1, (blasint)b_a->size[1], 1.0, &temp_data[0], (blasint)
                b_a->size[0], &b_iidx_data[0], (blasint)b_iidx->size[0], 0.0,
                &C_data[0], (blasint)b_a->size[0]);
  }

  emxInit_real_T(&ord, 1);

  /* 'gkmPWM:423' ord = zeros(m, 1); */
  i = ord->size[0];
  ord->size[0] = PWM->size[0];
  emxEnsureCapacity_real_T(ord, i);
  ord_data = ord->data;
  loop_ub = PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    ord_data[i] = 0.0;
  }

  /* 'gkmPWM:424' tic */
  tic();

  /* 'gkmPWM:425' while i < n */
  emxInit_real_T(&res, 1);
  emxInit_cell_wrap_14(&new_loc);
  emxInit_cell_wrap_2(&new_PWM);
  emxInit_cell_wrap_14(&b_new_loc);
  emxInit_cell_wrap_2(&b_new_PWM);
  emxInit_int32_T(&c_iidx, 2);
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&b_mat, 2);
  emxInit_real_T(&d_PWM, 2);
  emxInit_real_T(&b_loc, 1);
  emxInit_boolean_T(&d_ct, 1);
  do {
    exitg1 = 0;
    if (a < n) {
      /* 'gkmPWM:426' i = i+1; */
      a++;

      /* 'gkmPWM:427' if i~=1 */
      if (a != 1U) {
        /* 'gkmPWM:428' fprintf('Iteration %d: score %f\n', int32(i-1), scorevec(i-1)); */
        printf("Iteration %d: score %f\n", (int)(a - 1U), scorevec_data[(int)a -
               2]);
        fflush(stdout);
      }

      /* 'gkmPWM:430' if mod(i,10) == 0 */
      if (b_mod(a, 10.0) == 0.0) {
        /* 'gkmPWM:431' toc */
        toc();

        /* 'gkmPWM:432' tic */
        tic();

        /* 'gkmPWM:433' fprintf('%d iterations done...\n', int32(i)); */
        printf("%d iterations done...\n", (int)a);
        fflush(stdout);
      }

      /* 'gkmPWM:435' scount = scount + 1; */
      scount++;

      /* 'gkmPWM:436' if  i >= 10 && scount >= 5 && max(-1*diff(scorevec(i-5:i-1))./scorevec(i-4:i-1)) < 0.001 && acount < 5 && i ~= n */
      if ((a >= 10U) && (scount >= 5.0)) {
        b_diff(*(double (*)[5])&scorevec_data[(int)((double)a + -5.0) - 1], b);
        for (i = 0; i < 4; i++) {
          varargin_1[i] = -b[i] / scorevec_data[(int)((double)a + ((double)i +
            -4.0)) - 1];
        }

        if ((e_maximum(varargin_1) < 0.001) && (acount < 5) && (a != n)) {
          /* 'gkmPWM:437' acount = acount + 1; */
          acount++;

          /* 'gkmPWM:438' fprintf('adjusting PWMs after %d iterations (%d)\n', int32(i), int32(acount)); */
          printf("adjusting PWMs after %d iterations (%d)\n", (int)a, acount);
          fflush(stdout);

          /* 'gkmPWM:439' scount = 0; */
          scount = 0.0;

          /* 'gkmPWM:440' for ii = 1:m */
          for (j = 0; j <= varargin_2; j++) {
            /* 'gkmPWM:441' if i/n <= 0.8 && C(ord(ii)) > 0 */
            if ((double)a / n <= 0.8) {
              i = (int)ord_data[j] - 1;
              if (C_data[i] > 0.0) {
                /* 'gkmPWM:442' [PWM{ord(ii)}, lenvec(ord(ii))] = adjust_PWM(PWM{ord(ii)}(l_svm:(length(PWM{ord(ii)})-l_svm+1),:),GC); */
                u0 = b_PWM_data[(int)ord_data[j] - 1].f1->size[0];
                if (u0 < 4) {
                  u0 = 4;
                }

                if (b_PWM_data[(int)ord_data[j] - 1].f1->size[0] == 0) {
                  u0 = 0;
                }

                S = ((double)u0 - l_svm) + 1.0;
                if (l_svm > S) {
                  i1 = -1;
                  i2 = -1;
                } else {
                  i1 = (int)l_svm - 2;
                  i2 = (int)S - 1;
                }

                nx = i2 - i1;
                for (i2 = 0; i2 < 4; i2++) {
                  for (i3 = 0; i3 < nx; i3++) {
                    b_PWM_data[(int)ord_data[j] - 1].f1->data[i3 + nx * i2] =
                      b_PWM_data[i].f1->data[((i1 + i3) + b_PWM_data[i].f1->
                      size[0] * i2) + 1];
                  }
                }

                i1 = b_PWM_data[(int)ord_data[j] - 1].f1->size[0] * b_PWM_data
                  [(int)ord_data[j] - 1].f1->size[1];
                b_PWM_data[(int)ord_data[j] - 1].f1->size[0] = nx;
                b_PWM_data[(int)ord_data[j] - 1].f1->size[1] = 4;
                emxEnsureCapacity_real_T(b_PWM_data[(int)ord_data[j] - 1].f1, i1);
                for (i1 = 0; i1 < 4; i1++) {
                  varargin_1[i1] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] *
                    i1];
                }

                lenvec_data[i] = adjust_PWM(b_PWM_data[(int)ord_data[j] - 1].f1,
                  varargin_1);

                /* 'gkmPWM:443' PWM{ord(ii)} = extendPWM(PWM{ord(ii)}, l_svm-1, GC); */
                /* 'gkmPWM:971' mat = repmat(GCmat, n,1); */
                for (i1 = 0; i1 < 4; i1++) {
                  varargin_1[i1] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] *
                    i1];
                }

                b_repmat(varargin_1, l_svm - 1.0, mat);
                mat_data = mat->data;

                /* 'gkmPWM:972' ext_pwm = [mat;pwm;mat]; */
                i1 = b_mat->size[0] * b_mat->size[1];
                b_mat->size[0] = (mat->size[0] + b_PWM_data[(int)ord_data[j] - 1]
                                  .f1->size[0]) + mat->size[0];
                b_mat->size[1] = 4;
                emxEnsureCapacity_real_T(b_mat, i1);
                diffC_data = b_mat->data;
                loop_ub = mat->size[0];
                b_loop_ub = b_PWM_data[(int)ord_data[j] - 1].f1->size[0];
                for (i1 = 0; i1 < 4; i1++) {
                  for (i2 = 0; i2 < loop_ub; i2++) {
                    diffC_data[i2 + b_mat->size[0] * i1] = mat_data[i2 +
                      mat->size[0] * i1];
                  }

                  for (i2 = 0; i2 < b_loop_ub; i2++) {
                    diffC_data[(i2 + mat->size[0]) + b_mat->size[0] * i1] =
                      b_PWM_data[i].f1->data[i2 + b_PWM_data[i].f1->size[0] * i1];
                  }
                }

                loop_ub = mat->size[0];
                for (i1 = 0; i1 < 4; i1++) {
                  for (i2 = 0; i2 < loop_ub; i2++) {
                    diffC_data[((i2 + mat->size[0]) + b_PWM_data[(int)ord_data[j]
                                - 1].f1->size[0]) + b_mat->size[0] * i1] =
                      mat_data[i2 + mat->size[0] * i1];
                  }
                }

                i1 = b_PWM_data[(int)ord_data[j] - 1].f1->size[0] * b_PWM_data
                  [(int)ord_data[j] - 1].f1->size[1];
                b_PWM_data[(int)ord_data[j] - 1].f1->size[0] = b_mat->size[0];
                b_PWM_data[(int)ord_data[j] - 1].f1->size[1] = 4;
                emxEnsureCapacity_real_T(b_PWM_data[(int)ord_data[j] - 1].f1, i1);
                loop_ub = b_mat->size[0] * 4;
                for (i1 = 0; i1 < loop_ub; i1++) {
                  b_PWM_data[(int)ord_data[j] - 1].f1->data[i1] = diffC_data[i1];
                }

                /* 'gkmPWM:444' loc{ord(ii)} = zeros(lenvec(ord(ii))+2*l_svm-2, 1); */
                S = lenvec_data[i];
                loop_ub = (int)((S + 2.0 * l_svm) - 2.0);
                i = loc_data[(int)ord_data[j] - 1].f1->size[0];
                loc_data[(int)ord_data[j] - 1].f1->size[0] = loop_ub;
                emxEnsureCapacity_real_T(loc_data[(int)ord_data[j] - 1].f1, i);
                for (i = 0; i < loop_ub; i++) {
                  loc_data[(int)ord_data[j] - 1].f1->data[i] = 0.0;
                }

                /* 'gkmPWM:445' loc{ord(ii)}(l_svm:lenvec(ord(ii))+l_svm-1) = 1; */
                S = (S + l_svm) - 1.0;
                if (l_svm > S) {
                  i = -1;
                  i1 = 0;
                } else {
                  i = (int)l_svm - 2;
                  i1 = (int)S;
                }

                loop_ub = (i1 - i) - 1;
                for (i1 = 0; i1 < loop_ub; i1++) {
                  loc_data[(int)ord_data[j] - 1].f1->data[(i + i1) + 1] = 1.0;
                }
              }
            }
          }
        }
      }

      /* 'gkmPWM:449' C = (kmat'*kmat)^(-1)*(kmat'*kweig); */
      b_mtimes(kmat, kmat, c_ct);
      if ((kmat->size[0] == 0) || (kmat->size[1] == 0) || (kweig->size[0] == 0))
      {
        i = b_iidx->size[0];
        b_iidx->size[0] = kmat->size[1];
        emxEnsureCapacity_real_T(b_iidx, i);
        b_iidx_data = b_iidx->data;
        loop_ub = kmat->size[1];
        for (i = 0; i < loop_ub; i++) {
          b_iidx_data[i] = 0.0;
        }
      } else {
        i = b_iidx->size[0];
        b_iidx->size[0] = kmat->size[1];
        emxEnsureCapacity_real_T(b_iidx, i);
        b_iidx_data = b_iidx->data;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)kmat->
                    size[1], (blasint)1, (blasint)kmat->size[0], 1.0,
                    &kmat_data[0], (blasint)kmat->size[0], &kweig_data[0],
                    (blasint)kweig->size[0], 0.0, &b_iidx_data[0], (blasint)
                    kmat->size[1]);
      }

      mpower(c_ct, b_a);
      temp_data = b_a->data;
      if ((b_a->size[0] == 0) || (b_a->size[1] == 0) || (b_iidx->size[0] == 0))
      {
        i = C->size[0];
        C->size[0] = b_a->size[0];
        emxEnsureCapacity_real_T(C, i);
        C_data = C->data;
        loop_ub = b_a->size[0];
        for (i = 0; i < loop_ub; i++) {
          C_data[i] = 0.0;
        }
      } else {
        i = C->size[0];
        C->size[0] = b_a->size[0];
        emxEnsureCapacity_real_T(C, i);
        C_data = C->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                    b_a->size[0], (blasint)1, (blasint)b_a->size[1], 1.0,
                    &temp_data[0], (blasint)b_a->size[0], &b_iidx_data[0],
                    (blasint)b_iidx->size[0], 0.0, &C_data[0], (blasint)
                    b_a->size[0]);
      }

      /* 'gkmPWM:450' res = kweig-kmat*C; */
      if ((kmat->size[0] == 0) || (kmat->size[1] == 0) || (C->size[0] == 0)) {
        i = res->size[0];
        res->size[0] = kmat->size[0];
        emxEnsureCapacity_real_T(res, i);
        res_data = res->data;
        loop_ub = kmat->size[0];
        for (i = 0; i < loop_ub; i++) {
          res_data[i] = 0.0;
        }
      } else {
        i = res->size[0];
        res->size[0] = kmat->size[0];
        emxEnsureCapacity_real_T(res, i);
        res_data = res->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                    kmat->size[0], (blasint)1, (blasint)kmat->size[1], 1.0,
                    &kmat_data[0], (blasint)kmat->size[0], &C_data[0], (blasint)
                    C->size[0], 0.0, &res_data[0], (blasint)kmat->size[0]);
      }

      if (kweig->size[0] == res->size[0]) {
        i = res->size[0];
        res->size[0] = kweig->size[0];
        emxEnsureCapacity_real_T(res, i);
        res_data = res->data;
        loop_ub = kweig->size[0];
        for (i = 0; i < loop_ub; i++) {
          res_data[i] = kweig_data[i] - res_data[i];
        }
      } else {
        minus(res, kweig);
        res_data = res->data;
      }

      /* 'gkmPWM:451' corrvec = zeros(m,1); */
      /* 'gkmPWM:452' for ii = 1:m */
      i = ord->size[0];
      ord->size[0] = PWM->size[0];
      emxEnsureCapacity_real_T(ord, i);
      ord_data = ord->data;
      for (j = 0; j <= varargin_2; j++) {
        /* 'gkmPWM:453' [~,ind] = sort(kmat(:,ii), 'descend'); */
        loop_ub = kmat->size[0];
        i = f->size[0];
        f->size[0] = kmat->size[0];
        emxEnsureCapacity_real_T(f, i);
        f_data = f->data;
        for (i = 0; i < loop_ub; i++) {
          f_data[i] = kmat_data[i + kmat->size[0] * j];
        }

        sort(f, iidx);
        iidx_data = iidx->data;
        i = f->size[0];
        f->size[0] = iidx->size[0];
        emxEnsureCapacity_real_T(f, i);
        f_data = f->data;
        loop_ub = iidx->size[0];
        for (i = 0; i < loop_ub; i++) {
          f_data[i] = iidx_data[i];
        }

        /* 'gkmPWM:454' corrvec(ii) = sum(res(ind(1:lcomb)).^2); */
        if (1 > lcomb) {
          loop_ub = 0;
        } else {
          loop_ub = lcomb;
        }

        i = b_iidx->size[0];
        b_iidx->size[0] = loop_ub;
        emxEnsureCapacity_real_T(b_iidx, i);
        b_iidx_data = b_iidx->data;
        for (i = 0; i < loop_ub; i++) {
          M = res_data[(int)f_data[i] - 1];
          b_iidx_data[i] = pow(M, 2.0);
        }

        ord_data[j] = blockedSummation(b_iidx, b_iidx->size[0]);

        /* 'gkmPWM:455' if mod(i,20) == 0 */
        if (b_mod(a, 20.0) == 0.0) {
          /* 'gkmPWM:456' if RC */
          if (RC != 0.0) {
            /* 'gkmPWM:457' kmat(:,ii) = PWM2kmers(PWM{ii},negmat,rc,diffc,indc,loc{ii},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(ii)); */
            PWM2kmers(b_PWM_data[j].f1, negmat, rc, diffc, indc, loc_data[j].f1,
                      xc, l_svm, k_svm, rcnum, b_iidx);
            b_iidx_data = b_iidx->data;
            S = (l_svm - 1.0) + lenvec_data[j];
            if (b_iidx->size[0] == negvec->size[0]) {
              loop_ub = b_iidx->size[0];
              for (i = 0; i < loop_ub; i++) {
                kmat_data[i + kmat->size[0] * j] = b_iidx_data[i] -
                  negvec_data[i] * S;
              }
            } else {
              m_binary_expand_op(kmat, j, b_iidx, negvec, S);
              kmat_data = kmat->data;
            }
          } else {
            /* 'gkmPWM:458' else */
            /* 'gkmPWM:459' kmat(:,ii) = PWM2kmers_norc(PWM{ii},negmat,rc,diffc,indc,loc{ii},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(ii)); */
            PWM2kmers_norc(b_PWM_data[j].f1, negmat, rc, diffc, indc, loc_data[j]
                           .f1, xc, l_svm, k_svm, b_iidx);
            b_iidx_data = b_iidx->data;
            S = (l_svm - 1.0) + lenvec_data[j];
            if (b_iidx->size[0] == negvec->size[0]) {
              loop_ub = b_iidx->size[0];
              for (i = 0; i < loop_ub; i++) {
                kmat_data[i + kmat->size[0] * j] = b_iidx_data[i] -
                  negvec_data[i] * S;
              }
            } else {
              m_binary_expand_op(kmat, j, b_iidx, negvec, S);
              kmat_data = kmat->data;
            }
          }
        }
      }

      /* The order of PWM optimization is determined by the correlation of its top 110 kmers with the gapped kmer weight vector */
      /* 'gkmPWM:464' [~,ord] = sort(corrvec, 'descend'); */
      sort(ord, iidx);
      iidx_data = iidx->data;
      i = ord->size[0];
      ord->size[0] = iidx->size[0];
      emxEnsureCapacity_real_T(ord, i);
      ord_data = ord->data;
      loop_ub = iidx->size[0];
      for (i = 0; i < loop_ub; i++) {
        ord_data[i] = iidx_data[i];
      }

      /* 'gkmPWM:465' for ii = 1:m */
      for (j = 0; j <= varargin_2; j++) {
        /* The order of the column optimization is determined by the max probability in each column */
        /* 'gkmPWM:467' v = max(PWM{ord(ii)}(l_svm:lenvec(ord(ii))+l_svm-1,:)'); */
        i = (int)ord_data[j] - 1;
        S = (lenvec_data[i] + l_svm) - 1.0;
        if (l_svm > S) {
          i1 = 0;
          i2 = 0;
        } else {
          i1 = (int)l_svm - 1;
          i2 = (int)S;
        }

        /* 'gkmPWM:468' [~,c] = sort(v, 'ascend'); */
        i3 = d_PWM->size[0] * d_PWM->size[1];
        d_PWM->size[0] = 4;
        loop_ub = i2 - i1;
        d_PWM->size[1] = loop_ub;
        emxEnsureCapacity_real_T(d_PWM, i3);
        diffC_data = d_PWM->data;
        for (i2 = 0; i2 < loop_ub; i2++) {
          for (i3 = 0; i3 < 4; i3++) {
            diffC_data[i3 + 4 * i2] = b_PWM_data[i].f1->data[(i1 + i2) +
              b_PWM_data[i].f1->size[0] * i3];
          }
        }

        f_maximum(d_PWM, temp);
        e_sort(temp, c_iidx);
        iidx_data = c_iidx->data;
        i1 = temp->size[0] * temp->size[1];
        temp->size[0] = 1;
        temp->size[1] = c_iidx->size[1];
        emxEnsureCapacity_real_T(temp, i1);
        temp_data = temp->data;
        loop_ub = c_iidx->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          temp_data[i1] = iidx_data[i1];
        }

        /* C = (kmat'*kmat)^(-1)*(kmat'*kweig); */
        /* res = kweig-kmat*C; */
        /* 'gkmPWM:471' for iii = 1:length(c) */
        i1 = temp->size[1];
        for (jj = 0; jj < i1; jj++) {
          /* 'gkmPWM:472' PWMtemp = PWM{ord(ii)}(c(iii):c(iii)+l_svm*2-2,:); */
          S = temp_data[jj];
          M = ((double)(int)S + l_svm * 2.0) - 2.0;
          if ((int)S > M) {
            i2 = 0;
            i3 = 0;
            nx = 0;
            u0 = 0;
          } else {
            i2 = (int)S - 1;
            i3 = (int)M;
            nx = (int)S - 1;
            u0 = (int)M;
          }

          /* 'gkmPWM:473' [kweigdiff,PWM{ord(ii)}(c(iii)+l_svm-1,:)] = getEMprob_v3(PWMtemp,res/C(ord(ii)),negmat,poscell,rc,diffC,indc,loc{ord(ii)}(c(iii):c(iii)+2*l_svm-2),xc,reg,l_svm,k_svm,rcnum,RC); */
          loop_ub = i3 - i2;
          i3 = mat->size[0] * mat->size[1];
          mat->size[0] = loop_ub;
          mat->size[1] = 4;
          emxEnsureCapacity_real_T(mat, i3);
          mat_data = mat->data;
          for (i3 = 0; i3 < 4; i3++) {
            for (b_loop_ub = 0; b_loop_ub < loop_ub; b_loop_ub++) {
              mat_data[b_loop_ub + mat->size[0] * i3] = b_PWM_data[i].f1->data
                [(i2 + b_loop_ub) + b_PWM_data[i].f1->size[0] * i3];
            }
          }

          loop_ub = res->size[0];
          i2 = b_iidx->size[0];
          b_iidx->size[0] = res->size[0];
          emxEnsureCapacity_real_T(b_iidx, i2);
          b_iidx_data = b_iidx->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            b_iidx_data[i2] = res_data[i2] / C_data[i];
          }

          loop_ub = u0 - nx;
          i2 = b_loc->size[0];
          b_loc->size[0] = loop_ub;
          emxEnsureCapacity_real_T(b_loc, i2);
          diffC_data = b_loc->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            diffC_data[i2] = loc_data[i].f1->data[nx + i2];
          }

          getEMprob_v3(mat, b_iidx, negmat, poscell, rc, diffC, indc, b_loc, xc,
                       reg, l_svm, k_svm, rcnum, RC, f, tmp_data, &nx);
          f_data = f->data;
          i2 = (int)(((double)(int)temp_data[jj] + l_svm) - 1.0) - 1;
          for (i3 = 0; i3 < 4; i3++) {
            b_PWM_data[(int)ord_data[j] - 1].f1->data[i2 + b_PWM_data[(int)
              ord_data[j] - 1].f1->size[0] * i3] = tmp_data[i3];
          }

          /* 'gkmPWM:474' kmat(:,ord(ii)) = kmat(:,ord(ii)) + kweigdiff; */
          if (kmat->size[0] == f->size[0]) {
            nx = kmat->size[0] - 1;
            S = ord_data[j];
            i2 = b_iidx->size[0];
            b_iidx->size[0] = kmat->size[0];
            emxEnsureCapacity_real_T(b_iidx, i2);
            b_iidx_data = b_iidx->data;
            for (i2 = 0; i2 <= nx; i2++) {
              b_iidx_data[i2] = kmat_data[i2 + kmat->size[0] * ((int)S - 1)] +
                f_data[i2];
            }

            loop_ub = b_iidx->size[0];
            for (i2 = 0; i2 < loop_ub; i2++) {
              kmat_data[i2 + kmat->size[0] * ((int)S - 1)] = b_iidx_data[i2];
            }
          } else {
            o_binary_expand_op(kmat, ord, j, f);
            kmat_data = kmat->data;
          }

          /* 'gkmPWM:475' res = res-kweigdiff*C(ord(ii)); */
          loop_ub = res->size[0];
          if (res->size[0] == f->size[0]) {
            for (i2 = 0; i2 < loop_ub; i2++) {
              res_data[i2] -= f_data[i2] * C_data[i];
            }
          } else {
            n_binary_expand_op(res, f, C, ord, j);
            res_data = res->data;
          }
        }
      }

      /* Reseed PWMs if two or more of them are too highly correlated */
      /* 'gkmPWM:479' if i/n <= 0.80 */
      if ((double)a / n <= 0.8) {
        /* 'gkmPWM:480' info = avg_info(PWM,l_svm); */
        avg_info(c_PWM, l_svm, f);

        /* 'gkmPWM:481' [~,ord] = sort(info,'descend'); */
        sort(f, iidx);
        iidx_data = iidx->data;
        i = ord->size[0];
        ord->size[0] = iidx->size[0];
        emxEnsureCapacity_real_T(ord, i);
        ord_data = ord->data;
        loop_ub = iidx->size[0];
        for (i = 0; i < loop_ub; i++) {
          ord_data[i] = iidx_data[i];
        }

        /* 'gkmPWM:482' kmat = kmat(:,ord); */
        nx = kmat->size[0] - 1;
        i = c_ct->size[0] * c_ct->size[1];
        c_ct->size[0] = kmat->size[0];
        c_ct->size[1] = ord->size[0];
        emxEnsureCapacity_real_T(c_ct, i);
        mat_data = c_ct->data;
        loop_ub = ord->size[0];
        for (i = 0; i < loop_ub; i++) {
          for (i1 = 0; i1 <= nx; i1++) {
            mat_data[i1 + c_ct->size[0] * i] = kmat_data[i1 + kmat->size[0] *
              ((int)ord_data[i] - 1)];
          }
        }

        i = kmat->size[0] * kmat->size[1];
        kmat->size[0] = c_ct->size[0];
        kmat->size[1] = c_ct->size[1];
        emxEnsureCapacity_real_T(kmat, i);
        kmat_data = kmat->data;
        loop_ub = c_ct->size[0] * c_ct->size[1];
        for (i = 0; i < loop_ub; i++) {
          kmat_data[i] = mat_data[i];
        }

        /* 'gkmPWM:483' lenvec = lenvec(ord); */
        i = b_iidx->size[0];
        b_iidx->size[0] = ord->size[0];
        emxEnsureCapacity_real_T(b_iidx, i);
        b_iidx_data = b_iidx->data;
        loop_ub = ord->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_iidx_data[i] = lenvec_data[(int)ord_data[i] - 1];
        }

        i = lenvec->size[0];
        lenvec->size[0] = b_iidx->size[0];
        emxEnsureCapacity_real_T(lenvec, i);
        lenvec_data = lenvec->data;
        loop_ub = b_iidx->size[0];
        for (i = 0; i < loop_ub; i++) {
          lenvec_data[i] = b_iidx_data[i];
        }

        /*  loc = loc(ord); */
        /* 'gkmPWM:486' ord_len = length(ord); */
        /* 'gkmPWM:487' new_loc = cell(ord_len,1); */
        nx = ord->size[0];
        i = new_loc->size[0];
        new_loc->size[0] = ord->size[0];
        emxEnsureCapacity_cell_wrap_14(new_loc, i);
        poscell_data = new_loc->data;
        for (i = 0; i < nx; i++) {
          poscell_data[i].f1->size[0] = 0;
        }

        /* 'gkmPWM:488' new_loc = coder.nullcopy(new_loc); */
        i = b_new_loc->size[0];
        b_new_loc->size[0] = new_loc->size[0];
        emxEnsureCapacity_cell_wrap_14(b_new_loc, i);
        poscell_data = b_new_loc->data;

        /* 'gkmPWM:489' for cur_idx=1:ord_len */
        i = ord->size[0];
        for (u0 = 0; u0 < i; u0++) {
          /* 'gkmPWM:490' new_loc{cur_idx} = loc{ord(cur_idx)}; */
          i1 = poscell_data[u0].f1->size[0];
          poscell_data[u0].f1->size[0] = loc_data[(int)ord_data[u0] - 1]
            .f1->size[0];
          emxEnsureCapacity_real_T(poscell_data[u0].f1, i1);
          loop_ub = loc_data[(int)ord_data[u0] - 1].f1->size[0];
          for (i1 = 0; i1 < loop_ub; i1++) {
            poscell_data[u0].f1->data[i1] = loc_data[(int)ord_data[u0] - 1]
              .f1->data[i1];
          }
        }

        /* 'gkmPWM:492' loc = coder.nullcopy(loc); */
        /* 'gkmPWM:493' for cur_idx=1:ord_len */
        i = ord->size[0];
        for (u0 = 0; u0 < i; u0++) {
          /* 'gkmPWM:494' loc{cur_idx} = new_loc{cur_idx}; */
          loop_ub = poscell_data[u0].f1->size[0];
          i1 = loc_data[u0].f1->size[0];
          loc_data[u0].f1->size[0] = poscell_data[u0].f1->size[0];
          emxEnsureCapacity_real_T(loc_data[u0].f1, i1);
          for (i1 = 0; i1 < loop_ub; i1++) {
            loc_data[u0].f1->data[i1] = poscell_data[u0].f1->data[i1];
          }
        }

        /*  PWM = PWM(ord); */
        /* 'gkmPWM:499' new_PWM = cell(ord_len,1); */
        i = new_PWM->size[0];
        new_PWM->size[0] = ord->size[0];
        emxEnsureCapacity_cell_wrap_2(new_PWM, i);
        new_PWM_data = new_PWM->data;
        for (i = 0; i < nx; i++) {
          new_PWM_data[i].f1->size[0] = 0;
          new_PWM_data[i].f1->size[1] = 4;
        }

        /* 'gkmPWM:500' new_PWM = coder.nullcopy(new_PWM); */
        i = b_new_PWM->size[0];
        b_new_PWM->size[0] = new_PWM->size[0];
        emxEnsureCapacity_cell_wrap_2(b_new_PWM, i);
        new_PWM_data = b_new_PWM->data;

        /* 'gkmPWM:501' for cur_idx=1:ord_len */
        i = ord->size[0];
        for (u0 = 0; u0 < i; u0++) {
          /* 'gkmPWM:502' new_PWM{cur_idx} = PWM{ord(cur_idx)}; */
          i1 = new_PWM_data[u0].f1->size[0] * new_PWM_data[u0].f1->size[1];
          new_PWM_data[u0].f1->size[0] = b_PWM_data[(int)ord_data[u0] - 1]
            .f1->size[0];
          new_PWM_data[u0].f1->size[1] = 4;
          emxEnsureCapacity_real_T(new_PWM_data[u0].f1, i1);
          loop_ub = b_PWM_data[(int)ord_data[u0] - 1].f1->size[0] * 4;
          for (i1 = 0; i1 < loop_ub; i1++) {
            new_PWM_data[u0].f1->data[i1] = b_PWM_data[(int)ord_data[u0] - 1].
              f1->data[i1];
          }
        }

        /* 'gkmPWM:504' PWM = coder.nullcopy(PWM); */
        /* 'gkmPWM:505' for cur_idx=1:ord_len */
        i = ord->size[0];
        for (u0 = 0; u0 < i; u0++) {
          /* 'gkmPWM:506' PWM{cur_idx} = new_PWM{cur_idx}; */
          i1 = b_PWM_data[u0].f1->size[0] * b_PWM_data[u0].f1->size[1];
          b_PWM_data[u0].f1->size[0] = new_PWM_data[u0].f1->size[0];
          b_PWM_data[u0].f1->size[1] = 4;
          emxEnsureCapacity_real_T(b_PWM_data[u0].f1, i1);
          loop_ub = new_PWM_data[u0].f1->size[0] * 4;
          for (i1 = 0; i1 < loop_ub; i1++) {
            b_PWM_data[u0].f1->data[i1] = new_PWM_data[u0].f1->data[i1];
          }
        }

        /* 'gkmPWM:509' C = C(ord); */
        i = b_iidx->size[0];
        b_iidx->size[0] = ord->size[0];
        emxEnsureCapacity_real_T(b_iidx, i);
        b_iidx_data = b_iidx->data;
        loop_ub = ord->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_iidx_data[i] = C_data[(int)ord_data[i] - 1];
        }

        i = C->size[0];
        C->size[0] = b_iidx->size[0];
        emxEnsureCapacity_real_T(C, i);
        C_data = C->data;
        loop_ub = b_iidx->size[0];
        for (i = 0; i < loop_ub; i++) {
          C_data[i] = b_iidx_data[i];
        }

        /* 'gkmPWM:510' for j = 1:m */
        for (j = 0; j <= varargin_2; j++) {
          /* 'gkmPWM:511' KMAT(:,j) = kmat(:,j)/sqrt(kmat(:,j)'*kmat(:,j)); */
          loop_ub = kmat->size[0];
          i = f->size[0];
          f->size[0] = kmat->size[0];
          emxEnsureCapacity_real_T(f, i);
          f_data = f->data;
          for (i = 0; i < loop_ub; i++) {
            f_data[i] = kmat_data[i + kmat->size[0] * j];
          }

          loop_ub = kmat->size[0];
          i = b_iidx->size[0];
          b_iidx->size[0] = kmat->size[0];
          emxEnsureCapacity_real_T(b_iidx, i);
          b_iidx_data = b_iidx->data;
          for (i = 0; i < loop_ub; i++) {
            b_iidx_data[i] = kmat_data[i + kmat->size[0] * j];
          }

          loop_ub = kmat->size[0];
          if (kmat->size[0] < 1) {
            c = 0.0;
          } else {
            c = cblas_ddot((blasint)kmat->size[0], &f_data[0], (blasint)1,
                           &b_iidx_data[0], (blasint)1);
          }

          M = sqrt(c);
          for (i = 0; i < loop_ub; i++) {
            CT_data[i + CT->size[0] * j] = kmat_data[i + kmat->size[0] * j] / M;
          }
        }

        /* 'gkmPWM:513' MAT = KMAT'*KMAT; */
        b_mtimes(CT, CT, ct);
        ct_data = ct->data;

        /* 'gkmPWM:514' for j = 1:m-1 */
        i = PWM->size[0];
        for (j = 0; j <= i - 2; j++) {
          /* 'gkmPWM:515' vec = MAT(j+1:end,j); */
          if (j + 2 > ct->size[0]) {
            i1 = 0;
            i2 = 0;
          } else {
            i1 = j + 1;
            i2 = ct->size[0];
          }

          /* 'gkmPWM:516' [a b] = max(vec); */
          loop_ub = i2 - i1;
          i2 = b_iidx->size[0];
          b_iidx->size[0] = loop_ub;
          emxEnsureCapacity_real_T(b_iidx, i2);
          b_iidx_data = b_iidx->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            b_iidx_data[i2] = ct_data[(i1 + i2) + ct->size[0] * j];
          }

          d_maximum(b_iidx, &M, &nx);

          /* 'gkmPWM:517' if a > rcorr && C(j) > 0 */
          if ((M > rcorr) && (C_data[j] > 0.0)) {
            /* 'gkmPWM:518' scount = 0; */
            scount = 0.0;

            /* 'gkmPWM:519' fprintf('reseeding\n'); */
            printf("reseeding\n");
            fflush(stdout);

            /* 'gkmPWM:520' f = j+find(vec > rcorr); */
            i2 = d_ct->size[0];
            d_ct->size[0] = loop_ub;
            emxEnsureCapacity_boolean_T(d_ct, i2);
            b_ct_data = d_ct->data;
            for (i2 = 0; i2 < loop_ub; i2++) {
              b_ct_data[i2] = (ct_data[(i1 + i2) + ct->size[0] * j] > rcorr);
            }

            b_eml_find(d_ct, iidx);
            iidx_data = iidx->data;
            i1 = f->size[0];
            f->size[0] = iidx->size[0];
            emxEnsureCapacity_real_T(f, i1);
            f_data = f->data;
            loop_ub = iidx->size[0];
            for (i1 = 0; i1 < loop_ub; i1++) {
              f_data[i1] = ((double)j + 1.0) + (double)iidx_data[i1];
            }

            /* 'gkmPWM:521' for jj = 1:length(f) */
            i1 = f->size[0];
            if (0 <= f->size[0] - 1) {
              b_lenvec[1] = 12.0;
            }

            for (jj = 0; jj < i1; jj++) {
              /* 'gkmPWM:522' for jjj = 1:min([lenvec(f(jj)) 12]) */
              u0 = (int)f_data[jj] - 1;
              M = lenvec_data[u0];
              b_lenvec[0] = M;
              i2 = (int)b_minimum(b_lenvec);
              for (b_loop_ub = 0; b_loop_ub < i2; b_loop_ub++) {
                /* 'gkmPWM:523' PWM{f(jj)}(jj+9,:) =  PWM{f(jj)}(jjj+l_svm-1,randperm(4)); */
                randperm(b);
                nx = (int)((((double)b_loop_ub + 1.0) + l_svm) - 1.0);
                for (i3 = 0; i3 < 4; i3++) {
                  varargin_1[i3] = b_PWM_data[u0].f1->data[(nx + b_PWM_data[u0].
                    f1->size[0] * ((int)b[i3] - 1)) - 1];
                }

                for (i3 = 0; i3 < 4; i3++) {
                  b_PWM_data[(int)f_data[jj] - 1].f1->data[(jj + b_PWM_data[(int)
                    f_data[jj] - 1].f1->size[0] * i3) + 9] = varargin_1[i3];
                }
              }

              /* 'gkmPWM:525' if lenvec(f(jj)) >= 12 */
              if (M >= 12.0) {
                /* 'gkmPWM:526' PWM{f(jj)} = PWM{f(jj)}(l_svm:l_svm+11,:); */
                i2 = mat->size[0] * mat->size[1];
                mat->size[0] = 12;
                mat->size[1] = 4;
                emxEnsureCapacity_real_T(mat, i2);
                mat_data = mat->data;
                for (i2 = 0; i2 < 4; i2++) {
                  for (i3 = 0; i3 < 12; i3++) {
                    mat_data[i3 + mat->size[0] * i2] = b_PWM_data[u0].f1->data
                      [((int)(l_svm + (double)i3) + b_PWM_data[u0].f1->size[0] *
                        i2) - 1];
                  }
                }

                i2 = b_PWM_data[(int)f_data[jj] - 1].f1->size[0] * b_PWM_data
                  [(int)f_data[jj] - 1].f1->size[1];
                b_PWM_data[(int)f_data[jj] - 1].f1->size[0] = 12;
                b_PWM_data[(int)f_data[jj] - 1].f1->size[1] = 4;
                emxEnsureCapacity_real_T(b_PWM_data[(int)f_data[jj] - 1].f1, i2);
                for (i2 = 0; i2 < 48; i2++) {
                  b_PWM_data[(int)f_data[jj] - 1].f1->data[i2] = mat_data[i2];
                }

                /* 'gkmPWM:527' PWM{f(jj)} = extendPWM(PWM{f(jj)}, l_svm-1, GC); */
                /* 'gkmPWM:971' mat = repmat(GCmat, n,1); */
                for (i2 = 0; i2 < 4; i2++) {
                  varargin_1[i2] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] *
                    i2];
                }

                b_repmat(varargin_1, l_svm - 1.0, mat);
                mat_data = mat->data;

                /* 'gkmPWM:972' ext_pwm = [mat;pwm;mat]; */
                i2 = b_mat->size[0] * b_mat->size[1];
                b_mat->size[0] = (mat->size[0] + b_PWM_data[(int)f_data[jj] - 1]
                                  .f1->size[0]) + mat->size[0];
                b_mat->size[1] = 4;
                emxEnsureCapacity_real_T(b_mat, i2);
                diffC_data = b_mat->data;
                loop_ub = mat->size[0];
                b_loop_ub = b_PWM_data[(int)f_data[jj] - 1].f1->size[0];
                for (i2 = 0; i2 < 4; i2++) {
                  for (i3 = 0; i3 < loop_ub; i3++) {
                    diffC_data[i3 + b_mat->size[0] * i2] = mat_data[i3 +
                      mat->size[0] * i2];
                  }

                  for (i3 = 0; i3 < b_loop_ub; i3++) {
                    diffC_data[(i3 + mat->size[0]) + b_mat->size[0] * i2] =
                      b_PWM_data[u0].f1->data[i3 + b_PWM_data[u0].f1->size[0] *
                      i2];
                  }
                }

                loop_ub = mat->size[0];
                for (i2 = 0; i2 < 4; i2++) {
                  for (i3 = 0; i3 < loop_ub; i3++) {
                    diffC_data[((i3 + mat->size[0]) + b_PWM_data[(int)f_data[jj]
                                - 1].f1->size[0]) + b_mat->size[0] * i2] =
                      mat_data[i3 + mat->size[0] * i2];
                  }
                }

                i2 = b_PWM_data[(int)f_data[jj] - 1].f1->size[0] * b_PWM_data
                  [(int)f_data[jj] - 1].f1->size[1];
                b_PWM_data[(int)f_data[jj] - 1].f1->size[0] = b_mat->size[0];
                b_PWM_data[(int)f_data[jj] - 1].f1->size[1] = 4;
                emxEnsureCapacity_real_T(b_PWM_data[(int)f_data[jj] - 1].f1, i2);
                loop_ub = b_mat->size[0] * 4;
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_PWM_data[(int)f_data[jj] - 1].f1->data[i2] = diffC_data[i2];
                }

                /* 'gkmPWM:528' lenvec(f(jj)) = 12; */
                lenvec_data[u0] = 12.0;

                /* 'gkmPWM:529' loc{f(jj)} = zeros(lenvec(f(jj))+l_svm*2-2, 1); */
                loop_ub = (int)((lenvec_data[u0] + l_svm * 2.0) - 2.0);
                i2 = loc_data[(int)f_data[jj] - 1].f1->size[0];
                loc_data[(int)f_data[jj] - 1].f1->size[0] = loop_ub;
                emxEnsureCapacity_real_T(loc_data[(int)f_data[jj] - 1].f1, i2);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  loc_data[(int)f_data[jj] - 1].f1->data[i2] = 0.0;
                }

                /* 'gkmPWM:530' loc{f(jj)}(l_svm:lenvec(f(jj))+l_svm-1) = 1; */
                S = (lenvec_data[u0] + l_svm) - 1.0;
                if (l_svm > S) {
                  i2 = -1;
                  i3 = 0;
                } else {
                  i2 = (int)l_svm - 2;
                  i3 = (int)S;
                }

                loop_ub = (i3 - i2) - 1;
                for (i3 = 0; i3 < loop_ub; i3++) {
                  loc_data[(int)f_data[jj] - 1].f1->data[(i2 + i3) + 1] = 1.0;
                }
              }

              /* 'gkmPWM:532' if RC */
              if (RC != 0.0) {
                /* 'gkmPWM:533' kmat(:,f(jj)) = PWM2kmers(PWM{f(jj)},negmat,rc,diffc,indc,loc{f(jj)},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(f(jj))); */
                PWM2kmers(b_PWM_data[u0].f1, negmat, rc, diffc, indc,
                          loc_data[u0].f1, xc, l_svm, k_svm, rcnum, b_iidx);
                b_iidx_data = b_iidx->data;
                S = (l_svm - 1.0) + lenvec_data[u0];
                if (b_iidx->size[0] == negvec->size[0]) {
                  loop_ub = b_iidx->size[0];
                  for (i2 = 0; i2 < loop_ub; i2++) {
                    kmat_data[i2 + kmat->size[0] * u0] = b_iidx_data[i2] -
                      negvec_data[i2] * S;
                  }
                } else {
                  p_binary_expand_op(kmat, f, jj, b_iidx, negvec, S);
                  kmat_data = kmat->data;
                }
              } else {
                /* 'gkmPWM:534' else */
                /* 'gkmPWM:535' kmat(:,f(jj)) = PWM2kmers_norc(PWM{f(jj)},negmat,rc,diffc,indc,loc{f(jj)},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(f(jj))); */
                PWM2kmers_norc(b_PWM_data[u0].f1, negmat, rc, diffc, indc,
                               loc_data[u0].f1, xc, l_svm, k_svm, b_iidx);
                b_iidx_data = b_iidx->data;
                S = (l_svm - 1.0) + lenvec_data[u0];
                if (b_iidx->size[0] == negvec->size[0]) {
                  loop_ub = b_iidx->size[0];
                  for (i2 = 0; i2 < loop_ub; i2++) {
                    kmat_data[i2 + kmat->size[0] * u0] = b_iidx_data[i2] -
                      negvec_data[i2] * S;
                  }
                } else {
                  p_binary_expand_op(kmat, f, jj, b_iidx, negvec, S);
                  kmat_data = kmat->data;
                }
              }
            }
          }
        }
      }

      /* Breaks the loop if it looks like it converged */
      /* 'gkmPWM:542' C = (kmat'*kmat)^(-1)*(kmat'*kweig); */
      b_mtimes(kmat, kmat, c_ct);
      if ((kmat->size[0] == 0) || (kmat->size[1] == 0) || (kweig->size[0] == 0))
      {
        i = b_iidx->size[0];
        b_iidx->size[0] = kmat->size[1];
        emxEnsureCapacity_real_T(b_iidx, i);
        b_iidx_data = b_iidx->data;
        loop_ub = kmat->size[1];
        for (i = 0; i < loop_ub; i++) {
          b_iidx_data[i] = 0.0;
        }
      } else {
        i = b_iidx->size[0];
        b_iidx->size[0] = kmat->size[1];
        emxEnsureCapacity_real_T(b_iidx, i);
        b_iidx_data = b_iidx->data;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)kmat->
                    size[1], (blasint)1, (blasint)kmat->size[0], 1.0,
                    &kmat_data[0], (blasint)kmat->size[0], &kweig_data[0],
                    (blasint)kweig->size[0], 0.0, &b_iidx_data[0], (blasint)
                    kmat->size[1]);
      }

      mpower(c_ct, b_a);
      temp_data = b_a->data;
      if ((b_a->size[0] == 0) || (b_a->size[1] == 0) || (b_iidx->size[0] == 0))
      {
        i = C->size[0];
        C->size[0] = b_a->size[0];
        emxEnsureCapacity_real_T(C, i);
        C_data = C->data;
        loop_ub = b_a->size[0];
        for (i = 0; i < loop_ub; i++) {
          C_data[i] = 0.0;
        }
      } else {
        i = C->size[0];
        C->size[0] = b_a->size[0];
        emxEnsureCapacity_real_T(C, i);
        C_data = C->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)
                    b_a->size[0], (blasint)1, (blasint)b_a->size[1], 1.0,
                    &temp_data[0], (blasint)b_a->size[0], &b_iidx_data[0],
                    (blasint)b_iidx->size[0], 0.0, &C_data[0], (blasint)
                    b_a->size[0]);
      }

      /* 'gkmPWM:543' scorevec(i) = sqrt(res'*res); */
      if (res->size[0] < 1) {
        c = 0.0;
      } else {
        c = cblas_ddot((blasint)res->size[0], &res_data[0], (blasint)1,
                       &res_data[0], (blasint)1);
      }

      scorevec_data[(int)a - 1] = sqrt(c);

      /* 'gkmPWM:544' if i >= 10 && acount == 5 && scount >= 10 && max(abs(diff(scorevec(i-9:i))./scorevec(i-8:i))) < 0.001 */
      guard1 = false;
      if (((int)a >= 10) && (acount == 5) && (scount >= 10.0)) {
        c_diff(*(double (*)[10])&scorevec_data[(int)a - 10], b_varargin_1);
        for (u0 = 0; u0 < 9; u0++) {
          b_varargin_1[u0] = fabs(b_varargin_1[u0] / scorevec_data[(u0 + (int)a)
            - 9]);
        }

        if (i_maximum(b_varargin_1) < 0.001) {
          /* 'gkmPWM:545' scorevec = scorevec(1:i); */
          i = scorevec->size[0] * scorevec->size[1];
          scorevec->size[1] = (int)a;
          emxEnsureCapacity_real_T(scorevec, i);
          scorevec_data = scorevec->data;
          exitg1 = 1;
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1 && (((int)a > 10) && (acount == 5) && (scount >= 10.0))) {
        /* 'gkmPWM:548' if i > 10 && acount == 5 && scount >= 10 && sum(diff(scorevec(i-9:i))>0) > 7 */
        c_diff(*(double (*)[10])&scorevec_data[(int)a - 10], b_varargin_1);
        for (i = 0; i < 9; i++) {
          x[i] = (b_varargin_1[i] > 0.0);
        }

        nx = x[0];
        for (u0 = 0; u0 < 8; u0++) {
          nx += x[u0 + 1];
        }

        if (nx > 7) {
          /* 'gkmPWM:549' scorevec = scorevec(1:i); */
          i = scorevec->size[0] * scorevec->size[1];
          scorevec->size[1] = (int)a;
          emxEnsureCapacity_real_T(scorevec, i);
          scorevec_data = scorevec->data;
          exitg1 = 1;
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_boolean_T(&d_ct);
  emxFree_real_T(&b_loc);
  emxFree_real_T(&d_PWM);
  emxFree_real_T(&b_mat);
  emxFree_real_T(&mat);
  emxFree_int32_T(&c_iidx);
  emxFree_cell_wrap_2(&b_new_PWM);
  emxFree_cell_wrap_14(&b_new_loc);
  emxFree_cell_wrap_2(&new_PWM);
  emxFree_cell_wrap_14(&new_loc);
  emxFree_real_T(&ord);
  emxFree_cell_wrap_14(&poscell);
  emxFree_cell_wrap_14(&loc);
  emxFree_real_T(&lenvec);
  emxFree_real_T(&diffC);

  /* 'gkmPWM:553' toc */
  toc();

  /* 'gkmPWM:554' fprintf('gkmPWM completed after %d iterations\n', int32(i)); */
  printf("gkmPWM completed after %d iterations\n", (int)a);
  fflush(stdout);

  /* 'gkmPWM:556' for i = 1:length(PWM) */
  i = c_PWM->size[0];
  for (acount = 0; acount < i; acount++) {
    /* 'gkmPWM:557' PWM{i} = PWM{i}(l_svm:(length(PWM{i})-l_svm+1),:); */
    u0 = b_PWM_data[acount].f1->size[0];
    if (u0 < 4) {
      u0 = 4;
    }

    if (b_PWM_data[acount].f1->size[0] == 0) {
      u0 = 0;
    }

    S = ((double)u0 - l_svm) + 1.0;
    if (l_svm > S) {
      i1 = -1;
      i2 = -1;
    } else {
      i1 = (int)l_svm - 2;
      i2 = (int)S - 1;
    }

    nx = i2 - i1;
    for (i2 = 0; i2 < 4; i2++) {
      for (i3 = 0; i3 < nx; i3++) {
        b_PWM_data[acount].f1->data[i3 + nx * i2] = b_PWM_data[acount].f1->data
          [((i1 + i3) + b_PWM_data[acount].f1->size[0] * i2) + 1];
      }
    }

    i1 = b_PWM_data[acount].f1->size[0] * b_PWM_data[acount].f1->size[1];
    b_PWM_data[acount].f1->size[0] = nx;
    b_PWM_data[acount].f1->size[1] = 4;
    emxEnsureCapacity_real_T(b_PWM_data[acount].f1, i1);
  }

  /* the following just calculates a few interesting quantities */
  /*  r = corr(kweig, kmat*C); */
  /* 'gkmPWM:561' r = corrcoef(kweig, kmat*C); */
  if ((kmat->size[0] == 0) || (kmat->size[1] == 0) || (C->size[0] == 0)) {
    i = b_iidx->size[0];
    b_iidx->size[0] = kmat->size[0];
    emxEnsureCapacity_real_T(b_iidx, i);
    b_iidx_data = b_iidx->data;
    loop_ub = kmat->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_iidx_data[i] = 0.0;
    }
  } else {
    i = b_iidx->size[0];
    b_iidx->size[0] = kmat->size[0];
    emxEnsureCapacity_real_T(b_iidx, i);
    b_iidx_data = b_iidx->data;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)kmat->size[0],
                (blasint)1, (blasint)kmat->size[1], 1.0, &kmat_data[0], (blasint)
                kmat->size[0], &C_data[0], (blasint)C->size[0], 0.0,
                &b_iidx_data[0], (blasint)kmat->size[0]);
  }

  /* 'gkmPWM:562' r = r(1,2); */
  corrcoef(kweig, b_iidx, varargin_1);
  for (i = 0; i < 4; i++) {
    b[i] = varargin_1[i];
  }

  /* 'gkmPWM:563' M = mean(kweig); */
  M = blockedSummation(kweig, kweig->size[0]) / (double)kweig->size[0];

  /* 'gkmPWM:564' S = std(kweig); */
  S = b_std(kweig);

  /* 'gkmPWM:565' R = zeros(m,1); */
  /* 'gkmPWM:566' E = zeros(m,1); */
  /* 'gkmPWM:567' CM = corrcoef(kmat)-eye(m); */
  /* 'gkmPWM:568' for i = 1:m */
  i = PWM->size[0];
  i1 = R->size[0];
  R->size[0] = PWM->size[0];
  emxEnsureCapacity_real_T(R, i1);
  diffC_data = R->data;
  i1 = E->size[0];
  E->size[0] = PWM->size[0];
  emxEnsureCapacity_real_T(E, i1);
  mat_data = E->data;
  if (0 <= PWM->size[0] - 1) {
    i4 = kmat->size[0];
    c_loop_ub = kmat->size[0];
    if (1 > diffc->size[0]) {
      d_loop_ub = 0;
    } else {
      d_loop_ub = diffc->size[0];
    }

    e_loop_ub = kmat->size[0] * kmat->size[1];
  }

  for (acount = 0; acount < i; acount++) {
    /* 'gkmPWM:569' [~,a] = sort(kmat(:,i),'descend'); */
    i1 = f->size[0];
    f->size[0] = i4;
    emxEnsureCapacity_real_T(f, i1);
    f_data = f->data;
    for (i1 = 0; i1 < c_loop_ub; i1++) {
      f_data[i1] = kmat_data[i1 + kmat->size[0] * acount];
    }

    sort(f, iidx);
    iidx_data = iidx->data;
    i1 = f->size[0];
    f->size[0] = iidx->size[0];
    emxEnsureCapacity_real_T(f, i1);
    f_data = f->data;
    loop_ub = iidx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[i1] = iidx_data[i1];
    }

    /* 'gkmPWM:570' R(i) = (mean(kweig(a(1:lcomb)))-M)/S; */
    i1 = b_iidx->size[0];
    b_iidx->size[0] = d_loop_ub;
    emxEnsureCapacity_real_T(b_iidx, i1);
    b_iidx_data = b_iidx->data;
    for (i1 = 0; i1 < d_loop_ub; i1++) {
      b_iidx_data[i1] = kweig_data[(int)f_data[i1] - 1];
    }

    diffC_data[acount] = (blockedSummation(b_iidx, d_loop_ub) / (double)
                          d_loop_ub - M) / S;

    /* 'gkmPWM:571' Kmat = kmat; */
    /* 'gkmPWM:572' Kmat(:,i) = []; */
    i1 = ct->size[0] * ct->size[1];
    ct->size[0] = kmat->size[0];
    ct->size[1] = kmat->size[1];
    emxEnsureCapacity_real_T(ct, i1);
    ct_data = ct->data;
    for (i1 = 0; i1 < e_loop_ub; i1++) {
      ct_data[i1] = kmat_data[i1];
    }

    b_nullAssignment(ct, acount + 1);
    ct_data = ct->data;

    /* 'gkmPWM:573' c = (Kmat'*Kmat)^(-1)*(Kmat'*kweig); */
    b_mtimes(ct, ct, c_ct);
    if ((ct->size[0] == 0) || (ct->size[1] == 0) || (kweig->size[0] == 0)) {
      i1 = b_iidx->size[0];
      b_iidx->size[0] = ct->size[1];
      emxEnsureCapacity_real_T(b_iidx, i1);
      b_iidx_data = b_iidx->data;
      loop_ub = ct->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_iidx_data[i1] = 0.0;
      }
    } else {
      i1 = b_iidx->size[0];
      b_iidx->size[0] = ct->size[1];
      emxEnsureCapacity_real_T(b_iidx, i1);
      b_iidx_data = b_iidx->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)ct->size[1],
                  (blasint)1, (blasint)ct->size[0], 1.0, &ct_data[0], (blasint)
                  ct->size[0], &kweig_data[0], (blasint)kweig->size[0], 0.0,
                  &b_iidx_data[0], (blasint)ct->size[1]);
    }

    mpower(c_ct, b_a);
    temp_data = b_a->data;
    if ((b_a->size[0] == 0) || (b_a->size[1] == 0) || (b_iidx->size[0] == 0)) {
      i1 = f->size[0];
      f->size[0] = b_a->size[0];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      loop_ub = b_a->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = 0.0;
      }
    } else {
      i1 = f->size[0];
      f->size[0] = b_a->size[0];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)b_a->size
                  [0], (blasint)1, (blasint)b_a->size[1], 1.0, &temp_data[0],
                  (blasint)b_a->size[0], &b_iidx_data[0], (blasint)b_iidx->size
                  [0], 0.0, &f_data[0], (blasint)b_a->size[0]);
    }

    /* 'gkmPWM:574' res = kweig-Kmat*c; */
    if ((ct->size[0] == 0) || (ct->size[1] == 0) || (f->size[0] == 0)) {
      i1 = res->size[0];
      res->size[0] = ct->size[0];
      emxEnsureCapacity_real_T(res, i1);
      res_data = res->data;
      loop_ub = ct->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        res_data[i1] = 0.0;
      }
    } else {
      i1 = res->size[0];
      res->size[0] = ct->size[0];
      emxEnsureCapacity_real_T(res, i1);
      res_data = res->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)ct->size[0],
                  (blasint)1, (blasint)ct->size[1], 1.0, &ct_data[0], (blasint)
                  ct->size[0], &f_data[0], (blasint)f->size[0], 0.0, &res_data[0],
                  (blasint)ct->size[0]);
    }

    loop_ub = kweig->size[0];
    if (kweig->size[0] == res->size[0]) {
      i1 = res->size[0];
      res->size[0] = kweig->size[0];
      emxEnsureCapacity_real_T(res, i1);
      res_data = res->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        res_data[i1] = kweig_data[i1] - res_data[i1];
      }
    } else {
      minus(res, kweig);
      res_data = res->data;
    }

    /* 'gkmPWM:575' E(i) = (sqrt(res'*res)-scorevec(end))/scorevec(end); */
    if (res->size[0] < 1) {
      c = 0.0;
    } else {
      c = cblas_ddot((blasint)res->size[0], &res_data[0], (blasint)1, &res_data
                     [0], (blasint)1);
    }

    mat_data[acount] = (sqrt(c) - scorevec_data[scorevec->size[1] - 1]) /
      scorevec_data[scorevec->size[1] - 1];
  }

  emxFree_real_T(&b_a);
  emxFree_real_T(&res);

  /* 'gkmPWM:577' [R,a] = sort(R, 'descend'); */
  sort(R, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_real_T(f, i);
  f_data = f->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = iidx_data[i];
  }

  emxFree_int32_T(&iidx);

  /*  PWM = PWM(a); */
  /* 'gkmPWM:579' a_len = length(a); */
  /* 'gkmPWM:580' new_PWM = cell(a_len,1); */
  /* 'gkmPWM:581' for cur_idx=1:a_len */
  i = f->size[0];
  i1 = b_PWM->size[0];
  b_PWM->size[0] = f->size[0];
  emxEnsureCapacity_cell_wrap_2(b_PWM, i1);
  new_PWM_data = b_PWM->data;
  for (u0 = 0; u0 < i; u0++) {
    /* 'gkmPWM:582' new_PWM{cur_idx} = PWM{a(cur_idx)}; */
    i1 = new_PWM_data[u0].f1->size[0] * new_PWM_data[u0].f1->size[1];
    new_PWM_data[u0].f1->size[0] = b_PWM_data[(int)f_data[u0] - 1].f1->size[0];
    new_PWM_data[u0].f1->size[1] = 4;
    emxEnsureCapacity_real_T(new_PWM_data[u0].f1, i1);
    loop_ub = b_PWM_data[(int)f_data[u0] - 1].f1->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      new_PWM_data[u0].f1->data[i1] = b_PWM_data[(int)f_data[u0] - 1].f1->
        data[i1];
    }
  }

  emxFree_cell_wrap_2(&c_PWM);

  /* 'gkmPWM:584' PWM = new_PWM; */
  /* 'gkmPWM:586' C = C(a); */
  i = b_iidx->size[0];
  b_iidx->size[0] = f->size[0];
  emxEnsureCapacity_real_T(b_iidx, i);
  b_iidx_data = b_iidx->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_iidx_data[i] = C_data[(int)f_data[i] - 1];
  }

  i = C->size[0];
  C->size[0] = b_iidx->size[0];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  loop_ub = b_iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    C_data[i] = b_iidx_data[i];
  }

  /* 'gkmPWM:587' C = C/max(C); */
  M = maximum(C);
  loop_ub = C->size[0];
  for (i = 0; i < loop_ub; i++) {
    C_data[i] /= M;
  }

  /* 'gkmPWM:588' E = E(a); */
  i = b_iidx->size[0];
  b_iidx->size[0] = f->size[0];
  emxEnsureCapacity_real_T(b_iidx, i);
  b_iidx_data = b_iidx->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_iidx_data[i] = mat_data[(int)f_data[i] - 1];
  }

  i = E->size[0];
  E->size[0] = b_iidx->size[0];
  emxEnsureCapacity_real_T(E, i);
  mat_data = E->data;
  loop_ub = b_iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    mat_data[i] = b_iidx_data[i];
  }

  /* 'gkmPWM:589' E = E/max(abs(E)); */
  nx = E->size[0];
  i = b_iidx->size[0];
  b_iidx->size[0] = E->size[0];
  emxEnsureCapacity_real_T(b_iidx, i);
  b_iidx_data = b_iidx->data;
  for (u0 = 0; u0 < nx; u0++) {
    b_iidx_data[u0] = fabs(mat_data[u0]);
  }

  M = maximum(b_iidx);
  loop_ub = E->size[0];
  emxFree_real_T(&b_iidx);
  for (i = 0; i < loop_ub; i++) {
    mat_data[i] /= M;
  }

  /* 'gkmPWM:590' Rd = max(CM); */
  nx = PWM->size[0];
  i = ct->size[0] * ct->size[1];
  ct->size[0] = PWM->size[0];
  ct->size[1] = PWM->size[0];
  emxEnsureCapacity_real_T(ct, i);
  ct_data = ct->data;
  loop_ub = PWM->size[0] * PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    ct_data[i] = 0.0;
  }

  for (u0 = 0; u0 < nx; u0++) {
    ct_data[u0 + ct->size[0] * u0] = 1.0;
  }

  b_corrcoef(kmat, CT);
  CT_data = CT->data;
  emxFree_real_T(&kmat);
  if ((CT->size[0] == ct->size[0]) && (CT->size[1] == ct->size[1])) {
    i = c_ct->size[0] * c_ct->size[1];
    c_ct->size[0] = CT->size[0];
    c_ct->size[1] = CT->size[1];
    emxEnsureCapacity_real_T(c_ct, i);
    mat_data = c_ct->data;
    loop_ub = CT->size[0] * CT->size[1];
    for (i = 0; i < loop_ub; i++) {
      mat_data[i] = CT_data[i] - ct_data[i];
    }

    g_maximum(c_ct, Rd);
    diffC_data = Rd->data;
  } else {
    q_binary_expand_op(Rd, CT, ct);
    diffC_data = Rd->data;
  }

  emxFree_real_T(&c_ct);
  emxFree_real_T(&CT);
  emxFree_real_T(&ct);

  /* 'gkmPWM:591' Rd = Rd(a); */
  i = temp->size[0] * temp->size[1];
  temp->size[0] = 1;
  temp->size[1] = f->size[0];
  emxEnsureCapacity_real_T(temp, i);
  temp_data = temp->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    temp_data[i] = diffC_data[(int)f_data[i] - 1];
  }

  emxFree_real_T(&f);
  i = Rd->size[0] * Rd->size[1];
  Rd->size[0] = 1;
  Rd->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(Rd, i);
  diffC_data = Rd->data;
  loop_ub = temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    diffC_data[i] = temp_data[i];
  }

  emxFree_real_T(&temp);
  *r = b[2];
}

static void h_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3)
{
  const cell_wrap_1 *mat_data;
  const double *rs_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  bool *x_data;
  rs_data = rs->data;
  mat_data = mat->data;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  if ((i3 - i2) + 1 == 1) {
    x->size[1] = i1 + 3;
  } else {
    x->size[1] = (i3 - i2) + 1;
  }

  emxEnsureCapacity_boolean_T(x, i);
  x_data = x->data;
  stride_0_1 = (i1 + 3 != 1);
  stride_1_1 = ((i3 - i2) + 1 != 1);
  if ((i3 - i2) + 1 == 1) {
    loop_ub = i1 + 3;
  } else {
    loop_ub = (i3 - i2) + 1;
  }

  for (i = 0; i < loop_ub; i++) {
    x_data[i] = (mat_data[j].f1->data[i * stride_0_1] == rs_data[i2 + i *
                 stride_1_1]);
  }
}

static void i_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3)
{
  const cell_wrap_1 *mat_data;
  const double *rs_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  bool *x_data;
  rs_data = rs->data;
  mat_data = mat->data;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (i3 + 1 == 1) {
    x->size[1] = (i2 - i1) + 1;
  } else {
    x->size[1] = i3 + 1;
  }

  emxEnsureCapacity_boolean_T(x, i);
  x_data = x->data;
  stride_0_1 = ((i2 - i1) + 1 != 1);
  stride_1_1 = (i3 + 1 != 1);
  if (i3 + 1 == 1) {
    loop_ub = (i2 - i1) + 1;
  } else {
    loop_ub = i3 + 1;
  }

  for (i = 0; i < loop_ub; i++) {
    x_data[i] = (mat_data[j].f1->data[i1 + i * stride_0_1] == rs_data[i *
                 stride_1_1]);
  }
}

static void j_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3)
{
  const cell_wrap_1 *mat_data;
  const double *rs_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  bool *x_data;
  rs_data = rs->data;
  mat_data = mat->data;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (i3 + 2 == 1) {
    x->size[1] = (i2 - i1) + 1;
  } else {
    x->size[1] = i3 + 2;
  }

  emxEnsureCapacity_boolean_T(x, i);
  x_data = x->data;
  stride_0_1 = ((i2 - i1) + 1 != 1);
  stride_1_1 = (i3 + 2 != 1);
  if (i3 + 2 == 1) {
    loop_ub = (i2 - i1) + 1;
  } else {
    loop_ub = i3 + 2;
  }

  for (i = 0; i < loop_ub; i++) {
    x_data[i] = (mat_data[j].f1->data[i1 + i * stride_0_1] == rs_data[i *
                 stride_1_1]);
  }
}

static void k_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3)
{
  const cell_wrap_1 *mat_data;
  const double *rs_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  bool *x_data;
  rs_data = rs->data;
  mat_data = mat->data;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (i3 + 3 == 1) {
    x->size[1] = (i2 - i1) + 1;
  } else {
    x->size[1] = i3 + 3;
  }

  emxEnsureCapacity_boolean_T(x, i);
  x_data = x->data;
  stride_0_1 = ((i2 - i1) + 1 != 1);
  stride_1_1 = (i3 + 3 != 1);
  if (i3 + 3 == 1) {
    loop_ub = (i2 - i1) + 1;
  } else {
    loop_ub = i3 + 3;
  }

  for (i = 0; i < loop_ub; i++) {
    x_data[i] = (mat_data[j].f1->data[i1 + i * stride_0_1] == rs_data[i *
                 stride_1_1]);
  }
}

static void l_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, const emxArray_real_T *rs)
{
  const cell_wrap_1 *mat_data;
  const double *rs_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  bool *x_data;
  rs_data = rs->data;
  mat_data = mat->data;
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (rs->size[1] == 1) {
    x->size[1] = mat_data[j].f1->size[1];
  } else {
    x->size[1] = rs->size[1];
  }

  emxEnsureCapacity_boolean_T(x, i);
  x_data = x->data;
  stride_0_1 = (mat_data[j].f1->size[1] != 1);
  stride_1_1 = (rs->size[1] != 1);
  if (rs->size[1] == 1) {
    loop_ub = mat_data[j].f1->size[1];
  } else {
    loop_ub = rs->size[1];
  }

  for (i = 0; i < loop_ub; i++) {
    x_data[i] = (mat_data[j].f1->data[i * stride_0_1] == rs_data[i * stride_1_1]);
  }
}

/*
 * function kweig = ls_kweigtree(mat,negmat,poscell,c,s,ind,indloc,x,l,k,rcnum)
 */
static void ls_kweigtree(const emxArray_real_T *mat, const double negmat[16],
  const emxArray_cell_wrap_14 *poscell, const emxArray_real_T *c, const
  emxArray_real_T *s, const emxArray_real_T *ind, const emxArray_real_T *indloc,
  const emxArray_real_T *x, double l, double k, double rcnum, emxArray_real_T
  *kweig)
{
  cell_wrap_13 *p_data;
  const cell_wrap_14 *poscell_data;
  cell_wrap_14 *ktree2_data;
  cell_wrap_14 *ktree_data;
  emxArray_boolean_T *b_x;
  emxArray_cell_wrap_13 *p;
  emxArray_cell_wrap_14 *ktree;
  emxArray_cell_wrap_14 *ktree2;
  emxArray_int32_T *f1;
  emxArray_real_T *a;
  emxArray_real_T *b_indvec;
  emxArray_real_T *b_kweig;
  emxArray_real_T *indloc2;
  emxArray_real_T *indvec;
  emxArray_real_T *indvec2;
  emxArray_real_T *mat2;
  emxArray_real_T *sPWM;
  emxArray_real_T *sPWM2;
  emxArray_uint32_T *X;
  double b_p[16];
  double c_x[4];
  const double *c_data;
  const double *ind_data;
  const double *indloc_data;
  const double *mat_data;
  const double *s_data;
  const double *x_data;
  double b_c;
  double b_m;
  double b_n_tmp;
  double c_tmp;
  double d;
  double n;
  double n_tmp;
  double *a_data;
  double *b_indvec_data;
  double *indloc2_data;
  double *indvec2_data;
  double *indvec_data;
  double *kweig_data;
  double *mat2_data;
  double *sPWM2_data;
  double *sPWM_data;
  int b_i;
  int b_k;
  int f;
  int i;
  int i1;
  int i2;
  int ii;
  int iii;
  int indloc2_tmp;
  int j;
  int loop_ub_tmp;
  int m;
  int md2;
  int rx;
  unsigned int *X_data;
  int *f1_data;
  bool *b_x_data;
  x_data = x->data;
  indloc_data = indloc->data;
  ind_data = ind->data;
  s_data = s->data;
  c_data = c->data;
  poscell_data = poscell->data;
  mat_data = mat->data;
  emxInit_cell_wrap_13(&p);

  /* uses dynamic programming to find the matrix A (kweig in this code) such that A*p = expected gapped kmer vector, where p is the middle column of a 2*l-1 PWM. */
  /* I use a 1st order Markov model to model the flanking bases for a PWM */
  /* 'gkmPWM:717' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_13(p, i);
  p_data = p->data;

  /* 'gkmPWM:718' p = coder.nullcopy(p); */
  /* 'gkmPWM:719' p{1} = eye(4); */
  for (i = 0; i < 16; i++) {
    p_data[0].f1[i] = 0.0;
  }

  for (b_k = 0; b_k < 4; b_k++) {
    p_data[0].f1[b_k + (b_k << 2)] = 1.0;
  }

  /* 'gkmPWM:720' for i = 1:l-1 */
  i = (int)(l - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:721' p{i+1} = p{i}*negmat; */
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < 4; i2++) {
        d = 0.0;
        for (m = 0; m < 4; m++) {
          d += p_data[b_i].f1[i1 + (m << 2)] * negmat[m + (i2 << 2)];
        }

        b_p[i1 + (i2 << 2)] = d;
      }
    }

    for (i1 = 0; i1 < 16; i1++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i1] = b_p[i1];
    }
  }

  emxInit_real_T(&mat2, 2);
  emxInit_real_T(&indvec, 2);

  /* 'gkmPWM:723' n = 4^k*max(max(x)); */
  g_maximum(x, indvec);
  n_tmp = h_maximum(indvec);
  b_n_tmp = pow(4.0, k);
  n = b_n_tmp * n_tmp;

  /* number of possible k-mers */
  /* 'gkmPWM:724' mat2 = rot90(mat,2); */
  d_rot90(mat, mat2);
  mat2_data = mat2->data;

  /* 'gkmPWM:725' kweig = zeros(n, 4); */
  i = kweig->size[0] * kweig->size[1];
  kweig->size[0] = (int)n;
  kweig->size[1] = 4;
  emxEnsureCapacity_real_T(kweig, i);
  kweig_data = kweig->data;
  md2 = (int)n << 2;
  for (i = 0; i < md2; i++) {
    kweig_data[i] = 0.0;
  }

  emxInit_cell_wrap_14(&ktree);
  emxInit_cell_wrap_14(&ktree2);

  /* 'gkmPWM:726' ktree = cell(k-1,1); */
  m = (int)(k - 1.0);
  i = ktree->size[0];
  ktree->size[0] = (int)(k - 1.0);
  emxEnsureCapacity_cell_wrap_14(ktree, i);
  ktree_data = ktree->data;

  /* 'gkmPWM:727' ktree = coder.nullcopy(ktree); */
  /* 'gkmPWM:728' ktree2 = cell(k-1,1); */
  i = ktree2->size[0];
  ktree2->size[0] = (int)(k - 1.0);
  emxEnsureCapacity_cell_wrap_14(ktree2, i);
  ktree2_data = ktree2->data;
  for (i = 0; i < m; i++) {
    ktree_data[i].f1->size[0] = 0;
    ktree2_data[i].f1->size[0] = 0;
  }

  emxInit_uint32_T(&X);

  /* 'gkmPWM:729' ktree2 = coder.nullcopy(ktree2); */
  /* 'gkmPWM:730' [rx,cx] = size(x); */
  rx = x->size[0];

  /* 'gkmPWM:731' m=rx; */
  b_m = x->size[0];

  /* 'gkmPWM:732' M = l-1; */
  /* 'gkmPWM:733' X = cx*ones(M+1,1); */
  loop_ub_tmp = (int)((l - 1.0) + 1.0);
  i = X->size[0];
  X->size[0] = (int)((l - 1.0) + 1.0);
  emxEnsureCapacity_uint32_T(X, i);
  X_data = X->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    X_data[i] = (unsigned int)x->size[1];
  }

  /* 'gkmPWM:734' for i = 1:cx */
  i = x->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:735' X(i) = i; */
    X_data[b_i] = (unsigned int)(b_i + 1);
  }

  emxInit_real_T(&indloc2, 1);

  /* 'gkmPWM:737' indloc2 = flipud(indloc); */
  i = indloc2->size[0];
  indloc2->size[0] = indloc->size[0];
  emxEnsureCapacity_real_T(indloc2, i);
  indloc2_data = indloc2->data;
  md2 = indloc->size[0];
  for (i = 0; i < md2; i++) {
    indloc2_data[i] = indloc_data[i];
  }

  m = indloc->size[0] - 1;
  md2 = indloc->size[0] >> 1;
  for (b_i = 0; b_i < md2; b_i++) {
    n = indloc2_data[b_i];
    indloc2_tmp = m - b_i;
    indloc2_data[b_i] = indloc2_data[indloc2_tmp];
    indloc2_data[indloc2_tmp] = n;
  }

  /* 'gkmPWM:738' for i = 2:5 */
  for (b_i = 0; b_i < 4; b_i++) {
    /* 'gkmPWM:739' ktree{i} = zeros(4^i,1); */
    c_tmp = pow(4.0, (double)b_i + 2.0);
    m = (int)pow(4.0, (double)b_i + 2.0);
    i = ktree_data[b_i + 1].f1->size[0];
    ktree_data[b_i + 1].f1->size[0] = (int)c_tmp;
    emxEnsureCapacity_real_T(ktree_data[b_i + 1].f1, i);
    for (i = 0; i < m; i++) {
      ktree_data[b_i + 1].f1->data[i] = 0.0;
    }

    /* 'gkmPWM:740' ktree2{i} = zeros(4^i,1); */
    i = ktree2_data[b_i + 1].f1->size[0];
    ktree2_data[b_i + 1].f1->size[0] = (int)c_tmp;
    emxEnsureCapacity_real_T(ktree2_data[b_i + 1].f1, i);
    for (i = 0; i < m; i++) {
      ktree2_data[b_i + 1].f1->data[i] = 0.0;
    }
  }

  /* 'gkmPWM:742' for i = 0:M */
  emxInit_int32_T(&f1, 2);
  emxInit_real_T(&sPWM, 2);
  emxInit_real_T(&sPWM2, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&indvec2, 1);
  emxInit_real_T(&b_indvec, 1);
  emxInit_boolean_T(&b_x, 2);
  emxInit_real_T(&b_kweig, 1);
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    /* 'gkmPWM:743' if i > M-cx+1 */
    if (b_i > ((l - 1.0) - (double)x->size[1]) + 1.0) {
      /* 'gkmPWM:744' m = numel(c)/k; */
      b_m = (double)(c->size[0] * c->size[1]) / k;
    }

    /* the following loops is basically dynamic programming for tensor multiplication.  there are multiple cases to consider, hence the if statements. */
    /* 'gkmPWM:747' for ii = 1:m */
    i = (int)b_m;
    for (ii = 0; ii < i; ii++) {
      /* 'gkmPWM:748' if sum((c(ii,:)+i)==l) > 0 && ~(i == M-1 && ii > rx && ii ~= m) */
      md2 = c->size[1];
      i1 = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 1;
      b_x->size[1] = c->size[1];
      emxEnsureCapacity_boolean_T(b_x, i1);
      b_x_data = b_x->data;
      for (i1 = 0; i1 < md2; i1++) {
        b_x_data[i1] = (c_data[ii + c->size[0] * i1] + (double)b_i == l);
      }

      m = b_x->size[1];
      if (b_x->size[1] == 0) {
        md2 = 0;
      } else {
        md2 = b_x_data[0];
        for (b_k = 2; b_k <= m; b_k++) {
          md2 += b_x_data[b_k - 1];
        }
      }

      if ((md2 > 0) && ((b_i != (l - 1.0) - 1.0) || (ii + 1U <= (unsigned int)rx)
                        || ((double)ii + 1.0 == b_m))) {
        /* 'gkmPWM:749' indvec = c(ii,:)+i; */
        md2 = c->size[1];
        i1 = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i1);
        indvec_data = indvec->data;
        for (i1 = 0; i1 < md2; i1++) {
          indvec_data[i1] = c_data[ii + c->size[0] * i1] + (double)b_i;
        }

        /*  f = find(indvec == l); */
        /* 'gkmPWM:751' coder.varsize('f', [1 1]); */
        /* 'gkmPWM:752' f1 = find(indvec == l); */
        i1 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = indvec->size[1];
        emxEnsureCapacity_boolean_T(b_x, i1);
        b_x_data = b_x->data;
        md2 = indvec->size[1];
        for (i1 = 0; i1 < md2; i1++) {
          b_x_data[i1] = (indvec_data[i1] == l);
        }

        eml_find(b_x, f1);
        f1_data = f1->data;

        /* 'gkmPWM:753' if length(f1) == 0 */
        if (f1->size[1] == 0) {
          /* 'gkmPWM:754' f = zeros(1,1); */
          f = -1;

          /* 'gkmPWM:755' fprintf("Resizing f\n"); */
          printf("Resizing f\n");
          fflush(stdout);
        } else {
          /* 'gkmPWM:756' else */
          /* 'gkmPWM:757' f = f1(1); */
          f = f1_data[0] - 1;
        }

        /* 'gkmPWM:759' indvec(f) = []; */
        indloc2_tmp = f + 1;
        m = indvec->size[1];
        md2 = indvec->size[1] - 1;
        for (b_k = indloc2_tmp; b_k <= md2; b_k++) {
          indvec_data[b_k - 1] = indvec_data[b_k];
        }

        i1 = indvec->size[0] * indvec->size[1];
        if (1 > md2) {
          indvec->size[1] = 0;
        } else {
          indvec->size[1] = m - 1;
        }

        emxEnsureCapacity_real_T(indvec, i1);
        indvec_data = indvec->data;

        /* 'gkmPWM:760' loc = indloc(indvec); */
        /* 'gkmPWM:761' loc2 = indloc2(indvec); */
        /* 'gkmPWM:762' sPWM = mat(indvec,:).'; */
        i1 = b_indvec->size[0];
        b_indvec->size[0] = indvec->size[1];
        emxEnsureCapacity_real_T(b_indvec, i1);
        b_indvec_data = b_indvec->data;
        md2 = indvec->size[1];
        for (i1 = 0; i1 < md2; i1++) {
          b_indvec_data[i1] = indvec_data[i1];
        }

        i1 = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = 4;
        sPWM->size[1] = b_indvec->size[0];
        emxEnsureCapacity_real_T(sPWM, i1);
        sPWM_data = sPWM->data;
        md2 = b_indvec->size[0];
        for (i1 = 0; i1 < md2; i1++) {
          for (i2 = 0; i2 < 4; i2++) {
            sPWM_data[i2 + 4 * i1] = mat_data[((int)b_indvec_data[i1] +
              mat->size[0] * i2) - 1];
          }
        }

        /* 'gkmPWM:763' sPWM2 = mat2(indvec,:).'; */
        i1 = sPWM2->size[0] * sPWM2->size[1];
        sPWM2->size[0] = 4;
        sPWM2->size[1] = b_indvec->size[0];
        emxEnsureCapacity_real_T(sPWM2, i1);
        sPWM2_data = sPWM2->data;
        md2 = b_indvec->size[0];
        for (i1 = 0; i1 < md2; i1++) {
          for (i2 = 0; i2 < 4; i2++) {
            sPWM2_data[i2 + 4 * i1] = mat2_data[((int)b_indvec_data[i1] +
              mat2->size[0] * i2) - 1];
          }
        }

        /* 'gkmPWM:764' ktree{1} = sPWM(:,1); */
        i1 = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree_data[0].f1, i1);

        /* 'gkmPWM:765' ktree2{1} = sPWM2(:,1); */
        i1 = ktree2_data[0].f1->size[0];
        ktree2_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree2_data[0].f1, i1);
        for (i1 = 0; i1 < 4; i1++) {
          ktree_data[0].f1->data[i1] = sPWM_data[i1];
          ktree2_data[0].f1->data[i1] = sPWM2_data[i1];
        }

        /* 'gkmPWM:766' for iii = s(ii,i+1):k-1 */
        d = s_data[ii + s->size[0] * b_i];
        i1 = (int)((k - 1.0) + (1.0 - d));
        for (iii = 0; iii < i1; iii++) {
          n = d + (double)iii;

          /* 'gkmPWM:767' if loc(iii)==0 */
          i2 = (int)indvec_data[(int)n - 1] - 1;
          if (indloc_data[i2] == 0.0) {
            /* 'gkmPWM:768' if loc(iii-1)==1 && indvec(iii-1) < l */
            if ((indloc_data[(int)indvec_data[(int)(n - 1.0) - 1] - 1] == 1.0) &&
                ((unsigned int)indvec_data[(int)(n - 1.0) - 1] < l)) {
              /* 'gkmPWM:769' matt = sPWM(:,iii).'*p{indvec(iii)-l}; */
              for (i2 = 0; i2 < 16; i2++) {
                b_p[i2] = p_data[(int)((double)(unsigned int)indvec_data[(int)n
                  - 1] - l) - 1].f1[i2];
              }

              /* 'gkmPWM:770' a = ktree2{iii-1}.*sPWM2(:,iii).'; */
              i2 = a->size[0] * a->size[1];
              a->size[0] = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, i2);
              a_data = a->data;
              md2 = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
              for (i2 = 0; i2 < 4; i2++) {
                for (m = 0; m < md2; m++) {
                  a_data[m + a->size[0] * i2] = ktree2_data[(int)(n - 1.0) - 1].
                    f1->data[m] * sPWM2_data[i2 + 4 * ((int)n - 1)];
                }
              }

              /* 'gkmPWM:771' ktree2{iii} = a(:); */
              i2 = ktree2_data[(int)n - 1].f1->size[0];
              ktree2_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, i2);
              md2 = a->size[0] << 2;
              for (i2 = 0; i2 < md2; i2++) {
                ktree2_data[(int)n - 1].f1->data[i2] = a_data[i2];
              }

              /* 'gkmPWM:772' ktree{iii} = repmat(ktree{iii-1},4,1).*repelem(matt', 4^(iii-1)); */
              c_tmp = pow(4.0, n - 1.0);
              i2 = b_indvec->size[0];
              b_indvec->size[0] = (int)c_tmp << 2;
              emxEnsureCapacity_real_T(b_indvec, i2);
              b_indvec_data = b_indvec->data;
              indloc2_tmp = -1;
              m = (int)c_tmp;
              for (b_k = 0; b_k < 4; b_k++) {
                b_c = 0.0;
                for (i2 = 0; i2 < 4; i2++) {
                  b_c += sPWM_data[i2 + 4 * ((int)n - 1)] * b_p[i2 + (b_k << 2)];
                }

                c_x[b_k] = b_c;
                for (j = 0; j < m; j++) {
                  b_indvec_data[(indloc2_tmp + j) + 1] = c_x[b_k];
                }

                indloc2_tmp += (int)c_tmp;
              }

              c_repmat(ktree_data[(int)(n - 1.0) - 1].f1, indvec2);
              indvec2_data = indvec2->data;
              if (indvec2->size[0] == b_indvec->size[0]) {
                i2 = ktree_data[(int)n - 1].f1->size[0];
                ktree_data[(int)n - 1].f1->size[0] = indvec2->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, i2);
                md2 = indvec2->size[0];
                for (i2 = 0; i2 < md2; i2++) {
                  ktree_data[(int)n - 1].f1->data[i2] = indvec2_data[i2] *
                    b_indvec_data[i2];
                }
              } else {
                r_binary_expand_op(ktree, n, indvec2, b_indvec);
                ktree_data = ktree->data;
              }
            } else {
              /* 'gkmPWM:773' else */
              /* 'gkmPWM:774' matt = p{indvec(iii)-indvec(iii-1)+1}; */
              /* 'gkmPWM:775' ktree{iii} = repmat(ktree{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2)); */
              c_tmp = pow(4.0, n - 2.0);
              i2 = b_indvec->size[0];
              b_indvec->size[0] = (int)c_tmp << 4;
              emxEnsureCapacity_real_T(b_indvec, i2);
              b_indvec_data = b_indvec->data;
              indloc2_tmp = -1;
              m = (int)c_tmp;
              for (b_k = 0; b_k < 16; b_k++) {
                for (j = 0; j < m; j++) {
                  b_indvec_data[(indloc2_tmp + j) + 1] = p_data[(int)(unsigned
                    int)indvec_data[(int)n - 1] - (int)(unsigned int)
                    indvec_data[(int)(n - 1.0) - 1]].f1[b_k];
                }

                indloc2_tmp += (int)c_tmp;
              }

              c_repmat(ktree_data[(int)(n - 1.0) - 1].f1, indvec2);
              indvec2_data = indvec2->data;
              if (indvec2->size[0] == b_indvec->size[0]) {
                i2 = ktree_data[(int)n - 1].f1->size[0];
                ktree_data[(int)n - 1].f1->size[0] = indvec2->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, i2);
                md2 = indvec2->size[0];
                for (i2 = 0; i2 < md2; i2++) {
                  ktree_data[(int)n - 1].f1->data[i2] = indvec2_data[i2] *
                    b_indvec_data[i2];
                }
              } else {
                r_binary_expand_op(ktree, n, indvec2, b_indvec);
                ktree_data = ktree->data;
              }

              /* 'gkmPWM:776' a = ktree2{iii-1}.*sPWM2(:,iii).'; */
              i2 = a->size[0] * a->size[1];
              a->size[0] = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, i2);
              a_data = a->data;
              md2 = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
              for (i2 = 0; i2 < 4; i2++) {
                for (m = 0; m < md2; m++) {
                  a_data[m + a->size[0] * i2] = ktree2_data[(int)(n - 1.0) - 1].
                    f1->data[m] * sPWM2_data[i2 + 4 * ((int)n - 1)];
                }
              }

              /* 'gkmPWM:777' ktree2{iii} = a(:); */
              i2 = ktree2_data[(int)n - 1].f1->size[0];
              ktree2_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, i2);
              md2 = a->size[0] << 2;
              for (i2 = 0; i2 < md2; i2++) {
                ktree2_data[(int)n - 1].f1->data[i2] = a_data[i2];
              }
            }
          } else if (indloc2_data[i2] == 0.0) {
            /* 'gkmPWM:779' elseif loc2(iii)==0 */
            /* 'gkmPWM:780' if loc2(iii-1)==1 && indvec(iii-1) < l */
            if ((indloc2_data[(int)indvec_data[(int)(n - 1.0) - 1] - 1] == 1.0) &&
                ((unsigned int)indvec_data[(int)(n - 1.0) - 1] < l)) {
              /* 'gkmPWM:781' matt = sPWM2(:,iii).'*p{indvec(iii)-l}; */
              for (i2 = 0; i2 < 16; i2++) {
                b_p[i2] = p_data[(int)((double)(unsigned int)indvec_data[(int)n
                  - 1] - l) - 1].f1[i2];
              }

              /* 'gkmPWM:782' a = ktree{iii-1}.*sPWM(:,iii).'; */
              i2 = a->size[0] * a->size[1];
              a->size[0] = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, i2);
              a_data = a->data;
              md2 = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
              for (i2 = 0; i2 < 4; i2++) {
                for (m = 0; m < md2; m++) {
                  a_data[m + a->size[0] * i2] = ktree_data[(int)(n - 1.0) - 1].
                    f1->data[m] * sPWM_data[i2 + 4 * ((int)n - 1)];
                }
              }

              /* 'gkmPWM:783' ktree{iii} = a(:); */
              i2 = ktree_data[(int)n - 1].f1->size[0];
              ktree_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, i2);
              md2 = a->size[0] << 2;
              for (i2 = 0; i2 < md2; i2++) {
                ktree_data[(int)n - 1].f1->data[i2] = a_data[i2];
              }

              /* 'gkmPWM:784' ktree2{iii} = repmat(ktree2{iii-1},4,1).*repelem(matt', 4^(iii-1)); */
              c_tmp = pow(4.0, n - 1.0);
              i2 = b_indvec->size[0];
              b_indvec->size[0] = (int)c_tmp << 2;
              emxEnsureCapacity_real_T(b_indvec, i2);
              b_indvec_data = b_indvec->data;
              indloc2_tmp = -1;
              m = (int)c_tmp;
              for (b_k = 0; b_k < 4; b_k++) {
                b_c = 0.0;
                for (i2 = 0; i2 < 4; i2++) {
                  b_c += sPWM2_data[i2 + 4 * ((int)n - 1)] * b_p[i2 + (b_k << 2)];
                }

                c_x[b_k] = b_c;
                for (j = 0; j < m; j++) {
                  b_indvec_data[(indloc2_tmp + j) + 1] = c_x[b_k];
                }

                indloc2_tmp += (int)c_tmp;
              }

              c_repmat(ktree2_data[(int)(n - 1.0) - 1].f1, indvec2);
              indvec2_data = indvec2->data;
              if (indvec2->size[0] == b_indvec->size[0]) {
                i2 = ktree2_data[(int)n - 1].f1->size[0];
                ktree2_data[(int)n - 1].f1->size[0] = indvec2->size[0];
                emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, i2);
                md2 = indvec2->size[0];
                for (i2 = 0; i2 < md2; i2++) {
                  ktree2_data[(int)n - 1].f1->data[i2] = indvec2_data[i2] *
                    b_indvec_data[i2];
                }
              } else {
                r_binary_expand_op(ktree2, n, indvec2, b_indvec);
                ktree2_data = ktree2->data;
              }
            } else {
              /* 'gkmPWM:785' else */
              /* 'gkmPWM:786' matt = p{indvec(iii)-indvec(iii-1)+1}; */
              /* 'gkmPWM:787' ktree2{iii} = repmat(ktree2{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2)); */
              c_tmp = pow(4.0, n - 2.0);
              i2 = b_indvec->size[0];
              b_indvec->size[0] = (int)c_tmp << 4;
              emxEnsureCapacity_real_T(b_indvec, i2);
              b_indvec_data = b_indvec->data;
              indloc2_tmp = -1;
              m = (int)c_tmp;
              for (b_k = 0; b_k < 16; b_k++) {
                for (j = 0; j < m; j++) {
                  b_indvec_data[(indloc2_tmp + j) + 1] = p_data[(int)(unsigned
                    int)indvec_data[(int)n - 1] - (int)(unsigned int)
                    indvec_data[(int)(n - 1.0) - 1]].f1[b_k];
                }

                indloc2_tmp += (int)c_tmp;
              }

              c_repmat(ktree2_data[(int)(n - 1.0) - 1].f1, indvec2);
              indvec2_data = indvec2->data;
              if (indvec2->size[0] == b_indvec->size[0]) {
                i2 = ktree2_data[(int)n - 1].f1->size[0];
                ktree2_data[(int)n - 1].f1->size[0] = indvec2->size[0];
                emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, i2);
                md2 = indvec2->size[0];
                for (i2 = 0; i2 < md2; i2++) {
                  ktree2_data[(int)n - 1].f1->data[i2] = indvec2_data[i2] *
                    b_indvec_data[i2];
                }
              } else {
                r_binary_expand_op(ktree2, n, indvec2, b_indvec);
                ktree2_data = ktree2->data;
              }

              /* 'gkmPWM:788' a = ktree{iii-1}.*sPWM(:,iii).'; */
              i2 = a->size[0] * a->size[1];
              a->size[0] = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, i2);
              a_data = a->data;
              md2 = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
              for (i2 = 0; i2 < 4; i2++) {
                for (m = 0; m < md2; m++) {
                  a_data[m + a->size[0] * i2] = ktree_data[(int)(n - 1.0) - 1].
                    f1->data[m] * sPWM_data[i2 + 4 * ((int)n - 1)];
                }
              }

              /* 'gkmPWM:789' ktree{iii} = a(:); */
              i2 = ktree_data[(int)n - 1].f1->size[0];
              ktree_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, i2);
              md2 = a->size[0] << 2;
              for (i2 = 0; i2 < md2; i2++) {
                ktree_data[(int)n - 1].f1->data[i2] = a_data[i2];
              }
            }
          } else {
            /* 'gkmPWM:791' else */
            /* 'gkmPWM:792' a = ktree{iii-1}.*sPWM(:,iii).'; */
            i2 = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, i2);
            a_data = a->data;
            md2 = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
            for (i2 = 0; i2 < 4; i2++) {
              for (m = 0; m < md2; m++) {
                a_data[m + a->size[0] * i2] = ktree_data[(int)(n - 1.0) - 1].
                  f1->data[m] * sPWM_data[i2 + 4 * ((int)n - 1)];
              }
            }

            /* 'gkmPWM:793' ktree{iii} = a(:); */
            i2 = ktree_data[(int)n - 1].f1->size[0];
            ktree_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, i2);
            md2 = a->size[0] << 2;
            for (i2 = 0; i2 < md2; i2++) {
              ktree_data[(int)n - 1].f1->data[i2] = a_data[i2];
            }

            /* 'gkmPWM:794' a = ktree2{iii-1}.*sPWM2(:,iii).'; */
            i2 = a->size[0] * a->size[1];
            a->size[0] = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, i2);
            a_data = a->data;
            md2 = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
            for (i2 = 0; i2 < 4; i2++) {
              for (m = 0; m < md2; m++) {
                a_data[m + a->size[0] * i2] = ktree2_data[(int)(n - 1.0) - 1].
                  f1->data[m] * sPWM2_data[i2 + 4 * ((int)n - 1)];
              }
            }

            /* 'gkmPWM:795' ktree2{iii} = a(:); */
            i2 = ktree2_data[(int)n - 1].f1->size[0];
            ktree2_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, i2);
            md2 = a->size[0] << 2;
            for (i2 = 0; i2 < md2; i2++) {
              ktree2_data[(int)n - 1].f1->data[i2] = a_data[i2];
            }
          }
        }

        /* the weird indexing that I did early in the code comes to fruition.  It is critical to do so to make this computation as fast as possible. */
        /* 'gkmPWM:799' if ii <= rx */
        if (ii + 1U <= (unsigned int)rx) {
          /* 'gkmPWM:800' for j = 1:X(i+1) */
          i1 = (int)X_data[b_i];
          for (j = 0; j < i1; j++) {
            /* 'gkmPWM:801' if x(ii,j) ~= 0 */
            d = x_data[ii + x->size[0] * j];
            if (d != 0.0) {
              /* 'gkmPWM:802' for iii = 1:2 */
              c_tmp = pow(4.0, (double)(f + 1) - 1.0);
              b_c = b_n_tmp * (ind_data[(int)d - 1] - 1.0);
              md2 = poscell_data[f].f1->size[0];
              for (iii = 0; iii < 2; iii++) {
                /* 'gkmPWM:803' indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(x(ii,j))-1); */
                n = c_tmp * (((double)iii + 1.0) - 1.0);
                i2 = b_indvec->size[0];
                b_indvec->size[0] = poscell_data[f].f1->size[0];
                emxEnsureCapacity_real_T(b_indvec, i2);
                b_indvec_data = b_indvec->data;
                for (i2 = 0; i2 < md2; i2++) {
                  b_indvec_data[i2] = (poscell_data[f].f1->data[i2] + n) + b_c;
                }

                /* 'gkmPWM:804' indvec2 = indvec+(5-2*iii)*4^(f-1); */
                d = (5.0 - 2.0 * ((double)iii + 1.0)) * c_tmp;
                i2 = indvec2->size[0];
                indvec2->size[0] = b_indvec->size[0];
                emxEnsureCapacity_real_T(indvec2, i2);
                indvec2_data = indvec2->data;
                m = b_indvec->size[0];
                for (i2 = 0; i2 < m; i2++) {
                  indvec2_data[i2] = b_indvec_data[i2] + d;
                }

                /* 'gkmPWM:805' kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1}; */
                if (b_indvec->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->
                    size[0]) {
                  i2 = b_kweig->size[0];
                  b_kweig->size[0] = b_indvec->size[0];
                  emxEnsureCapacity_real_T(b_kweig, i2);
                  a_data = b_kweig->data;
                  m = b_indvec->size[0];
                  for (i2 = 0; i2 < m; i2++) {
                    a_data[i2] = kweig_data[((int)b_indvec_data[i2] +
                      kweig->size[0] * iii) - 1] + ktree_data[(int)(k - 1.0) - 1]
                      .f1->data[i2];
                  }

                  m = b_kweig->size[0];
                  for (i2 = 0; i2 < m; i2++) {
                    kweig_data[((int)b_indvec_data[i2] + kweig->size[0] * iii) -
                      1] = a_data[i2];
                  }
                } else {
                  bb_binary_expand_op(kweig, b_indvec, iii, ktree, k);
                  kweig_data = kweig->data;
                }

                /* 'gkmPWM:806' kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1}; */
                if (indvec2->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->size
                    [0]) {
                  i2 = b_kweig->size[0];
                  b_kweig->size[0] = indvec2->size[0];
                  emxEnsureCapacity_real_T(b_kweig, i2);
                  a_data = b_kweig->data;
                  m = indvec2->size[0];
                  for (i2 = 0; i2 < m; i2++) {
                    a_data[i2] = kweig_data[((int)indvec2_data[i2] + kweig->
                      size[0] * (3 - iii)) - 1] + ktree_data[(int)(k - 1.0) - 1]
                      .f1->data[i2];
                  }

                  m = b_kweig->size[0];
                  for (i2 = 0; i2 < m; i2++) {
                    kweig_data[((int)indvec2_data[i2] + kweig->size[0] * (3 -
                      iii)) - 1] = a_data[i2];
                  }
                } else {
                  ab_binary_expand_op(kweig, indvec2, iii, ktree, k);
                  kweig_data = kweig->data;
                }

                /* 'gkmPWM:807' kweig(indvec2,iii) = kweig(indvec2,iii) + ktree2{k-1}; */
                if (indvec2->size[0] == ktree2_data[(int)(k - 1.0) - 1].f1->
                    size[0]) {
                  i2 = b_kweig->size[0];
                  b_kweig->size[0] = indvec2->size[0];
                  emxEnsureCapacity_real_T(b_kweig, i2);
                  a_data = b_kweig->data;
                  m = indvec2->size[0];
                  for (i2 = 0; i2 < m; i2++) {
                    a_data[i2] = kweig_data[((int)indvec2_data[i2] + kweig->
                      size[0] * iii) - 1] + ktree2_data[(int)(k - 1.0) - 1]
                      .f1->data[i2];
                  }

                  m = b_kweig->size[0];
                  for (i2 = 0; i2 < m; i2++) {
                    kweig_data[((int)indvec2_data[i2] + kweig->size[0] * iii) -
                      1] = a_data[i2];
                  }
                } else {
                  bb_binary_expand_op(kweig, indvec2, iii, ktree2, k);
                  kweig_data = kweig->data;
                }

                /* 'gkmPWM:808' kweig(indvec,5-iii) = kweig(indvec,5-iii) + ktree2{k-1}; */
                if (b_indvec->size[0] == ktree2_data[(int)(k - 1.0) - 1]
                    .f1->size[0]) {
                  i2 = b_kweig->size[0];
                  b_kweig->size[0] = b_indvec->size[0];
                  emxEnsureCapacity_real_T(b_kweig, i2);
                  a_data = b_kweig->data;
                  m = b_indvec->size[0];
                  for (i2 = 0; i2 < m; i2++) {
                    a_data[i2] = kweig_data[((int)b_indvec_data[i2] +
                      kweig->size[0] * (3 - iii)) - 1] + ktree2_data[(int)(k -
                      1.0) - 1].f1->data[i2];
                  }

                  m = b_kweig->size[0];
                  for (i2 = 0; i2 < m; i2++) {
                    kweig_data[((int)b_indvec_data[i2] + kweig->size[0] * (3 -
                      iii)) - 1] = a_data[i2];
                  }
                } else {
                  ab_binary_expand_op(kweig, b_indvec, iii, ktree2, k);
                  kweig_data = kweig->data;
                }
              }
            }
          }
        } else {
          /* 'gkmPWM:812' else */
          /* 'gkmPWM:813' for iii = 1:2 */
          c_tmp = pow(4.0, (double)(f + 1) - 1.0);
          b_c = b_n_tmp * (ind_data[ii] - 1.0);
          md2 = poscell_data[f].f1->size[0];
          for (iii = 0; iii < 2; iii++) {
            /* 'gkmPWM:814' indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(ii)-1); */
            n = c_tmp * (((double)iii + 1.0) - 1.0);
            i1 = b_indvec->size[0];
            b_indvec->size[0] = poscell_data[f].f1->size[0];
            emxEnsureCapacity_real_T(b_indvec, i1);
            b_indvec_data = b_indvec->data;
            for (i1 = 0; i1 < md2; i1++) {
              b_indvec_data[i1] = (poscell_data[f].f1->data[i1] + n) + b_c;
            }

            /* 'gkmPWM:815' indvec2 = indvec+(5-2*iii)*4^(f-1); */
            d = (5.0 - 2.0 * ((double)iii + 1.0)) * c_tmp;
            i1 = indvec2->size[0];
            indvec2->size[0] = b_indvec->size[0];
            emxEnsureCapacity_real_T(indvec2, i1);
            indvec2_data = indvec2->data;
            m = b_indvec->size[0];
            for (i1 = 0; i1 < m; i1++) {
              indvec2_data[i1] = b_indvec_data[i1] + d;
            }

            /* 'gkmPWM:816' kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1}; */
            if (b_indvec->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->size[0])
            {
              i1 = b_kweig->size[0];
              b_kweig->size[0] = b_indvec->size[0];
              emxEnsureCapacity_real_T(b_kweig, i1);
              a_data = b_kweig->data;
              m = b_indvec->size[0];
              for (i1 = 0; i1 < m; i1++) {
                a_data[i1] = kweig_data[((int)b_indvec_data[i1] + kweig->size[0]
                  * iii) - 1] + ktree_data[(int)(k - 1.0) - 1].f1->data[i1];
              }

              m = b_kweig->size[0];
              for (i1 = 0; i1 < m; i1++) {
                kweig_data[((int)b_indvec_data[i1] + kweig->size[0] * iii) - 1] =
                  a_data[i1];
              }
            } else {
              bb_binary_expand_op(kweig, b_indvec, iii, ktree, k);
              kweig_data = kweig->data;
            }

            /* 'gkmPWM:817' kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1}; */
            if (indvec2->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->size[0])
            {
              i1 = b_kweig->size[0];
              b_kweig->size[0] = indvec2->size[0];
              emxEnsureCapacity_real_T(b_kweig, i1);
              a_data = b_kweig->data;
              m = indvec2->size[0];
              for (i1 = 0; i1 < m; i1++) {
                a_data[i1] = kweig_data[((int)indvec2_data[i1] + kweig->size[0] *
                  (3 - iii)) - 1] + ktree_data[(int)(k - 1.0) - 1].f1->data[i1];
              }

              m = b_kweig->size[0];
              for (i1 = 0; i1 < m; i1++) {
                kweig_data[((int)indvec2_data[i1] + kweig->size[0] * (3 - iii))
                  - 1] = a_data[i1];
              }
            } else {
              ab_binary_expand_op(kweig, indvec2, iii, ktree, k);
              kweig_data = kweig->data;
            }

            /* 'gkmPWM:818' kweig(indvec2,iii) = kweig(indvec2,iii) + ktree2{k-1}; */
            if (indvec2->size[0] == ktree2_data[(int)(k - 1.0) - 1].f1->size[0])
            {
              i1 = b_kweig->size[0];
              b_kweig->size[0] = indvec2->size[0];
              emxEnsureCapacity_real_T(b_kweig, i1);
              a_data = b_kweig->data;
              m = indvec2->size[0];
              for (i1 = 0; i1 < m; i1++) {
                a_data[i1] = kweig_data[((int)indvec2_data[i1] + kweig->size[0] *
                  iii) - 1] + ktree2_data[(int)(k - 1.0) - 1].f1->data[i1];
              }

              m = b_kweig->size[0];
              for (i1 = 0; i1 < m; i1++) {
                kweig_data[((int)indvec2_data[i1] + kweig->size[0] * iii) - 1] =
                  a_data[i1];
              }
            } else {
              bb_binary_expand_op(kweig, indvec2, iii, ktree2, k);
              kweig_data = kweig->data;
            }

            /* 'gkmPWM:819' kweig(indvec,5-iii) = kweig(indvec,5-iii) + ktree2{k-1}; */
            if (b_indvec->size[0] == ktree2_data[(int)(k - 1.0) - 1].f1->size[0])
            {
              i1 = b_kweig->size[0];
              b_kweig->size[0] = b_indvec->size[0];
              emxEnsureCapacity_real_T(b_kweig, i1);
              a_data = b_kweig->data;
              m = b_indvec->size[0];
              for (i1 = 0; i1 < m; i1++) {
                a_data[i1] = kweig_data[((int)b_indvec_data[i1] + kweig->size[0]
                  * (3 - iii)) - 1] + ktree2_data[(int)(k - 1.0) - 1].f1->
                  data[i1];
              }

              m = b_kweig->size[0];
              for (i1 = 0; i1 < m; i1++) {
                kweig_data[((int)b_indvec_data[i1] + kweig->size[0] * (3 - iii))
                  - 1] = a_data[i1];
              }
            } else {
              ab_binary_expand_op(kweig, b_indvec, iii, ktree2, k);
              kweig_data = kweig->data;
            }
          }
        }
      }
    }
  }

  emxFree_real_T(&b_kweig);
  emxFree_boolean_T(&b_x);
  emxFree_real_T(&b_indvec);
  emxFree_cell_wrap_13(&p);
  emxFree_real_T(&indvec2);
  emxFree_real_T(&a);
  emxFree_real_T(&sPWM2);
  emxFree_real_T(&sPWM);
  emxFree_int32_T(&f1);
  emxFree_real_T(&indvec);
  emxFree_real_T(&indloc2);
  emxFree_uint32_T(&X);
  emxFree_cell_wrap_14(&ktree2);
  emxFree_cell_wrap_14(&ktree);

  /* 'gkmPWM:825' kweig(4^k*(max(max(x))-rcnum)+1:end,:) = kweig(4^k*(max(max(x))-rcnum)+1:end,:)/sqrt(2); */
  d = b_n_tmp * (n_tmp - rcnum) + 1.0;
  if (d > kweig->size[0]) {
    i = 0;
    i1 = 0;
    i2 = 1;
  } else {
    i = (int)d - 1;
    i1 = kweig->size[0];
    i2 = (int)d;
  }

  md2 = i1 - i;
  i1 = mat2->size[0] * mat2->size[1];
  mat2->size[0] = md2;
  mat2->size[1] = 4;
  emxEnsureCapacity_real_T(mat2, i1);
  mat2_data = mat2->data;
  for (i1 = 0; i1 < 4; i1++) {
    for (m = 0; m < md2; m++) {
      mat2_data[m + mat2->size[0] * i1] = kweig_data[(i + m) + kweig->size[0] *
        i1] / 1.4142135623730951;
    }
  }

  md2 = mat2->size[0];
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < md2; i1++) {
      kweig_data[((i2 + i1) + kweig->size[0] * i) - 1] = mat2_data[i1 +
        mat2->size[0] * i];
    }
  }

  emxFree_real_T(&mat2);
}

/*
 * function kweig = ls_kweigtree_norc(mat,negmat,poscell,c,s,ind,indloc,x,l,k,rcnum)
 */
static void ls_kweigtree_norc(const emxArray_real_T *mat, const double negmat[16],
  const emxArray_cell_wrap_14 *poscell, const emxArray_real_T *c, const
  emxArray_real_T *s, const emxArray_real_T *ind, const emxArray_real_T *indloc,
  const emxArray_real_T *x, double l, double k, emxArray_real_T *kweig)
{
  cell_wrap_13 *p_data;
  const cell_wrap_14 *poscell_data;
  cell_wrap_14 *ktree_data;
  emxArray_boolean_T *b_x;
  emxArray_cell_wrap_13 *p;
  emxArray_cell_wrap_14 *ktree;
  emxArray_int32_T *f1;
  emxArray_real_T *a;
  emxArray_real_T *b_indvec;
  emxArray_real_T *b_kweig;
  emxArray_real_T *indloc2;
  emxArray_real_T *indvec;
  emxArray_real_T *indvec2;
  emxArray_real_T *sPWM;
  emxArray_uint32_T *X;
  double b_p[16];
  double c_x[4];
  const double *c_data;
  const double *ind_data;
  const double *indloc_data;
  const double *mat_data;
  const double *s_data;
  const double *x_data;
  double b_c;
  double b_m;
  double b_n;
  double d;
  double n_tmp;
  double varargin_1;
  double *a_data;
  double *b_indvec_data;
  double *indloc2_data;
  double *indvec2_data;
  double *indvec_data;
  double *kweig_data;
  double *sPWM_data;
  int b_i;
  int b_k;
  int f;
  int i;
  int i1;
  int idx;
  int ii;
  int iii;
  int itilerow;
  int j;
  int loop_ub_tmp;
  int m;
  int n;
  int nz;
  int rx;
  unsigned int *X_data;
  int *f1_data;
  bool *b_x_data;
  x_data = x->data;
  indloc_data = indloc->data;
  ind_data = ind->data;
  s_data = s->data;
  c_data = c->data;
  poscell_data = poscell->data;
  mat_data = mat->data;
  emxInit_cell_wrap_13(&p);

  /* uses dynamic programming to find the matrix A (kweig in this code) such that A*p = expected gapped kmer vector, where p is the middle column of a 2*l-1 PWM. */
  /* I use a 1st order Markov model to model the flanking bases for a PWM */
  /* 'gkmPWM:832' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_13(p, i);
  p_data = p->data;

  /* 'gkmPWM:833' p = coder.nullcopy(p); */
  /* 'gkmPWM:834' p{1} = eye(4); */
  for (i = 0; i < 16; i++) {
    p_data[0].f1[i] = 0.0;
  }

  for (b_k = 0; b_k < 4; b_k++) {
    p_data[0].f1[b_k + (b_k << 2)] = 1.0;
  }

  /* 'gkmPWM:835' for i = 1:l-1 */
  i = (int)(l - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:836' p{i+1} = p{i}*negmat; */
    for (i1 = 0; i1 < 4; i1++) {
      for (b_k = 0; b_k < 4; b_k++) {
        d = 0.0;
        for (n = 0; n < 4; n++) {
          d += p_data[b_i].f1[i1 + (n << 2)] * negmat[n + (b_k << 2)];
        }

        b_p[i1 + (b_k << 2)] = d;
      }
    }

    for (i1 = 0; i1 < 16; i1++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i1] = b_p[i1];
    }
  }

  emxInit_real_T(&indvec, 2);

  /* 'gkmPWM:838' n = 4^k*max(max(x)); */
  m = x->size[0];
  n = x->size[1];
  i = indvec->size[0] * indvec->size[1];
  indvec->size[0] = 1;
  indvec->size[1] = x->size[1];
  emxEnsureCapacity_real_T(indvec, i);
  indvec_data = indvec->data;
  if (x->size[1] >= 1) {
    for (j = 0; j < n; j++) {
      indvec_data[j] = x_data[x->size[0] * j];
      for (b_i = 2; b_i <= m; b_i++) {
        d = x_data[(b_i + x->size[0] * j) - 1];
        if (indvec_data[j] < d) {
          indvec_data[j] = d;
        }
      }
    }
  }

  n = indvec->size[1];
  if (indvec->size[1] <= 2) {
    if (indvec->size[1] == 1) {
      b_n = indvec_data[0];
    } else if (indvec_data[0] < indvec_data[indvec->size[1] - 1]) {
      b_n = indvec_data[indvec->size[1] - 1];
    } else {
      b_n = indvec_data[0];
    }
  } else {
    b_n = indvec_data[0];
    for (b_k = 2; b_k <= n; b_k++) {
      d = indvec_data[b_k - 1];
      if (b_n < d) {
        b_n = d;
      }
    }
  }

  n_tmp = pow(4.0, k);
  b_n *= n_tmp;

  /* number of possible k-mers */
  /* 'gkmPWM:839' mat2 = rot90(mat,2); */
  /* 'gkmPWM:840' kweig = zeros(n, 4); */
  i = kweig->size[0] * kweig->size[1];
  kweig->size[0] = (int)b_n;
  kweig->size[1] = 4;
  emxEnsureCapacity_real_T(kweig, i);
  kweig_data = kweig->data;
  nz = (int)b_n << 2;
  for (i = 0; i < nz; i++) {
    kweig_data[i] = 0.0;
  }

  emxInit_cell_wrap_14(&ktree);

  /* 'gkmPWM:841' ktree = cell(k-1,1); */
  n = (int)(k - 1.0);
  i = ktree->size[0];
  ktree->size[0] = (int)(k - 1.0);
  emxEnsureCapacity_cell_wrap_14(ktree, i);
  ktree_data = ktree->data;
  for (i = 0; i < n; i++) {
    ktree_data[i].f1->size[0] = 0;
  }

  emxInit_uint32_T(&X);

  /* 'gkmPWM:842' ktree = coder.nullcopy(ktree); */
  /* 'gkmPWM:843' [rx,cx] = size(x); */
  rx = x->size[0];

  /* 'gkmPWM:844' m=rx; */
  b_m = x->size[0];

  /* 'gkmPWM:845' M = l-1; */
  /* 'gkmPWM:846' X = cx*ones(M+1,1); */
  loop_ub_tmp = (int)((l - 1.0) + 1.0);
  i = X->size[0];
  X->size[0] = (int)((l - 1.0) + 1.0);
  emxEnsureCapacity_uint32_T(X, i);
  X_data = X->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    X_data[i] = (unsigned int)x->size[1];
  }

  /* 'gkmPWM:847' for i = 1:cx */
  i = x->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:848' X(i) = i; */
    X_data[b_i] = (unsigned int)(b_i + 1);
  }

  emxInit_real_T(&indloc2, 1);

  /* 'gkmPWM:850' indloc2 = flipud(indloc); */
  i = indloc2->size[0];
  indloc2->size[0] = indloc->size[0];
  emxEnsureCapacity_real_T(indloc2, i);
  indloc2_data = indloc2->data;
  nz = indloc->size[0];
  for (i = 0; i < nz; i++) {
    indloc2_data[i] = indloc_data[i];
  }

  m = indloc->size[0] - 1;
  n = indloc->size[0] >> 1;
  for (b_i = 0; b_i < n; b_i++) {
    b_n = indloc2_data[b_i];
    nz = m - b_i;
    indloc2_data[b_i] = indloc2_data[nz];
    indloc2_data[nz] = b_n;
  }

  /* 'gkmPWM:851' for i = 2:5 */
  for (b_i = 0; b_i < 4; b_i++) {
    /* 'gkmPWM:852' ktree{i} = zeros(4^i,1); */
    n = (int)pow(4.0, (double)b_i + 2.0);
    i = ktree_data[b_i + 1].f1->size[0];
    ktree_data[b_i + 1].f1->size[0] = n;
    emxEnsureCapacity_real_T(ktree_data[b_i + 1].f1, i);
    for (i = 0; i < n; i++) {
      ktree_data[b_i + 1].f1->data[i] = 0.0;
    }
  }

  /* 'gkmPWM:854' for i = 0:M */
  emxInit_int32_T(&f1, 2);
  emxInit_real_T(&sPWM, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&indvec2, 1);
  emxInit_real_T(&b_indvec, 1);
  emxInit_boolean_T(&b_x, 2);
  emxInit_real_T(&b_kweig, 1);
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    /* 'gkmPWM:855' if i > M-cx+1 */
    if (b_i > ((l - 1.0) - (double)x->size[1]) + 1.0) {
      /* 'gkmPWM:856' m = numel(c)/k; */
      b_m = (double)(c->size[0] * c->size[1]) / k;
    }

    /* the following loops is basically dynamic programming for tensor multiplication.  there are multiple cases to consider, hence the if statements. */
    /* 'gkmPWM:859' for ii = 1:m */
    i = (int)b_m;
    for (ii = 0; ii < i; ii++) {
      /* 'gkmPWM:860' if sum((c(ii,:)+i)==l) > 0 && ~(i == M-1 && ii > rx && ii ~= m) */
      nz = c->size[1];
      i1 = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 1;
      b_x->size[1] = c->size[1];
      emxEnsureCapacity_boolean_T(b_x, i1);
      b_x_data = b_x->data;
      for (i1 = 0; i1 < nz; i1++) {
        b_x_data[i1] = (c_data[ii + c->size[0] * i1] + (double)b_i == l);
      }

      n = b_x->size[1];
      if (b_x->size[1] == 0) {
        nz = 0;
      } else {
        nz = b_x_data[0];
        for (b_k = 2; b_k <= n; b_k++) {
          nz += b_x_data[b_k - 1];
        }
      }

      if ((nz > 0) && ((b_i != (l - 1.0) - 1.0) || (ii + 1U <= (unsigned int)rx)
                       || ((double)ii + 1.0 == b_m))) {
        /* 'gkmPWM:861' indvec = c(ii,:)+i; */
        nz = c->size[1];
        i1 = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i1);
        indvec_data = indvec->data;
        for (i1 = 0; i1 < nz; i1++) {
          indvec_data[i1] = c_data[ii + c->size[0] * i1] + (double)b_i;
        }

        /*  f = find(indvec == l); */
        /* 'gkmPWM:863' coder.varsize('f', [1 1]); */
        /* 'gkmPWM:864' f1 = find(indvec == l); */
        i1 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = indvec->size[1];
        emxEnsureCapacity_boolean_T(b_x, i1);
        b_x_data = b_x->data;
        nz = indvec->size[1];
        for (i1 = 0; i1 < nz; i1++) {
          b_x_data[i1] = (indvec_data[i1] == l);
        }

        eml_find(b_x, f1);
        f1_data = f1->data;

        /* 'gkmPWM:865' if length(f1) == 0 */
        if (f1->size[1] == 0) {
          /* 'gkmPWM:866' f = zeros(1,1); */
          f = -1;

          /* 'gkmPWM:867' fprintf("Resizing f\n"); */
          printf("Resizing f\n");
          fflush(stdout);
        } else {
          /* 'gkmPWM:868' else */
          /* 'gkmPWM:869' f = f1(1); */
          f = f1_data[0] - 1;
        }

        /* 'gkmPWM:871' indvec(f) = []; */
        idx = f + 1;
        n = indvec->size[1];
        nz = indvec->size[1] - 1;
        for (b_k = idx; b_k <= nz; b_k++) {
          indvec_data[b_k - 1] = indvec_data[b_k];
        }

        i1 = indvec->size[0] * indvec->size[1];
        if (1 > nz) {
          indvec->size[1] = 0;
        } else {
          indvec->size[1] = n - 1;
        }

        emxEnsureCapacity_real_T(indvec, i1);
        indvec_data = indvec->data;

        /* 'gkmPWM:872' loc = indloc(indvec); */
        /* 'gkmPWM:873' loc2 = indloc2(indvec); */
        /* 'gkmPWM:874' sPWM = mat(indvec,:).'; */
        i1 = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = 4;
        sPWM->size[1] = indvec->size[1];
        emxEnsureCapacity_real_T(sPWM, i1);
        sPWM_data = sPWM->data;
        nz = indvec->size[1];
        for (i1 = 0; i1 < nz; i1++) {
          for (b_k = 0; b_k < 4; b_k++) {
            sPWM_data[b_k + 4 * i1] = mat_data[((int)indvec_data[i1] + mat->
              size[0] * b_k) - 1];
          }
        }

        /* 'gkmPWM:875' ktree{1} = sPWM(:,1); */
        i1 = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree_data[0].f1, i1);
        for (i1 = 0; i1 < 4; i1++) {
          ktree_data[0].f1->data[i1] = sPWM_data[i1];
        }

        /* 'gkmPWM:876' for iii = s(ii,i+1):k-1 */
        d = s_data[ii + s->size[0] * b_i];
        i1 = (int)((k - 1.0) + (1.0 - d));
        for (iii = 0; iii < i1; iii++) {
          b_n = d + (double)iii;

          /* 'gkmPWM:877' if loc(iii)==0 */
          b_k = (int)indvec_data[(int)b_n - 1] - 1;
          if (indloc_data[b_k] == 0.0) {
            /* 'gkmPWM:878' if loc(iii-1)==1 && indvec(iii-1) < l */
            if ((indloc_data[(int)indvec_data[(int)(b_n - 1.0) - 1] - 1] == 1.0)
                && ((unsigned int)indvec_data[(int)(b_n - 1.0) - 1] < l)) {
              /* 'gkmPWM:879' matt = sPWM(:,iii).'*p{indvec(iii)-l}; */
              for (b_k = 0; b_k < 16; b_k++) {
                b_p[b_k] = p_data[(int)((double)(unsigned int)indvec_data[(int)
                  b_n - 1] - l) - 1].f1[b_k];
              }

              /* 'gkmPWM:880' ktree{iii} = repmat(ktree{iii-1},4,1).*repelem(matt', 4^(iii-1)); */
              b_k = b_indvec->size[0];
              b_indvec->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0] <<
                2;
              emxEnsureCapacity_real_T(b_indvec, b_k);
              b_indvec_data = b_indvec->data;
              nz = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              varargin_1 = pow(4.0, b_n - 1.0);
              b_k = indvec2->size[0];
              indvec2->size[0] = (int)varargin_1 << 2;
              emxEnsureCapacity_real_T(indvec2, b_k);
              indvec2_data = indvec2->data;
              idx = -1;
              m = (int)varargin_1;
              for (itilerow = 0; itilerow < 4; itilerow++) {
                n = itilerow * nz;
                for (b_k = 0; b_k < nz; b_k++) {
                  b_indvec_data[n + b_k] = ktree_data[(int)(b_n - 1.0) - 1]
                    .f1->data[b_k];
                }

                b_c = 0.0;
                for (b_k = 0; b_k < 4; b_k++) {
                  b_c += sPWM_data[b_k + 4 * ((int)b_n - 1)] * b_p[b_k +
                    (itilerow << 2)];
                }

                c_x[itilerow] = b_c;
                for (j = 0; j < m; j++) {
                  indvec2_data[(idx + j) + 1] = c_x[itilerow];
                }

                idx += (int)varargin_1;
              }

              if (b_indvec->size[0] == indvec2->size[0]) {
                b_k = ktree_data[(int)b_n - 1].f1->size[0];
                ktree_data[(int)b_n - 1].f1->size[0] = b_indvec->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
                nz = b_indvec->size[0];
                for (b_k = 0; b_k < nz; b_k++) {
                  ktree_data[(int)b_n - 1].f1->data[b_k] = b_indvec_data[b_k] *
                    indvec2_data[b_k];
                }
              } else {
                r_binary_expand_op(ktree, b_n, b_indvec, indvec2);
                ktree_data = ktree->data;
              }
            } else {
              /* 'gkmPWM:881' else */
              /* 'gkmPWM:882' matt = p{indvec(iii)-indvec(iii-1)+1}; */
              /* 'gkmPWM:883' ktree{iii} = repmat(ktree{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2)); */
              b_k = b_indvec->size[0];
              b_indvec->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0] <<
                2;
              emxEnsureCapacity_real_T(b_indvec, b_k);
              b_indvec_data = b_indvec->data;
              nz = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              for (itilerow = 0; itilerow < 4; itilerow++) {
                n = itilerow * nz;
                for (b_k = 0; b_k < nz; b_k++) {
                  b_indvec_data[n + b_k] = ktree_data[(int)(b_n - 1.0) - 1]
                    .f1->data[b_k];
                }
              }

              varargin_1 = pow(4.0, b_n - 2.0);
              b_k = indvec2->size[0];
              indvec2->size[0] = (int)varargin_1 << 4;
              emxEnsureCapacity_real_T(indvec2, b_k);
              indvec2_data = indvec2->data;
              idx = -1;
              m = (int)varargin_1;
              for (b_k = 0; b_k < 16; b_k++) {
                for (j = 0; j < m; j++) {
                  indvec2_data[(idx + j) + 1] = p_data[(int)(unsigned int)
                    indvec_data[(int)b_n - 1] - (int)(unsigned int)indvec_data
                    [(int)(b_n - 1.0) - 1]].f1[b_k];
                }

                idx += (int)varargin_1;
              }

              if (b_indvec->size[0] == indvec2->size[0]) {
                b_k = ktree_data[(int)b_n - 1].f1->size[0];
                ktree_data[(int)b_n - 1].f1->size[0] = b_indvec->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
                nz = b_indvec->size[0];
                for (b_k = 0; b_k < nz; b_k++) {
                  ktree_data[(int)b_n - 1].f1->data[b_k] = b_indvec_data[b_k] *
                    indvec2_data[b_k];
                }
              } else {
                r_binary_expand_op(ktree, b_n, b_indvec, indvec2);
                ktree_data = ktree->data;
              }
            }
          } else if (indloc2_data[b_k] == 0.0) {
            /* 'gkmPWM:885' elseif loc2(iii)==0 */
            /* 'gkmPWM:886' if loc2(iii-1)==1 && indvec(iii-1) < l */
            if ((indloc2_data[(int)indvec_data[(int)(b_n - 1.0) - 1] - 1] == 1.0)
                && ((unsigned int)indvec_data[(int)(b_n - 1.0) - 1] < l)) {
              /* 'gkmPWM:887' a = ktree{iii-1}.*sPWM(:,iii).'; */
              b_k = a->size[0] * a->size[1];
              a->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, b_k);
              a_data = a->data;
              nz = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              for (b_k = 0; b_k < 4; b_k++) {
                for (n = 0; n < nz; n++) {
                  a_data[n + a->size[0] * b_k] = ktree_data[(int)(b_n - 1.0) - 1]
                    .f1->data[n] * sPWM_data[b_k + 4 * ((int)b_n - 1)];
                }
              }

              /* 'gkmPWM:888' ktree{iii} = a(:); */
              b_k = ktree_data[(int)b_n - 1].f1->size[0];
              ktree_data[(int)b_n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
              nz = a->size[0] << 2;
              for (b_k = 0; b_k < nz; b_k++) {
                ktree_data[(int)b_n - 1].f1->data[b_k] = a_data[b_k];
              }
            } else {
              /* 'gkmPWM:889' else */
              /* 'gkmPWM:890' a = ktree{iii-1}.*sPWM(:,iii).'; */
              b_k = a->size[0] * a->size[1];
              a->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, b_k);
              a_data = a->data;
              nz = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              for (b_k = 0; b_k < 4; b_k++) {
                for (n = 0; n < nz; n++) {
                  a_data[n + a->size[0] * b_k] = ktree_data[(int)(b_n - 1.0) - 1]
                    .f1->data[n] * sPWM_data[b_k + 4 * ((int)b_n - 1)];
                }
              }

              /* 'gkmPWM:891' ktree{iii} = a(:); */
              b_k = ktree_data[(int)b_n - 1].f1->size[0];
              ktree_data[(int)b_n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
              nz = a->size[0] << 2;
              for (b_k = 0; b_k < nz; b_k++) {
                ktree_data[(int)b_n - 1].f1->data[b_k] = a_data[b_k];
              }
            }
          } else {
            /* 'gkmPWM:893' else */
            /* 'gkmPWM:894' a = ktree{iii-1}.*sPWM(:,iii).'; */
            b_k = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, b_k);
            a_data = a->data;
            nz = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
            for (b_k = 0; b_k < 4; b_k++) {
              for (n = 0; n < nz; n++) {
                a_data[n + a->size[0] * b_k] = ktree_data[(int)(b_n - 1.0) - 1].
                  f1->data[n] * sPWM_data[b_k + 4 * ((int)b_n - 1)];
              }
            }

            /* 'gkmPWM:895' ktree{iii} = a(:); */
            b_k = ktree_data[(int)b_n - 1].f1->size[0];
            ktree_data[(int)b_n - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
            nz = a->size[0] << 2;
            for (b_k = 0; b_k < nz; b_k++) {
              ktree_data[(int)b_n - 1].f1->data[b_k] = a_data[b_k];
            }
          }
        }

        /* the weird indexing that I did early in the code comes to fruition.  It is critical to do so to make this computation as fast as possible. */
        /* 'gkmPWM:899' if ii <= rx */
        if (ii + 1U <= (unsigned int)rx) {
          /* 'gkmPWM:900' for j = 1:X(i+1) */
          i1 = (int)X_data[b_i];
          for (j = 0; j < i1; j++) {
            /* 'gkmPWM:901' if x(ii,j) ~= 0 */
            d = x_data[ii + x->size[0] * j];
            if (d != 0.0) {
              /* 'gkmPWM:902' for iii = 1:2 */
              varargin_1 = pow(4.0, (double)(f + 1) - 1.0);
              b_c = n_tmp * (ind_data[(int)d - 1] - 1.0);
              nz = poscell_data[f].f1->size[0];
              for (iii = 0; iii < 2; iii++) {
                /* 'gkmPWM:903' indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(x(ii,j))-1); */
                b_n = varargin_1 * (((double)iii + 1.0) - 1.0);
                b_k = b_indvec->size[0];
                b_indvec->size[0] = poscell_data[f].f1->size[0];
                emxEnsureCapacity_real_T(b_indvec, b_k);
                b_indvec_data = b_indvec->data;
                for (b_k = 0; b_k < nz; b_k++) {
                  b_indvec_data[b_k] = (poscell_data[f].f1->data[b_k] + b_n) +
                    b_c;
                }

                /* 'gkmPWM:904' indvec2 = indvec+(5-2*iii)*4^(f-1); */
                d = (5.0 - 2.0 * ((double)iii + 1.0)) * varargin_1;
                b_k = indvec2->size[0];
                indvec2->size[0] = b_indvec->size[0];
                emxEnsureCapacity_real_T(indvec2, b_k);
                indvec2_data = indvec2->data;
                n = b_indvec->size[0];
                for (b_k = 0; b_k < n; b_k++) {
                  indvec2_data[b_k] = b_indvec_data[b_k] + d;
                }

                /* 'gkmPWM:905' kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1}; */
                if (b_indvec->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->
                    size[0]) {
                  b_k = b_kweig->size[0];
                  b_kweig->size[0] = b_indvec->size[0];
                  emxEnsureCapacity_real_T(b_kweig, b_k);
                  a_data = b_kweig->data;
                  n = b_indvec->size[0];
                  for (b_k = 0; b_k < n; b_k++) {
                    a_data[b_k] = kweig_data[((int)b_indvec_data[b_k] +
                      kweig->size[0] * iii) - 1] + ktree_data[(int)(k - 1.0) - 1]
                      .f1->data[b_k];
                  }

                  n = b_kweig->size[0];
                  for (b_k = 0; b_k < n; b_k++) {
                    kweig_data[((int)b_indvec_data[b_k] + kweig->size[0] * iii)
                      - 1] = a_data[b_k];
                  }
                } else {
                  bb_binary_expand_op(kweig, b_indvec, iii, ktree, k);
                  kweig_data = kweig->data;
                }

                /* 'gkmPWM:906' kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1}; */
                if (indvec2->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->size
                    [0]) {
                  b_k = b_kweig->size[0];
                  b_kweig->size[0] = indvec2->size[0];
                  emxEnsureCapacity_real_T(b_kweig, b_k);
                  a_data = b_kweig->data;
                  n = indvec2->size[0];
                  for (b_k = 0; b_k < n; b_k++) {
                    a_data[b_k] = kweig_data[((int)indvec2_data[b_k] +
                      kweig->size[0] * (3 - iii)) - 1] + ktree_data[(int)(k -
                      1.0) - 1].f1->data[b_k];
                  }

                  n = b_kweig->size[0];
                  for (b_k = 0; b_k < n; b_k++) {
                    kweig_data[((int)indvec2_data[b_k] + kweig->size[0] * (3 -
                      iii)) - 1] = a_data[b_k];
                  }
                } else {
                  ab_binary_expand_op(kweig, indvec2, iii, ktree, k);
                  kweig_data = kweig->data;
                }
              }
            }
          }
        } else {
          /* 'gkmPWM:910' else */
          /* 'gkmPWM:911' for iii = 1:2 */
          varargin_1 = pow(4.0, (double)(f + 1) - 1.0);
          b_c = n_tmp * (ind_data[ii] - 1.0);
          nz = poscell_data[f].f1->size[0];
          for (iii = 0; iii < 2; iii++) {
            /* 'gkmPWM:912' indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(ii)-1); */
            b_n = varargin_1 * (((double)iii + 1.0) - 1.0);
            i1 = b_indvec->size[0];
            b_indvec->size[0] = poscell_data[f].f1->size[0];
            emxEnsureCapacity_real_T(b_indvec, i1);
            b_indvec_data = b_indvec->data;
            for (i1 = 0; i1 < nz; i1++) {
              b_indvec_data[i1] = (poscell_data[f].f1->data[i1] + b_n) + b_c;
            }

            /* 'gkmPWM:913' indvec2 = indvec+(5-2*iii)*4^(f-1); */
            d = (5.0 - 2.0 * ((double)iii + 1.0)) * varargin_1;
            i1 = indvec2->size[0];
            indvec2->size[0] = b_indvec->size[0];
            emxEnsureCapacity_real_T(indvec2, i1);
            indvec2_data = indvec2->data;
            n = b_indvec->size[0];
            for (i1 = 0; i1 < n; i1++) {
              indvec2_data[i1] = b_indvec_data[i1] + d;
            }

            /* 'gkmPWM:914' kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1}; */
            if (b_indvec->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->size[0])
            {
              i1 = b_kweig->size[0];
              b_kweig->size[0] = b_indvec->size[0];
              emxEnsureCapacity_real_T(b_kweig, i1);
              a_data = b_kweig->data;
              n = b_indvec->size[0];
              for (i1 = 0; i1 < n; i1++) {
                a_data[i1] = kweig_data[((int)b_indvec_data[i1] + kweig->size[0]
                  * iii) - 1] + ktree_data[(int)(k - 1.0) - 1].f1->data[i1];
              }

              n = b_kweig->size[0];
              for (i1 = 0; i1 < n; i1++) {
                kweig_data[((int)b_indvec_data[i1] + kweig->size[0] * iii) - 1] =
                  a_data[i1];
              }
            } else {
              bb_binary_expand_op(kweig, b_indvec, iii, ktree, k);
              kweig_data = kweig->data;
            }

            /* 'gkmPWM:915' kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1}; */
            if (indvec2->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->size[0])
            {
              i1 = b_kweig->size[0];
              b_kweig->size[0] = indvec2->size[0];
              emxEnsureCapacity_real_T(b_kweig, i1);
              a_data = b_kweig->data;
              n = indvec2->size[0];
              for (i1 = 0; i1 < n; i1++) {
                a_data[i1] = kweig_data[((int)indvec2_data[i1] + kweig->size[0] *
                  (3 - iii)) - 1] + ktree_data[(int)(k - 1.0) - 1].f1->data[i1];
              }

              n = b_kweig->size[0];
              for (i1 = 0; i1 < n; i1++) {
                kweig_data[((int)indvec2_data[i1] + kweig->size[0] * (3 - iii))
                  - 1] = a_data[i1];
              }
            } else {
              ab_binary_expand_op(kweig, indvec2, iii, ktree, k);
              kweig_data = kweig->data;
            }
          }
        }
      }
    }
  }

  emxFree_real_T(&b_kweig);
  emxFree_boolean_T(&b_x);
  emxFree_real_T(&b_indvec);
  emxFree_cell_wrap_13(&p);
  emxFree_real_T(&indvec2);
  emxFree_real_T(&a);
  emxFree_real_T(&sPWM);
  emxFree_int32_T(&f1);
  emxFree_real_T(&indvec);
  emxFree_real_T(&indloc2);
  emxFree_uint32_T(&X);
  emxFree_cell_wrap_14(&ktree);
}

static void m_binary_expand_op(emxArray_real_T *kmat, int i, const
  emxArray_real_T *y, const emxArray_real_T *negvec, double b)
{
  const double *negvec_data;
  const double *y_data;
  double *kmat_data;
  int b_i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  negvec_data = negvec->data;
  y_data = y->data;
  kmat_data = kmat->data;
  stride_0_0 = (y->size[0] != 1);
  stride_1_0 = (negvec->size[0] != 1);
  if (negvec->size[0] == 1) {
    loop_ub = y->size[0];
  } else {
    loop_ub = negvec->size[0];
  }

  for (b_i = 0; b_i < loop_ub; b_i++) {
    kmat_data[b_i + kmat->size[0] * i] = y_data[b_i * stride_0_0] -
      negvec_data[b_i * stride_1_0] * b;
  }
}

static void minus(emxArray_real_T *res, const emxArray_real_T *kweig)
{
  emxArray_real_T *b_kweig;
  const double *kweig_data;
  double *b_kweig_data;
  double *res_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  kweig_data = kweig->data;
  res_data = res->data;
  emxInit_real_T(&b_kweig, 1);
  i = b_kweig->size[0];
  if (res->size[0] == 1) {
    b_kweig->size[0] = kweig->size[0];
  } else {
    b_kweig->size[0] = res->size[0];
  }

  emxEnsureCapacity_real_T(b_kweig, i);
  b_kweig_data = b_kweig->data;
  stride_0_0 = (kweig->size[0] != 1);
  stride_1_0 = (res->size[0] != 1);
  if (res->size[0] == 1) {
    loop_ub = kweig->size[0];
  } else {
    loop_ub = res->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_kweig_data[i] = kweig_data[i * stride_0_0] - res_data[i * stride_1_0];
  }

  i = res->size[0];
  res->size[0] = b_kweig->size[0];
  emxEnsureCapacity_real_T(res, i);
  res_data = res->data;
  loop_ub = b_kweig->size[0];
  for (i = 0; i < loop_ub; i++) {
    res_data[i] = b_kweig_data[i];
  }

  emxFree_real_T(&b_kweig);
}

static void n_binary_expand_op(emxArray_real_T *res, const emxArray_real_T *f,
  const emxArray_real_T *C, const emxArray_real_T *ord, int ii)
{
  emxArray_real_T *b_res;
  const double *C_data;
  const double *f_data;
  const double *ord_data;
  double b_C;
  double *b_res_data;
  double *res_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  ord_data = ord->data;
  C_data = C->data;
  f_data = f->data;
  res_data = res->data;
  emxInit_real_T(&b_res, 1);
  b_C = C_data[(int)ord_data[ii] - 1];
  i = b_res->size[0];
  if (f->size[0] == 1) {
    b_res->size[0] = res->size[0];
  } else {
    b_res->size[0] = f->size[0];
  }

  emxEnsureCapacity_real_T(b_res, i);
  b_res_data = b_res->data;
  stride_0_0 = (res->size[0] != 1);
  stride_1_0 = (f->size[0] != 1);
  if (f->size[0] == 1) {
    loop_ub = res->size[0];
  } else {
    loop_ub = f->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_res_data[i] = res_data[i * stride_0_0] - f_data[i * stride_1_0] * b_C;
  }

  i = res->size[0];
  res->size[0] = b_res->size[0];
  emxEnsureCapacity_real_T(res, i);
  res_data = res->data;
  loop_ub = b_res->size[0];
  for (i = 0; i < loop_ub; i++) {
    res_data[i] = b_res_data[i];
  }

  emxFree_real_T(&b_res);
}

static void o_binary_expand_op(emxArray_real_T *kmat, const emxArray_real_T *ord,
  int ii, const emxArray_real_T *f)
{
  emxArray_real_T *b_kmat;
  const double *f_data;
  const double *ord_data;
  double *b_kmat_data;
  double *kmat_data;
  int c_kmat;
  int i;
  int ord_tmp;
  int stride_0_0;
  int stride_1_0;
  f_data = f->data;
  ord_data = ord->data;
  kmat_data = kmat->data;
  emxInit_real_T(&b_kmat, 1);
  c_kmat = kmat->size[0] - 1;
  ord_tmp = (int)ord_data[ii];
  i = b_kmat->size[0];
  if (f->size[0] == 1) {
    b_kmat->size[0] = c_kmat + 1;
  } else {
    b_kmat->size[0] = f->size[0];
  }

  emxEnsureCapacity_real_T(b_kmat, i);
  b_kmat_data = b_kmat->data;
  stride_0_0 = (c_kmat + 1 != 1);
  stride_1_0 = (f->size[0] != 1);
  if (f->size[0] == 1) {
    c_kmat++;
  } else {
    c_kmat = f->size[0];
  }

  for (i = 0; i < c_kmat; i++) {
    b_kmat_data[i] = kmat_data[i * stride_0_0 + kmat->size[0] * (ord_tmp - 1)] +
      f_data[i * stride_1_0];
  }

  c_kmat = b_kmat->size[0];
  for (i = 0; i < c_kmat; i++) {
    kmat_data[i + kmat->size[0] * (ord_tmp - 1)] = b_kmat_data[i];
  }

  emxFree_real_T(&b_kmat);
}

static void p_binary_expand_op(emxArray_real_T *kmat, const emxArray_real_T *f,
  int jj, const emxArray_real_T *y, const emxArray_real_T *negvec, double b)
{
  const double *f_data;
  const double *negvec_data;
  const double *y_data;
  double *kmat_data;
  int b_f;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  negvec_data = negvec->data;
  y_data = y->data;
  f_data = f->data;
  kmat_data = kmat->data;
  b_f = (int)f_data[jj];
  stride_0_0 = (y->size[0] != 1);
  stride_1_0 = (negvec->size[0] != 1);
  if (negvec->size[0] == 1) {
    loop_ub = y->size[0];
  } else {
    loop_ub = negvec->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    kmat_data[i + kmat->size[0] * (b_f - 1)] = y_data[i * stride_0_0] -
      negvec_data[i * stride_1_0] * b;
  }
}

/*
 * function [p, mat,pwms, c] = seed_kmers(fn, num, pn, ik)
 */
static void seed_kmers(const emxArray_char_T *fn, double num,
  emxArray_cell_wrap_0 *p, emxArray_cell_wrap_1 *mat, emxArray_cell_wrap_2 *pwms,
  double *c)
{
  cell_wrap_0 *p_data;
  cell_wrap_0 *s_data;
  cell_wrap_0 *sequences_data;
  cell_wrap_1 *b_mat_data;
  cell_wrap_1 *mat_data;
  cell_wrap_2 *b_pwms_data;
  cell_wrap_2 *pwms_data;
  emxArray_boolean_T *b_x;
  emxArray_cell_wrap_0 *b_p;
  emxArray_cell_wrap_0 *s;
  emxArray_cell_wrap_0 *sequences;
  emxArray_cell_wrap_1 *b_mat;
  emxArray_cell_wrap_2 *b_pwms;
  emxArray_char_T *b_fileid;
  emxArray_char_T *cur_alpha;
  emxArray_char_T *cur_line;
  emxArray_char_T *cur_seq;
  emxArray_int32_T *D;
  emxArray_int8_T *DD;
  emxArray_real_T *M;
  emxArray_real_T *alpha;
  emxArray_real_T *rs;
  emxArray_real_T *ss;
  creal_T dc;
  double h_y[9];
  double m[2];
  double curr_pos;
  double d;
  double idx;
  double *M_data;
  double *alpha_data;
  double *rs_data;
  double *ss_data;
  int CC[9];
  int x_size[2];
  int b_i;
  int b_loop_ub;
  int b_y;
  int c_loop_ub;
  int c_y;
  int cur_idx;
  int d_loop_ub;
  int d_y;
  int e_loop_ub;
  int e_y;
  int exitg1;
  int f_y;
  int g_y;
  int i;
  int i1;
  int i10;
  int i11;
  int i12;
  int i13;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int i8;
  int i9;
  int iindx;
  int j;
  int k;
  int l;
  int loop_ub;
  int loop_ub_tmp;
  int unnamed_idx_0_tmp_tmp;
  int varargin_2;
  int vlen;
  int x;
  int *D_data;
  signed char B[9];
  signed char BB[9];
  signed char fileid;
  signed char *DD_data;
  bool exitg2;
  bool y;
  bool *x_data;

  /* 'gkmPWM:172' fid = fopen(fn, 'r'); */
  fileid = cfopen(fn, "rb");

  /* 'gkmPWM:173' if fid == -1 */
  if (fileid == -1) {
    /* 'gkmPWM:174' fprintf("ERROR: Weight file cannot be opened.\n") */
    printf("ERROR: Weight file cannot be opened.\n");
    fflush(stdout);
    exit(1);
  }

  /*  a = textscan(fid, '%s\t%f\n'); */
  /* 'gkmPWM:177' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'gkmPWM:178' idx=0; */
  idx = 0.0;

  /* 'gkmPWM:179' while ~feof(fid) */
  emxInit_char_T(&b_fileid, 2);
  do {
    exitg1 = 0;
    d = b_feof(fileid);
    if (d == 0.0) {
      /* 'gkmPWM:180' idx=idx+1; */
      idx++;

      /* 'gkmPWM:181' fgetl(fid); */
      b_fgets(fileid, b_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_fileid);
  emxInit_cell_wrap_0(&sequences, 1);

  /* 'gkmPWM:183' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'gkmPWM:184' sequences = cell(idx, 1); */
  unnamed_idx_0_tmp_tmp = (int)idx;
  i = sequences->size[0];
  sequences->size[0] = (int)idx;
  emxEnsureCapacity_cell_wrap_0(sequences, i);
  sequences_data = sequences->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    sequences_data[i].f1->size[0] = 1;
    sequences_data[i].f1->size[1] = 0;
  }

  emxInit_real_T(&alpha, 1);

  /* 'gkmPWM:185' sequences = coder.nullcopy(sequences); */
  /* 'gkmPWM:186' alpha = zeros(idx, 1); */
  i = alpha->size[0];
  alpha->size[0] = (int)idx;
  emxEnsureCapacity_real_T(alpha, i);
  alpha_data = alpha->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    alpha_data[i] = 0.0;
  }

  /* 'gkmPWM:187' for cur_idx=1:idx */
  cur_idx = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&cur_seq, 2);
  emxInit_char_T(&cur_alpha, 2);
  exitg2 = false;
  while ((!exitg2) && (cur_idx <= (int)idx - 1)) {
    /* 'gkmPWM:188' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'gkmPWM:189' if cur_line == -1 */
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
      /* 'gkmPWM:192' [cur_seq, cur_alpha] = strtok(cur_line, char(9)); */
      b_strtok(cur_line, cur_seq, cur_alpha);

      /* 'gkmPWM:193' alpha(cur_idx,1) = real(str2double(cur_alpha)); */
      dc = str2double(cur_alpha);
      alpha_data[cur_idx] = dc.re;

      /* 'gkmPWM:194' sequences{cur_idx} = (strip(cur_seq)); */
      strip(cur_seq, sequences_data[cur_idx].f1);
      cur_idx++;
    }
  }

  emxFree_char_T(&cur_alpha);
  emxFree_char_T(&cur_seq);
  emxFree_char_T(&cur_line);
  emxInit_real_T(&M, 1);
  emxInit_int32_T(&D, 1);

  /* 'gkmPWM:196' fclose(fid); */
  cfclose(fileid);

  /*  [w, ind] = sort(a{2}, pn); */
  /*  s = a{1}(ind(1:min([100000 length(a{1})]))); */
  /* 'gkmPWM:201' [w, ind] = sort(alpha, pn); */
  sort(alpha, D);
  D_data = D->data;
  alpha_data = alpha->data;
  i = M->size[0];
  M->size[0] = D->size[0];
  emxEnsureCapacity_real_T(M, i);
  M_data = M->data;
  loop_ub = D->size[0];
  for (i = 0; i < loop_ub; i++) {
    M_data[i] = D_data[i];
  }

  emxInit_cell_wrap_0(&s, 1);

  /* 'gkmPWM:202' s_len = min([100000 length(sequences)]); */
  m[0] = 100000.0;
  m[1] = sequences->size[0];

  /* 'gkmPWM:203' s = cell(s_len, 1); */
  unnamed_idx_0_tmp_tmp = (int)b_minimum(m);
  i = s->size[0];
  s->size[0] = unnamed_idx_0_tmp_tmp;
  emxEnsureCapacity_cell_wrap_0(s, i);
  s_data = s->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    s_data[i].f1->size[0] = 1;
    s_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:204' s = coder.nullcopy(s); */
  /* 'gkmPWM:205' for cur_idx=1:s_len */
  for (cur_idx = 0; cur_idx < unnamed_idx_0_tmp_tmp; cur_idx++) {
    /* 'gkmPWM:206' s{cur_idx} = sequences{ind(cur_idx)}; */
    i = s_data[cur_idx].f1->size[0] * s_data[cur_idx].f1->size[1];
    s_data[cur_idx].f1->size[0] = 1;
    s_data[cur_idx].f1->size[1] = sequences_data[(int)M_data[cur_idx] - 1]
      .f1->size[1];
    emxEnsureCapacity_char_T(s_data[cur_idx].f1, i);
    loop_ub = sequences_data[(int)M_data[cur_idx] - 1].f1->size[1];
    for (i = 0; i < loop_ub; i++) {
      s_data[cur_idx].f1->data[i] = sequences_data[(int)M_data[cur_idx] - 1].
        f1->data[i];
    }
  }

  emxInit_cell_wrap_0(&b_p, 1);

  /* 'gkmPWM:210' l = length(s{1}); */
  varargin_2 = s_data[0].f1->size[1] - 1;
  l = s_data[0].f1->size[1] - 4;

  /* 'gkmPWM:211' k = round(l/2)+1; */
  x = (int)rt_roundd((double)s_data[0].f1->size[1] / 2.0);

  /* 'gkmPWM:212' ikl = length(ik); */
  /* 'gkmPWM:213' p = cell(num,1); */
  unnamed_idx_0_tmp_tmp = (int)num;
  i = b_p->size[0];
  b_p->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_0(b_p, i);
  p_data = b_p->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    p_data[i].f1->size[0] = 1;
    p_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:214' p = coder.nullcopy(p); */
  i = sequences->size[0];
  sequences->size[0] = b_p->size[0];
  emxEnsureCapacity_cell_wrap_0(sequences, i);
  sequences_data = sequences->data;

  /* 'gkmPWM:215' c = ikl+1; */
  *c = 1.0;

  /* 'gkmPWM:216' p{1} = s{1}; */
  i = sequences_data[0].f1->size[0] * sequences_data[0].f1->size[1];
  sequences_data[0].f1->size[0] = 1;
  sequences_data[0].f1->size[1] = s_data[0].f1->size[1];
  emxEnsureCapacity_char_T(sequences_data[0].f1, i);
  loop_ub = s_data[0].f1->size[1];
  emxFree_cell_wrap_0(&b_p);
  for (i = 0; i < loop_ub; i++) {
    sequences_data[0].f1->data[i] = s_data[0].f1->data[i];
  }

  emxInit_cell_wrap_1(&b_mat);
  emxInit_cell_wrap_2(&b_pwms);

  /* 'gkmPWM:217' mat = cell(ikl+num,1); */
  i = b_mat->size[0];
  b_mat->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_1(b_mat, i);
  mat_data = b_mat->data;

  /* 'gkmPWM:218' mat = coder.nullcopy(mat); */
  /*  mat(1:ikl) = ik; */
  /* 'gkmPWM:220' for cur_idx=1:length(ik) */
  /* 'gkmPWM:223' mat{c} = letterconvert(s{1}); */
  /* 'gkmPWM:224' pwms = cell(num,1); */
  i = b_pwms->size[0];
  b_pwms->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_2(b_pwms, i);
  pwms_data = b_pwms->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    mat_data[i].f1->size[0] = 1;
    mat_data[i].f1->size[1] = 0;
    pwms_data[i].f1->size[0] = 0;
    pwms_data[i].f1->size[1] = 4;
  }

  letterconvert(s_data[0].f1, mat_data[0].f1);

  /* 'gkmPWM:225' pwms = coder.nullcopy(pwms); */
  /* 'gkmPWM:226' for i = 1:num */
  for (b_i = 0; b_i < unnamed_idx_0_tmp_tmp; b_i++) {
    /* 'gkmPWM:227' pwms{i} = zeros(l,4); */
    i = pwms_data[b_i].f1->size[0] * pwms_data[b_i].f1->size[1];
    pwms_data[b_i].f1->size[0] = varargin_2 + 1;
    pwms_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(pwms_data[b_i].f1, i);
    loop_ub = (varargin_2 + 1) << 2;
    for (i = 0; i < loop_ub; i++) {
      pwms_data[b_i].f1->data[i] = 0.0;
    }
  }

  /* 'gkmPWM:229' for i = 1:l */
  i = s_data[0].f1->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:230' pwms{1}(i,mat{c}(i)+1) = pwms{1}(i,mat{c}(i)+1)+w(i); */
    d = mat_data[0].f1->data[b_i];
    pwms_data[0].f1->data[b_i + pwms_data[0].f1->size[0] * ((int)(d + 1.0) - 1)]
      += alpha_data[b_i];
  }

  /* 'gkmPWM:232' B = zeros(9,1); */
  /* 'gkmPWM:233' BB = zeros(9,1); */
  for (b_i = 0; b_i < 9; b_i++) {
    B[b_i] = 0;
    BB[b_i] = 0;
  }

  /* 'gkmPWM:234' B(1:5) = (0:4)'; */
  for (i = 0; i < 5; i++) {
    B[i] = (signed char)i;
  }

  /* 'gkmPWM:235' B(6:9) = 0; */
  for (b_i = 0; b_i < 4; b_i++) {
    B[b_i + 5] = 0;
  }

  /* 'gkmPWM:236' BB(1:5) = 0; */
  for (b_i = 0; b_i < 5; b_i++) {
    BB[b_i] = 0;
  }

  /* 'gkmPWM:237' BB(6:9) = (1:4)'; */
  for (i = 0; i < 4; i++) {
    BB[i + 5] = (signed char)(i + 1);
  }

  /* 'gkmPWM:238' CC = [l l-1 l-2 l-3 l-4 l-1 l-2 l-3 l-4]; */
  CC[0] = s_data[0].f1->size[1];
  CC[1] = s_data[0].f1->size[1] - 1;
  CC[2] = s_data[0].f1->size[1] - 2;
  CC[3] = s_data[0].f1->size[1] - 3;
  CC[4] = s_data[0].f1->size[1] - 4;
  CC[5] = s_data[0].f1->size[1] - 1;
  CC[6] = s_data[0].f1->size[1] - 2;
  CC[7] = s_data[0].f1->size[1] - 3;
  CC[8] = s_data[0].f1->size[1] - 4;

  /* this process picks kmers to seed the PWMs.  kmers that match one of the seeds by round(l/2)+1 or more are added that particular seed.  Otherwise, it becomes another seed. */
  /* 'gkmPWM:240' for i = 2:100000 */
  b_i = 1;
  emxInit_real_T(&ss, 2);
  emxInit_real_T(&rs, 2);
  emxInit_int8_T(&DD, 1);
  emxInit_boolean_T(&b_x, 2);
  exitg2 = false;
  while ((!exitg2) && (b_i - 1 < 99999)) {
    /* 'gkmPWM:241' ss = letterconvert(s{i}); */
    letterconvert(s_data[b_i].f1, ss);
    ss_data = ss->data;

    /* 'gkmPWM:242' rs = 3-fliplr(ss); */
    i = rs->size[0] * rs->size[1];
    rs->size[0] = 1;
    rs->size[1] = ss->size[1];
    emxEnsureCapacity_real_T(rs, i);
    rs_data = rs->data;
    loop_ub = ss->size[1];
    for (i = 0; i < loop_ub; i++) {
      rs_data[i] = ss_data[i];
    }

    fliplr(rs);
    i = rs->size[0] * rs->size[1];
    rs->size[0] = 1;
    emxEnsureCapacity_real_T(rs, i);
    rs_data = rs->data;
    loop_ub = rs->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      rs_data[i] = 3.0 - rs_data[i];
    }

    /* 'gkmPWM:243' M = zeros(c,1); */
    loop_ub_tmp = (int)*c;
    i = M->size[0];
    M->size[0] = (int)*c;
    emxEnsureCapacity_real_T(M, i);
    M_data = M->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      M_data[i] = 0.0;
    }

    /* 'gkmPWM:244' D = zeros(c,1); */
    i = D->size[0];
    D->size[0] = (int)*c;
    emxEnsureCapacity_int32_T(D, i);
    D_data = D->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      D_data[i] = 0;
    }

    /* 'gkmPWM:245' DD = zeros(c,1); */
    i = DD->size[0];
    DD->size[0] = (int)*c;
    emxEnsureCapacity_int8_T(DD, i);
    DD_data = DD->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      DD_data[i] = 0;
    }

    /* 'gkmPWM:246' for j = 1:c */
    for (j = 0; j < loop_ub_tmp; j++) {
      /* 'gkmPWM:247' [m,d] = max([sum(mat{j}==ss) sum(mat{j}(2:end)==ss(1:l-1)) sum(mat{j}(3:end)==ss(1:l-2)) sum(mat{j}(4:end)==ss(1:l-3)) sum(mat{j}(5:end)==ss(1:l-4)) sum(mat{j}(1:l-1)==ss(2:end)) sum(mat{j}(1:l-2)==ss(3:end)) sum(mat{j}(1:l-3)==ss(4:end)) sum(mat{j}(1:l-4)==ss(5:end))]); */
      if (2 > mat_data[j].f1->size[1]) {
        i = 0;
        b_y = 0;
      } else {
        i = 1;
        b_y = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2) {
        unnamed_idx_0_tmp_tmp = -1;
      } else {
        unnamed_idx_0_tmp_tmp = l + 2;
      }

      if (3 > mat_data[j].f1->size[1]) {
        cur_idx = 0;
        c_y = 0;
      } else {
        cur_idx = 2;
        c_y = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 1) {
        d_y = -1;
      } else {
        d_y = l + 1;
      }

      if (4 > mat_data[j].f1->size[1]) {
        e_y = 0;
        i1 = 0;
      } else {
        e_y = 3;
        i1 = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 2) {
        i2 = -1;
      } else {
        i2 = l;
      }

      if (5 > mat_data[j].f1->size[1]) {
        i3 = 0;
        i4 = 0;
      } else {
        i3 = 4;
        i4 = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 3) {
        i5 = 0;
      } else {
        i5 = l;
      }

      if (1 > varargin_2) {
        loop_ub = -1;
      } else {
        loop_ub = l + 2;
      }

      if (2 > ss->size[1]) {
        i6 = 0;
        i7 = 0;
      } else {
        i6 = 1;
        i7 = ss->size[1];
      }

      if (1 > varargin_2 - 1) {
        b_loop_ub = -1;
      } else {
        b_loop_ub = l + 1;
      }

      if (3 > ss->size[1]) {
        i8 = 0;
        i9 = 0;
      } else {
        i8 = 2;
        i9 = ss->size[1];
      }

      if (1 > varargin_2 - 2) {
        c_loop_ub = -1;
      } else {
        c_loop_ub = l;
      }

      if (4 > ss->size[1]) {
        i10 = 0;
        i11 = 0;
      } else {
        i10 = 3;
        i11 = ss->size[1];
      }

      if (1 > varargin_2 - 3) {
        d_loop_ub = 0;
      } else {
        d_loop_ub = l;
      }

      if (5 > ss->size[1]) {
        i12 = 0;
        i13 = 0;
      } else {
        i12 = 4;
        i13 = ss->size[1];
      }

      if (mat_data[j].f1->size[1] == ss->size[1]) {
        f_y = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = mat_data[j].f1->size[1];
        emxEnsureCapacity_boolean_T(b_x, f_y);
        x_data = b_x->data;
        e_loop_ub = mat_data[j].f1->size[1];
        for (f_y = 0; f_y < e_loop_ub; f_y++) {
          x_data[f_y] = (mat_data[j].f1->data[f_y] == ss_data[f_y]);
        }
      } else {
        l_binary_expand_op(b_x, b_mat, j, ss);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        g_y = 0;
      } else {
        g_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          g_y += x_data[k - 1];
        }
      }

      e_loop_ub = b_y - i;
      if (e_loop_ub == unnamed_idx_0_tmp_tmp + 1) {
        b_y = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, b_y);
        x_data = b_x->data;
        for (b_y = 0; b_y < e_loop_ub; b_y++) {
          x_data[b_y] = (mat_data[j].f1->data[i + b_y] == ss_data[b_y]);
        }
      } else {
        k_binary_expand_op(b_x, b_mat, j, i, b_y - 1, ss, unnamed_idx_0_tmp_tmp
                           - 2);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        b_y = 0;
      } else {
        b_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          b_y += x_data[k - 1];
        }
      }

      e_loop_ub = c_y - cur_idx;
      if (e_loop_ub == d_y + 1) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[cur_idx + i] == ss_data[i]);
        }
      } else {
        j_binary_expand_op(b_x, b_mat, j, cur_idx, c_y - 1, ss, d_y - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        f_y = 0;
      } else {
        f_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          f_y += x_data[k - 1];
        }
      }

      e_loop_ub = i1 - e_y;
      if (e_loop_ub == i2 + 1) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[e_y + i] == ss_data[i]);
        }
      } else {
        i_binary_expand_op(b_x, b_mat, j, e_y, i1 - 1, ss, i2);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        e_y = 0;
      } else {
        e_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          e_y += x_data[k - 1];
        }
      }

      e_loop_ub = i4 - i3;
      if (e_loop_ub == i5) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i3 + i] == ss_data[i]);
        }
      } else {
        i_binary_expand_op(b_x, b_mat, j, i3, i4 - 1, ss, i5 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        e_loop_ub = 0;
      } else {
        e_loop_ub = x_data[0];
        for (k = 2; k <= vlen; k++) {
          e_loop_ub += x_data[k - 1];
        }
      }

      if (loop_ub + 1 == i7 - i6) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == ss_data[i6 + i]);
        }
      } else {
        h_binary_expand_op(b_x, b_mat, j, loop_ub - 2, ss, i6, i7 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        d_y = 0;
      } else {
        d_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          d_y += x_data[k - 1];
        }
      }

      if (b_loop_ub + 1 == i9 - i8) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = b_loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= b_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == ss_data[i8 + i]);
        }
      } else {
        g_binary_expand_op(b_x, b_mat, j, b_loop_ub - 1, ss, i8, i9 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        c_y = 0;
      } else {
        c_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          c_y += x_data[k - 1];
        }
      }

      if (c_loop_ub + 1 == i11 - i10) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = c_loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= c_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == ss_data[i10 + i]);
        }
      } else {
        f_binary_expand_op(b_x, b_mat, j, c_loop_ub, ss, i10, i11 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        cur_idx = 0;
      } else {
        cur_idx = x_data[0];
        for (k = 2; k <= vlen; k++) {
          cur_idx += x_data[k - 1];
        }
      }

      if (d_loop_ub == i13 - i12) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = d_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < d_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == ss_data[i12 + i]);
        }
      } else {
        f_binary_expand_op(b_x, b_mat, j, d_loop_ub - 1, ss, i12, i13 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        unnamed_idx_0_tmp_tmp = 0;
      } else {
        unnamed_idx_0_tmp_tmp = x_data[0];
        for (k = 2; k <= vlen; k++) {
          unnamed_idx_0_tmp_tmp += x_data[k - 1];
        }
      }

      h_y[0] = g_y;
      h_y[1] = b_y;
      h_y[2] = f_y;
      h_y[3] = e_y;
      h_y[4] = e_loop_ub;
      h_y[5] = d_y;
      h_y[6] = c_y;
      h_y[7] = cur_idx;
      h_y[8] = unnamed_idx_0_tmp_tmp;
      b_maximum(h_y, &curr_pos, &iindx);

      /* 'gkmPWM:248' [mm,dd] = max([sum(mat{j}==rs) sum(mat{j}(2:end)==rs(1:l-1)) sum(mat{j}(3:end)==rs(1:l-2)) sum(mat{j}(4:end)==rs(1:l-3)) sum(mat{j}(5:end)==rs(1:l-4)) sum(mat{j}(1:l-1)==rs(2:end)) sum(mat{j}(1:l-2)==rs(3:end)) sum(mat{j}(1:l-3)==rs(4:end)) sum(mat{j}(1:l-4)==rs(5:end))]); */
      if (2 > mat_data[j].f1->size[1]) {
        i = 0;
        b_y = 0;
      } else {
        i = 1;
        b_y = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2) {
        unnamed_idx_0_tmp_tmp = -1;
      } else {
        unnamed_idx_0_tmp_tmp = l + 2;
      }

      if (3 > mat_data[j].f1->size[1]) {
        cur_idx = 0;
        c_y = 0;
      } else {
        cur_idx = 2;
        c_y = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 1) {
        d_y = -1;
      } else {
        d_y = l + 1;
      }

      if (4 > mat_data[j].f1->size[1]) {
        e_y = 0;
        i1 = 0;
      } else {
        e_y = 3;
        i1 = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 2) {
        i2 = -1;
      } else {
        i2 = l;
      }

      if (5 > mat_data[j].f1->size[1]) {
        i3 = 0;
        i4 = 0;
      } else {
        i3 = 4;
        i4 = mat_data[j].f1->size[1];
      }

      if (1 > varargin_2 - 3) {
        i5 = 0;
      } else {
        i5 = l;
      }

      if (1 > varargin_2) {
        loop_ub = -1;
      } else {
        loop_ub = l + 2;
      }

      if (2 > rs->size[1]) {
        i6 = 0;
        i7 = 0;
      } else {
        i6 = 1;
        i7 = rs->size[1];
      }

      if (1 > varargin_2 - 1) {
        b_loop_ub = -1;
      } else {
        b_loop_ub = l + 1;
      }

      if (3 > rs->size[1]) {
        i8 = 0;
        i9 = 0;
      } else {
        i8 = 2;
        i9 = rs->size[1];
      }

      if (1 > varargin_2 - 2) {
        c_loop_ub = -1;
      } else {
        c_loop_ub = l;
      }

      if (4 > rs->size[1]) {
        i10 = 0;
        i11 = 0;
      } else {
        i10 = 3;
        i11 = rs->size[1];
      }

      if (1 > varargin_2 - 3) {
        d_loop_ub = 0;
      } else {
        d_loop_ub = l;
      }

      if (5 > rs->size[1]) {
        i12 = 0;
        i13 = 0;
      } else {
        i12 = 4;
        i13 = rs->size[1];
      }

      if (mat_data[j].f1->size[1] == rs->size[1]) {
        f_y = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = mat_data[j].f1->size[1];
        emxEnsureCapacity_boolean_T(b_x, f_y);
        x_data = b_x->data;
        e_loop_ub = mat_data[j].f1->size[1];
        for (f_y = 0; f_y < e_loop_ub; f_y++) {
          x_data[f_y] = (mat_data[j].f1->data[f_y] == rs_data[f_y]);
        }
      } else {
        l_binary_expand_op(b_x, b_mat, j, rs);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        g_y = 0;
      } else {
        g_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          g_y += x_data[k - 1];
        }
      }

      e_loop_ub = b_y - i;
      if (e_loop_ub == unnamed_idx_0_tmp_tmp + 1) {
        b_y = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, b_y);
        x_data = b_x->data;
        for (b_y = 0; b_y < e_loop_ub; b_y++) {
          x_data[b_y] = (mat_data[j].f1->data[i + b_y] == rs_data[b_y]);
        }
      } else {
        k_binary_expand_op(b_x, b_mat, j, i, b_y - 1, rs, unnamed_idx_0_tmp_tmp
                           - 2);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        b_y = 0;
      } else {
        b_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          b_y += x_data[k - 1];
        }
      }

      e_loop_ub = c_y - cur_idx;
      if (e_loop_ub == d_y + 1) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[cur_idx + i] == rs_data[i]);
        }
      } else {
        j_binary_expand_op(b_x, b_mat, j, cur_idx, c_y - 1, rs, d_y - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        f_y = 0;
      } else {
        f_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          f_y += x_data[k - 1];
        }
      }

      e_loop_ub = i1 - e_y;
      if (e_loop_ub == i2 + 1) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[e_y + i] == rs_data[i]);
        }
      } else {
        i_binary_expand_op(b_x, b_mat, j, e_y, i1 - 1, rs, i2);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        e_y = 0;
      } else {
        e_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          e_y += x_data[k - 1];
        }
      }

      e_loop_ub = i4 - i3;
      if (e_loop_ub == i5) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = e_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < e_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i3 + i] == rs_data[i]);
        }
      } else {
        i_binary_expand_op(b_x, b_mat, j, i3, i4 - 1, rs, i5 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        e_loop_ub = 0;
      } else {
        e_loop_ub = x_data[0];
        for (k = 2; k <= vlen; k++) {
          e_loop_ub += x_data[k - 1];
        }
      }

      if (loop_ub + 1 == i7 - i6) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == rs_data[i6 + i]);
        }
      } else {
        h_binary_expand_op(b_x, b_mat, j, loop_ub - 2, rs, i6, i7 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        d_y = 0;
      } else {
        d_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          d_y += x_data[k - 1];
        }
      }

      if (b_loop_ub + 1 == i9 - i8) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = b_loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= b_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == rs_data[i8 + i]);
        }
      } else {
        g_binary_expand_op(b_x, b_mat, j, b_loop_ub - 1, rs, i8, i9 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        c_y = 0;
      } else {
        c_y = x_data[0];
        for (k = 2; k <= vlen; k++) {
          c_y += x_data[k - 1];
        }
      }

      if (c_loop_ub + 1 == i11 - i10) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = c_loop_ub + 1;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i <= c_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == rs_data[i10 + i]);
        }
      } else {
        f_binary_expand_op(b_x, b_mat, j, c_loop_ub, rs, i10, i11 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        cur_idx = 0;
      } else {
        cur_idx = x_data[0];
        for (k = 2; k <= vlen; k++) {
          cur_idx += x_data[k - 1];
        }
      }

      if (d_loop_ub == i13 - i12) {
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = d_loop_ub;
        emxEnsureCapacity_boolean_T(b_x, i);
        x_data = b_x->data;
        for (i = 0; i < d_loop_ub; i++) {
          x_data[i] = (mat_data[j].f1->data[i] == rs_data[i12 + i]);
        }
      } else {
        f_binary_expand_op(b_x, b_mat, j, d_loop_ub - 1, rs, i12, i13 - 1);
        x_data = b_x->data;
      }

      vlen = b_x->size[1];
      if (b_x->size[1] == 0) {
        unnamed_idx_0_tmp_tmp = 0;
      } else {
        unnamed_idx_0_tmp_tmp = x_data[0];
        for (k = 2; k <= vlen; k++) {
          unnamed_idx_0_tmp_tmp += x_data[k - 1];
        }
      }

      h_y[0] = g_y;
      h_y[1] = b_y;
      h_y[2] = f_y;
      h_y[3] = e_y;
      h_y[4] = e_loop_ub;
      h_y[5] = d_y;
      h_y[6] = c_y;
      h_y[7] = cur_idx;
      h_y[8] = unnamed_idx_0_tmp_tmp;
      b_maximum(h_y, &idx, &unnamed_idx_0_tmp_tmp);

      /* 'gkmPWM:249' [M(j),ddd] = max([m mm]); */
      m[0] = curr_pos;
      m[1] = idx;
      c_maximum(m, &M_data[j], &cur_idx);

      /* 'gkmPWM:250' if ddd == 1 */
      if (cur_idx == 1) {
        /* 'gkmPWM:251' D(j) = d; */
        D_data[j] = iindx;

        /* 'gkmPWM:252' DD(j) = 1; */
        DD_data[j] = 1;
      } else {
        /* 'gkmPWM:253' else */
        /* 'gkmPWM:254' D(j) = dd; */
        D_data[j] = unnamed_idx_0_tmp_tmp;

        /* 'gkmPWM:255' DD(j) = 2; */
        DD_data[j] = 2;
      }
    }

    /* 'gkmPWM:258' if max(M) < k */
    if (maximum(M) < (double)x + 1.0) {
      /* 'gkmPWM:259' c = c+1; */
      (*c)++;

      /* 'gkmPWM:260' p{c-ikl} = s{i}; */
      i = sequences_data[(int)*c - 1].f1->size[0] * sequences_data[(int)*c - 1].
        f1->size[1];
      sequences_data[(int)*c - 1].f1->size[0] = 1;
      sequences_data[(int)*c - 1].f1->size[1] = s_data[b_i].f1->size[1];
      emxEnsureCapacity_char_T(sequences_data[(int)*c - 1].f1, i);
      loop_ub = s_data[b_i].f1->size[1];
      for (i = 0; i < loop_ub; i++) {
        sequences_data[(int)*c - 1].f1->data[i] = s_data[b_i].f1->data[i];
      }

      /* 'gkmPWM:261' mat{c} = ss; */
      i = mat_data[(int)*c - 1].f1->size[0] * mat_data[(int)*c - 1].f1->size[1];
      mat_data[(int)*c - 1].f1->size[0] = 1;
      mat_data[(int)*c - 1].f1->size[1] = ss->size[1];
      emxEnsureCapacity_real_T(mat_data[(int)*c - 1].f1, i);
      loop_ub = ss->size[1];
      for (i = 0; i < loop_ub; i++) {
        mat_data[(int)*c - 1].f1->data[i] = ss_data[i];
      }

      /* 'gkmPWM:262' ss = ss+1; */
      i = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity_real_T(ss, i);
      ss_data = ss->data;
      loop_ub = ss->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        ss_data[i]++;
      }

      /* 'gkmPWM:263' for j = 1:l */
      for (j = 0; j <= varargin_2; j++) {
        /* 'gkmPWM:264' pwms{c-ikl}(j,ss(j)) = pwms{c-ikl}(j,ss(j))+w(i); */
        i = (int)ss_data[j] - 1;
        pwms_data[(int)*c - 1].f1->data[j + pwms_data[(int)*c - 1].f1->size[0] *
          i] += alpha_data[b_i];
      }
    } else {
      /* 'gkmPWM:266' else */
      /* 'gkmPWM:267' [~,d] = max(M); */
      d_maximum(M, &curr_pos, &iindx);

      /* 'gkmPWM:268' if DD(d) == 1 && d > ikl */
      i = DD_data[iindx - 1];
      if (i == 1) {
        /* 'gkmPWM:269' ss = ss+1; */
        i = ss->size[0] * ss->size[1];
        ss->size[0] = 1;
        emxEnsureCapacity_real_T(ss, i);
        ss_data = ss->data;
        loop_ub = ss->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          ss_data[i]++;
        }

        /* 'gkmPWM:270' d = d-ikl; */
        /* 'gkmPWM:271' for j = 1:CC(D(d)) */
        i = D_data[iindx - 1] - 1;
        b_y = CC[i];
        for (j = 0; j < b_y; j++) {
          /* 'gkmPWM:272' pwms{d}(j+B(D(d)),ss(j+BB(D(d)))) = pwms{d}(j+B(D(d)),ss(j+BB(D(d))))+w(i); */
          unnamed_idx_0_tmp_tmp = (int)((unsigned int)j + B[i]);
          cur_idx = (int)ss_data[(int)((unsigned int)j + BB[i])] - 1;
          pwms_data[iindx - 1].f1->data[unnamed_idx_0_tmp_tmp + pwms_data[iindx
            - 1].f1->size[0] * cur_idx] += alpha_data[b_i];
        }
      } else if (i == 2) {
        /* 'gkmPWM:274' elseif DD(d) == 2 && d > ikl */
        /* 'gkmPWM:275' rs = rs+1; */
        i = rs->size[0] * rs->size[1];
        rs->size[0] = 1;
        emxEnsureCapacity_real_T(rs, i);
        rs_data = rs->data;
        loop_ub = rs->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          rs_data[i]++;
        }

        /* 'gkmPWM:276' d = d-ikl; */
        /* 'gkmPWM:277' for j = 1:CC(D(d)) */
        i = D_data[iindx - 1] - 1;
        b_y = CC[i];
        for (j = 0; j < b_y; j++) {
          /* 'gkmPWM:278' pwms{d}(j+B(D(d)),rs(j+BB(D(d)))) = pwms{d}(j+B(D(d)),rs(j+BB(D(d))))+w(i); */
          unnamed_idx_0_tmp_tmp = (int)((unsigned int)j + B[i]);
          cur_idx = (int)rs_data[(int)((unsigned int)j + BB[i])] - 1;
          pwms_data[iindx - 1].f1->data[unnamed_idx_0_tmp_tmp + pwms_data[iindx
            - 1].f1->size[0] * cur_idx] += alpha_data[b_i];
        }
      }
    }

    /* 'gkmPWM:282' if c == num+ikl */
    if (*c == num) {
      exitg2 = true;
    } else {
      b_i++;
    }
  }

  emxFree_boolean_T(&b_x);
  emxFree_int8_T(&DD);
  emxFree_int32_T(&D);
  emxFree_real_T(&M);
  emxFree_real_T(&rs);
  emxFree_real_T(&ss);
  emxFree_cell_wrap_0(&s);
  emxFree_real_T(&alpha);

  /*  mat = mat(1:c); */
  /* 'gkmPWM:287' new_mat = cell(c, 1); */
  /* 'gkmPWM:288' for cur_idx=1:c */
  i = (int)*c;
  b_y = mat->size[0];
  mat->size[0] = (int)*c;
  emxEnsureCapacity_cell_wrap_1(mat, b_y);
  b_mat_data = mat->data;
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'gkmPWM:289' new_mat{cur_idx} = mat{cur_idx}; */
    b_y = b_mat_data[cur_idx].f1->size[0] * b_mat_data[cur_idx].f1->size[1];
    b_mat_data[cur_idx].f1->size[0] = 1;
    b_mat_data[cur_idx].f1->size[1] = mat_data[cur_idx].f1->size[1];
    emxEnsureCapacity_real_T(b_mat_data[cur_idx].f1, b_y);
    loop_ub = mat_data[cur_idx].f1->size[1];
    for (b_y = 0; b_y < loop_ub; b_y++) {
      b_mat_data[cur_idx].f1->data[b_y] = mat_data[cur_idx].f1->data[b_y];
    }
  }

  emxFree_cell_wrap_1(&b_mat);

  /* 'gkmPWM:291' mat = new_mat; */
  /*  p = p(1:c-ikl); */
  /* 'gkmPWM:294' p_len = c-ikl; */
  /* 'gkmPWM:295' new_p = cell(p_len, 1); */
  /* 'gkmPWM:296' for cur_idx=1:p_len */
  b_y = p->size[0] * p->size[1];
  p->size[0] = (int)*c;
  p->size[1] = 1;
  emxEnsureCapacity_cell_wrap_0(p, b_y);
  p_data = p->data;
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'gkmPWM:297' new_p{cur_idx} = p{cur_idx}; */
    b_y = p_data[cur_idx].f1->size[0] * p_data[cur_idx].f1->size[1];
    p_data[cur_idx].f1->size[0] = 1;
    p_data[cur_idx].f1->size[1] = sequences_data[cur_idx].f1->size[1];
    emxEnsureCapacity_char_T(p_data[cur_idx].f1, b_y);
    loop_ub = sequences_data[cur_idx].f1->size[1];
    for (b_y = 0; b_y < loop_ub; b_y++) {
      p_data[cur_idx].f1->data[b_y] = sequences_data[cur_idx].f1->data[b_y];
    }
  }

  emxFree_cell_wrap_0(&sequences);

  /* 'gkmPWM:299' p = new_p; */
  /*  pwms = pwms(1:c-ikl); */
  /* 'gkmPWM:302' pwms_len = c-ikl; */
  /* 'gkmPWM:303' new_pwms = cell(pwms_len, 1); */
  /* 'gkmPWM:304' for cur_idx=1:pwms_len */
  b_y = pwms->size[0];
  pwms->size[0] = (int)*c;
  emxEnsureCapacity_cell_wrap_2(pwms, b_y);
  b_pwms_data = pwms->data;
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'gkmPWM:305' new_pwms{cur_idx} = pwms{cur_idx}; */
    b_y = b_pwms_data[cur_idx].f1->size[0] * b_pwms_data[cur_idx].f1->size[1];
    b_pwms_data[cur_idx].f1->size[0] = pwms_data[cur_idx].f1->size[0];
    b_pwms_data[cur_idx].f1->size[1] = 4;
    emxEnsureCapacity_real_T(b_pwms_data[cur_idx].f1, b_y);
    loop_ub = pwms_data[cur_idx].f1->size[0] * 4;
    for (b_y = 0; b_y < loop_ub; b_y++) {
      b_pwms_data[cur_idx].f1->data[b_y] = pwms_data[cur_idx].f1->data[b_y];
    }
  }

  emxFree_cell_wrap_2(&b_pwms);

  /* 'gkmPWM:307' pwms = new_pwms; */
  /* 'gkmPWM:309' for i = 1:c-ikl */
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:310' for j = 1:l */
    for (j = 0; j <= varargin_2; j++) {
      /* 'gkmPWM:311' pwms{i}(j,:) = pwms{i}(j,:)/sum(pwms{i}(j,:)); */
      curr_pos = b_pwms_data[b_i].f1->data[j];
      for (k = 0; k < 3; k++) {
        curr_pos += b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0] *
          (k + 1)];
      }

      for (b_y = 0; b_y < 4; b_y++) {
        b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0] * b_y] /=
          curr_pos;
      }
    }
  }
}

static void w_binary_expand_op(double p_data[], int *p_size, const creal_T
  ps_data[], const int *ps_size, const creal_T t, const creal_T B_data[], const
  int B_size[2], const creal_T E)
{
  creal_T t_data[9];
  creal_T dcv[3];
  double b_t_re_tmp;
  double brm;
  double t_im;
  double t_re;
  double t_re_tmp;
  int i;
  int i1;
  int loop_ub;
  int t_re_tmp_tmp;
  loop_ub = B_size[0] * B_size[1];
  for (i = 0; i < loop_ub; i++) {
    t_re_tmp = B_data[i].im;
    b_t_re_tmp = B_data[i].re;
    t_re = t.re * b_t_re_tmp - t.im * t_re_tmp;
    t_im = t.re * t_re_tmp + t.im * b_t_re_tmp;
    if (E.im == 0.0) {
      if (t_im == 0.0) {
        t_data[i].re = t_re / E.re;
        t_data[i].im = 0.0;
      } else if (t_re == 0.0) {
        t_data[i].re = 0.0;
        t_data[i].im = t_im / E.re;
      } else {
        t_data[i].re = t_re / E.re;
        t_data[i].im = t_im / E.re;
      }
    } else if (E.re == 0.0) {
      if (t_re == 0.0) {
        t_data[i].re = t_im / E.im;
        t_data[i].im = 0.0;
      } else if (t_im == 0.0) {
        t_data[i].re = 0.0;
        t_data[i].im = -(t_re / E.im);
      } else {
        t_data[i].re = t_im / E.im;
        t_data[i].im = -(t_re / E.im);
      }
    } else {
      brm = fabs(E.re);
      t_re_tmp = fabs(E.im);
      if (brm > t_re_tmp) {
        t_re_tmp = E.im / E.re;
        b_t_re_tmp = E.re + t_re_tmp * E.im;
        t_data[i].re = (t_re + t_re_tmp * t_im) / b_t_re_tmp;
        t_data[i].im = (t_im - t_re_tmp * t_re) / b_t_re_tmp;
      } else if (t_re_tmp == brm) {
        if (E.re > 0.0) {
          t_re_tmp = 0.5;
        } else {
          t_re_tmp = -0.5;
        }

        if (E.im > 0.0) {
          b_t_re_tmp = 0.5;
        } else {
          b_t_re_tmp = -0.5;
        }

        t_data[i].re = (t_re * t_re_tmp + t_im * b_t_re_tmp) / brm;
        t_data[i].im = (t_im * t_re_tmp - t_re * b_t_re_tmp) / brm;
      } else {
        t_re_tmp = E.re / E.im;
        b_t_re_tmp = E.im + t_re_tmp * E.re;
        t_data[i].re = (t_re_tmp * t_re + t_im) / b_t_re_tmp;
        t_data[i].im = (t_re_tmp * t_im - t_re) / b_t_re_tmp;
      }
    }
  }

  for (i = 0; i < 3; i++) {
    dcv[i].re = 1.0;
    dcv[i].im = 0.0;
  }

  *p_size = 3;
  loop_ub = (*ps_size != 1);
  for (i = 0; i < 3; i++) {
    t_re_tmp = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      t_re_tmp_tmp = i + 3 * i1;
      t_re_tmp += t_data[t_re_tmp_tmp].re - t_data[t_re_tmp_tmp].im * dcv[i1].im;
    }

    p_data[i] = ps_data[i * loop_ub].re + t_re_tmp;
  }
}

static void x_binary_expand_op(creal_T x[4], const double MAT_data[], const int
  MAT_size[2], const creal_T E, const signed char b[4])
{
  creal_T MAT[4];
  int aux_0_1;
  int i;
  int i1;
  int i2;
  int i3;
  int stride_0_0;
  int stride_0_1;
  stride_0_0 = (MAT_size[0] != 1);
  stride_0_1 = (MAT_size[1] != 1);
  aux_0_1 = 0;
  for (i = 0; i < 2; i++) {
    for (i1 = 0; i1 < 2; i1++) {
      i2 = i1 + (i << 1);
      i3 = b[i2];
      MAT[i2].re = MAT_data[i1 * stride_0_0 + MAT_size[0] * aux_0_1] - E.re *
        (double)i3;
      MAT[i2].im = 0.0 - E.im * (double)i3;
    }

    aux_0_1 += stride_0_1;
  }

  memcpy(&x[0], &MAT[0], 4U * sizeof(creal_T));
}

static void y_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *b,
  const emxArray_real_T *A, int iindx)
{
  const double *A_data;
  const double *b_data;
  double *vec_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  A_data = A->data;
  b_data = b->data;
  i = A->size[0];
  stride_0_0 = vec->size[0];
  if (i == 1) {
    vec->size[0] = b->size[0];
  } else {
    vec->size[0] = i;
  }

  emxEnsureCapacity_real_T(vec, stride_0_0);
  vec_data = vec->data;
  stride_0_0 = (b->size[0] != 1);
  stride_1_0 = (i != 1);
  if (i == 1) {
    loop_ub = b->size[0];
  } else {
    loop_ub = i;
  }

  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = b_data[i * stride_0_0] - A_data[i * stride_1_0 + A->size[0] *
      (iindx - 1)];
  }
}

/*
 * function gkmPWM(varargin)
 */
void gkmPWM(const emxArray_char_T *varargin_1, const emxArray_char_T *varargin_2,
            const emxArray_char_T *varargin_3, double varargin_4, double
            varargin_5, double varargin_6, double varargin_7, double varargin_8,
            double varargin_9, double varargin_10, double varargin_11, double
            varargin_12, double varargin_13, double varargin_14)
{
  FILE* b_NULL;
  FILE* filestar;
  cell_wrap_0 *kmers2_data;
  cell_wrap_0 *kmers_data;
  cell_wrap_13 *p_data;
  cell_wrap_2 *b_p_data;
  cell_wrap_2 *pp_data;
  emxArray_cell_wrap_0 *kmers;
  emxArray_cell_wrap_0 *kmers2;
  emxArray_cell_wrap_1 *seed;
  emxArray_cell_wrap_1 *seed2;
  emxArray_cell_wrap_13 *p;
  emxArray_cell_wrap_2 *b_p;
  emxArray_cell_wrap_2 *pp;
  emxArray_char_T *b_charStr;
  emxArray_char_T *b_varargin_1;
  emxArray_char_T *c_varargin_1;
  emxArray_char_T *charStr;
  emxArray_real_T *C;
  emxArray_real_T *E;
  emxArray_real_T *R;
  emxArray_real_T *b_comb;
  emxArray_real_T *b_mat;
  emxArray_real_T *c2;
  emxArray_real_T *c_mat;
  emxArray_real_T *cfile;
  emxArray_real_T *comb;
  emxArray_real_T *dc;
  emxArray_real_T *dc2;
  emxArray_real_T *diffc;
  emxArray_real_T *indc;
  emxArray_real_T *negvec;
  emxArray_real_T *rc;
  emxArray_real_T *seqvec;
  emxArray_real_T *seqvec2;
  emxArray_real_T *xc;
  double mat[16];
  double mat2[16];
  double GCmat[4];
  double startvec[4];
  double lk_data[2];
  double GCmat_tmp;
  double GCneg1;
  double c_tmp;
  double len;
  double nfrac;
  double pnr;
  double rcnum;
  double tmp;
  double *b_comb_data;
  double *c2_data;
  double *cfile_data;
  double *comb_data;
  double *dc2_data;
  double *dc_data;
  double *negvec_data;
  double *seqvec2_data;
  double *seqvec_data;
  int lk_size[2];
  unsigned int a;
  int b_i;
  int b_loop_ub;
  int c_loop_ub;
  int i;
  int i1;
  int i2;
  int j2;
  int loop_ub;
  int m;
  int nd2;
  const char *varargin_1_data;
  signed char fileid;
  char *b_varargin_1_data;
  char *charStr_data;
  bool ipnr;
  if (!isInitialized_gkmPWM) {
    gkmPWM_initialize();
  }

  varargin_1_data = varargin_1->data;

  /*  gkmPWM finds predictive motifs from sequence-based models of regulatory */
  /*      elements de novo   */
  /*   */
  /*      gkmPWM(fileprefix, wfile, memefile, motifNum, ...) */
  /*   */
  /*      Works conveniently with the output of the gkmSVM R package or the lsgkm */
  /*      package (https://github.com/Dongwon-Lee/lsgkm).  Can be leveraged to  */
  /*      extract motifs from other sequence based models as long as you have scores  */
  /*      paired with sequences.  The sequences should be in fasta format with the  */
  /*      suffix *_svseq.fa.  The scores should be in a 2 column tab delimited  */
  /*      file with the suffix *svalpha.out with the same prefix as the *svseq.fa  */
  /*      file.  The first column containing the labels of each sequence and the  */
  /*      scores in the second.   */
  /*   */
  /* 'gkmPWM:18' fileprefix = varargin{1}; */
  /* 'gkmPWM:19' wfile = varargin{2}; */
  /* 'gkmPWM:20' memefile = varargin{3}; */
  /* 'gkmPWM:21' mnum = varargin{4}; */
  /* 'gkmPWM:22' num = varargin{5}; */
  /* 'gkmPWM:23' rcorr = varargin{6}; */
  /* 'gkmPWM:24' reg = varargin{7}; */
  /* 'gkmPWM:25' l_svm = varargin{8}; */
  /* 'gkmPWM:26' k_svm = varargin{9}; */
  /* 'gkmPWM:27' BG_GC = varargin{10}; */
  /* 'gkmPWM:28' RC = varargin{11}; */
  /* 'gkmPWM:29' pnr = varargin{12}; */
  pnr = varargin_12;

  /* 'gkmPWM:30' nfrac = varargin{13}; */
  /* 'gkmPWM:31' nfracLim = varargin{14}; */
  /* 'gkmPWM:33' if pnr == 0 */
  if (varargin_12 == 0.0) {
    /* 'gkmPWM:34' ipnr = true; */
    ipnr = true;
  } else {
    /* 'gkmPWM:35' else */
    /* 'gkmPWM:36' ipnr = false; */
    ipnr = false;
  }

  /* 'gkmPWM:38' lk = 1; */
  lk_size[0] = 1;
  lk_size[1] = 1;
  lk_data[0] = 1.0;

  /* 'gkmPWM:39' if nfrac ~= 1 */
  if (varargin_13 != 1.0) {
    /* 'gkmPWM:40' lk = 0; */
    lk_size[0] = 1;
    lk_size[1] = 1;
    lk_data[0] = 0.0;
  }

  emxInit_real_T(&comb, 2);
  emxInit_real_T(&rc, 2);
  emxInit_real_T(&diffc, 1);
  emxInit_real_T(&indc, 1);
  emxInit_real_T(&xc, 2);

  /* 'gkmPWM:44' [comb,rc,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm,nfrac); */
  genIndex(varargin_8, varargin_9, varargin_13, comb, rc, diffc, indc, xc,
           &rcnum);
  comb_data = comb->data;

  /* generate gapped positions, adjusted for reverse complements */
  /* 'gkmPWM:45' if nfracLim && numel(comb)/k_svm*4^k_svm > 5*10^5 */
  if (varargin_14 != 0.0) {
    nfrac = pow(4.0, varargin_9);
    if ((double)(comb->size[0] * comb->size[1]) / varargin_9 * nfrac > 500000.0)
    {
      /* 'gkmPWM:46' nfrac = round(5*10^7/4^k_svm/numel(comb)*k_svm)/100; */
      nfrac = rt_roundd(5.0E+7 / nfrac / (double)(comb->size[0] * comb->size[1])
                        * varargin_9) / 100.0;

      /* 'gkmPWM:48' lk = [l_svm k_svm]; */
      lk_size[0] = 1;
      lk_size[1] = 2;
      lk_data[0] = varargin_8;
      lk_data[1] = varargin_9;

      /* 'gkmPWM:49' [comb,rc,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm,nfrac); */
      genIndex(varargin_8, varargin_9, nfrac, comb, rc, diffc, indc, xc, &rcnum);
      comb_data = comb->data;
        
      printf("WARNING: Using %d gapped kmers\n", (int)((double)(comb->size[0] * comb->size[1]) / varargin_9 * pow(4.0, varargin_9)));
      fflush(stdout); 
    }
  }

  emxInit_char_T(&b_varargin_1, 2);

  /* 'gkmPWM:52' fprintf('Running gkmPWM on %s for %d motifs and %d iterations\n', fileprefix, int32(mnum), int32(num)); */
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  b_varargin_1_data = b_varargin_1->data;
  loop_ub = varargin_1->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_varargin_1_data[i] = varargin_1_data[i];
  }

  emxInit_real_T(&cfile, 1);
  b_varargin_1_data[varargin_1->size[1]] = '\x00';
  printf("Running gkmPWM on %s for %d motifs and %d iterations\n",
         &b_varargin_1_data[0], (int)rt_roundd(varargin_4), (int)rt_roundd
         (varargin_5));
  fflush(stdout);

  /* 'gkmPWM:53' fprintf('Counting gapped k-mers\n'); */
  printf("Counting gapped k-mers\n");
  fflush(stdout);

  /* 'gkmPWM:55' [A, GCpos1, GCneg1,mat,mat2] = getgkmcounts(fileprefix,l_svm,k_svm,lk,RC,comb,rcnum); */
  getgkmcounts(varargin_1, varargin_8, varargin_9, lk_data, lk_size, varargin_11,
               comb, rcnum, cfile, &nfrac, &GCneg1, mat, mat2);
  cfile_data = cfile->data;

  /* count gapped k-mers, get GC content, and dinucleotide distribution */
  /* 'gkmPWM:56' if BG_GC == 1 */
  if (varargin_10 == 1.0) {
    /* 'gkmPWM:57' mat = (mat+mat2)/2; */
    for (i = 0; i < 16; i++) {
      mat[i] = (mat[i] + mat2[i]) / 2.0;
    }

    /* 'gkmPWM:58' GCpos1 = (GCpos1+GCneg1)/2; */
    GCneg1 = (nfrac + GCneg1) / 2.0;

    /* 'gkmPWM:59' GCneg1 = GCpos1; */
  }

  emxInit_real_T(&negvec, 1);

  /* 'gkmPWM:61' negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC); */
  /* 'BGkmer:2' len = numel(c)/k; */
  len = (double)(comb->size[0] * comb->size[1]) / varargin_9;

  /* 'BGkmer:3' alen = len-rcnum; */
  /* 'BGkmer:4' negweights = zeros(len*4^k,1); */
  c_tmp = pow(4.0, varargin_9);
  loop_ub = (int)(len * c_tmp);
  i = negvec->size[0];
  negvec->size[0] = loop_ub;
  emxEnsureCapacity_real_T(negvec, i);
  negvec_data = negvec->data;
  for (i = 0; i < loop_ub; i++) {
    negvec_data[i] = 0.0;
  }

  emxInit_real_T(&c2, 2);
  emxInit_real_T(&seqvec2, 2);
  emxInit_real_T(&dc2, 2);
  emxInit_cell_wrap_13(&p);

  /* 'BGkmer:5' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  GCmat_tmp = 0.5 - GCneg1 / 2.0;
  GCmat[0] = GCmat_tmp;
  GCmat[1] = GCneg1 / 2.0;
  GCmat[2] = GCneg1 / 2.0;
  GCmat[3] = GCmat_tmp;

  /* 'BGkmer:6' c2 = 0; */
  i = c2->size[0] * c2->size[1];
  c2->size[0] = 1;
  c2->size[1] = 1;
  emxEnsureCapacity_real_T(c2, i);
  c2_data = c2->data;
  c2_data[0] = 0.0;

  /*  To appease Coder */
  /* 'BGkmer:7' seqvec2 = 0; */
  i = seqvec2->size[0] * seqvec2->size[1];
  seqvec2->size[0] = 1;
  seqvec2->size[1] = 1;
  emxEnsureCapacity_real_T(seqvec2, i);
  seqvec2_data = seqvec2->data;
  seqvec2_data[0] = 0.0;

  /*  To appease Coder */
  /* 'BGkmer:8' dc2 = 0; */
  i = dc2->size[0] * dc2->size[1];
  dc2->size[0] = 1;
  dc2->size[1] = 1;
  emxEnsureCapacity_real_T(dc2, i);
  dc2_data = dc2->data;
  dc2_data[0] = 0.0;

  /*  To appease Coder */
  /* 'BGkmer:9' tmp = 0; */
  tmp = 0.0;

  /*  To appease Coder */
  /* 'BGkmer:10' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)varargin_8;
  emxEnsureCapacity_cell_wrap_13(p, i);
  p_data = p->data;

  /* 'BGkmer:11' p = coder.nullcopy(p); */
  /* 'BGkmer:12' p{1} = eye(4); */
  for (i = 0; i < 16; i++) {
    p_data[0].f1[i] = 0.0;
  }

  for (m = 0; m < 4; m++) {
    p_data[0].f1[m + (m << 2)] = 1.0;
  }

  /* 'BGkmer:13' for i = 1:l-1 */
  i = (int)(varargin_8 - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'BGkmer:14' p{i+1} = p{i}*mat; */
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < 4; i2++) {
        nfrac = 0.0;
        for (m = 0; m < 4; m++) {
          nfrac += p_data[b_i].f1[i1 + (m << 2)] * mat[m + (i2 << 2)];
        }

        mat2[i1 + (i2 << 2)] = nfrac;
      }
    }

    for (i1 = 0; i1 < 16; i1++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i1] = mat2[i1];
    }
  }

  emxInit_real_T(&seqvec, 2);

  /* 'BGkmer:16' seqvec = zeros(4^k, k); */
  i = seqvec->size[0] * seqvec->size[1];
  i1 = (int)c_tmp;
  seqvec->size[0] = (int)c_tmp;
  i2 = (int)varargin_9;
  seqvec->size[1] = (int)varargin_9;
  emxEnsureCapacity_real_T(seqvec, i);
  seqvec_data = seqvec->data;
  loop_ub = (int)c_tmp * (int)varargin_9;
  for (i = 0; i < loop_ub; i++) {
    seqvec_data[i] = 0.0;
  }

  /* 'BGkmer:17' for i = 1:k */
  for (b_i = 0; b_i < i2; b_i++) {
    /* 'BGkmer:18' for j = 1:4^k */
    for (m = 0; m < i1; m++) {
      /* 'BGkmer:19' seqvec(j,i) = mod(floor((j-1)/4^(i-1)), 4)+1; */
      seqvec_data[m + seqvec->size[0] * b_i] = b_mod(floor((((double)m + 1.0) -
        1.0) / pow(4.0, ((double)b_i + 1.0) - 1.0)), 4.0) + 1.0;
    }
  }

  /* 'BGkmer:22' if RC */
  if (varargin_11 != 0.0) {
    /* 'BGkmer:23' seqvec2 = 5-fliplr(seqvec); */
    i = seqvec2->size[0] * seqvec2->size[1];
    seqvec2->size[0] = seqvec->size[0];
    seqvec2->size[1] = seqvec->size[1];
    emxEnsureCapacity_real_T(seqvec2, i);
    seqvec2_data = seqvec2->data;
    loop_ub = seqvec->size[0] * seqvec->size[1];
    for (i = 0; i < loop_ub; i++) {
      seqvec2_data[i] = seqvec_data[i];
    }

    m = seqvec->size[0];
    nd2 = seqvec->size[1] >> 1;
    for (loop_ub = 0; loop_ub < nd2; loop_ub++) {
      j2 = (seqvec->size[1] - loop_ub) - 1;
      for (b_i = 0; b_i < m; b_i++) {
        nfrac = seqvec2_data[b_i + seqvec2->size[0] * loop_ub];
        seqvec2_data[b_i + seqvec2->size[0] * loop_ub] = seqvec2_data[b_i +
          seqvec2->size[0] * j2];
        seqvec2_data[b_i + seqvec2->size[0] * j2] = nfrac;
      }
    }

    loop_ub = seqvec2->size[0] * seqvec2->size[1];
    for (i = 0; i < loop_ub; i++) {
      seqvec2_data[i] = 5.0 - seqvec2_data[i];
    }

    /* 'BGkmer:24' c2 = l+1-fliplr(c); */
    i = c2->size[0] * c2->size[1];
    c2->size[0] = comb->size[0];
    c2->size[1] = comb->size[1];
    emxEnsureCapacity_real_T(c2, i);
    c2_data = c2->data;
    loop_ub = comb->size[0] * comb->size[1];
    for (i = 0; i < loop_ub; i++) {
      c2_data[i] = comb_data[i];
    }

    m = comb->size[0];
    nd2 = comb->size[1] >> 1;
    for (loop_ub = 0; loop_ub < nd2; loop_ub++) {
      j2 = (comb->size[1] - loop_ub) - 1;
      for (b_i = 0; b_i < m; b_i++) {
        nfrac = c2_data[b_i + c2->size[0] * loop_ub];
        c2_data[b_i + c2->size[0] * loop_ub] = c2_data[b_i + c2->size[0] * j2];
        c2_data[b_i + c2->size[0] * j2] = nfrac;
      }
    }

    loop_ub = c2->size[0] * c2->size[1];
    for (i = 0; i < loop_ub; i++) {
      c2_data[i] = (varargin_8 + 1.0) - c2_data[i];
    }
  }

  /* 'BGkmer:26' a = 1; */
  a = 1U;

  /* 'BGkmer:27' for i = 1:len */
  i = (int)len;
  emxInit_real_T(&dc, 2);
  emxInit_real_T(&b_comb, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'BGkmer:28' dc = diff(c(i,:)); */
    loop_ub = comb->size[1];
    i1 = b_comb->size[0] * b_comb->size[1];
    b_comb->size[0] = 1;
    b_comb->size[1] = comb->size[1];
    emxEnsureCapacity_real_T(b_comb, i1);
    b_comb_data = b_comb->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_comb_data[i1] = comb_data[b_i + comb->size[0] * i1];
    }

    diff(b_comb, dc);
    dc_data = dc->data;

    /* 'BGkmer:29' startvec = GCmat*p{c(i,1)}; */
    for (i1 = 0; i1 < 16; i1++) {
      mat2[i1] = p_data[(int)comb_data[b_i] - 1].f1[i1];
    }

    for (i1 = 0; i1 < 4; i1++) {
      nfrac = 0.0;
      for (i2 = 0; i2 < 4; i2++) {
        nfrac += GCmat[i2] * mat2[i2 + (i1 << 2)];
      }

      startvec[i1] = nfrac;
    }

    /* 'BGkmer:30' if RC */
    if (varargin_11 != 0.0) {
      /* 'BGkmer:31' dc2 = diff(c2(i,:)); */
      loop_ub = c2->size[1];
      i1 = b_comb->size[0] * b_comb->size[1];
      b_comb->size[0] = 1;
      b_comb->size[1] = c2->size[1];
      emxEnsureCapacity_real_T(b_comb, i1);
      b_comb_data = b_comb->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1] = c2_data[b_i + c2->size[0] * i1];
      }

      diff(b_comb, dc2);
      dc2_data = dc2->data;

      /* 'BGkmer:32' startvec2 = GCmat*p{c2(i,1)}; */
    }

    /* 'BGkmer:34' for ii = 1:4^k */
    i1 = (int)pow(4.0, varargin_9);
    for (nd2 = 0; nd2 < i1; nd2++) {
      /* 'BGkmer:35' negweights(a) = startvec(seqvec(ii,1)); */
      m = (int)(a + nd2) - 1;
      negvec_data[m] = startvec[(int)seqvec_data[nd2] - 1];

      /* 'BGkmer:36' if RC */
      if (varargin_11 != 0.0) {
        /* 'BGkmer:37' tmp = startvec(seqvec2(ii,1)); */
        tmp = startvec[(int)seqvec2_data[nd2] - 1];
      }

      /* 'BGkmer:39' for iii = 1:k-1 */
      i2 = (int)(varargin_9 - 1.0);
      for (j2 = 0; j2 < i2; j2++) {
        /* 'BGkmer:40' matt = p{dc(iii)+1}; */
        /* 'BGkmer:41' negweights(a) = negweights(a)*matt(seqvec(ii,iii), seqvec(ii,iii+1)); */
        negvec_data[m] *= p_data[(int)(dc_data[j2] + 1.0) - 1].f1[((int)
          seqvec_data[nd2 + seqvec->size[0] * j2] + (((int)seqvec_data[nd2 +
          seqvec->size[0] * (j2 + 1)] - 1) << 2)) - 1];

        /* 'BGkmer:42' if RC */
        if (varargin_11 != 0.0) {
          /* 'BGkmer:43' matt2 = p{dc2(iii)+1}; */
          /* 'BGkmer:44' tmp = tmp*matt2(seqvec2(ii,iii), seqvec2(ii,iii+1)); */
          tmp *= p_data[(int)(dc2_data[j2] + 1.0) - 1].f1[((int)seqvec2_data[nd2
            + seqvec2->size[0] * j2] + (((int)seqvec2_data[nd2 + seqvec2->size[0]
            * (j2 + 1)] - 1) << 2)) - 1];
        }
      }

      /* 'BGkmer:47' if RC */
      if (varargin_11 != 0.0) {
        /* 'BGkmer:48' negweights(a) = (negweights(a)+tmp); */
        negvec_data[m] += tmp;
      }

      /* 'BGkmer:50' a = a+1; */
    }

    a += i1;
  }

  emxFree_cell_wrap_13(&p);
  emxFree_real_T(&seqvec);
  emxFree_real_T(&seqvec2);
  emxFree_real_T(&c2);
  emxFree_real_T(&comb);

  /* 'BGkmer:53' if rcnum > 0 && RC */
  if ((rcnum > 0.0) && (varargin_11 != 0.0)) {
    /* 'BGkmer:54' negweights(4^k*alen+1:end) = negweights(4^k*alen+1:end)/sqrt(2); */
    nfrac = c_tmp * (len - rcnum) + 1.0;
    if (nfrac > negvec->size[0]) {
      i = 1;
      i1 = -1;
      i2 = 0;
    } else {
      i = (int)nfrac;
      i1 = (int)nfrac - 2;
      i2 = negvec->size[0];
    }

    m = (i2 - i1) - 1;
    i2 = b_comb->size[0] * b_comb->size[1];
    b_comb->size[0] = 1;
    b_comb->size[1] = m;
    emxEnsureCapacity_real_T(b_comb, i2);
    b_comb_data = b_comb->data;
    for (i2 = 0; i2 < m; i2++) {
      b_comb_data[i2] = negvec_data[(i + i2) - 1] / 1.4142135623730951;
    }

    loop_ub = b_comb->size[1];
    for (i = 0; i < loop_ub; i++) {
      negvec_data[(i1 + i) + 1] = b_comb_data[i];
    }
  }

  emxFree_real_T(&b_comb);

  /* generate expected gapped kmer distribution of background */
  /* 'gkmPWM:62' GC=[0.5-GCneg1/2 GCneg1/2 GCneg1/2 0.5-GCneg1/2]; */
  GCmat[0] = GCmat_tmp;
  GCmat[1] = GCneg1 / 2.0;
  GCmat[2] = GCneg1 / 2.0;
  GCmat[3] = GCmat_tmp;

  /* GC content vector */
  /* 'gkmPWM:63' if ipnr */
  if (ipnr) {
    /* 'gkmPWM:64' pnr = abs(max(A)/min(A)); */
    pnr = fabs(maximum(cfile) / minimum(cfile));
  }

  emxInit_cell_wrap_0(&kmers, 2);
  emxInit_cell_wrap_2(&b_p);
  emxInit_cell_wrap_1(&seed);

  /* 'gkmPWM:66' fprintf('Finding PWM seeds\n'); */
  printf("Finding PWM seeds\n");
  fflush(stdout);

  /* get PWM seeds using the kmer weight vectors */
  /* 'gkmPWM:68' coder.varsize('p'); */
  /* 'gkmPWM:69' [kmers, seed,p,c] = seed_kmers(wfile, mnum,'descend', {}); */
  seed_kmers(varargin_2, varargin_4, kmers, seed, b_p, &len);
  b_p_data = b_p->data;
  kmers_data = kmers->data;

  /* 'gkmPWM:70' if pnr ~= 0 */
  emxInit_cell_wrap_2(&pp);
  if (pnr != 0.0) {
    emxInit_cell_wrap_0(&kmers2, 1);
    emxInit_cell_wrap_1(&seed2);

    /* 'gkmPWM:71' [kmers2, seed2,pp, c2] = seed_kmers(wfile, max([floor(mnum/pnr) 2]),'ascend',seed); */
    nfrac = floor(varargin_4 / pnr);
    if (nfrac < 2.0) {
      nfrac = 2.0;
    }

    b_seed_kmers(varargin_2, nfrac, seed, kmers2, seed2, pp, &len);
    pp_data = pp->data;
    kmers2_data = kmers2->data;

    /*      kmers = [kmers;kmers2]; */
    /*      p = [p;pp]; */
    /* 'gkmPWM:74' for cur_idx=1:length(kmers2) */
    i = kmers2->size[0];
    emxFree_cell_wrap_1(&seed2);
    for (nd2 = 0; nd2 < i; nd2++) {
      /* 'gkmPWM:75' kmers{end + 1} = kmers2{cur_idx}; */
      m = kmers->size[0] + 1;
      i1 = kmers->size[0] * kmers->size[1];
      kmers->size[0]++;
      kmers->size[1] = 1;
      emxEnsureCapacity_cell_wrap_0(kmers, i1);
      kmers_data = kmers->data;
      loop_ub = kmers2_data[nd2].f1->size[1] - 1;
      i1 = kmers_data[m - 1].f1->size[0] * kmers_data[m - 1].f1->size[1];
      kmers_data[m - 1].f1->size[0] = 1;
      emxEnsureCapacity_char_T(kmers_data[m - 1].f1, i1);
      i1 = kmers_data[kmers->size[0] - 1].f1->size[0] * kmers_data[kmers->size[0]
        - 1].f1->size[1];
      kmers_data[kmers->size[0] - 1].f1->size[1] = kmers2_data[nd2].f1->size[1];
      emxEnsureCapacity_char_T(kmers_data[kmers->size[0] - 1].f1, i1);
      m = kmers->size[0] - 1;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        kmers_data[m].f1->data[i1] = kmers2_data[nd2].f1->data[i1];
      }
    }

    emxFree_cell_wrap_0(&kmers2);

    /* 'gkmPWM:77' for cur_idx=1:length(pp) */
    i = pp->size[0];
    for (nd2 = 0; nd2 < i; nd2++) {
      /* 'gkmPWM:78' p{end + 1} = pp{cur_idx}; */
      m = b_p->size[0] + 1;
      i1 = b_p->size[0];
      b_p->size[0]++;
      emxEnsureCapacity_cell_wrap_2(b_p, i1);
      b_p_data = b_p->data;
      loop_ub = pp_data[nd2].f1->size[0] * 4;
      i1 = b_p_data[m - 1].f1->size[0] * b_p_data[m - 1].f1->size[1];
      b_p_data[m - 1].f1->size[0] = pp_data[nd2].f1->size[0];
      emxEnsureCapacity_real_T(b_p_data[m - 1].f1, i1);
      i1 = b_p_data[b_p->size[0] - 1].f1->size[0] * b_p_data[b_p->size[0] - 1].
        f1->size[1];
      b_p_data[b_p->size[0] - 1].f1->size[1] = 4;
      emxEnsureCapacity_real_T(b_p_data[b_p->size[0] - 1].f1, i1);
      m = b_p->size[0] - 1;
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_p_data[m].f1->data[i1] = pp_data[nd2].f1->data[i1];
      }
    }

    /* 'gkmPWM:80' tot = c2; */
  } else {
    /* 'gkmPWM:81' else */
    /* 'gkmPWM:82' tot = c; */
  }

  emxFree_cell_wrap_1(&seed);

  /* 'gkmPWM:84' fprintf('Seeding PWMs at the following kmers\n'); */
  printf("Seeding PWMs at the following kmers\n");
  fflush(stdout);

  /* 'gkmPWM:85' for i = 1:tot */
  i = (int)len;
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:86' fprintf("%s\n", kmers{i}); */
    i1 = b_varargin_1->size[0] * b_varargin_1->size[1];
    b_varargin_1->size[0] = 1;
    b_varargin_1->size[1] = kmers_data[b_i].f1->size[1] + 1;
    emxEnsureCapacity_char_T(b_varargin_1, i1);
    b_varargin_1_data = b_varargin_1->data;
    loop_ub = kmers_data[b_i].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_varargin_1_data[i1] = kmers_data[b_i].f1->data[i1];
    }

    b_varargin_1_data[kmers_data[b_i].f1->size[1]] = '\x00';
    printf("%s\n", &b_varargin_1_data[0]);
    fflush(stdout);
  }

  emxFree_cell_wrap_0(&kmers);

  /* 'gkmPWM:89' for i = 1:length(p) */
  i = b_p->size[0];
  emxInit_real_T(&b_mat, 2);
  comb_data = b_mat->data;
  if (0 <= b_p->size[0] - 1) {
    b_repmat(GCmat, varargin_8 + 1.0, b_mat);
    comb_data = b_mat->data;
    b_loop_ub = b_mat->size[0];
    c_loop_ub = b_mat->size[0];
  }

  emxInit_real_T(&c_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:90' p{i} = extendPWM(p{i},l_svm+1,GC); */
    /* 'gkmPWM:971' mat = repmat(GCmat, n,1); */
    /* 'gkmPWM:972' ext_pwm = [mat;pwm;mat]; */
    loop_ub = b_p_data[b_i].f1->size[0];
    i1 = c_mat->size[0] * c_mat->size[1];
    c_mat->size[0] = (b_mat->size[0] + b_p_data[b_i].f1->size[0]) + b_mat->size
      [0];
    c_mat->size[1] = 4;
    emxEnsureCapacity_real_T(c_mat, i1);
    c2_data = c_mat->data;
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        c2_data[i2 + c_mat->size[0] * i1] = comb_data[i2 + b_mat->size[0] * i1];
      }

      for (i2 = 0; i2 < loop_ub; i2++) {
        c2_data[(i2 + b_mat->size[0]) + c_mat->size[0] * i1] = b_p_data[b_i].
          f1->data[i2 + b_p_data[b_i].f1->size[0] * i1];
      }
    }

    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < c_loop_ub; i2++) {
        c2_data[((i2 + b_mat->size[0]) + b_p_data[b_i].f1->size[0]) +
          c_mat->size[0] * i1] = comb_data[i2 + b_mat->size[0] * i1];
      }
    }

    i1 = b_p_data[b_i].f1->size[0] * b_p_data[b_i].f1->size[1];
    b_p_data[b_i].f1->size[0] = c_mat->size[0];
    b_p_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(b_p_data[b_i].f1, i1);
    loop_ub = c_mat->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_p_data[b_i].f1->data[i1] = c2_data[i1];
    }
  }

  emxFree_real_T(&c_mat);
  emxFree_real_T(&b_mat);

  /* 'gkmPWM:93' fprintf('Running de novo motif discovery\n'); */
  printf("Running de novo motif discovery\n");
  fflush(stdout);

  /* 'gkmPWM:94' m = mean(A); */
  /* 'gkmPWM:95' s = std(A); */
  /* 'gkmPWM:96' cfile = A-negvec/sum(negvec)*sum(A); */
  nfrac = blockedSummation(negvec, negvec->size[0]);
  len = blockedSummation(cfile, cfile->size[0]);
  if (cfile->size[0] == negvec->size[0]) {
    loop_ub = cfile->size[0];
    for (i = 0; i < loop_ub; i++) {
      cfile_data[i] -= negvec_data[i] / nfrac * len;
    }
  } else {
    binary_expand_op(cfile, negvec, nfrac, len);
    cfile_data = cfile->data;
  }

  /* cfile = cfile/max(abs(cfile));% normalize to speed up computation */
  /* 'gkmPWM:98' cfile = cfile/std(cfile); */
  nfrac = b_std(cfile);
  loop_ub = cfile->size[0];
  for (i = 0; i < loop_ub; i++) {
    cfile_data[i] /= nfrac;
  }

  emxInit_real_T(&C, 1);
  emxInit_real_T(&R, 1);
  emxInit_real_T(&E, 1);

  /*  clear A */
  /* 'gkmPWM:100' [pp, scorevec, C, r, R, E, Rd] = gkmPWM_lagrange(cfile,mat,p,negvec,num,rcorr,reg,l_svm,k_svm,RC,rc,diffc,indc,xc,rcnum); */
  gkmPWM_lagrange(cfile, mat, b_p, negvec, varargin_5, varargin_6, varargin_7,
                  varargin_8, varargin_9, varargin_11, rc, diffc, indc, xc,
                  rcnum, pp, dc2, C, &nfrac, R, E, dc);
  dc2_data = dc2->data;

  /*  createMEME([fileprefix '_' num2str(l_svm) '_' num2str(k_svm) '_' num2str(reg) '_' num2str(mnum)], pp, memefile,GCneg1, C, r, R, rcorr, E, Rd); */
  /* 'gkmPWM:103' createMEME(sprintf("%s_%d_%d_%d_%d", fileprefix, int32(l_svm), int32(k_svm), int32(reg), int32(mnum)), pp, memefile, GCneg1, C, r, R, rcorr, E, Rd); */
  nd2 = (int)rt_roundd(varargin_8);
  j2 = (int)rt_roundd(varargin_9);
  b_loop_ub = (int)rt_roundd(varargin_7);
  c_loop_ub = (int)rt_roundd(varargin_4);
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  b_varargin_1_data = b_varargin_1->data;
  loop_ub = varargin_1->size[1];
  emxFree_real_T(&cfile);
  emxFree_cell_wrap_2(&b_p);
  emxFree_real_T(&negvec);
  emxFree_real_T(&xc);
  emxFree_real_T(&indc);
  emxFree_real_T(&diffc);
  emxFree_real_T(&rc);
  for (i = 0; i < loop_ub; i++) {
    b_varargin_1_data[i] = varargin_1_data[i];
  }

  emxInit_char_T(&c_varargin_1, 2);
  b_varargin_1_data[varargin_1->size[1]] = '\x00';
  i = c_varargin_1->size[0] * c_varargin_1->size[1];
  c_varargin_1->size[0] = 1;
  c_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_1, i);
  charStr_data = c_varargin_1->data;
  loop_ub = varargin_1->size[1];
  for (i = 0; i < loop_ub; i++) {
    charStr_data[i] = varargin_1_data[i];
  }

  emxInit_char_T(&charStr, 2);
  charStr_data[varargin_1->size[1]] = '\x00';
  m = snprintf(NULL, 0, "%s_%d_%d_%d_%d", &charStr_data[0], nd2, j2, b_loop_ub,
               c_loop_ub);
  i = charStr->size[0] * charStr->size[1];
  charStr->size[0] = 1;
  charStr->size[1] = m + 1;
  emxEnsureCapacity_char_T(charStr, i);
  charStr_data = charStr->data;
  snprintf(&charStr_data[0], (size_t)(m + 1), "%s_%d_%d_%d_%d",
           &b_varargin_1_data[0], nd2, j2, b_loop_ub, c_loop_ub);
  i = charStr->size[0] * charStr->size[1];
  if (1 > m) {
    charStr->size[1] = 0;
  } else {
    charStr->size[1] = m;
  }

  emxEnsureCapacity_char_T(charStr, i);
  createMEME(charStr, pp, varargin_3, GCneg1, C, nfrac, R, varargin_6, E, dc);

  /*  dlmwrite([fileprefix '_' num2str(l_svm) '_' num2str(k_svm) '_' num2str(reg) '_' num2str(mnum) '_error.out'],scorevec'); */
  /* 'gkmPWM:106' fidw = fopen(sprintf("%s_%d_%d_%d_%d_error.out", fileprefix, int32(l_svm), int32(k_svm), int32(reg), int32(mnum)), 'w'); */
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  b_varargin_1_data = b_varargin_1->data;
  loop_ub = varargin_1->size[1];
  emxFree_real_T(&dc);
  emxFree_char_T(&charStr);
  emxFree_real_T(&E);
  emxFree_real_T(&R);
  emxFree_real_T(&C);
  emxFree_cell_wrap_2(&pp);
  for (i = 0; i < loop_ub; i++) {
    b_varargin_1_data[i] = varargin_1_data[i];
  }

  b_varargin_1_data[varargin_1->size[1]] = '\x00';
  i = c_varargin_1->size[0] * c_varargin_1->size[1];
  c_varargin_1->size[0] = 1;
  c_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(c_varargin_1, i);
  charStr_data = c_varargin_1->data;
  loop_ub = varargin_1->size[1];
  for (i = 0; i < loop_ub; i++) {
    charStr_data[i] = varargin_1_data[i];
  }

  emxInit_char_T(&b_charStr, 2);
  charStr_data[varargin_1->size[1]] = '\x00';
  m = snprintf(NULL, 0, "%s_%d_%d_%d_%d_error.out", &charStr_data[0], nd2, j2,
               b_loop_ub, c_loop_ub);
  i = b_charStr->size[0] * b_charStr->size[1];
  b_charStr->size[0] = 1;
  b_charStr->size[1] = m + 1;
  emxEnsureCapacity_char_T(b_charStr, i);
  charStr_data = b_charStr->data;
  snprintf(&charStr_data[0], (size_t)(m + 1), "%s_%d_%d_%d_%d_error.out",
           &b_varargin_1_data[0], nd2, j2, b_loop_ub, c_loop_ub);
  i = b_charStr->size[0] * b_charStr->size[1];
  if (1 > m) {
    b_charStr->size[1] = 0;
  } else {
    b_charStr->size[1] = m;
  }

  emxEnsureCapacity_char_T(b_charStr, i);
  fileid = cfopen(b_charStr, "wb");

  /* 'gkmPWM:107' if fidw == -1 */
  emxFree_char_T(&b_charStr);
  emxFree_char_T(&c_varargin_1);
  emxFree_char_T(&b_varargin_1);
  if (fileid == -1) {
    /* 'gkmPWM:108' fprintf("ERROR: Cannot create training error log.\n"); */
    printf("ERROR: Cannot create training error log.\n");
    fflush(stdout);
    exit(1);
  }

  /* 'gkmPWM:110' for idx=1:length(scorevec) */
  i = dc2->size[1];
  for (m = 0; m < i; m++) {
    /* 'gkmPWM:111' fprintf(fidw, "%f\n", scorevec(idx)); */
    b_NULL = NULL;
    getfilestar(fileid, &filestar, &ipnr);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "%f\n", dc2_data[m]);
      if (ipnr) {
        fflush(filestar);
      }
    }
  }

  emxFree_real_T(&dc2);
}

/* End of code generation (gkmPWM.c) */
