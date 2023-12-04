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
#include "mpower.h"
#include "mtimes.h"
#include "nullAssignment.h"
#include "randperm.h"
#include "relop.h"
#include "repelem.h"
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
  *indvec_loop2, int iii, const emxArray_cell_wrap_14 *ktree2, double k);
static double adjust_PWM(emxArray_real_T *p, const double GC[4]);
static void avg_info(const emxArray_cell_wrap_2 *p, double l_svm,
                     emxArray_real_T *info);
static void b_minus(emxArray_real_T *CT, const emxArray_real_T *ct);
static void b_seed_kmers(const emxArray_char_T *fn, double num, const
  emxArray_cell_wrap_1 *ik, emxArray_cell_wrap_0 *p, emxArray_cell_wrap_1 *mat,
  emxArray_cell_wrap_2 *pwms, double *c);
static void bb_binary_expand_op(emxArray_real_T *kweig, const emxArray_real_T
  *indvec2_loop2, int iii, const emxArray_cell_wrap_14 *ktree2, double k);
static void binary_expand_op(emxArray_real_T *cfile, const emxArray_real_T
  *negvec, double y, double b);
static void cb_binary_expand_op(emxArray_real_T *mat, const emxArray_real_T *x);
static void createMEME(const emxArray_char_T *fileh_Value, const
  emxArray_cell_wrap_2 *PWM, const emxArray_char_T *memefile, double GC, const
  emxArray_real_T *C, double r, const emxArray_real_T *R, double rcorr, const
  emxArray_real_T *E, const emxArray_real_T *Rd);
static void db_binary_expand_op(emxArray_real_T *mat, const emxArray_cell_wrap_2
  *p, int i, int i2, int i3, int i4, int i5);
static void eb_binary_expand_op(emxArray_real_T *kmat, int i, const
  emxArray_real_T *y, const emxArray_real_T *negvec, double b);
static void fb_binary_expand_op(emxArray_real_T *res, const emxArray_real_T *f,
  const emxArray_real_T *C, const emxArray_real_T *ord, int ii);
static void gb_binary_expand_op(emxArray_real_T *kmat, const emxArray_real_T
  *ord, int ii, const emxArray_real_T *f);
static void getEMprob_v3(const emxArray_real_T *PWM, const emxArray_real_T *res,
  const double negmat[16], const emxArray_cell_wrap_14 *poscell, const
  emxArray_real_T *rc, const emxArray_real_T *diffc, const emxArray_real_T *indc,
  const emxArray_real_T *indloc, const emxArray_real_T *xc, double reg, double
  l_svm, double k_svm, double rcnum, double RC, emxArray_real_T *kweig, double
  P[4]);
static void gkmPWM_lagrange(const emxArray_real_T *kweig, const double negmat[16],
  emxArray_cell_wrap_2 *PWM, const emxArray_real_T *negvec, double n, double
  rcorr, double reg, double l_svm, double k_svm, double RC, emxArray_real_T
  *scorevec, emxArray_real_T *C, double *r, emxArray_real_T *R, emxArray_real_T *
  E, emxArray_real_T *Rd);
static void hb_binary_expand_op(emxArray_real_T *kmat, const emxArray_real_T *f,
  int jj, const emxArray_real_T *y, const emxArray_real_T *negvec, double b);
static void i_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3);
static void j_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3);
static void k_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, const emxArray_real_T *rs, int i2, int i3);
static void l_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3);
static void ls_kweigtree(const emxArray_real_T *mat, const double negmat[16],
  const emxArray_cell_wrap_14 *poscell, const emxArray_real_T *c, const
  emxArray_real_T *s, const emxArray_real_T *ind, const emxArray_real_T *indloc,
  const emxArray_real_T *x, double l, double k, double rcnum, emxArray_real_T
  *kweig);
static void ls_kweigtree_norc(const emxArray_real_T *mat, const double negmat[16],
  const emxArray_cell_wrap_14 *poscell, const emxArray_real_T *c, const
  emxArray_real_T *s, const emxArray_real_T *ind, const emxArray_real_T *indloc,
  const emxArray_real_T *x, double l, double k, emxArray_real_T *kweig);
static void m_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3);
static void minus(emxArray_real_T *res, const emxArray_real_T *kweig);
static void n_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, int i1, int i2, const emxArray_real_T *rs, int i3);
static void o_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
  *mat, int j, const emxArray_real_T *rs);
static void plus(emxArray_real_T *b, const emxArray_real_T *res);
static void seed_kmers(const emxArray_char_T *fn, double num,
  emxArray_cell_wrap_0 *p, emxArray_cell_wrap_1 *mat, emxArray_cell_wrap_2 *pwms,
  double *c);
static void u_binary_expand_op(emxArray_creal_T *vec, const emxArray_real_T *b,
  const emxArray_real_T *A, const emxArray_int8_T *x, const creal_T p3_data[]);
static void v_binary_expand_op(creal_T p3_data[], int *p3_size, const creal_T
  ps_data[], const int *ps_size, const creal_T E, const creal_T B_data[], const
  int B_size[2], const creal_T y);
static void w_binary_expand_op(emxArray_creal_T *vec, const emxArray_real_T *b,
  const emxArray_real_T *A, const signed char Posvec_data[], const int
  Posvec_size[2], const creal_T p2_data[]);
static void x_binary_expand_op(creal_T p[4], const double MAT_data[], const int
  MAT_size[2], const creal_T y, const signed char b[4]);
static void y_binary_expand_op(emxArray_real_T *b, const emxArray_real_T *A, int
  iindx);

/* Function Definitions */
static void ab_binary_expand_op(emxArray_real_T *kweig, const emxArray_real_T
  *indvec_loop2, int iii, const emxArray_cell_wrap_14 *ktree2, double k)
{
  const cell_wrap_14 *ktree2_data;
  emxArray_real_T *b_kweig;
  const double *indvec_loop2_data;
  double *b_kweig_data;
  double *kweig_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  ktree2_data = ktree2->data;
  indvec_loop2_data = indvec_loop2->data;
  kweig_data = kweig->data;
  emxInit_real_T(&b_kweig, 1);
  i = b_kweig->size[0];
  if (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] == 1) {
    b_kweig->size[0] = indvec_loop2->size[0];
  } else {
    b_kweig->size[0] = ktree2_data[(int)(k - 1.0) - 1].f1->size[0];
  }

  emxEnsureCapacity_real_T(b_kweig, i);
  b_kweig_data = b_kweig->data;
  stride_0_0 = (indvec_loop2->size[0] != 1);
  stride_1_0 = (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] != 1);
  if (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] == 1) {
    loop_ub = indvec_loop2->size[0];
  } else {
    loop_ub = ktree2_data[(int)(k - 1.0) - 1].f1->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_kweig_data[i] = kweig_data[((int)indvec_loop2_data[i * stride_0_0] +
      kweig->size[0] * (3 - iii)) - 1] + ktree2_data[(int)(k - 1.0) - 1]
      .f1->data[i * stride_1_0];
  }

  loop_ub = b_kweig->size[0];
  for (i = 0; i < loop_ub; i++) {
    kweig_data[((int)indvec_loop2_data[i] + kweig->size[0] * (3 - iii)) - 1] =
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
  bool b;
  bool b2;
  bool guard1 = false;
  bool guard2 = false;
  bool guard3 = false;
  bool guard4 = false;
  p_data = p->data;
  emxInit_real_T(&mat, 2);

  /* extends or truncates the PWM based on information.  I try to do it intelligently, so you may disagree on the condition required for adjustment. */
  /* 'gkmPWM:955' [len,~] = size(p); */
  len = p->size[0];

  /* 'gkmPWM:956' info = zeros(len, 1); */
  /* 'gkmPWM:957' cut = 0.2; */
  /* maximum information of a column to truncate */
  /* 'gkmPWM:958' ext = 0.7; */
  /* minimum information needed to extend */
  /* The rest of the code is easy enough to read through quickly */
  /* 'gkmPWM:960' mat = p+(p==0); */
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

  /* 'gkmPWM:961' vec = 2+sum(mat.*log(mat)/log(2),2); */
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

  /* 'gkmPWM:962' b = true; */
  b = true;

  /* 'gkmPWM:963' b2 = true; */
  b2 = true;

  /* 'gkmPWM:964' while ((vec(1) < cut && max(vec(2:3)) <= ext) || mean(vec(1:3) < cut)) && len > 12 */
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
      if ((double)(((vec_data[0] < 0.2) + (vec_data[1] < 0.2)) + (vec_data[2] <
            0.2)) / 3.0 != 0.0) {
        guard3 = true;
      } else {
        exitg1 = 1;
      }
    }

    if (guard3) {
      if (len > 12.0) {
        /* 'gkmPWM:965' p(1,:) = []; */
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

        /* 'gkmPWM:966' vec(1) = []; */
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

        /* 'gkmPWM:967' len = len-1; */
        len--;

        /* 'gkmPWM:968' b = false; */
        b = false;

        /* 'gkmPWM:969' if ((vec(end) < cut && max(vec(end-2:end-1)) <= ext) || mean(vec(end-2:end) < cut)) && len > 12 */
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

    if (guard2 && ((double)(((vec_data[vec->size[0] - 3] < 0.2) + (vec_data
            [vec->size[0] - 2] < 0.2)) + (vec_data[vec->size[0] - 1] < 0.2)) /
                   3.0 != 0.0)) {
      guard1 = true;
    }

    if (guard1 && (len > 12.0)) {
      /* 'gkmPWM:970' vec(end) = []; */
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

      /* 'gkmPWM:971' p(end,:) = []; */
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

      /* 'gkmPWM:972' len = len-1; */
      len--;

      /* 'gkmPWM:973' b2 = false; */
      b2 = false;
    }
  } while (exitg1 == 0);

  /* 'gkmPWM:976' if b && min(vec(1:2)) > ext && len < 20 */
  if (b) {
    if (vec_data[0] > vec_data[1]) {
      b_vec_data = vec_data[1];
    } else {
      b_vec_data = vec_data[0];
    }

    if ((b_vec_data > 0.7) && (len < 20.0)) {
      /* 'gkmPWM:977' mat = [GC ; mat]; */
      /* 'gkmPWM:978' p = [GC;p]; */
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

      /* 'gkmPWM:979' len = len+1; */
      len++;
    }
  }

  /* 'gkmPWM:981' while ((vec(end) < cut && max(vec(end-2:end-1)) <= ext) || mean(vec(end-2:end) < cut)) && len > 12 */
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
      if ((double)(((vec_data[vec->size[0] - 3] < 0.2) + (vec_data[vec->size[0]
             - 2] < 0.2)) + (vec_data[vec->size[0] - 1] < 0.2)) / 3.0 != 0.0) {
        guard1 = true;
      } else {
        exitg1 = 1;
      }
    }

    if (guard1) {
      if (len > 12.0) {
        /* 'gkmPWM:982' vec(end) = []; */
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

        /* 'gkmPWM:983' p(end,:) = []; */
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

        /* 'gkmPWM:984' len = len-1; */
        len--;

        /* 'gkmPWM:985' b2 = false; */
        b2 = false;
      } else {
        exitg1 = 1;
      }
    }
  } while (exitg1 == 0);

  /* 'gkmPWM:987' if b2 && min(vec(end-1:end)) > ext && len < 20 */
  if (b2) {
    if (vec_data[vec->size[0] - 2] > vec_data[vec->size[0] - 1]) {
      b_vec_data = vec_data[vec->size[0] - 1];
    } else {
      b_vec_data = vec_data[vec->size[0] - 2];
    }

    if ((b_vec_data > 0.7) && (len < 20.0)) {
      /* 'gkmPWM:988' p = [p;GC]; */
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

      mat_data[p->size[0]] = GC[0];
      mat_data[p->size[0] + mat->size[0]] = GC[1];
      mat_data[p->size[0] + mat->size[0] * 2] = GC[2];
      mat_data[p->size[0] + mat->size[0] * 3] = GC[3];
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

      /* 'gkmPWM:989' len = len+1; */
      len++;
    }
  }

  emxFree_real_T(&vec);
  emxFree_real_T(&mat);

  /* 'gkmPWM:991' pp = p; */
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

  /* 'gkmPWM:994' info = zeros(length(p),1); */
  i = info->size[0];
  info->size[0] = p->size[0];
  emxEnsureCapacity_real_T(info, i);
  info_data = info->data;
  nx = p->size[0];
  for (i = 0; i < nx; i++) {
    info_data[i] = 0.0;
  }

  /* 'gkmPWM:995' for i = 1:length(p) */
  i = p->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&x, 2);
  emxInit_real_T(&y, 1);
  y_data = y->data;
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:996' [l,~] = size(p{i}); */
    /* 'gkmPWM:997' mat = p{i}(l_svm:end-l_svm+1,:)+(p{i}(l_svm:end-l_svm+1,:)==0); */
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

    /* 'gkmPWM:998' info(i) = sum(2+sum(mat.*log(max(mat,0))/log(2),2)); */
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

static void b_minus(emxArray_real_T *CT, const emxArray_real_T *ct)
{
  emxArray_real_T *b_CT;
  const double *ct_data;
  double *CT_data;
  double *b_CT_data;
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
  ct_data = ct->data;
  CT_data = CT->data;
  emxInit_real_T(&b_CT, 2);
  i = b_CT->size[0] * b_CT->size[1];
  if (ct->size[0] == 1) {
    b_CT->size[0] = CT->size[0];
  } else {
    b_CT->size[0] = ct->size[0];
  }

  if (ct->size[1] == 1) {
    b_CT->size[1] = CT->size[1];
  } else {
    b_CT->size[1] = ct->size[1];
  }

  emxEnsureCapacity_real_T(b_CT, i);
  b_CT_data = b_CT->data;
  stride_0_0 = (CT->size[0] != 1);
  stride_0_1 = (CT->size[1] != 1);
  stride_1_0 = (ct->size[0] != 1);
  stride_1_1 = (ct->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (ct->size[1] == 1) {
    loop_ub = CT->size[1];
  } else {
    loop_ub = ct->size[1];
  }

  for (i = 0; i < loop_ub; i++) {
    if (ct->size[0] == 1) {
      b_loop_ub = CT->size[0];
    } else {
      b_loop_ub = ct->size[0];
    }

    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_CT_data[i1 + b_CT->size[0] * i] = CT_data[i1 * stride_0_0 + CT->size[0] *
        aux_0_1] - ct_data[i1 * stride_1_0 + ct->size[0] * aux_1_1];
    }

    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }

  i = CT->size[0] * CT->size[1];
  CT->size[0] = b_CT->size[0];
  CT->size[1] = b_CT->size[1];
  emxEnsureCapacity_real_T(CT, i);
  CT_data = CT->data;
  loop_ub = b_CT->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = b_CT->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      CT_data[i1 + CT->size[0] * i] = b_CT_data[i1 + b_CT->size[0] * i];
    }
  }

  emxFree_real_T(&b_CT);
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

  /* 'gkmPWM:208' fid = fopen(fn, 'r'); */
  fileid = cfopen(fn, "rb");

  /*  a = textscan(fid, '%s\t%f\n'); */
  /* 'gkmPWM:211' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'gkmPWM:212' idx=0; */
  idx = 0.0;

  /* 'gkmPWM:213' while ~feof(fid) */
  emxInit_char_T(&b_fileid, 2);
  do {
    exitg1 = 0;
    d = b_feof(fileid);
    if (d == 0.0) {
      /* 'gkmPWM:214' idx=idx+1; */
      idx++;

      /* 'gkmPWM:215' fgetl(fid); */
      b_fgets(fileid, b_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_fileid);
  emxInit_cell_wrap_0(&sequences, 1);

  /* 'gkmPWM:217' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'gkmPWM:218' sequences = cell(idx, 1); */
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

  /* 'gkmPWM:219' sequences = coder.nullcopy(sequences); */
  /* 'gkmPWM:220' alpha = zeros(idx, 1); */
  i = alpha->size[0];
  alpha->size[0] = (int)idx;
  emxEnsureCapacity_real_T(alpha, i);
  alpha_data = alpha->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    alpha_data[i] = 0.0;
  }

  /* 'gkmPWM:221' for cur_idx=1:idx */
  b_d = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&cur_seq, 2);
  emxInit_char_T(&cur_alpha, 2);
  exitg2 = false;
  while ((!exitg2) && (b_d <= (int)idx - 1)) {
    /* 'gkmPWM:222' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'gkmPWM:223' if cur_line == -1 */
    y = (cur_line->size[1] != 0);
    if (y) {
      y = (0 > cur_line->size[1] - 1);
    }

    if (y) {
      exitg2 = true;
    } else {
      /* 'gkmPWM:226' [cur_seq, cur_alpha] = strtok(cur_line, char(9)); */
      b_strtok(cur_line, cur_seq, cur_alpha);

      /* 'gkmPWM:227' alpha(cur_idx,1) = real(str2double(cur_alpha)); */
      dc = str2double(cur_alpha);
      alpha_data[b_d] = dc.re;

      /* 'gkmPWM:228' sequences{cur_idx} = (strip(cur_seq)); */
      strip(cur_seq, sequences_data[b_d].f1);
      b_d++;
    }
  }

  emxFree_char_T(&cur_alpha);
  emxFree_char_T(&cur_seq);
  emxFree_char_T(&cur_line);
  emxInit_real_T(&M, 1);
  emxInit_int32_T(&D, 1);

  /* 'gkmPWM:230' fclose(fid); */
  cfclose(fileid);

  /*  [w, ind] = sort(a{2}, pn); */
  /*  s = a{1}(ind(1:min([100000 length(a{1})]))); */
  /* 'gkmPWM:235' [w, ind] = sort(alpha, pn); */
  b_sort(alpha, D);
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

  /* 'gkmPWM:236' s_len = min([100000 length(sequences)]); */
  m[0] = 100000.0;
  m[1] = sequences->size[0];

  /* 'gkmPWM:237' s = cell(s_len, 1); */
  unnamed_idx_0_tmp_tmp = (int)b_minimum(m);
  i = s->size[0];
  s->size[0] = unnamed_idx_0_tmp_tmp;
  emxEnsureCapacity_cell_wrap_0(s, i);
  s_data = s->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    s_data[i].f1->size[0] = 1;
    s_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:238' s = coder.nullcopy(s); */
  /* 'gkmPWM:239' for cur_idx=1:s_len */
  for (b_d = 0; b_d < unnamed_idx_0_tmp_tmp; b_d++) {
    /* 'gkmPWM:240' s{cur_idx} = sequences{ind(cur_idx)}; */
    i = s_data[b_d].f1->size[0] * s_data[b_d].f1->size[1];
    s_data[b_d].f1->size[0] = 1;
    s_data[b_d].f1->size[1] = sequences_data[(int)M_data[b_d] - 1].f1->size[1];
    emxEnsureCapacity_char_T(s_data[b_d].f1, i);
    loop_ub = sequences_data[(int)M_data[b_d] - 1].f1->size[1];
    for (i = 0; i < loop_ub; i++) {
      s_data[b_d].f1->data[i] = sequences_data[(int)M_data[b_d] - 1].f1->data[i];
    }
  }

  /* 'gkmPWM:244' l = length(s{1}); */
  varargin_2 = s_data[0].f1->size[1] - 1;
  l = s_data[0].f1->size[1] - 4;

  /* 'gkmPWM:245' k = round(l/2)+1; */
  x = (int)rt_roundd((double)s_data[0].f1->size[1] / 2.0);

  /* 'gkmPWM:246' ikl = length(ik); */
  b_varargin_2 = ik->size[0];

  /* 'gkmPWM:247' p = cell(num,1); */
  unnamed_idx_0_tmp_tmp = (int)num;
  i = p->size[0];
  p->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_0(p, i);
  p_data = p->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    p_data[i].f1->size[0] = 1;
    p_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:248' p = coder.nullcopy(p); */
  i = sequences->size[0];
  sequences->size[0] = p->size[0];
  emxEnsureCapacity_cell_wrap_0(sequences, i);
  sequences_data = sequences->data;

  /* 'gkmPWM:249' c = ikl+1; */
  *c = (double)ik->size[0] + 1.0;

  /* 'gkmPWM:250' p{1} = s{1}; */
  i = sequences_data[0].f1->size[0] * sequences_data[0].f1->size[1];
  sequences_data[0].f1->size[0] = 1;
  sequences_data[0].f1->size[1] = s_data[0].f1->size[1];
  emxEnsureCapacity_char_T(sequences_data[0].f1, i);
  loop_ub = s_data[0].f1->size[1];
  for (i = 0; i < loop_ub; i++) {
    sequences_data[0].f1->data[i] = s_data[0].f1->data[i];
  }

  emxInit_cell_wrap_1(&b_mat);

  /* 'gkmPWM:251' mat = cell(ikl+num,1); */
  b_d = (int)((double)ik->size[0] + num);
  i = b_mat->size[0];
  b_mat->size[0] = (int)((double)ik->size[0] + num);
  emxEnsureCapacity_cell_wrap_1(b_mat, i);
  mat_data = b_mat->data;
  for (i = 0; i < b_d; i++) {
    mat_data[i].f1->size[0] = 1;
    mat_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:252' mat = coder.nullcopy(mat); */
  /*  mat(1:ikl) = ik; */
  /* 'gkmPWM:254' for cur_idx=1:length(ik) */
  i = ik->size[0];
  for (b_d = 0; b_d < i; b_d++) {
    /* 'gkmPWM:255' mat{cur_idx} = ik{cur_idx}; */
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

  /* 'gkmPWM:257' mat{c} = letterconvert(s{1}); */
  letterconvert(s_data[0].f1, mat_data[ik->size[0]].f1);

  /* 'gkmPWM:258' pwms = cell(num,1); */
  i = b_pwms->size[0];
  b_pwms->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_2(b_pwms, i);
  pwms_data = b_pwms->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    pwms_data[i].f1->size[0] = 0;
    pwms_data[i].f1->size[1] = 4;
  }

  /* 'gkmPWM:259' pwms = coder.nullcopy(pwms); */
  /* 'gkmPWM:260' for i = 1:num */
  for (b_i = 0; b_i < unnamed_idx_0_tmp_tmp; b_i++) {
    /* 'gkmPWM:261' pwms{i} = zeros(l,4); */
    i = pwms_data[b_i].f1->size[0] * pwms_data[b_i].f1->size[1];
    pwms_data[b_i].f1->size[0] = varargin_2 + 1;
    pwms_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(pwms_data[b_i].f1, i);
    loop_ub = (varargin_2 + 1) << 2;
    for (i = 0; i < loop_ub; i++) {
      pwms_data[b_i].f1->data[i] = 0.0;
    }
  }

  /* 'gkmPWM:263' for i = 1:l */
  i = s_data[0].f1->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:264' pwms{1}(i,mat{c}(i)+1) = pwms{1}(i,mat{c}(i)+1)+w(i); */
    d = mat_data[b_varargin_2].f1->data[b_i];
    pwms_data[0].f1->data[b_i + pwms_data[0].f1->size[0] * ((int)(d + 1.0) - 1)]
      += alpha_data[b_i];
  }

  /* 'gkmPWM:266' B = zeros(9,1); */
  /* 'gkmPWM:267' BB = zeros(9,1); */
  for (b_i = 0; b_i < 9; b_i++) {
    B[b_i] = 0;
    BB[b_i] = 0;
  }

  /* 'gkmPWM:268' B(1:5) = (0:4)'; */
  /* 'gkmPWM:269' B(6:9) = 0; */
  B[5] = 0;
  B[6] = 0;
  B[7] = 0;
  B[8] = 0;

  /* 'gkmPWM:270' BB(1:5) = 0; */
  for (b_i = 0; b_i < 5; b_i++) {
    B[b_i] = (signed char)b_i;
    BB[b_i] = 0;
  }

  /* 'gkmPWM:271' BB(6:9) = (1:4)'; */
  BB[5] = 1;
  BB[6] = 2;
  BB[7] = 3;
  BB[8] = 4;

  /* 'gkmPWM:272' CC = [l l-1 l-2 l-3 l-4 l-1 l-2 l-3 l-4]; */
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
  /* 'gkmPWM:274' for i = 2:100000 */
  b_i = 1;
  emxInit_real_T(&ss, 2);
  emxInit_real_T(&rs, 2);
  emxInit_real_T(&DD, 1);
  emxInit_boolean_T(&b_x, 2);
  exitg2 = false;
  while ((!exitg2) && (b_i - 1 < 99999)) {
    /* 'gkmPWM:275' ss = letterconvert(s{i}); */
    letterconvert(s_data[b_i].f1, ss);
    ss_data = ss->data;

    /* 'gkmPWM:276' rs = 3-fliplr(ss); */
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

    /* 'gkmPWM:277' M = zeros(c,1); */
    loop_ub_tmp = (int)*c;
    i = M->size[0];
    M->size[0] = (int)*c;
    emxEnsureCapacity_real_T(M, i);
    M_data = M->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      M_data[i] = 0.0;
    }

    /* 'gkmPWM:278' D = zeros(c,1); */
    i = D->size[0];
    D->size[0] = (int)*c;
    emxEnsureCapacity_int32_T(D, i);
    D_data = D->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      D_data[i] = 0;
    }

    /* 'gkmPWM:279' DD = zeros(c,1); */
    i = DD->size[0];
    DD->size[0] = (int)*c;
    emxEnsureCapacity_real_T(DD, i);
    DD_data = DD->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      DD_data[i] = 0.0;
    }

    /* 'gkmPWM:280' for j = 1:c */
    for (j = 0; j < loop_ub_tmp; j++) {
      /* 'gkmPWM:281' [m,d] = max([sum(mat{j}==ss) sum(mat{j}(2:end)==ss(1:l-1)) sum(mat{j}(3:end)==ss(1:l-2)) sum(mat{j}(4:end)==ss(1:l-3)) sum(mat{j}(5:end)==ss(1:l-4)) sum(mat{j}(1:l-1)==ss(2:end)) sum(mat{j}(1:l-2)==ss(3:end)) sum(mat{j}(1:l-3)==ss(4:end)) sum(mat{j}(1:l-4)==ss(5:end))]); */
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
        o_binary_expand_op(b_x, b_mat, j, ss);
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
        n_binary_expand_op(b_x, b_mat, j, i, b_y - 1, ss, unnamed_idx_0_tmp_tmp
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
        m_binary_expand_op(b_x, b_mat, j, b_d, c_y - 1, ss, d_y - 1);
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
        l_binary_expand_op(b_x, b_mat, j, e_y, i1 - 1, ss, i2);
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
        l_binary_expand_op(b_x, b_mat, j, i3, i4 - 1, ss, i5 - 1);
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
        k_binary_expand_op(b_x, b_mat, j, loop_ub - 2, ss, i6, i7 - 1);
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
        j_binary_expand_op(b_x, b_mat, j, b_loop_ub - 1, ss, i8, i9 - 1);
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
        i_binary_expand_op(b_x, b_mat, j, c_loop_ub, ss, i10, i11 - 1);
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
        i_binary_expand_op(b_x, b_mat, j, d_loop_ub - 1, ss, i12, i13 - 1);
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

      /* 'gkmPWM:282' [mm,dd] = max([sum(mat{j}==rs) sum(mat{j}(2:end)==rs(1:l-1)) sum(mat{j}(3:end)==rs(1:l-2)) sum(mat{j}(4:end)==rs(1:l-3)) sum(mat{j}(5:end)==rs(1:l-4)) sum(mat{j}(1:l-1)==rs(2:end)) sum(mat{j}(1:l-2)==rs(3:end)) sum(mat{j}(1:l-3)==rs(4:end)) sum(mat{j}(1:l-4)==rs(5:end))]); */
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
        o_binary_expand_op(b_x, b_mat, j, rs);
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
        n_binary_expand_op(b_x, b_mat, j, i, b_y - 1, rs, unnamed_idx_0_tmp_tmp
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
        m_binary_expand_op(b_x, b_mat, j, b_d, c_y - 1, rs, d_y - 1);
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
        l_binary_expand_op(b_x, b_mat, j, e_y, i1 - 1, rs, i2);
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
        l_binary_expand_op(b_x, b_mat, j, i3, i4 - 1, rs, i5 - 1);
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
        k_binary_expand_op(b_x, b_mat, j, loop_ub - 2, rs, i6, i7 - 1);
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
        j_binary_expand_op(b_x, b_mat, j, b_loop_ub - 1, rs, i8, i9 - 1);
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
        i_binary_expand_op(b_x, b_mat, j, c_loop_ub, rs, i10, i11 - 1);
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
        i_binary_expand_op(b_x, b_mat, j, d_loop_ub - 1, rs, i12, i13 - 1);
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

      /* 'gkmPWM:283' [M(j),ddd] = max([m mm]); */
      m[0] = curr_pos;
      m[1] = idx;
      c_maximum(m, &M_data[j], &unnamed_idx_0_tmp_tmp);

      /* 'gkmPWM:284' if ddd == 1 */
      if (unnamed_idx_0_tmp_tmp == 1) {
        /* 'gkmPWM:285' D(j) = d; */
        D_data[j] = iindx;

        /* 'gkmPWM:286' DD(j) = 1; */
        DD_data[j] = 1.0;
      } else {
        /* 'gkmPWM:287' else */
        /* 'gkmPWM:288' D(j) = dd; */
        D_data[j] = b_d;

        /* 'gkmPWM:289' DD(j) = 2; */
        DD_data[j] = 2.0;
      }
    }

    /* 'gkmPWM:292' if max(M) < k */
    if (maximum(M) < (double)x + 1.0) {
      /* 'gkmPWM:293' c = c+1; */
      (*c)++;

      /* 'gkmPWM:294' p{c-ikl} = s{i}; */
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

      /* 'gkmPWM:295' mat{c} = ss; */
      i = mat_data[(int)*c - 1].f1->size[0] * mat_data[(int)*c - 1].f1->size[1];
      mat_data[(int)*c - 1].f1->size[0] = 1;
      mat_data[(int)*c - 1].f1->size[1] = ss->size[1];
      emxEnsureCapacity_real_T(mat_data[(int)*c - 1].f1, i);
      loop_ub = ss->size[1];
      for (i = 0; i < loop_ub; i++) {
        mat_data[(int)*c - 1].f1->data[i] = ss_data[i];
      }

      /* 'gkmPWM:296' ss = ss+1; */
      i = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity_real_T(ss, i);
      ss_data = ss->data;
      loop_ub = ss->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        ss_data[i]++;
      }

      /* 'gkmPWM:297' for j = 1:l */
      for (j = 0; j <= varargin_2; j++) {
        /* 'gkmPWM:298' pwms{c-ikl}(j,ss(j)) = pwms{c-ikl}(j,ss(j))+w(i); */
        i = (int)ss_data[j] - 1;
        pwms_data[(int)(*c - (double)b_varargin_2) - 1].f1->data[j + pwms_data
          [(int)(*c - (double)b_varargin_2) - 1].f1->size[0] * i] = pwms_data
          [(int)(*c - (double)ik->size[0]) - 1].f1->data[j + pwms_data[(int)(*c
          - (double)ik->size[0]) - 1].f1->size[0] * i] + alpha_data[b_i];
      }
    } else {
      /* 'gkmPWM:300' else */
      /* 'gkmPWM:301' [~,d] = max(M); */
      d_maximum(M, &curr_pos, &iindx);

      /* 'gkmPWM:302' if DD(d) == 1 && d > ikl */
      d = DD_data[iindx - 1];
      if ((d == 1.0) && (iindx > b_varargin_2)) {
        /* 'gkmPWM:303' ss = ss+1; */
        i = ss->size[0] * ss->size[1];
        ss->size[0] = 1;
        emxEnsureCapacity_real_T(ss, i);
        ss_data = ss->data;
        loop_ub = ss->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          ss_data[i]++;
        }

        /* 'gkmPWM:304' d = d-ikl; */
        b_d = (iindx - b_varargin_2) - 1;

        /* 'gkmPWM:305' for j = 1:CC(D(d)) */
        i = CC[D_data[b_d] - 1];
        for (j = 0; j < i; j++) {
          /* 'gkmPWM:306' pwms{d}(j+B(D(d)),ss(j+BB(D(d)))) = pwms{d}(j+B(D(d)),ss(j+BB(D(d))))+w(i); */
          b_y = (int)((unsigned int)j + B[D_data[b_d] - 1]);
          unnamed_idx_0_tmp_tmp = (int)ss_data[(int)((unsigned int)j +
            BB[D_data[b_d] - 1])] - 1;
          pwms_data[b_d].f1->data[b_y + pwms_data[b_d].f1->size[0] *
            unnamed_idx_0_tmp_tmp] += alpha_data[b_i];
        }
      } else if ((d == 2.0) && (iindx > b_varargin_2)) {
        /* 'gkmPWM:308' elseif DD(d) == 2 && d > ikl */
        /* 'gkmPWM:309' rs = rs+1; */
        i = rs->size[0] * rs->size[1];
        rs->size[0] = 1;
        emxEnsureCapacity_real_T(rs, i);
        rs_data = rs->data;
        loop_ub = rs->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          rs_data[i]++;
        }

        /* 'gkmPWM:310' d = d-ikl; */
        b_d = (iindx - b_varargin_2) - 1;

        /* 'gkmPWM:311' for j = 1:CC(D(d)) */
        i = CC[D_data[b_d] - 1];
        for (j = 0; j < i; j++) {
          /* 'gkmPWM:312' pwms{d}(j+B(D(d)),rs(j+BB(D(d)))) = pwms{d}(j+B(D(d)),rs(j+BB(D(d))))+w(i); */
          b_y = (int)((unsigned int)j + B[D_data[b_d] - 1]);
          unnamed_idx_0_tmp_tmp = (int)rs_data[(int)((unsigned int)j +
            BB[D_data[b_d] - 1])] - 1;
          pwms_data[b_d].f1->data[b_y + pwms_data[b_d].f1->size[0] *
            unnamed_idx_0_tmp_tmp] += alpha_data[b_i];
        }
      }
    }

    /* 'gkmPWM:316' if c == num+ikl */
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
  /* 'gkmPWM:321' new_mat = cell(c, 1); */
  /* 'gkmPWM:322' for cur_idx=1:c */
  i = (int)*c;
  b_y = mat->size[0];
  mat->size[0] = (int)*c;
  emxEnsureCapacity_cell_wrap_1(mat, b_y);
  b_mat_data = mat->data;
  for (b_d = 0; b_d < i; b_d++) {
    /* 'gkmPWM:323' new_mat{cur_idx} = mat{cur_idx}; */
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

  /* 'gkmPWM:325' mat = new_mat; */
  /*  p = p(1:c-ikl); */
  /* 'gkmPWM:328' p_len = c-ikl; */
  /* 'gkmPWM:329' new_p = cell(p_len, 1); */
  /* 'gkmPWM:330' for cur_idx=1:p_len */
  i = (int)(*c - (double)ik->size[0]);
  b_y = p->size[0];
  p->size[0] = (int)(*c - (double)ik->size[0]);
  emxEnsureCapacity_cell_wrap_0(p, b_y);
  p_data = p->data;
  for (b_d = 0; b_d < i; b_d++) {
    /* 'gkmPWM:331' new_p{cur_idx} = p{cur_idx}; */
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

  /* 'gkmPWM:333' p = new_p; */
  /*  pwms = pwms(1:c-ikl); */
  /* 'gkmPWM:336' pwms_len = c-ikl; */
  /* 'gkmPWM:337' new_pwms = cell(pwms_len, 1); */
  /* 'gkmPWM:338' for cur_idx=1:pwms_len */
  i = (int)(*c - (double)ik->size[0]);
  b_y = pwms->size[0];
  pwms->size[0] = (int)(*c - (double)ik->size[0]);
  emxEnsureCapacity_cell_wrap_2(pwms, b_y);
  b_pwms_data = pwms->data;
  for (b_d = 0; b_d < i; b_d++) {
    /* 'gkmPWM:339' new_pwms{cur_idx} = pwms{cur_idx}; */
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

  /* 'gkmPWM:341' pwms = new_pwms; */
  /* 'gkmPWM:343' for i = 1:c-ikl */
  i = (int)(*c - (double)ik->size[0]);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:344' for j = 1:l */
    for (j = 0; j <= varargin_2; j++) {
      /* 'gkmPWM:345' pwms{i}(j,:) = pwms{i}(j,:)/sum(pwms{i}(j,:)); */
      curr_pos = ((b_pwms_data[b_i].f1->data[j] + b_pwms_data[b_i].f1->data[j +
                   b_pwms_data[b_i].f1->size[0]]) + b_pwms_data[b_i].f1->data[j
                  + b_pwms_data[b_i].f1->size[0] * 2]) + b_pwms_data[b_i]
        .f1->data[j + b_pwms_data[b_i].f1->size[0] * 3];
      b_pwms_data[b_i].f1->data[j] /= curr_pos;
      b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0]] /= curr_pos;
      b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0] * 2] /=
        curr_pos;
      b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0] * 3] /=
        curr_pos;
    }
  }
}

static void bb_binary_expand_op(emxArray_real_T *kweig, const emxArray_real_T
  *indvec2_loop2, int iii, const emxArray_cell_wrap_14 *ktree2, double k)
{
  const cell_wrap_14 *ktree2_data;
  emxArray_real_T *b_kweig;
  const double *indvec2_loop2_data;
  double *b_kweig_data;
  double *kweig_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  ktree2_data = ktree2->data;
  indvec2_loop2_data = indvec2_loop2->data;
  kweig_data = kweig->data;
  emxInit_real_T(&b_kweig, 1);
  i = b_kweig->size[0];
  if (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] == 1) {
    b_kweig->size[0] = indvec2_loop2->size[0];
  } else {
    b_kweig->size[0] = ktree2_data[(int)(k - 1.0) - 1].f1->size[0];
  }

  emxEnsureCapacity_real_T(b_kweig, i);
  b_kweig_data = b_kweig->data;
  stride_0_0 = (indvec2_loop2->size[0] != 1);
  stride_1_0 = (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] != 1);
  if (ktree2_data[(int)(k - 1.0) - 1].f1->size[0] == 1) {
    loop_ub = indvec2_loop2->size[0];
  } else {
    loop_ub = ktree2_data[(int)(k - 1.0) - 1].f1->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    b_kweig_data[i] = kweig_data[((int)indvec2_loop2_data[i * stride_0_0] +
      kweig->size[0] * iii) - 1] + ktree2_data[(int)(k - 1.0) - 1].f1->data[i *
      stride_1_0];
  }

  loop_ub = b_kweig->size[0];
  for (i = 0; i < loop_ub; i++) {
    kweig_data[((int)indvec2_loop2_data[i] + kweig->size[0] * iii) - 1] =
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
  cell_wrap_3 *cur_PWM_data;
  cell_wrap_3 *p_data;
  emxArray_cell_wrap_0 *names;
  emxArray_cell_wrap_3 *b_cur_PWM;
  emxArray_cell_wrap_3 *cur_PWM;
  emxArray_cell_wrap_3 *p;
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
  double a;
  double b_MM;
  double b_varargin_1;
  double d;
  double *mat_data;
  double *rmat_data;
  double *varargin_1_data;
  int MM;
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

  /* 'gkmPWM:133' num = numel(C); */
  /* 'gkmPWM:134' GC = round(GC*100)/100; */
  GC = rt_roundd(GC * 100.0) / 100.0;

  /* 'gkmPWM:135' GCvec = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  /* 'gkmPWM:137' f = strfind(fileread(memefile),'MOTIF'); */
  fileread(memefile, b_memefile);

  /* 'gkmPWM:138' num2 = length(strfind(fileread(memefile),'MOTIF')); */
  fileread(memefile, text);
  text_data = text->data;
  emxFree_char_T(&b_memefile);
  if (text->size[1] == 0) {
    match_idx = 0;
  } else {
    emxInit_int32_T(&match_out, 2);
    emxInit_int32_T(&matches, 2);
    i = matches->size[0] * matches->size[1];
    matches->size[0] = 1;
    matches->size[1] = text->size[1];
    emxEnsureCapacity_int32_T(matches, i);
    matches_data = matches->data;
    match_idx = 0;
    i = text->size[1];
    for (b_i = 0; b_i <= i - 5; b_i++) {
      j = 1;
      while ((j <= 5) && (text_data[(b_i + j) - 1] == cv[j - 1])) {
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

  /* 'gkmPWM:139' [p,names] = getmotif(memefile,1:num2); */
  emxInit_real_T(&y, 2);
  if (match_idx < 1) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = match_idx;
    emxEnsureCapacity_real_T(y, i);
    mat_data = y->data;
    loop_ub = match_idx - 1;
    for (i = 0; i <= loop_ub; i++) {
      mat_data[i] = (double)i + 1.0;
    }
  }

  emxInit_uint32_T(&lenvec);
  emxInit_cell_wrap_3(&p);
  emxInit_cell_wrap_0(&names, 1);
  getmotif(memefile, y, p, names);
  names_data = names->data;
  p_data = p->data;

  /* 'gkmPWM:141' fid = fopen(sprintf('%s_denovo.meme', fileh), 'w'); */
  b_sprintf(fileh_Value, text);
  fileid = cfopen(text, "wb");

  /* 'gkmPWM:142' fid2 = fopen(sprintf('%s_gkmPWM.out', fileh), 'w'); */
  c_sprintf(fileh_Value, text);
  b_fileid = cfopen(text, "wb");

  /* 'gkmPWM:143' lenvec = zeros(num2,1); */
  i = lenvec->size[0];
  lenvec->size[0] = match_idx;
  emxEnsureCapacity_uint32_T(lenvec, i);
  lenvec_data = lenvec->data;

  /* 'gkmPWM:144' for i = 1:num2 */
  for (b_i = 0; b_i < match_idx; b_i++) {
    /* 'gkmPWM:145' [lenvec(i),~] = size(p{i}); */
    lenvec_data[b_i] = (unsigned int)p_data[b_i].f1->size[0];
  }

  /*  mf = fopen('motif_logos/short_name', 'r'); */
  /*  % names = textscan(mf,'%*d\t%s\t%*s\n'); */
  /*   */
  /*  curr_pos = ftell(mf); */
  /*  idx=0; */
  /*  while ~feof(mf) */
  /*      idx=idx+1; */
  /*      fgetl(mf); */
  /*  end */
  /*  fseek(mf, curr_pos, 'bof'); */
  /*  names = cell(idx, 1); */
  /*  names = coder.nullcopy(names); */
  /*  for cur_idx=1:idx */
  /*      cur_line = fgetl(mf); */
  /*      if cur_line == -1 */
  /*          break */
  /*      end */
  /*      [~, cur_names] = strtok(cur_line, char(9)); */
  /*      [cur_short_name, ~] = strtok(cur_names, char(9)); */
  /*      names{cur_idx} = string(strip(cur_short_name)); */
  /*  end */
  /*  fclose(mf); */
  /*  names = names{1}; */
  /* 'gkmPWM:171' fprintf(fid, 'MEME\n\n'); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "MEME\n\n");
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:172' fprintf(fid, 'ALPHABET= ACGT\n\n'); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "ALPHABET= ACGT\n\n");
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:173' fprintf(fid, 'Correlation with SVM weight vector: %0.3f\n\n', r); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Correlation with SVM weight vector: %0.3f\n\n", r);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:174' fprintf(fid, 'Max PWM Correlation: %0.3f\n\n', rcorr); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Max PWM Correlation: %0.3f\n\n", rcorr);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:175' fprintf(fid, 'Background letter frequencies (from negative set)\n'); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Background letter frequencies (from negative set)\n");
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:176' fprintf(fid, 'A %0.2f C %0.2f G %0.2f T %0.2f\n\n', GCvec(1), GCvec(2), GCvec(3), GCvec(4)); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "A %0.2f C %0.2f G %0.2f T %0.2f\n\n", 0.5 - GC / 2.0, GC /
            2.0, GC / 2.0, 0.5 - GC / 2.0);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:177' fprintf(fid2, 'Correlation with SVM weight vector:\t%0.3f\n', r); */
  b_NULL = NULL;
  getfilestar(b_fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Correlation with SVM weight vector:\t%0.3f\n", r);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:178' fprintf(fid2, 'Max PWM Correlation:\t%0.3f\n', rcorr); */
  b_NULL = NULL;
  getfilestar(b_fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Max PWM Correlation:\t%0.3f\n", rcorr);
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:179' fprintf(fid2, 'MOTIF\tID\tSimilarity\tRedundancy\tWeight\tZ\tError\n'); */
  b_NULL = NULL;
  getfilestar(b_fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "MOTIF\tID\tSimilarity\tRedundancy\tWeight\tZ\tError\n");
    if (autoflush) {
      fflush(filestar);
    }
  }

  /* 'gkmPWM:180' a = 1; */
  a = 1.0;

  /* 'gkmPWM:181' b = 1; */
  /* 'gkmPWM:182' for i = 1:num */
  i = C->size[0];
  emxInit_cell_wrap_3(&cur_PWM);
  emxInit_cell_wrap_3(&b_cur_PWM);
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&rmat, 2);
  emxInit_real_T(&varargin_1, 2);
  emxInit_real_T(&A, 2);
  emxInit_real_T(&b_r, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:183' [len,~] = size(PWM{i}); */
    /* 'gkmPWM:184' p_len = length(p); */
    /* 'gkmPWM:185' cur_PWM = cell(p_len+1,1); */
    match_idx = p->size[0] + 1;
    i1 = cur_PWM->size[0];
    cur_PWM->size[0] = p->size[0] + 1;
    emxEnsureCapacity_cell_wrap_3(cur_PWM, i1);
    cur_PWM_data = cur_PWM->data;
    for (i1 = 0; i1 < match_idx; i1++) {
      cur_PWM_data[i1].f1->size[0] = 0;
      cur_PWM_data[i1].f1->size[1] = 0;
    }

    /* 'gkmPWM:186' cur_PWM = coder.nullcopy(cur_PWM); */
    i1 = b_cur_PWM->size[0];
    b_cur_PWM->size[0] = cur_PWM->size[0];
    emxEnsureCapacity_cell_wrap_3(b_cur_PWM, i1);
    cur_PWM_data = b_cur_PWM->data;

    /* 'gkmPWM:187' cur_PWM{1} = PWM{i}; */
    i1 = cur_PWM_data[0].f1->size[0] * cur_PWM_data[0].f1->size[1];
    cur_PWM_data[0].f1->size[0] = PWM_data[b_i].f1->size[0];
    cur_PWM_data[0].f1->size[1] = 4;
    emxEnsureCapacity_real_T(cur_PWM_data[0].f1, i1);
    loop_ub = PWM_data[b_i].f1->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      cur_PWM_data[0].f1->data[i1] = PWM_data[b_i].f1->data[i1];
    }

    /* 'gkmPWM:188' for cur_idx=1:p_len */
    i1 = p->size[0];
    for (match_idx = 0; match_idx < i1; match_idx++) {
      /* 'gkmPWM:189' cur_PWM{cur_idx+1} = p{cur_idx}; */
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

    /* 'gkmPWM:191' [M, MM] = matchMotif(cur_PWM, [len;lenvec]); */
    /* 'gkmPWM:367' n = length(lenvec)-1; */
    /* 'gkmPWM:368' simmat = ones(n-1,1); */
    /* 'gkmPWM:369' for i = 1:n+1 */
    i1 = lenvec->size[0] + 1;
    for (match_idx = 0; match_idx < i1; match_idx++) {
      /* 'gkmPWM:370' mot{i} = mot{i}-1/4; */
      loop_ub = cur_PWM_data[match_idx].f1->size[0] * cur_PWM_data[match_idx].
        f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        cur_PWM_data[match_idx].f1->data[i2] -= 0.25;
      }

      /* 'gkmPWM:371' mot{i} = mot{i}/sqrt(sum(sum(mot{i}.^2))); */
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

      f_sum(A, b_r);
      b_varargin_1 = sqrt(sum(b_r));
      loop_ub = cur_PWM_data[match_idx].f1->size[0] * cur_PWM_data[match_idx].
        f1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        cur_PWM_data[match_idx].f1->data[i2] /= b_varargin_1;
      }
    }

    /* 'gkmPWM:373' M = 0; */
    b_varargin_1 = 0.0;

    /* 'gkmPWM:374' ind = 1; */
    MM = 0;

    /* 'gkmPWM:375' for j = 2:n+1 */
    i1 = lenvec->size[0] + 1;
    for (j = 0; j <= i1 - 2; j++) {
      /* 'gkmPWM:376' mat = mot{1}*mot{j}'; */
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

      /* 'gkmPWM:377' rmat = rot90(mot{1},2)*mot{j}'; */
      d_rot90(cur_PWM_data[0].f1, A);
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

      /* 'gkmPWM:378' MM = max([sum(spdiags(mat)) sum(spdiags(rmat))]); */
      spdiags(mat, A);
      f_sum(A, y);
      mat_data = y->data;
      spdiags(rmat, A);
      f_sum(A, b_r);
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

      match_idx = varargin_1->size[1];
      if (varargin_1->size[1] <= 2) {
        if (varargin_1->size[1] == 1) {
          b_MM = varargin_1_data[0];
        } else if (varargin_1_data[0] < varargin_1_data[varargin_1->size[1] - 1])
        {
          b_MM = varargin_1_data[varargin_1->size[1] - 1];
        } else {
          b_MM = varargin_1_data[0];
        }
      } else {
        b_MM = varargin_1_data[0];
        for (loop_ub = 2; loop_ub <= match_idx; loop_ub++) {
          d = varargin_1_data[loop_ub - 1];
          if (b_MM < d) {
            b_MM = d;
          }
        }
      }

      /* 'gkmPWM:379' if MM > M */
      if (b_MM > b_varargin_1) {
        /* 'gkmPWM:380' M = MM; */
        b_varargin_1 = b_MM;

        /* 'gkmPWM:381' ind = j-1; */
        MM = j;
      }
    }

    /*  matches motifs to the best motif in our database */
    /*  [M, MM] = matchMotif([PWM{i}; p], [len;lenvec]);% matches motifs to the best motif in our database */
    /* 'gkmPWM:195' fprintf(fid, 'MOTIF %d %s\n', int32(a), names{MM}); */
    i1 = text->size[0] * text->size[1];
    text->size[0] = 1;
    text->size[1] = names_data[MM].f1->size[1] + 1;
    emxEnsureCapacity_char_T(text, i1);
    text_data = text->data;
    loop_ub = names_data[MM].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      text_data[i1] = names_data[MM].f1->data[i1];
    }

    text_data[names_data[MM].f1->size[1]] = '\x00';
    b_NULL = NULL;
    getfilestar(fileid, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "MOTIF %d %s\n", (int)a, &text_data[0]);
      if (autoflush) {
        fflush(filestar);
      }
    }

    /* 'gkmPWM:196' fprintf(fid2,'%s\t%d\t%0.3f\t%0.3f\t%0.2f\t%0.3f\t%0.3f\n',names{MM},int32(MM),M,Rd(i),C(i),R(i),E(i)); */
    i1 = text->size[0] * text->size[1];
    text->size[0] = 1;
    text->size[1] = names_data[MM].f1->size[1] + 1;
    emxEnsureCapacity_char_T(text, i1);
    text_data = text->data;
    loop_ub = names_data[MM].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      text_data[i1] = names_data[MM].f1->data[i1];
    }

    text_data[names_data[MM].f1->size[1]] = '\x00';
    b_NULL = NULL;
    getfilestar(b_fileid, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "%s\t%d\t%0.3f\t%0.3f\t%0.2f\t%0.3f\t%0.3f\n",
              &text_data[0], MM + 1, b_varargin_1, Rd_data[b_i], C_data[b_i],
              R_data[b_i], E_data[b_i]);
      if (autoflush) {
        fflush(filestar);
      }
    }

    /* 'gkmPWM:197' fprintf(fid, 'weight= %0.3f l= 4 w= %d z-score= %0.2f motifsim= %0.3f\n', C(i), int32(len), R(i), M); */
    b_NULL = NULL;
    getfilestar(fileid, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar,
              "weight= %0.3f l= 4 w= %d z-score= %0.2f motifsim= %0.3f\n",
              C_data[b_i], PWM_data[b_i].f1->size[0], R_data[b_i], b_varargin_1);
      if (autoflush) {
        fflush(filestar);
      }
    }

    /* 'gkmPWM:198' for j = 1:len */
    i1 = PWM_data[b_i].f1->size[0];
    for (j = 0; j < i1; j++) {
      /* 'gkmPWM:199' fprintf(fid, '%0.3f %0.3f %0.3f %0.3f\n',PWM{i}(j,1),PWM{i}(j,2),PWM{i}(j,3),PWM{i}(j,4)); */
      b_NULL = NULL;
      getfilestar(fileid, &filestar, &autoflush);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "%0.3f %0.3f %0.3f %0.3f\n", PWM_data[b_i].f1->data[j],
                PWM_data[b_i].f1->data[j + PWM_data[b_i].f1->size[0]],
                PWM_data[b_i].f1->data[j + PWM_data[b_i].f1->size[0] * 2],
                PWM_data[b_i].f1->data[j + PWM_data[b_i].f1->size[0] * 3]);
        if (autoflush) {
          fflush(filestar);
        }
      }
    }

    /* 'gkmPWM:201' fprintf(fid, '\n'); */
    b_NULL = NULL;
    getfilestar(fileid, &filestar, &autoflush);
    if (!(filestar == b_NULL)) {
      fprintf(filestar, "\n");
      if (autoflush) {
        fflush(filestar);
      }
    }

    /* 'gkmPWM:202' a = a+1; */
    a++;
  }

  emxFree_real_T(&b_r);
  emxFree_real_T(&A);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&rmat);
  emxFree_real_T(&mat);
  emxFree_real_T(&y);
  emxFree_char_T(&text);
  emxFree_cell_wrap_3(&b_cur_PWM);
  emxFree_cell_wrap_0(&names);
  emxFree_cell_wrap_3(&p);
  emxFree_cell_wrap_3(&cur_PWM);
  emxFree_uint32_T(&lenvec);

  /* 'gkmPWM:204' fclose(fid); */
  cfclose(fileid);

  /* 'gkmPWM:205' fclose(fid2); */
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

static void eb_binary_expand_op(emxArray_real_T *kmat, int i, const
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

static void fb_binary_expand_op(emxArray_real_T *res, const emxArray_real_T *f,
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

static void gb_binary_expand_op(emxArray_real_T *kmat, const emxArray_real_T
  *ord, int ii, const emxArray_real_T *f)
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

/*
 * function [kweig,P] = getEMprob_v3(PWM,res,negmat,poscell,rc,diffc,indc,indloc,xc,reg,l_svm,k_svm,rcnum,RC)
 */
static void getEMprob_v3(const emxArray_real_T *PWM, const emxArray_real_T *res,
  const double negmat[16], const emxArray_cell_wrap_14 *poscell, const
  emxArray_real_T *rc, const emxArray_real_T *diffc, const emxArray_real_T *indc,
  const emxArray_real_T *indloc, const emxArray_real_T *xc, double reg, double
  l_svm, double k_svm, double rcnum, double RC, emxArray_real_T *kweig, double
  P[4])
{
  static const signed char b_b[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 1 };

  static const signed char c_b[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  static const signed char d_b[4] = { 1, 0, 0, 1 };

  emxArray_creal_T *b_A;
  emxArray_creal_T *c_A;
  emxArray_creal_T *r;
  emxArray_creal_T *vec;
  emxArray_int8_T *b_x;
  emxArray_int8_T *x;
  emxArray_real_T b_MAT_data;
  emxArray_real_T d_x_data;
  emxArray_real_T *A;
  emxArray_real_T *b;
  emxArray_real_T *c_x;
  creal_T B[16];
  creal_T b_mat[16];
  creal_T d_x[9];
  creal_T dcv[9];
  creal_T E_data[8];
  creal_T e_x_data[4];
  creal_T p[4];
  creal_T p2_data[4];
  creal_T ps[4];
  creal_T ps_data[4];
  creal_T y_data[4];
  creal_T p3_data[3];
  creal_T E;
  creal_T M;
  creal_T b_y;
  creal_T *b_A_data;
  creal_T *vec_data;
  double MAT_data[16];
  double mat[16];
  double tmp_data[16];
  double c_x_data[12];
  double y[4];
  const double *PWM_data;
  const double *res_data;
  double A_re_tmp;
  double E_im;
  double E_re;
  double E_re_tmp;
  double b_E_re_tmp;
  double b_y_tmp;
  double brm;
  double c_y_tmp;
  double d_y_tmp;
  double e_re;
  double y_tmp;
  double *A_data;
  double *b_data;
  double *kweig_data;
  int B_size[2];
  int MAT_size[2];
  int Posvec_size[2];
  int b_x_size[2];
  int tmp_size[2];
  int x_size[2];
  int b_I;
  int b_i;
  int i;
  int idx_tmp;
  int iindx;
  int ixlast;
  int k;
  int loop_ub;
  signed char ind[6];
  signed char Posvec_data[4];
  signed char b_tmp_data[3];
  signed char *b_x_data;
  signed char *x_data;
  bool b_b_data[4];
  res_data = res->data;
  PWM_data = PWM->data;

  /* Lagrange optimization (see paper for derivation) */
  /* 'gkmPWM:628' a = true; */
  /* 'gkmPWM:629' posvec = 1:4; */
  /* 'gkmPWM:630' if RC */
  emxInit_real_T(&A, 2);
  if (RC != 0.0) {
    /* 'gkmPWM:631' A =  ls_kweigtree(PWM,negmat,poscell,rc,diffc,indc,indloc,xc,l_svm,k_svm,rcnum); */
    ls_kweigtree(PWM, negmat, poscell, rc, diffc, indc, indloc, xc, l_svm, k_svm,
                 rcnum, A);
    A_data = A->data;
  } else {
    /* 'gkmPWM:632' else */
    /* 'gkmPWM:633' A =  ls_kweigtree_norc(PWM,negmat,poscell,rc,diffc,indc,indloc,xc,l_svm,k_svm,rcnum); */
    ls_kweigtree_norc(PWM, negmat, poscell, rc, diffc, indc, indloc, xc, l_svm,
                      k_svm, A);
    A_data = A->data;
  }

  /* 'gkmPWM:635' b = res+A*PWM(l_svm,:)'; */
  y_tmp = PWM_data[(int)l_svm - 1];
  y[0] = y_tmp;
  b_y_tmp = PWM_data[((int)l_svm + PWM->size[0]) - 1];
  y[1] = b_y_tmp;
  c_y_tmp = PWM_data[((int)l_svm + PWM->size[0] * 2) - 1];
  y[2] = c_y_tmp;
  d_y_tmp = PWM_data[((int)l_svm + PWM->size[0] * 3) - 1];
  y[3] = d_y_tmp;
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
    plus(b, res);
    b_data = b->data;
  }

  /* 'gkmPWM:636' mat = A'*A; */
  if (A->size[0] == 0) {
    memset(&mat[0], 0, 16U * sizeof(double));
  } else {
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)4, (blasint)4,
                (blasint)A->size[0], 1.0, &A_data[0], (blasint)A->size[0],
                &A_data[0], (blasint)A->size[0], 0.0, &mat[0], (blasint)4);
  }

  /* 'gkmPWM:637' y = A'*b; */
  b_mtimes(A, b, y);

  /* 'gkmPWM:638' if reg > 0 */
  if (reg > 0.0) {
    /* 'gkmPWM:639' M = min(eig(mat)); */
    eig(mat, p);
    M = c_minimum(p);

    /* 'gkmPWM:640' B=(mat-reg*M*eye(4))^-1; */
    b_y.re = reg * M.re;
    b_y.im = reg * M.im;
    for (i = 0; i < 16; i++) {
      k = b_b[i];
      b_mat[i].re = mat[i] - b_y.re * (double)k;
      b_mat[i].im = 0.0 - b_y.im * (double)k;
    }

    b_mpower(b_mat, B);
  } else {
    /* 'gkmPWM:641' else */
    /* 'gkmPWM:642' M = 0; */
    M.re = 0.0;
    M.im = 0.0;

    /* 'gkmPWM:643' B = mat^-1; */
    c_mpower(mat, tmp_data);
    for (i = 0; i < 16; i++) {
      B[i].re = tmp_data[i];
      B[i].im = 0.0;
    }
  }

  /* 'gkmPWM:645' ps = B*y; */
  brm = y[0];
  E_im = y[1];
  E_re_tmp = y[2];
  b_E_re_tmp = y[3];
  for (i = 0; i < 4; i++) {
    ps[i].re = ((B[i].re * brm + B[i + 4].re * E_im) + B[i + 8].re * E_re_tmp) +
      B[i + 12].re * b_E_re_tmp;
    ps[i].im = ((B[i].im * brm + B[i + 4].im * E_im) + B[i + 8].im * E_re_tmp) +
      B[i + 12].im * b_E_re_tmp;
  }

  /*  coder.varsize('p'); */
  /* 'gkmPWM:647' p = ps+(1-sum(ps))*B/sum(sum(B))*ones(4,1); */
  E.re = 1.0 - (((ps[0].re + ps[1].re) + ps[2].re) + ps[3].re);
  E.im = 0.0 - (((ps[0].im + ps[1].im) + ps[2].im) + ps[3].im);
  d_sum(B, ps_data);
  b_y.re = ((ps_data[0].re + ps_data[1].re) + ps_data[2].re) + ps_data[3].re;
  b_y.im = ((ps_data[0].im + ps_data[1].im) + ps_data[2].im) + ps_data[3].im;
  for (i = 0; i < 16; i++) {
    E_re_tmp = B[i].im;
    b_E_re_tmp = B[i].re;
    E_re = E.re * b_E_re_tmp - E.im * E_re_tmp;
    E_im = E.re * E_re_tmp + E.im * b_E_re_tmp;
    if (b_y.im == 0.0) {
      if (E_im == 0.0) {
        B[i].re = E_re / b_y.re;
        B[i].im = 0.0;
      } else if (E_re == 0.0) {
        B[i].re = 0.0;
        B[i].im = E_im / b_y.re;
      } else {
        B[i].re = E_re / b_y.re;
        B[i].im = E_im / b_y.re;
      }
    } else if (b_y.re == 0.0) {
      if (E_re == 0.0) {
        B[i].re = E_im / b_y.im;
        B[i].im = 0.0;
      } else if (E_im == 0.0) {
        B[i].re = 0.0;
        B[i].im = -(E_re / b_y.im);
      } else {
        B[i].re = E_im / b_y.im;
        B[i].im = -(E_re / b_y.im);
      }
    } else {
      brm = fabs(b_y.re);
      E_re_tmp = fabs(b_y.im);
      if (brm > E_re_tmp) {
        E_re_tmp = b_y.im / b_y.re;
        b_E_re_tmp = b_y.re + E_re_tmp * b_y.im;
        B[i].re = (E_re + E_re_tmp * E_im) / b_E_re_tmp;
        B[i].im = (E_im - E_re_tmp * E_re) / b_E_re_tmp;
      } else if (E_re_tmp == brm) {
        if (b_y.re > 0.0) {
          E_re_tmp = 0.5;
        } else {
          E_re_tmp = -0.5;
        }

        if (b_y.im > 0.0) {
          b_E_re_tmp = 0.5;
        } else {
          b_E_re_tmp = -0.5;
        }

        B[i].re = (E_re * E_re_tmp + E_im * b_E_re_tmp) / brm;
        B[i].im = (E_im * E_re_tmp - E_re * b_E_re_tmp) / brm;
      } else {
        E_re_tmp = b_y.re / b_y.im;
        b_E_re_tmp = b_y.im + E_re_tmp * b_y.re;
        B[i].re = (E_re_tmp * E_re + E_im) / b_E_re_tmp;
        B[i].im = (E_re_tmp * E_im - E_re) / b_E_re_tmp;
      }
    }
  }

  /* 'gkmPWM:648' P = zeros(4,1); */
  for (b_i = 0; b_i < 4; b_i++) {
    p[b_i].re = ps[b_i].re + (((B[b_i].re + B[b_i + 4].re) + B[b_i + 8].re) +
      B[b_i + 12].re);
    p[b_i].im = ps[b_i].im + (((B[b_i].im + B[b_i + 4].im) + B[b_i + 8].im) +
      B[b_i + 12].im);
    ps[b_i].re = 0.0;
    ps[b_i].im = 0.0;
  }

  /* 'gkmPWM:649' e = 0; */
  e_re = 0.0;

  /* The solution to the lagrange optimization problem */
  /* The following deals with the case if the optimal solution has a negative probability.  I cheat a little by only considering the cases where the base with the maximum solution is non-zero.  This works just fine in practice since the sum(p) = 1 constraint forces one of the bases to be positive.  It speeds up computation almost two fold. */
  /* 'gkmPWM:652' if min(p) < 0 */
  E = c_minimum(p);
  if (E.re < 0.0) {
    emxInit_int8_T(&x, 2);

    /* 'gkmPWM:653' I = 0; */
    b_I = 0;

    /* 'gkmPWM:654' [~,a] = max(p); */
    i_maximum(p, &E, &iindx);

    /* 'gkmPWM:655' nvec = posvec; */
    /* 'gkmPWM:656' nvec(a) = []; */
    i = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = 4;
    emxEnsureCapacity_int8_T(x, i);
    x_data = x->data;
    x_data[0] = 1;
    x_data[1] = 2;
    x_data[2] = 3;
    x_data[3] = 4;
    for (k = iindx; k < 4; k++) {
      x_data[k - 1] = x_data[k];
    }

    i = x->size[0] * x->size[1];
    x->size[1] = 3;
    emxEnsureCapacity_int8_T(x, i);
    x_data = x->data;

    /* Check cases where one of the probabilities is zero */
    /* 'gkmPWM:658' for i = 1:3 */
    emxInit_creal_T(&vec, 1);
    emxInit_int8_T(&b_x, 2);
    emxInit_real_T(&c_x, 2);
    emxInit_creal_T(&r, 1);
    emxInit_creal_T(&c_A, 2);
    Posvec_size[1] = 3;
    for (b_i = 0; b_i < 3; b_i++) {
      /* 'gkmPWM:659' Posvec = posvec; */
      /* 'gkmPWM:660' Posvec(nvec(i)) = []; */
      i = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 1;
      b_x->size[1] = 4;
      emxEnsureCapacity_int8_T(b_x, i);
      b_x_data = b_x->data;
      b_x_data[0] = 1;
      b_x_data[1] = 2;
      b_x_data[2] = 3;
      b_x_data[3] = 4;
      idx_tmp = x_data[b_i];
      for (k = idx_tmp; k < 4; k++) {
        b_x_data[k - 1] = b_x_data[k];
      }

      i = b_x->size[0] * b_x->size[1];
      b_x->size[1] = 3;
      emxEnsureCapacity_int8_T(b_x, i);
      b_x_data = b_x->data;

      /* 'gkmPWM:661' MAT = mat; */
      /* 'gkmPWM:662' MAT(:,nvec(i)) =[]; */
      /* 'gkmPWM:663' MAT(nvec(i),:) =[]; */
      i = c_x->size[0] * c_x->size[1];
      c_x->size[0] = 4;
      c_x->size[1] = 4;
      emxEnsureCapacity_real_T(c_x, i);
      kweig_data = c_x->data;
      for (i = 0; i < 16; i++) {
        kweig_data[i] = mat[i];
      }

      for (ixlast = idx_tmp; ixlast < 4; ixlast++) {
        kweig_data[c_x->size[0] * (ixlast - 1)] = kweig_data[c_x->size[0] *
          ixlast];
        kweig_data[c_x->size[0] * (ixlast - 1) + 1] = kweig_data[c_x->size[0] *
          ixlast + 1];
        kweig_data[c_x->size[0] * (ixlast - 1) + 2] = kweig_data[c_x->size[0] *
          ixlast + 2];
        kweig_data[c_x->size[0] * (ixlast - 1) + 3] = kweig_data[c_x->size[0] *
          ixlast + 3];
      }

      for (i = 0; i < 3; i++) {
        ixlast = i << 2;
        kweig_data[ixlast] = kweig_data[c_x->size[0] * i];
        kweig_data[ixlast + 1] = kweig_data[c_x->size[0] * i + 1];
        kweig_data[ixlast + 2] = kweig_data[c_x->size[0] * i + 2];
        kweig_data[ixlast + 3] = kweig_data[c_x->size[0] * i + 3];
      }

      i = c_x->size[0] * c_x->size[1];
      c_x->size[0] = 4;
      c_x->size[1] = 3;
      emxEnsureCapacity_real_T(c_x, i);
      kweig_data = c_x->data;
      b_x_size[0] = 4;
      b_x_size[1] = 3;
      for (i = 0; i < 12; i++) {
        c_x_data[i] = kweig_data[i];
      }

      for (ixlast = 0; ixlast < 3; ixlast++) {
        for (loop_ub = idx_tmp; loop_ub < 4; loop_ub++) {
          c_x_data[(loop_ub + b_x_size[0] * ixlast) - 1] = c_x_data[loop_ub +
            b_x_size[0] * ixlast];
        }
      }

      for (i = 0; i < 3; i++) {
        c_x_data[3 * i] = c_x_data[b_x_size[0] * i];
        c_x_data[3 * i + 1] = c_x_data[b_x_size[0] * i + 1];
        c_x_data[3 * i + 2] = c_x_data[b_x_size[0] * i + 2];
      }

      b_x_size[0] = 3;
      b_x_size[1] = 3;

      /* 'gkmPWM:664' Y = y(Posvec); */
      /* 'gkmPWM:665' if reg > 0 */
      if (reg > 0.0) {
        /* 'gkmPWM:666' B=(MAT-reg*M*eye(3))^-1; */
        b_y.re = reg * M.re;
        b_y.im = reg * M.im;
        for (i = 0; i < 3; i++) {
          for (k = 0; k < 3; k++) {
            idx_tmp = k + 3 * i;
            ixlast = c_b[idx_tmp];
            d_x[idx_tmp].re = c_x_data[k + b_x_size[0] * i] - b_y.re * (double)
              ixlast;
            d_x[idx_tmp].im = 0.0 - b_y.im * (double)ixlast;
          }
        }

        d_mpower(d_x, dcv);
        B_size[0] = 3;
        B_size[1] = 3;
        memcpy(&B[0], &dcv[0], 9U * sizeof(creal_T));
      } else {
        /* 'gkmPWM:667' else */
        /* 'gkmPWM:668' B=(MAT)^-1; */
        d_x_data.data = &c_x_data[0];
        d_x_data.size = &b_x_size[0];
        d_x_data.allocatedSize = 12;
        d_x_data.numDimensions = 2;
        d_x_data.canFreeData = false;
        mpower(&d_x_data, c_x);
        kweig_data = c_x->data;
        B_size[0] = c_x->size[0];
        B_size[1] = c_x->size[1];
        loop_ub = c_x->size[0] * c_x->size[1];
        for (i = 0; i < loop_ub; i++) {
          B[i].re = kweig_data[i];
          B[i].im = 0.0;
        }
      }

      /* 'gkmPWM:670' ps = B*Y; */
      for (i = 0; i < 3; i++) {
        y_data[i].re = y[b_x_data[i] - 1];
        y_data[i].im = 0.0;
      }

      loop_ub = B_size[0];
      for (i = 0; i < loop_ub; i++) {
        ps_data[i].re = 0.0;
        ps_data[i].im = 0.0;
        ixlast = B_size[1];
        for (k = 0; k < ixlast; k++) {
          idx_tmp = i + B_size[0] * k;
          E_re_tmp = B[idx_tmp].re;
          b_E_re_tmp = y_data[k].im;
          brm = B[idx_tmp].im;
          E_im = y_data[k].re;
          ps_data[i].re += E_re_tmp * E_im - brm * b_E_re_tmp;
          ps_data[i].im += E_re_tmp * b_E_re_tmp + brm * E_im;
        }
      }

      /* 'gkmPWM:671' p3 = ps+(1-sum(ps))*B/sum(sum(B))*ones(3,1); */
      E.re = 1.0 - ((ps_data[0].re + ps_data[1].re) + ps_data[2].re);
      E.im = 0.0 - ((ps_data[0].im + ps_data[1].im) + ps_data[2].im);
      e_sum(B, B_size, e_x_data, x_size);
      ixlast = x_size[1];
      if (x_size[1] == 0) {
        b_y.re = 0.0;
        b_y.im = 0.0;
      } else {
        b_y = e_x_data[0];
        for (k = 2; k <= ixlast; k++) {
          b_y.re += e_x_data[k - 1].re;
          b_y.im += e_x_data[k - 1].im;
        }
      }

      if (B_size[0] == 3) {
        loop_ub = B_size[0] * B_size[1];
        for (i = 0; i < loop_ub; i++) {
          brm = B[i].re;
          E_im = B[i].im;
          E_re = E.re * brm - E.im * E_im;
          E_im = E.re * E_im + E.im * brm;
          if (b_y.im == 0.0) {
            if (E_im == 0.0) {
              d_x[i].re = E_re / b_y.re;
              d_x[i].im = 0.0;
            } else if (E_re == 0.0) {
              d_x[i].re = 0.0;
              d_x[i].im = E_im / b_y.re;
            } else {
              d_x[i].re = E_re / b_y.re;
              d_x[i].im = E_im / b_y.re;
            }
          } else if (b_y.re == 0.0) {
            if (E_re == 0.0) {
              d_x[i].re = E_im / b_y.im;
              d_x[i].im = 0.0;
            } else if (E_im == 0.0) {
              d_x[i].re = 0.0;
              d_x[i].im = -(E_re / b_y.im);
            } else {
              d_x[i].re = E_im / b_y.im;
              d_x[i].im = -(E_re / b_y.im);
            }
          } else {
            brm = fabs(b_y.re);
            E_re_tmp = fabs(b_y.im);
            if (brm > E_re_tmp) {
              E_re_tmp = b_y.im / b_y.re;
              b_E_re_tmp = b_y.re + E_re_tmp * b_y.im;
              d_x[i].re = (E_re + E_re_tmp * E_im) / b_E_re_tmp;
              d_x[i].im = (E_im - E_re_tmp * E_re) / b_E_re_tmp;
            } else if (E_re_tmp == brm) {
              if (b_y.re > 0.0) {
                E_re_tmp = 0.5;
              } else {
                E_re_tmp = -0.5;
              }

              if (b_y.im > 0.0) {
                b_E_re_tmp = 0.5;
              } else {
                b_E_re_tmp = -0.5;
              }

              d_x[i].re = (E_re * E_re_tmp + E_im * b_E_re_tmp) / brm;
              d_x[i].im = (E_im * E_re_tmp - E_re * b_E_re_tmp) / brm;
            } else {
              E_re_tmp = b_y.re / b_y.im;
              b_E_re_tmp = b_y.im + E_re_tmp * b_y.re;
              d_x[i].re = (E_re_tmp * E_re + E_im) / b_E_re_tmp;
              d_x[i].im = (E_re_tmp * E_im - E_re) / b_E_re_tmp;
            }
          }
        }

        for (i = 0; i < 3; i++) {
          p3_data[i].re = ps_data[i].re + ((d_x[i].re + d_x[i + 3].re) + d_x[i +
            6].re);
          p3_data[i].im = ps_data[i].im + ((d_x[i].im + d_x[i + 3].im) + d_x[i +
            6].im);
        }
      } else {
        v_binary_expand_op(p3_data, &ixlast, ps_data, &B_size[0], E, B, B_size,
                           b_y);
      }

      /* solution */
      /* Checks if solution is permitted.  If so, makes sure that it creates a smaller error than other cases */
      /* 'gkmPWM:673' if min(p3) >= 0 */
      E = p3_data[0];
      absRelopProxies(p3_data[0], p3_data[1], &E_re_tmp, &b_E_re_tmp);
      if (E_re_tmp > b_E_re_tmp) {
        E = p3_data[1];
      }

      absRelopProxies(E, p3_data[2], &E_re_tmp, &b_E_re_tmp);
      if (E_re_tmp > b_E_re_tmp) {
        E = p3_data[2];
      }

      if (E.re >= 0.0) {
        /* 'gkmPWM:674' vec = b-A(:,Posvec)*p3; */
        loop_ub = A->size[0];
        i = c_A->size[0] * c_A->size[1];
        c_A->size[0] = A->size[0];
        c_A->size[1] = 3;
        emxEnsureCapacity_creal_T(c_A, i);
        b_A_data = c_A->data;
        for (i = 0; i < 3; i++) {
          for (k = 0; k < loop_ub; k++) {
            b_A_data[k + c_A->size[0] * i].re = A_data[k + A->size[0] *
              (b_x_data[i] - 1)];
            b_A_data[k + c_A->size[0] * i].im = 0.0;
          }
        }

        i = r->size[0];
        r->size[0] = c_A->size[0];
        emxEnsureCapacity_creal_T(r, i);
        vec_data = r->data;
        loop_ub = c_A->size[0];
        for (i = 0; i < loop_ub; i++) {
          vec_data[i].re = 0.0;
          vec_data[i].im = 0.0;
          for (k = 0; k < 3; k++) {
            E_re_tmp = b_A_data[i + c_A->size[0] * k].re;
            vec_data[i].re += E_re_tmp * p3_data[k].re;
            vec_data[i].im += E_re_tmp * p3_data[k].im;
          }
        }

        if (b->size[0] == r->size[0]) {
          loop_ub = A->size[0];
          i = c_A->size[0] * c_A->size[1];
          c_A->size[0] = A->size[0];
          c_A->size[1] = 3;
          emxEnsureCapacity_creal_T(c_A, i);
          b_A_data = c_A->data;
          for (i = 0; i < 3; i++) {
            for (k = 0; k < loop_ub; k++) {
              b_A_data[k + c_A->size[0] * i].re = A_data[k + A->size[0] *
                (b_x_data[i] - 1)];
              b_A_data[k + c_A->size[0] * i].im = 0.0;
            }
          }

          i = vec->size[0];
          vec->size[0] = c_A->size[0];
          emxEnsureCapacity_creal_T(vec, i);
          vec_data = vec->data;
          loop_ub = c_A->size[0];
          for (i = 0; i < loop_ub; i++) {
            b_E_re_tmp = 0.0;
            E_im = 0.0;
            for (k = 0; k < 3; k++) {
              E_re_tmp = b_A_data[i + c_A->size[0] * k].re;
              b_E_re_tmp += E_re_tmp * p3_data[k].re;
              E_im += E_re_tmp * p3_data[k].im;
            }

            vec_data[i].re = b_data[i] - b_E_re_tmp;
            vec_data[i].im = 0.0 - E_im;
          }
        } else {
          u_binary_expand_op(vec, b, A, b_x, p3_data);
          vec_data = vec->data;
        }

        /* 'gkmPWM:675' if I == 0 */
        if (b_I == 0) {
          /* 'gkmPWM:676' I = 1; */
          b_I = 1;

          /* 'gkmPWM:677' P = zeros(4,1); */
          memset(&ps[0], 0, 4U * sizeof(creal_T));

          /* 'gkmPWM:678' P(Posvec) = p3; */
          loop_ub = b_x->size[1];
          for (i = 0; i < loop_ub; i++) {
            b_tmp_data[i] = b_x_data[i];
          }

          ps[b_tmp_data[0] - 1] = p3_data[0];
          ps[b_tmp_data[1] - 1] = p3_data[1];
          ps[b_tmp_data[2] - 1] = p3_data[2];

          /* 'gkmPWM:679' e = vec'*vec-reg*M*p3'*p3; */
          b_y.re = 0.0;
          if (vec->size[0] >= 1) {
            ixlast = vec->size[0];
            for (k = 0; k < ixlast; k++) {
              E_re_tmp = vec_data[k].re;
              b_E_re_tmp = vec_data[k].im;
              b_y.re += E_re_tmp * E_re_tmp + b_E_re_tmp * b_E_re_tmp;
            }
          }

          E.re = reg * M.re;
          E.im = reg * M.im;
          e_x_data[0].re = E.re * p3_data[0].re - E.im * -p3_data[0].im;
          e_x_data[0].im = E.re * -p3_data[0].im + E.im * p3_data[0].re;
          e_x_data[1].re = E.re * p3_data[1].re - E.im * -p3_data[1].im;
          e_x_data[1].im = E.re * -p3_data[1].im + E.im * p3_data[1].re;
          e_x_data[2].re = E.re * p3_data[2].re - E.im * -p3_data[2].im;
          e_x_data[2].im = E.re * -p3_data[2].im + E.im * p3_data[2].re;
          e_re = b_y.re - (((e_x_data[0].re * p3_data[0].re - e_x_data[0].im *
                             p3_data[0].im) + (e_x_data[1].re * p3_data[1].re -
            e_x_data[1].im * p3_data[1].im)) + (e_x_data[2].re * p3_data[2].re -
            e_x_data[2].im * p3_data[2].im));
        } else {
          /* 'gkmPWM:680' else */
          /* 'gkmPWM:681' E = vec'*vec-reg*M*p3'*p3; */
          b_y.re = 0.0;
          if (vec->size[0] >= 1) {
            ixlast = vec->size[0];
            for (k = 0; k < ixlast; k++) {
              E_re_tmp = vec_data[k].re;
              b_E_re_tmp = vec_data[k].im;
              b_y.re += E_re_tmp * E_re_tmp + b_E_re_tmp * b_E_re_tmp;
            }
          }

          E.re = reg * M.re;
          E.im = reg * M.im;
          e_x_data[0].re = E.re * p3_data[0].re - E.im * -p3_data[0].im;
          e_x_data[0].im = E.re * -p3_data[0].im + E.im * p3_data[0].re;
          e_x_data[1].re = E.re * p3_data[1].re - E.im * -p3_data[1].im;
          e_x_data[1].im = E.re * -p3_data[1].im + E.im * p3_data[1].re;
          e_x_data[2].re = E.re * p3_data[2].re - E.im * -p3_data[2].im;
          e_x_data[2].im = E.re * -p3_data[2].im + E.im * p3_data[2].re;
          E.re = b_y.re - (((e_x_data[0].re * p3_data[0].re - e_x_data[0].im *
                             p3_data[0].im) + (e_x_data[1].re * p3_data[1].re -
            e_x_data[1].im * p3_data[1].im)) + (e_x_data[2].re * p3_data[2].re -
            e_x_data[2].im * p3_data[2].im));

          /* 'gkmPWM:682' if E < e */
          if (E.re < e_re) {
            /* 'gkmPWM:683' e = E; */
            e_re = E.re;

            /* 'gkmPWM:684' P = zeros(4,1); */
            memset(&ps[0], 0, 4U * sizeof(creal_T));

            /* 'gkmPWM:685' P(Posvec) = p3; */
            loop_ub = b_x->size[1];
            for (i = 0; i < loop_ub; i++) {
              b_tmp_data[i] = b_x_data[i];
            }

            ps[b_tmp_data[0] - 1] = p3_data[0];
            ps[b_tmp_data[1] - 1] = p3_data[1];
            ps[b_tmp_data[2] - 1] = p3_data[2];
          }
        }
      }
    }

    emxFree_int8_T(&b_x);

    /* Check cases where two of the probabilities are zero */
    /* 'gkmPWM:691' ind = [nvec(1) nvec(2);nvec(1) nvec(3);nvec(2) nvec(3)]; */
    ind[0] = x_data[0];
    ind[3] = x_data[1];
    ind[1] = x_data[0];
    ind[4] = x_data[2];
    ind[2] = x_data[1];
    ind[5] = x_data[2];

    /* 'gkmPWM:692' for i = 1:3 */
    emxFree_int8_T(&x);
    Posvec_size[0] = 1;
    for (b_i = 0; b_i < 3; b_i++) {
      /* 'gkmPWM:693' Posvec = posvec; */
      /* 'gkmPWM:694' Posvec(ind(i,:)) = []; */
      x_size[0] = ind[b_i];
      x_size[1] = ind[b_i + 3];
      Posvec_data[0] = 1;
      b_b_data[0] = false;
      Posvec_data[1] = 2;
      b_b_data[1] = false;
      Posvec_data[2] = 3;
      b_b_data[2] = false;
      Posvec_data[3] = 4;
      b_b_data[3] = false;
      b_b_data[x_size[0] - 1] = true;
      b_b_data[x_size[1] - 1] = true;
      ixlast = -1;
      if (!b_b_data[0]) {
        ixlast = 0;
      }

      if (!b_b_data[1]) {
        ixlast++;
        Posvec_data[ixlast] = 2;
      }

      if (!b_b_data[2]) {
        ixlast++;
        Posvec_data[ixlast] = 3;
      }

      idx_tmp = ((b_b_data[0] + b_b_data[1]) + b_b_data[2]) + b_b_data[3];
      if (!b_b_data[3]) {
        ixlast++;
        Posvec_data[ixlast] = 4;
      }

      if (1 > 4 - idx_tmp) {
        Posvec_size[1] = 0;
      } else {
        Posvec_size[1] = 4 - idx_tmp;
      }

      /* 'gkmPWM:695' MAT = mat; */
      /* 'gkmPWM:696' MAT(:,ind(i,:)) =[]; */
      /* 'gkmPWM:697' MAT(ind(i,:),:) =[]; */
      tmp_size[0] = 4;
      tmp_size[1] = 4;
      memcpy(&tmp_data[0], &mat[0], 16U * sizeof(double));
      b_nullAssignment(tmp_data, tmp_size, x_size);
      MAT_size[0] = 4;
      MAT_size[1] = tmp_size[1];
      loop_ub = 4 * tmp_size[1];
      if (0 <= loop_ub - 1) {
        memcpy(&MAT_data[0], &tmp_data[0], loop_ub * sizeof(double));
      }

      c_nullAssignment(MAT_data, MAT_size, x_size);

      /* 'gkmPWM:698' Y = y(Posvec); */
      /* 'gkmPWM:699' if reg > 0 */
      if (reg > 0.0) {
        /* 'gkmPWM:700' B=(MAT-reg*M*eye(2))^-1; */
        b_y.re = reg * M.re;
        b_y.im = reg * M.im;
        if ((MAT_size[0] == 2) && (MAT_size[1] == 2)) {
          for (i = 0; i < 2; i++) {
            for (k = 0; k < 2; k++) {
              idx_tmp = k + (i << 1);
              ixlast = d_b[idx_tmp];
              ps_data[idx_tmp].re = MAT_data[k + MAT_size[0] * i] - b_y.re *
                (double)ixlast;
              ps_data[idx_tmp].im = 0.0 - b_y.im * (double)ixlast;
            }
          }

          e_mpower(ps_data, p);
        } else {
          x_binary_expand_op(p, MAT_data, MAT_size, b_y, d_b);
        }

        B_size[0] = 2;
        B_size[1] = 2;
        memcpy(&B[0], &p[0], 4U * sizeof(creal_T));
      } else {
        /* 'gkmPWM:701' else */
        /* 'gkmPWM:702' B=MAT^-1; */
        b_MAT_data.data = &MAT_data[0];
        b_MAT_data.size = &MAT_size[0];
        b_MAT_data.allocatedSize = 16;
        b_MAT_data.numDimensions = 2;
        b_MAT_data.canFreeData = false;
        mpower(&b_MAT_data, c_x);
        kweig_data = c_x->data;
        B_size[0] = c_x->size[0];
        B_size[1] = c_x->size[1];
        loop_ub = c_x->size[0] * c_x->size[1];
        for (i = 0; i < loop_ub; i++) {
          B[i].re = kweig_data[i];
          B[i].im = 0.0;
        }
      }

      /* 'gkmPWM:704' ps = B*Y; */
      loop_ub = Posvec_size[1];
      for (i = 0; i < loop_ub; i++) {
        y_data[i].re = y[Posvec_data[i] - 1];
        y_data[i].im = 0.0;
      }

      loop_ub = B_size[0];
      for (i = 0; i < loop_ub; i++) {
        ps_data[i].re = 0.0;
        ps_data[i].im = 0.0;
        ixlast = B_size[1];
        for (k = 0; k < ixlast; k++) {
          idx_tmp = i + B_size[0] * k;
          E_re_tmp = B[idx_tmp].re;
          b_E_re_tmp = y_data[k].im;
          brm = B[idx_tmp].im;
          E_im = y_data[k].re;
          ps_data[i].re += E_re_tmp * E_im - brm * b_E_re_tmp;
          ps_data[i].im += E_re_tmp * b_E_re_tmp + brm * E_im;
        }
      }

      /* 'gkmPWM:705' p2 = ps+(1-sum(ps))*B/sum(sum(B))*ones(2,1); */
      if (B_size[0] == 0) {
        b_y.re = 0.0;
        b_y.im = 0.0;
      } else {
        b_y = ps_data[0];
        for (k = 2; k <= loop_ub; k++) {
          b_y.re += ps_data[k - 1].re;
          b_y.im += ps_data[k - 1].im;
        }
      }

      E.re = 1.0 - b_y.re;
      E.im = 0.0 - b_y.im;
      e_sum(B, B_size, e_x_data, x_size);
      ixlast = x_size[1];
      if (x_size[1] == 0) {
        b_y.re = 0.0;
        b_y.im = 0.0;
      } else {
        b_y = e_x_data[0];
        for (k = 2; k <= ixlast; k++) {
          b_y.re += e_x_data[k - 1].re;
          b_y.im += e_x_data[k - 1].im;
        }
      }

      ixlast = B_size[0];
      loop_ub = B_size[0] * B_size[1];
      for (i = 0; i < loop_ub; i++) {
        brm = B[i].re;
        E_im = B[i].im;
        E_re = E.re * brm - E.im * E_im;
        E_im = E.re * E_im + E.im * brm;
        if (b_y.im == 0.0) {
          if (E_im == 0.0) {
            E_data[i].re = E_re / b_y.re;
            E_data[i].im = 0.0;
          } else if (E_re == 0.0) {
            E_data[i].re = 0.0;
            E_data[i].im = E_im / b_y.re;
          } else {
            E_data[i].re = E_re / b_y.re;
            E_data[i].im = E_im / b_y.re;
          }
        } else if (b_y.re == 0.0) {
          if (E_re == 0.0) {
            E_data[i].re = E_im / b_y.im;
            E_data[i].im = 0.0;
          } else if (E_im == 0.0) {
            E_data[i].re = 0.0;
            E_data[i].im = -(E_re / b_y.im);
          } else {
            E_data[i].re = E_im / b_y.im;
            E_data[i].im = -(E_re / b_y.im);
          }
        } else {
          brm = fabs(b_y.re);
          E_re_tmp = fabs(b_y.im);
          if (brm > E_re_tmp) {
            E_re_tmp = b_y.im / b_y.re;
            b_E_re_tmp = b_y.re + E_re_tmp * b_y.im;
            E_data[i].re = (E_re + E_re_tmp * E_im) / b_E_re_tmp;
            E_data[i].im = (E_im - E_re_tmp * E_re) / b_E_re_tmp;
          } else if (E_re_tmp == brm) {
            if (b_y.re > 0.0) {
              E_re_tmp = 0.5;
            } else {
              E_re_tmp = -0.5;
            }

            if (b_y.im > 0.0) {
              b_E_re_tmp = 0.5;
            } else {
              b_E_re_tmp = -0.5;
            }

            E_data[i].re = (E_re * E_re_tmp + E_im * b_E_re_tmp) / brm;
            E_data[i].im = (E_im * E_re_tmp - E_re * b_E_re_tmp) / brm;
          } else {
            E_re_tmp = b_y.re / b_y.im;
            b_E_re_tmp = b_y.im + E_re_tmp * b_y.re;
            E_data[i].re = (E_re_tmp * E_re + E_im) / b_E_re_tmp;
            E_data[i].im = (E_re_tmp * E_im - E_re) / b_E_re_tmp;
          }
        }
      }

      idx_tmp = B_size[0];
      for (i = 0; i < ixlast; i++) {
        k = i + ixlast;
        p2_data[i].re = ps_data[i].re + (E_data[i].re + E_data[k].re);
        p2_data[i].im = ps_data[i].im + (E_data[i].im + E_data[k].im);
      }

      /* solution */
      /* Checks if solution is permitted.  If so, makes sure that it creates a smaller error than other cases */
      /* 'gkmPWM:707' if min(p2) >= 0 */
      E = p2_data[0];
      for (k = 2; k <= idx_tmp; k++) {
        b_y = p2_data[k - 1];
        absRelopProxies(E, b_y, &E_re_tmp, &b_E_re_tmp);
        if (E_re_tmp > b_E_re_tmp) {
          E = b_y;
        }
      }

      if (E.re >= 0.0) {
        /* 'gkmPWM:708' vec = b-A(:,Posvec)*p2; */
        loop_ub = A->size[0];
        i = c_A->size[0] * c_A->size[1];
        c_A->size[0] = A->size[0];
        ixlast = Posvec_size[1];
        c_A->size[1] = Posvec_size[1];
        emxEnsureCapacity_creal_T(c_A, i);
        b_A_data = c_A->data;
        for (i = 0; i < ixlast; i++) {
          for (k = 0; k < loop_ub; k++) {
            b_A_data[k + c_A->size[0] * i].re = A_data[k + A->size[0] *
              (Posvec_data[i] - 1)];
            b_A_data[k + c_A->size[0] * i].im = 0.0;
          }
        }

        i = r->size[0];
        r->size[0] = c_A->size[0];
        emxEnsureCapacity_creal_T(r, i);
        vec_data = r->data;
        loop_ub = c_A->size[0];
        for (i = 0; i < loop_ub; i++) {
          vec_data[i].re = 0.0;
          vec_data[i].im = 0.0;
          ixlast = c_A->size[1];
          for (k = 0; k < ixlast; k++) {
            E_re_tmp = b_A_data[i + c_A->size[0] * k].re;
            brm = p2_data[k].im;
            E_re = b_A_data[i + c_A->size[0] * k].im;
            A_re_tmp = p2_data[k].re;
            vec_data[i].re += E_re_tmp * A_re_tmp - E_re * brm;
            vec_data[i].im += E_re_tmp * brm + E_re * A_re_tmp;
          }
        }

        if (b->size[0] == r->size[0]) {
          loop_ub = A->size[0];
          i = c_A->size[0] * c_A->size[1];
          c_A->size[0] = A->size[0];
          ixlast = Posvec_size[1];
          c_A->size[1] = Posvec_size[1];
          emxEnsureCapacity_creal_T(c_A, i);
          b_A_data = c_A->data;
          for (i = 0; i < ixlast; i++) {
            for (k = 0; k < loop_ub; k++) {
              b_A_data[k + c_A->size[0] * i].re = A_data[k + A->size[0] *
                (Posvec_data[i] - 1)];
              b_A_data[k + c_A->size[0] * i].im = 0.0;
            }
          }

          i = vec->size[0];
          vec->size[0] = c_A->size[0];
          emxEnsureCapacity_creal_T(vec, i);
          vec_data = vec->data;
          loop_ub = c_A->size[0];
          for (i = 0; i < loop_ub; i++) {
            b_E_re_tmp = 0.0;
            E_im = 0.0;
            ixlast = c_A->size[1];
            for (k = 0; k < ixlast; k++) {
              E_re_tmp = b_A_data[i + c_A->size[0] * k].re;
              brm = p2_data[k].im;
              E_re = b_A_data[i + c_A->size[0] * k].im;
              A_re_tmp = p2_data[k].re;
              b_E_re_tmp += E_re_tmp * A_re_tmp - E_re * brm;
              E_im += E_re_tmp * brm + E_re * A_re_tmp;
            }

            vec_data[i].re = b_data[i] - b_E_re_tmp;
            vec_data[i].im = 0.0 - E_im;
          }
        } else {
          w_binary_expand_op(vec, b, A, Posvec_data, Posvec_size, p2_data);
          vec_data = vec->data;
        }

        /* 'gkmPWM:709' if I == 0 */
        if (b_I == 0) {
          /* 'gkmPWM:710' I = 1; */
          b_I = 1;

          /* 'gkmPWM:711' P = zeros(4,1); */
          memset(&ps[0], 0, 4U * sizeof(creal_T));

          /* 'gkmPWM:712' P(Posvec) = p2; */
          loop_ub = Posvec_size[1];
          for (i = 0; i < loop_ub; i++) {
            ps[Posvec_data[i] - 1] = p2_data[i];
          }

          /* 'gkmPWM:713' e = vec'*vec-reg*M*p2'*p2; */
          b_y.re = 0.0;
          if (vec->size[0] >= 1) {
            ixlast = vec->size[0];
            for (k = 0; k < ixlast; k++) {
              E_re_tmp = vec_data[k].re;
              b_E_re_tmp = vec_data[k].im;
              b_y.re += E_re_tmp * E_re_tmp + b_E_re_tmp * b_E_re_tmp;
            }
          }

          E.re = reg * M.re;
          E.im = reg * M.im;
          for (i = 0; i < idx_tmp; i++) {
            E_re_tmp = p2_data[i].re;
            b_E_re_tmp = -p2_data[i].im;
            e_x_data[i].re = E.re * E_re_tmp - E.im * b_E_re_tmp;
            e_x_data[i].im = E.re * b_E_re_tmp + E.im * E_re_tmp;
          }

          E.re = 0.0;
          if (idx_tmp >= 1) {
            for (k = 0; k < idx_tmp; k++) {
              E.re += e_x_data[k].re * p2_data[k].re - e_x_data[k].im *
                p2_data[k].im;
            }
          }

          e_re = b_y.re - E.re;
        } else {
          /* 'gkmPWM:714' else */
          /* 'gkmPWM:715' E = vec'*vec-reg*M*p2'*p2; */
          b_y.re = 0.0;
          if (vec->size[0] >= 1) {
            ixlast = vec->size[0];
            for (k = 0; k < ixlast; k++) {
              E_re_tmp = vec_data[k].re;
              b_E_re_tmp = vec_data[k].im;
              b_y.re += E_re_tmp * E_re_tmp + b_E_re_tmp * b_E_re_tmp;
            }
          }

          E.re = reg * M.re;
          E.im = reg * M.im;
          for (i = 0; i < idx_tmp; i++) {
            E_re_tmp = p2_data[i].re;
            b_E_re_tmp = -p2_data[i].im;
            e_x_data[i].re = E.re * E_re_tmp - E.im * b_E_re_tmp;
            e_x_data[i].im = E.re * b_E_re_tmp + E.im * E_re_tmp;
          }

          E.re = 0.0;
          if (idx_tmp >= 1) {
            for (k = 0; k < idx_tmp; k++) {
              E.re += e_x_data[k].re * p2_data[k].re - e_x_data[k].im *
                p2_data[k].im;
            }
          }

          E.re = b_y.re - E.re;

          /* 'gkmPWM:716' if E < e */
          if (E.re < e_re) {
            /* 'gkmPWM:717' e = E; */
            e_re = E.re;

            /* 'gkmPWM:718' P = zeros(4,1); */
            memset(&ps[0], 0, 4U * sizeof(creal_T));

            /* 'gkmPWM:719' P(Posvec) = p2; */
            loop_ub = Posvec_size[1];
            for (i = 0; i < loop_ub; i++) {
              ps[Posvec_data[i] - 1] = p2_data[i];
            }
          }
        }
      }
    }

    emxFree_creal_T(&c_A);
    emxFree_creal_T(&r);
    emxFree_real_T(&c_x);
    emxFree_creal_T(&vec);

    /* Checks to see if one non-zero case is better than the other cases */
    /* 'gkmPWM:725' if I == 0 */
    if (b_I == 0) {
      /* 'gkmPWM:726' P = zeros(4,1); */
      memset(&ps[0], 0, 4U * sizeof(creal_T));

      /* 'gkmPWM:727' P(a) = 1; */
      ps[iindx - 1].re = 1.0;
      ps[iindx - 1].im = 0.0;
    } else {
      /* 'gkmPWM:728' else */
      /* 'gkmPWM:729' vec = b-A(:,a); */
      if (b->size[0] == A->size[0]) {
        loop_ub = b->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_data[i] -= A_data[i + A->size[0] * (iindx - 1)];
        }
      } else {
        y_binary_expand_op(b, A, iindx);
        b_data = b->data;
      }

      /* 'gkmPWM:730' E = vec'*vec-reg*M; */
      /* 'gkmPWM:731' if E < e */
      if (b->size[0] < 1) {
        E_re_tmp = 0.0;
      } else {
        E_re_tmp = cblas_ddot((blasint)b->size[0], &b_data[0], (blasint)1,
                              &b_data[0], (blasint)1);
      }

      if (E_re_tmp - reg * M.re < e_re) {
        /* 'gkmPWM:732' e = E; */
        /* 'gkmPWM:733' P = zeros(4,1); */
        memset(&ps[0], 0, 4U * sizeof(creal_T));

        /* 'gkmPWM:734' P(Posvec) = p2; */
        loop_ub = Posvec_size[1];
        for (i = 0; i < loop_ub; i++) {
          ps[Posvec_data[i] - 1] = p2_data[i];
        }
      }
    }
  } else {
    /* 'gkmPWM:737' else */
    /* 'gkmPWM:738' P = p; */
    memcpy(&ps[0], &p[0], 4U * sizeof(creal_T));
  }

  emxFree_real_T(&b);
  emxInit_creal_T(&b_A, 2);

  /* 'gkmPWM:740' kweig = A*(P-PWM(l_svm,:)'); */
  /* 'gkmPWM:741' P = real(P); */
  P[0] = ps[0].re;
  P[1] = ps[1].re;
  P[2] = ps[2].re;
  P[3] = ps[3].re;

  /* 'gkmPWM:742' kweig = real(kweig); */
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = 4;
  emxEnsureCapacity_creal_T(b_A, i);
  b_A_data = b_A->data;
  loop_ub = A->size[0] * 4;
  for (i = 0; i < loop_ub; i++) {
    b_A_data[i].re = A_data[i];
    b_A_data[i].im = 0.0;
  }

  emxFree_real_T(&A);
  ps[0].re -= y_tmp;
  ps[1].re -= b_y_tmp;
  ps[2].re -= c_y_tmp;
  ps[3].re -= d_y_tmp;
  i = kweig->size[0];
  kweig->size[0] = b_A->size[0];
  emxEnsureCapacity_real_T(kweig, i);
  kweig_data = kweig->data;
  loop_ub = b_A->size[0];
  for (i = 0; i < loop_ub; i++) {
    kweig_data[i] = (((b_A_data[i].re * ps[0].re - b_A_data[i].im * ps[0].im) +
                      (b_A_data[i + b_A->size[0]].re * ps[1].re - b_A_data[i +
                       b_A->size[0]].im * ps[1].im)) + (b_A_data[i + b_A->size[0]
      * 2].re * ps[2].re - b_A_data[i + b_A->size[0] * 2].im * ps[2].im)) +
      (b_A_data[i + b_A->size[0] * 3].re * ps[3].re - b_A_data[i + b_A->size[0] *
       3].im * ps[3].im);
  }

  emxFree_creal_T(&b_A);
}

/*
 * function [PWM, scorevec, C, r, R, E, Rd] = gkmPWM_lagrange(kweig,negmat,PWM,negvec,n,rcorr,reg,l_svm,k_svm,RC)
 */
static void gkmPWM_lagrange(const emxArray_real_T *kweig, const double negmat[16],
  emxArray_cell_wrap_2 *PWM, const emxArray_real_T *negvec, double n, double
  rcorr, double reg, double l_svm, double k_svm, double RC, emxArray_real_T
  *scorevec, emxArray_real_T *C, double *r, emxArray_real_T *R, emxArray_real_T *
  E, emxArray_real_T *Rd)
{
  cell_wrap_14 *loc_data;
  cell_wrap_14 *poscell_data;
  cell_wrap_2 *PWM_data;
  cell_wrap_2 *b_PWM_data;
  cell_wrap_2 *new_PWM_data;
  emxArray_boolean_T *b_ct;
  emxArray_boolean_T *d_ct;
  emxArray_cell_wrap_14 *b_new_loc;
  emxArray_cell_wrap_14 *loc;
  emxArray_cell_wrap_14 *new_loc;
  emxArray_cell_wrap_14 *poscell;
  emxArray_cell_wrap_2 *b_PWM;
  emxArray_cell_wrap_2 *b_new_PWM;
  emxArray_cell_wrap_2 *new_PWM;
  emxArray_int32_T *c_iidx;
  emxArray_int32_T *iidx;
  emxArray_real_T *CT;
  emxArray_real_T *b_a;
  emxArray_real_T *b_iidx;
  emxArray_real_T *b_loc;
  emxArray_real_T *b_mat;
  emxArray_real_T *c_PWM;
  emxArray_real_T *c_ct;
  emxArray_real_T *ct;
  emxArray_real_T *diffC;
  emxArray_real_T *diffc;
  emxArray_real_T *f;
  emxArray_real_T *indc;
  emxArray_real_T *kmat;
  emxArray_real_T *lenvec;
  emxArray_real_T *mat;
  emxArray_real_T *ord;
  emxArray_real_T *rc;
  emxArray_real_T *res;
  emxArray_real_T *temp;
  emxArray_real_T *vec;
  emxArray_real_T *xc;
  double b_varargin_1[9];
  double b[4];
  double varargin_1[4];
  double b_lenvec[2];
  const double *kweig_data;
  const double *negvec_data;
  double M;
  double S;
  double c;
  double rcnum;
  double scount;
  double *CT_data;
  double *C_data;
  double *b_iidx_data;
  double *ct_data;
  double *diffC_data;
  double *f_data;
  double *kmat_data;
  double *lenvec_data;
  double *ord_data;
  double *res_data;
  double *scorevec_data;
  double *temp_data;
  double *vec_data;
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
  int m;
  int u0;
  int unnamed_idx_0;
  int varargin_2;
  int *iidx_data;
  bool x[9];
  bool guard1 = false;
  bool *b_ct_data;
  negvec_data = negvec->data;
  PWM_data = PWM->data;
  kweig_data = kweig->data;
  emxInit_cell_wrap_2(&b_PWM);
  i = b_PWM->size[0];
  b_PWM->size[0] = PWM->size[0];
  emxEnsureCapacity_cell_wrap_2(b_PWM, i);
  b_PWM_data = b_PWM->data;
  loop_ub = PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    emxCopyStruct_cell_wrap_2(&b_PWM_data[i], &PWM_data[i]);
  }

  emxInit_real_T(&diffC, 2);
  emxInit_real_T(&ct, 2);
  emxInit_real_T(&rc, 2);
  emxInit_real_T(&diffc, 1);
  emxInit_real_T(&indc, 1);
  emxInit_real_T(&xc, 2);

  /* Note: This code is rather messy.  I block commmented to the best of my ability, so hopefully this makes enough sense.  If something seems non-trivial, then I probably found a mathematical trick to speed up computation time (in particular dynamic programming). */
  /* 'gkmPWM:388' GC = PWM{1}(1,:); */
  /* 'gkmPWM:389' [~,rc,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm); */
  genIndex(l_svm, k_svm, ct, rc, diffc, indc, xc, &rcnum);
  vec_data = rc->data;

  /* generate index information */
  /* 'gkmPWM:390' lcomb = length(diffc); */
  lcomb = diffc->size[0];

  /* 'gkmPWM:391' diffC = zeros(lcomb,l_svm); */
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

  /* 'gkmPWM:392' for i = 1:l_svm */
  emxInit_real_T(&f, 1);
  emxInit_real_T(&CT, 2);
  emxInit_int32_T(&iidx, 1);
  emxInit_boolean_T(&b_ct, 2);
  emxInit_real_T(&b_iidx, 1);
  emxInit_real_T(&c_ct, 2);
  for (acount = 0; acount < i1; acount++) {
    /* 'gkmPWM:393' ct = rc+i-1; */
    i = ct->size[0] * ct->size[1];
    ct->size[0] = rc->size[0];
    ct->size[1] = rc->size[1];
    emxEnsureCapacity_real_T(ct, i);
    ct_data = ct->data;
    loop_ub = rc->size[0] * rc->size[1];
    for (i = 0; i < loop_ub; i++) {
      ct_data[i] = (vec_data[i] + ((double)acount + 1.0)) - 1.0;
    }

    /* 'gkmPWM:394' f = find(sum(ct==l_svm,2)); */
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

    /* 'gkmPWM:395' ct = ct(f,:); */
    m = ct->size[1] - 1;
    i = c_ct->size[0] * c_ct->size[1];
    c_ct->size[0] = f->size[0];
    c_ct->size[1] = ct->size[1];
    emxEnsureCapacity_real_T(c_ct, i);
    temp_data = c_ct->data;
    for (i = 0; i <= m; i++) {
      loop_ub = f->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        temp_data[i2 + c_ct->size[0] * i] = ct_data[((int)f_data[i2] + ct->size
          [0] * i) - 1];
      }
    }

    i = ct->size[0] * ct->size[1];
    ct->size[0] = c_ct->size[0];
    ct->size[1] = c_ct->size[1];
    emxEnsureCapacity_real_T(ct, i);
    ct_data = ct->data;
    loop_ub = c_ct->size[0] * c_ct->size[1];
    for (i = 0; i < loop_ub; i++) {
      ct_data[i] = temp_data[i];
    }

    /* 'gkmPWM:396' CT = zeros(length(ct),k_svm-1); */
    if ((ct->size[0] == 0) || (ct->size[1] == 0)) {
      m = 0;
    } else {
      u0 = ct->size[0];
      m = ct->size[1];
      if (u0 >= m) {
        m = u0;
      }
    }

    i = CT->size[0] * CT->size[1];
    CT->size[0] = m;
    CT->size[1] = (int)(k_svm - 1.0);
    emxEnsureCapacity_real_T(CT, i);
    CT_data = CT->data;
    loop_ub = m * (int)(k_svm - 1.0);
    for (i = 0; i < loop_ub; i++) {
      CT_data[i] = 0.0;
    }

    /* 'gkmPWM:397' for j = 1:length(ct) */
    if ((ct->size[0] == 0) || (ct->size[1] == 0)) {
      m = 0;
    } else {
      u0 = ct->size[0];
      m = ct->size[1];
      if (u0 >= m) {
        m = u0;
      }
    }

    for (j = 0; j < m; j++) {
      /* 'gkmPWM:398' a = 1; */
      a = 1U;

      /* 'gkmPWM:399' for jj = 1:k_svm */
      i = (int)k_svm;
      for (jj = 0; jj < i; jj++) {
        /* 'gkmPWM:400' if ct(j,jj) ~= l_svm */
        S = ct_data[j + ct->size[0] * jj];
        if (S != l_svm) {
          /* 'gkmPWM:401' CT(j,a) = ct(j,jj); */
          CT_data[j + CT->size[0] * ((int)a - 1)] = S;

          /* 'gkmPWM:402' a = a+1; */
          a++;
        }
      }
    }

    /* 'gkmPWM:406' for j = 2:length(f) */
    i = f->size[0];
    for (j = 0; j <= i - 2; j++) {
      /* 'gkmPWM:407' a = 1; */
      /* 'gkmPWM:408' while CT(j,a)==CT(j-1,a) */
      for (a = 1U; CT_data[(j + CT->size[0] * ((int)a - 1)) + 1] == CT_data[j +
           CT->size[0] * ((int)a - 1)]; a++) {
        /* 'gkmPWM:409' a = a+1; */
      }

      /* 'gkmPWM:411' if a < 2 */
      if ((int)a < 2) {
        /* 'gkmPWM:412' a = 2; */
        a = 2U;
      }

      /* 'gkmPWM:414' diffC(f(j),i)=a; */
      diffC_data[((int)f_data[j + 1] + diffC->size[0] * acount) - 1] = a;
    }

    /* 'gkmPWM:416' diffC(f(1),i)=2; */
    diffC_data[((int)f_data[0] + diffC->size[0] * acount) - 1] = 2.0;
  }

  emxFree_boolean_T(&b_ct);

  /* 'gkmPWM:418' m = length(PWM); */
  varargin_2 = PWM->size[0] - 1;

  /* 'gkmPWM:419' scorevec = zeros(1,n); */
  i = scorevec->size[0] * scorevec->size[1];
  scorevec->size[0] = 1;
  m = (int)n;
  scorevec->size[1] = (int)n;
  emxEnsureCapacity_real_T(scorevec, i);
  scorevec_data = scorevec->data;
  for (i = 0; i < m; i++) {
    scorevec_data[i] = 0.0;
  }

  emxInit_real_T(&lenvec, 1);
  emxInit_cell_wrap_14(&loc);

  /* 'gkmPWM:420' lenvec = zeros(m,1); */
  unnamed_idx_0 = PWM->size[0];

  /* 'gkmPWM:421' loc = cell(m, 1); */
  /* 'gkmPWM:422' for i = 1:m */
  i = PWM->size[0];
  i1 = lenvec->size[0];
  lenvec->size[0] = unnamed_idx_0;
  emxEnsureCapacity_real_T(lenvec, i1);
  lenvec_data = lenvec->data;
  i1 = loc->size[0];
  loc->size[0] = PWM->size[0];
  emxEnsureCapacity_cell_wrap_14(loc, i1);
  loc_data = loc->data;
  for (acount = 0; acount < i; acount++) {
    /* 'gkmPWM:423' lenvec(i) = length(PWM{i})-l_svm*2+2; */
    u0 = PWM_data[acount].f1->size[0];
    if (u0 < 4) {
      u0 = 4;
    }

    if (PWM_data[acount].f1->size[0] == 0) {
      u0 = 0;
    }

    lenvec_data[acount] = ((double)u0 - l_svm * 2.0) + 2.0;

    /* 'gkmPWM:424' loc{i} = zeros(length(PWM{i}), 1); */
    u0 = PWM_data[acount].f1->size[0];
    if (u0 < 4) {
      u0 = 4;
    }

    if (PWM_data[acount].f1->size[0] == 0) {
      m = 0;
    } else {
      m = u0;
    }

    i1 = loc_data[acount].f1->size[0];
    loc_data[acount].f1->size[0] = m;
    emxEnsureCapacity_real_T(loc_data[acount].f1, i1);
    for (i1 = 0; i1 < m; i1++) {
      loc_data[acount].f1->data[i1] = 0.0;
    }

    /* 'gkmPWM:425' loc{i}(l_svm:lenvec(i)+l_svm-1) = 1; */
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

  /* 'gkmPWM:428' kmat = zeros(lcomb*4^k_svm, m); */
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

  /* 'gkmPWM:429' KMAT = zeros(lcomb*4^k_svm, m); */
  i = CT->size[0] * CT->size[1];
  CT->size[0] = (int)((double)diffc->size[0] * M);
  CT->size[1] = PWM->size[0];
  emxEnsureCapacity_real_T(CT, i);
  CT_data = CT->data;
  loop_ub = (int)((double)diffc->size[0] * M) * PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    CT_data[i] = 0.0;
  }

  /* 'gkmPWM:431' fprintf('Mapping PWMs to gkm space\n'); */
  printf("Mapping PWMs to gkm space\n");
  fflush(stdout);

  /* 'gkmPWM:432' if RC */
  if (RC != 0.0) {
    /* 'gkmPWM:433' for i = 1:m */
    i = PWM->size[0];
    for (acount = 0; acount < i; acount++) {
      /* 'gkmPWM:434' kmat(:,i) = PWM2kmers(PWM{i},negmat,rc,diffc,indc,loc{i},xc,l_svm,k_svm,rcnum)-negvec*(lenvec(i)+l_svm-1); */
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
        eb_binary_expand_op(kmat, acount, b_iidx, negvec, S);
        kmat_data = kmat->data;
      }

      /* map PWMs to gapped kmers */
    }
  } else {
    /* 'gkmPWM:436' else */
    /* 'gkmPWM:437' for i = 1:m */
    i = PWM->size[0];
    for (acount = 0; acount < i; acount++) {
      /* 'gkmPWM:438' kmat(:,i) = PWM2kmers_norc(PWM{i},negmat,rc,diffc,indc,loc{i},xc,l_svm,k_svm,rcnum)-negvec*(lenvec(i)+l_svm-1); */
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
        eb_binary_expand_op(kmat, acount, b_iidx, negvec, S);
        kmat_data = kmat->data;
      }

      /* map PWMs to gapped kmers */
    }
  }

  emxInit_cell_wrap_14(&poscell);

  /* the following loop creates indices for the PWM column optimize to utilize dynamic programming. */
  /* 'gkmPWM:443' poscell = cell(k_svm,1); */
  /* 'gkmPWM:444' for i = 1:k_svm */
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
    /* 'gkmPWM:445' temp = zeros(4^(k_svm-1),1)'; */
    i1 = temp->size[0] * temp->size[1];
    temp->size[0] = 1;
    m = (int)c;
    temp->size[1] = (int)c;
    emxEnsureCapacity_real_T(temp, i1);
    temp_data = temp->data;
    for (i1 = 0; i1 < m; i1++) {
      temp_data[i1] = 0.0;
    }

    /* 'gkmPWM:446' vec = 1:4^(i):4^(k_svm); */
    M = pow(4.0, (double)acount + 1.0);
    if (floor(M) == M) {
      i1 = vec->size[0] * vec->size[1];
      vec->size[0] = 1;
      loop_ub = (int)floor((scount - 1.0) / M);
      vec->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(vec, i1);
      vec_data = vec->data;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        vec_data[i1] = M * (double)i1 + 1.0;
      }
    } else {
      eml_float_colon(M, scount, vec);
      vec_data = vec->data;
    }

    /* 'gkmPWM:447' for ii = 1:4^(i-1) */
    i1 = (int)pow(4.0, ((double)acount + 1.0) - 1.0);
    if (0 <= i1 - 1) {
      b_loop_ub = vec->size[1];
    }

    for (j = 0; j < i1; j++) {
      /* 'gkmPWM:448' t = length(vec); */
      /* 'gkmPWM:449' temp(1+(ii-1)*t:t+(ii-1)*t) = vec+ii-1; */
      S = (((double)j + 1.0) - 1.0) * (double)vec->size[1] + 1.0;
      if (S > (double)vec->size[1] + (((double)j + 1.0) - 1.0) * (double)
          vec->size[1]) {
        i2 = 1;
      } else {
        i2 = (int)S;
      }

      for (i3 = 0; i3 < b_loop_ub; i3++) {
        temp_data[(i2 + i3) - 1] = (vec_data[i3] + ((double)j + 1.0)) - 1.0;
      }
    }

    /* 'gkmPWM:451' poscell{i} = sort(temp)'; */
    c_sort(temp);
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

  /* 'gkmPWM:453' fprintf('Running Recursion\n'); */
  printf("Running Recursion\n");
  fflush(stdout);

  /* 'gkmPWM:454' acount = 0; */
  acount = 0;

  /* 'gkmPWM:455' i = 0; */
  a = 0U;

  /* 'gkmPWM:456' scount = 0; */
  scount = 0.0;

  /* 'gkmPWM:457' C = (kmat'*kmat)^(-1)*(kmat'*kweig); */
  mtimes(kmat, kmat, c_ct);
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

  /* 'gkmPWM:458' ord = zeros(m, 1); */
  i = ord->size[0];
  ord->size[0] = PWM->size[0];
  emxEnsureCapacity_real_T(ord, i);
  ord_data = ord->data;
  loop_ub = PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    ord_data[i] = 0.0;
  }

  /* 'gkmPWM:459' tic */
  tic();

  /* 'gkmPWM:460' while i < n */
  emxInit_real_T(&res, 1);
  emxInit_cell_wrap_14(&new_loc);
  emxInit_cell_wrap_2(&new_PWM);
  emxInit_cell_wrap_14(&b_new_loc);
  emxInit_cell_wrap_2(&b_new_PWM);
  emxInit_int32_T(&c_iidx, 2);
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&b_mat, 2);
  emxInit_real_T(&c_PWM, 2);
  emxInit_real_T(&b_loc, 1);
  emxInit_boolean_T(&d_ct, 1);
  do {
    exitg1 = 0;
    if (a < n) {
      /* 'gkmPWM:461' i = i+1; */
      a++;

      /* 'gkmPWM:462' if mod(i,10) == 0 */
      if (fmod(a, 10.0) == 0.0) {
        /* 'gkmPWM:463' toc */
        toc();

        /* 'gkmPWM:464' tic */
        tic();

        /* 'gkmPWM:465' fprintf('%d iterations done...\n', int32(i)); */
        printf("%d iterations done...\n", (int)a);
        fflush(stdout);
      }

      /* 'gkmPWM:467' scount = scount + 1; */
      scount++;

      /* 'gkmPWM:468' if  i >= 10 && scount >= 5 && max(-1*diff(scorevec(i-5:i-1))./scorevec(i-4:i-1)) < 0.001 && acount < 5 && i ~= n */
      if ((a >= 10U) && (scount >= 5.0)) {
        b_diff(*(double (*)[5])&scorevec_data[(int)((double)a + -5.0) - 1], b);
        varargin_1[0] = -b[0] / scorevec_data[(int)((double)a + -4.0) - 1];
        varargin_1[1] = -b[1] / scorevec_data[(int)((double)a + -3.0) - 1];
        varargin_1[2] = -b[2] / scorevec_data[(int)((double)a + -2.0) - 1];
        varargin_1[3] = -b[3] / scorevec_data[(int)((double)a + -1.0) - 1];
        if ((e_maximum(varargin_1) < 0.001) && (acount < 5) && (a != n)) {
          /* 'gkmPWM:469' acount = acount + 1; */
          acount++;

          /* 'gkmPWM:470' fprintf('adjusting PWMs after %d iterations (%d)\n', int32(i), int32(acount)); */
          printf("adjusting PWMs after %d iterations (%d)\n", (int)a, acount);
          fflush(stdout);

          /* 'gkmPWM:471' scount = 0; */
          scount = 0.0;

          /* 'gkmPWM:472' for ii = 1:m */
          for (j = 0; j <= varargin_2; j++) {
            /* 'gkmPWM:473' if i/n <= 0.8 && C(ord(ii)) > 0 */
            if ((double)a / n <= 0.8) {
              i = (int)ord_data[j] - 1;
              if (C_data[i] > 0.0) {
                /* 'gkmPWM:474' [PWM{ord(ii)}, lenvec(ord(ii))] = adjust_PWM(PWM{ord(ii)}(l_svm:(length(PWM{ord(ii)})-l_svm+1),:),GC); */
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

                m = i2 - i1;
                for (i2 = 0; i2 < 4; i2++) {
                  for (i3 = 0; i3 < m; i3++) {
                    b_PWM_data[(int)ord_data[j] - 1].f1->data[i3 + m * i2] =
                      b_PWM_data[i].f1->data[((i1 + i3) + b_PWM_data[i].f1->
                      size[0] * i2) + 1];
                  }
                }

                i1 = b_PWM_data[(int)ord_data[j] - 1].f1->size[0] * b_PWM_data
                  [(int)ord_data[j] - 1].f1->size[1];
                b_PWM_data[(int)ord_data[j] - 1].f1->size[0] = m;
                b_PWM_data[(int)ord_data[j] - 1].f1->size[1] = 4;
                emxEnsureCapacity_real_T(b_PWM_data[(int)ord_data[j] - 1].f1, i1);
                varargin_1[0] = PWM_data[0].f1->data[0];
                varargin_1[1] = PWM_data[0].f1->data[PWM_data[0].f1->size[0]];
                varargin_1[2] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] * 2];
                varargin_1[3] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] * 3];
                lenvec_data[i] = adjust_PWM(b_PWM_data[(int)ord_data[j] - 1].f1,
                  varargin_1);

                /* 'gkmPWM:475' PWM{ord(ii)} = extendPWM(PWM{ord(ii)}, l_svm-1, GC); */
                /* 'gkmPWM:1002' mat = repmat(GCmat, n,1); */
                varargin_1[0] = PWM_data[0].f1->data[0];
                varargin_1[1] = PWM_data[0].f1->data[PWM_data[0].f1->size[0]];
                varargin_1[2] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] * 2];
                varargin_1[3] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] * 3];
                b_repmat(varargin_1, l_svm - 1.0, mat);
                diffC_data = mat->data;

                /* 'gkmPWM:1003' ext_pwm = [mat;pwm;mat]; */
                i1 = b_mat->size[0] * b_mat->size[1];
                b_mat->size[0] = (mat->size[0] + b_PWM_data[(int)ord_data[j] - 1]
                                  .f1->size[0]) + mat->size[0];
                b_mat->size[1] = 4;
                emxEnsureCapacity_real_T(b_mat, i1);
                vec_data = b_mat->data;
                loop_ub = mat->size[0];
                b_loop_ub = b_PWM_data[(int)ord_data[j] - 1].f1->size[0];
                for (i1 = 0; i1 < 4; i1++) {
                  for (i2 = 0; i2 < loop_ub; i2++) {
                    vec_data[i2 + b_mat->size[0] * i1] = diffC_data[i2 +
                      mat->size[0] * i1];
                  }

                  for (i2 = 0; i2 < b_loop_ub; i2++) {
                    vec_data[(i2 + mat->size[0]) + b_mat->size[0] * i1] =
                      b_PWM_data[i].f1->data[i2 + b_PWM_data[i].f1->size[0] * i1];
                  }
                }

                loop_ub = mat->size[0];
                for (i1 = 0; i1 < 4; i1++) {
                  for (i2 = 0; i2 < loop_ub; i2++) {
                    vec_data[((i2 + mat->size[0]) + b_PWM_data[(int)ord_data[j]
                              - 1].f1->size[0]) + b_mat->size[0] * i1] =
                      diffC_data[i2 + mat->size[0] * i1];
                  }
                }

                i1 = b_PWM_data[(int)ord_data[j] - 1].f1->size[0] * b_PWM_data
                  [(int)ord_data[j] - 1].f1->size[1];
                b_PWM_data[(int)ord_data[j] - 1].f1->size[0] = b_mat->size[0];
                b_PWM_data[(int)ord_data[j] - 1].f1->size[1] = 4;
                emxEnsureCapacity_real_T(b_PWM_data[(int)ord_data[j] - 1].f1, i1);
                loop_ub = b_mat->size[0] * 4;
                for (i1 = 0; i1 < loop_ub; i1++) {
                  b_PWM_data[(int)ord_data[j] - 1].f1->data[i1] = vec_data[i1];
                }

                /* 'gkmPWM:476' loc{ord(ii)} = zeros(lenvec(ord(ii))+2*l_svm-2, 1); */
                S = lenvec_data[i];
                loop_ub = (int)((S + 2.0 * l_svm) - 2.0);
                i = loc_data[(int)ord_data[j] - 1].f1->size[0];
                loc_data[(int)ord_data[j] - 1].f1->size[0] = loop_ub;
                emxEnsureCapacity_real_T(loc_data[(int)ord_data[j] - 1].f1, i);
                for (i = 0; i < loop_ub; i++) {
                  loc_data[(int)ord_data[j] - 1].f1->data[i] = 0.0;
                }

                /* 'gkmPWM:477' loc{ord(ii)}(l_svm:lenvec(ord(ii))+l_svm-1) = 1; */
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

      /* 'gkmPWM:481' C = (kmat'*kmat)^(-1)*(kmat'*kweig); */
      mtimes(kmat, kmat, c_ct);
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

      /* 'gkmPWM:482' res = kweig-kmat*C; */
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

      /* 'gkmPWM:483' corrvec = zeros(m,1); */
      /* 'gkmPWM:484' for ii = 1:m */
      i = ord->size[0];
      ord->size[0] = unnamed_idx_0;
      emxEnsureCapacity_real_T(ord, i);
      ord_data = ord->data;
      for (j = 0; j <= varargin_2; j++) {
        /* 'gkmPWM:485' [~,ind] = sort(kmat(:,ii), 'descend'); */
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

        /* 'gkmPWM:486' corrvec(ii) = sum(res(ind(1:lcomb)).^2); */
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

        /* 'gkmPWM:487' if mod(i,20) == 0 */
        if (fmod(a, 20.0) == 0.0) {
          /* 'gkmPWM:488' if RC */
          if (RC != 0.0) {
            /* 'gkmPWM:489' kmat(:,ii) = PWM2kmers(PWM{ii},negmat,rc,diffc,indc,loc{ii},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(ii)); */
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
              eb_binary_expand_op(kmat, j, b_iidx, negvec, S);
              kmat_data = kmat->data;
            }
          } else {
            /* 'gkmPWM:490' else */
            /* 'gkmPWM:491' kmat(:,ii) = PWM2kmers_norc(PWM{ii},negmat,rc,diffc,indc,loc{ii},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(ii)); */
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
              eb_binary_expand_op(kmat, j, b_iidx, negvec, S);
              kmat_data = kmat->data;
            }
          }
        }
      }

      /* The order of PWM optimization is determined by the correlation of its top 110 kmers with the gapped kmer weight vector */
      /* 'gkmPWM:496' [~,ord] = sort(corrvec, 'descend'); */
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

      /* 'gkmPWM:497' for ii = 1:m */
      for (j = 0; j <= varargin_2; j++) {
        /* The order of the column optimization is determined by the max probability in each column */
        /* 'gkmPWM:499' v = max(PWM{ord(ii)}(l_svm:lenvec(ord(ii))+l_svm-1,:)'); */
        i = (int)ord_data[j] - 1;
        S = (lenvec_data[i] + l_svm) - 1.0;
        if (l_svm > S) {
          i1 = 0;
          i2 = 0;
        } else {
          i1 = (int)l_svm - 1;
          i2 = (int)S;
        }

        /* 'gkmPWM:500' [~,c] = sort(v, 'ascend'); */
        i3 = c_PWM->size[0] * c_PWM->size[1];
        c_PWM->size[0] = 4;
        loop_ub = i2 - i1;
        c_PWM->size[1] = loop_ub;
        emxEnsureCapacity_real_T(c_PWM, i3);
        vec_data = c_PWM->data;
        for (i2 = 0; i2 < loop_ub; i2++) {
          m = i1 + i2;
          vec_data[4 * i2] = b_PWM_data[i].f1->data[m];
          vec_data[4 * i2 + 1] = b_PWM_data[i].f1->data[m + b_PWM_data[i]
            .f1->size[0]];
          vec_data[4 * i2 + 2] = b_PWM_data[i].f1->data[m + b_PWM_data[i]
            .f1->size[0] * 2];
          vec_data[4 * i2 + 3] = b_PWM_data[i].f1->data[m + b_PWM_data[i]
            .f1->size[0] * 3];
        }

        f_maximum(c_PWM, temp);
        d_sort(temp, c_iidx);
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
        /* 'gkmPWM:503' for iii = 1:length(c) */
        i1 = temp->size[1];
        for (jj = 0; jj < i1; jj++) {
          /* 'gkmPWM:504' PWMtemp = PWM{ord(ii)}(c(iii):c(iii)+l_svm*2-2,:); */
          S = temp_data[jj];
          M = ((double)(int)S + l_svm * 2.0) - 2.0;
          if ((int)S > M) {
            i2 = 0;
            i3 = 0;
            m = 0;
            u0 = 0;
          } else {
            i2 = (int)S - 1;
            i3 = (int)M;
            m = (int)S - 1;
            u0 = (int)M;
          }

          /* 'gkmPWM:505' [kweigdiff,PWM{ord(ii)}(c(iii)+l_svm-1,:)] = getEMprob_v3(PWMtemp,res/C(ord(ii)),negmat,poscell,rc,diffC,indc,loc{ord(ii)}(c(iii):c(iii)+2*l_svm-2),xc,reg,l_svm,k_svm,rcnum,RC); */
          loop_ub = i3 - i2;
          i3 = mat->size[0] * mat->size[1];
          mat->size[0] = loop_ub;
          mat->size[1] = 4;
          emxEnsureCapacity_real_T(mat, i3);
          diffC_data = mat->data;
          for (i3 = 0; i3 < 4; i3++) {
            for (b_loop_ub = 0; b_loop_ub < loop_ub; b_loop_ub++) {
              diffC_data[b_loop_ub + mat->size[0] * i3] = b_PWM_data[i].f1->
                data[(i2 + b_loop_ub) + b_PWM_data[i].f1->size[0] * i3];
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

          loop_ub = u0 - m;
          i2 = b_loc->size[0];
          b_loc->size[0] = loop_ub;
          emxEnsureCapacity_real_T(b_loc, i2);
          vec_data = b_loc->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            vec_data[i2] = loc_data[i].f1->data[m + i2];
          }

          getEMprob_v3(mat, b_iidx, negmat, poscell, rc, diffC, indc, b_loc, xc,
                       reg, l_svm, k_svm, rcnum, RC, f, b);
          f_data = f->data;
          m = (int)(((double)(int)temp_data[jj] + l_svm) - 1.0);
          b_PWM_data[(int)ord_data[j] - 1].f1->data[m - 1] = b[0];
          b_PWM_data[(int)ord_data[j] - 1].f1->data[(m + b_PWM_data[(int)
            ord_data[j] - 1].f1->size[0]) - 1] = b[1];
          b_PWM_data[(int)ord_data[j] - 1].f1->data[(m + b_PWM_data[(int)
            ord_data[j] - 1].f1->size[0] * 2) - 1] = b[2];
          b_PWM_data[(int)ord_data[j] - 1].f1->data[(m + b_PWM_data[(int)
            ord_data[j] - 1].f1->size[0] * 3) - 1] = b[3];

          /* 'gkmPWM:506' kmat(:,ord(ii)) = kmat(:,ord(ii)) + kweigdiff; */
          if (kmat->size[0] == f->size[0]) {
            m = kmat->size[0] - 1;
            S = ord_data[j];
            i2 = b_iidx->size[0];
            b_iidx->size[0] = kmat->size[0];
            emxEnsureCapacity_real_T(b_iidx, i2);
            b_iidx_data = b_iidx->data;
            for (i2 = 0; i2 <= m; i2++) {
              b_iidx_data[i2] = kmat_data[i2 + kmat->size[0] * ((int)S - 1)] +
                f_data[i2];
            }

            loop_ub = b_iidx->size[0];
            for (i2 = 0; i2 < loop_ub; i2++) {
              kmat_data[i2 + kmat->size[0] * ((int)S - 1)] = b_iidx_data[i2];
            }
          } else {
            gb_binary_expand_op(kmat, ord, j, f);
            kmat_data = kmat->data;
          }

          /* 'gkmPWM:507' res = res-kweigdiff*C(ord(ii)); */
          loop_ub = res->size[0];
          if (res->size[0] == f->size[0]) {
            for (i2 = 0; i2 < loop_ub; i2++) {
              res_data[i2] -= f_data[i2] * C_data[i];
            }
          } else {
            fb_binary_expand_op(res, f, C, ord, j);
            res_data = res->data;
          }
        }
      }

      /* Reseed PWMs if two or more of them are too highly correlated */
      /* 'gkmPWM:511' if i/n <= 0.80 */
      if ((double)a / n <= 0.8) {
        /* 'gkmPWM:512' info = avg_info(PWM,l_svm); */
        avg_info(b_PWM, l_svm, f);

        /* 'gkmPWM:513' [~,ord] = sort(info,'descend'); */
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

        /* 'gkmPWM:514' kmat = kmat(:,ord); */
        m = kmat->size[0] - 1;
        i = c_ct->size[0] * c_ct->size[1];
        c_ct->size[0] = kmat->size[0];
        c_ct->size[1] = ord->size[0];
        emxEnsureCapacity_real_T(c_ct, i);
        temp_data = c_ct->data;
        loop_ub = ord->size[0];
        for (i = 0; i < loop_ub; i++) {
          for (i1 = 0; i1 <= m; i1++) {
            temp_data[i1 + c_ct->size[0] * i] = kmat_data[i1 + kmat->size[0] *
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
          kmat_data[i] = temp_data[i];
        }

        /* 'gkmPWM:515' lenvec = lenvec(ord); */
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
        /* 'gkmPWM:517' ord_len = length(ord); */
        /* 'gkmPWM:518' new_loc = cell(ord_len,1); */
        m = ord->size[0];
        i = new_loc->size[0];
        new_loc->size[0] = ord->size[0];
        emxEnsureCapacity_cell_wrap_14(new_loc, i);
        poscell_data = new_loc->data;
        for (i = 0; i < m; i++) {
          poscell_data[i].f1->size[0] = 0;
        }

        /* 'gkmPWM:519' new_loc = coder.nullcopy(new_loc); */
        i = b_new_loc->size[0];
        b_new_loc->size[0] = new_loc->size[0];
        emxEnsureCapacity_cell_wrap_14(b_new_loc, i);
        poscell_data = b_new_loc->data;

        /* 'gkmPWM:520' for cur_idx=1:ord_len */
        i = ord->size[0];
        for (u0 = 0; u0 < i; u0++) {
          /* 'gkmPWM:521' new_loc{cur_idx} = loc{ord(cur_idx)}; */
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

        /* 'gkmPWM:523' loc = coder.nullcopy(loc); */
        /* 'gkmPWM:524' for cur_idx=1:ord_len */
        i = ord->size[0];
        for (u0 = 0; u0 < i; u0++) {
          /* 'gkmPWM:525' loc{cur_idx} = new_loc{cur_idx}; */
          loop_ub = poscell_data[u0].f1->size[0];
          i1 = loc_data[u0].f1->size[0];
          loc_data[u0].f1->size[0] = poscell_data[u0].f1->size[0];
          emxEnsureCapacity_real_T(loc_data[u0].f1, i1);
          for (i1 = 0; i1 < loop_ub; i1++) {
            loc_data[u0].f1->data[i1] = poscell_data[u0].f1->data[i1];
          }
        }

        /*  PWM = PWM(ord); */
        /* 'gkmPWM:530' new_PWM = cell(ord_len,1); */
        i = new_PWM->size[0];
        new_PWM->size[0] = ord->size[0];
        emxEnsureCapacity_cell_wrap_2(new_PWM, i);
        new_PWM_data = new_PWM->data;
        for (i = 0; i < m; i++) {
          new_PWM_data[i].f1->size[0] = 0;
          new_PWM_data[i].f1->size[1] = 4;
        }

        /* 'gkmPWM:531' new_PWM = coder.nullcopy(new_PWM); */
        i = b_new_PWM->size[0];
        b_new_PWM->size[0] = new_PWM->size[0];
        emxEnsureCapacity_cell_wrap_2(b_new_PWM, i);
        new_PWM_data = b_new_PWM->data;

        /* 'gkmPWM:532' for cur_idx=1:ord_len */
        i = ord->size[0];
        for (u0 = 0; u0 < i; u0++) {
          /* 'gkmPWM:533' new_PWM{cur_idx} = PWM{ord(cur_idx)}; */
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

        /* 'gkmPWM:535' PWM = coder.nullcopy(PWM); */
        /* 'gkmPWM:536' for cur_idx=1:ord_len */
        i = ord->size[0];
        for (u0 = 0; u0 < i; u0++) {
          /* 'gkmPWM:537' PWM{cur_idx} = new_PWM{cur_idx}; */
          i1 = b_PWM_data[u0].f1->size[0] * b_PWM_data[u0].f1->size[1];
          b_PWM_data[u0].f1->size[0] = new_PWM_data[u0].f1->size[0];
          b_PWM_data[u0].f1->size[1] = 4;
          emxEnsureCapacity_real_T(b_PWM_data[u0].f1, i1);
          loop_ub = new_PWM_data[u0].f1->size[0] * 4;
          for (i1 = 0; i1 < loop_ub; i1++) {
            b_PWM_data[u0].f1->data[i1] = new_PWM_data[u0].f1->data[i1];
          }
        }

        /* 'gkmPWM:540' C = C(ord); */
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

        /* 'gkmPWM:541' for j = 1:m */
        for (j = 0; j <= varargin_2; j++) {
          /* 'gkmPWM:542' KMAT(:,j) = kmat(:,j)/(kmat(:,j)'*kmat(:,j)); */
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
            M = 0.0;
          } else {
            M = cblas_ddot((blasint)kmat->size[0], &f_data[0], (blasint)1,
                           &b_iidx_data[0], (blasint)1);
          }

          for (i = 0; i < loop_ub; i++) {
            CT_data[i + CT->size[0] * j] = kmat_data[i + kmat->size[0] * j] / M;
          }
        }

        /* 'gkmPWM:544' MAT = KMAT'*KMAT; */
        mtimes(CT, CT, ct);
        ct_data = ct->data;

        /* 'gkmPWM:545' for j = 1:m-1 */
        i = PWM->size[0];
        for (j = 0; j <= i - 2; j++) {
          /* 'gkmPWM:546' vec = MAT(j+1:end,j); */
          if (j + 2U > (unsigned int)ct->size[0]) {
            i1 = 0;
            i2 = 0;
          } else {
            i1 = j + 1;
            i2 = ct->size[0];
          }

          /* 'gkmPWM:547' [a b] = max(vec); */
          loop_ub = i2 - i1;
          i2 = b_iidx->size[0];
          b_iidx->size[0] = loop_ub;
          emxEnsureCapacity_real_T(b_iidx, i2);
          b_iidx_data = b_iidx->data;
          for (i2 = 0; i2 < loop_ub; i2++) {
            b_iidx_data[i2] = ct_data[(i1 + i2) + ct->size[0] * j];
          }

          d_maximum(b_iidx, &M, &m);

          /* 'gkmPWM:548' if a > rcorr && C(j) > 0 */
          if ((M > rcorr) && (C_data[j] > 0.0)) {
            /* 'gkmPWM:549' scount = 0; */
            scount = 0.0;

            /* 'gkmPWM:550' fprintf('reseeding\n'); */
            printf("reseeding\n");
            fflush(stdout);

            /* 'gkmPWM:551' f = j+find(vec > rcorr); */
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

            /* 'gkmPWM:552' for jj = 1:length(f) */
            i1 = f->size[0];
            if (0 <= f->size[0] - 1) {
              b_lenvec[1] = 12.0;
            }

            for (jj = 0; jj < i1; jj++) {
              /* 'gkmPWM:553' for jjj = 1:min([lenvec(f(jj)) 12]) */
              u0 = (int)f_data[jj] - 1;
              M = lenvec_data[u0];
              b_lenvec[0] = M;
              i2 = (int)b_minimum(b_lenvec);
              for (b_loop_ub = 0; b_loop_ub < i2; b_loop_ub++) {
                /* 'gkmPWM:554' PWM{f(jj)}(jj+9,:) =  PWM{f(jj)}(jjj+l_svm-1,randperm(4)); */
                randperm(b);
                m = (int)((((double)b_loop_ub + 1.0) + l_svm) - 1.0);
                varargin_1[1] = b_PWM_data[u0].f1->data[(m + b_PWM_data[u0]
                  .f1->size[0] * ((int)b[1] - 1)) - 1];
                varargin_1[2] = b_PWM_data[u0].f1->data[(m + b_PWM_data[u0]
                  .f1->size[0] * ((int)b[2] - 1)) - 1];
                varargin_1[3] = b_PWM_data[u0].f1->data[(m + b_PWM_data[u0]
                  .f1->size[0] * ((int)b[3] - 1)) - 1];
                b_PWM_data[(int)f_data[jj] - 1].f1->data[jj + 9] = b_PWM_data[u0]
                  .f1->data[(m + b_PWM_data[u0].f1->size[0] * ((int)b[0] - 1)) -
                  1];
                b_PWM_data[(int)f_data[jj] - 1].f1->data[(jj + b_PWM_data[(int)
                  f_data[jj] - 1].f1->size[0]) + 9] = varargin_1[1];
                b_PWM_data[(int)f_data[jj] - 1].f1->data[(jj + b_PWM_data[(int)
                  f_data[jj] - 1].f1->size[0] * 2) + 9] = varargin_1[2];
                b_PWM_data[(int)f_data[jj] - 1].f1->data[(jj + b_PWM_data[(int)
                  f_data[jj] - 1].f1->size[0] * 3) + 9] = varargin_1[3];
              }

              /* 'gkmPWM:556' if lenvec(f(jj)) >= 12 */
              if (M >= 12.0) {
                /* 'gkmPWM:557' PWM{f(jj)} = PWM{f(jj)}(l_svm:l_svm+11,:); */
                i2 = mat->size[0] * mat->size[1];
                mat->size[0] = 12;
                mat->size[1] = 4;
                emxEnsureCapacity_real_T(mat, i2);
                diffC_data = mat->data;
                for (i2 = 0; i2 < 4; i2++) {
                  for (i3 = 0; i3 < 12; i3++) {
                    diffC_data[i3 + mat->size[0] * i2] = b_PWM_data[u0].f1->
                      data[((int)(l_svm + (double)i3) + b_PWM_data[u0].f1->size
                            [0] * i2) - 1];
                  }
                }

                i2 = b_PWM_data[(int)f_data[jj] - 1].f1->size[0] * b_PWM_data
                  [(int)f_data[jj] - 1].f1->size[1];
                b_PWM_data[(int)f_data[jj] - 1].f1->size[0] = 12;
                b_PWM_data[(int)f_data[jj] - 1].f1->size[1] = 4;
                emxEnsureCapacity_real_T(b_PWM_data[(int)f_data[jj] - 1].f1, i2);
                for (i2 = 0; i2 < 48; i2++) {
                  b_PWM_data[(int)f_data[jj] - 1].f1->data[i2] = diffC_data[i2];
                }

                /* 'gkmPWM:558' PWM{f(jj)} = extendPWM(PWM{f(jj)}, l_svm-1, GC); */
                /* 'gkmPWM:1002' mat = repmat(GCmat, n,1); */
                varargin_1[0] = PWM_data[0].f1->data[0];
                varargin_1[1] = PWM_data[0].f1->data[PWM_data[0].f1->size[0]];
                varargin_1[2] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] * 2];
                varargin_1[3] = PWM_data[0].f1->data[PWM_data[0].f1->size[0] * 3];
                b_repmat(varargin_1, l_svm - 1.0, mat);
                diffC_data = mat->data;

                /* 'gkmPWM:1003' ext_pwm = [mat;pwm;mat]; */
                i2 = b_mat->size[0] * b_mat->size[1];
                b_mat->size[0] = (mat->size[0] + b_PWM_data[(int)f_data[jj] - 1]
                                  .f1->size[0]) + mat->size[0];
                b_mat->size[1] = 4;
                emxEnsureCapacity_real_T(b_mat, i2);
                vec_data = b_mat->data;
                loop_ub = mat->size[0];
                b_loop_ub = b_PWM_data[(int)f_data[jj] - 1].f1->size[0];
                for (i2 = 0; i2 < 4; i2++) {
                  for (i3 = 0; i3 < loop_ub; i3++) {
                    vec_data[i3 + b_mat->size[0] * i2] = diffC_data[i3 +
                      mat->size[0] * i2];
                  }

                  for (i3 = 0; i3 < b_loop_ub; i3++) {
                    vec_data[(i3 + mat->size[0]) + b_mat->size[0] * i2] =
                      b_PWM_data[u0].f1->data[i3 + b_PWM_data[u0].f1->size[0] *
                      i2];
                  }
                }

                loop_ub = mat->size[0];
                for (i2 = 0; i2 < 4; i2++) {
                  for (i3 = 0; i3 < loop_ub; i3++) {
                    vec_data[((i3 + mat->size[0]) + b_PWM_data[(int)f_data[jj] -
                              1].f1->size[0]) + b_mat->size[0] * i2] =
                      diffC_data[i3 + mat->size[0] * i2];
                  }
                }

                i2 = b_PWM_data[(int)f_data[jj] - 1].f1->size[0] * b_PWM_data
                  [(int)f_data[jj] - 1].f1->size[1];
                b_PWM_data[(int)f_data[jj] - 1].f1->size[0] = b_mat->size[0];
                b_PWM_data[(int)f_data[jj] - 1].f1->size[1] = 4;
                emxEnsureCapacity_real_T(b_PWM_data[(int)f_data[jj] - 1].f1, i2);
                loop_ub = b_mat->size[0] * 4;
                for (i2 = 0; i2 < loop_ub; i2++) {
                  b_PWM_data[(int)f_data[jj] - 1].f1->data[i2] = vec_data[i2];
                }

                /* 'gkmPWM:559' lenvec(f(jj)) = 12; */
                lenvec_data[u0] = 12.0;

                /* 'gkmPWM:560' loc{f(jj)} = zeros(lenvec(f(jj))+l_svm*2-2, 1); */
                loop_ub = (int)((lenvec_data[u0] + l_svm * 2.0) - 2.0);
                i2 = loc_data[(int)f_data[jj] - 1].f1->size[0];
                loc_data[(int)f_data[jj] - 1].f1->size[0] = loop_ub;
                emxEnsureCapacity_real_T(loc_data[(int)f_data[jj] - 1].f1, i2);
                for (i2 = 0; i2 < loop_ub; i2++) {
                  loc_data[(int)f_data[jj] - 1].f1->data[i2] = 0.0;
                }

                /* 'gkmPWM:561' loc{f(jj)}(l_svm:lenvec(f(jj))+l_svm-1) = 1; */
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

              /* 'gkmPWM:563' if RC */
              if (RC != 0.0) {
                /* 'gkmPWM:564' kmat(:,f(jj)) = PWM2kmers(PWM{f(jj)},negmat,rc,diffc,indc,loc{f(jj)},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(f(jj))); */
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
                  hb_binary_expand_op(kmat, f, jj, b_iidx, negvec, S);
                  kmat_data = kmat->data;
                }
              } else {
                /* 'gkmPWM:565' else */
                /* 'gkmPWM:566' kmat(:,f(jj)) = PWM2kmers_norc(PWM{f(jj)},negmat,rc,diffc,indc,loc{f(jj)},xc,l_svm,k_svm,rcnum)-negvec*(l_svm-1+lenvec(f(jj))); */
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
                  hb_binary_expand_op(kmat, f, jj, b_iidx, negvec, S);
                  kmat_data = kmat->data;
                }
              }
            }
          }
        }
      }

      /* Breaks the loop if it looks like it converged */
      /* 'gkmPWM:573' C = (kmat'*kmat)^(-1)*(kmat'*kweig); */
      mtimes(kmat, kmat, c_ct);
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

      /* 'gkmPWM:574' scorevec(i) = sqrt(res'*res); */
      if (res->size[0] < 1) {
        c = 0.0;
      } else {
        c = cblas_ddot((blasint)res->size[0], &res_data[0], (blasint)1,
                       &res_data[0], (blasint)1);
      }

      scorevec_data[(int)a - 1] = sqrt(c);

      /* 'gkmPWM:575' if i >= 10 && acount == 5 && scount >= 10 && max(abs(diff(scorevec(i-9:i))./scorevec(i-8:i))) < 0.0001 */
      guard1 = false;
      if (((int)a >= 10) && (acount == 5) && (scount >= 10.0)) {
        c_diff(*(double (*)[10])&scorevec_data[(int)a - 10], b_varargin_1);
        for (u0 = 0; u0 < 9; u0++) {
          b_varargin_1[u0] = fabs(b_varargin_1[u0] / scorevec_data[(u0 + (int)a)
            - 9]);
        }

        if (j_maximum(b_varargin_1) < 0.0001) {
          /* 'gkmPWM:576' scorevec = scorevec(1:i); */
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
        /* 'gkmPWM:579' if i > 10 && acount == 5 && scount >= 10 && sum(diff(scorevec(i-9:i))>0) > 7 */
        c_diff(*(double (*)[10])&scorevec_data[(int)a - 10], b_varargin_1);
        for (i = 0; i < 9; i++) {
          x[i] = (b_varargin_1[i] > 0.0);
        }

        m = x[0];
        for (u0 = 0; u0 < 8; u0++) {
          m += x[u0 + 1];
        }

        if (m > 7) {
          /* 'gkmPWM:580' scorevec = scorevec(1:i); */
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
  emxFree_real_T(&c_PWM);
  emxFree_real_T(&b_mat);
  emxFree_real_T(&mat);
  emxFree_int32_T(&c_iidx);
  emxFree_cell_wrap_2(&b_new_PWM);
  emxFree_cell_wrap_14(&b_new_loc);
  emxFree_real_T(&xc);
  emxFree_real_T(&indc);
  emxFree_real_T(&rc);
  emxFree_cell_wrap_2(&new_PWM);
  emxFree_cell_wrap_14(&new_loc);
  emxFree_real_T(&ord);
  emxFree_cell_wrap_14(&poscell);
  emxFree_cell_wrap_14(&loc);
  emxFree_real_T(&lenvec);
  emxFree_real_T(&diffC);

  /* 'gkmPWM:584' toc */
  toc();

  /* 'gkmPWM:585' fprintf('gkmPWM completed after %d iterations\n', int32(i)); */
  printf("gkmPWM completed after %d iterations\n", (int)a);
  fflush(stdout);

  /* 'gkmPWM:587' for i = 1:length(PWM) */
  i = b_PWM->size[0];
  for (acount = 0; acount < i; acount++) {
    /* 'gkmPWM:588' PWM{i} = PWM{i}(l_svm:(length(PWM{i})-l_svm+1),:); */
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

    m = i2 - i1;
    for (i2 = 0; i2 < 4; i2++) {
      for (i3 = 0; i3 < m; i3++) {
        b_PWM_data[acount].f1->data[i3 + m * i2] = b_PWM_data[acount].f1->data
          [((i1 + i3) + b_PWM_data[acount].f1->size[0] * i2) + 1];
      }
    }

    i1 = b_PWM_data[acount].f1->size[0] * b_PWM_data[acount].f1->size[1];
    b_PWM_data[acount].f1->size[0] = m;
    b_PWM_data[acount].f1->size[1] = 4;
    emxEnsureCapacity_real_T(b_PWM_data[acount].f1, i1);
  }

  /* the following just calculates a few interesting quantities */
  /*  r = corr(kweig, kmat*C); */
  /* 'gkmPWM:592' r = corrcoef(kweig, kmat*C); */
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

  /* 'gkmPWM:593' r = r(1,2); */
  corrcoef(kweig, b_iidx, b);

  /* 'gkmPWM:594' M = mean(kweig); */
  M = blockedSummation(kweig, kweig->size[0]) / (double)kweig->size[0];

  /* 'gkmPWM:595' S = std(kweig); */
  S = b_std(kweig);

  /* 'gkmPWM:596' R = zeros(m,1); */
  /* 'gkmPWM:597' E = zeros(m,1); */
  /* 'gkmPWM:598' CM = corrcoef(kmat)-eye(m); */
  m = PWM->size[0];
  i = ct->size[0] * ct->size[1];
  ct->size[0] = PWM->size[0];
  ct->size[1] = PWM->size[0];
  emxEnsureCapacity_real_T(ct, i);
  ct_data = ct->data;
  loop_ub = PWM->size[0] * PWM->size[0];
  for (i = 0; i < loop_ub; i++) {
    ct_data[i] = 0.0;
  }

  for (u0 = 0; u0 < m; u0++) {
    ct_data[u0 + ct->size[0] * u0] = 1.0;
  }

  b_corrcoef(kmat, CT);
  CT_data = CT->data;
  if ((CT->size[0] == ct->size[0]) && (CT->size[1] == ct->size[1])) {
    loop_ub = CT->size[0] * CT->size[1];
    for (i = 0; i < loop_ub; i++) {
      CT_data[i] -= ct_data[i];
    }
  } else {
    b_minus(CT, ct);
    CT_data = CT->data;
  }

  /* 'gkmPWM:599' for i = 1:m */
  i = PWM->size[0];
  i1 = R->size[0];
  R->size[0] = unnamed_idx_0;
  emxEnsureCapacity_real_T(R, i1);
  vec_data = R->data;
  i1 = E->size[0];
  E->size[0] = unnamed_idx_0;
  emxEnsureCapacity_real_T(E, i1);
  diffC_data = E->data;
  if (0 <= i - 1) {
    i4 = kmat->size[0];
    c_loop_ub = kmat->size[0];
    if (1 > diffc->size[0]) {
      d_loop_ub = 0;
    } else {
      d_loop_ub = diffc->size[0];
    }

    e_loop_ub = kmat->size[0] * kmat->size[1];
  }

  emxFree_real_T(&diffc);
  for (acount = 0; acount < i; acount++) {
    /* 'gkmPWM:600' [~,a] = sort(kmat(:,i),'descend'); */
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

    /* 'gkmPWM:601' R(i) = (mean(kweig(a(1:lcomb)))-M)/S; */
    i1 = b_iidx->size[0];
    b_iidx->size[0] = d_loop_ub;
    emxEnsureCapacity_real_T(b_iidx, i1);
    b_iidx_data = b_iidx->data;
    for (i1 = 0; i1 < d_loop_ub; i1++) {
      b_iidx_data[i1] = kweig_data[(int)f_data[i1] - 1];
    }

    vec_data[acount] = (blockedSummation(b_iidx, d_loop_ub) / (double)d_loop_ub
                        - M) / S;

    /* 'gkmPWM:602' Kmat = kmat; */
    /* 'gkmPWM:603' Kmat(:,i) = []; */
    i1 = ct->size[0] * ct->size[1];
    ct->size[0] = kmat->size[0];
    ct->size[1] = kmat->size[1];
    emxEnsureCapacity_real_T(ct, i1);
    ct_data = ct->data;
    for (i1 = 0; i1 < e_loop_ub; i1++) {
      ct_data[i1] = kmat_data[i1];
    }

    nullAssignment(ct, acount + 1);
    ct_data = ct->data;

    /* 'gkmPWM:604' c = (Kmat'*Kmat)^(-1)*(Kmat'*kweig); */
    mtimes(ct, ct, c_ct);
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

    /* 'gkmPWM:605' res = kweig-Kmat*c; */
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

    /* 'gkmPWM:606' E(i) = (sqrt(res'*res)-scorevec(end))/scorevec(end); */
    if (res->size[0] < 1) {
      c = 0.0;
    } else {
      c = cblas_ddot((blasint)res->size[0], &res_data[0], (blasint)1, &res_data
                     [0], (blasint)1);
    }

    diffC_data[acount] = (sqrt(c) - scorevec_data[scorevec->size[1] - 1]) /
      scorevec_data[scorevec->size[1] - 1];

    /* 'gkmPWM:607' if C(i) < 0 */
    if (C_data[acount] < 0.0) {
      /* 'gkmPWM:608' CM(i,:) = 0; */
      loop_ub = CT->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        CT_data[acount + CT->size[0] * i1] = 0.0;
      }

      /* 'gkmPWM:609' CM(:,i) = 0; */
      loop_ub = CT->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        CT_data[i1 + CT->size[0] * acount] = 0.0;
      }
    }
  }

  emxFree_real_T(&c_ct);
  emxFree_real_T(&b_a);
  emxFree_real_T(&res);
  emxFree_real_T(&kmat);
  emxFree_real_T(&ct);

  /* 'gkmPWM:612' [R,a] = sort(R, 'descend'); */
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
  /* 'gkmPWM:614' a_len = length(a); */
  /* 'gkmPWM:615' new_PWM = cell(a_len,1); */
  unnamed_idx_0 = f->size[0];
  i = PWM->size[0];
  PWM->size[0] = f->size[0];
  emxEnsureCapacity_cell_wrap_2(PWM, i);
  PWM_data = PWM->data;
  for (i = 0; i < unnamed_idx_0; i++) {
    PWM_data[i].f1->size[0] = 0;
    PWM_data[i].f1->size[1] = 4;
  }

  /* 'gkmPWM:616' for cur_idx=1:a_len */
  i = f->size[0];
  i1 = PWM->size[0];
  PWM->size[0] = f->size[0];
  emxEnsureCapacity_cell_wrap_2(PWM, i1);
  PWM_data = PWM->data;
  for (u0 = 0; u0 < i; u0++) {
    /* 'gkmPWM:617' new_PWM{cur_idx} = PWM{a(cur_idx)}; */
    i1 = PWM_data[u0].f1->size[0] * PWM_data[u0].f1->size[1];
    PWM_data[u0].f1->size[0] = b_PWM_data[(int)f_data[u0] - 1].f1->size[0];
    PWM_data[u0].f1->size[1] = 4;
    emxEnsureCapacity_real_T(PWM_data[u0].f1, i1);
    loop_ub = b_PWM_data[(int)f_data[u0] - 1].f1->size[0] * 4;
    for (i1 = 0; i1 < loop_ub; i1++) {
      PWM_data[u0].f1->data[i1] = b_PWM_data[(int)f_data[u0] - 1].f1->data[i1];
    }
  }

  emxFree_cell_wrap_2(&b_PWM);

  /* 'gkmPWM:619' PWM = new_PWM; */
  /* 'gkmPWM:621' C = C(a); */
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

  /* 'gkmPWM:622' E = E(a); */
  i = b_iidx->size[0];
  b_iidx->size[0] = f->size[0];
  emxEnsureCapacity_real_T(b_iidx, i);
  b_iidx_data = b_iidx->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_iidx_data[i] = diffC_data[(int)f_data[i] - 1];
  }

  i = E->size[0];
  E->size[0] = b_iidx->size[0];
  emxEnsureCapacity_real_T(E, i);
  diffC_data = E->data;
  loop_ub = b_iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    diffC_data[i] = b_iidx_data[i];
  }

  emxFree_real_T(&b_iidx);

  /* 'gkmPWM:623' Rd = max(CM); */
  g_maximum(CT, Rd);
  vec_data = Rd->data;

  /* 'gkmPWM:624' Rd = Rd(a); */
  i = temp->size[0] * temp->size[1];
  temp->size[0] = 1;
  temp->size[1] = f->size[0];
  emxEnsureCapacity_real_T(temp, i);
  temp_data = temp->data;
  loop_ub = f->size[0];
  emxFree_real_T(&CT);
  for (i = 0; i < loop_ub; i++) {
    temp_data[i] = vec_data[(int)f_data[i] - 1];
  }

  emxFree_real_T(&f);
  i = Rd->size[0] * Rd->size[1];
  Rd->size[0] = 1;
  Rd->size[1] = temp->size[1];
  emxEnsureCapacity_real_T(Rd, i);
  vec_data = Rd->data;
  loop_ub = temp->size[1];
  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = temp_data[i];
  }

  emxFree_real_T(&temp);
  *r = b[2];
}

static void hb_binary_expand_op(emxArray_real_T *kmat, const emxArray_real_T *f,
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

static void i_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
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

static void j_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
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

static void k_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
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

static void l_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
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
  emxArray_real_T *b_kweig;
  emxArray_real_T *indloc2;
  emxArray_real_T *indvec;
  emxArray_real_T *indvec2_loop2;
  emxArray_real_T *indvec_loop2;
  emxArray_real_T *mat2;
  emxArray_real_T *sPWM;
  emxArray_real_T *sPWM2;
  emxArray_uint32_T *X;
  double b_p[16];
  double b_sPWM2[4];
  const double *c_data;
  const double *ind_data;
  const double *indloc_data;
  const double *mat_data;
  const double *s_data;
  const double *x_data;
  double b_c;
  double b_n_tmp;
  double c_tmp;
  double d;
  double n;
  double n_tmp;
  double *a_data;
  double *indloc2_data;
  double *indvec2_loop2_data;
  double *indvec_data;
  double *indvec_loop2_data;
  double *kweig_data;
  double *mat2_data;
  double *sPWM2_data;
  double *sPWM_data;
  int b_i;
  int b_k;
  int b_m;
  int f;
  int i;
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
  /* 'gkmPWM:750' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_13(p, i);
  p_data = p->data;

  /* 'gkmPWM:751' p = coder.nullcopy(p); */
  /* 'gkmPWM:752' p{1} = eye(4); */
  for (i = 0; i < 16; i++) {
    p_data[0].f1[i] = 0.0;
  }

  p_data[0].f1[0] = 1.0;
  p_data[0].f1[5] = 1.0;
  p_data[0].f1[10] = 1.0;
  p_data[0].f1[15] = 1.0;

  /* 'gkmPWM:753' for i = 1:l-1 */
  i = (int)(l - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:754' p{i+1} = p{i}*negmat; */
    for (b_k = 0; b_k < 4; b_k++) {
      for (indloc2_tmp = 0; indloc2_tmp < 4; indloc2_tmp++) {
        m = indloc2_tmp << 2;
        b_p[b_k + m] = ((p_data[b_i].f1[b_k] * negmat[m] + p_data[b_i].f1[b_k +
                         4] * negmat[m + 1]) + p_data[b_i].f1[b_k + 8] *
                        negmat[m + 2]) + p_data[b_i].f1[b_k + 12] * negmat[m + 3];
      }
    }

    for (b_k = 0; b_k < 16; b_k++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[b_k] = b_p[b_k];
    }
  }

  emxInit_real_T(&mat2, 2);
  emxInit_real_T(&indvec, 2);

  /* 'gkmPWM:756' n = 4^k*max(max(x)); */
  g_maximum(x, indvec);
  n_tmp = h_maximum(indvec);
  b_n_tmp = pow(4.0, k);
  n = b_n_tmp * n_tmp;

  /* number of possible k-mers */
  /* 'gkmPWM:757' mat2 = rot90(mat,2); */
  c_rot90(mat, mat2);
  mat2_data = mat2->data;

  /* 'gkmPWM:758' kweig = zeros(n, 4); */
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

  /* 'gkmPWM:759' ktree = cell(k-1,1); */
  m = (int)(k - 1.0);
  i = ktree->size[0];
  ktree->size[0] = (int)(k - 1.0);
  emxEnsureCapacity_cell_wrap_14(ktree, i);
  ktree_data = ktree->data;

  /* 'gkmPWM:760' ktree = coder.nullcopy(ktree); */
  /* 'gkmPWM:761' ktree2 = cell(k-1,1); */
  i = ktree2->size[0];
  ktree2->size[0] = (int)(k - 1.0);
  emxEnsureCapacity_cell_wrap_14(ktree2, i);
  ktree2_data = ktree2->data;
  for (i = 0; i < m; i++) {
    ktree_data[i].f1->size[0] = 0;
    ktree2_data[i].f1->size[0] = 0;
  }

  emxInit_uint32_T(&X);

  /* 'gkmPWM:762' ktree2 = coder.nullcopy(ktree2); */
  /* 'gkmPWM:763' [rx,cx] = size(x); */
  rx = x->size[0];

  /* 'gkmPWM:764' m=rx; */
  b_m = x->size[0];

  /* 'gkmPWM:765' M = l-1; */
  /* 'gkmPWM:766' X = cx*ones(M+1,1); */
  loop_ub_tmp = (int)((l - 1.0) + 1.0);
  i = X->size[0];
  X->size[0] = (int)((l - 1.0) + 1.0);
  emxEnsureCapacity_uint32_T(X, i);
  X_data = X->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    X_data[i] = (unsigned int)x->size[1];
  }

  /* 'gkmPWM:767' for i = 1:cx */
  i = x->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:768' X(i) = i; */
    X_data[b_i] = (unsigned int)(b_i + 1);
  }

  emxInit_real_T(&indloc2, 1);

  /* 'gkmPWM:770' indloc2 = flipud(indloc); */
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

  /* 'gkmPWM:771' for i = 2:5 */
  for (b_i = 0; b_i < 4; b_i++) {
    /* 'gkmPWM:772' ktree{i} = zeros(4^i,1); */
    c_tmp = pow(4.0, (double)b_i + 2.0);
    m = (int)pow(4.0, (double)b_i + 2.0);
    i = ktree_data[b_i + 1].f1->size[0];
    ktree_data[b_i + 1].f1->size[0] = (int)c_tmp;
    emxEnsureCapacity_real_T(ktree_data[b_i + 1].f1, i);
    for (i = 0; i < m; i++) {
      ktree_data[b_i + 1].f1->data[i] = 0.0;
    }

    /* 'gkmPWM:773' ktree2{i} = zeros(4^i,1); */
    i = ktree2_data[b_i + 1].f1->size[0];
    ktree2_data[b_i + 1].f1->size[0] = (int)c_tmp;
    emxEnsureCapacity_real_T(ktree2_data[b_i + 1].f1, i);
    for (i = 0; i < m; i++) {
      ktree2_data[b_i + 1].f1->data[i] = 0.0;
    }
  }

  /* 'gkmPWM:775' for i = 0:M */
  emxInit_int32_T(&f1, 2);
  emxInit_real_T(&sPWM, 2);
  emxInit_real_T(&sPWM2, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&indvec_loop2, 1);
  emxInit_real_T(&indvec2_loop2, 1);
  emxInit_boolean_T(&b_x, 2);
  emxInit_real_T(&b_kweig, 1);
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    /* 'gkmPWM:776' if i > M-cx+1 */
    if (b_i > ((l - 1.0) - (double)x->size[1]) + 1.0) {
      /* 'gkmPWM:777' m = length(c); */
      if ((c->size[0] == 0) || (c->size[1] == 0)) {
        b_m = 0;
      } else {
        m = c->size[0];
        b_m = c->size[1];
        if (m >= b_m) {
          b_m = m;
        }
      }
    }

    /* the following loops is basically dynamic programming for tensor multiplication.  there are multiple cases to consider, hence the if statements. */
    /* 'gkmPWM:780' for ii = 1:m */
    for (ii = 0; ii < b_m; ii++) {
      /* 'gkmPWM:781' if sum((c(ii,:)+i)==l) > 0 && ~(i == M-1 && ii > rx && ii ~= m) */
      md2 = c->size[1];
      i = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 1;
      b_x->size[1] = c->size[1];
      emxEnsureCapacity_boolean_T(b_x, i);
      b_x_data = b_x->data;
      for (i = 0; i < md2; i++) {
        b_x_data[i] = (c_data[ii + c->size[0] * i] + (double)b_i == l);
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

      if ((md2 > 0) && ((b_i != (l - 1.0) - 1.0) || (ii + 1 <= rx) || (ii + 1 ==
            b_m))) {
        /* 'gkmPWM:782' indvec = c(ii,:)+i; */
        md2 = c->size[1];
        i = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i);
        indvec_data = indvec->data;
        for (i = 0; i < md2; i++) {
          indvec_data[i] = c_data[ii + c->size[0] * i] + (double)b_i;
        }

        /* 'gkmPWM:783' coder.varsize('f', [1 1]); */
        /* 'gkmPWM:784' f1 = find(indvec == l); */
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = indvec->size[1];
        emxEnsureCapacity_boolean_T(b_x, i);
        b_x_data = b_x->data;
        md2 = indvec->size[1];
        for (i = 0; i < md2; i++) {
          b_x_data[i] = (indvec_data[i] == l);
        }

        eml_find(b_x, f1);
        f1_data = f1->data;

        /* 'gkmPWM:785' if length(f1) == 0 */
        if (f1->size[1] == 0) {
          /* 'gkmPWM:786' f = zeros(1,1); */
          f = -1;

          /* 'gkmPWM:787' fprintf("Resizing f\n"); */
          printf("Resizing f\n");
          fflush(stdout);
        } else {
          /* 'gkmPWM:788' else */
          /* 'gkmPWM:789' f = f1(1); */
          f = f1_data[0] - 1;
        }

        /* 'gkmPWM:791' indvec(f) = []; */
        indloc2_tmp = f + 1;
        m = indvec->size[1];
        md2 = indvec->size[1] - 1;
        for (b_k = indloc2_tmp; b_k <= md2; b_k++) {
          indvec_data[b_k - 1] = indvec_data[b_k];
        }

        i = indvec->size[0] * indvec->size[1];
        if (1 > md2) {
          indvec->size[1] = 0;
        } else {
          indvec->size[1] = m - 1;
        }

        emxEnsureCapacity_real_T(indvec, i);
        indvec_data = indvec->data;

        /* 'gkmPWM:792' loc = indloc(indvec); */
        /* 'gkmPWM:793' loc2 = indloc2(indvec); */
        /* 'gkmPWM:794' sPWM = mat(indvec,:).'; */
        i = indvec_loop2->size[0];
        indvec_loop2->size[0] = indvec->size[1];
        emxEnsureCapacity_real_T(indvec_loop2, i);
        indvec_loop2_data = indvec_loop2->data;
        md2 = indvec->size[1];
        for (i = 0; i < md2; i++) {
          indvec_loop2_data[i] = indvec_data[i];
        }

        i = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = 4;
        sPWM->size[1] = indvec_loop2->size[0];
        emxEnsureCapacity_real_T(sPWM, i);
        sPWM_data = sPWM->data;
        md2 = indvec_loop2->size[0];
        for (i = 0; i < md2; i++) {
          m = (int)indvec_loop2_data[i] - 1;
          sPWM_data[4 * i] = mat_data[m];
          sPWM_data[4 * i + 1] = mat_data[m + mat->size[0]];
          sPWM_data[4 * i + 2] = mat_data[m + mat->size[0] * 2];
          sPWM_data[4 * i + 3] = mat_data[m + mat->size[0] * 3];
        }

        /* 'gkmPWM:795' sPWM2 = mat2(indvec,:).'; */
        i = sPWM2->size[0] * sPWM2->size[1];
        sPWM2->size[0] = 4;
        sPWM2->size[1] = indvec_loop2->size[0];
        emxEnsureCapacity_real_T(sPWM2, i);
        sPWM2_data = sPWM2->data;
        md2 = indvec_loop2->size[0];
        for (i = 0; i < md2; i++) {
          m = (int)indvec_loop2_data[i] - 1;
          sPWM2_data[4 * i] = mat2_data[m];
          sPWM2_data[4 * i + 1] = mat2_data[m + mat2->size[0]];
          sPWM2_data[4 * i + 2] = mat2_data[m + mat2->size[0] * 2];
          sPWM2_data[4 * i + 3] = mat2_data[m + mat2->size[0] * 3];
        }

        /* 'gkmPWM:796' ktree{1} = sPWM(:,1); */
        i = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree_data[0].f1, i);

        /* 'gkmPWM:797' ktree2{1} = sPWM2(:,1); */
        i = ktree2_data[0].f1->size[0];
        ktree2_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree2_data[0].f1, i);
        ktree_data[0].f1->data[0] = sPWM_data[0];
        ktree2_data[0].f1->data[0] = sPWM2_data[0];
        ktree_data[0].f1->data[1] = sPWM_data[1];
        ktree2_data[0].f1->data[1] = sPWM2_data[1];
        ktree_data[0].f1->data[2] = sPWM_data[2];
        ktree2_data[0].f1->data[2] = sPWM2_data[2];
        ktree_data[0].f1->data[3] = sPWM_data[3];
        ktree2_data[0].f1->data[3] = sPWM2_data[3];

        /* 'gkmPWM:798' for iii = s(ii,i+1):k-1 */
        d = s_data[ii + s->size[0] * b_i];
        i = (int)((k - 1.0) + (1.0 - d));
        for (iii = 0; iii < i; iii++) {
          n = d + (double)iii;

          /* 'gkmPWM:799' if loc(iii)==0 */
          b_k = (int)indvec_data[(int)n - 1] - 1;
          if (indloc_data[b_k] == 0.0) {
            /* 'gkmPWM:800' if loc(iii-1)==1 && indvec(iii-1) < l */
            if ((indloc_data[(int)indvec_data[(int)(n - 1.0) - 1] - 1] == 1.0) &&
                ((unsigned int)indvec_data[(int)(n - 1.0) - 1] < l)) {
              /* 'gkmPWM:801' matt = sPWM(:,iii).'*p{indvec(iii)-l}; */
              for (b_k = 0; b_k < 16; b_k++) {
                b_p[b_k] = p_data[(int)((double)(unsigned int)indvec_data[(int)n
                  - 1] - l) - 1].f1[b_k];
              }

              /* 'gkmPWM:802' a = ktree2{iii-1}.*sPWM2(:,iii).'; */
              b_k = a->size[0] * a->size[1];
              a->size[0] = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, b_k);
              a_data = a->data;
              md2 = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
              for (b_k = 0; b_k < 4; b_k++) {
                for (indloc2_tmp = 0; indloc2_tmp < md2; indloc2_tmp++) {
                  a_data[indloc2_tmp + a->size[0] * b_k] = ktree2_data[(int)(n -
                    1.0) - 1].f1->data[indloc2_tmp] * sPWM2_data[b_k + 4 * ((int)
                    n - 1)];
                }
              }

              /* 'gkmPWM:803' ktree2{iii} = a(:); */
              b_k = ktree2_data[(int)n - 1].f1->size[0];
              ktree2_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, b_k);
              md2 = a->size[0] << 2;
              for (b_k = 0; b_k < md2; b_k++) {
                ktree2_data[(int)n - 1].f1->data[b_k] = a_data[b_k];
              }

              /* 'gkmPWM:804' ktree{iii} = repmat(ktree{iii-1},4,1).*repelem(matt', 4^(iii-1)); */
              c_repmat(ktree_data[(int)(n - 1.0) - 1].f1, indvec2_loop2);
              indvec2_loop2_data = indvec2_loop2->data;
              m = 4 * ((int)n - 1);
              for (b_k = 0; b_k < 4; b_k++) {
                indloc2_tmp = b_k << 2;
                b_sPWM2[b_k] = ((sPWM_data[m] * b_p[indloc2_tmp] + sPWM_data[m +
                                 1] * b_p[indloc2_tmp + 1]) + sPWM_data[m + 2] *
                                b_p[indloc2_tmp + 2]) + sPWM_data[m + 3] *
                  b_p[indloc2_tmp + 3];
              }

              repelem(b_sPWM2, pow(4.0, n - 1.0), indvec_loop2);
              indvec_loop2_data = indvec_loop2->data;
              if (indvec2_loop2->size[0] == indvec_loop2->size[0]) {
                b_k = ktree_data[(int)n - 1].f1->size[0];
                ktree_data[(int)n - 1].f1->size[0] = indvec2_loop2->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, b_k);
                md2 = indvec2_loop2->size[0];
                for (b_k = 0; b_k < md2; b_k++) {
                  ktree_data[(int)n - 1].f1->data[b_k] = indvec2_loop2_data[b_k]
                    * indvec_loop2_data[b_k];
                }
              } else {
                p_binary_expand_op(ktree, n, indvec2_loop2, indvec_loop2);
                ktree_data = ktree->data;
              }
            } else {
              /* 'gkmPWM:805' else */
              /* 'gkmPWM:806' matt = p{indvec(iii)-indvec(iii-1)+1}; */
              /* 'gkmPWM:807' ktree{iii} = repmat(ktree{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2)); */
              c_tmp = pow(4.0, n - 2.0);
              b_k = indvec_loop2->size[0];
              indvec_loop2->size[0] = (int)c_tmp << 4;
              emxEnsureCapacity_real_T(indvec_loop2, b_k);
              indvec_loop2_data = indvec_loop2->data;
              indloc2_tmp = -1;
              m = (int)c_tmp;
              for (b_k = 0; b_k < 16; b_k++) {
                for (j = 0; j < m; j++) {
                  indvec_loop2_data[(indloc2_tmp + j) + 1] = p_data[(int)
                    (unsigned int)indvec_data[(int)n - 1] - (int)(unsigned int)
                    indvec_data[(int)(n - 1.0) - 1]].f1[b_k];
                }

                indloc2_tmp += (int)c_tmp;
              }

              c_repmat(ktree_data[(int)(n - 1.0) - 1].f1, indvec2_loop2);
              indvec2_loop2_data = indvec2_loop2->data;
              if (indvec2_loop2->size[0] == indvec_loop2->size[0]) {
                b_k = ktree_data[(int)n - 1].f1->size[0];
                ktree_data[(int)n - 1].f1->size[0] = indvec2_loop2->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, b_k);
                md2 = indvec2_loop2->size[0];
                for (b_k = 0; b_k < md2; b_k++) {
                  ktree_data[(int)n - 1].f1->data[b_k] = indvec2_loop2_data[b_k]
                    * indvec_loop2_data[b_k];
                }
              } else {
                p_binary_expand_op(ktree, n, indvec2_loop2, indvec_loop2);
                ktree_data = ktree->data;
              }

              /* 'gkmPWM:808' a = ktree2{iii-1}.*sPWM2(:,iii).'; */
              b_k = a->size[0] * a->size[1];
              a->size[0] = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, b_k);
              a_data = a->data;
              md2 = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
              for (b_k = 0; b_k < 4; b_k++) {
                for (indloc2_tmp = 0; indloc2_tmp < md2; indloc2_tmp++) {
                  a_data[indloc2_tmp + a->size[0] * b_k] = ktree2_data[(int)(n -
                    1.0) - 1].f1->data[indloc2_tmp] * sPWM2_data[b_k + 4 * ((int)
                    n - 1)];
                }
              }

              /* 'gkmPWM:809' ktree2{iii} = a(:); */
              b_k = ktree2_data[(int)n - 1].f1->size[0];
              ktree2_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, b_k);
              md2 = a->size[0] << 2;
              for (b_k = 0; b_k < md2; b_k++) {
                ktree2_data[(int)n - 1].f1->data[b_k] = a_data[b_k];
              }
            }
          } else if (indloc2_data[b_k] == 0.0) {
            /* 'gkmPWM:811' elseif loc2(iii)==0 */
            /* 'gkmPWM:812' if loc2(iii-1)==1 && indvec(iii-1) < l */
            if ((indloc2_data[(int)indvec_data[(int)(n - 1.0) - 1] - 1] == 1.0) &&
                ((unsigned int)indvec_data[(int)(n - 1.0) - 1] < l)) {
              /* 'gkmPWM:813' matt = sPWM2(:,iii).'*p{indvec(iii)-l}; */
              for (b_k = 0; b_k < 16; b_k++) {
                b_p[b_k] = p_data[(int)((double)(unsigned int)indvec_data[(int)n
                  - 1] - l) - 1].f1[b_k];
              }

              /* 'gkmPWM:814' a = ktree{iii-1}.*sPWM(:,iii).'; */
              b_k = a->size[0] * a->size[1];
              a->size[0] = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, b_k);
              a_data = a->data;
              md2 = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
              for (b_k = 0; b_k < 4; b_k++) {
                for (indloc2_tmp = 0; indloc2_tmp < md2; indloc2_tmp++) {
                  a_data[indloc2_tmp + a->size[0] * b_k] = ktree_data[(int)(n -
                    1.0) - 1].f1->data[indloc2_tmp] * sPWM_data[b_k + 4 * ((int)
                    n - 1)];
                }
              }

              /* 'gkmPWM:815' ktree{iii} = a(:); */
              b_k = ktree_data[(int)n - 1].f1->size[0];
              ktree_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, b_k);
              md2 = a->size[0] << 2;
              for (b_k = 0; b_k < md2; b_k++) {
                ktree_data[(int)n - 1].f1->data[b_k] = a_data[b_k];
              }

              /* 'gkmPWM:816' ktree2{iii} = repmat(ktree2{iii-1},4,1).*repelem(matt', 4^(iii-1)); */
              c_repmat(ktree2_data[(int)(n - 1.0) - 1].f1, indvec2_loop2);
              indvec2_loop2_data = indvec2_loop2->data;
              for (b_k = 0; b_k < 4; b_k++) {
                indloc2_tmp = b_k << 2;
                b_sPWM2[b_k] = ((sPWM2_data[4 * ((int)n - 1)] * b_p[indloc2_tmp]
                                 + sPWM2_data[4 * ((int)n - 1) + 1] *
                                 b_p[indloc2_tmp + 1]) + sPWM2_data[4 * ((int)n
                  - 1) + 2] * b_p[indloc2_tmp + 2]) + sPWM2_data[4 * ((int)n - 1)
                  + 3] * b_p[indloc2_tmp + 3];
              }

              repelem(b_sPWM2, pow(4.0, n - 1.0), indvec_loop2);
              indvec_loop2_data = indvec_loop2->data;
              if (indvec2_loop2->size[0] == indvec_loop2->size[0]) {
                b_k = ktree2_data[(int)n - 1].f1->size[0];
                ktree2_data[(int)n - 1].f1->size[0] = indvec2_loop2->size[0];
                emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, b_k);
                md2 = indvec2_loop2->size[0];
                for (b_k = 0; b_k < md2; b_k++) {
                  ktree2_data[(int)n - 1].f1->data[b_k] = indvec2_loop2_data[b_k]
                    * indvec_loop2_data[b_k];
                }
              } else {
                p_binary_expand_op(ktree2, n, indvec2_loop2, indvec_loop2);
                ktree2_data = ktree2->data;
              }
            } else {
              /* 'gkmPWM:817' else */
              /* 'gkmPWM:818' matt = p{indvec(iii)-indvec(iii-1)+1}; */
              /* 'gkmPWM:819' ktree2{iii} = repmat(ktree2{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2)); */
              c_tmp = pow(4.0, n - 2.0);
              b_k = indvec_loop2->size[0];
              indvec_loop2->size[0] = (int)c_tmp << 4;
              emxEnsureCapacity_real_T(indvec_loop2, b_k);
              indvec_loop2_data = indvec_loop2->data;
              indloc2_tmp = -1;
              m = (int)c_tmp;
              for (b_k = 0; b_k < 16; b_k++) {
                for (j = 0; j < m; j++) {
                  indvec_loop2_data[(indloc2_tmp + j) + 1] = p_data[(int)
                    (unsigned int)indvec_data[(int)n - 1] - (int)(unsigned int)
                    indvec_data[(int)(n - 1.0) - 1]].f1[b_k];
                }

                indloc2_tmp += (int)c_tmp;
              }

              c_repmat(ktree2_data[(int)(n - 1.0) - 1].f1, indvec2_loop2);
              indvec2_loop2_data = indvec2_loop2->data;
              if (indvec2_loop2->size[0] == indvec_loop2->size[0]) {
                b_k = ktree2_data[(int)n - 1].f1->size[0];
                ktree2_data[(int)n - 1].f1->size[0] = indvec2_loop2->size[0];
                emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, b_k);
                md2 = indvec2_loop2->size[0];
                for (b_k = 0; b_k < md2; b_k++) {
                  ktree2_data[(int)n - 1].f1->data[b_k] = indvec2_loop2_data[b_k]
                    * indvec_loop2_data[b_k];
                }
              } else {
                p_binary_expand_op(ktree2, n, indvec2_loop2, indvec_loop2);
                ktree2_data = ktree2->data;
              }

              /* 'gkmPWM:820' a = ktree{iii-1}.*sPWM(:,iii).'; */
              b_k = a->size[0] * a->size[1];
              a->size[0] = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, b_k);
              a_data = a->data;
              md2 = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
              for (b_k = 0; b_k < 4; b_k++) {
                for (indloc2_tmp = 0; indloc2_tmp < md2; indloc2_tmp++) {
                  a_data[indloc2_tmp + a->size[0] * b_k] = ktree_data[(int)(n -
                    1.0) - 1].f1->data[indloc2_tmp] * sPWM_data[b_k + 4 * ((int)
                    n - 1)];
                }
              }

              /* 'gkmPWM:821' ktree{iii} = a(:); */
              b_k = ktree_data[(int)n - 1].f1->size[0];
              ktree_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, b_k);
              md2 = a->size[0] << 2;
              for (b_k = 0; b_k < md2; b_k++) {
                ktree_data[(int)n - 1].f1->data[b_k] = a_data[b_k];
              }
            }
          } else {
            /* 'gkmPWM:823' else */
            /* 'gkmPWM:824' a = ktree{iii-1}.*sPWM(:,iii).'; */
            b_k = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, b_k);
            a_data = a->data;
            md2 = ktree_data[(int)(n - 1.0) - 1].f1->size[0];
            for (b_k = 0; b_k < 4; b_k++) {
              for (indloc2_tmp = 0; indloc2_tmp < md2; indloc2_tmp++) {
                a_data[indloc2_tmp + a->size[0] * b_k] = ktree_data[(int)(n -
                  1.0) - 1].f1->data[indloc2_tmp] * sPWM_data[b_k + 4 * ((int)n
                  - 1)];
              }
            }

            /* 'gkmPWM:825' ktree{iii} = a(:); */
            b_k = ktree_data[(int)n - 1].f1->size[0];
            ktree_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree_data[(int)n - 1].f1, b_k);
            md2 = a->size[0] << 2;
            for (b_k = 0; b_k < md2; b_k++) {
              ktree_data[(int)n - 1].f1->data[b_k] = a_data[b_k];
            }

            /* 'gkmPWM:826' a = ktree2{iii-1}.*sPWM2(:,iii).'; */
            b_k = a->size[0] * a->size[1];
            a->size[0] = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, b_k);
            a_data = a->data;
            md2 = ktree2_data[(int)(n - 1.0) - 1].f1->size[0];
            for (b_k = 0; b_k < 4; b_k++) {
              for (indloc2_tmp = 0; indloc2_tmp < md2; indloc2_tmp++) {
                a_data[indloc2_tmp + a->size[0] * b_k] = ktree2_data[(int)(n -
                  1.0) - 1].f1->data[indloc2_tmp] * sPWM2_data[b_k + 4 * ((int)n
                  - 1)];
              }
            }

            /* 'gkmPWM:827' ktree2{iii} = a(:); */
            b_k = ktree2_data[(int)n - 1].f1->size[0];
            ktree2_data[(int)n - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree2_data[(int)n - 1].f1, b_k);
            md2 = a->size[0] << 2;
            for (b_k = 0; b_k < md2; b_k++) {
              ktree2_data[(int)n - 1].f1->data[b_k] = a_data[b_k];
            }
          }
        }

        /* the weird indexing that I did early in the code comes to fruition.  It is critical to do so to make this computation as fast as possible. */
        /* 'gkmPWM:831' if ii <= rx */
        if (ii + 1 <= rx) {
          /* 'gkmPWM:832' for j = 1:X(i+1) */
          i = (int)X_data[b_i];
          for (j = 0; j < i; j++) {
            /* 'gkmPWM:833' if x(ii,j) ~= 0 */
            d = x_data[ii + x->size[0] * j];
            if (d != 0.0) {
              /* 'gkmPWM:834' for iii = 1:2 */
              c_tmp = pow(4.0, (double)(f + 1) - 1.0);
              b_c = b_n_tmp * (ind_data[(int)d - 1] - 1.0);
              md2 = poscell_data[f].f1->size[0];
              for (iii = 0; iii < 2; iii++) {
                /* 'gkmPWM:835' indvec_loop1 = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(x(ii,j))-1); */
                n = c_tmp * (((double)iii + 1.0) - 1.0);
                b_k = indvec_loop2->size[0];
                indvec_loop2->size[0] = poscell_data[f].f1->size[0];
                emxEnsureCapacity_real_T(indvec_loop2, b_k);
                indvec_loop2_data = indvec_loop2->data;
                for (b_k = 0; b_k < md2; b_k++) {
                  indvec_loop2_data[b_k] = (poscell_data[f].f1->data[b_k] + n) +
                    b_c;
                }

                /* 'gkmPWM:836' indvec2_loop1 = indvec_loop1+(5-2*iii)*4^(f-1); */
                d = (5.0 - 2.0 * ((double)iii + 1.0)) * c_tmp;
                b_k = indvec2_loop2->size[0];
                indvec2_loop2->size[0] = indvec_loop2->size[0];
                emxEnsureCapacity_real_T(indvec2_loop2, b_k);
                indvec2_loop2_data = indvec2_loop2->data;
                m = indvec_loop2->size[0];
                for (b_k = 0; b_k < m; b_k++) {
                  indvec2_loop2_data[b_k] = indvec_loop2_data[b_k] + d;
                }

                /* 'gkmPWM:837' kweig(indvec_loop1,iii) = kweig(indvec_loop1,iii) + ktree{k-1}; */
                if (indvec_loop2->size[0] == ktree_data[(int)(k - 1.0) - 1]
                    .f1->size[0]) {
                  b_k = b_kweig->size[0];
                  b_kweig->size[0] = indvec_loop2->size[0];
                  emxEnsureCapacity_real_T(b_kweig, b_k);
                  a_data = b_kweig->data;
                  m = indvec_loop2->size[0];
                  for (b_k = 0; b_k < m; b_k++) {
                    a_data[b_k] = kweig_data[((int)indvec_loop2_data[b_k] +
                      kweig->size[0] * iii) - 1] + ktree_data[(int)(k - 1.0) - 1]
                      .f1->data[b_k];
                  }

                  m = b_kweig->size[0];
                  for (b_k = 0; b_k < m; b_k++) {
                    kweig_data[((int)indvec_loop2_data[b_k] + kweig->size[0] *
                                iii) - 1] = a_data[b_k];
                  }
                } else {
                  bb_binary_expand_op(kweig, indvec_loop2, iii, ktree, k);
                  kweig_data = kweig->data;
                }

                /* 'gkmPWM:838' kweig(indvec2_loop1,5-iii) = kweig(indvec2_loop1,5-iii) + ktree{k-1}; */
                if (indvec2_loop2->size[0] == ktree_data[(int)(k - 1.0) - 1].
                    f1->size[0]) {
                  b_k = b_kweig->size[0];
                  b_kweig->size[0] = indvec2_loop2->size[0];
                  emxEnsureCapacity_real_T(b_kweig, b_k);
                  a_data = b_kweig->data;
                  m = indvec2_loop2->size[0];
                  for (b_k = 0; b_k < m; b_k++) {
                    a_data[b_k] = kweig_data[((int)indvec2_loop2_data[b_k] +
                      kweig->size[0] * (3 - iii)) - 1] + ktree_data[(int)(k -
                      1.0) - 1].f1->data[b_k];
                  }

                  m = b_kweig->size[0];
                  for (b_k = 0; b_k < m; b_k++) {
                    kweig_data[((int)indvec2_loop2_data[b_k] + kweig->size[0] *
                                (3 - iii)) - 1] = a_data[b_k];
                  }
                } else {
                  ab_binary_expand_op(kweig, indvec2_loop2, iii, ktree, k);
                  kweig_data = kweig->data;
                }

                /* 'gkmPWM:839' kweig(indvec2_loop1,iii) = kweig(indvec2_loop1,iii) + ktree2{k-1}; */
                if (indvec2_loop2->size[0] == ktree2_data[(int)(k - 1.0) - 1].
                    f1->size[0]) {
                  b_k = b_kweig->size[0];
                  b_kweig->size[0] = indvec2_loop2->size[0];
                  emxEnsureCapacity_real_T(b_kweig, b_k);
                  a_data = b_kweig->data;
                  m = indvec2_loop2->size[0];
                  for (b_k = 0; b_k < m; b_k++) {
                    a_data[b_k] = kweig_data[((int)indvec2_loop2_data[b_k] +
                      kweig->size[0] * iii) - 1] + ktree2_data[(int)(k - 1.0) -
                      1].f1->data[b_k];
                  }

                  m = b_kweig->size[0];
                  for (b_k = 0; b_k < m; b_k++) {
                    kweig_data[((int)indvec2_loop2_data[b_k] + kweig->size[0] *
                                iii) - 1] = a_data[b_k];
                  }
                } else {
                  bb_binary_expand_op(kweig, indvec2_loop2, iii, ktree2, k);
                  kweig_data = kweig->data;
                }

                /* 'gkmPWM:840' kweig(indvec_loop1,5-iii) = kweig(indvec_loop1,5-iii) + ktree2{k-1}; */
                if (indvec_loop2->size[0] == ktree2_data[(int)(k - 1.0) - 1].
                    f1->size[0]) {
                  b_k = b_kweig->size[0];
                  b_kweig->size[0] = indvec_loop2->size[0];
                  emxEnsureCapacity_real_T(b_kweig, b_k);
                  a_data = b_kweig->data;
                  m = indvec_loop2->size[0];
                  for (b_k = 0; b_k < m; b_k++) {
                    a_data[b_k] = kweig_data[((int)indvec_loop2_data[b_k] +
                      kweig->size[0] * (3 - iii)) - 1] + ktree2_data[(int)(k -
                      1.0) - 1].f1->data[b_k];
                  }

                  m = b_kweig->size[0];
                  for (b_k = 0; b_k < m; b_k++) {
                    kweig_data[((int)indvec_loop2_data[b_k] + kweig->size[0] *
                                (3 - iii)) - 1] = a_data[b_k];
                  }
                } else {
                  ab_binary_expand_op(kweig, indvec_loop2, iii, ktree2, k);
                  kweig_data = kweig->data;
                }
              }
            }
          }
        } else {
          /* 'gkmPWM:844' else */
          /* 'gkmPWM:845' for iii = 1:2 */
          c_tmp = pow(4.0, (double)(f + 1) - 1.0);
          b_c = b_n_tmp * (ind_data[ii] - 1.0);
          md2 = poscell_data[f].f1->size[0];
          for (iii = 0; iii < 2; iii++) {
            /* 'gkmPWM:846' indvec_loop2 = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(ii)-1); */
            n = c_tmp * (((double)iii + 1.0) - 1.0);
            i = indvec_loop2->size[0];
            indvec_loop2->size[0] = poscell_data[f].f1->size[0];
            emxEnsureCapacity_real_T(indvec_loop2, i);
            indvec_loop2_data = indvec_loop2->data;
            for (i = 0; i < md2; i++) {
              indvec_loop2_data[i] = (poscell_data[f].f1->data[i] + n) + b_c;
            }

            /* 'gkmPWM:847' indvec2_loop2 = indvec_loop2+(5-2*iii)*4^(f-1); */
            d = (5.0 - 2.0 * ((double)iii + 1.0)) * c_tmp;
            i = indvec2_loop2->size[0];
            indvec2_loop2->size[0] = indvec_loop2->size[0];
            emxEnsureCapacity_real_T(indvec2_loop2, i);
            indvec2_loop2_data = indvec2_loop2->data;
            m = indvec_loop2->size[0];
            for (i = 0; i < m; i++) {
              indvec2_loop2_data[i] = indvec_loop2_data[i] + d;
            }

            /* 'gkmPWM:848' kweig(indvec_loop2,iii) = kweig(indvec_loop2,iii) + ktree{k-1}; */
            if (indvec_loop2->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->
                size[0]) {
              i = b_kweig->size[0];
              b_kweig->size[0] = indvec_loop2->size[0];
              emxEnsureCapacity_real_T(b_kweig, i);
              a_data = b_kweig->data;
              m = indvec_loop2->size[0];
              for (i = 0; i < m; i++) {
                a_data[i] = kweig_data[((int)indvec_loop2_data[i] + kweig->size
                  [0] * iii) - 1] + ktree_data[(int)(k - 1.0) - 1].f1->data[i];
              }

              m = b_kweig->size[0];
              for (i = 0; i < m; i++) {
                kweig_data[((int)indvec_loop2_data[i] + kweig->size[0] * iii) -
                  1] = a_data[i];
              }
            } else {
              bb_binary_expand_op(kweig, indvec_loop2, iii, ktree, k);
              kweig_data = kweig->data;
            }

            /* 'gkmPWM:849' kweig(indvec2_loop2,5-iii) = kweig(indvec2_loop2,5-iii) + ktree{k-1}; */
            if (indvec2_loop2->size[0] == ktree_data[(int)(k - 1.0) - 1]
                .f1->size[0]) {
              i = b_kweig->size[0];
              b_kweig->size[0] = indvec2_loop2->size[0];
              emxEnsureCapacity_real_T(b_kweig, i);
              a_data = b_kweig->data;
              m = indvec2_loop2->size[0];
              for (i = 0; i < m; i++) {
                a_data[i] = kweig_data[((int)indvec2_loop2_data[i] + kweig->
                  size[0] * (3 - iii)) - 1] + ktree_data[(int)(k - 1.0) - 1].
                  f1->data[i];
              }

              m = b_kweig->size[0];
              for (i = 0; i < m; i++) {
                kweig_data[((int)indvec2_loop2_data[i] + kweig->size[0] * (3 -
                  iii)) - 1] = a_data[i];
              }
            } else {
              ab_binary_expand_op(kweig, indvec2_loop2, iii, ktree, k);
              kweig_data = kweig->data;
            }

            /* 'gkmPWM:850' kweig(indvec2_loop2,iii) = kweig(indvec2_loop2,iii) + ktree2{k-1}; */
            if (indvec2_loop2->size[0] == ktree2_data[(int)(k - 1.0) - 1]
                .f1->size[0]) {
              i = b_kweig->size[0];
              b_kweig->size[0] = indvec2_loop2->size[0];
              emxEnsureCapacity_real_T(b_kweig, i);
              a_data = b_kweig->data;
              m = indvec2_loop2->size[0];
              for (i = 0; i < m; i++) {
                a_data[i] = kweig_data[((int)indvec2_loop2_data[i] + kweig->
                  size[0] * iii) - 1] + ktree2_data[(int)(k - 1.0) - 1].f1->
                  data[i];
              }

              m = b_kweig->size[0];
              for (i = 0; i < m; i++) {
                kweig_data[((int)indvec2_loop2_data[i] + kweig->size[0] * iii) -
                  1] = a_data[i];
              }
            } else {
              bb_binary_expand_op(kweig, indvec2_loop2, iii, ktree2, k);
              kweig_data = kweig->data;
            }

            /* 'gkmPWM:851' kweig(indvec_loop2,5-iii) = kweig(indvec_loop2,5-iii) + ktree2{k-1}; */
            if (indvec_loop2->size[0] == ktree2_data[(int)(k - 1.0) - 1]
                .f1->size[0]) {
              i = b_kweig->size[0];
              b_kweig->size[0] = indvec_loop2->size[0];
              emxEnsureCapacity_real_T(b_kweig, i);
              a_data = b_kweig->data;
              m = indvec_loop2->size[0];
              for (i = 0; i < m; i++) {
                a_data[i] = kweig_data[((int)indvec_loop2_data[i] + kweig->size
                  [0] * (3 - iii)) - 1] + ktree2_data[(int)(k - 1.0) - 1]
                  .f1->data[i];
              }

              m = b_kweig->size[0];
              for (i = 0; i < m; i++) {
                kweig_data[((int)indvec_loop2_data[i] + kweig->size[0] * (3 -
                  iii)) - 1] = a_data[i];
              }
            } else {
              ab_binary_expand_op(kweig, indvec_loop2, iii, ktree2, k);
              kweig_data = kweig->data;
            }
          }
        }
      }
    }
  }

  emxFree_real_T(&b_kweig);
  emxFree_boolean_T(&b_x);
  emxFree_cell_wrap_13(&p);
  emxFree_real_T(&indvec2_loop2);
  emxFree_real_T(&indvec_loop2);
  emxFree_real_T(&a);
  emxFree_real_T(&sPWM2);
  emxFree_real_T(&sPWM);
  emxFree_int32_T(&f1);
  emxFree_real_T(&indvec);
  emxFree_real_T(&indloc2);
  emxFree_uint32_T(&X);
  emxFree_cell_wrap_14(&ktree2);
  emxFree_cell_wrap_14(&ktree);

  /* 'gkmPWM:857' kweig(4^k*(max(max(x))-rcnum)+1:end,:) = kweig(4^k*(max(max(x))-rcnum)+1:end,:)/sqrt(2); */
  d = b_n_tmp * (n_tmp - rcnum) + 1.0;
  if (d > kweig->size[0]) {
    i = 0;
    b_k = 0;
    indloc2_tmp = 1;
  } else {
    i = (int)d - 1;
    b_k = kweig->size[0];
    indloc2_tmp = (int)d;
  }

  md2 = b_k - i;
  b_k = mat2->size[0] * mat2->size[1];
  mat2->size[0] = md2;
  mat2->size[1] = 4;
  emxEnsureCapacity_real_T(mat2, b_k);
  mat2_data = mat2->data;
  for (b_k = 0; b_k < 4; b_k++) {
    for (m = 0; m < md2; m++) {
      mat2_data[m + mat2->size[0] * b_k] = kweig_data[(i + m) + kweig->size[0] *
        b_k] / 1.4142135623730951;
    }
  }

  md2 = mat2->size[0];
  for (i = 0; i < 4; i++) {
    for (b_k = 0; b_k < md2; b_k++) {
      kweig_data[((indloc2_tmp + b_k) + kweig->size[0] * i) - 1] = mat2_data[b_k
        + mat2->size[0] * i];
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
  int idx;
  int ii;
  int iii;
  int itilerow;
  int j;
  int loop_ub_tmp;
  int m;
  int md2;
  int n;
  int ni;
  int rx;
  int x_tmp;
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
  /* 'gkmPWM:864' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_13(p, i);
  p_data = p->data;

  /* 'gkmPWM:865' p = coder.nullcopy(p); */
  /* 'gkmPWM:866' p{1} = eye(4); */
  for (i = 0; i < 16; i++) {
    p_data[0].f1[i] = 0.0;
  }

  p_data[0].f1[0] = 1.0;
  p_data[0].f1[5] = 1.0;
  p_data[0].f1[10] = 1.0;
  p_data[0].f1[15] = 1.0;

  /* 'gkmPWM:867' for i = 1:l-1 */
  i = (int)(l - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:868' p{i+1} = p{i}*negmat; */
    for (b_k = 0; b_k < 4; b_k++) {
      for (ni = 0; ni < 4; ni++) {
        n = ni << 2;
        b_p[b_k + n] = ((p_data[b_i].f1[b_k] * negmat[n] + p_data[b_i].f1[b_k +
                         4] * negmat[n + 1]) + p_data[b_i].f1[b_k + 8] *
                        negmat[n + 2]) + p_data[b_i].f1[b_k + 12] * negmat[n + 3];
      }
    }

    for (b_k = 0; b_k < 16; b_k++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[b_k] = b_p[b_k];
    }
  }

  emxInit_real_T(&indvec, 2);

  /* 'gkmPWM:870' n = 4^k*max(max(x)); */
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
  /* 'gkmPWM:871' mat2 = rot90(mat,2); */
  /* 'gkmPWM:872' kweig = zeros(n, 4); */
  i = kweig->size[0] * kweig->size[1];
  kweig->size[0] = (int)b_n;
  kweig->size[1] = 4;
  emxEnsureCapacity_real_T(kweig, i);
  kweig_data = kweig->data;
  md2 = (int)b_n << 2;
  for (i = 0; i < md2; i++) {
    kweig_data[i] = 0.0;
  }

  emxInit_cell_wrap_14(&ktree);

  /* 'gkmPWM:873' ktree = cell(k-1,1); */
  n = (int)(k - 1.0);
  i = ktree->size[0];
  ktree->size[0] = (int)(k - 1.0);
  emxEnsureCapacity_cell_wrap_14(ktree, i);
  ktree_data = ktree->data;
  for (i = 0; i < n; i++) {
    ktree_data[i].f1->size[0] = 0;
  }

  emxInit_uint32_T(&X);

  /* 'gkmPWM:874' ktree = coder.nullcopy(ktree); */
  /* 'gkmPWM:875' [rx,cx] = size(x); */
  rx = x->size[0];

  /* 'gkmPWM:876' m=rx; */
  m = x->size[0];

  /* 'gkmPWM:877' M = l-1; */
  /* 'gkmPWM:878' X = cx*ones(M+1,1); */
  loop_ub_tmp = (int)((l - 1.0) + 1.0);
  i = X->size[0];
  X->size[0] = (int)((l - 1.0) + 1.0);
  emxEnsureCapacity_uint32_T(X, i);
  X_data = X->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    X_data[i] = (unsigned int)x->size[1];
  }

  /* 'gkmPWM:879' for i = 1:cx */
  i = x->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:880' X(i) = i; */
    X_data[b_i] = (unsigned int)(b_i + 1);
  }

  emxInit_real_T(&indloc2, 1);

  /* 'gkmPWM:882' indloc2 = flipud(indloc); */
  i = indloc2->size[0];
  indloc2->size[0] = indloc->size[0];
  emxEnsureCapacity_real_T(indloc2, i);
  indloc2_data = indloc2->data;
  md2 = indloc->size[0];
  for (i = 0; i < md2; i++) {
    indloc2_data[i] = indloc_data[i];
  }

  n = indloc->size[0] - 1;
  md2 = indloc->size[0] >> 1;
  for (b_i = 0; b_i < md2; b_i++) {
    b_n = indloc2_data[b_i];
    ni = n - b_i;
    indloc2_data[b_i] = indloc2_data[ni];
    indloc2_data[ni] = b_n;
  }

  /* 'gkmPWM:883' for i = 2:5 */
  for (b_i = 0; b_i < 4; b_i++) {
    /* 'gkmPWM:884' ktree{i} = zeros(4^i,1); */
    n = (int)pow(4.0, (double)b_i + 2.0);
    i = ktree_data[b_i + 1].f1->size[0];
    ktree_data[b_i + 1].f1->size[0] = n;
    emxEnsureCapacity_real_T(ktree_data[b_i + 1].f1, i);
    for (i = 0; i < n; i++) {
      ktree_data[b_i + 1].f1->data[i] = 0.0;
    }
  }

  /* 'gkmPWM:886' for i = 0:M */
  emxInit_int32_T(&f1, 2);
  emxInit_real_T(&sPWM, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&indvec2, 1);
  emxInit_real_T(&b_indvec, 1);
  emxInit_boolean_T(&b_x, 2);
  emxInit_real_T(&b_kweig, 1);
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    /* 'gkmPWM:887' if i > M-cx+1 */
    if (b_i > ((l - 1.0) - (double)x->size[1]) + 1.0) {
      /* 'gkmPWM:888' m = length(c); */
      if ((c->size[0] == 0) || (c->size[1] == 0)) {
        m = 0;
      } else {
        n = c->size[0];
        m = c->size[1];
        if (n >= m) {
          m = n;
        }
      }
    }

    /* the following loops is basically dynamic programming for tensor multiplication.  there are multiple cases to consider, hence the if statements. */
    /* 'gkmPWM:891' for ii = 1:m */
    for (ii = 0; ii < m; ii++) {
      /* 'gkmPWM:892' if sum((c(ii,:)+i)==l) > 0 && ~(i == M-1 && ii > rx && ii ~= m) */
      md2 = c->size[1];
      i = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 1;
      b_x->size[1] = c->size[1];
      emxEnsureCapacity_boolean_T(b_x, i);
      b_x_data = b_x->data;
      for (i = 0; i < md2; i++) {
        b_x_data[i] = (c_data[ii + c->size[0] * i] + (double)b_i == l);
      }

      n = b_x->size[1];
      if (b_x->size[1] == 0) {
        md2 = 0;
      } else {
        md2 = b_x_data[0];
        for (b_k = 2; b_k <= n; b_k++) {
          md2 += b_x_data[b_k - 1];
        }
      }

      if ((md2 > 0) && ((b_i != (l - 1.0) - 1.0) || (ii + 1 <= rx) || (ii + 1 ==
            m))) {
        /* 'gkmPWM:893' indvec = c(ii,:)+i; */
        md2 = c->size[1];
        i = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i);
        indvec_data = indvec->data;
        for (i = 0; i < md2; i++) {
          indvec_data[i] = c_data[ii + c->size[0] * i] + (double)b_i;
        }

        /* 'gkmPWM:894' coder.varsize('f', [1 1]); */
        /* 'gkmPWM:895' f1 = find(indvec == l); */
        i = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = indvec->size[1];
        emxEnsureCapacity_boolean_T(b_x, i);
        b_x_data = b_x->data;
        md2 = indvec->size[1];
        for (i = 0; i < md2; i++) {
          b_x_data[i] = (indvec_data[i] == l);
        }

        eml_find(b_x, f1);
        f1_data = f1->data;

        /* 'gkmPWM:896' if length(f1) == 0 */
        if (f1->size[1] == 0) {
          /* 'gkmPWM:897' f = zeros(1,1); */
          f = -1;

          /* 'gkmPWM:898' fprintf("Resizing f\n"); */
          printf("Resizing f\n");
          fflush(stdout);
        } else {
          /* 'gkmPWM:899' else */
          /* 'gkmPWM:900' f = f1(1); */
          f = f1_data[0] - 1;
        }

        /* 'gkmPWM:902' indvec(f) = []; */
        idx = f + 1;
        n = indvec->size[1];
        md2 = indvec->size[1] - 1;
        for (b_k = idx; b_k <= md2; b_k++) {
          indvec_data[b_k - 1] = indvec_data[b_k];
        }

        i = indvec->size[0] * indvec->size[1];
        if (1 > md2) {
          indvec->size[1] = 0;
        } else {
          indvec->size[1] = n - 1;
        }

        emxEnsureCapacity_real_T(indvec, i);
        indvec_data = indvec->data;

        /* 'gkmPWM:903' loc = indloc(indvec); */
        /* 'gkmPWM:904' loc2 = indloc2(indvec); */
        /* 'gkmPWM:905' sPWM = mat(indvec,:).'; */
        i = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = 4;
        sPWM->size[1] = indvec->size[1];
        emxEnsureCapacity_real_T(sPWM, i);
        sPWM_data = sPWM->data;
        md2 = indvec->size[1];
        for (i = 0; i < md2; i++) {
          n = (int)indvec_data[i] - 1;
          sPWM_data[4 * i] = mat_data[n];
          sPWM_data[4 * i + 1] = mat_data[n + mat->size[0]];
          sPWM_data[4 * i + 2] = mat_data[n + mat->size[0] * 2];
          sPWM_data[4 * i + 3] = mat_data[n + mat->size[0] * 3];
        }

        /* 'gkmPWM:906' ktree{1} = sPWM(:,1); */
        i = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree_data[0].f1, i);
        ktree_data[0].f1->data[0] = sPWM_data[0];
        ktree_data[0].f1->data[1] = sPWM_data[1];
        ktree_data[0].f1->data[2] = sPWM_data[2];
        ktree_data[0].f1->data[3] = sPWM_data[3];

        /* 'gkmPWM:907' for iii = s(ii,i+1):k-1 */
        d = s_data[ii + s->size[0] * b_i];
        i = (int)((k - 1.0) + (1.0 - d));
        for (iii = 0; iii < i; iii++) {
          b_n = d + (double)iii;

          /* 'gkmPWM:908' if loc(iii)==0 */
          b_k = (int)indvec_data[(int)b_n - 1] - 1;
          if (indloc_data[b_k] == 0.0) {
            /* 'gkmPWM:909' if loc(iii-1)==1 && indvec(iii-1) < l */
            if ((indloc_data[(int)indvec_data[(int)(b_n - 1.0) - 1] - 1] == 1.0)
                && ((unsigned int)indvec_data[(int)(b_n - 1.0) - 1] < l)) {
              /* 'gkmPWM:910' matt = sPWM(:,iii).'*p{indvec(iii)-l}; */
              for (b_k = 0; b_k < 16; b_k++) {
                b_p[b_k] = p_data[(int)((double)(unsigned int)indvec_data[(int)
                  b_n - 1] - l) - 1].f1[b_k];
              }

              /* 'gkmPWM:911' ktree{iii} = repmat(ktree{iii-1},4,1).*repelem(matt', 4^(iii-1)); */
              b_k = b_indvec->size[0];
              b_indvec->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0] <<
                2;
              emxEnsureCapacity_real_T(b_indvec, b_k);
              b_indvec_data = b_indvec->data;
              md2 = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              varargin_1 = pow(4.0, b_n - 1.0);
              b_k = indvec2->size[0];
              indvec2->size[0] = (int)varargin_1 << 2;
              emxEnsureCapacity_real_T(indvec2, b_k);
              indvec2_data = indvec2->data;
              idx = -1;
              ni = (int)varargin_1;
              x_tmp = 4 * ((int)b_n - 1);
              for (itilerow = 0; itilerow < 4; itilerow++) {
                n = itilerow * md2;
                for (b_k = 0; b_k < md2; b_k++) {
                  b_indvec_data[n + b_k] = ktree_data[(int)(b_n - 1.0) - 1]
                    .f1->data[b_k];
                }

                b_k = itilerow << 2;
                c_x[itilerow] = ((sPWM_data[x_tmp] * b_p[b_k] + sPWM_data[x_tmp
                                  + 1] * b_p[b_k + 1]) + sPWM_data[x_tmp + 2] *
                                 b_p[b_k + 2]) + sPWM_data[x_tmp + 3] * b_p[b_k
                  + 3];
                for (j = 0; j < ni; j++) {
                  indvec2_data[(idx + j) + 1] = c_x[itilerow];
                }

                idx += (int)varargin_1;
              }

              if (b_indvec->size[0] == indvec2->size[0]) {
                b_k = ktree_data[(int)b_n - 1].f1->size[0];
                ktree_data[(int)b_n - 1].f1->size[0] = b_indvec->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
                md2 = b_indvec->size[0];
                for (b_k = 0; b_k < md2; b_k++) {
                  ktree_data[(int)b_n - 1].f1->data[b_k] = b_indvec_data[b_k] *
                    indvec2_data[b_k];
                }
              } else {
                p_binary_expand_op(ktree, b_n, b_indvec, indvec2);
                ktree_data = ktree->data;
              }
            } else {
              /* 'gkmPWM:912' else */
              /* 'gkmPWM:913' matt = p{indvec(iii)-indvec(iii-1)+1}; */
              /* 'gkmPWM:914' ktree{iii} = repmat(ktree{iii-1}, 4, 1).*repelem(matt(:), 4^(iii-2)); */
              b_k = b_indvec->size[0];
              b_indvec->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0] <<
                2;
              emxEnsureCapacity_real_T(b_indvec, b_k);
              b_indvec_data = b_indvec->data;
              md2 = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              for (itilerow = 0; itilerow < 4; itilerow++) {
                n = itilerow * md2;
                for (b_k = 0; b_k < md2; b_k++) {
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
              ni = (int)varargin_1;
              for (b_k = 0; b_k < 16; b_k++) {
                for (j = 0; j < ni; j++) {
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
                md2 = b_indvec->size[0];
                for (b_k = 0; b_k < md2; b_k++) {
                  ktree_data[(int)b_n - 1].f1->data[b_k] = b_indvec_data[b_k] *
                    indvec2_data[b_k];
                }
              } else {
                p_binary_expand_op(ktree, b_n, b_indvec, indvec2);
                ktree_data = ktree->data;
              }
            }
          } else if (indloc2_data[b_k] == 0.0) {
            /* 'gkmPWM:916' elseif loc2(iii)==0 */
            /* 'gkmPWM:917' if loc2(iii-1)==1 && indvec(iii-1) < l */
            if ((indloc2_data[(int)indvec_data[(int)(b_n - 1.0) - 1] - 1] == 1.0)
                && ((unsigned int)indvec_data[(int)(b_n - 1.0) - 1] < l)) {
              /* 'gkmPWM:918' a = ktree{iii-1}.*sPWM(:,iii).'; */
              b_k = a->size[0] * a->size[1];
              a->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, b_k);
              a_data = a->data;
              md2 = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              for (b_k = 0; b_k < 4; b_k++) {
                for (ni = 0; ni < md2; ni++) {
                  a_data[ni + a->size[0] * b_k] = ktree_data[(int)(b_n - 1.0) -
                    1].f1->data[ni] * sPWM_data[b_k + 4 * ((int)b_n - 1)];
                }
              }

              /* 'gkmPWM:919' ktree{iii} = a(:); */
              b_k = ktree_data[(int)b_n - 1].f1->size[0];
              ktree_data[(int)b_n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
              md2 = a->size[0] << 2;
              for (b_k = 0; b_k < md2; b_k++) {
                ktree_data[(int)b_n - 1].f1->data[b_k] = a_data[b_k];
              }
            } else {
              /* 'gkmPWM:920' else */
              /* 'gkmPWM:921' a = ktree{iii-1}.*sPWM(:,iii).'; */
              b_k = a->size[0] * a->size[1];
              a->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              a->size[1] = 4;
              emxEnsureCapacity_real_T(a, b_k);
              a_data = a->data;
              md2 = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
              for (b_k = 0; b_k < 4; b_k++) {
                for (ni = 0; ni < md2; ni++) {
                  a_data[ni + a->size[0] * b_k] = ktree_data[(int)(b_n - 1.0) -
                    1].f1->data[ni] * sPWM_data[b_k + 4 * ((int)b_n - 1)];
                }
              }

              /* 'gkmPWM:922' ktree{iii} = a(:); */
              b_k = ktree_data[(int)b_n - 1].f1->size[0];
              ktree_data[(int)b_n - 1].f1->size[0] = a->size[0] << 2;
              emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
              md2 = a->size[0] << 2;
              for (b_k = 0; b_k < md2; b_k++) {
                ktree_data[(int)b_n - 1].f1->data[b_k] = a_data[b_k];
              }
            }
          } else {
            /* 'gkmPWM:924' else */
            /* 'gkmPWM:925' a = ktree{iii-1}.*sPWM(:,iii).'; */
            b_k = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, b_k);
            a_data = a->data;
            md2 = ktree_data[(int)(b_n - 1.0) - 1].f1->size[0];
            for (b_k = 0; b_k < 4; b_k++) {
              for (ni = 0; ni < md2; ni++) {
                a_data[ni + a->size[0] * b_k] = ktree_data[(int)(b_n - 1.0) - 1]
                  .f1->data[ni] * sPWM_data[b_k + 4 * ((int)b_n - 1)];
              }
            }

            /* 'gkmPWM:926' ktree{iii} = a(:); */
            b_k = ktree_data[(int)b_n - 1].f1->size[0];
            ktree_data[(int)b_n - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree_data[(int)b_n - 1].f1, b_k);
            md2 = a->size[0] << 2;
            for (b_k = 0; b_k < md2; b_k++) {
              ktree_data[(int)b_n - 1].f1->data[b_k] = a_data[b_k];
            }
          }
        }

        /* the weird indexing that I did early in the code comes to fruition.  It is critical to do so to make this computation as fast as possible. */
        /* 'gkmPWM:930' if ii <= rx */
        if (ii + 1 <= rx) {
          /* 'gkmPWM:931' for j = 1:X(i+1) */
          i = (int)X_data[b_i];
          for (j = 0; j < i; j++) {
            /* 'gkmPWM:932' if x(ii,j) ~= 0 */
            d = x_data[ii + x->size[0] * j];
            if (d != 0.0) {
              /* 'gkmPWM:933' for iii = 1:2 */
              varargin_1 = pow(4.0, (double)(f + 1) - 1.0);
              b_c = n_tmp * (ind_data[(int)d - 1] - 1.0);
              md2 = poscell_data[f].f1->size[0];
              for (iii = 0; iii < 2; iii++) {
                /* 'gkmPWM:934' indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(x(ii,j))-1); */
                b_n = varargin_1 * (((double)iii + 1.0) - 1.0);
                b_k = b_indvec->size[0];
                b_indvec->size[0] = poscell_data[f].f1->size[0];
                emxEnsureCapacity_real_T(b_indvec, b_k);
                b_indvec_data = b_indvec->data;
                for (b_k = 0; b_k < md2; b_k++) {
                  b_indvec_data[b_k] = (poscell_data[f].f1->data[b_k] + b_n) +
                    b_c;
                }

                /* 'gkmPWM:935' indvec2 = indvec+(5-2*iii)*4^(f-1); */
                d = (5.0 - 2.0 * ((double)iii + 1.0)) * varargin_1;
                b_k = indvec2->size[0];
                indvec2->size[0] = b_indvec->size[0];
                emxEnsureCapacity_real_T(indvec2, b_k);
                indvec2_data = indvec2->data;
                n = b_indvec->size[0];
                for (b_k = 0; b_k < n; b_k++) {
                  indvec2_data[b_k] = b_indvec_data[b_k] + d;
                }

                /* 'gkmPWM:936' kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1}; */
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

                /* 'gkmPWM:937' kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1}; */
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
          /* 'gkmPWM:941' else */
          /* 'gkmPWM:942' for iii = 1:2 */
          varargin_1 = pow(4.0, (double)(f + 1) - 1.0);
          b_c = n_tmp * (ind_data[ii] - 1.0);
          md2 = poscell_data[f].f1->size[0];
          for (iii = 0; iii < 2; iii++) {
            /* 'gkmPWM:943' indvec = poscell{f}+4^(f-1)*(iii-1)+4^k*(ind(ii)-1); */
            b_n = varargin_1 * (((double)iii + 1.0) - 1.0);
            i = b_indvec->size[0];
            b_indvec->size[0] = poscell_data[f].f1->size[0];
            emxEnsureCapacity_real_T(b_indvec, i);
            b_indvec_data = b_indvec->data;
            for (i = 0; i < md2; i++) {
              b_indvec_data[i] = (poscell_data[f].f1->data[i] + b_n) + b_c;
            }

            /* 'gkmPWM:944' indvec2 = indvec+(5-2*iii)*4^(f-1); */
            d = (5.0 - 2.0 * ((double)iii + 1.0)) * varargin_1;
            i = indvec2->size[0];
            indvec2->size[0] = b_indvec->size[0];
            emxEnsureCapacity_real_T(indvec2, i);
            indvec2_data = indvec2->data;
            n = b_indvec->size[0];
            for (i = 0; i < n; i++) {
              indvec2_data[i] = b_indvec_data[i] + d;
            }

            /* 'gkmPWM:945' kweig(indvec,iii) = kweig(indvec,iii) + ktree{k-1}; */
            if (b_indvec->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->size[0])
            {
              i = b_kweig->size[0];
              b_kweig->size[0] = b_indvec->size[0];
              emxEnsureCapacity_real_T(b_kweig, i);
              a_data = b_kweig->data;
              n = b_indvec->size[0];
              for (i = 0; i < n; i++) {
                a_data[i] = kweig_data[((int)b_indvec_data[i] + kweig->size[0] *
                  iii) - 1] + ktree_data[(int)(k - 1.0) - 1].f1->data[i];
              }

              n = b_kweig->size[0];
              for (i = 0; i < n; i++) {
                kweig_data[((int)b_indvec_data[i] + kweig->size[0] * iii) - 1] =
                  a_data[i];
              }
            } else {
              bb_binary_expand_op(kweig, b_indvec, iii, ktree, k);
              kweig_data = kweig->data;
            }

            /* 'gkmPWM:946' kweig(indvec2,5-iii) = kweig(indvec2,5-iii) + ktree{k-1}; */
            if (indvec2->size[0] == ktree_data[(int)(k - 1.0) - 1].f1->size[0])
            {
              i = b_kweig->size[0];
              b_kweig->size[0] = indvec2->size[0];
              emxEnsureCapacity_real_T(b_kweig, i);
              a_data = b_kweig->data;
              n = indvec2->size[0];
              for (i = 0; i < n; i++) {
                a_data[i] = kweig_data[((int)indvec2_data[i] + kweig->size[0] *
                  (3 - iii)) - 1] + ktree_data[(int)(k - 1.0) - 1].f1->data[i];
              }

              n = b_kweig->size[0];
              for (i = 0; i < n; i++) {
                kweig_data[((int)indvec2_data[i] + kweig->size[0] * (3 - iii)) -
                  1] = a_data[i];
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

static void m_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
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

static void n_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
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

static void o_binary_expand_op(emxArray_boolean_T *x, const emxArray_cell_wrap_1
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

static void plus(emxArray_real_T *b, const emxArray_real_T *res)
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

  /* 'gkmPWM:208' fid = fopen(fn, 'r'); */
  fileid = cfopen(fn, "rb");

  /*  a = textscan(fid, '%s\t%f\n'); */
  /* 'gkmPWM:211' curr_pos = ftell(fid); */
  curr_pos = b_ftell(fileid);

  /* 'gkmPWM:212' idx=0; */
  idx = 0.0;

  /* 'gkmPWM:213' while ~feof(fid) */
  emxInit_char_T(&b_fileid, 2);
  do {
    exitg1 = 0;
    d = b_feof(fileid);
    if (d == 0.0) {
      /* 'gkmPWM:214' idx=idx+1; */
      idx++;

      /* 'gkmPWM:215' fgetl(fid); */
      b_fgets(fileid, b_fileid);
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  emxFree_char_T(&b_fileid);
  emxInit_cell_wrap_0(&sequences, 1);

  /* 'gkmPWM:217' fseek(fid, curr_pos, 'bof'); */
  b_fseek(fileid, curr_pos);

  /* 'gkmPWM:218' sequences = cell(idx, 1); */
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

  /* 'gkmPWM:219' sequences = coder.nullcopy(sequences); */
  /* 'gkmPWM:220' alpha = zeros(idx, 1); */
  i = alpha->size[0];
  alpha->size[0] = (int)idx;
  emxEnsureCapacity_real_T(alpha, i);
  alpha_data = alpha->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    alpha_data[i] = 0.0;
  }

  /* 'gkmPWM:221' for cur_idx=1:idx */
  cur_idx = 0;
  emxInit_char_T(&cur_line, 2);
  emxInit_char_T(&cur_seq, 2);
  emxInit_char_T(&cur_alpha, 2);
  exitg2 = false;
  while ((!exitg2) && (cur_idx <= (int)idx - 1)) {
    /* 'gkmPWM:222' cur_line = fgetl(fid); */
    fgetl(fileid, cur_line);

    /* 'gkmPWM:223' if cur_line == -1 */
    y = (cur_line->size[1] != 0);
    if (y) {
      y = (0 > cur_line->size[1] - 1);
    }

    if (y) {
      exitg2 = true;
    } else {
      /* 'gkmPWM:226' [cur_seq, cur_alpha] = strtok(cur_line, char(9)); */
      b_strtok(cur_line, cur_seq, cur_alpha);

      /* 'gkmPWM:227' alpha(cur_idx,1) = real(str2double(cur_alpha)); */
      dc = str2double(cur_alpha);
      alpha_data[cur_idx] = dc.re;

      /* 'gkmPWM:228' sequences{cur_idx} = (strip(cur_seq)); */
      strip(cur_seq, sequences_data[cur_idx].f1);
      cur_idx++;
    }
  }

  emxFree_char_T(&cur_alpha);
  emxFree_char_T(&cur_seq);
  emxFree_char_T(&cur_line);
  emxInit_real_T(&M, 1);
  emxInit_int32_T(&D, 1);

  /* 'gkmPWM:230' fclose(fid); */
  cfclose(fileid);

  /*  [w, ind] = sort(a{2}, pn); */
  /*  s = a{1}(ind(1:min([100000 length(a{1})]))); */
  /* 'gkmPWM:235' [w, ind] = sort(alpha, pn); */
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

  /* 'gkmPWM:236' s_len = min([100000 length(sequences)]); */
  m[0] = 100000.0;
  m[1] = sequences->size[0];

  /* 'gkmPWM:237' s = cell(s_len, 1); */
  unnamed_idx_0_tmp_tmp = (int)b_minimum(m);
  i = s->size[0];
  s->size[0] = unnamed_idx_0_tmp_tmp;
  emxEnsureCapacity_cell_wrap_0(s, i);
  s_data = s->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    s_data[i].f1->size[0] = 1;
    s_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:238' s = coder.nullcopy(s); */
  /* 'gkmPWM:239' for cur_idx=1:s_len */
  for (cur_idx = 0; cur_idx < unnamed_idx_0_tmp_tmp; cur_idx++) {
    /* 'gkmPWM:240' s{cur_idx} = sequences{ind(cur_idx)}; */
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

  /* 'gkmPWM:244' l = length(s{1}); */
  varargin_2 = s_data[0].f1->size[1] - 1;
  l = s_data[0].f1->size[1] - 4;

  /* 'gkmPWM:245' k = round(l/2)+1; */
  x = (int)rt_roundd((double)s_data[0].f1->size[1] / 2.0);

  /* 'gkmPWM:246' ikl = length(ik); */
  /* 'gkmPWM:247' p = cell(num,1); */
  unnamed_idx_0_tmp_tmp = (int)num;
  i = b_p->size[0];
  b_p->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_0(b_p, i);
  p_data = b_p->data;
  for (i = 0; i < unnamed_idx_0_tmp_tmp; i++) {
    p_data[i].f1->size[0] = 1;
    p_data[i].f1->size[1] = 0;
  }

  /* 'gkmPWM:248' p = coder.nullcopy(p); */
  i = sequences->size[0];
  sequences->size[0] = b_p->size[0];
  emxEnsureCapacity_cell_wrap_0(sequences, i);
  sequences_data = sequences->data;

  /* 'gkmPWM:249' c = ikl+1; */
  *c = 1.0;

  /* 'gkmPWM:250' p{1} = s{1}; */
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

  /* 'gkmPWM:251' mat = cell(ikl+num,1); */
  i = b_mat->size[0];
  b_mat->size[0] = (int)num;
  emxEnsureCapacity_cell_wrap_1(b_mat, i);
  mat_data = b_mat->data;

  /* 'gkmPWM:252' mat = coder.nullcopy(mat); */
  /*  mat(1:ikl) = ik; */
  /* 'gkmPWM:254' for cur_idx=1:length(ik) */
  /* 'gkmPWM:257' mat{c} = letterconvert(s{1}); */
  /* 'gkmPWM:258' pwms = cell(num,1); */
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

  /* 'gkmPWM:259' pwms = coder.nullcopy(pwms); */
  /* 'gkmPWM:260' for i = 1:num */
  for (b_i = 0; b_i < unnamed_idx_0_tmp_tmp; b_i++) {
    /* 'gkmPWM:261' pwms{i} = zeros(l,4); */
    i = pwms_data[b_i].f1->size[0] * pwms_data[b_i].f1->size[1];
    pwms_data[b_i].f1->size[0] = varargin_2 + 1;
    pwms_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(pwms_data[b_i].f1, i);
    loop_ub = (varargin_2 + 1) << 2;
    for (i = 0; i < loop_ub; i++) {
      pwms_data[b_i].f1->data[i] = 0.0;
    }
  }

  /* 'gkmPWM:263' for i = 1:l */
  i = s_data[0].f1->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:264' pwms{1}(i,mat{c}(i)+1) = pwms{1}(i,mat{c}(i)+1)+w(i); */
    d = mat_data[0].f1->data[b_i];
    pwms_data[0].f1->data[b_i + pwms_data[0].f1->size[0] * ((int)(d + 1.0) - 1)]
      += alpha_data[b_i];
  }

  /* 'gkmPWM:266' B = zeros(9,1); */
  /* 'gkmPWM:267' BB = zeros(9,1); */
  for (b_i = 0; b_i < 9; b_i++) {
    B[b_i] = 0;
    BB[b_i] = 0;
  }

  /* 'gkmPWM:268' B(1:5) = (0:4)'; */
  /* 'gkmPWM:269' B(6:9) = 0; */
  B[5] = 0;
  B[6] = 0;
  B[7] = 0;
  B[8] = 0;

  /* 'gkmPWM:270' BB(1:5) = 0; */
  for (b_i = 0; b_i < 5; b_i++) {
    B[b_i] = (signed char)b_i;
    BB[b_i] = 0;
  }

  /* 'gkmPWM:271' BB(6:9) = (1:4)'; */
  BB[5] = 1;
  BB[6] = 2;
  BB[7] = 3;
  BB[8] = 4;

  /* 'gkmPWM:272' CC = [l l-1 l-2 l-3 l-4 l-1 l-2 l-3 l-4]; */
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
  /* 'gkmPWM:274' for i = 2:100000 */
  b_i = 1;
  emxInit_real_T(&ss, 2);
  emxInit_real_T(&rs, 2);
  emxInit_int8_T(&DD, 1);
  emxInit_boolean_T(&b_x, 2);
  exitg2 = false;
  while ((!exitg2) && (b_i - 1 < 99999)) {
    /* 'gkmPWM:275' ss = letterconvert(s{i}); */
    letterconvert(s_data[b_i].f1, ss);
    ss_data = ss->data;

    /* 'gkmPWM:276' rs = 3-fliplr(ss); */
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

    /* 'gkmPWM:277' M = zeros(c,1); */
    loop_ub_tmp = (int)*c;
    i = M->size[0];
    M->size[0] = (int)*c;
    emxEnsureCapacity_real_T(M, i);
    M_data = M->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      M_data[i] = 0.0;
    }

    /* 'gkmPWM:278' D = zeros(c,1); */
    i = D->size[0];
    D->size[0] = (int)*c;
    emxEnsureCapacity_int32_T(D, i);
    D_data = D->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      D_data[i] = 0;
    }

    /* 'gkmPWM:279' DD = zeros(c,1); */
    i = DD->size[0];
    DD->size[0] = (int)*c;
    emxEnsureCapacity_int8_T(DD, i);
    DD_data = DD->data;
    for (i = 0; i < loop_ub_tmp; i++) {
      DD_data[i] = 0;
    }

    /* 'gkmPWM:280' for j = 1:c */
    for (j = 0; j < loop_ub_tmp; j++) {
      /* 'gkmPWM:281' [m,d] = max([sum(mat{j}==ss) sum(mat{j}(2:end)==ss(1:l-1)) sum(mat{j}(3:end)==ss(1:l-2)) sum(mat{j}(4:end)==ss(1:l-3)) sum(mat{j}(5:end)==ss(1:l-4)) sum(mat{j}(1:l-1)==ss(2:end)) sum(mat{j}(1:l-2)==ss(3:end)) sum(mat{j}(1:l-3)==ss(4:end)) sum(mat{j}(1:l-4)==ss(5:end))]); */
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
        o_binary_expand_op(b_x, b_mat, j, ss);
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
        n_binary_expand_op(b_x, b_mat, j, i, b_y - 1, ss, unnamed_idx_0_tmp_tmp
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
        m_binary_expand_op(b_x, b_mat, j, cur_idx, c_y - 1, ss, d_y - 1);
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
        l_binary_expand_op(b_x, b_mat, j, e_y, i1 - 1, ss, i2);
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
        l_binary_expand_op(b_x, b_mat, j, i3, i4 - 1, ss, i5 - 1);
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
        k_binary_expand_op(b_x, b_mat, j, loop_ub - 2, ss, i6, i7 - 1);
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
        j_binary_expand_op(b_x, b_mat, j, b_loop_ub - 1, ss, i8, i9 - 1);
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
        i_binary_expand_op(b_x, b_mat, j, c_loop_ub, ss, i10, i11 - 1);
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
        i_binary_expand_op(b_x, b_mat, j, d_loop_ub - 1, ss, i12, i13 - 1);
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

      /* 'gkmPWM:282' [mm,dd] = max([sum(mat{j}==rs) sum(mat{j}(2:end)==rs(1:l-1)) sum(mat{j}(3:end)==rs(1:l-2)) sum(mat{j}(4:end)==rs(1:l-3)) sum(mat{j}(5:end)==rs(1:l-4)) sum(mat{j}(1:l-1)==rs(2:end)) sum(mat{j}(1:l-2)==rs(3:end)) sum(mat{j}(1:l-3)==rs(4:end)) sum(mat{j}(1:l-4)==rs(5:end))]); */
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
        o_binary_expand_op(b_x, b_mat, j, rs);
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
        n_binary_expand_op(b_x, b_mat, j, i, b_y - 1, rs, unnamed_idx_0_tmp_tmp
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
        m_binary_expand_op(b_x, b_mat, j, cur_idx, c_y - 1, rs, d_y - 1);
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
        l_binary_expand_op(b_x, b_mat, j, e_y, i1 - 1, rs, i2);
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
        l_binary_expand_op(b_x, b_mat, j, i3, i4 - 1, rs, i5 - 1);
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
        k_binary_expand_op(b_x, b_mat, j, loop_ub - 2, rs, i6, i7 - 1);
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
        j_binary_expand_op(b_x, b_mat, j, b_loop_ub - 1, rs, i8, i9 - 1);
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
        i_binary_expand_op(b_x, b_mat, j, c_loop_ub, rs, i10, i11 - 1);
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
        i_binary_expand_op(b_x, b_mat, j, d_loop_ub - 1, rs, i12, i13 - 1);
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

      /* 'gkmPWM:283' [M(j),ddd] = max([m mm]); */
      m[0] = curr_pos;
      m[1] = idx;
      c_maximum(m, &M_data[j], &cur_idx);

      /* 'gkmPWM:284' if ddd == 1 */
      if (cur_idx == 1) {
        /* 'gkmPWM:285' D(j) = d; */
        D_data[j] = iindx;

        /* 'gkmPWM:286' DD(j) = 1; */
        DD_data[j] = 1;
      } else {
        /* 'gkmPWM:287' else */
        /* 'gkmPWM:288' D(j) = dd; */
        D_data[j] = unnamed_idx_0_tmp_tmp;

        /* 'gkmPWM:289' DD(j) = 2; */
        DD_data[j] = 2;
      }
    }

    /* 'gkmPWM:292' if max(M) < k */
    if (maximum(M) < (double)x + 1.0) {
      /* 'gkmPWM:293' c = c+1; */
      (*c)++;

      /* 'gkmPWM:294' p{c-ikl} = s{i}; */
      i = sequences_data[(int)*c - 1].f1->size[0] * sequences_data[(int)*c - 1].
        f1->size[1];
      sequences_data[(int)*c - 1].f1->size[0] = 1;
      sequences_data[(int)*c - 1].f1->size[1] = s_data[b_i].f1->size[1];
      emxEnsureCapacity_char_T(sequences_data[(int)*c - 1].f1, i);
      loop_ub = s_data[b_i].f1->size[1];
      for (i = 0; i < loop_ub; i++) {
        sequences_data[(int)*c - 1].f1->data[i] = s_data[b_i].f1->data[i];
      }

      /* 'gkmPWM:295' mat{c} = ss; */
      i = mat_data[(int)*c - 1].f1->size[0] * mat_data[(int)*c - 1].f1->size[1];
      mat_data[(int)*c - 1].f1->size[0] = 1;
      mat_data[(int)*c - 1].f1->size[1] = ss->size[1];
      emxEnsureCapacity_real_T(mat_data[(int)*c - 1].f1, i);
      loop_ub = ss->size[1];
      for (i = 0; i < loop_ub; i++) {
        mat_data[(int)*c - 1].f1->data[i] = ss_data[i];
      }

      /* 'gkmPWM:296' ss = ss+1; */
      i = ss->size[0] * ss->size[1];
      ss->size[0] = 1;
      emxEnsureCapacity_real_T(ss, i);
      ss_data = ss->data;
      loop_ub = ss->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        ss_data[i]++;
      }

      /* 'gkmPWM:297' for j = 1:l */
      for (j = 0; j <= varargin_2; j++) {
        /* 'gkmPWM:298' pwms{c-ikl}(j,ss(j)) = pwms{c-ikl}(j,ss(j))+w(i); */
        i = (int)ss_data[j] - 1;
        pwms_data[(int)*c - 1].f1->data[j + pwms_data[(int)*c - 1].f1->size[0] *
          i] += alpha_data[b_i];
      }
    } else {
      /* 'gkmPWM:300' else */
      /* 'gkmPWM:301' [~,d] = max(M); */
      d_maximum(M, &curr_pos, &iindx);

      /* 'gkmPWM:302' if DD(d) == 1 && d > ikl */
      i = DD_data[iindx - 1];
      if (i == 1) {
        /* 'gkmPWM:303' ss = ss+1; */
        i = ss->size[0] * ss->size[1];
        ss->size[0] = 1;
        emxEnsureCapacity_real_T(ss, i);
        ss_data = ss->data;
        loop_ub = ss->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          ss_data[i]++;
        }

        /* 'gkmPWM:304' d = d-ikl; */
        /* 'gkmPWM:305' for j = 1:CC(D(d)) */
        i = D_data[iindx - 1] - 1;
        b_y = CC[i];
        for (j = 0; j < b_y; j++) {
          /* 'gkmPWM:306' pwms{d}(j+B(D(d)),ss(j+BB(D(d)))) = pwms{d}(j+B(D(d)),ss(j+BB(D(d))))+w(i); */
          unnamed_idx_0_tmp_tmp = (int)((unsigned int)j + B[i]);
          cur_idx = (int)ss_data[(int)((unsigned int)j + BB[i])] - 1;
          pwms_data[iindx - 1].f1->data[unnamed_idx_0_tmp_tmp + pwms_data[iindx
            - 1].f1->size[0] * cur_idx] += alpha_data[b_i];
        }
      } else if (i == 2) {
        /* 'gkmPWM:308' elseif DD(d) == 2 && d > ikl */
        /* 'gkmPWM:309' rs = rs+1; */
        i = rs->size[0] * rs->size[1];
        rs->size[0] = 1;
        emxEnsureCapacity_real_T(rs, i);
        rs_data = rs->data;
        loop_ub = rs->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          rs_data[i]++;
        }

        /* 'gkmPWM:310' d = d-ikl; */
        /* 'gkmPWM:311' for j = 1:CC(D(d)) */
        i = D_data[iindx - 1] - 1;
        b_y = CC[i];
        for (j = 0; j < b_y; j++) {
          /* 'gkmPWM:312' pwms{d}(j+B(D(d)),rs(j+BB(D(d)))) = pwms{d}(j+B(D(d)),rs(j+BB(D(d))))+w(i); */
          unnamed_idx_0_tmp_tmp = (int)((unsigned int)j + B[i]);
          cur_idx = (int)rs_data[(int)((unsigned int)j + BB[i])] - 1;
          pwms_data[iindx - 1].f1->data[unnamed_idx_0_tmp_tmp + pwms_data[iindx
            - 1].f1->size[0] * cur_idx] += alpha_data[b_i];
        }
      }
    }

    /* 'gkmPWM:316' if c == num+ikl */
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
  /* 'gkmPWM:321' new_mat = cell(c, 1); */
  /* 'gkmPWM:322' for cur_idx=1:c */
  i = (int)*c;
  b_y = mat->size[0];
  mat->size[0] = (int)*c;
  emxEnsureCapacity_cell_wrap_1(mat, b_y);
  b_mat_data = mat->data;
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'gkmPWM:323' new_mat{cur_idx} = mat{cur_idx}; */
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

  /* 'gkmPWM:325' mat = new_mat; */
  /*  p = p(1:c-ikl); */
  /* 'gkmPWM:328' p_len = c-ikl; */
  /* 'gkmPWM:329' new_p = cell(p_len, 1); */
  /* 'gkmPWM:330' for cur_idx=1:p_len */
  b_y = p->size[0] * p->size[1];
  p->size[0] = (int)*c;
  p->size[1] = 1;
  emxEnsureCapacity_cell_wrap_0(p, b_y);
  p_data = p->data;
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'gkmPWM:331' new_p{cur_idx} = p{cur_idx}; */
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

  /* 'gkmPWM:333' p = new_p; */
  /*  pwms = pwms(1:c-ikl); */
  /* 'gkmPWM:336' pwms_len = c-ikl; */
  /* 'gkmPWM:337' new_pwms = cell(pwms_len, 1); */
  /* 'gkmPWM:338' for cur_idx=1:pwms_len */
  b_y = pwms->size[0];
  pwms->size[0] = (int)*c;
  emxEnsureCapacity_cell_wrap_2(pwms, b_y);
  b_pwms_data = pwms->data;
  for (cur_idx = 0; cur_idx < i; cur_idx++) {
    /* 'gkmPWM:339' new_pwms{cur_idx} = pwms{cur_idx}; */
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

  /* 'gkmPWM:341' pwms = new_pwms; */
  /* 'gkmPWM:343' for i = 1:c-ikl */
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:344' for j = 1:l */
    for (j = 0; j <= varargin_2; j++) {
      /* 'gkmPWM:345' pwms{i}(j,:) = pwms{i}(j,:)/sum(pwms{i}(j,:)); */
      curr_pos = ((b_pwms_data[b_i].f1->data[j] + b_pwms_data[b_i].f1->data[j +
                   b_pwms_data[b_i].f1->size[0]]) + b_pwms_data[b_i].f1->data[j
                  + b_pwms_data[b_i].f1->size[0] * 2]) + b_pwms_data[b_i]
        .f1->data[j + b_pwms_data[b_i].f1->size[0] * 3];
      b_pwms_data[b_i].f1->data[j] /= curr_pos;
      b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0]] /= curr_pos;
      b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0] * 2] /=
        curr_pos;
      b_pwms_data[b_i].f1->data[j + b_pwms_data[b_i].f1->size[0] * 3] /=
        curr_pos;
    }
  }
}

static void u_binary_expand_op(emxArray_creal_T *vec, const emxArray_real_T *b,
  const emxArray_real_T *A, const emxArray_int8_T *x, const creal_T p3_data[])
{
  emxArray_creal_T *b_A;
  emxArray_creal_T *b_b;
  creal_T *b_A_data;
  creal_T *vec_data;
  const double *A_data;
  const double *b_data;
  double A_re_tmp;
  double b_A_re_tmp;
  double c_A_re_tmp;
  double d_A_re_tmp;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_1_0;
  const signed char *x_data;
  x_data = x->data;
  A_data = A->data;
  b_data = b->data;
  emxInit_creal_T(&b_A, 2);
  loop_ub = A->size[0];
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = loop_ub;
  b_A->size[1] = x->size[1];
  emxEnsureCapacity_creal_T(b_A, i);
  b_A_data = b_A->data;
  b_loop_ub = x->size[1];
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_A_data[i1 + b_A->size[0] * i].re = A_data[i1 + A->size[0] * (x_data[i] -
        1)];
      b_A_data[i1 + b_A->size[0] * i].im = 0.0;
    }
  }

  i = vec->size[0];
  vec->size[0] = b_A->size[0];
  emxEnsureCapacity_creal_T(vec, i);
  vec_data = vec->data;
  loop_ub = b_A->size[0];
  for (i = 0; i < loop_ub; i++) {
    vec_data[i].re = 0.0;
    vec_data[i].im = 0.0;
    b_loop_ub = b_A->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      A_re_tmp = b_A_data[i + b_A->size[0] * i1].re;
      b_A_re_tmp = p3_data[i1].im;
      c_A_re_tmp = b_A_data[i + b_A->size[0] * i1].im;
      d_A_re_tmp = p3_data[i1].re;
      vec_data[i].re += A_re_tmp * d_A_re_tmp - c_A_re_tmp * b_A_re_tmp;
      vec_data[i].im += A_re_tmp * b_A_re_tmp + c_A_re_tmp * d_A_re_tmp;
    }
  }

  emxFree_creal_T(&b_A);
  emxInit_creal_T(&b_b, 1);
  i = b_b->size[0];
  if (vec->size[0] == 1) {
    b_b->size[0] = b->size[0];
  } else {
    b_b->size[0] = vec->size[0];
  }

  emxEnsureCapacity_creal_T(b_b, i);
  b_A_data = b_b->data;
  b_loop_ub = (b->size[0] != 1);
  stride_1_0 = (vec->size[0] != 1);
  if (vec->size[0] == 1) {
    loop_ub = b->size[0];
  } else {
    loop_ub = vec->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    i1 = i * stride_1_0;
    b_A_data[i].re = b_data[i * b_loop_ub] - vec_data[i1].re;
    b_A_data[i].im = 0.0 - vec_data[i1].im;
  }

  i = vec->size[0];
  vec->size[0] = b_b->size[0];
  emxEnsureCapacity_creal_T(vec, i);
  vec_data = vec->data;
  loop_ub = b_b->size[0];
  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = b_A_data[i];
  }

  emxFree_creal_T(&b_b);
}

static void v_binary_expand_op(creal_T p3_data[], int *p3_size, const creal_T
  ps_data[], const int *ps_size, const creal_T E, const creal_T B_data[], const
  int B_size[2], const creal_T y)
{
  creal_T E_data[9];
  creal_T b_ps_data[3];
  double E_im;
  double E_re;
  double E_re_tmp;
  double b_E_re_tmp;
  double brm;
  int i;
  int loop_ub;
  loop_ub = B_size[0] * B_size[1];
  for (i = 0; i < loop_ub; i++) {
    E_re_tmp = B_data[i].im;
    b_E_re_tmp = B_data[i].re;
    E_re = E.re * b_E_re_tmp - E.im * E_re_tmp;
    E_im = E.re * E_re_tmp + E.im * b_E_re_tmp;
    if (y.im == 0.0) {
      if (E_im == 0.0) {
        E_data[i].re = E_re / y.re;
        E_data[i].im = 0.0;
      } else if (E_re == 0.0) {
        E_data[i].re = 0.0;
        E_data[i].im = E_im / y.re;
      } else {
        E_data[i].re = E_re / y.re;
        E_data[i].im = E_im / y.re;
      }
    } else if (y.re == 0.0) {
      if (E_re == 0.0) {
        E_data[i].re = E_im / y.im;
        E_data[i].im = 0.0;
      } else if (E_im == 0.0) {
        E_data[i].re = 0.0;
        E_data[i].im = -(E_re / y.im);
      } else {
        E_data[i].re = E_im / y.im;
        E_data[i].im = -(E_re / y.im);
      }
    } else {
      brm = fabs(y.re);
      E_re_tmp = fabs(y.im);
      if (brm > E_re_tmp) {
        E_re_tmp = y.im / y.re;
        b_E_re_tmp = y.re + E_re_tmp * y.im;
        E_data[i].re = (E_re + E_re_tmp * E_im) / b_E_re_tmp;
        E_data[i].im = (E_im - E_re_tmp * E_re) / b_E_re_tmp;
      } else if (E_re_tmp == brm) {
        if (y.re > 0.0) {
          E_re_tmp = 0.5;
        } else {
          E_re_tmp = -0.5;
        }

        if (y.im > 0.0) {
          b_E_re_tmp = 0.5;
        } else {
          b_E_re_tmp = -0.5;
        }

        E_data[i].re = (E_re * E_re_tmp + E_im * b_E_re_tmp) / brm;
        E_data[i].im = (E_im * E_re_tmp - E_re * b_E_re_tmp) / brm;
      } else {
        E_re_tmp = y.re / y.im;
        b_E_re_tmp = y.im + E_re_tmp * y.re;
        E_data[i].re = (E_re_tmp * E_re + E_im) / b_E_re_tmp;
        E_data[i].im = (E_re_tmp * E_im - E_re) / b_E_re_tmp;
      }
    }
  }

  *p3_size = 3;
  for (i = 0; i < 3; i++) {
    p3_data[i].re = 0.0;
    p3_data[i].im = 0.0;
    p3_data[i].re += E_data[i].re;
    p3_data[i].im += E_data[i].im;
    p3_data[i].re += E_data[i + 3].re;
    p3_data[i].im += E_data[i + 3].im;
    p3_data[i].re += E_data[i + 6].re;
    p3_data[i].im += E_data[i + 6].im;
  }

  loop_ub = (*ps_size != 1);
  b_ps_data[0].re = ps_data[0].re + p3_data[0].re;
  b_ps_data[0].im = ps_data[0].im + p3_data[0].im;
  b_ps_data[1].re = ps_data[loop_ub].re + p3_data[1].re;
  b_ps_data[1].im = ps_data[loop_ub].im + p3_data[1].im;
  i = loop_ub << 1;
  b_ps_data[2].re = ps_data[i].re + p3_data[2].re;
  b_ps_data[2].im = ps_data[i].im + p3_data[2].im;
  *p3_size = 3;
  p3_data[0] = b_ps_data[0];
  p3_data[1] = b_ps_data[1];
  p3_data[2] = b_ps_data[2];
}

static void w_binary_expand_op(emxArray_creal_T *vec, const emxArray_real_T *b,
  const emxArray_real_T *A, const signed char Posvec_data[], const int
  Posvec_size[2], const creal_T p2_data[])
{
  emxArray_creal_T *b_A;
  emxArray_creal_T *b_b;
  creal_T *b_A_data;
  creal_T *vec_data;
  const double *A_data;
  const double *b_data;
  double A_re_tmp;
  double b_A_re_tmp;
  double c_A_re_tmp;
  double d_A_re_tmp;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_1_0;
  A_data = A->data;
  b_data = b->data;
  emxInit_creal_T(&b_A, 2);
  loop_ub = A->size[0];
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = loop_ub;
  b_A->size[1] = Posvec_size[1];
  emxEnsureCapacity_creal_T(b_A, i);
  b_A_data = b_A->data;
  b_loop_ub = Posvec_size[1];
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_A_data[i1 + b_A->size[0] * i].re = A_data[i1 + A->size[0] *
        (Posvec_data[i] - 1)];
      b_A_data[i1 + b_A->size[0] * i].im = 0.0;
    }
  }

  i = vec->size[0];
  vec->size[0] = b_A->size[0];
  emxEnsureCapacity_creal_T(vec, i);
  vec_data = vec->data;
  loop_ub = b_A->size[0];
  for (i = 0; i < loop_ub; i++) {
    vec_data[i].re = 0.0;
    vec_data[i].im = 0.0;
    b_loop_ub = b_A->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      A_re_tmp = b_A_data[i + b_A->size[0] * i1].re;
      b_A_re_tmp = p2_data[i1].im;
      c_A_re_tmp = b_A_data[i + b_A->size[0] * i1].im;
      d_A_re_tmp = p2_data[i1].re;
      vec_data[i].re += A_re_tmp * d_A_re_tmp - c_A_re_tmp * b_A_re_tmp;
      vec_data[i].im += A_re_tmp * b_A_re_tmp + c_A_re_tmp * d_A_re_tmp;
    }
  }

  emxFree_creal_T(&b_A);
  emxInit_creal_T(&b_b, 1);
  i = b_b->size[0];
  if (vec->size[0] == 1) {
    b_b->size[0] = b->size[0];
  } else {
    b_b->size[0] = vec->size[0];
  }

  emxEnsureCapacity_creal_T(b_b, i);
  b_A_data = b_b->data;
  b_loop_ub = (b->size[0] != 1);
  stride_1_0 = (vec->size[0] != 1);
  if (vec->size[0] == 1) {
    loop_ub = b->size[0];
  } else {
    loop_ub = vec->size[0];
  }

  for (i = 0; i < loop_ub; i++) {
    i1 = i * stride_1_0;
    b_A_data[i].re = b_data[i * b_loop_ub] - vec_data[i1].re;
    b_A_data[i].im = 0.0 - vec_data[i1].im;
  }

  i = vec->size[0];
  vec->size[0] = b_b->size[0];
  emxEnsureCapacity_creal_T(vec, i);
  vec_data = vec->data;
  loop_ub = b_b->size[0];
  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = b_A_data[i];
  }

  emxFree_creal_T(&b_b);
}

static void x_binary_expand_op(creal_T p[4], const double MAT_data[], const int
  MAT_size[2], const creal_T y, const signed char b[4])
{
  creal_T MAT[4];
  creal_T dcv[4];
  int i;
  int stride_0_0;
  stride_0_0 = (MAT_size[0] != 1);
  MAT[0].re = MAT_data[0] - y.re * (double)b[0];
  MAT[0].im = 0.0 - y.im * (double)b[0];
  MAT[1].re = MAT_data[stride_0_0] - y.re * (double)b[1];
  MAT[1].im = 0.0 - y.im * (double)b[1];
  i = MAT_size[0] * (MAT_size[1] != 1);
  MAT[2].re = MAT_data[i] - y.re * (double)b[2];
  MAT[2].im = 0.0 - y.im * (double)b[2];
  MAT[3].re = MAT_data[stride_0_0 + i] - y.re * (double)b[3];
  MAT[3].im = 0.0 - y.im * (double)b[3];
  e_mpower(MAT, dcv);
  memcpy(&p[0], &dcv[0], 4U * sizeof(creal_T));
}

static void y_binary_expand_op(emxArray_real_T *b, const emxArray_real_T *A, int
  iindx)
{
  emxArray_real_T *b_b;
  const double *A_data;
  double *b_b_data;
  double *b_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  A_data = A->data;
  b_data = b->data;
  emxInit_real_T(&b_b, 1);
  i = A->size[0];
  stride_0_0 = b_b->size[0];
  if (i == 1) {
    b_b->size[0] = b->size[0];
  } else {
    b_b->size[0] = i;
  }

  emxEnsureCapacity_real_T(b_b, stride_0_0);
  b_b_data = b_b->data;
  stride_0_0 = (b->size[0] != 1);
  stride_1_0 = (i != 1);
  if (i == 1) {
    loop_ub = b->size[0];
  } else {
    loop_ub = i;
  }

  for (i = 0; i < loop_ub; i++) {
    b_b_data[i] = b_data[i * stride_0_0] - A_data[i * stride_1_0 + A->size[0] *
      (iindx - 1)];
  }

  i = b->size[0];
  b->size[0] = b_b->size[0];
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  loop_ub = b_b->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_data[i] = b_b_data[i];
  }

  emxFree_real_T(&b_b);
}

/*
 * function gkmPWM(varargin)
 */
void gkmPWM(const emxArray_char_T *varargin_1, const emxArray_char_T *varargin_2,
            const emxArray_char_T *varargin_3, double varargin_4, double
            varargin_5, double varargin_6, double varargin_7, double varargin_8,
            double varargin_9, double varargin_10, double varargin_11, double
            varargin_12)
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
  emxArray_real_T *E;
  emxArray_real_T *a__2;
  emxArray_real_T *a__3;
  emxArray_real_T *b_comb;
  emxArray_real_T *b_mat;
  emxArray_real_T *c2;
  emxArray_real_T *c_mat;
  emxArray_real_T *cfile;
  emxArray_real_T *comb;
  emxArray_real_T *dc;
  emxArray_real_T *dc2;
  emxArray_real_T *negvec;
  emxArray_real_T *seqvec;
  emxArray_real_T *seqvec2;
  double mat[16];
  double mat2[16];
  double GCmat[4];
  double startvec[4];
  double GCmat_tmp;
  double GCneg1;
  double GCpos1;
  double c;
  double pnr;
  double rcnum;
  double tmp;
  double xtmp;
  double *c2_data;
  double *cfile_data;
  double *comb_data;
  double *dc2_data;
  double *dc_data;
  double *mat_data;
  double *negvec_data;
  double *seqvec2_data;
  double *seqvec_data;
  unsigned int a;
  int b_i;
  int b_j1;
  int i;
  int j2;
  int loop_ub;
  int m;
  int nd2;
  int ntilerows;
  int outsize_idx_0;
  int u1;
  const char *varargin_1_data;
  signed char fileid;
  char *b_varargin_1_data;
  char *charStr_data;
  bool ipnr;
  if (!isInitialized_gkmPWM) {
    gkmPWM_initialize();
  }

  varargin_1_data = varargin_1->data;

  /*  usage  gkmPWM(fileheader, wfile, mnum, num) */
  /*  fileheader: everything before *svseq.fa */
  /*  wfile: a list of kmer weights */
  /*  memefile: memefile to identify motifs */
  /*  mmun: the number of motifs to identify */
  /*  num: number of iterations (I recommend at least 10*mnum) */
  /* 'gkmPWM:8' if nargin < 4 */
  /* 'gkmPWM:12' fileheader = varargin{1}; */
  /* 'gkmPWM:13' wfile = varargin{2}; */
  /* 'gkmPWM:14' memefile = varargin{3}; */
  /* 'gkmPWM:15' mnum = varargin{4}; */
  /* 'gkmPWM:16' num = varargin{5}; */
  /* 'gkmPWM:17' rcorr = varargin{6}; */
  /* 'gkmPWM:18' reg = varargin{7}; */
  /* 'gkmPWM:19' l_svm = varargin{8}; */
  /* 'gkmPWM:20' k_svm = varargin{9}; */
  /* 'gkmPWM:21' BG_GC = varargin{10}; */
  /* 'gkmPWM:22' RC = varargin{11}; */
  /* 'gkmPWM:23' pnr = varargin{12}; */
  pnr = varargin_12;

  /* 'gkmPWM:25' if pnr == 0 */
  if (varargin_12 == 0.0) {
    /* 'gkmPWM:26' ipnr = true; */
    ipnr = true;
  } else {
    /* 'gkmPWM:27' else */
    /* 'gkmPWM:28' ipnr = false; */
    ipnr = false;
  }

  emxInit_real_T(&a__2, 1);
  emxInit_real_T(&a__3, 1);
  emxInit_real_T(&comb, 2);
  emxInit_real_T(&c2, 2);
  emxInit_real_T(&seqvec2, 2);

  /*  if nargin > 4 */
  /*      f = find(strcmp('MaxCorr', varargin)); */
  /*      if ~isempty(f); */
  /*          rcorr = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('PNratio', varargin)); */
  /*      if ~isempty(f); */
  /*          pnr = varargin{f+1}; */
  /*          ipnr=false; */
  /*      else */
  /*          ipnr=true; */
  /*      end */
  /*      f = find(strcmp('RegFrac', varargin)); */
  /*      if ~isempty(f); */
  /*          reg = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('l', varargin)); */
  /*      if ~isempty(f); */
  /*          l_svm = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('k', varargin)); */
  /*      if ~isempty(f); */
  /*          k_svm = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('Mode', varargin)); */
  /*      if ~isempty(f) && strcmp('Compare',varargin{f+1}); */
  /*          BG_GC = 1; */
  /*      end */
  /*      f = find(strcmp('RC', varargin)); */
  /*      if ~isempty(f); */
  /*          RC = varargin{f+1}; */
  /*      end */
  /*  end */
  /* 'gkmPWM:65' [comb,~,~,~,~,rcnum] = genIndex(l_svm,k_svm); */
  genIndex(varargin_8, varargin_9, comb, c2, a__2, a__3, seqvec2, &rcnum);
  mat_data = comb->data;

  /* generate gapped positions, adjusted for reverse complements */
  /* 'gkmPWM:67' if length(comb)*4^k_svm > 10^6 */
  if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
    u1 = 0;
  } else {
    m = comb->size[0];
    u1 = comb->size[1];
    if (m >= u1) {
      u1 = m;
    }
  }

  c = pow(4.0, varargin_9);
  if ((double)u1 * c > 1.0E+6) {
    /*  error([num2str(length(comb)*4^k_svm) ' exceeds the maximum number of gapped kmers allowed (10^6)']) */
    /* 'gkmPWM:69' fprintf('The choice of l and k combination exceeds the maximum number of gapped kmers allowed (10^6)\n'); */
    printf("The choice of l and k combination exceeds the maximum number of gapped kmers allowed (10^6)\n");
    fflush(stdout);

    /* 'gkmPWM:70' fprintf('l = 10, k = 6\n'); */
    printf("l = 10, k = 6\n");
    fflush(stdout);
  }

  emxInit_char_T(&b_varargin_1, 2);

  /*  disp(['Running gkmPWM on ' fileheader ' for ' num2str(mnum) ' motifs and ' num2str(num) ' iterations']) */
  /* 'gkmPWM:74' fprintf('Running gkmPWM on %s for %d motifs and %d iterations\n', fileheader, int32(mnum), int32(num)); */
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

  /* 'gkmPWM:75' fprintf('Counting gapped k-mers\n'); */
  printf("Counting gapped k-mers\n");
  fflush(stdout);

  /* 'gkmPWM:77' [A, GCpos1, GCneg1,mat,mat2] = getgkmcounts(fileheader,l_svm,k_svm,1,RC); */
  getgkmcounts(varargin_1, varargin_8, varargin_9, varargin_11, cfile, &GCpos1,
               &GCneg1, mat, mat2);
  cfile_data = cfile->data;

  /* count gapped k-mers, get GC content, and dinucleotide distribution */
  /* 'gkmPWM:78' if BG_GC == 1 */
  if (varargin_10 == 1.0) {
    /* 'gkmPWM:79' mat = (mat+mat2)/2; */
    for (i = 0; i < 16; i++) {
      mat[i] = (mat[i] + mat2[i]) / 2.0;
    }

    /* 'gkmPWM:80' GCpos1 = (GCpos1+GCneg1)/2; */
    GCneg1 = (GCpos1 + GCneg1) / 2.0;

    /* 'gkmPWM:81' GCneg1 = GCpos1; */
  }

  /* 'gkmPWM:83' negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC); */
  /* 'BGkmer:2' len = length(c); */
  if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
    u1 = 0;
  } else {
    m = comb->size[0];
    u1 = comb->size[1];
    if (m >= u1) {
      u1 = m;
    }
  }

  emxInit_real_T(&negvec, 1);

  /* 'BGkmer:3' alen = len-rcnum; */
  /* 'BGkmer:4' negweights = zeros(len*4^k,1); */
  loop_ub = (int)((double)u1 * c);
  i = negvec->size[0];
  negvec->size[0] = loop_ub;
  emxEnsureCapacity_real_T(negvec, i);
  negvec_data = negvec->data;
  for (i = 0; i < loop_ub; i++) {
    negvec_data[i] = 0.0;
  }

  emxInit_real_T(&dc2, 2);
  emxInit_cell_wrap_13(&p);

  /* 'BGkmer:5' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  GCmat_tmp = 0.5 - GCneg1 / 2.0;
  GCmat[1] = GCneg1 / 2.0;
  GCmat[2] = GCneg1 / 2.0;

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

  p_data[0].f1[0] = 1.0;
  p_data[0].f1[5] = 1.0;
  p_data[0].f1[10] = 1.0;
  p_data[0].f1[15] = 1.0;

  /* 'BGkmer:13' for i = 1:l-1 */
  i = (int)(varargin_8 - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'BGkmer:14' p{i+1} = p{i}*mat; */
    for (j2 = 0; j2 < 4; j2++) {
      for (nd2 = 0; nd2 < 4; nd2++) {
        m = nd2 << 2;
        mat2[j2 + m] = ((p_data[b_i].f1[j2] * mat[m] + p_data[b_i].f1[j2 + 4] *
                         mat[m + 1]) + p_data[b_i].f1[j2 + 8] * mat[m + 2]) +
          p_data[b_i].f1[j2 + 12] * mat[m + 3];
      }
    }

    for (j2 = 0; j2 < 16; j2++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[j2] = mat2[j2];
    }
  }

  emxInit_real_T(&seqvec, 2);

  /* 'BGkmer:16' seqvec = zeros(4^k, k); */
  i = seqvec->size[0] * seqvec->size[1];
  j2 = (int)c;
  seqvec->size[0] = (int)c;
  nd2 = (int)varargin_9;
  seqvec->size[1] = (int)varargin_9;
  emxEnsureCapacity_real_T(seqvec, i);
  seqvec_data = seqvec->data;
  loop_ub = (int)c * (int)varargin_9;
  for (i = 0; i < loop_ub; i++) {
    seqvec_data[i] = 0.0;
  }

  /* 'BGkmer:17' for i = 1:k */
  for (b_i = 0; b_i < nd2; b_i++) {
    /* 'BGkmer:18' for j = 1:4^k */
    for (m = 0; m < j2; m++) {
      /* 'BGkmer:19' seqvec(j,i) = mod(floor((j-1)/4^(i-1)), 4)+1; */
      xtmp = floor((((double)m + 1.0) - 1.0) / pow(4.0, ((double)b_i + 1.0) -
        1.0));
      if (xtmp == 0.0) {
        xtmp = 0.0;
      } else {
        xtmp = fmod(xtmp, 4.0);
        if (xtmp == 0.0) {
          xtmp = 0.0;
        }
      }

      seqvec_data[m + seqvec->size[0] * b_i] = xtmp + 1.0;
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
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (seqvec->size[1] - b_j1) - 1;
      for (b_i = 0; b_i < m; b_i++) {
        xtmp = seqvec2_data[b_i + seqvec2->size[0] * b_j1];
        seqvec2_data[b_i + seqvec2->size[0] * b_j1] = seqvec2_data[b_i +
          seqvec2->size[0] * j2];
        seqvec2_data[b_i + seqvec2->size[0] * j2] = xtmp;
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
      c2_data[i] = mat_data[i];
    }

    m = comb->size[0];
    nd2 = comb->size[1] >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (comb->size[1] - b_j1) - 1;
      for (b_i = 0; b_i < m; b_i++) {
        xtmp = c2_data[b_i + c2->size[0] * b_j1];
        c2_data[b_i + c2->size[0] * b_j1] = c2_data[b_i + c2->size[0] * j2];
        c2_data[b_i + c2->size[0] * j2] = xtmp;
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
  emxInit_real_T(&dc, 2);
  emxInit_real_T(&b_comb, 2);
  for (b_i = 0; b_i < u1; b_i++) {
    /* 'BGkmer:28' dc = diff(c(i,:)); */
    loop_ub = comb->size[1];
    i = b_comb->size[0] * b_comb->size[1];
    b_comb->size[0] = 1;
    b_comb->size[1] = comb->size[1];
    emxEnsureCapacity_real_T(b_comb, i);
    comb_data = b_comb->data;
    for (i = 0; i < loop_ub; i++) {
      comb_data[i] = mat_data[b_i + comb->size[0] * i];
    }

    diff(b_comb, dc);
    dc_data = dc->data;

    /* 'BGkmer:29' startvec = GCmat*p{c(i,1)}; */
    for (i = 0; i < 16; i++) {
      mat2[i] = p_data[(int)mat_data[b_i] - 1].f1[i];
    }

    GCpos1 = GCmat[1];
    xtmp = GCmat[2];
    for (i = 0; i < 4; i++) {
      j2 = i << 2;
      startvec[i] = ((GCmat_tmp * mat2[j2] + GCpos1 * mat2[j2 + 1]) + xtmp *
                     mat2[j2 + 2]) + GCmat_tmp * mat2[j2 + 3];
    }

    /* 'BGkmer:30' if RC */
    if (varargin_11 != 0.0) {
      /* 'BGkmer:31' dc2 = diff(c2(i,:)); */
      loop_ub = c2->size[1];
      i = b_comb->size[0] * b_comb->size[1];
      b_comb->size[0] = 1;
      b_comb->size[1] = c2->size[1];
      emxEnsureCapacity_real_T(b_comb, i);
      comb_data = b_comb->data;
      for (i = 0; i < loop_ub; i++) {
        comb_data[i] = c2_data[b_i + c2->size[0] * i];
      }

      diff(b_comb, dc2);
      dc2_data = dc2->data;

      /* 'BGkmer:32' startvec2 = GCmat*p{c2(i,1)}; */
    }

    /* 'BGkmer:34' for ii = 1:4^k */
    i = (int)pow(4.0, varargin_9);
    for (m = 0; m < i; m++) {
      /* 'BGkmer:35' negweights(a) = startvec(seqvec(ii,1)); */
      b_j1 = (int)(a + m) - 1;
      negvec_data[b_j1] = startvec[(int)seqvec_data[m] - 1];

      /* 'BGkmer:36' if RC */
      if (varargin_11 != 0.0) {
        /* 'BGkmer:37' tmp = startvec(seqvec2(ii,1)); */
        tmp = startvec[(int)seqvec2_data[m] - 1];
      }

      /* 'BGkmer:39' for iii = 1:k-1 */
      j2 = (int)(varargin_9 - 1.0);
      for (nd2 = 0; nd2 < j2; nd2++) {
        /* 'BGkmer:40' matt = p{dc(iii)+1}; */
        /* 'BGkmer:41' negweights(a) = negweights(a)*matt(seqvec(ii,iii), seqvec(ii,iii+1)); */
        negvec_data[b_j1] *= p_data[(int)(dc_data[nd2] + 1.0) - 1].f1[((int)
          seqvec_data[m + seqvec->size[0] * nd2] + (((int)seqvec_data[m +
          seqvec->size[0] * (nd2 + 1)] - 1) << 2)) - 1];

        /* 'BGkmer:42' if RC */
        if (varargin_11 != 0.0) {
          /* 'BGkmer:43' matt2 = p{dc2(iii)+1}; */
          /* 'BGkmer:44' tmp = tmp*matt2(seqvec2(ii,iii), seqvec2(ii,iii+1)); */
          tmp *= p_data[(int)(dc2_data[nd2] + 1.0) - 1].f1[((int)seqvec2_data[m
            + seqvec2->size[0] * nd2] + (((int)seqvec2_data[m + seqvec2->size[0]
            * (nd2 + 1)] - 1) << 2)) - 1];
        }
      }

      /* 'BGkmer:47' if RC */
      if (varargin_11 != 0.0) {
        /* 'BGkmer:48' negweights(a) = (negweights(a)+tmp); */
        negvec_data[b_j1] += tmp;
      }

      /* 'BGkmer:50' a = a+1; */
    }

    a += i;
  }

  emxFree_cell_wrap_13(&p);
  emxFree_real_T(&seqvec);
  emxFree_real_T(&seqvec2);
  emxFree_real_T(&c2);
  emxFree_real_T(&comb);

  /* 'BGkmer:53' if rcnum > 0 && RC */
  if ((rcnum > 0.0) && (varargin_11 != 0.0)) {
    /* 'BGkmer:54' negweights(4^k*alen+1:end) = negweights(4^k*alen+1:end)/sqrt(2); */
    GCpos1 = c * ((double)u1 - rcnum) + 1.0;
    if (GCpos1 > negvec->size[0]) {
      i = 1;
      j2 = -1;
      nd2 = 0;
    } else {
      i = (int)GCpos1;
      j2 = (int)GCpos1 - 2;
      nd2 = negvec->size[0];
    }

    m = (nd2 - j2) - 1;
    nd2 = b_comb->size[0] * b_comb->size[1];
    b_comb->size[0] = 1;
    b_comb->size[1] = m;
    emxEnsureCapacity_real_T(b_comb, nd2);
    comb_data = b_comb->data;
    for (nd2 = 0; nd2 < m; nd2++) {
      comb_data[nd2] = negvec_data[(i + nd2) - 1] / 1.4142135623730951;
    }

    loop_ub = b_comb->size[1];
    for (i = 0; i < loop_ub; i++) {
      negvec_data[(j2 + i) + 1] = comb_data[i];
    }
  }

  emxFree_real_T(&b_comb);

  /* generate expected gapped kmer distribution of background */
  /* 'gkmPWM:84' GC=[0.5-GCneg1/2 GCneg1/2 GCneg1/2 0.5-GCneg1/2]; */
  GCmat[0] = GCmat_tmp;
  GCmat[1] = GCneg1 / 2.0;
  GCmat[2] = GCneg1 / 2.0;
  GCmat[3] = GCmat_tmp;

  /* GC content vector */
  /* 'gkmPWM:85' if ipnr */
  if (ipnr) {
    /* 'gkmPWM:86' pnr = abs(min(A)/max(A)); */
    pnr = fabs(minimum(cfile) / maximum(cfile));
  }

  emxInit_cell_wrap_0(&kmers, 2);
  emxInit_cell_wrap_2(&b_p);
  emxInit_cell_wrap_1(&seed);

  /* 'gkmPWM:88' fprintf('Finding PWM seeds\n'); */
  printf("Finding PWM seeds\n");
  fflush(stdout);

  /* get PWM seeds using the kmer weight vectors */
  /* 'gkmPWM:90' coder.varsize('p'); */
  /* 'gkmPWM:91' [kmers, seed,p,c] = seed_kmers(wfile, mnum,'descend', {}); */
  seed_kmers(varargin_2, varargin_4, kmers, seed, b_p, &xtmp);
  b_p_data = b_p->data;
  kmers_data = kmers->data;

  /* 'gkmPWM:92' if pnr ~= 0 */
  if (pnr != 0.0) {
    emxInit_cell_wrap_0(&kmers2, 1);
    emxInit_cell_wrap_1(&seed2);
    emxInit_cell_wrap_2(&pp);

    /* 'gkmPWM:93' [kmers2, seed2,pp, c2] = seed_kmers(wfile, max([floor(mnum/pnr) 2]),'ascend',seed); */
    GCpos1 = floor(varargin_4 / pnr);
    if (GCpos1 < 2.0) {
      GCpos1 = 2.0;
    }

    b_seed_kmers(varargin_2, GCpos1, seed, kmers2, seed2, pp, &xtmp);
    pp_data = pp->data;
    kmers2_data = kmers2->data;

    /*      kmers = [kmers;kmers2]; */
    /*      p = [p;pp]; */
    /* 'gkmPWM:96' for cur_idx=1:length(kmers2) */
    i = kmers2->size[0];
    emxFree_cell_wrap_1(&seed2);
    for (nd2 = 0; nd2 < i; nd2++) {
      /* 'gkmPWM:97' kmers{end + 1} = kmers2{cur_idx}; */
      m = kmers->size[0] + 1;
      j2 = kmers->size[0] * kmers->size[1];
      kmers->size[0]++;
      kmers->size[1] = 1;
      emxEnsureCapacity_cell_wrap_0(kmers, j2);
      kmers_data = kmers->data;
      loop_ub = kmers2_data[nd2].f1->size[1] - 1;
      j2 = kmers_data[m - 1].f1->size[0] * kmers_data[m - 1].f1->size[1];
      kmers_data[m - 1].f1->size[0] = 1;
      emxEnsureCapacity_char_T(kmers_data[m - 1].f1, j2);
      j2 = kmers_data[kmers->size[0] - 1].f1->size[0] * kmers_data[kmers->size[0]
        - 1].f1->size[1];
      kmers_data[kmers->size[0] - 1].f1->size[1] = kmers2_data[nd2].f1->size[1];
      emxEnsureCapacity_char_T(kmers_data[kmers->size[0] - 1].f1, j2);
      m = kmers->size[0] - 1;
      for (j2 = 0; j2 <= loop_ub; j2++) {
        kmers_data[m].f1->data[j2] = kmers2_data[nd2].f1->data[j2];
      }
    }

    emxFree_cell_wrap_0(&kmers2);

    /* 'gkmPWM:99' for cur_idx=1:length(pp) */
    i = pp->size[0];
    for (nd2 = 0; nd2 < i; nd2++) {
      /* 'gkmPWM:100' p{end + 1} = pp{cur_idx}; */
      m = b_p->size[0] + 1;
      j2 = b_p->size[0];
      b_p->size[0]++;
      emxEnsureCapacity_cell_wrap_2(b_p, j2);
      b_p_data = b_p->data;
      loop_ub = pp_data[nd2].f1->size[0] * 4;
      j2 = b_p_data[m - 1].f1->size[0] * b_p_data[m - 1].f1->size[1];
      b_p_data[m - 1].f1->size[0] = pp_data[nd2].f1->size[0];
      emxEnsureCapacity_real_T(b_p_data[m - 1].f1, j2);
      j2 = b_p_data[b_p->size[0] - 1].f1->size[0] * b_p_data[b_p->size[0] - 1].
        f1->size[1];
      b_p_data[b_p->size[0] - 1].f1->size[1] = 4;
      emxEnsureCapacity_real_T(b_p_data[b_p->size[0] - 1].f1, j2);
      m = b_p->size[0] - 1;
      for (j2 = 0; j2 < loop_ub; j2++) {
        b_p_data[m].f1->data[j2] = pp_data[nd2].f1->data[j2];
      }
    }

    emxFree_cell_wrap_2(&pp);

    /* 'gkmPWM:102' tot = c2; */
  } else {
    /* 'gkmPWM:103' else */
    /* 'gkmPWM:104' tot = c; */
  }

  emxFree_cell_wrap_1(&seed);

  /* 'gkmPWM:106' fprintf('Seeding %d PWMs at the following kmers\n', int32(tot)); */
  printf("Seeding %d PWMs at the following kmers\n", (int)rt_roundd(xtmp));
  fflush(stdout);

  /* 'gkmPWM:107' for i = 1:tot */
  i = (int)xtmp;
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:108' fprintf("%s\n", kmers{i}); */
    j2 = b_varargin_1->size[0] * b_varargin_1->size[1];
    b_varargin_1->size[0] = 1;
    b_varargin_1->size[1] = kmers_data[b_i].f1->size[1] + 1;
    emxEnsureCapacity_char_T(b_varargin_1, j2);
    b_varargin_1_data = b_varargin_1->data;
    loop_ub = kmers_data[b_i].f1->size[1];
    for (j2 = 0; j2 < loop_ub; j2++) {
      b_varargin_1_data[j2] = kmers_data[b_i].f1->data[j2];
    }

    b_varargin_1_data[kmers_data[b_i].f1->size[1]] = '\x00';
    printf("%s\n", &b_varargin_1_data[0]);
    fflush(stdout);
  }

  emxFree_cell_wrap_0(&kmers);

  /* 'gkmPWM:111' for i = 1:length(p) */
  i = b_p->size[0];
  if (0 <= b_p->size[0] - 1) {
    outsize_idx_0 = (int)(varargin_8 + 1.0);
    ntilerows = (int)(varargin_8 + 1.0);
  }

  emxInit_real_T(&b_mat, 2);
  emxInit_real_T(&c_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWM:112' p{i} = extendPWM(p{i},l_svm+1,GC); */
    /* 'gkmPWM:1002' mat = repmat(GCmat, n,1); */
    j2 = b_mat->size[0] * b_mat->size[1];
    b_mat->size[0] = outsize_idx_0;
    b_mat->size[1] = 4;
    emxEnsureCapacity_real_T(b_mat, j2);
    mat_data = b_mat->data;
    for (nd2 = 0; nd2 < 4; nd2++) {
      m = nd2 * ntilerows;
      for (j2 = 0; j2 < ntilerows; j2++) {
        mat_data[m + j2] = GCmat[nd2];
      }
    }

    /* 'gkmPWM:1003' ext_pwm = [mat;pwm;mat]; */
    j2 = c_mat->size[0] * c_mat->size[1];
    c_mat->size[0] = (b_mat->size[0] + b_p_data[b_i].f1->size[0]) + b_mat->size
      [0];
    c_mat->size[1] = 4;
    emxEnsureCapacity_real_T(c_mat, j2);
    c2_data = c_mat->data;
    loop_ub = b_mat->size[0];
    for (j2 = 0; j2 < 4; j2++) {
      for (nd2 = 0; nd2 < loop_ub; nd2++) {
        c2_data[nd2 + c_mat->size[0] * j2] = mat_data[nd2 + b_mat->size[0] * j2];
      }
    }

    loop_ub = b_p_data[b_i].f1->size[0];
    for (j2 = 0; j2 < 4; j2++) {
      for (nd2 = 0; nd2 < loop_ub; nd2++) {
        c2_data[(nd2 + b_mat->size[0]) + c_mat->size[0] * j2] = b_p_data[b_i].
          f1->data[nd2 + b_p_data[b_i].f1->size[0] * j2];
      }
    }

    loop_ub = b_mat->size[0];
    for (j2 = 0; j2 < 4; j2++) {
      for (nd2 = 0; nd2 < loop_ub; nd2++) {
        c2_data[((nd2 + b_mat->size[0]) + b_p_data[b_i].f1->size[0]) +
          c_mat->size[0] * j2] = mat_data[nd2 + b_mat->size[0] * j2];
      }
    }

    j2 = b_p_data[b_i].f1->size[0] * b_p_data[b_i].f1->size[1];
    b_p_data[b_i].f1->size[0] = c_mat->size[0];
    b_p_data[b_i].f1->size[1] = 4;
    emxEnsureCapacity_real_T(b_p_data[b_i].f1, j2);
    loop_ub = c_mat->size[0] * 4;
    for (j2 = 0; j2 < loop_ub; j2++) {
      b_p_data[b_i].f1->data[j2] = c2_data[j2];
    }
  }

  emxFree_real_T(&c_mat);
  emxFree_real_T(&b_mat);

  /* 'gkmPWM:115' fprintf('Running de novo motif discovery\n'); */
  printf("Running de novo motif discovery\n");
  fflush(stdout);

  /* 'gkmPWM:116' m = mean(A); */
  /* 'gkmPWM:117' s = std(A); */
  /* 'gkmPWM:118' cfile = A-negvec/sum(negvec)*sum(A); */
  GCpos1 = blockedSummation(cfile, cfile->size[0]);
  xtmp = blockedSummation(negvec, negvec->size[0]);
  if (cfile->size[0] == negvec->size[0]) {
    loop_ub = cfile->size[0];
    for (i = 0; i < loop_ub; i++) {
      cfile_data[i] -= negvec_data[i] / xtmp * GCpos1;
    }
  } else {
    binary_expand_op(cfile, negvec, xtmp, GCpos1);
    cfile_data = cfile->data;
  }

  /* 'gkmPWM:119' cfile = cfile/max(abs(cfile)); */
  m = cfile->size[0];
  i = a__2->size[0];
  a__2->size[0] = cfile->size[0];
  emxEnsureCapacity_real_T(a__2, i);
  mat_data = a__2->data;
  for (nd2 = 0; nd2 < m; nd2++) {
    mat_data[nd2] = fabs(cfile_data[nd2]);
  }

  GCpos1 = maximum(a__2);
  loop_ub = cfile->size[0];
  for (i = 0; i < loop_ub; i++) {
    cfile_data[i] /= GCpos1;
  }

  emxInit_real_T(&E, 1);

  /*  normalize to speed up computation */
  /*  clear A */
  /* 'gkmPWM:121' [pp, scorevec, C, r, R, E, Rd] = gkmPWM_lagrange(cfile,mat,p,negvec,num,rcorr,reg,l_svm,k_svm,RC); */
  gkmPWM_lagrange(cfile, mat, b_p, negvec, varargin_5, varargin_6, varargin_7,
                  varargin_8, varargin_9, varargin_11, dc2, a__2, &xtmp, a__3, E,
                  dc);
  dc2_data = dc2->data;

  /*  createMEME([fileheader '_' num2str(l_svm) '_' num2str(k_svm) '_' num2str(reg) '_' num2str(mnum)], pp, GCneg1, C, r, R, rcorr, E, Rd); */
  /* 'gkmPWM:124' createMEME(sprintf("%s_%d_%d_%d_%d", fileheader, int32(l_svm), int32(k_svm), int32(reg), int32(mnum)), pp, memefile, GCneg1, C, r, R, rcorr, E, Rd); */
  nd2 = (int)rt_roundd(varargin_8);
  j2 = (int)rt_roundd(varargin_9);
  b_j1 = (int)rt_roundd(varargin_7);
  u1 = (int)rt_roundd(varargin_4);
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  b_varargin_1_data = b_varargin_1->data;
  loop_ub = varargin_1->size[1];
  emxFree_real_T(&cfile);
  emxFree_real_T(&negvec);
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
  m = snprintf(NULL, 0, "%s_%d_%d_%d_%d", &charStr_data[0], nd2, j2, b_j1, u1);
  i = charStr->size[0] * charStr->size[1];
  charStr->size[0] = 1;
  charStr->size[1] = m + 1;
  emxEnsureCapacity_char_T(charStr, i);
  charStr_data = charStr->data;
  snprintf(&charStr_data[0], (size_t)(m + 1), "%s_%d_%d_%d_%d",
           &b_varargin_1_data[0], nd2, j2, b_j1, u1);
  i = charStr->size[0] * charStr->size[1];
  if (1 > m) {
    charStr->size[1] = 0;
  } else {
    charStr->size[1] = m;
  }

  emxEnsureCapacity_char_T(charStr, i);
  createMEME(charStr, b_p, varargin_3, GCneg1, a__2, xtmp, a__3, varargin_6, E,
             dc);

  /*  dlmwrite([fileheader '_' num2str(l_svm) '_' num2str(k_svm) '_' num2str(reg) '_' num2str(mnum) '_error.out'],scorevec'); */
  /* 'gkmPWM:126' fidw = fopen(sprintf("%s_%d_%d_%d_%d_error.out", fileheader, int32(l_svm), int32(k_svm), int32(reg), int32(mnum)), 'w'); */
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  b_varargin_1_data = b_varargin_1->data;
  loop_ub = varargin_1->size[1];
  emxFree_char_T(&charStr);
  emxFree_real_T(&dc);
  emxFree_real_T(&E);
  emxFree_cell_wrap_2(&b_p);
  emxFree_real_T(&a__3);
  emxFree_real_T(&a__2);
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
               b_j1, u1);
  i = b_charStr->size[0] * b_charStr->size[1];
  b_charStr->size[0] = 1;
  b_charStr->size[1] = m + 1;
  emxEnsureCapacity_char_T(b_charStr, i);
  charStr_data = b_charStr->data;
  snprintf(&charStr_data[0], (size_t)(m + 1), "%s_%d_%d_%d_%d_error.out",
           &b_varargin_1_data[0], nd2, j2, b_j1, u1);
  i = b_charStr->size[0] * b_charStr->size[1];
  if (1 > m) {
    b_charStr->size[1] = 0;
  } else {
    b_charStr->size[1] = m;
  }

  emxEnsureCapacity_char_T(b_charStr, i);
  fileid = cfopen(b_charStr, "wb");

  /* 'gkmPWM:127' for idx=1:length(scorevec) */
  i = dc2->size[1];
  emxFree_char_T(&b_charStr);
  emxFree_char_T(&c_varargin_1);
  emxFree_char_T(&b_varargin_1);
  for (m = 0; m < i; m++) {
    /* 'gkmPWM:128' fprintf(fidw, "%f\n", scorevec(idx)); */
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
