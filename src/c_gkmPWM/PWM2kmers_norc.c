/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * PWM2kmers_norc.c
 *
 * Code generation for function 'PWM2kmers_norc'
 *
 */

/* Include files */
#include "PWM2kmers_norc.h"
#include "PWM2kmers.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void q_binary_expand_op(emxArray_real_T *kweig, int i1, int i2, int i3,
                               const emxArray_cell_wrap_14 *ktree, double k,
                               int i4, int i5);

/* Function Definitions */
static void q_binary_expand_op(emxArray_real_T *kweig, int i1, int i2, int i3,
                               const emxArray_cell_wrap_14 *ktree, double k,
                               int i4, int i5)
{
  const cell_wrap_14 *ktree_data;
  emxArray_real_T *b_kweig;
  double *b_kweig_data;
  double *kweig_data;
  int i;
  int stride_0_1;
  int stride_1_1;
  int unnamed_idx_1;
  ktree_data = ktree->data;
  kweig_data = kweig->data;
  emxInit_real_T(&b_kweig, 2);
  unnamed_idx_1 = i4 - i5;
  i = b_kweig->size[0] * b_kweig->size[1];
  b_kweig->size[0] = 1;
  b_kweig->size[1] = unnamed_idx_1;
  emxEnsureCapacity_real_T(b_kweig, i);
  b_kweig_data = b_kweig->data;
  stride_0_1 = ((i3 - i2) + 1 != 1);
  stride_1_1 = (ktree_data[(int)k - 1].f1->size[0] != 1);
  for (i = 0; i < unnamed_idx_1; i++) {
    b_kweig_data[i] = kweig_data[i2 + i * stride_0_1] +
                      ktree_data[(int)k - 1].f1->data[i * stride_1_1];
  }
  unnamed_idx_1 = b_kweig->size[1];
  for (i = 0; i < unnamed_idx_1; i++) {
    kweig_data[i1 + i] = b_kweig_data[i];
  }
  emxFree_real_T(&b_kweig);
}

/*
 * function [kweig] = PWM2kmers_norc(mat,negmat,c,s,ind,indloc,x,l,k,rcnum)
 */
void PWM2kmers_norc(const emxArray_real_T *mat, const double negmat[16],
                    const emxArray_real_T *c, const emxArray_real_T *s,
                    const emxArray_real_T *ind, const emxArray_real_T *indloc,
                    const emxArray_real_T *x, double l, double k,
                    emxArray_real_T *kweig)
{
  cell_wrap_13 *p_data;
  cell_wrap_14 *ktree_data;
  emxArray_cell_wrap_13 *p;
  emxArray_cell_wrap_14 *ktree;
  emxArray_real_T *a;
  emxArray_real_T *b_ktree;
  emxArray_real_T *indvec;
  emxArray_real_T *r;
  emxArray_real_T *sPWM;
  emxArray_real_T *y;
  emxArray_uint32_T *X;
  double b_p[16];
  const double *c_data;
  const double *ind_data;
  const double *indloc_data;
  const double *mat_data;
  const double *s_data;
  const double *x_data;
  double M;
  double b_i3;
  double d;
  double d1;
  double matt;
  double n_tmp;
  double varargin_1;
  double *a_data;
  double *indvec_data;
  double *kweig_data;
  double *sPWM_data;
  double *y_data;
  int b_i;
  int b_i2;
  int b_k;
  int i;
  int i1;
  int i2;
  int i3;
  int ibcol;
  int j;
  int loop_ub;
  int m;
  int nrows;
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
  /* 'PWM2kmers_norc:3' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_13(p, i);
  p_data = p->data;
  /* 'PWM2kmers_norc:4' p = coder.nullcopy(p); */
  /* 'PWM2kmers_norc:5' p{1} = eye(4); */
  for (i = 0; i < 16; i++) {
    p_data[0].f1[i] = 0.0;
  }
  p_data[0].f1[0] = 1.0;
  p_data[0].f1[5] = 1.0;
  p_data[0].f1[10] = 1.0;
  p_data[0].f1[15] = 1.0;
  /* 'PWM2kmers_norc:6' for i = 1:l-1 */
  i = (int)(l - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers_norc:7' p{i+1} = p{i}*negmat; */
    for (i1 = 0; i1 < 4; i1++) {
      for (b_k = 0; b_k < 4; b_k++) {
        i2 = b_k << 2;
        b_p[i1 + i2] = ((p_data[b_i].f1[i1] * negmat[i2] +
                         p_data[b_i].f1[i1 + 4] * negmat[i2 + 1]) +
                        p_data[b_i].f1[i1 + 8] * negmat[i2 + 2]) +
                       p_data[b_i].f1[i1 + 12] * negmat[i2 + 3];
      }
    }
    for (i1 = 0; i1 < 16; i1++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i1] = b_p[i1];
    }
  }
  /* 'PWM2kmers_norc:9' m = length(x); */
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    m = 0;
  } else {
    ibcol = x->size[0];
    m = x->size[1];
    if (ibcol >= m) {
      m = ibcol;
    }
  }
  /* 'PWM2kmers_norc:10' M = length(mat)-l; */
  ibcol = mat->size[0];
  if (ibcol < 4) {
    ibcol = 4;
  }
  if (mat->size[0] == 0) {
    ibcol = 0;
  }
  M = (double)ibcol - l;
  /* 'PWM2kmers_norc:11' n = 4^k*length(c); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    nrows = 0;
  } else {
    ibcol = c->size[0];
    nrows = c->size[1];
    if (ibcol >= nrows) {
      nrows = ibcol;
    }
  }
  n_tmp = pow(4.0, k);
  /* number of possible k-mers */
  /* 'PWM2kmers_norc:12' mat2 = rot90(mat,2); */
  /* 'PWM2kmers_norc:13' kweig = zeros(n,1); */
  nrows = (int)(n_tmp * (double)nrows);
  i = kweig->size[0];
  kweig->size[0] = nrows;
  emxEnsureCapacity_real_T(kweig, i);
  kweig_data = kweig->data;
  for (i = 0; i < nrows; i++) {
    kweig_data[i] = 0.0;
  }
  emxInit_cell_wrap_14(&ktree);
  /* 'PWM2kmers_norc:14' ktree = cell(k,1); */
  nrows = (int)k;
  i = ktree->size[0];
  ktree->size[0] = (int)k;
  emxEnsureCapacity_cell_wrap_14(ktree, i);
  ktree_data = ktree->data;
  for (i = 0; i < nrows; i++) {
    ktree_data[i].f1->size[0] = 0;
  }
  emxInit_uint32_T(&X);
  /* 'PWM2kmers_norc:15' ktree = coder.nullcopy(ktree); */
  /* 'PWM2kmers_norc:16' [rx,cx] = size(x); */
  rx = x->size[0];
  /* 'PWM2kmers_norc:17' X = cx*ones(length(mat)-l+1,1); */
  ibcol = mat->size[0];
  if (ibcol < 4) {
    ibcol = 4;
  }
  if (mat->size[0] == 0) {
    ibcol = 0;
  }
  loop_ub = (int)(((double)ibcol - l) + 1.0);
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
  /* 'PWM2kmers_norc:18' for i = 1:cx */
  i = x->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers_norc:19' X(i) = i; */
    X_data[b_i] = (unsigned int)(b_i + 1);
  }
  /* 'PWM2kmers_norc:21' for i = 2:k */
  i = (int)(k + -1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers_norc:22' ktree{i} = zeros(4^i,1); */
    nrows = (int)pow(4.0, (double)b_i + 2.0);
    i1 = ktree_data[b_i + 1].f1->size[0];
    ktree_data[b_i + 1].f1->size[0] = nrows;
    emxEnsureCapacity_real_T(ktree_data[b_i + 1].f1, i1);
    for (i1 = 0; i1 < nrows; i1++) {
      ktree_data[b_i + 1].f1->data[i1] = 0.0;
    }
  }
  /* 'PWM2kmers_norc:24' a = 0; */
  /* 'PWM2kmers_norc:25' for i = 0:M */
  i = (int)(M + 1.0);
  emxInit_real_T(&indvec, 2);
  emxInit_real_T(&sPWM, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&r, 1);
  emxInit_real_T(&y, 1);
  emxInit_real_T(&b_ktree, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers_norc:26' if i == M-1 */
    if (b_i == M - 1.0) {
      /* 'PWM2kmers_norc:27' m = length(c); */
      if ((c->size[0] == 0) || (c->size[1] == 0)) {
        m = 0;
      } else {
        ibcol = c->size[0];
        m = c->size[1];
        if (ibcol >= m) {
          m = ibcol;
        }
      }
    }
    /* 'PWM2kmers_norc:29' for i2 = 1:m */
    for (b_i2 = 0; b_i2 < m; b_i2++) {
      /* 'PWM2kmers_norc:30' if ~(i == M-1 && i2 > rx && i2 ~= m) */
      if ((b_i != M - 1.0) || (b_i2 + 1 <= rx) || (b_i2 + 1 == m)) {
        /* 'PWM2kmers_norc:31' indvec = c(i2,:)+i; */
        loop_ub = c->size[1];
        i1 = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i1);
        indvec_data = indvec->data;
        for (i1 = 0; i1 < loop_ub; i1++) {
          indvec_data[i1] = c_data[b_i2 + c->size[0] * i1] + (double)b_i;
        }
        /* 'PWM2kmers_norc:32' loc = indloc(indvec); */
        /* 'PWM2kmers_norc:33' sPWM = mat(indvec,:).'; */
        i1 = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = 4;
        sPWM->size[1] = indvec->size[1];
        emxEnsureCapacity_real_T(sPWM, i1);
        sPWM_data = sPWM->data;
        loop_ub = indvec->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          nrows = (int)indvec_data[i1] - 1;
          sPWM_data[4 * i1] = mat_data[nrows];
          sPWM_data[4 * i1 + 1] = mat_data[nrows + mat->size[0]];
          sPWM_data[4 * i1 + 2] = mat_data[nrows + mat->size[0] * 2];
          sPWM_data[4 * i1 + 3] = mat_data[nrows + mat->size[0] * 3];
        }
        /* 'PWM2kmers_norc:34' ktree{1} = sPWM(:,1); */
        i1 = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree_data[0].f1, i1);
        ktree_data[0].f1->data[0] = sPWM_data[0];
        ktree_data[0].f1->data[1] = sPWM_data[1];
        ktree_data[0].f1->data[2] = sPWM_data[2];
        ktree_data[0].f1->data[3] = sPWM_data[3];
        /* 'PWM2kmers_norc:35' for i3 = s(i2):k */
        i1 = (int)(k + (1.0 - s_data[b_i2]));
        for (i3 = 0; i3 < i1; i3++) {
          b_i3 = s_data[b_i2] + (double)i3;
          /* 'PWM2kmers_norc:36' if loc(i3)==0 */
          if (indloc_data[(int)indvec_data[(int)b_i3 - 1] - 1] == 0.0) {
            /* 'PWM2kmers_norc:37' if loc(i3-1) == 1 */
            if (indloc_data[(int)indvec_data[(int)(b_i3 - 1.0) - 1] - 1] ==
                1.0) {
              /* 'PWM2kmers_norc:38' matt =
               * mat(1,:)*p{indvec(i3)-indvec(i3-1)+1}; */
              for (b_k = 0; b_k < 16; b_k++) {
                b_p[b_k] =
                    p_data[(int)(unsigned int)indvec_data[(int)b_i3 - 1] -
                           (int)(unsigned int)
                               indvec_data[(int)(b_i3 - 1.0) - 1]]
                        .f1[b_k];
              }
              /* 'PWM2kmers_norc:39' for i4 = 1:4 */
              varargin_1 = pow(4.0, b_i3 - 1.0);
              for (ibcol = 0; ibcol < 4; ibcol++) {
                nrows = ibcol << 2;
                matt = ((mat_data[0] * b_p[nrows] +
                         mat_data[mat->size[0]] * b_p[nrows + 1]) +
                        mat_data[mat->size[0] * 2] * b_p[nrows + 2]) +
                       mat_data[mat->size[0] * 3] * b_p[nrows + 3];
                /* 'PWM2kmers_norc:40'
                 * ktree{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4)) =
                 * ktree{i3-1}*matt(i4); */
                d = (((double)ibcol + 1.0) - 1.0) * varargin_1 + 1.0;
                d1 = varargin_1 * ((double)ibcol + 1.0);
                if (d > d1) {
                  b_k = -1;
                  i2 = 0;
                } else {
                  b_k = (int)d - 2;
                  i2 = (int)d1;
                }
                loop_ub = (i2 - b_k) - 1;
                i2 = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = loop_ub;
                emxEnsureCapacity_real_T(b_ktree, i2);
                a_data = b_ktree->data;
                for (i2 = 0; i2 < loop_ub; i2++) {
                  a_data[i2] =
                      ktree_data[(int)(b_i3 - 1.0) - 1].f1->data[i2] * matt;
                }
                loop_ub = b_ktree->size[1];
                for (i2 = 0; i2 < loop_ub; i2++) {
                  ktree_data[(int)b_i3 - 1].f1->data[(b_k + i2) + 1] =
                      a_data[i2];
                }
              }
            } else {
              /* 'PWM2kmers_norc:42' else */
              /* 'PWM2kmers_norc:43' matt = p{indvec(i3)-indvec(i3-1)+1}; */
              /* 'PWM2kmers_norc:44' ktree{i3} = repmat(ktree{i3-1}, 4,
               * 1).*repelem(matt(:), 4^(i3-2)); */
              b_k = r->size[0];
              r->size[0] = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0] << 2;
              emxEnsureCapacity_real_T(r, b_k);
              a_data = r->data;
              nrows = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
              for (loop_ub = 0; loop_ub < 4; loop_ub++) {
                ibcol = loop_ub * nrows;
                for (b_k = 0; b_k < nrows; b_k++) {
                  a_data[ibcol + b_k] =
                      ktree_data[(int)(b_i3 - 1.0) - 1].f1->data[b_k];
                }
              }
              varargin_1 = pow(4.0, b_i3 - 2.0);
              b_k = y->size[0];
              y->size[0] = (int)varargin_1 << 4;
              emxEnsureCapacity_real_T(y, b_k);
              y_data = y->data;
              nrows = -1;
              ibcol = (int)varargin_1;
              for (b_k = 0; b_k < 16; b_k++) {
                for (j = 0; j < ibcol; j++) {
                  y_data[(nrows + j) + 1] =
                      p_data[(int)(unsigned int)indvec_data[(int)b_i3 - 1] -
                             (int)(unsigned int)
                                 indvec_data[(int)(b_i3 - 1.0) - 1]]
                          .f1[b_k];
                }
                nrows += (int)varargin_1;
              }
              if (r->size[0] == y->size[0]) {
                b_k = ktree_data[(int)b_i3 - 1].f1->size[0];
                ktree_data[(int)b_i3 - 1].f1->size[0] = r->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)b_i3 - 1].f1, b_k);
                loop_ub = r->size[0];
                for (b_k = 0; b_k < loop_ub; b_k++) {
                  ktree_data[(int)b_i3 - 1].f1->data[b_k] =
                      a_data[b_k] * y_data[b_k];
                }
              } else {
                m_binary_expand_op(ktree, b_i3, r, y);
                ktree_data = ktree->data;
              }
            }
          } else {
            /* 'PWM2kmers_norc:46' else */
            /* 'PWM2kmers_norc:47' a = ktree{i3-1}.*sPWM(:,i3).'; */
            b_k = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, b_k);
            a_data = a->data;
            loop_ub = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
            for (b_k = 0; b_k < 4; b_k++) {
              for (i2 = 0; i2 < loop_ub; i2++) {
                a_data[i2 + a->size[0] * b_k] =
                    ktree_data[(int)(b_i3 - 1.0) - 1].f1->data[i2] *
                    sPWM_data[b_k + 4 * ((int)b_i3 - 1)];
              }
            }
            /* 'PWM2kmers_norc:48' ktree{i3} = a(:); */
            b_k = ktree_data[(int)b_i3 - 1].f1->size[0];
            ktree_data[(int)b_i3 - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree_data[(int)b_i3 - 1].f1, b_k);
            loop_ub = a->size[0] << 2;
            for (b_k = 0; b_k < loop_ub; b_k++) {
              ktree_data[(int)b_i3 - 1].f1->data[b_k] = a_data[b_k];
            }
          }
        }
        /* 'PWM2kmers_norc:51' if i2 <= rx */
        if (b_i2 + 1 <= rx) {
          /* 'PWM2kmers_norc:52' for j = 1:X(i+1) */
          i1 = (int)X_data[b_i];
          for (j = 0; j < i1; j++) {
            /* 'PWM2kmers_norc:53' if x(i2,j) ~= 0 */
            varargin_1 = x_data[b_i2 + x->size[0] * j];
            if (varargin_1 != 0.0) {
              /* 'PWM2kmers_norc:54'
               * kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) =
               * kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) + ktree{k};
               */
              varargin_1 = ind_data[(int)varargin_1 - 1];
              d = n_tmp * (varargin_1 - 1.0) + 1.0;
              varargin_1 *= n_tmp;
              if (d > varargin_1) {
                b_k = 0;
                i2 = 0;
                ibcol = 0;
                nrows = -1;
              } else {
                b_k = (int)d - 1;
                i2 = (int)varargin_1;
                ibcol = (int)d - 1;
                nrows = (int)varargin_1 - 1;
              }
              if (i2 - b_k == ktree_data[(int)k - 1].f1->size[0]) {
                nrows = (nrows - ibcol) + 1;
                i2 = indvec->size[0] * indvec->size[1];
                indvec->size[0] = 1;
                indvec->size[1] = nrows;
                emxEnsureCapacity_real_T(indvec, i2);
                indvec_data = indvec->data;
                for (i2 = 0; i2 < nrows; i2++) {
                  indvec_data[i2] = kweig_data[b_k + i2] +
                                    ktree_data[(int)k - 1].f1->data[i2];
                }
                loop_ub = indvec->size[1];
                for (b_k = 0; b_k < loop_ub; b_k++) {
                  kweig_data[ibcol + b_k] = indvec_data[b_k];
                }
              } else {
                q_binary_expand_op(kweig, ibcol, b_k, i2 - 1, ktree, k, nrows,
                                   ibcol - 1);
                kweig_data = kweig->data;
              }
            }
          }
        } else {
          /* 'PWM2kmers_norc:57' else */
          /* 'PWM2kmers_norc:58' kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) =
           * kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) + ktree{k}; */
          varargin_1 = n_tmp * (ind_data[b_i2] - 1.0) + 1.0;
          d = n_tmp * ind_data[b_i2];
          if (varargin_1 > d) {
            i1 = 0;
            b_k = 0;
            i2 = 0;
            ibcol = -1;
          } else {
            i1 = (int)varargin_1 - 1;
            b_k = (int)d;
            i2 = (int)varargin_1 - 1;
            ibcol = (int)d - 1;
          }
          if (b_k - i1 == ktree_data[(int)k - 1].f1->size[0]) {
            nrows = (ibcol - i2) + 1;
            b_k = indvec->size[0] * indvec->size[1];
            indvec->size[0] = 1;
            indvec->size[1] = nrows;
            emxEnsureCapacity_real_T(indvec, b_k);
            indvec_data = indvec->data;
            for (b_k = 0; b_k < nrows; b_k++) {
              indvec_data[b_k] =
                  kweig_data[i1 + b_k] + ktree_data[(int)k - 1].f1->data[b_k];
            }
            loop_ub = indvec->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              kweig_data[i2 + i1] = indvec_data[i1];
            }
          } else {
            q_binary_expand_op(kweig, i2, i1, b_k - 1, ktree, k, ibcol, i2 - 1);
            kweig_data = kweig->data;
          }
        }
      }
    }
  }
  emxFree_real_T(&b_ktree);
  emxFree_real_T(&y);
  emxFree_real_T(&r);
  emxFree_real_T(&a);
  emxFree_cell_wrap_13(&p);
  emxFree_real_T(&sPWM);
  emxFree_real_T(&indvec);
  emxFree_uint32_T(&X);
  emxFree_cell_wrap_14(&ktree);
  /* alen = length(c)-rcnum; */
  /* kweig(4^k*alen+1:end) = kweig(4^k*alen+1:end)/sqrt(2); */
}

/* End of code generation (PWM2kmers_norc.c) */