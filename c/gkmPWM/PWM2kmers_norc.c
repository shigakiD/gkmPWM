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
static void v_binary_expand_op(emxArray_real_T *kweig, int i1, int i2, int i3,
                               const emxArray_cell_wrap_14 *ktree, double k,
                               int i4, int i5);

/* Function Definitions */
static void v_binary_expand_op(emxArray_real_T *kweig, int i1, int i2, int i3,
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
  double m;
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
  int itilerow;
  int j;
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
  for (b_k = 0; b_k < 4; b_k++) {
    p_data[0].f1[b_k + (b_k << 2)] = 1.0;
  }
  /* 'PWM2kmers_norc:6' for i = 1:l-1 */
  i = (int)(l - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers_norc:7' p{i+1} = p{i}*negmat; */
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < 4; i2++) {
        varargin_1 = 0.0;
        for (itilerow = 0; itilerow < 4; itilerow++) {
          varargin_1 += p_data[b_i].f1[i1 + (itilerow << 2)] *
                        negmat[itilerow + (i2 << 2)];
        }
        b_p[i1 + (i2 << 2)] = varargin_1;
      }
    }
    for (i1 = 0; i1 < 16; i1++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i1] = b_p[i1];
    }
  }
  /* 'PWM2kmers_norc:9' m = length(x); */
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    nrows = 0;
  } else {
    ibcol = x->size[0];
    nrows = x->size[1];
    if (ibcol >= nrows) {
      nrows = ibcol;
    }
  }
  m = nrows;
  /* 'PWM2kmers_norc:10' M = length(mat)-l; */
  ibcol = mat->size[0];
  if (ibcol < 4) {
    ibcol = 4;
  }
  if (mat->size[0] == 0) {
    ibcol = 0;
  }
  M = (double)ibcol - l;
  /* 'PWM2kmers_norc:11' n = 4^k*numel(c)/k; */
  n_tmp = pow(4.0, k);
  /* number of possible k-mers */
  /* 'PWM2kmers_norc:12' mat2 = rot90(mat,2); */
  /* 'PWM2kmers_norc:13' kweig = zeros(n,1); */
  nrows = (int)(n_tmp * (double)(c->size[0] * c->size[1]) / k);
  i = kweig->size[0];
  kweig->size[0] = (int)(n_tmp * (double)(c->size[0] * c->size[1]) / k);
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
  nrows = (int)(((double)ibcol - l) + 1.0);
  i = X->size[0];
  X->size[0] = nrows;
  emxEnsureCapacity_uint32_T(X, i);
  X_data = X->data;
  for (i = 0; i < nrows; i++) {
    X_data[i] = 1U;
  }
  nrows = X->size[0];
  for (i = 0; i < nrows; i++) {
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
      /* 'PWM2kmers_norc:27' m = numel(c)/k; */
      m = (double)(c->size[0] * c->size[1]) / k;
    }
    /* 'PWM2kmers_norc:29' for i2 = 1:m */
    i1 = (int)m;
    for (b_i2 = 0; b_i2 < i1; b_i2++) {
      /* 'PWM2kmers_norc:30' if ~(i == M-1 && i2 > rx && i2 ~= m) */
      if ((b_i != M - 1.0) || (b_i2 + 1U <= (unsigned int)rx) ||
          ((double)b_i2 + 1.0 == m)) {
        /* 'PWM2kmers_norc:31' indvec = c(i2,:)+i; */
        nrows = c->size[1];
        i2 = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i2);
        indvec_data = indvec->data;
        for (i2 = 0; i2 < nrows; i2++) {
          indvec_data[i2] = c_data[b_i2 + c->size[0] * i2] + (double)b_i;
        }
        /* 'PWM2kmers_norc:32' loc = indloc(indvec); */
        /* 'PWM2kmers_norc:33' sPWM = mat(indvec,:).'; */
        i2 = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = 4;
        sPWM->size[1] = indvec->size[1];
        emxEnsureCapacity_real_T(sPWM, i2);
        sPWM_data = sPWM->data;
        nrows = indvec->size[1];
        for (i2 = 0; i2 < nrows; i2++) {
          for (itilerow = 0; itilerow < 4; itilerow++) {
            sPWM_data[itilerow + 4 * i2] =
                mat_data[((int)indvec_data[i2] + mat->size[0] * itilerow) - 1];
          }
        }
        /* 'PWM2kmers_norc:34' ktree{1} = sPWM(:,1); */
        i2 = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = 4;
        emxEnsureCapacity_real_T(ktree_data[0].f1, i2);
        for (i2 = 0; i2 < 4; i2++) {
          ktree_data[0].f1->data[i2] = sPWM_data[i2];
        }
        /* 'PWM2kmers_norc:35' for i3 = s(i2):k */
        i2 = (int)(k + (1.0 - s_data[b_i2]));
        for (i3 = 0; i3 < i2; i3++) {
          b_i3 = s_data[b_i2] + (double)i3;
          /* 'PWM2kmers_norc:36' if loc(i3)==0 */
          if (indloc_data[(int)indvec_data[(int)b_i3 - 1] - 1] == 0.0) {
            /* 'PWM2kmers_norc:37' if loc(i3-1) == 1 */
            if (indloc_data[(int)indvec_data[(int)(b_i3 - 1.0) - 1] - 1] ==
                1.0) {
              /* 'PWM2kmers_norc:38' matt =
               * mat(1,:)*p{indvec(i3)-indvec(i3-1)+1}; */
              for (itilerow = 0; itilerow < 16; itilerow++) {
                b_p[itilerow] =
                    p_data[(int)(unsigned int)indvec_data[(int)b_i3 - 1] -
                           (int)(unsigned int)
                               indvec_data[(int)(b_i3 - 1.0) - 1]]
                        .f1[itilerow];
              }
              /* 'PWM2kmers_norc:39' for i4 = 1:4 */
              varargin_1 = pow(4.0, b_i3 - 1.0);
              for (ibcol = 0; ibcol < 4; ibcol++) {
                matt = 0.0;
                for (itilerow = 0; itilerow < 4; itilerow++) {
                  matt += mat_data[mat->size[0] * itilerow] *
                          b_p[itilerow + (ibcol << 2)];
                }
                /* 'PWM2kmers_norc:40'
                 * ktree{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4)) =
                 * ktree{i3-1}*matt(i4); */
                d = (((double)ibcol + 1.0) - 1.0) * varargin_1 + 1.0;
                d1 = varargin_1 * ((double)ibcol + 1.0);
                if (d > d1) {
                  itilerow = -1;
                  b_k = 0;
                } else {
                  itilerow = (int)d - 2;
                  b_k = (int)d1;
                }
                nrows = (b_k - itilerow) - 1;
                b_k = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = nrows;
                emxEnsureCapacity_real_T(b_ktree, b_k);
                a_data = b_ktree->data;
                for (b_k = 0; b_k < nrows; b_k++) {
                  a_data[b_k] =
                      ktree_data[(int)(b_i3 - 1.0) - 1].f1->data[b_k] * matt;
                }
                nrows = b_ktree->size[1];
                for (b_k = 0; b_k < nrows; b_k++) {
                  ktree_data[(int)b_i3 - 1].f1->data[(itilerow + b_k) + 1] =
                      a_data[b_k];
                }
              }
            } else {
              /* 'PWM2kmers_norc:42' else */
              /* 'PWM2kmers_norc:43' matt = p{indvec(i3)-indvec(i3-1)+1}; */
              /* 'PWM2kmers_norc:44' ktree{i3} = repmat(ktree{i3-1}, 4,
               * 1).*repelem(matt(:), 4^(i3-2)); */
              itilerow = r->size[0];
              r->size[0] = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0] << 2;
              emxEnsureCapacity_real_T(r, itilerow);
              a_data = r->data;
              nrows = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
              for (itilerow = 0; itilerow < 4; itilerow++) {
                ibcol = itilerow * nrows;
                for (b_k = 0; b_k < nrows; b_k++) {
                  a_data[ibcol + b_k] =
                      ktree_data[(int)(b_i3 - 1.0) - 1].f1->data[b_k];
                }
              }
              varargin_1 = pow(4.0, b_i3 - 2.0);
              itilerow = y->size[0];
              y->size[0] = (int)varargin_1 << 4;
              emxEnsureCapacity_real_T(y, itilerow);
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
                itilerow = ktree_data[(int)b_i3 - 1].f1->size[0];
                ktree_data[(int)b_i3 - 1].f1->size[0] = r->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)b_i3 - 1].f1,
                                         itilerow);
                nrows = r->size[0];
                for (itilerow = 0; itilerow < nrows; itilerow++) {
                  ktree_data[(int)b_i3 - 1].f1->data[itilerow] =
                      a_data[itilerow] * y_data[itilerow];
                }
              } else {
                r_binary_expand_op(ktree, b_i3, r, y);
                ktree_data = ktree->data;
              }
            }
          } else {
            /* 'PWM2kmers_norc:46' else */
            /* 'PWM2kmers_norc:47' a = ktree{i3-1}.*sPWM(:,i3).'; */
            itilerow = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
            a->size[1] = 4;
            emxEnsureCapacity_real_T(a, itilerow);
            a_data = a->data;
            nrows = ktree_data[(int)(b_i3 - 1.0) - 1].f1->size[0];
            for (itilerow = 0; itilerow < 4; itilerow++) {
              for (b_k = 0; b_k < nrows; b_k++) {
                a_data[b_k + a->size[0] * itilerow] =
                    ktree_data[(int)(b_i3 - 1.0) - 1].f1->data[b_k] *
                    sPWM_data[itilerow + 4 * ((int)b_i3 - 1)];
              }
            }
            /* 'PWM2kmers_norc:48' ktree{i3} = a(:); */
            itilerow = ktree_data[(int)b_i3 - 1].f1->size[0];
            ktree_data[(int)b_i3 - 1].f1->size[0] = a->size[0] << 2;
            emxEnsureCapacity_real_T(ktree_data[(int)b_i3 - 1].f1, itilerow);
            nrows = a->size[0] << 2;
            for (itilerow = 0; itilerow < nrows; itilerow++) {
              ktree_data[(int)b_i3 - 1].f1->data[itilerow] = a_data[itilerow];
            }
          }
        }
        /* 'PWM2kmers_norc:51' if i2 <= rx */
        if (b_i2 + 1U <= (unsigned int)rx) {
          /* 'PWM2kmers_norc:52' for j = 1:X(i+1) */
          i2 = (int)X_data[b_i];
          for (j = 0; j < i2; j++) {
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
                itilerow = 0;
                b_k = 0;
                ibcol = 0;
                nrows = -1;
              } else {
                itilerow = (int)d - 1;
                b_k = (int)varargin_1;
                ibcol = (int)d - 1;
                nrows = (int)varargin_1 - 1;
              }
              if (b_k - itilerow == ktree_data[(int)k - 1].f1->size[0]) {
                nrows = (nrows - ibcol) + 1;
                b_k = indvec->size[0] * indvec->size[1];
                indvec->size[0] = 1;
                indvec->size[1] = nrows;
                emxEnsureCapacity_real_T(indvec, b_k);
                indvec_data = indvec->data;
                for (b_k = 0; b_k < nrows; b_k++) {
                  indvec_data[b_k] = kweig_data[itilerow + b_k] +
                                     ktree_data[(int)k - 1].f1->data[b_k];
                }
                nrows = indvec->size[1];
                for (itilerow = 0; itilerow < nrows; itilerow++) {
                  kweig_data[ibcol + itilerow] = indvec_data[itilerow];
                }
              } else {
                v_binary_expand_op(kweig, ibcol, itilerow, b_k - 1, ktree, k,
                                   nrows, ibcol - 1);
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
            i2 = 0;
            itilerow = 0;
            b_k = 0;
            ibcol = -1;
          } else {
            i2 = (int)varargin_1 - 1;
            itilerow = (int)d;
            b_k = (int)varargin_1 - 1;
            ibcol = (int)d - 1;
          }
          if (itilerow - i2 == ktree_data[(int)k - 1].f1->size[0]) {
            nrows = (ibcol - b_k) + 1;
            itilerow = indvec->size[0] * indvec->size[1];
            indvec->size[0] = 1;
            indvec->size[1] = nrows;
            emxEnsureCapacity_real_T(indvec, itilerow);
            indvec_data = indvec->data;
            for (itilerow = 0; itilerow < nrows; itilerow++) {
              indvec_data[itilerow] = kweig_data[i2 + itilerow] +
                                      ktree_data[(int)k - 1].f1->data[itilerow];
            }
            nrows = indvec->size[1];
            for (i2 = 0; i2 < nrows; i2++) {
              kweig_data[b_k + i2] = indvec_data[i2];
            }
          } else {
            v_binary_expand_op(kweig, b_k, i2, itilerow - 1, ktree, k, ibcol,
                               b_k - 1);
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
