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
#include "eye.h"
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include "repelem.h"
#include "repmat.h"
#include "rot90.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void b_plus(emxArray_real_T *kweig, const emxArray_real_T *kweig2);

/* Function Definitions */
static void b_plus(emxArray_real_T *kweig, const emxArray_real_T *kweig2)
{
  emxArray_real_T *b_kweig;
  const double *kweig2_data;
  double *b_kweig_data;
  double *kweig_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  kweig2_data = kweig2->data;
  kweig_data = kweig->data;
  emxInit_real_T(&b_kweig, 1);
  i = b_kweig->size[0];
  if (kweig2->size[0] == 1) {
    b_kweig->size[0] = kweig->size[0];
  } else {
    b_kweig->size[0] = kweig2->size[0];
  }
  emxEnsureCapacity_real_T(b_kweig, i);
  b_kweig_data = b_kweig->data;
  stride_0_0 = (kweig->size[0] != 1);
  stride_1_0 = (kweig2->size[0] != 1);
  if (kweig2->size[0] == 1) {
    loop_ub = kweig->size[0];
  } else {
    loop_ub = kweig2->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_kweig_data[i] = kweig_data[i * stride_0_0] + kweig2_data[i * stride_1_0];
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
 * function [kweig,C] = PWM2kmers(mat,negmat,c,s,ind,indloc,x,l,k,rcnum)
 */
void PWM2kmers(const emxArray_real_T *mat, const double negmat[16],
               const emxArray_real_T *c, const emxArray_real_T *s,
               const emxArray_real_T *ind, const emxArray_real_T *indloc,
               const emxArray_real_T *x, double l, double k, double rcnum,
               emxArray_real_T *kweig)
{
  cell_wrap_12 *p_data;
  cell_wrap_3 *ktree2_data;
  cell_wrap_3 *ktree_data;
  emxArray_cell_wrap_12 *p;
  emxArray_cell_wrap_3 *ktree;
  emxArray_cell_wrap_3 *ktree2;
  emxArray_real_T *a;
  emxArray_real_T *b_ktree;
  emxArray_real_T *b_r;
  emxArray_real_T *indvec;
  emxArray_real_T *kweig2;
  emxArray_real_T *mat2;
  emxArray_real_T *sPWM;
  emxArray_real_T *sPWM2;
  emxArray_real_T *sPWM_tmp;
  emxArray_uint32_T *X;
  double matt[16];
  double c_matt[4];
  const double *c_data;
  const double *ind_data;
  const double *indloc_data;
  const double *mat_data;
  const double *s_data;
  const double *x_data;
  double M;
  double alen;
  double b_matt;
  double d;
  double d1;
  double d2;
  double n_tmp;
  double *a_data;
  double *indvec_data;
  double *kweig2_data;
  double *kweig_data;
  double *mat2_data;
  double *r1;
  double *sPWM2_data;
  double *sPWM_data;
  int b_i;
  int b_i2;
  int b_i3;
  int i;
  int i1;
  int i2;
  int i3;
  int loop_ub;
  int m;
  int rx;
  int u0;
  int u1;
  unsigned int *X_data;
  x_data = x->data;
  indloc_data = indloc->data;
  ind_data = ind->data;
  s_data = s->data;
  c_data = c->data;
  mat_data = mat->data;
  emxInit_cell_wrap_12(&p);
  /* makes gkm-pwm faster */
  /* 'PWM2kmers:3' p = cell(l,1); */
  i = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_12(p, i);
  p_data = p->data;
  /* 'PWM2kmers:4' p = coder.nullcopy(p); */
  /* 'PWM2kmers:5' p{1} = eye(4); */
  eye(p_data[0].f1);
  /* 'PWM2kmers:6' for i = 1:l-1 */
  i = (int)(l - 1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:7' p{i+1} = p{i}*negmat; */
    for (i1 = 0; i1 < 4; i1++) {
      for (i2 = 0; i2 < 4; i2++) {
        i3 = i2 << 2;
        matt[i1 + i3] = ((p_data[b_i].f1[i1] * negmat[i3] +
                          p_data[b_i].f1[i1 + 4] * negmat[i3 + 1]) +
                         p_data[b_i].f1[i1 + 8] * negmat[i3 + 2]) +
                        p_data[b_i].f1[i1 + 12] * negmat[i3 + 3];
      }
    }
    for (i1 = 0; i1 < 16; i1++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i1] = matt[i1];
    }
  }
  /* 'PWM2kmers:9' m = length(x); */
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    m = 0;
  } else {
    u0 = x->size[0];
    m = x->size[1];
    if (u0 >= m) {
      m = u0;
    }
  }
  /* 'PWM2kmers:10' M = length(mat)-l; */
  if ((mat->size[0] == 0) || (mat->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = mat->size[0];
    u1 = mat->size[1];
    if (u0 >= u1) {
      u1 = u0;
    }
  }
  M = (double)u1 - l;
  /* 'PWM2kmers:11' n = 4^k*length(c); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = c->size[0];
    u1 = c->size[1];
    if (u0 >= u1) {
      u1 = u0;
    }
  }
  emxInit_real_T(&mat2, 2);
  n_tmp = pow(4.0, k);
  alen = n_tmp * (double)u1;
  /* number of possible k-mers */
  /* 'PWM2kmers:12' mat2 = rot90(mat,2); */
  rot90(mat, mat2);
  mat2_data = mat2->data;
  /* 'PWM2kmers:13' kweig = zeros(n,1); */
  u1 = (int)alen;
  i = kweig->size[0];
  kweig->size[0] = (int)alen;
  emxEnsureCapacity_real_T(kweig, i);
  kweig_data = kweig->data;
  for (i = 0; i < u1; i++) {
    kweig_data[i] = 0.0;
  }
  emxInit_real_T(&kweig2, 1);
  /* 'PWM2kmers:14' kweig2= zeros(n,1); */
  i = kweig2->size[0];
  kweig2->size[0] = (int)alen;
  emxEnsureCapacity_real_T(kweig2, i);
  kweig2_data = kweig2->data;
  for (i = 0; i < u1; i++) {
    kweig2_data[i] = 0.0;
  }
  emxInit_cell_wrap_3(&ktree);
  /* 'PWM2kmers:15' ktree = cell(k,1); */
  u1 = (int)k;
  i = ktree->size[0];
  ktree->size[0] = (int)k;
  emxEnsureCapacity_cell_wrap_32(ktree, i);
  ktree_data = ktree->data;
  for (i = 0; i < u1; i++) {
    ktree_data[i].f1->size[0] = 0;
  }
  emxInit_cell_wrap_3(&ktree2);
  /* 'PWM2kmers:16' ktree = coder.nullcopy(ktree); */
  /* 'PWM2kmers:17' ktree2 = cell(k,1); */
  i = ktree2->size[0];
  ktree2->size[0] = (int)k;
  emxEnsureCapacity_cell_wrap_32(ktree2, i);
  ktree2_data = ktree2->data;
  for (i = 0; i < u1; i++) {
    ktree2_data[i].f1->size[0] = 0;
  }
  /* 'PWM2kmers:18' ktree2 = coder.nullcopy(ktree2); */
  /* 'PWM2kmers:19' [rx,cx] = size(x); */
  rx = x->size[0];
  /* 'PWM2kmers:20' X = cx*ones(length(mat)-l+1,1); */
  if ((mat->size[0] == 0) || (mat->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = mat->size[0];
    u1 = mat->size[1];
    if (u0 >= u1) {
      u1 = u0;
    }
  }
  emxInit_uint32_T(&X, 1);
  loop_ub = (int)(((double)u1 - l) + 1.0);
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
  /* 'PWM2kmers:21' for i = 1:cx */
  i = x->size[1];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:22' X(i) = i; */
    X_data[b_i] = (unsigned int)(b_i + 1);
  }
  /* 'PWM2kmers:24' for i = 2:k */
  i = (int)(k + -1.0);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:25' ktree{i} = zeros(4^i,1); */
    alen = pow(4.0, (double)b_i + 2.0);
    u1 = (int)pow(4.0, (double)b_i + 2.0);
    i1 = ktree_data[b_i + 1].f1->size[0];
    ktree_data[b_i + 1].f1->size[0] = (int)alen;
    emxEnsureCapacity_real_T(ktree_data[b_i + 1].f1, i1);
    for (i1 = 0; i1 < u1; i1++) {
      ktree_data[b_i + 1].f1->data[i1] = 0.0;
    }
    /* 'PWM2kmers:26' ktree2{i} = zeros(4^i,1); */
    i1 = ktree2_data[b_i + 1].f1->size[0];
    ktree2_data[b_i + 1].f1->size[0] = (int)alen;
    emxEnsureCapacity_real_T(ktree2_data[b_i + 1].f1, i1);
    for (i1 = 0; i1 < u1; i1++) {
      ktree2_data[b_i + 1].f1->data[i1] = 0.0;
    }
  }
  /* 'PWM2kmers:28' a = 0; */
  /* 'PWM2kmers:29' for i = 0:M */
  i = (int)(M + 1.0);
  emxInit_real_T(&indvec, 2);
  emxInit_real_T(&sPWM, 2);
  emxInit_real_T(&sPWM2, 2);
  emxInit_real_T(&a, 2);
  emxInit_real_T(&sPWM_tmp, 1);
  emxInit_real_T(&b_r, 1);
  emxInit_real_T(&b_ktree, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'PWM2kmers:30' if i == M-1 */
    if (b_i == M - 1.0) {
      /* 'PWM2kmers:31' m = length(c); */
      if ((c->size[0] == 0) || (c->size[1] == 0)) {
        m = 0;
      } else {
        u0 = c->size[0];
        m = c->size[1];
        if (u0 >= m) {
          m = u0;
        }
      }
    }
    /* 'PWM2kmers:33' for i2 = 1:m */
    for (b_i2 = 0; b_i2 < m; b_i2++) {
      /* 'PWM2kmers:34' if ~(i == M-1 && i2 > rx && i2 ~= m) */
      if ((b_i != M - 1.0) || (b_i2 + 1 <= rx) || (b_i2 + 1 == m)) {
        /* 'PWM2kmers:35' indvec = c(i2,:)+i; */
        loop_ub = c->size[1];
        i1 = indvec->size[0] * indvec->size[1];
        indvec->size[0] = 1;
        indvec->size[1] = c->size[1];
        emxEnsureCapacity_real_T(indvec, i1);
        indvec_data = indvec->data;
        for (i1 = 0; i1 < loop_ub; i1++) {
          indvec_data[i1] = c_data[b_i2 + c->size[0] * i1] + (double)b_i;
        }
        /* 'PWM2kmers:36' loc = indloc(indvec); */
        /* 'PWM2kmers:37' sPWM = mat(indvec,:).'; */
        i1 = sPWM_tmp->size[0];
        sPWM_tmp->size[0] = indvec->size[1];
        emxEnsureCapacity_real_T(sPWM_tmp, i1);
        a_data = sPWM_tmp->data;
        loop_ub = indvec->size[1];
        for (i1 = 0; i1 < loop_ub; i1++) {
          a_data[i1] = indvec_data[i1];
        }
        loop_ub = mat->size[1];
        i1 = sPWM->size[0] * sPWM->size[1];
        sPWM->size[0] = mat->size[1];
        sPWM->size[1] = sPWM_tmp->size[0];
        emxEnsureCapacity_real_T(sPWM, i1);
        sPWM_data = sPWM->data;
        u1 = sPWM_tmp->size[0];
        for (i1 = 0; i1 < u1; i1++) {
          for (i2 = 0; i2 < loop_ub; i2++) {
            sPWM_data[i2 + sPWM->size[0] * i1] =
                mat_data[((int)a_data[i1] + mat->size[0] * i2) - 1];
          }
        }
        /* 'PWM2kmers:38' sPWM2 = mat2(indvec,:).'; */
        loop_ub = mat2->size[1];
        i1 = sPWM2->size[0] * sPWM2->size[1];
        sPWM2->size[0] = mat2->size[1];
        sPWM2->size[1] = sPWM_tmp->size[0];
        emxEnsureCapacity_real_T(sPWM2, i1);
        sPWM2_data = sPWM2->data;
        u1 = sPWM_tmp->size[0];
        for (i1 = 0; i1 < u1; i1++) {
          for (i2 = 0; i2 < loop_ub; i2++) {
            sPWM2_data[i2 + sPWM2->size[0] * i1] =
                mat2_data[((int)a_data[i1] + mat2->size[0] * i2) - 1];
          }
        }
        /* 'PWM2kmers:39' ktree{1} = sPWM(:,1); */
        loop_ub = sPWM->size[0];
        i1 = ktree_data[0].f1->size[0];
        ktree_data[0].f1->size[0] = sPWM->size[0];
        emxEnsureCapacity_real_T(ktree_data[0].f1, i1);
        for (i1 = 0; i1 < loop_ub; i1++) {
          ktree_data[0].f1->data[i1] = sPWM_data[i1];
        }
        /* 'PWM2kmers:40' ktree2{1} = sPWM2(:,1); */
        loop_ub = sPWM2->size[0];
        i1 = ktree2_data[0].f1->size[0];
        ktree2_data[0].f1->size[0] = sPWM2->size[0];
        emxEnsureCapacity_real_T(ktree2_data[0].f1, i1);
        for (i1 = 0; i1 < loop_ub; i1++) {
          ktree2_data[0].f1->data[i1] = sPWM2_data[i1];
        }
        /* 'PWM2kmers:41' for i3 = s(i2):k */
        i1 = (int)(k + (1.0 - s_data[b_i2]));
        for (b_i3 = 0; b_i3 < i1; b_i3++) {
          alen = s_data[b_i2] + (double)b_i3;
          /* 'PWM2kmers:42' if loc(i3)==0 */
          if (indloc_data[(int)indvec_data[(int)alen - 1] - 1] == 0.0) {
            /* 'PWM2kmers:43' if loc(i3-1) == 1 */
            if (indloc_data[(int)indvec_data[(int)(alen - 1.0) - 1] - 1] ==
                1.0) {
              /* 'PWM2kmers:44' matt = mat(1,:)*p{indvec(i3)-indvec(i3-1)+1}; */
              for (i2 = 0; i2 < 16; i2++) {
                matt[i2] =
                    p_data[(int)(unsigned int)indvec_data[(int)alen - 1] -
                           (int)(unsigned int)
                               indvec_data[(int)(alen - 1.0) - 1]]
                        .f1[i2];
              }
              /* 'PWM2kmers:45' for i4 = 1:4 */
              d = pow(4.0, alen - 1.0);
              for (u0 = 0; u0 < 4; u0++) {
                u1 = u0 << 2;
                b_matt = ((mat_data[0] * matt[u1] +
                           mat_data[mat->size[0]] * matt[u1 + 1]) +
                          mat_data[mat->size[0] * 2] * matt[u1 + 2]) +
                         mat_data[mat->size[0] * 3] * matt[u1 + 3];
                c_matt[u0] = b_matt;
                /* 'PWM2kmers:46' ktree{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4)) =
                 * ktree{i3-1}*matt(i4); */
                d1 = (((double)u0 + 1.0) - 1.0) * d + 1.0;
                d2 = d * ((double)u0 + 1.0);
                if (d1 > d2) {
                  i2 = -1;
                  i3 = 0;
                } else {
                  i2 = (int)d1 - 2;
                  i3 = (int)d2;
                }
                u1 = (i3 - i2) - 1;
                i3 = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = u1;
                emxEnsureCapacity_real_T(b_ktree, i3);
                a_data = b_ktree->data;
                for (i3 = 0; i3 < u1; i3++) {
                  a_data[i3] =
                      ktree_data[(int)(alen - 1.0) - 1].f1->data[i3] * b_matt;
                }
                loop_ub = b_ktree->size[1];
                for (i3 = 0; i3 < loop_ub; i3++) {
                  ktree_data[(int)alen - 1].f1->data[(i2 + i3) + 1] =
                      a_data[i3];
                }
                /* 'PWM2kmers:47' ktree2{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4))
                 * = ktree2{i3-1}*matt(i4); */
                d1 = (((double)u0 + 1.0) - 1.0) * pow(4.0, alen - 1.0) + 1.0;
                d2 = pow(4.0, alen - 1.0) * ((double)u0 + 1.0);
                if (d1 > d2) {
                  i2 = -1;
                  i3 = 0;
                } else {
                  i2 = (int)d1 - 2;
                  i3 = (int)d2;
                }
                u1 = (i3 - i2) - 1;
                i3 = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = u1;
                emxEnsureCapacity_real_T(b_ktree, i3);
                a_data = b_ktree->data;
                for (i3 = 0; i3 < u1; i3++) {
                  a_data[i3] = ktree2_data[(int)(alen - 1.0) - 1].f1->data[i3] *
                               c_matt[u0];
                }
                loop_ub = b_ktree->size[1];
                for (i3 = 0; i3 < loop_ub; i3++) {
                  ktree2_data[(int)alen - 1].f1->data[(i2 + i3) + 1] =
                      a_data[i3];
                }
              }
            } else {
              /* 'PWM2kmers:49' else */
              /* 'PWM2kmers:50' matt = p{indvec(i3)-indvec(i3-1)+1}; */
              for (i2 = 0; i2 < 16; i2++) {
                matt[i2] =
                    p_data[(int)(unsigned int)indvec_data[(int)alen - 1] -
                           (int)(unsigned int)
                               indvec_data[(int)(alen - 1.0) - 1]]
                        .f1[i2];
              }
              /* 'PWM2kmers:51' ktree{i3} = repmat(ktree{i3-1}, 4,
               * 1).*repelem(matt(:), 4^(i3-2)); */
              repelem(matt, pow(4.0, alen - 2.0), sPWM_tmp);
              a_data = sPWM_tmp->data;
              c_repmat(ktree_data[(int)(alen - 1.0) - 1].f1, b_r);
              r1 = b_r->data;
              if (b_r->size[0] == sPWM_tmp->size[0]) {
                i2 = ktree_data[(int)alen - 1].f1->size[0];
                ktree_data[(int)alen - 1].f1->size[0] = b_r->size[0];
                emxEnsureCapacity_real_T(ktree_data[(int)alen - 1].f1, i2);
                loop_ub = b_r->size[0];
                for (i2 = 0; i2 < loop_ub; i2++) {
                  ktree_data[(int)alen - 1].f1->data[i2] = r1[i2] * a_data[i2];
                }
              } else {
                l_binary_expand_op(ktree, alen, b_r, sPWM_tmp);
                ktree_data = ktree->data;
              }
              /* 'PWM2kmers:52' ktree2{i3} = repmat(ktree2{i3-1}, 4,
               * 1).*repelem(matt(:), 4^(i3-2)); */
              c_repmat(ktree2_data[(int)(alen - 1.0) - 1].f1, b_r);
              r1 = b_r->data;
              if (b_r->size[0] == sPWM_tmp->size[0]) {
                i2 = ktree2_data[(int)alen - 1].f1->size[0];
                ktree2_data[(int)alen - 1].f1->size[0] = b_r->size[0];
                emxEnsureCapacity_real_T(ktree2_data[(int)alen - 1].f1, i2);
                loop_ub = b_r->size[0];
                for (i2 = 0; i2 < loop_ub; i2++) {
                  ktree2_data[(int)alen - 1].f1->data[i2] = r1[i2] * a_data[i2];
                }
              } else {
                l_binary_expand_op(ktree2, alen, b_r, sPWM_tmp);
                ktree2_data = ktree2->data;
              }
              /* for i4 = 1:4 */
              /*     for i5 = 1:4 */
              /*         ktree{i3}(((i5-1)*4^(i3-2)+(i4-1)*4^(i3-1)+1):(4^(i3-2)*i5+(i4-1)*4^(i3-1)))=ktree{i3-1}(((i5-1)*4^(i3-2)+1):(4^(i3-2)*i5))*matt(i5,i4);
               */
              /*         ktree2{i3}(((i5-1)*4^(i3-2)+(i4-1)*4^(i3-1)+1):(4^(i3-2)*i5+(i4-1)*4^(i3-1)))=ktree2{i3-1}(((i5-1)*4^(i3-2)+1):(4^(i3-2)*i5))*matt(i5,i4);
               */
              /*     end */
              /* end */
            }
          } else {
            /* 'PWM2kmers:60' else */
            /* 'PWM2kmers:61' a = ktree{i3-1}.*sPWM(:,i3).'; */
            loop_ub = sPWM->size[0];
            i2 = a->size[0] * a->size[1];
            a->size[0] = ktree_data[(int)(alen - 1.0) - 1].f1->size[0];
            a->size[1] = sPWM->size[0];
            emxEnsureCapacity_real_T(a, i2);
            a_data = a->data;
            for (i2 = 0; i2 < loop_ub; i2++) {
              u1 = ktree_data[(int)(alen - 1.0) - 1].f1->size[0];
              for (i3 = 0; i3 < u1; i3++) {
                a_data[i3 + a->size[0] * i2] =
                    ktree_data[(int)(alen - 1.0) - 1].f1->data[i3] *
                    sPWM_data[i2 + sPWM->size[0] * ((int)alen - 1)];
              }
            }
            /* 'PWM2kmers:62' ktree{i3} = a(:); */
            i2 = ktree_data[(int)alen - 1].f1->size[0];
            ktree_data[(int)alen - 1].f1->size[0] = a->size[0] * a->size[1];
            emxEnsureCapacity_real_T(ktree_data[(int)alen - 1].f1, i2);
            loop_ub = a->size[0] * a->size[1];
            for (i2 = 0; i2 < loop_ub; i2++) {
              ktree_data[(int)alen - 1].f1->data[i2] = a_data[i2];
            }
            /* 'PWM2kmers:63' a = ktree2{i3-1}.*sPWM2(:,i3).'; */
            loop_ub = sPWM2->size[0];
            i2 = a->size[0] * a->size[1];
            a->size[0] = ktree2_data[(int)(alen - 1.0) - 1].f1->size[0];
            a->size[1] = sPWM2->size[0];
            emxEnsureCapacity_real_T(a, i2);
            a_data = a->data;
            for (i2 = 0; i2 < loop_ub; i2++) {
              u1 = ktree2_data[(int)(alen - 1.0) - 1].f1->size[0];
              for (i3 = 0; i3 < u1; i3++) {
                a_data[i3 + a->size[0] * i2] =
                    ktree2_data[(int)(alen - 1.0) - 1].f1->data[i3] *
                    sPWM2_data[i2 + sPWM2->size[0] * ((int)alen - 1)];
              }
            }
            /* 'PWM2kmers:64' ktree2{i3} = a(:); */
            i2 = ktree2_data[(int)alen - 1].f1->size[0];
            ktree2_data[(int)alen - 1].f1->size[0] = a->size[0] * a->size[1];
            emxEnsureCapacity_real_T(ktree2_data[(int)alen - 1].f1, i2);
            loop_ub = a->size[0] * a->size[1];
            for (i2 = 0; i2 < loop_ub; i2++) {
              ktree2_data[(int)alen - 1].f1->data[i2] = a_data[i2];
            }
            /* for i4 = 1:4 */
            /*     ktree{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4)) =
             * ktree{i3-1}*sPWM(i4,i3); */
            /*     ktree2{i3}(((i4-1)*4^(i3-1)+1):(4^(i3-1)*i4)) =
             * ktree2{i3-1}*sPWM2(i4,i3); */
            /* end */
          }
        }
        /* 'PWM2kmers:71' if i2 <= rx */
        if (b_i2 + 1 <= rx) {
          /* 'PWM2kmers:72' for j = 1:X(i+1) */
          i1 = (int)X_data[b_i];
          for (b_i3 = 0; b_i3 < i1; b_i3++) {
            /* 'PWM2kmers:73' if x(i2,j) ~= 0 */
            d = x_data[b_i2 + x->size[0] * b_i3];
            if (d != 0.0) {
              /* 'PWM2kmers:74'
               * kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) =
               * kweig((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) + ktree{k};
               */
              d = ind_data[(int)d - 1];
              d1 = n_tmp * (d - 1.0) + 1.0;
              d *= n_tmp;
              if (d1 > d) {
                i2 = 0;
                i3 = 0;
                u0 = 0;
                u1 = -1;
              } else {
                i2 = (int)d1 - 1;
                i3 = (int)d;
                u0 = (int)d1 - 1;
                u1 = (int)d - 1;
              }
              if (i3 - i2 == ktree_data[(int)k - 1].f1->size[0]) {
                u1 = (u1 - u0) + 1;
                i3 = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = u1;
                emxEnsureCapacity_real_T(b_ktree, i3);
                a_data = b_ktree->data;
                for (i3 = 0; i3 < u1; i3++) {
                  a_data[i3] =
                      kweig_data[i2 + i3] + ktree_data[(int)k - 1].f1->data[i3];
                }
                loop_ub = b_ktree->size[1];
                for (i2 = 0; i2 < loop_ub; i2++) {
                  kweig_data[u0 + i2] = a_data[i2];
                }
              } else {
                m_binary_expand_op(kweig, u0, i2, i3 - 1, ktree, k, u1, u0 - 1);
                kweig_data = kweig->data;
              }
              /* 'PWM2kmers:75'
               * kweig2((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) =
               * kweig2((4^k*(ind(x(i2,j))-1)+1):4^k*(ind(x(i2,j)))) +
               * ktree2{k}; */
              d = ind_data[(int)x_data[b_i2 + x->size[0] * b_i3] - 1];
              d1 = n_tmp * (d - 1.0) + 1.0;
              d *= n_tmp;
              if (d1 > d) {
                i2 = 0;
                i3 = 0;
                u0 = 0;
                u1 = -1;
              } else {
                i2 = (int)d1 - 1;
                i3 = (int)d;
                u0 = (int)d1 - 1;
                u1 = (int)d - 1;
              }
              if (i3 - i2 == ktree2_data[(int)k - 1].f1->size[0]) {
                u1 = (u1 - u0) + 1;
                i3 = b_ktree->size[0] * b_ktree->size[1];
                b_ktree->size[0] = 1;
                b_ktree->size[1] = u1;
                emxEnsureCapacity_real_T(b_ktree, i3);
                a_data = b_ktree->data;
                for (i3 = 0; i3 < u1; i3++) {
                  a_data[i3] = kweig2_data[i2 + i3] +
                               ktree2_data[(int)k - 1].f1->data[i3];
                }
                loop_ub = b_ktree->size[1];
                for (i2 = 0; i2 < loop_ub; i2++) {
                  kweig2_data[u0 + i2] = a_data[i2];
                }
              } else {
                m_binary_expand_op(kweig2, u0, i2, i3 - 1, ktree2, k, u1,
                                   u0 - 1);
                kweig2_data = kweig2->data;
              }
            }
          }
        } else {
          /* 'PWM2kmers:78' else */
          /* 'PWM2kmers:79' kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) =
           * kweig((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) + ktree{k}; */
          d = n_tmp * (ind_data[b_i2] - 1.0) + 1.0;
          d1 = n_tmp * ind_data[b_i2];
          if (d > d1) {
            i1 = 0;
            i2 = 0;
            i3 = 0;
            u0 = -1;
          } else {
            i1 = (int)d - 1;
            i2 = (int)d1;
            i3 = (int)d - 1;
            u0 = (int)d1 - 1;
          }
          if (i2 - i1 == ktree_data[(int)k - 1].f1->size[0]) {
            u1 = (u0 - i3) + 1;
            i2 = b_ktree->size[0] * b_ktree->size[1];
            b_ktree->size[0] = 1;
            b_ktree->size[1] = u1;
            emxEnsureCapacity_real_T(b_ktree, i2);
            a_data = b_ktree->data;
            for (i2 = 0; i2 < u1; i2++) {
              a_data[i2] =
                  kweig_data[i1 + i2] + ktree_data[(int)k - 1].f1->data[i2];
            }
            loop_ub = b_ktree->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              kweig_data[i3 + i1] = a_data[i1];
            }
          } else {
            m_binary_expand_op(kweig, i3, i1, i2 - 1, ktree, k, u0, i3 - 1);
            kweig_data = kweig->data;
          }
          /* 'PWM2kmers:80' kweig2((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) =
           * kweig2((4^k*(ind(i2)-1)+1):4^k*(ind(i2))) + ktree2{k}; */
          if (d > d1) {
            i1 = 0;
            i2 = 0;
            i3 = 0;
            u0 = -1;
          } else {
            i1 = (int)d - 1;
            i2 = (int)d1;
            i3 = (int)d - 1;
            u0 = (int)d1 - 1;
          }
          if (i2 - i1 == ktree2_data[(int)k - 1].f1->size[0]) {
            u1 = (u0 - i3) + 1;
            i2 = b_ktree->size[0] * b_ktree->size[1];
            b_ktree->size[0] = 1;
            b_ktree->size[1] = u1;
            emxEnsureCapacity_real_T(b_ktree, i2);
            a_data = b_ktree->data;
            for (i2 = 0; i2 < u1; i2++) {
              a_data[i2] =
                  kweig2_data[i1 + i2] + ktree2_data[(int)k - 1].f1->data[i2];
            }
            loop_ub = b_ktree->size[1];
            for (i1 = 0; i1 < loop_ub; i1++) {
              kweig2_data[i3 + i1] = a_data[i1];
            }
          } else {
            m_binary_expand_op(kweig2, i3, i1, i2 - 1, ktree2, k, u0, i3 - 1);
            kweig2_data = kweig2->data;
          }
        }
      }
    }
  }
  emxFree_real_T(&b_r);
  emxFree_real_T(&sPWM_tmp);
  emxFree_real_T(&a);
  emxFree_cell_wrap_12(&p);
  emxFree_real_T(&sPWM2);
  emxFree_real_T(&sPWM);
  emxFree_real_T(&indvec);
  emxFree_uint32_T(&X);
  emxFree_cell_wrap_3(&ktree2);
  emxFree_cell_wrap_3(&ktree);
  emxFree_real_T(&mat2);
  /* 'PWM2kmers:85' alen = length(c)-rcnum; */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    u1 = 0;
  } else {
    u0 = c->size[0];
    u1 = c->size[1];
    if (u0 >= u1) {
      u1 = u0;
    }
  }
  alen = (double)u1 - rcnum;
  /* 'PWM2kmers:86' kweig(4^k*alen+1:end) = kweig(4^k*alen+1:end)/sqrt(2); */
  d = n_tmp * alen + 1.0;
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
  u1 = (i2 - i1) - 1;
  i2 = b_ktree->size[0] * b_ktree->size[1];
  b_ktree->size[0] = 1;
  b_ktree->size[1] = u1;
  emxEnsureCapacity_real_T(b_ktree, i2);
  a_data = b_ktree->data;
  for (i2 = 0; i2 < u1; i2++) {
    a_data[i2] = kweig_data[(i + i2) - 1] / 1.4142135623730951;
  }
  loop_ub = b_ktree->size[1];
  for (i = 0; i < loop_ub; i++) {
    kweig_data[(i1 + i) + 1] = a_data[i];
  }
  /* 'PWM2kmers:87' kweig2(4^k*alen+1:end) = kweig2(4^k*alen+1:end)/sqrt(2); */
  d = pow(4.0, k) * alen + 1.0;
  if (d > kweig2->size[0]) {
    i = 1;
  } else {
    i = (int)d;
  }
  d = pow(4.0, k) * alen + 1.0;
  if (d > kweig2->size[0]) {
    i1 = -1;
    i2 = 0;
  } else {
    i1 = (int)d - 2;
    i2 = kweig2->size[0];
  }
  u1 = (i2 - i1) - 1;
  i2 = b_ktree->size[0] * b_ktree->size[1];
  b_ktree->size[0] = 1;
  b_ktree->size[1] = u1;
  emxEnsureCapacity_real_T(b_ktree, i2);
  a_data = b_ktree->data;
  for (i2 = 0; i2 < u1; i2++) {
    a_data[i2] = kweig2_data[(i + i2) - 1] / 1.4142135623730951;
  }
  loop_ub = b_ktree->size[1];
  for (i = 0; i < loop_ub; i++) {
    kweig2_data[(i1 + i) + 1] = a_data[i];
  }
  emxFree_real_T(&b_ktree);
  /* 'PWM2kmers:88' C = kweig'*kweig2/sqrt(kweig'*kweig)/sqrt(kweig2'*kweig2);
   */
  /* 'PWM2kmers:89' kweig = kweig+kweig2; */
  if (kweig->size[0] == kweig2->size[0]) {
    loop_ub = kweig->size[0];
    for (i = 0; i < loop_ub; i++) {
      kweig_data[i] += kweig2_data[i];
    }
  } else {
    b_plus(kweig, kweig2);
  }
  emxFree_real_T(&kweig2);
}

void l_binary_expand_op(emxArray_cell_wrap_3 *ktree2, double i3,
                        const emxArray_real_T *r1,
                        const emxArray_real_T *sPWM_tmp)
{
  cell_wrap_3 *ktree2_data;
  const double *b_r;
  const double *sPWM_tmp_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  sPWM_tmp_data = sPWM_tmp->data;
  b_r = r1->data;
  ktree2_data = ktree2->data;
  i = ktree2_data[(int)i3 - 1].f1->size[0];
  if (sPWM_tmp->size[0] == 1) {
    ktree2_data[(int)i3 - 1].f1->size[0] = r1->size[0];
  } else {
    ktree2_data[(int)i3 - 1].f1->size[0] = sPWM_tmp->size[0];
  }
  emxEnsureCapacity_real_T(ktree2_data[(int)i3 - 1].f1, i);
  stride_0_0 = (r1->size[0] != 1);
  stride_1_0 = (sPWM_tmp->size[0] != 1);
  if (sPWM_tmp->size[0] == 1) {
    loop_ub = r1->size[0];
  } else {
    loop_ub = sPWM_tmp->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    ktree2_data[(int)i3 - 1].f1->data[i] =
        b_r[i * stride_0_0] * sPWM_tmp_data[i * stride_1_0];
  }
}

void m_binary_expand_op(emxArray_real_T *kweig2, int i1, int i2, int i3,
                        const emxArray_cell_wrap_3 *ktree2, double k, int i4,
                        int i5)
{
  const cell_wrap_3 *ktree2_data;
  emxArray_real_T *b_kweig2;
  double *b_kweig2_data;
  double *kweig2_data;
  int i;
  int stride_0_1;
  int stride_1_1;
  int unnamed_idx_1;
  ktree2_data = ktree2->data;
  kweig2_data = kweig2->data;
  emxInit_real_T(&b_kweig2, 2);
  unnamed_idx_1 = i4 - i5;
  i = b_kweig2->size[0] * b_kweig2->size[1];
  b_kweig2->size[0] = 1;
  b_kweig2->size[1] = unnamed_idx_1;
  emxEnsureCapacity_real_T(b_kweig2, i);
  b_kweig2_data = b_kweig2->data;
  stride_0_1 = ((i3 - i2) + 1 != 1);
  stride_1_1 = (ktree2_data[(int)k - 1].f1->size[0] != 1);
  for (i = 0; i < unnamed_idx_1; i++) {
    b_kweig2_data[i] = kweig2_data[i2 + i * stride_0_1] +
                       ktree2_data[(int)k - 1].f1->data[i * stride_1_1];
  }
  unnamed_idx_1 = b_kweig2->size[1];
  for (i = 0; i < unnamed_idx_1; i++) {
    kweig2_data[i1 + i] = b_kweig2_data[i];
  }
  emxFree_real_T(&b_kweig2);
}

/* End of code generation (PWM2kmers.c) */
