/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * genIndex.c
 *
 * Code generation for function 'genIndex'
 *
 */

/* Include files */
#include "genIndex.h"
#include "find.h"
#include "flip.h"
#include "getgkmweights_data.h"
#include "getgkmweights_emxutil.h"
#include "getgkmweights_rtwutil.h"
#include "getgkmweights_types.h"
#include "minOrMax.h"
#include "mod.h"
#include "nchoosek.h"
#include "nullAssignment.h"
#include "rand.h"
#include "sort.h"
#include "sum.h"
#include <math.h>
#include <stdio.h>

/* Function Declarations */
static void b_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               const emxArray_real_T *r1,
                               const emxArray_real_T *mat);

static void f_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               unsigned int j, const emxArray_real_T *vec);

static void g_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               int i, const emxArray_real_T *vec);

static void genIndex_frac(double l, double k, const emxArray_real_T *c,
                          double Lfrac, emxArray_real_T *mat);

static void h_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *C,
                               const emxArray_int32_T *f,
                               const emxArray_real_T *d, int i);

static void i_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               const emxArray_int32_T *ind,
                               const emxArray_int32_T *f, int i,
                               const emxArray_real_T *C);

static void sort_comb(emxArray_real_T *c, double l, double k, double n_frac,
                      emxArray_real_T *C, emxArray_real_T *b_I,
                      emxArray_real_T *ind, emxArray_real_T *mat,
                      double *rcnum);

/* Function Definitions */
static void b_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               const emxArray_real_T *r1,
                               const emxArray_real_T *mat)
{
  const double *c_data;
  const double *mat_data;
  const double *r;
  double *x_data;
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
  mat_data = mat->data;
  r = r1->data;
  c_data = c->data;
  i = c->size[1];
  i1 = x->size[0] * x->size[1];
  if (mat->size[0] == 1) {
    x->size[0] = r1->size[0];
  } else {
    x->size[0] = mat->size[0];
  }
  if (mat->size[1] == 1) {
    x->size[1] = i;
  } else {
    x->size[1] = mat->size[1];
  }
  emxEnsureCapacity_real_T(x, i1);
  x_data = x->data;
  stride_0_0 = (r1->size[0] != 1);
  stride_0_1 = (i != 1);
  stride_1_0 = (mat->size[0] != 1);
  stride_1_1 = (mat->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (mat->size[1] == 1) {
    loop_ub = i;
  } else {
    loop_ub = mat->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    if (mat->size[0] == 1) {
      b_loop_ub = r1->size[0];
    } else {
      b_loop_ub = mat->size[0];
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      x_data[i1 + x->size[0] * i] =
          c_data[((int)r[i1 * stride_0_0] + c->size[0] * aux_0_1) - 1] -
          mat_data[i1 * stride_1_0 + mat->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
}

static void f_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               unsigned int j, const emxArray_real_T *vec)
{
  const double *c_data;
  const double *vec_data;
  double *x_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  vec_data = vec->data;
  c_data = c->data;
  i = c->size[1];
  stride_0_1 = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (vec->size[1] == 1) {
    x->size[1] = i;
  } else {
    x->size[1] = vec->size[1];
  }
  emxEnsureCapacity_real_T(x, stride_0_1);
  x_data = x->data;
  stride_0_1 = (i != 1);
  stride_1_1 = (vec->size[1] != 1);
  if (vec->size[1] == 1) {
    loop_ub = i;
  } else {
    loop_ub = vec->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    x_data[i] = c_data[((int)j + c->size[0] * (i * stride_0_1)) - 1] -
                vec_data[i * stride_1_1];
  }
}

static void g_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               int i, const emxArray_real_T *vec)
{
  const double *c_data;
  const double *vec_data;
  double *x_data;
  int b_i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  vec_data = vec->data;
  c_data = c->data;
  b_i = c->size[1];
  stride_0_1 = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (vec->size[1] == 1) {
    x->size[1] = b_i;
  } else {
    x->size[1] = vec->size[1];
  }
  emxEnsureCapacity_real_T(x, stride_0_1);
  x_data = x->data;
  stride_0_1 = (b_i != 1);
  stride_1_1 = (vec->size[1] != 1);
  if (vec->size[1] == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = vec->size[1];
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    x_data[b_i] = c_data[i + c->size[0] * (b_i * stride_0_1)] -
                  vec_data[b_i * stride_1_1];
  }
}

/*
 * function mat = genIndex_frac(l,k,c,rc,Lfrac)
 */
static void genIndex_frac(double l, double k, const emxArray_real_T *c,
                          double Lfrac, emxArray_real_T *mat)
{
  emxArray_boolean_T *b_f;
  emxArray_boolean_T *b_x;
  emxArray_int32_T *c_f;
  emxArray_real_T *C;
  emxArray_real_T *X;
  emxArray_real_T *cvec;
  emxArray_real_T *f;
  emxArray_real_T *ind;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  emxArray_real_T *x;
  emxArray_real_T *y;
  const double *c_data;
  double L;
  double M;
  double b_r;
  double xtmp;
  double *C_data;
  double *X_data;
  double *b_x_data;
  double *cvec_data;
  double *ind_data;
  double *mat_data;
  double *r2;
  int b_k;
  unsigned int count;
  int i;
  int j;
  int nxin;
  int nxout;
  int nz;
  int vlen;
  int *f_data;
  bool exitg1;
  bool *x_data;
  c_data = c->data;
  /* 'genIndex:158' L = numel(c)/k; */
  L = (double)(c->size[0] * c->size[1]) / k;
  /* 'genIndex:159' x = zeros(L, l); */
  /* 'genIndex:160' for i = 1:L */
  /* 'genIndex:164' M = l/2; */
  M = l / 2.0;
  /* 'genIndex:165' y = sum(x==2,2); */
  /* 'genIndex:166' nind = true; */
  /* 'genIndex:167' ind = 1:L; */
  emxInit_real_T(&ind, 2);
  ind_data = ind->data;
  if (L < 1.0) {
    ind->size[0] = 1;
    ind->size[1] = 0;
  } else {
    i = ind->size[0] * ind->size[1];
    ind->size[0] = 1;
    nxin = (int)floor(L - 1.0);
    ind->size[1] = nxin + 1;
    emxEnsureCapacity_real_T(ind, i);
    ind_data = ind->data;
    for (i = 0; i <= nxin; i++) {
      ind_data[i] = (double)i + 1.0;
    }
  }
  emxInit_real_T(&X, 1);
  /* 'genIndex:168' X = zeros(L,1); */
  i = X->size[0];
  nxin = (int)L;
  X->size[0] = (int)L;
  emxEnsureCapacity_real_T(X, i);
  X_data = X->data;
  for (i = 0; i < nxin; i++) {
    X_data[i] = 0.0;
  }
  /* 'genIndex:169' count = 0; */
  count = 0U;
  /* 'genIndex:170' i = 0; */
  /* 'genIndex:171' while count == 0 || count < Lfrac */
  emxInit_real_T(&x, 2);
  emxInit_real_T(&cvec, 2);
  emxInit_real_T(&C, 2);
  emxInit_real_T(&f, 1);
  emxInit_real_T(&r, 2);
  emxInit_boolean_T(&b_x, 2);
  emxInit_real_T(&y, 2);
  emxInit_real_T(&r1, 1);
  emxInit_boolean_T(&b_f, 1);
  emxInit_int32_T(&c_f, 1);
  while ((((int)count == 0) || (count < Lfrac)) &&
         (Lfrac - (double)count >= ((double)count + M) - Lfrac)) {
    /* 'genIndex:172' if Lfrac-count < count+M-Lfrac */
    /* 'genIndex:175' count = count + 1; */
    count++;
    /* 'genIndex:176' r = randi(L-count); */
    b_r = b_rand();
    b_r = floor(b_r * (L - (double)count));
    /* 'genIndex:177' cvec = c(ind(r),:); */
    xtmp = ind_data[(int)(b_r + 1.0) - 1];
    nxin = c->size[1];
    i = cvec->size[0] * cvec->size[1];
    cvec->size[0] = 1;
    cvec->size[1] = c->size[1];
    emxEnsureCapacity_real_T(cvec, i);
    cvec_data = cvec->data;
    for (i = 0; i < nxin; i++) {
      cvec_data[i] = c_data[((int)xtmp + c->size[0] * i) - 1];
    }
    /* 'genIndex:178' X(count) = ind(r); */
    X_data[(int)count - 1] = xtmp;
    /* 'genIndex:179' ind(r) = []; */
    vlen = (int)(b_r + 1.0);
    nxin = ind->size[1];
    nxout = ind->size[1] - 1;
    for (b_k = vlen; b_k <= nxout; b_k++) {
      ind_data[b_k - 1] = ind_data[b_k];
    }
    i = ind->size[0] * ind->size[1];
    if (1 > nxout) {
      ind->size[1] = 0;
    } else {
      ind->size[1] = nxin - 1;
    }
    emxEnsureCapacity_real_T(ind, i);
    ind_data = ind->data;
    /* 'genIndex:180' num = sum(cvec<M+1); */
    i = b_x->size[0] * b_x->size[1];
    b_x->size[0] = 1;
    b_x->size[1] = cvec->size[1];
    emxEnsureCapacity_boolean_T(b_x, i);
    x_data = b_x->data;
    nxin = cvec->size[1];
    for (i = 0; i < nxin; i++) {
      x_data[i] = (cvec_data[i] < M + 1.0);
    }
    vlen = b_x->size[1];
    if (b_x->size[1] == 0) {
      nz = 0;
    } else {
      nz = x_data[0];
      for (b_k = 2; b_k <= vlen; b_k++) {
        nz += x_data[b_k - 1];
      }
    }
    /* 'genIndex:181' for j = 1:M-1 */
    j = 0;
    exitg1 = false;
    while ((!exitg1) && (j <= (int)(M - 1.0) - 1)) {
      /* 'genIndex:182' C = zeros(1,k); */
      i = C->size[0] * C->size[1];
      C->size[0] = 1;
      vlen = (int)k;
      C->size[1] = (int)k;
      emxEnsureCapacity_real_T(C, i);
      C_data = C->data;
      for (i = 0; i < vlen; i++) {
        C_data[i] = 0.0;
      }
      /* 'genIndex:183' C(1:num) = sort(mod(cvec(1:num)+(j-1),M)+1); */
      if (1 > nz) {
        nxin = 0;
      } else {
        nxin = nz;
      }
      i = r->size[0] * r->size[1];
      r->size[0] = 1;
      r->size[1] = nxin;
      emxEnsureCapacity_real_T(r, i);
      r2 = r->data;
      for (i = 0; i < nxin; i++) {
        b_r = cvec_data[i] + (((double)j + 1.0) - 1.0);
        r2[i] = b_mod(b_r, M);
      }
      i = r->size[0] * r->size[1];
      r->size[0] = 1;
      emxEnsureCapacity_real_T(r, i);
      r2 = r->data;
      nxin = r->size[1] - 1;
      for (i = 0; i <= nxin; i++) {
        r2[i]++;
      }
      b_sort(r);
      r2 = r->data;
      nxin = r->size[1];
      for (i = 0; i < nxin; i++) {
        C_data[i] = r2[i];
      }
      /* 'genIndex:184' C(num+1:end) = sort(l-mod(l+1-cvec(num+1:end)+(j-1),M));
       */
      if ((double)nz + 1.0 > cvec->size[1]) {
        i = 0;
        b_k = 0;
      } else {
        i = nz;
        b_k = cvec->size[1];
      }
      if ((double)nz + 1.0 > C->size[1]) {
        nxout = 0;
      } else {
        nxout = nz;
      }
      vlen = r->size[0] * r->size[1];
      r->size[0] = 1;
      nxin = b_k - i;
      r->size[1] = nxin;
      emxEnsureCapacity_real_T(r, vlen);
      r2 = r->data;
      for (b_k = 0; b_k < nxin; b_k++) {
        b_r = ((l + 1.0) - cvec_data[i + b_k]) + (((double)j + 1.0) - 1.0);
        r2[b_k] = b_mod(b_r, M);
      }
      i = r->size[0] * r->size[1];
      r->size[0] = 1;
      emxEnsureCapacity_real_T(r, i);
      r2 = r->data;
      nxin = r->size[1] - 1;
      for (i = 0; i <= nxin; i++) {
        r2[i] = l - r2[i];
      }
      b_sort(r);
      r2 = r->data;
      nxin = r->size[1];
      for (i = 0; i < nxin; i++) {
        C_data[nxout + i] = r2[i];
      }
      /* 'genIndex:185' mat = repmat(C, L-count,1); */
      b_r = L - (double)count;
      i = (int)b_r;
      b_k = mat->size[0] * mat->size[1];
      mat->size[0] = (int)b_r;
      mat->size[1] = C->size[1];
      emxEnsureCapacity_real_T(mat, b_k);
      mat_data = mat->data;
      vlen = C->size[1];
      for (nxout = 0; nxout < vlen; nxout++) {
        nxin = nxout * (int)b_r;
        for (b_k = 0; b_k < i; b_k++) {
          mat_data[nxin + b_k] = C_data[nxout];
        }
      }
      /* 'genIndex:186' f = find(sum(abs(c(ind,:)-mat),2)==0); */
      b_k = r1->size[0];
      r1->size[0] = ind->size[1];
      emxEnsureCapacity_real_T(r1, b_k);
      r2 = r1->data;
      nxin = ind->size[1];
      for (b_k = 0; b_k < nxin; b_k++) {
        r2[b_k] = ind_data[b_k];
      }
      if ((r1->size[0] == mat->size[0]) && (c->size[1] == mat->size[1])) {
        nxin = c->size[1];
        b_k = x->size[0] * x->size[1];
        x->size[0] = r1->size[0];
        x->size[1] = c->size[1];
        emxEnsureCapacity_real_T(x, b_k);
        b_x_data = x->data;
        for (b_k = 0; b_k < nxin; b_k++) {
          vlen = r1->size[0];
          for (nxout = 0; nxout < vlen; nxout++) {
            b_x_data[nxout + x->size[0] * b_k] =
                c_data[((int)r2[nxout] + c->size[0] * b_k) - 1] -
                mat_data[nxout + mat->size[0] * b_k];
          }
        }
      } else {
        b_binary_expand_op(x, c, r1, mat);
        b_x_data = x->data;
      }
      vlen = x->size[0] * x->size[1];
      b_k = y->size[0] * y->size[1];
      y->size[0] = x->size[0];
      y->size[1] = x->size[1];
      emxEnsureCapacity_real_T(y, b_k);
      mat_data = y->data;
      for (b_k = 0; b_k < vlen; b_k++) {
        mat_data[b_k] = fabs(b_x_data[b_k]);
      }
      b_sum(y, f);
      mat_data = f->data;
      b_k = b_f->size[0];
      b_f->size[0] = f->size[0];
      emxEnsureCapacity_boolean_T(b_f, b_k);
      x_data = b_f->data;
      nxin = f->size[0];
      for (b_k = 0; b_k < nxin; b_k++) {
        x_data[b_k] = (mat_data[b_k] == 0.0);
      }
      b_eml_find(b_f, c_f);
      f_data = c_f->data;
      b_k = f->size[0];
      f->size[0] = c_f->size[0];
      emxEnsureCapacity_real_T(f, b_k);
      mat_data = f->data;
      nxin = c_f->size[0];
      for (b_k = 0; b_k < nxin; b_k++) {
        mat_data[b_k] = f_data[b_k];
      }
      /* 'genIndex:187' if isempty(f) */
      if (f->size[0] == 0) {
        /* 'genIndex:188' R = l+1-fliplr(C); */
        vlen = C->size[1] - 1;
        nxin = C->size[1] >> 1;
        for (b_k = 0; b_k < nxin; b_k++) {
          nxout = vlen - b_k;
          xtmp = C_data[b_k];
          C_data[b_k] = C_data[nxout];
          C_data[nxout] = xtmp;
        }
        /* 'genIndex:189' mat = repmat(R,L-count,1); */
        b_k = C->size[0] * C->size[1];
        C->size[0] = 1;
        emxEnsureCapacity_real_T(C, b_k);
        C_data = C->data;
        nxin = C->size[1] - 1;
        for (b_k = 0; b_k <= nxin; b_k++) {
          C_data[b_k] = (l + 1.0) - C_data[b_k];
        }
        b_k = mat->size[0] * mat->size[1];
        mat->size[0] = (int)b_r;
        mat->size[1] = C->size[1];
        emxEnsureCapacity_real_T(mat, b_k);
        mat_data = mat->data;
        vlen = C->size[1];
        for (nxout = 0; nxout < vlen; nxout++) {
          nxin = nxout * (int)b_r;
          for (b_k = 0; b_k < i; b_k++) {
            mat_data[nxin + b_k] = C_data[nxout];
          }
        }
        /* 'genIndex:190' f = find(sum(abs(c(ind,:)-mat),2)==0); */
        if ((r1->size[0] == mat->size[0]) && (c->size[1] == mat->size[1])) {
          nxin = c->size[1];
          i = x->size[0] * x->size[1];
          x->size[0] = r1->size[0];
          x->size[1] = c->size[1];
          emxEnsureCapacity_real_T(x, i);
          b_x_data = x->data;
          for (i = 0; i < nxin; i++) {
            vlen = r1->size[0];
            for (b_k = 0; b_k < vlen; b_k++) {
              b_x_data[b_k + x->size[0] * i] =
                  c_data[((int)r2[b_k] + c->size[0] * i) - 1] -
                  mat_data[b_k + mat->size[0] * i];
            }
          }
        } else {
          b_binary_expand_op(x, c, r1, mat);
          b_x_data = x->data;
        }
        vlen = x->size[0] * x->size[1];
        i = y->size[0] * y->size[1];
        y->size[0] = x->size[0];
        y->size[1] = x->size[1];
        emxEnsureCapacity_real_T(y, i);
        mat_data = y->data;
        for (b_k = 0; b_k < vlen; b_k++) {
          mat_data[b_k] = fabs(b_x_data[b_k]);
        }
        b_sum(y, r1);
        r2 = r1->data;
        i = b_f->size[0];
        b_f->size[0] = r1->size[0];
        emxEnsureCapacity_boolean_T(b_f, i);
        x_data = b_f->data;
        nxin = r1->size[0];
        for (i = 0; i < nxin; i++) {
          x_data[i] = (r2[i] == 0.0);
        }
        b_eml_find(b_f, c_f);
        f_data = c_f->data;
        i = f->size[0];
        f->size[0] = c_f->size[0];
        emxEnsureCapacity_real_T(f, i);
        mat_data = f->data;
        nxin = c_f->size[0];
        for (i = 0; i < nxin; i++) {
          mat_data[i] = f_data[i];
        }
        /* 'genIndex:191' if isempty(f) */
        if (f->size[0] == 0) {
          exitg1 = true;
        } else {
          /* 'genIndex:193' else */
          /* 'genIndex:194' count = count + 1; */
          count++;
          /* 'genIndex:195' X(count)=ind(f); */
          X_data[(int)count - 1] = ind_data[(int)mat_data[0] - 1];
          /* 'genIndex:196' ind(f)=[]; */
          i = c_f->size[0];
          c_f->size[0] = f->size[0];
          emxEnsureCapacity_int32_T(c_f, i);
          f_data = c_f->data;
          nxin = f->size[0];
          for (i = 0; i < nxin; i++) {
            f_data[i] = (int)mat_data[i];
          }
          nullAssignment(ind, c_f);
          ind_data = ind->data;
          j++;
        }
      } else {
        /* 'genIndex:198' else */
        /* 'genIndex:199' count = count + 1; */
        count++;
        /* 'genIndex:200' X(count)=ind(f); */
        X_data[(int)count - 1] = ind_data[(int)mat_data[0] - 1];
        /* 'genIndex:201' ind(f)=[]; */
        i = c_f->size[0];
        c_f->size[0] = f->size[0];
        emxEnsureCapacity_int32_T(c_f, i);
        f_data = c_f->data;
        nxin = f->size[0];
        for (i = 0; i < nxin; i++) {
          f_data[i] = (int)mat_data[i];
        }
        nullAssignment(ind, c_f);
        ind_data = ind->data;
        j++;
      }
    }
  }
  emxFree_int32_T(&c_f);
  emxFree_boolean_T(&b_f);
  emxFree_real_T(&r1);
  emxFree_real_T(&y);
  emxFree_boolean_T(&b_x);
  emxFree_real_T(&r);
  emxFree_real_T(&f);
  emxFree_real_T(&C);
  emxFree_real_T(&cvec);
  emxFree_real_T(&ind);
  emxFree_real_T(&x);
  /* 'genIndex:205' X = X(1:count); */
  i = X->size[0];
  if (1 > (int)count) {
    X->size[0] = 0;
  } else {
    X->size[0] = (int)count;
  }
  emxEnsureCapacity_real_T(X, i);
  /* 'genIndex:206' X = sort(X); */
  c_sort(X);
  X_data = X->data;
  /* 'genIndex:207' mat = c(X,:); */
  nxin = c->size[1];
  i = mat->size[0] * mat->size[1];
  mat->size[0] = X->size[0];
  mat->size[1] = c->size[1];
  emxEnsureCapacity_real_T(mat, i);
  mat_data = mat->data;
  for (i = 0; i < nxin; i++) {
    vlen = X->size[0];
    for (b_k = 0; b_k < vlen; b_k++) {
      mat_data[b_k + mat->size[0] * i] =
          c_data[((int)X_data[b_k] + c->size[0] * i) - 1];
    }
  }
  emxFree_real_T(&X);
}

static void h_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *C,
                               const emxArray_int32_T *f,
                               const emxArray_real_T *d, int i)
{
  emxArray_real_T *b_C;
  const double *C_data;
  const double *d_data;
  double *b_C_data;
  double *x_data;
  const int *f_data;
  int b_f;
  int b_i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  d_data = d->data;
  f_data = f->data;
  C_data = C->data;
  x_data = x->data;
  emxInit_real_T(&b_C, 2);
  b_f = f_data[(int)d_data[i] - 1];
  b_i = C->size[1];
  stride_0_1 = b_C->size[0] * b_C->size[1];
  b_C->size[0] = 1;
  if (x->size[1] == 1) {
    b_C->size[1] = b_i;
  } else {
    b_C->size[1] = x->size[1];
  }
  emxEnsureCapacity_real_T(b_C, stride_0_1);
  b_C_data = b_C->data;
  stride_0_1 = (b_i != 1);
  stride_1_1 = (x->size[1] != 1);
  if (x->size[1] == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = x->size[1];
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    b_C_data[b_i] =
        (C_data[(b_f + C->size[0] * (b_i * stride_0_1)) - 1] + 1.0) -
        x_data[b_i * stride_1_1];
  }
  b_i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = b_C->size[1];
  emxEnsureCapacity_real_T(x, b_i);
  x_data = x->data;
  loop_ub = b_C->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    x_data[b_i] = b_C_data[b_i];
  }
  emxFree_real_T(&b_C);
}

static void i_binary_expand_op(emxArray_real_T *x, const emxArray_real_T *c,
                               const emxArray_int32_T *ind,
                               const emxArray_int32_T *f, int i,
                               const emxArray_real_T *C)
{
  const double *C_data;
  const double *c_data;
  double *x_data;
  const int *f_data;
  const int *ind_data;
  int b_i;
  int b_ind;
  int i1;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  C_data = C->data;
  f_data = f->data;
  ind_data = ind->data;
  c_data = c->data;
  b_ind = ind_data[f_data[i] - 1];
  b_i = c->size[1];
  i1 = C->size[1];
  loop_ub = x->size[0] * x->size[1];
  x->size[0] = 1;
  if (i1 == 1) {
    x->size[1] = b_i;
  } else {
    x->size[1] = i1;
  }
  emxEnsureCapacity_real_T(x, loop_ub);
  x_data = x->data;
  stride_0_1 = (b_i != 1);
  stride_1_1 = (i1 != 1);
  if (i1 == 1) {
    loop_ub = b_i;
  } else {
    loop_ub = i1;
  }
  for (b_i = 0; b_i < loop_ub; b_i++) {
    x_data[b_i] = c_data[(b_ind + c->size[0] * (b_i * stride_0_1)) - 1] -
                  C_data[i + C->size[0] * (b_i * stride_1_1)];
  }
}

/*
 * function [c,C,I,ind,mat,rcnum] = sort_comb(c,l,k,n_frac)
 */
static void sort_comb(emxArray_real_T *c, double l, double k, double n_frac,
                      emxArray_real_T *C, emxArray_real_T *b_I,
                      emxArray_real_T *ind, emxArray_real_T *mat, double *rcnum)
{
  emxArray_boolean_T *b_d;
  emxArray_boolean_T *b_e;
  emxArray_int32_T *b_ind;
  emxArray_int32_T *e;
  emxArray_int32_T *f;
  emxArray_int32_T *iidx;
  emxArray_int32_T *r;
  emxArray_real_T *d;
  emxArray_real_T *varargin_1;
  emxArray_real_T *varargin_2;
  emxArray_real_T *vec;
  emxArray_real_T *x;
  emxArray_real_T *y;
  double L;
  double a;
  double xtmp;
  double *C_data;
  double *c_data;
  double *d_data;
  double *vec_data;
  double *x_data;
  double *y_data;
  int b_i;
  unsigned int b_j;
  int b_j1;
  int exitg1;
  int i;
  int i1;
  int j;
  int j2;
  int loop_ub;
  int loop_ub_tmp;
  int nd2;
  int nx;
  int *e_data;
  int *f_data;
  int *iidx_data;
  signed char input_sizes_idx_1;
  bool empty_non_axis_sizes;
  bool *b_e_data;
  c_data = c->data;
  emxInit_real_T(&d, 1);
  /* 'genIndex:58' L = numel(c)/k; */
  L = (double)(c->size[0] * c->size[1]) / k;
  /* 'genIndex:59' d = ones(L,1); */
  i = d->size[0];
  loop_ub = (int)L;
  d->size[0] = (int)L;
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  for (i = 0; i < loop_ub; i++) {
    d_data[i] = 1.0;
  }
  emxInit_int32_T(&e, 2);
  /*  Alternative to continual growth */
  /*  e = []; */
  /* 'genIndex:63' e = ones(1,length(d)-1) * -1; */
  i = e->size[0] * e->size[1];
  e->size[0] = 1;
  nx = (int)L - 1;
  e->size[1] = (int)L - 1;
  emxEnsureCapacity_int32_T(e, i);
  e_data = e->data;
  for (i = 0; i < nx; i++) {
    e_data[i] = -1;
  }
  /* 'genIndex:65' a = 1; */
  a = 1.0;
  /* 'genIndex:66' for i = 1:length(d)-1 */
  emxInit_real_T(&vec, 2);
  emxInit_real_T(&y, 2);
  y_data = y->data;
  emxInit_real_T(&x, 2);
  for (b_i = 0; b_i <= loop_ub - 2; b_i++) {
    /* 'genIndex:67' vec = l+1-fliplr(c(i,:)); */
    nx = c->size[1];
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    vec->size[1] = nx;
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    for (i = 0; i < nx; i++) {
      vec_data[i] = c_data[b_i + c->size[0] * i];
    }
    nd2 = c->size[1] >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (c->size[1] - b_j1) - 1;
      xtmp = vec_data[b_j1];
      vec_data[b_j1] = vec_data[j2];
      vec_data[j2] = xtmp;
    }
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    nx = vec->size[1] - 1;
    for (i = 0; i <= nx; i++) {
      vec_data[i] = (l + 1.0) - vec_data[i];
    }
    /* 'genIndex:68' if d(i) ~= 0 */
    if (d_data[b_i] != 0.0) {
      /* 'genIndex:69' if sum(abs(c(i,:)-vec))==0 */
      if (c->size[1] == vec->size[1]) {
        nx = c->size[1];
        i = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = nx;
        emxEnsureCapacity_real_T(x, i);
        x_data = x->data;
        for (i = 0; i < nx; i++) {
          x_data[i] = c_data[b_i + c->size[0] * i] - vec_data[i];
        }
      } else {
        g_binary_expand_op(x, c, b_i, vec);
        x_data = x->data;
      }
      nx = x->size[1];
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity_real_T(y, i);
      y_data = y->data;
      for (j2 = 0; j2 < nx; j2++) {
        y_data[j2] = fabs(x_data[j2]);
      }
      if (sum(y) == 0.0) {
        /* 'genIndex:70' e(a) = i; */
        e_data[(int)a - 1] = b_i + 1;
        /* 'genIndex:71' d(i) = 0; */
        d_data[b_i] = 0.0;
        /* 'genIndex:72' a = a+1; */
        a++;
      } else {
        /* 'genIndex:73' else */
        /* 'genIndex:74' for j = i+1:length(d) */
        i = d->size[0] - b_i;
        for (j = 0; j <= i - 2; j++) {
          b_j = ((unsigned int)b_i + j) + 2U;
          /* 'genIndex:75' if sum(abs(c(j,:)-vec))==0 */
          if (c->size[1] == vec->size[1]) {
            nx = c->size[1];
            i1 = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = nx;
            emxEnsureCapacity_real_T(x, i1);
            x_data = x->data;
            for (i1 = 0; i1 < nx; i1++) {
              x_data[i1] =
                  c_data[((int)b_j + c->size[0] * i1) - 1] - vec_data[i1];
            }
          } else {
            f_binary_expand_op(x, c, b_j, vec);
            x_data = x->data;
          }
          nx = x->size[1];
          i1 = y->size[0] * y->size[1];
          y->size[0] = 1;
          y->size[1] = x->size[1];
          emxEnsureCapacity_real_T(y, i1);
          y_data = y->data;
          for (j2 = 0; j2 < nx; j2++) {
            y_data[j2] = fabs(x_data[j2]);
          }
          if (sum(y) == 0.0) {
            /* 'genIndex:76' d(j) = 0; */
            d_data[(int)b_j - 1] = 0.0;
          }
        }
      }
    }
  }
  emxInit_boolean_T(&b_e, 2);
  /*  Remove all trailing -1 */
  /* 'genIndex:83' first_occ = find(e==-1); */
  i = b_e->size[0] * b_e->size[1];
  b_e->size[0] = 1;
  b_e->size[1] = e->size[1];
  emxEnsureCapacity_boolean_T(b_e, i);
  b_e_data = b_e->data;
  loop_ub = e->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (e_data[i] == -1);
  }
  emxInit_int32_T(&r, 2);
  eml_find(b_e, r);
  e_data = r->data;
  i = vec->size[0] * vec->size[1];
  vec->size[0] = 1;
  vec->size[1] = r->size[1];
  emxEnsureCapacity_real_T(vec, i);
  vec_data = vec->data;
  loop_ub = r->size[1];
  emxFree_boolean_T(&b_e);
  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = e_data[i];
  }
  emxFree_int32_T(&r);
  emxInit_boolean_T(&b_d, 1);
  /* 'genIndex:84' e = e(1:first_occ(1)-1); */
  if (1.0 > (double)(int)vec_data[0] - 1.0) {
    nd2 = 0;
  } else {
    nd2 = (int)vec_data[0] - 1;
  }
  i = e->size[0] * e->size[1];
  if (1.0 > (double)(int)vec_data[0] - 1.0) {
    e->size[1] = 0;
  } else {
    e->size[1] = (int)vec_data[0] - 1;
  }
  emxEnsureCapacity_int32_T(e, i);
  e_data = e->data;
  /* 'genIndex:86' rcnum = a-1; */
  *rcnum = a - 1.0;
  /* 'genIndex:87' f = find(d==1); */
  i = b_d->size[0];
  b_d->size[0] = d->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  loop_ub = d->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (d_data[i] == 1.0);
  }
  emxInit_int32_T(&f, 1);
  emxInit_int32_T(&iidx, 1);
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = iidx_data[i];
  }
  /* 'genIndex:88' c=[c(f,:);c(e,:)]; */
  i = b_d->size[0];
  b_d->size[0] = d->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  loop_ub = d->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (d_data[i] == 1.0);
  }
  emxInit_real_T(&varargin_1, 2);
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  loop_ub = c->size[1];
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = iidx->size[0];
  varargin_1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(varargin_1, i);
  vec_data = varargin_1->data;
  for (i = 0; i < loop_ub; i++) {
    nx = iidx->size[0];
    for (i1 = 0; i1 < nx; i1++) {
      vec_data[i1 + varargin_1->size[0] * i] =
          c_data[(iidx_data[i1] + c->size[0] * i) - 1];
    }
  }
  emxInit_real_T(&varargin_2, 2);
  loop_ub = c->size[1];
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = nd2;
  varargin_2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(varargin_2, i);
  x_data = varargin_2->data;
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < nd2; i1++) {
      x_data[i1 + varargin_2->size[0] * i] =
          c_data[(e_data[i1] + c->size[0] * i) - 1];
    }
  }
  emxFree_int32_T(&e);
  if ((f->size[0] != 0) && (c->size[1] != 0)) {
    j2 = c->size[1];
  } else if ((nd2 != 0) && (c->size[1] != 0)) {
    j2 = c->size[1];
  } else {
    if (c->size[1] > 0) {
      j2 = c->size[1];
    } else {
      j2 = 0;
    }
    if (c->size[1] > j2) {
      j2 = c->size[1];
    }
  }
  empty_non_axis_sizes = (j2 == 0);
  if (empty_non_axis_sizes || ((f->size[0] != 0) && (c->size[1] != 0))) {
    nx = f->size[0];
  } else {
    nx = 0;
  }
  if ((!empty_non_axis_sizes) && ((nd2 == 0) || (c->size[1] == 0))) {
    nd2 = 0;
  }
  i = c->size[0] * c->size[1];
  c->size[0] = nx + nd2;
  c->size[1] = j2;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  for (i = 0; i < j2; i++) {
    for (i1 = 0; i1 < nx; i1++) {
      c_data[i1 + c->size[0] * i] = vec_data[i1 + nx * i];
    }
  }
  for (i = 0; i < j2; i++) {
    for (i1 = 0; i1 < nd2; i1++) {
      c_data[(i1 + nx) + c->size[0] * i] = x_data[i1 + nd2 * i];
    }
  }
  /* 'genIndex:89' C = c; */
  i = C->size[0] * C->size[1];
  C->size[0] = c->size[0];
  C->size[1] = c->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  loop_ub = c->size[0] * c->size[1];
  for (i = 0; i < loop_ub; i++) {
    C_data[i] = c_data[i];
  }
  /* 'genIndex:90' L = numel(c)/k; */
  L = (double)(c->size[0] * c->size[1]) / k;
  /* 'genIndex:91' for i = 1:L */
  loop_ub_tmp = (int)L;
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    /* 'genIndex:92' c2 = l+1-fliplr(c(i,:)); */
    loop_ub = c->size[1];
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    vec->size[1] = loop_ub;
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    for (i = 0; i < loop_ub; i++) {
      vec_data[i] = c_data[b_i + c->size[0] * i];
    }
    nd2 = c->size[1] >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (c->size[1] - b_j1) - 1;
      xtmp = vec_data[b_j1];
      vec_data[b_j1] = vec_data[j2];
      vec_data[j2] = xtmp;
    }
    i = vec->size[0] * vec->size[1];
    vec->size[0] = 1;
    emxEnsureCapacity_real_T(vec, i);
    vec_data = vec->data;
    loop_ub = vec->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      vec_data[i] = (l + 1.0) - vec_data[i];
    }
    /* 'genIndex:93' if sum(0.5.^c(i,:)) < sum(0.5.^c2) */
    loop_ub = c->size[1];
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = loop_ub;
    emxEnsureCapacity_real_T(y, i);
    y_data = y->data;
    for (i = 0; i < loop_ub; i++) {
      a = c_data[b_i + c->size[0] * i];
      y_data[i] = pow(0.5, a);
    }
    i = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = vec->size[1];
    emxEnsureCapacity_real_T(x, i);
    x_data = x->data;
    loop_ub = vec->size[1];
    for (i = 0; i < loop_ub; i++) {
      a = vec_data[i];
      x_data[i] = pow(0.5, a);
    }
    if (sum(y) < sum(x)) {
      /* 'genIndex:94' C(i,:) = c2; */
      loop_ub = vec->size[1];
      for (i = 0; i < loop_ub; i++) {
        C_data[b_i + C->size[0] * i] = vec_data[i];
      }
    }
  }
  emxFree_real_T(&vec);
  /* 'genIndex:97' c=C; */
  i = c->size[0] * c->size[1];
  c->size[0] = C->size[0];
  c->size[1] = C->size[1];
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  loop_ub = C->size[1];
  for (i = 0; i < loop_ub; i++) {
    nx = C->size[0];
    for (i1 = 0; i1 < nx; i1++) {
      c_data[i1 + c->size[0] * i] = C_data[i1 + C->size[0] * i];
    }
  }
  /* 'genIndex:98' [~,ind] = sort(sum(0.5.^C,2),'descend'); */
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = C->size[0];
  varargin_1->size[1] = C->size[1];
  emxEnsureCapacity_real_T(varargin_1, i);
  vec_data = varargin_1->data;
  loop_ub = C->size[0] * C->size[1];
  for (i = 0; i < loop_ub; i++) {
    a = C_data[i];
    vec_data[i] = pow(0.5, a);
  }
  emxInit_int32_T(&b_ind, 1);
  b_sum(varargin_1, d);
  sort(d, iidx);
  iidx_data = iidx->data;
  i = b_ind->size[0];
  b_ind->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(b_ind, i);
  e_data = b_ind->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    e_data[i] = iidx_data[i];
  }
  /* 'genIndex:99' C = C(ind,:); */
  nx = C->size[1] - 1;
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = b_ind->size[0];
  varargin_2->size[1] = C->size[1];
  emxEnsureCapacity_real_T(varargin_2, i);
  x_data = varargin_2->data;
  for (i = 0; i <= nx; i++) {
    loop_ub = b_ind->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      x_data[i1 + varargin_2->size[0] * i] =
          C_data[(e_data[i1] + C->size[0] * i) - 1];
    }
  }
  i = C->size[0] * C->size[1];
  C->size[0] = varargin_2->size[0];
  C->size[1] = varargin_2->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  loop_ub = varargin_2->size[0] * varargin_2->size[1];
  for (i = 0; i < loop_ub; i++) {
    C_data[i] = x_data[i];
  }
  /* 'genIndex:100' f = find(C(:,1)==1); */
  loop_ub = C->size[0];
  i = b_d->size[0];
  b_d->size[0] = C->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] = (C_data[i] == 1.0);
  }
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = iidx_data[i];
  }
  /* 'genIndex:101' S = length(f); */
  /* 'genIndex:102' ff = find(C(f,end)~=l); */
  i = b_d->size[0];
  b_d->size[0] = f->size[0];
  emxEnsureCapacity_boolean_T(b_d, i);
  b_e_data = b_d->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_e_data[i] =
        (C_data[(f_data[i] + C->size[0] * (C->size[1] - 1)) - 1] != l);
  }
  b_eml_find(b_d, iidx);
  iidx_data = iidx->data;
  i = d->size[0];
  d->size[0] = iidx->size[0];
  emxEnsureCapacity_real_T(d, i);
  d_data = d->data;
  loop_ub = iidx->size[0];
  for (i = 0; i < loop_ub; i++) {
    d_data[i] = iidx_data[i];
  }
  /* 'genIndex:103' for i = 1:length(ff)-1 */
  i = d->size[0];
  for (b_i = 0; b_i <= i - 2; b_i++) {
    /* 'genIndex:104' for j = i+1:length(ff) */
    i1 = d->size[0] - b_i;
    for (j = 0; j <= i1 - 2; j++) {
      b_j = ((unsigned int)b_i + j) + 2U;
      /* 'genIndex:105' if
       * sum(abs(C(f(ff(i)),:)+1-fliplr(l+1-C(f(ff(j)),:))))==0 */
      loop_ub = C->size[1];
      nd2 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = C->size[1];
      emxEnsureCapacity_real_T(x, nd2);
      x_data = x->data;
      for (nd2 = 0; nd2 < loop_ub; nd2++) {
        x_data[nd2] =
            (l + 1.0) -
            C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] + C->size[0] * nd2) -
                   1];
      }
      nx = x->size[1] - 1;
      nd2 = x->size[1] >> 1;
      for (b_j1 = 0; b_j1 < nd2; b_j1++) {
        j2 = nx - b_j1;
        xtmp = x_data[b_j1];
        x_data[b_j1] = x_data[j2];
        x_data[j2] = xtmp;
      }
      loop_ub = C->size[1];
      if (C->size[1] == x->size[1]) {
        nd2 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = C->size[1];
        emxEnsureCapacity_real_T(x, nd2);
        x_data = x->data;
        for (nd2 = 0; nd2 < loop_ub; nd2++) {
          x_data[nd2] =
              (C_data[(f_data[(int)d_data[b_i] - 1] + C->size[0] * nd2) - 1] +
               1.0) -
              x_data[nd2];
        }
      } else {
        h_binary_expand_op(x, C, f, d, b_i);
        x_data = x->data;
      }
      nx = x->size[1];
      nd2 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = x->size[1];
      emxEnsureCapacity_real_T(y, nd2);
      y_data = y->data;
      for (j2 = 0; j2 < nx; j2++) {
        y_data[j2] = fabs(x_data[j2]);
      }
      if (sum(y) == 0.0) {
        /* 'genIndex:106' C(f(ff(j)),:) = fliplr(l+1-C(f(ff(j)),:)); */
        loop_ub = C->size[1];
        nd2 = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = C->size[1];
        emxEnsureCapacity_real_T(x, nd2);
        x_data = x->data;
        for (nd2 = 0; nd2 < loop_ub; nd2++) {
          x_data[nd2] =
              (l + 1.0) - C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] +
                                  C->size[0] * nd2) -
                                 1];
        }
        nx = x->size[1] - 1;
        nd2 = x->size[1] >> 1;
        for (b_j1 = 0; b_j1 < nd2; b_j1++) {
          j2 = nx - b_j1;
          xtmp = x_data[b_j1];
          x_data[b_j1] = x_data[j2];
          x_data[j2] = xtmp;
        }
        loop_ub = x->size[1];
        for (nd2 = 0; nd2 < loop_ub; nd2++) {
          C_data[(f_data[(int)d_data[(int)b_j - 1] - 1] + C->size[0] * nd2) -
                 1] = x_data[nd2];
        }
      }
    }
  }
  /* 'genIndex:110' [~,ind2] = sort(sum(0.5.^C,2),'descend'); */
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = C->size[0];
  varargin_1->size[1] = C->size[1];
  emxEnsureCapacity_real_T(varargin_1, i);
  vec_data = varargin_1->data;
  loop_ub = C->size[0] * C->size[1];
  for (i = 0; i < loop_ub; i++) {
    a = C_data[i];
    vec_data[i] = pow(0.5, a);
  }
  b_sum(varargin_1, d);
  sort(d, iidx);
  iidx_data = iidx->data;
  i = f->size[0];
  f->size[0] = iidx->size[0];
  emxEnsureCapacity_int32_T(f, i);
  f_data = f->data;
  loop_ub = iidx->size[0];
  emxFree_real_T(&varargin_1);
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = iidx_data[i];
  }
  emxFree_int32_T(&iidx);
  /* 'genIndex:111' C = C(ind2,:); */
  nx = C->size[1] - 1;
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = f->size[0];
  varargin_2->size[1] = C->size[1];
  emxEnsureCapacity_real_T(varargin_2, i);
  x_data = varargin_2->data;
  for (i = 0; i <= nx; i++) {
    loop_ub = f->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      x_data[i1 + varargin_2->size[0] * i] =
          C_data[(f_data[i1] + C->size[0] * i) - 1];
    }
  }
  i = C->size[0] * C->size[1];
  C->size[0] = varargin_2->size[0];
  C->size[1] = varargin_2->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  loop_ub = varargin_2->size[0] * varargin_2->size[1];
  for (i = 0; i < loop_ub; i++) {
    C_data[i] = x_data[i];
  }
  /* 'genIndex:112' ind = ind(ind2); */
  i = ind->size[0];
  ind->size[0] = f->size[0];
  emxEnsureCapacity_real_T(ind, i);
  vec_data = ind->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    vec_data[i] = e_data[f_data[i] - 1];
  }
  /* 'genIndex:113' if n_frac == 1 */
  if (n_frac == 1.0) {
    /* 'genIndex:114' S = sum(C(:,1)==1); */
    loop_ub = C->size[0];
    i = b_d->size[0];
    b_d->size[0] = C->size[0];
    emxEnsureCapacity_boolean_T(b_d, i);
    b_e_data = b_d->data;
    for (i = 0; i < loop_ub; i++) {
      b_e_data[i] = (C_data[i] == 1.0);
    }
    nx = b_d->size[0];
    if (b_d->size[0] == 0) {
      nd2 = 0;
    } else {
      nd2 = b_e_data[0];
      for (j2 = 2; j2 <= nx; j2++) {
        nd2 += b_e_data[j2 - 1];
      }
    }
    /* 'genIndex:115' mat = zeros(S, max(c(:,1))-1); */
    loop_ub = c->size[0];
    i = d->size[0];
    d->size[0] = loop_ub;
    emxEnsureCapacity_real_T(d, i);
    d_data = d->data;
    for (i = 0; i < loop_ub; i++) {
      d_data[i] = c_data[i];
    }
    a = maximum(d);
    i = mat->size[0] * mat->size[1];
    mat->size[0] = nd2;
    mat->size[1] = (int)(a - 1.0);
    emxEnsureCapacity_real_T(mat, i);
    vec_data = mat->data;
    loop_ub = nd2 * (int)(a - 1.0);
    for (i = 0; i < loop_ub; i++) {
      vec_data[i] = 0.0;
    }
    /* 'genIndex:116' for i = 1:S */
    for (b_i = 0; b_i < nd2; b_i++) {
      /* 'genIndex:117' if C(i,end) ~= l */
      if (C_data[b_i + C->size[0] * (C->size[1] - 1)] != l) {
        /* 'genIndex:118' for j = 2:max(c(:,1)) */
        loop_ub = c->size[0];
        i = d->size[0];
        d->size[0] = loop_ub;
        emxEnsureCapacity_real_T(d, i);
        d_data = d->data;
        for (i = 0; i < loop_ub; i++) {
          d_data[i] = c_data[i];
        }
        i = (int)(maximum(d) + -1.0);
        for (j = 0; j < i; j++) {
          /* 'genIndex:119' a = S+1; */
          a = (double)nd2 + 1.0;
          /* 'genIndex:120' while sum(abs(C(i,:)+j-1-C(a,:))) ~= 0 */
          do {
            exitg1 = 0;
            loop_ub = C->size[1];
            i1 = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = C->size[1];
            emxEnsureCapacity_real_T(x, i1);
            x_data = x->data;
            for (i1 = 0; i1 < loop_ub; i1++) {
              x_data[i1] =
                  ((C_data[b_i + C->size[0] * i1] + ((double)j + 2.0)) - 1.0) -
                  C_data[((int)a + C->size[0] * i1) - 1];
            }
            nx = x->size[1];
            i1 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = x->size[1];
            emxEnsureCapacity_real_T(y, i1);
            y_data = y->data;
            for (j2 = 0; j2 < nx; j2++) {
              y_data[j2] = fabs(x_data[j2]);
            }
            if (sum(y) != 0.0) {
              /* 'genIndex:121' a = a + 1; */
              a++;
              /* 'genIndex:122' if a > L */
              if (a > L) {
                exitg1 = 1;
              }
            } else {
              exitg1 = 1;
            }
          } while (exitg1 == 0);
          /* 'genIndex:126' if a <= L */
          if (a <= L) {
            /* 'genIndex:127' mat(i,j-1) = a; */
            vec_data[b_i + mat->size[0] * j] = a;
          }
        }
      }
    }
    /* 'genIndex:132' mat = [(1:S)' mat]; */
    if (nd2 < 1) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = nd2;
      emxEnsureCapacity_real_T(y, i);
      y_data = y->data;
      loop_ub = nd2 - 1;
      for (i = 0; i <= loop_ub; i++) {
        y_data[i] = (double)i + 1.0;
      }
    }
    i = d->size[0];
    d->size[0] = y->size[1];
    emxEnsureCapacity_real_T(d, i);
    d_data = d->data;
    loop_ub = y->size[1];
    for (i = 0; i < loop_ub; i++) {
      d_data[i] = y_data[i];
    }
    if (d->size[0] != 0) {
      j2 = d->size[0];
    } else if ((mat->size[0] != 0) && (mat->size[1] != 0)) {
      j2 = mat->size[0];
    } else {
      j2 = 0;
      if (mat->size[0] > 0) {
        j2 = mat->size[0];
      }
    }
    empty_non_axis_sizes = (j2 == 0);
    if (empty_non_axis_sizes || (d->size[0] != 0)) {
      input_sizes_idx_1 = 1;
    } else {
      input_sizes_idx_1 = 0;
    }
    if (empty_non_axis_sizes || ((mat->size[0] != 0) && (mat->size[1] != 0))) {
      nx = mat->size[1];
    } else {
      nx = 0;
    }
    i = varargin_2->size[0] * varargin_2->size[1];
    varargin_2->size[0] = j2;
    varargin_2->size[1] = input_sizes_idx_1 + nx;
    emxEnsureCapacity_real_T(varargin_2, i);
    x_data = varargin_2->data;
    loop_ub = input_sizes_idx_1;
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 < j2; i1++) {
        x_data[i1] = d_data[i1];
      }
    }
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < j2; i1++) {
        x_data[i1 + varargin_2->size[0] * (i + input_sizes_idx_1)] =
            vec_data[i1 + j2 * i];
      }
    }
    i = mat->size[0] * mat->size[1];
    mat->size[0] = varargin_2->size[0];
    mat->size[1] = varargin_2->size[1];
    emxEnsureCapacity_real_T(mat, i);
    vec_data = mat->data;
    loop_ub = varargin_2->size[0] * varargin_2->size[1];
    for (i = 0; i < loop_ub; i++) {
      vec_data[i] = x_data[i];
    }
  } else {
    /* 'genIndex:133' else */
    /* 'genIndex:134' mat = (1:L)'; */
    if (L < 1.0) {
      y->size[0] = 1;
      y->size[1] = 0;
    } else {
      i = y->size[0] * y->size[1];
      y->size[0] = 1;
      loop_ub = (int)floor(L - 1.0);
      y->size[1] = loop_ub + 1;
      emxEnsureCapacity_real_T(y, i);
      y_data = y->data;
      for (i = 0; i <= loop_ub; i++) {
        y_data[i] = (double)i + 1.0;
      }
    }
    i = d->size[0];
    d->size[0] = y->size[1];
    emxEnsureCapacity_real_T(d, i);
    d_data = d->data;
    loop_ub = y->size[1];
    for (i = 0; i < loop_ub; i++) {
      d_data[i] = y_data[i];
    }
    i = mat->size[0] * mat->size[1];
    mat->size[0] = y->size[1];
    mat->size[1] = 1;
    emxEnsureCapacity_real_T(mat, i);
    vec_data = mat->data;
    loop_ub = y->size[1];
    for (i = 0; i < loop_ub; i++) {
      vec_data[i] = d_data[i];
    }
  }
  emxFree_boolean_T(&b_d);
  emxFree_real_T(&varargin_2);
  emxFree_real_T(&d);
  /* 'genIndex:136' I = zeros(L,1); */
  i = b_I->size[0];
  b_I->size[0] = (int)L;
  emxEnsureCapacity_real_T(b_I, i);
  vec_data = b_I->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    vec_data[i] = 0.0;
  }
  /* 'genIndex:137' I(1)=2; */
  vec_data[0] = 2.0;
  /* 'genIndex:138' for i=2:length(I) */
  for (b_i = 0; b_i <= loop_ub_tmp - 2; b_i++) {
    /* 'genIndex:139' a = 1; */
    /* 'genIndex:140' while C(i-1,a)==C(i,a) */
    for (a = 1.0; C_data[b_i + C->size[0] * ((int)a - 1)] ==
                  C_data[(b_i + C->size[0] * ((int)a - 1)) + 1];
         a++) {
      /* 'genIndex:141' a = a+1; */
    }
    /* 'genIndex:143' I(i) = a; */
    vec_data[b_i + 1] = a;
    /* 'genIndex:144' if I(i) < 2 */
    if (a < 2.0) {
      /* 'genIndex:145' I(i)=2; */
      vec_data[b_i + 1] = 2.0;
    }
  }
  /* 'genIndex:148' for i = 1:L */
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    /* 'genIndex:149' if sum(abs(c(ind(i),:)-C(i,:)))~=0 */
    if (c->size[1] == C->size[1]) {
      loop_ub = c->size[1];
      i = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = loop_ub;
      emxEnsureCapacity_real_T(x, i);
      x_data = x->data;
      for (i = 0; i < loop_ub; i++) {
        x_data[i] = c_data[(e_data[f_data[b_i] - 1] + c->size[0] * i) - 1] -
                    C_data[b_i + C->size[0] * i];
      }
    } else {
      i_binary_expand_op(x, c, b_ind, f, b_i, C);
      x_data = x->data;
    }
    nx = x->size[1];
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real_T(y, i);
    y_data = y->data;
    for (j2 = 0; j2 < nx; j2++) {
      y_data[j2] = fabs(x_data[j2]);
    }
    if (sum(y) != 0.0) {
      /* 'genIndex:150' c(ind(i),:) = l+1-fliplr(c(ind(i),:)); */
      loop_ub = c->size[1];
      i = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = loop_ub;
      emxEnsureCapacity_real_T(x, i);
      x_data = x->data;
      for (i = 0; i < loop_ub; i++) {
        x_data[i] = c_data[(e_data[f_data[b_i] - 1] + c->size[0] * i) - 1];
      }
      nd2 = c->size[1] >> 1;
      for (b_j1 = 0; b_j1 < nd2; b_j1++) {
        j2 = (c->size[1] - b_j1) - 1;
        xtmp = x_data[b_j1];
        x_data[b_j1] = x_data[j2];
        x_data[j2] = xtmp;
      }
      loop_ub = x->size[1];
      for (i = 0; i < loop_ub; i++) {
        c_data[(e_data[f_data[b_i] - 1] + c->size[0] * i) - 1] =
            (l + 1.0) - x_data[i];
      }
    }
  }
  emxFree_real_T(&x);
  emxFree_real_T(&y);
  emxFree_int32_T(&b_ind);
  emxFree_int32_T(&f);
}

/*
 * function [c,C,I,ind,mat,rcnum] = genIndex(l,k,n_frac)
 */
void genIndex(double l, double k, double n_frac, emxArray_real_T *c,
              emxArray_real_T *C, emxArray_real_T *b_I, emxArray_real_T *ind,
              emxArray_real_T *mat, double *rcnum)
{
  emxArray_boolean_T *c_c1;
  emxArray_int32_T *b_r;
  emxArray_int32_T *r4;
  emxArray_real_T *a__2;
  emxArray_real_T *a__3;
  emxArray_real_T *a__4;
  emxArray_real_T *b_c1;
  emxArray_real_T *b_c2;
  emxArray_real_T *c1;
  emxArray_real_T *c2;
  emxArray_real_T *cc1;
  emxArray_real_T *f;
  emxArray_real_T *r1;
  double L1;
  double xtmp;
  double *c1_data;
  double *c2_data;
  double *cc1_data;
  double *f_data;
  double *r3;
  int b_i;
  int i;
  int i1;
  int j2;
  int nd2;
  unsigned int r;
  int sizes_idx_0;
  int *r2;
  bool empty_non_axis_sizes;
  bool *b_c1_data;
  /* l = k-mer length */
  /* k = # of ungapped positions */
  /* 'genIndex:4' if n_frac > 1  || n_frac < 0 */
  if ((n_frac > 1.0) || (n_frac < 0.0)) {
    /* 'genIndex:5' fprintf('ERROR: n_frac must be a fraction in [0,1]\n'); */
    printf("ERROR: n_frac must be a fraction in [0,1]\n");
    fflush(stdout);
    exit(1);
  }
  /* 'genIndex:7' if n_frac ~= 1 */
  emxInit_real_T(&f, 2);
  if (n_frac != 1.0) {
    /* 'genIndex:8' rng(1); */
    r = 1U;
    state[0] = 1U;
    for (nd2 = 0; nd2 < 623; nd2++) {
      r = ((r ^ r >> 30U) * 1812433253U + nd2) + 1U;
      state[nd2 + 1] = r;
    }
    state[624] = 624U;
    /* 'genIndex:9' if mod(l,2) == 0 */
    if (l == 0.0) {
      xtmp = 0.0;
    } else {
      xtmp = fmod(l, 2.0);
      if (xtmp == 0.0) {
        xtmp = 0.0;
      } else if (l < 0.0) {
        xtmp += 2.0;
      }
    }
    emxInit_real_T(&c2, 2);
    if (xtmp == 0.0) {
      /*  c = combnk(1:l,k); */
      /* 'genIndex:11' c = flip(nchoosek(1:l,k)); */
      if (l < 1.0) {
        f->size[0] = 1;
        f->size[1] = 0;
      } else {
        i = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = (int)floor(l - 1.0) + 1;
        emxEnsureCapacity_real_T(f, i);
        f_data = f->data;
        nd2 = (int)floor(l - 1.0);
        for (i = 0; i <= nd2; i++) {
          f_data[i] = (double)i + 1.0;
        }
      }
      nchoosek(f, k, c2);
      flip(c2);
      /* 'genIndex:12' L = numel(c)/k; */
      /* 'genIndex:13' [c,C,I,ind,mat,rcnum] = sort_comb(c,l,k,n_frac); */
      sort_comb(c2, l, k, n_frac, C, b_I, ind, mat, rcnum);
      /* 'genIndex:14' L = numel(c)/k; */
      /* 'genIndex:15' c = genIndex_frac(l,k,c,rcnum,L*n_frac); */
      genIndex_frac(l, k, c2, (double)(c2->size[0] * c2->size[1]) / k * n_frac,
                    c);
      /* 'genIndex:16' [c,C,I,ind,mat,rcnum] = sort_comb(c,l,k,n_frac); */
      sort_comb(c, l, k, n_frac, C, b_I, ind, mat, rcnum);
    } else {
      /* 'genIndex:17' else */
      /*  c = combnk(1:l,k); */
      /* 'genIndex:19' c = flip(nchoosek(1:l,k)); */
      if (l < 1.0) {
        f->size[0] = 1;
        f->size[1] = 0;
      } else {
        i = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = (int)floor(l - 1.0) + 1;
        emxEnsureCapacity_real_T(f, i);
        f_data = f->data;
        nd2 = (int)floor(l - 1.0);
        for (i = 0; i <= nd2; i++) {
          f_data[i] = (double)i + 1.0;
        }
      }
      emxInit_real_T(&a__2, 1);
      emxInit_real_T(&a__3, 1);
      emxInit_real_T(&a__4, 2);
      nchoosek(f, k, c);
      flip(c);
      /* 'genIndex:20' [c,~,~,~,~,rcnum1] = sort_comb(c,l,k,n_frac); */
      sort_comb(c, l, k, n_frac, c2, a__2, a__3, a__4, &xtmp);
      /*  c1 = combnk(1:(l-1),k-1); */
      /* 'genIndex:22' c1 = flip(nchoosek(1:(l-1),(k-1))); */
      if (l - 1.0 < 1.0) {
        f->size[0] = 1;
        f->size[1] = 0;
      } else {
        i = f->size[0] * f->size[1];
        f->size[0] = 1;
        nd2 = (int)floor((l - 1.0) - 1.0);
        f->size[1] = nd2 + 1;
        emxEnsureCapacity_real_T(f, i);
        f_data = f->data;
        for (i = 0; i <= nd2; i++) {
          f_data[i] = (double)i + 1.0;
        }
      }
      emxInit_real_T(&c1, 2);
      emxInit_real_T(&b_c1, 2);
      nchoosek(f, k - 1.0, b_c1);
      flip(b_c1);
      /* 'genIndex:23' [c1,~,~,~,~,rcnum1] = sort_comb(c1,l-1,k-1,n_frac); */
      sort_comb(b_c1, l - 1.0, k - 1.0, n_frac, c2, a__2, a__3, a__4, &xtmp);
      /* 'genIndex:24' L1 = numel(c1)/(k-1); */
      /* 'genIndex:25' c1 = genIndex_frac(l-1,k-1,c1,rcnum1,L1*n_frac); */
      genIndex_frac(
          l - 1.0, k - 1.0, b_c1,
          (double)(b_c1->size[0] * b_c1->size[1]) / (k - 1.0) * n_frac, c1);
      c1_data = c1->data;
      /* 'genIndex:26' L1 = numel(c1)/(k-1); */
      L1 = (double)(c1->size[0] * c1->size[1]) / (k - 1.0);
      /* 'genIndex:27' x = zeros(L1, l-1); */
      i = (int)L1;
      i1 = c2->size[0] * c2->size[1];
      c2->size[0] = (int)L1;
      c2->size[1] = (int)(l - 1.0);
      emxEnsureCapacity_real_T(c2, i1);
      c2_data = c2->data;
      nd2 = (int)L1 * (int)(l - 1.0);
      for (i1 = 0; i1 < nd2; i1++) {
        c2_data[i1] = 0.0;
      }
      /* 'genIndex:28' for i = 1:L1 */
      emxInit_int32_T(&b_r, 1);
      emxInit_real_T(&r1, 2);
      for (b_i = 0; b_i < i; b_i++) {
        /* 'genIndex:29' x(i,c1(i,:)) = x(i,c1(i,:))+1; */
        nd2 = c1->size[1];
        i1 = b_r->size[0];
        b_r->size[0] = c1->size[1];
        emxEnsureCapacity_int32_T(b_r, i1);
        r2 = b_r->data;
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = 1;
        r1->size[1] = c1->size[1];
        emxEnsureCapacity_real_T(r1, i1);
        r3 = r1->data;
        for (i1 = 0; i1 < nd2; i1++) {
          xtmp = c1_data[b_i + c1->size[0] * i1];
          r2[i1] = (int)xtmp - 1;
          r3[i1] = c2_data[b_i + c2->size[0] * ((int)xtmp - 1)];
        }
        nd2 = r1->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          c2_data[b_i + c2->size[0] * r2[i1]] = r3[i1] + 1.0;
        }
        /* 'genIndex:30' x(i,l-fliplr(c1(i,:))) = x(i,l-fliplr(c1(i,:)))+1; */
        nd2 = c1->size[1];
        i1 = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = c1->size[1];
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        for (i1 = 0; i1 < nd2; i1++) {
          f_data[i1] = c1_data[b_i + c1->size[0] * i1];
        }
        nd2 = c1->size[1] >> 1;
        for (sizes_idx_0 = 0; sizes_idx_0 < nd2; sizes_idx_0++) {
          j2 = (c1->size[1] - sizes_idx_0) - 1;
          xtmp = f_data[sizes_idx_0];
          f_data[sizes_idx_0] = f_data[j2];
          f_data[j2] = xtmp;
        }
        i1 = b_r->size[0];
        b_r->size[0] = f->size[1];
        emxEnsureCapacity_int32_T(b_r, i1);
        r2 = b_r->data;
        nd2 = f->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          r2[i1] = (int)(l - f_data[i1]) - 1;
        }
        nd2 = c1->size[1];
        i1 = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = c1->size[1];
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        for (i1 = 0; i1 < nd2; i1++) {
          f_data[i1] = c1_data[b_i + c1->size[0] * i1];
        }
        nd2 = c1->size[1] >> 1;
        for (sizes_idx_0 = 0; sizes_idx_0 < nd2; sizes_idx_0++) {
          j2 = (c1->size[1] - sizes_idx_0) - 1;
          xtmp = f_data[sizes_idx_0];
          f_data[sizes_idx_0] = f_data[j2];
          f_data[j2] = xtmp;
        }
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = 1;
        r1->size[1] = f->size[1];
        emxEnsureCapacity_real_T(r1, i1);
        r3 = r1->data;
        nd2 = f->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          r3[i1] = c2_data[b_i + c2->size[0] * ((int)(l - f_data[i1]) - 1)];
        }
        nd2 = r1->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          c2_data[b_i + c2->size[0] * r2[i1]] = r3[i1] + 1.0;
        }
      }
      /* 'genIndex:32' L2 = round((L1*(l-1)*2-sum(sum(x)))/k/2); */
      /*  c2 = combnk(1:(l-1),k); */
      /* 'genIndex:34' c2 = flip(nchoosek(1:(l-1),k)); */
      if (l - 1.0 < 1.0) {
        f->size[0] = 1;
        f->size[1] = 0;
      } else {
        i1 = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = (int)floor((l - 1.0) - 1.0) + 1;
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        nd2 = (int)floor((l - 1.0) - 1.0);
        for (i1 = 0; i1 <= nd2; i1++) {
          f_data[i1] = (double)i1 + 1.0;
        }
      }
      emxInit_real_T(&cc1, 2);
      emxInit_real_T(&b_c2, 2);
      nchoosek(f, k, b_c2);
      flip(b_c2);
      /* 'genIndex:35' [c2,~,~,~,~,rcnum2] = sort_comb(c2,l-1,k,n_frac); */
      sort_comb(b_c2, l - 1.0, k, n_frac, a__4, a__2, a__3, b_c1, &xtmp);
      /* 'genIndex:36' c2 = genIndex_frac(l-1,k,c2,rcnum2,L2); */
      c_sum(c2, r1);
      genIndex_frac(l - 1.0, k, b_c2,
                    rt_roundd((L1 * (l - 1.0) * 2.0 - sum(r1)) / k / 2.0), c2);
      c2_data = c2->data;
      /* 'genIndex:37' L2 = numel(c2)/k; */
      /* 'genIndex:38' cc1 = zeros(L1,k); */
      i1 = cc1->size[0] * cc1->size[1];
      cc1->size[0] = (int)L1;
      cc1->size[1] = (int)k;
      emxEnsureCapacity_real_T(cc1, i1);
      cc1_data = cc1->data;
      /* 'genIndex:39' M = ceil(l/2); */
      xtmp = ceil(l / 2.0);
      /* 'genIndex:40' for i = 1:L1 */
      emxFree_real_T(&b_c2);
      emxFree_real_T(&b_c1);
      emxFree_real_T(&a__4);
      emxFree_real_T(&a__3);
      emxInit_int32_T(&r4, 2);
      emxInit_boolean_T(&c_c1, 2);
      for (b_i = 0; b_i < i; b_i++) {
        /* 'genIndex:41' f = find(c1(i,:)>=M); */
        nd2 = c1->size[1];
        i1 = c_c1->size[0] * c_c1->size[1];
        c_c1->size[0] = 1;
        c_c1->size[1] = c1->size[1];
        emxEnsureCapacity_boolean_T(c_c1, i1);
        b_c1_data = c_c1->data;
        for (i1 = 0; i1 < nd2; i1++) {
          b_c1_data[i1] = (c1_data[b_i + c1->size[0] * i1] >= xtmp);
        }
        eml_find(c_c1, r4);
        r2 = r4->data;
        i1 = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = r4->size[1];
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        nd2 = r4->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          f_data[i1] = r2[i1];
        }
        /* 'genIndex:42' c1(i,f) = c1(i,f)+1; */
        i1 = b_r->size[0];
        b_r->size[0] = f->size[1];
        emxEnsureCapacity_int32_T(b_r, i1);
        r2 = b_r->data;
        nd2 = f->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          r2[i1] = (int)f_data[i1] - 1;
        }
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = 1;
        r1->size[1] = f->size[1];
        emxEnsureCapacity_real_T(r1, i1);
        r3 = r1->data;
        nd2 = f->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          r3[i1] = c1_data[b_i + c1->size[0] * ((int)f_data[i1] - 1)];
        }
        nd2 = r1->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          c1_data[b_i + c1->size[0] * r2[i1]] = r3[i1] + 1.0;
        }
        /* 'genIndex:43' cc1(i,:) = sort([c1(i,:) M]); */
        nd2 = c1->size[1];
        i1 = r1->size[0] * r1->size[1];
        r1->size[0] = 1;
        r1->size[1] = c1->size[1] + 1;
        emxEnsureCapacity_real_T(r1, i1);
        r3 = r1->data;
        for (i1 = 0; i1 < nd2; i1++) {
          r3[i1] = c1_data[b_i + c1->size[0] * i1];
        }
        r3[c1->size[1]] = xtmp;
        b_sort(r1);
        r3 = r1->data;
        nd2 = r1->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          cc1_data[b_i + cc1->size[0] * i1] = r3[i1];
        }
      }
      emxFree_real_T(&r1);
      emxFree_int32_T(&b_r);
      emxFree_real_T(&c1);
      /* 'genIndex:45' for i = 1:L2 */
      i = (int)((double)(c2->size[0] * c2->size[1]) / k);
      for (b_i = 0; b_i < i; b_i++) {
        /* 'genIndex:46' f = find(c2(i,:)>=M); */
        /* 'genIndex:47' c2(i,f) = c2(i,f)+1; */
        nd2 = c2->size[1];
        i1 = c_c1->size[0] * c_c1->size[1];
        c_c1->size[0] = 1;
        c_c1->size[1] = c2->size[1];
        emxEnsureCapacity_boolean_T(c_c1, i1);
        b_c1_data = c_c1->data;
        for (i1 = 0; i1 < nd2; i1++) {
          b_c1_data[i1] = (c2_data[b_i + c2->size[0] * i1] >= xtmp);
        }
        eml_find(c_c1, r4);
        r2 = r4->data;
        i1 = a__2->size[0];
        a__2->size[0] = r4->size[1];
        emxEnsureCapacity_real_T(a__2, i1);
        c1_data = a__2->data;
        nd2 = r4->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          c1_data[i1] = r2[i1];
        }
        i1 = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = a__2->size[0];
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        nd2 = a__2->size[0];
        for (i1 = 0; i1 < nd2; i1++) {
          f_data[i1] =
              c2_data[b_i + c2->size[0] * ((int)c1_data[i1] - 1)] + 1.0;
        }
        nd2 = f->size[1];
        for (i1 = 0; i1 < nd2; i1++) {
          c2_data[b_i + c2->size[0] * ((int)c1_data[i1] - 1)] = f_data[i1];
        }
      }
      emxFree_boolean_T(&c_c1);
      emxFree_int32_T(&r4);
      emxFree_real_T(&a__2);
      /* 'genIndex:49' [c,C,I,ind,mat,rcnum] = sort_comb([cc1;c2],l,k,n_frac);
       */
      if ((cc1->size[0] != 0) && (cc1->size[1] != 0)) {
        nd2 = cc1->size[1];
      } else if ((c2->size[0] != 0) && (c2->size[1] != 0)) {
        nd2 = c2->size[1];
      } else {
        nd2 = cc1->size[1];
        if (c2->size[1] > cc1->size[1]) {
          nd2 = c2->size[1];
        }
      }
      empty_non_axis_sizes = (nd2 == 0);
      if (empty_non_axis_sizes ||
          ((cc1->size[0] != 0) && (cc1->size[1] != 0))) {
        j2 = cc1->size[0];
      } else {
        j2 = 0;
      }
      if (empty_non_axis_sizes || ((c2->size[0] != 0) && (c2->size[1] != 0))) {
        sizes_idx_0 = c2->size[0];
      } else {
        sizes_idx_0 = 0;
      }
      i = c->size[0] * c->size[1];
      c->size[0] = j2 + sizes_idx_0;
      c->size[1] = nd2;
      emxEnsureCapacity_real_T(c, i);
      c1_data = c->data;
      for (i = 0; i < nd2; i++) {
        for (i1 = 0; i1 < j2; i1++) {
          c1_data[i1 + c->size[0] * i] = cc1_data[i1 + j2 * i];
        }
      }
      emxFree_real_T(&cc1);
      for (i = 0; i < nd2; i++) {
        for (i1 = 0; i1 < sizes_idx_0; i1++) {
          c1_data[(i1 + j2) + c->size[0] * i] = c2_data[i1 + sizes_idx_0 * i];
        }
      }
      sort_comb(c, l, k, n_frac, C, b_I, ind, mat, rcnum);
    }
    emxFree_real_T(&c2);
  } else {
    /* 'genIndex:51' else */
    /*  c = combnk(1:l,k); */
    /* 'genIndex:53' c = flip(nchoosek(1:l,k)); */
    if (l < 1.0) {
      f->size[0] = 1;
      f->size[1] = 0;
    } else {
      i = f->size[0] * f->size[1];
      f->size[0] = 1;
      nd2 = (int)floor(l - 1.0);
      f->size[1] = nd2 + 1;
      emxEnsureCapacity_real_T(f, i);
      f_data = f->data;
      for (i = 0; i <= nd2; i++) {
        f_data[i] = (double)i + 1.0;
      }
    }
    nchoosek(f, k, c);
    flip(c);
    /* 'genIndex:54' [c,C,I,ind,mat,rcnum] = sort_comb(c,l,k,n_frac); */
    sort_comb(c, l, k, 1.0, C, b_I, ind, mat, rcnum);
  }
  emxFree_real_T(&f);
}

/* End of code generation (genIndex.c) */
