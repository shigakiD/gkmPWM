/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sort.c
 *
 * Code generation for function 'sort'
 *
 */

/* Include files */
#include "sort.h"
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include "sortIdx.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  b_sortIdx(x, idx);
}

/*
 *
 */
void c_sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_int32_T *iidx;
  emxArray_int32_T *iwork;
  emxArray_real_T *vwork;
  emxArray_real_T *xwork;
  double b_xwork[256];
  double x4[4];
  double d;
  double d1;
  double *vwork_data;
  double *x_data;
  double *xwork_data;
  int b_iwork[256];
  int idx4[4];
  int b;
  int b_b;
  int b_n;
  int dim;
  int exitg1;
  int i1;
  int i2;
  int i3;
  int i4;
  int iidx_tmp;
  int j;
  int k;
  int n;
  int nBlocks;
  int nNaNs;
  int nNonNaN;
  int nPairs;
  int vlen;
  int vstride;
  int *idx_data;
  int *iidx_data;
  int *iwork_data;
  signed char perm[4];
  x_data = x->data;
  dim = 0;
  if (x->size[0] != 1) {
    dim = -1;
  }
  emxInit_real_T(&vwork, 1);
  if (dim + 2 <= 1) {
    i2 = x->size[0];
  } else {
    i2 = 1;
  }
  vlen = i2 - 1;
  i1 = vwork->size[0];
  vwork->size[0] = i2;
  emxEnsureCapacity_real_T(vwork, i1);
  vwork_data = vwork->data;
  i2 = idx->size[0];
  idx->size[0] = x->size[0];
  emxEnsureCapacity_int32_T(idx, i2);
  idx_data = idx->data;
  vstride = 1;
  for (k = 0; k <= dim; k++) {
    vstride *= x->size[0];
  }
  emxInit_int32_T(&iidx, 1);
  emxInit_int32_T(&iwork, 1);
  emxInit_real_T(&xwork, 1);
  for (j = 0; j < vstride; j++) {
    for (k = 0; k <= vlen; k++) {
      vwork_data[k] = x_data[j + k * vstride];
    }
    i2 = iidx->size[0];
    iidx->size[0] = vwork->size[0];
    emxEnsureCapacity_int32_T(iidx, i2);
    iidx_data = iidx->data;
    i1 = vwork->size[0];
    for (i2 = 0; i2 < i1; i2++) {
      iidx_data[i2] = 0;
    }
    if (vwork->size[0] != 0) {
      n = vwork->size[0];
      b_n = vwork->size[0];
      x4[0] = 0.0;
      idx4[0] = 0;
      x4[1] = 0.0;
      idx4[1] = 0;
      x4[2] = 0.0;
      idx4[2] = 0;
      x4[3] = 0.0;
      idx4[3] = 0;
      i2 = iwork->size[0];
      iwork->size[0] = vwork->size[0];
      emxEnsureCapacity_int32_T(iwork, i2);
      iwork_data = iwork->data;
      i1 = vwork->size[0];
      for (i2 = 0; i2 < i1; i2++) {
        iwork_data[i2] = 0;
      }
      i2 = xwork->size[0];
      xwork->size[0] = vwork->size[0];
      emxEnsureCapacity_real_T(xwork, i2);
      xwork_data = xwork->data;
      i1 = vwork->size[0];
      for (i2 = 0; i2 < i1; i2++) {
        xwork_data[i2] = 0.0;
      }
      nNaNs = 0;
      dim = -1;
      for (k = 0; k < b_n; k++) {
        if (rtIsNaN(vwork_data[k])) {
          iidx_tmp = (b_n - nNaNs) - 1;
          iidx_data[iidx_tmp] = k + 1;
          xwork_data[iidx_tmp] = vwork_data[k];
          nNaNs++;
        } else {
          dim++;
          idx4[dim] = k + 1;
          x4[dim] = vwork_data[k];
          if (dim + 1 == 4) {
            dim = k - nNaNs;
            if (x4[0] >= x4[1]) {
              i1 = 1;
              i2 = 2;
            } else {
              i1 = 2;
              i2 = 1;
            }
            if (x4[2] >= x4[3]) {
              i3 = 3;
              i4 = 4;
            } else {
              i3 = 4;
              i4 = 3;
            }
            d = x4[i1 - 1];
            d1 = x4[i3 - 1];
            if (d >= d1) {
              d = x4[i2 - 1];
              if (d >= d1) {
                perm[0] = (signed char)i1;
                perm[1] = (signed char)i2;
                perm[2] = (signed char)i3;
                perm[3] = (signed char)i4;
              } else if (d >= x4[i4 - 1]) {
                perm[0] = (signed char)i1;
                perm[1] = (signed char)i3;
                perm[2] = (signed char)i2;
                perm[3] = (signed char)i4;
              } else {
                perm[0] = (signed char)i1;
                perm[1] = (signed char)i3;
                perm[2] = (signed char)i4;
                perm[3] = (signed char)i2;
              }
            } else {
              d1 = x4[i4 - 1];
              if (d >= d1) {
                if (x4[i2 - 1] >= d1) {
                  perm[0] = (signed char)i3;
                  perm[1] = (signed char)i1;
                  perm[2] = (signed char)i2;
                  perm[3] = (signed char)i4;
                } else {
                  perm[0] = (signed char)i3;
                  perm[1] = (signed char)i1;
                  perm[2] = (signed char)i4;
                  perm[3] = (signed char)i2;
                }
              } else {
                perm[0] = (signed char)i3;
                perm[1] = (signed char)i4;
                perm[2] = (signed char)i1;
                perm[3] = (signed char)i2;
              }
            }
            iidx_data[dim - 3] = idx4[perm[0] - 1];
            iidx_data[dim - 2] = idx4[perm[1] - 1];
            iidx_data[dim - 1] = idx4[perm[2] - 1];
            iidx_data[dim] = idx4[perm[3] - 1];
            vwork_data[dim - 3] = x4[perm[0] - 1];
            vwork_data[dim - 2] = x4[perm[1] - 1];
            vwork_data[dim - 1] = x4[perm[2] - 1];
            vwork_data[dim] = x4[perm[3] - 1];
            dim = -1;
          }
        }
      }
      i3 = (b_n - nNaNs) - 1;
      if (dim + 1 > 0) {
        perm[1] = 0;
        perm[2] = 0;
        perm[3] = 0;
        if (dim + 1 == 1) {
          perm[0] = 1;
        } else if (dim + 1 == 2) {
          if (x4[0] >= x4[1]) {
            perm[0] = 1;
            perm[1] = 2;
          } else {
            perm[0] = 2;
            perm[1] = 1;
          }
        } else if (x4[0] >= x4[1]) {
          if (x4[1] >= x4[2]) {
            perm[0] = 1;
            perm[1] = 2;
            perm[2] = 3;
          } else if (x4[0] >= x4[2]) {
            perm[0] = 1;
            perm[1] = 3;
            perm[2] = 2;
          } else {
            perm[0] = 3;
            perm[1] = 1;
            perm[2] = 2;
          }
        } else if (x4[0] >= x4[2]) {
          perm[0] = 2;
          perm[1] = 1;
          perm[2] = 3;
        } else if (x4[1] >= x4[2]) {
          perm[0] = 2;
          perm[1] = 3;
          perm[2] = 1;
        } else {
          perm[0] = 3;
          perm[1] = 2;
          perm[2] = 1;
        }
        for (k = 0; k <= dim; k++) {
          iidx_tmp = perm[k] - 1;
          i1 = (i3 - dim) + k;
          iidx_data[i1] = idx4[iidx_tmp];
          vwork_data[i1] = x4[iidx_tmp];
        }
      }
      dim = (nNaNs >> 1) + 1;
      for (k = 0; k <= dim - 2; k++) {
        i1 = (i3 + k) + 1;
        i2 = iidx_data[i1];
        iidx_tmp = (b_n - k) - 1;
        iidx_data[i1] = iidx_data[iidx_tmp];
        iidx_data[iidx_tmp] = i2;
        vwork_data[i1] = xwork_data[iidx_tmp];
        vwork_data[iidx_tmp] = xwork_data[i1];
      }
      if ((nNaNs & 1) != 0) {
        dim += i3;
        vwork_data[dim] = xwork_data[dim];
      }
      nNonNaN = n - nNaNs;
      dim = 2;
      if (nNonNaN > 1) {
        if (n >= 256) {
          nBlocks = nNonNaN >> 8;
          if (nBlocks > 0) {
            for (b = 0; b < nBlocks; b++) {
              i4 = (b << 8) - 1;
              for (b_b = 0; b_b < 6; b_b++) {
                b_n = 1 << (b_b + 2);
                n = b_n << 1;
                nPairs = 256 >> (b_b + 3);
                for (k = 0; k < nPairs; k++) {
                  i2 = (i4 + k * n) + 1;
                  for (i1 = 0; i1 < n; i1++) {
                    dim = i2 + i1;
                    b_iwork[i1] = iidx_data[dim];
                    b_xwork[i1] = vwork_data[dim];
                  }
                  i3 = 0;
                  i1 = b_n;
                  dim = i2 - 1;
                  do {
                    exitg1 = 0;
                    dim++;
                    if (b_xwork[i3] >= b_xwork[i1]) {
                      iidx_data[dim] = b_iwork[i3];
                      vwork_data[dim] = b_xwork[i3];
                      if (i3 + 1 < b_n) {
                        i3++;
                      } else {
                        exitg1 = 1;
                      }
                    } else {
                      iidx_data[dim] = b_iwork[i1];
                      vwork_data[dim] = b_xwork[i1];
                      if (i1 + 1 < n) {
                        i1++;
                      } else {
                        dim -= i3;
                        for (i1 = i3 + 1; i1 <= b_n; i1++) {
                          iidx_tmp = dim + i1;
                          iidx_data[iidx_tmp] = b_iwork[i1 - 1];
                          vwork_data[iidx_tmp] = b_xwork[i1 - 1];
                        }
                        exitg1 = 1;
                      }
                    }
                  } while (exitg1 == 0);
                }
              }
            }
            dim = nBlocks << 8;
            i1 = nNonNaN - dim;
            if (i1 > 0) {
              b_merge_block(iidx, vwork, dim, i1, 2, iwork, xwork);
            }
            dim = 8;
          }
        }
        b_merge_block(iidx, vwork, 0, nNonNaN, dim, iwork, xwork);
        xwork_data = xwork->data;
        iwork_data = iwork->data;
        vwork_data = vwork->data;
        iidx_data = iidx->data;
      }
      if ((nNaNs > 0) && (nNonNaN > 0)) {
        for (k = 0; k < nNaNs; k++) {
          dim = nNonNaN + k;
          xwork_data[k] = vwork_data[dim];
          iwork_data[k] = iidx_data[dim];
        }
        for (k = nNonNaN; k >= 1; k--) {
          dim = (nNaNs + k) - 1;
          vwork_data[dim] = vwork_data[k - 1];
          iidx_data[dim] = iidx_data[k - 1];
        }
        for (k = 0; k < nNaNs; k++) {
          vwork_data[k] = xwork_data[k];
          iidx_data[k] = iwork_data[k];
        }
      }
    }
    for (k = 0; k <= vlen; k++) {
      i2 = j + k * vstride;
      x_data[i2] = vwork_data[k];
      idx_data[i2] = iidx_data[k];
    }
  }
  emxFree_real_T(&xwork);
  emxFree_int32_T(&iwork);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

/*
 *
 */
void d_sort(emxArray_real_T *x)
{
  emxArray_int32_T *b_x;
  emxInit_int32_T(&b_x, 2);
  b_sortIdx(x, b_x);
  emxFree_int32_T(&b_x);
}

/*
 *
 */
void e_sort(emxArray_real_T *x)
{
  emxArray_int32_T *b_vwork;
  emxArray_real_T *vwork;
  double *vwork_data;
  double *x_data;
  int dim;
  int j;
  int k;
  int vlen;
  int vstride;
  x_data = x->data;
  dim = 0;
  if (x->size[0] != 1) {
    dim = -1;
  }
  emxInit_real_T(&vwork, 1);
  if (dim + 2 <= 1) {
    vstride = x->size[0];
  } else {
    vstride = 1;
  }
  vlen = vstride - 1;
  j = vwork->size[0];
  vwork->size[0] = vstride;
  emxEnsureCapacity_real_T(vwork, j);
  vwork_data = vwork->data;
  vstride = 1;
  for (k = 0; k <= dim; k++) {
    vstride *= x->size[0];
  }
  emxInit_int32_T(&b_vwork, 1);
  for (j = 0; j < vstride; j++) {
    for (k = 0; k <= vlen; k++) {
      vwork_data[k] = x_data[j + k * vstride];
    }
    sortIdx(vwork, b_vwork);
    vwork_data = vwork->data;
    for (k = 0; k <= vlen; k++) {
      x_data[j + k * vstride] = vwork_data[k];
    }
  }
  emxFree_int32_T(&b_vwork);
  emxFree_real_T(&vwork);
}

/*
 *
 */
void sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_int32_T *iidx;
  emxArray_real_T *vwork;
  double *vwork_data;
  double *x_data;
  int dim;
  int i;
  int k;
  int vlen;
  int vstride;
  int *idx_data;
  int *iidx_data;
  x_data = x->data;
  dim = 0;
  if (x->size[0] != 1) {
    dim = -1;
  }
  emxInit_real_T(&vwork, 1);
  if (dim + 2 <= 1) {
    i = x->size[0];
  } else {
    i = 1;
  }
  vlen = i - 1;
  vstride = vwork->size[0];
  vwork->size[0] = i;
  emxEnsureCapacity_real_T(vwork, vstride);
  vwork_data = vwork->data;
  i = idx->size[0];
  idx->size[0] = x->size[0];
  emxEnsureCapacity_int32_T(idx, i);
  idx_data = idx->data;
  vstride = 1;
  for (k = 0; k <= dim; k++) {
    vstride *= x->size[0];
  }
  emxInit_int32_T(&iidx, 1);
  for (dim = 0; dim < vstride; dim++) {
    for (k = 0; k <= vlen; k++) {
      vwork_data[k] = x_data[dim + k * vstride];
    }
    sortIdx(vwork, iidx);
    iidx_data = iidx->data;
    vwork_data = vwork->data;
    for (k = 0; k <= vlen; k++) {
      i = dim + k * vstride;
      x_data[i] = vwork_data[k];
      idx_data[i] = iidx_data[k];
    }
  }
  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

/* End of code generation (sort.c) */
