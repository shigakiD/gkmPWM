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
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include "sortIdx.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_sort(emxArray_real_T *x)
{
  emxArray_int32_T *b_x;
  emxInit_int32_T(&b_x, 2);
  sortIdx(x, b_x);
  emxFree_int32_T(&b_x);
}

/*
 *
 */
void c_sort(emxArray_real_T *x)
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
    b_sortIdx(vwork, b_vwork);
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
void d_sort(emxArray_real_T *x, emxArray_int32_T *idx)
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
    b_sortIdx(vwork, iidx);
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

/*
 *
 */
void e_sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  sortIdx(x, idx);
}

/*
 *
 */
void sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_int32_T *b_iwork;
  emxArray_int32_T *iidx;
  emxArray_int32_T *iwork;
  emxArray_real_T *b_xwork;
  emxArray_real_T *vwork;
  emxArray_real_T *xwork;
  double c_xwork[256];
  double x4[4];
  double d;
  double d1;
  double *vwork_data;
  double *x_data;
  double *xwork_data;
  int c_iwork[256];
  int idx4[4];
  int b;
  int bLen;
  int bLen2;
  int b_b;
  int dim;
  int exitg1;
  int i;
  int i3;
  int i4;
  int ib;
  int j;
  int k;
  int nBlocks;
  int nPairs;
  int offset;
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
    i = x->size[0];
  } else {
    i = 1;
  }
  vlen = i - 1;
  i3 = vwork->size[0];
  vwork->size[0] = i;
  emxEnsureCapacity_real_T(vwork, i3);
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
  emxInit_int32_T(&iwork, 1);
  emxInit_real_T(&xwork, 1);
  emxInit_int32_T(&b_iwork, 1);
  emxInit_real_T(&b_xwork, 1);
  for (j = 0; j < vstride; j++) {
    for (k = 0; k <= vlen; k++) {
      vwork_data[k] = x_data[j + k * vstride];
    }
    i = iidx->size[0];
    iidx->size[0] = vwork->size[0];
    emxEnsureCapacity_int32_T(iidx, i);
    iidx_data = iidx->data;
    dim = vwork->size[0];
    for (i = 0; i < dim; i++) {
      iidx_data[i] = 0;
    }
    i = vwork->size[0];
    if (vwork->size[0] != 0) {
      for (dim = 0; dim < 4; dim++) {
        x4[dim] = 0.0;
        idx4[dim] = 0;
      }
      i3 = iwork->size[0];
      iwork->size[0] = vwork->size[0];
      emxEnsureCapacity_int32_T(iwork, i3);
      iwork_data = iwork->data;
      dim = vwork->size[0];
      for (i3 = 0; i3 < dim; i3++) {
        iwork_data[i3] = 0;
      }
      i3 = xwork->size[0];
      xwork->size[0] = vwork->size[0];
      emxEnsureCapacity_real_T(xwork, i3);
      xwork_data = xwork->data;
      dim = vwork->size[0];
      for (i3 = 0; i3 < dim; i3++) {
        xwork_data[i3] = 0.0;
      }
      ib = 0;
      for (k = 0; k < i; k++) {
        ib++;
        idx4[ib - 1] = k + 1;
        x4[ib - 1] = vwork_data[k];
        if (ib == 4) {
          if (x4[0] >= x4[1]) {
            dim = 1;
            ib = 2;
          } else {
            dim = 2;
            ib = 1;
          }
          if (x4[2] >= x4[3]) {
            i3 = 3;
            i4 = 4;
          } else {
            i3 = 4;
            i4 = 3;
          }
          d = x4[dim - 1];
          d1 = x4[i3 - 1];
          if (d >= d1) {
            d = x4[ib - 1];
            if (d >= d1) {
              perm[0] = (signed char)dim;
              perm[1] = (signed char)ib;
              perm[2] = (signed char)i3;
              perm[3] = (signed char)i4;
            } else if (d >= x4[i4 - 1]) {
              perm[0] = (signed char)dim;
              perm[1] = (signed char)i3;
              perm[2] = (signed char)ib;
              perm[3] = (signed char)i4;
            } else {
              perm[0] = (signed char)dim;
              perm[1] = (signed char)i3;
              perm[2] = (signed char)i4;
              perm[3] = (signed char)ib;
            }
          } else {
            d1 = x4[i4 - 1];
            if (d >= d1) {
              if (x4[ib - 1] >= d1) {
                perm[0] = (signed char)i3;
                perm[1] = (signed char)dim;
                perm[2] = (signed char)ib;
                perm[3] = (signed char)i4;
              } else {
                perm[0] = (signed char)i3;
                perm[1] = (signed char)dim;
                perm[2] = (signed char)i4;
                perm[3] = (signed char)ib;
              }
            } else {
              perm[0] = (signed char)i3;
              perm[1] = (signed char)i4;
              perm[2] = (signed char)dim;
              perm[3] = (signed char)ib;
            }
          }
          iidx_data[k - 3] = idx4[perm[0] - 1];
          iidx_data[k - 2] = idx4[perm[1] - 1];
          iidx_data[k - 1] = idx4[perm[2] - 1];
          iidx_data[k] = idx4[perm[3] - 1];
          vwork_data[k - 3] = x4[perm[0] - 1];
          vwork_data[k - 2] = x4[perm[1] - 1];
          vwork_data[k - 1] = x4[perm[2] - 1];
          vwork_data[k] = x4[perm[3] - 1];
          ib = 0;
        }
      }
      if (ib > 0) {
        for (dim = 0; dim < 4; dim++) {
          perm[dim] = 0;
        }
        if (ib == 1) {
          perm[0] = 1;
        } else if (ib == 2) {
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
        for (k = 0; k < ib; k++) {
          i3 = perm[k] - 1;
          dim = (i - ib) + k;
          iidx_data[dim] = idx4[i3];
          vwork_data[dim] = x4[i3];
        }
      }
      ib = 2;
      if (i > 1) {
        if (i >= 256) {
          nBlocks = i >> 8;
          for (b = 0; b < nBlocks; b++) {
            offset = (b << 8) - 1;
            for (b_b = 0; b_b < 6; b_b++) {
              bLen = 1 << (b_b + 2);
              bLen2 = bLen << 1;
              nPairs = 256 >> (b_b + 3);
              for (k = 0; k < nPairs; k++) {
                i3 = (offset + k * bLen2) + 1;
                for (ib = 0; ib < bLen2; ib++) {
                  dim = i3 + ib;
                  c_iwork[ib] = iidx_data[dim];
                  c_xwork[ib] = vwork_data[dim];
                }
                i4 = 0;
                ib = bLen;
                dim = i3 - 1;
                do {
                  exitg1 = 0;
                  dim++;
                  if (c_xwork[i4] >= c_xwork[ib]) {
                    iidx_data[dim] = c_iwork[i4];
                    vwork_data[dim] = c_xwork[i4];
                    if (i4 + 1 < bLen) {
                      i4++;
                    } else {
                      exitg1 = 1;
                    }
                  } else {
                    iidx_data[dim] = c_iwork[ib];
                    vwork_data[dim] = c_xwork[ib];
                    if (ib + 1 < bLen2) {
                      ib++;
                    } else {
                      dim -= i4;
                      for (ib = i4 + 1; ib <= bLen; ib++) {
                        i3 = dim + ib;
                        iidx_data[i3] = c_iwork[ib - 1];
                        vwork_data[i3] = c_xwork[ib - 1];
                      }
                      exitg1 = 1;
                    }
                  }
                } while (exitg1 == 0);
              }
            }
          }
          dim = nBlocks << 8;
          ib = i - dim;
          if (ib > 0) {
            merge_block(iidx, vwork, dim, ib, 2, iwork, xwork);
            xwork_data = xwork->data;
            iwork_data = iwork->data;
          }
          ib = 8;
        }
        dim = iwork->size[0];
        i3 = b_iwork->size[0];
        b_iwork->size[0] = iwork->size[0];
        emxEnsureCapacity_int32_T(b_iwork, i3);
        iidx_data = b_iwork->data;
        for (i3 = 0; i3 < dim; i3++) {
          iidx_data[i3] = iwork_data[i3];
        }
        dim = xwork->size[0];
        i3 = b_xwork->size[0];
        b_xwork->size[0] = xwork->size[0];
        emxEnsureCapacity_real_T(b_xwork, i3);
        vwork_data = b_xwork->data;
        for (i3 = 0; i3 < dim; i3++) {
          vwork_data[i3] = xwork_data[i3];
        }
        merge_block(iidx, vwork, 0, i, ib, b_iwork, b_xwork);
        vwork_data = vwork->data;
        iidx_data = iidx->data;
      }
    }
    for (k = 0; k <= vlen; k++) {
      i = j + k * vstride;
      x_data[i] = vwork_data[k];
      idx_data[i] = iidx_data[k];
    }
  }
  emxFree_real_T(&b_xwork);
  emxFree_int32_T(&b_iwork);
  emxFree_real_T(&xwork);
  emxFree_int32_T(&iwork);
  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

/* End of code generation (sort.c) */
