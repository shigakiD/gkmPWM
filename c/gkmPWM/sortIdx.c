/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sortIdx.c
 *
 * Code generation for function 'sortIdx'
 *
 */

/* Include files */
#include "sortIdx.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_types.h"
#include <string.h>

/* Function Declarations */
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                    int np, int nq, emxArray_int32_T *iwork,
                    emxArray_real_T *xwork);

static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                          int n, int preSortLevel, emxArray_int32_T *iwork,
                          emxArray_real_T *xwork);

static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);

/* Function Definitions */
/*
 *
 */
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                    int np, int nq, emxArray_int32_T *iwork,
                    emxArray_real_T *xwork)
{
  double *x_data;
  double *xwork_data;
  int exitg1;
  int iout;
  int j;
  int n_tmp;
  int p;
  int q;
  int *idx_data;
  int *iwork_data;
  xwork_data = xwork->data;
  iwork_data = iwork->data;
  x_data = x->data;
  idx_data = idx->data;
  if (nq != 0) {
    n_tmp = np + nq;
    for (j = 0; j < n_tmp; j++) {
      iout = offset + j;
      iwork_data[j] = idx_data[iout];
      xwork_data[j] = x_data[iout];
    }
    p = 0;
    q = np;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork_data[p] <= xwork_data[q]) {
        idx_data[iout] = iwork_data[p];
        x_data[iout] = xwork_data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx_data[iout] = iwork_data[q];
        x_data[iout] = xwork_data[q];
        if (q + 1 < n_tmp) {
          q++;
        } else {
          q = iout - p;
          for (j = p + 1; j <= np; j++) {
            iout = q + j;
            idx_data[iout] = iwork_data[j - 1];
            x_data[iout] = xwork_data[j - 1];
          }
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

/*
 *
 */
static void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                          int n, int preSortLevel, emxArray_int32_T *iwork,
                          emxArray_real_T *xwork)
{
  int bLen;
  int nPairs;
  int nTail;
  int tailOffset;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        b_merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }
    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 0; nTail < nPairs; nTail++) {
      b_merge(idx, x, offset + nTail * tailOffset, bLen, bLen, iwork, xwork);
    }
    bLen = tailOffset;
  }
  if (n > bLen) {
    b_merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

/*
 *
 */
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  double *x_data;
  double *xwork_data;
  int exitg1;
  int iout;
  int j;
  int n_tmp;
  int p;
  int q;
  int *idx_data;
  int *iwork_data;
  xwork_data = xwork->data;
  iwork_data = iwork->data;
  x_data = x->data;
  idx_data = idx->data;
  if (nq != 0) {
    n_tmp = np + nq;
    for (j = 0; j < n_tmp; j++) {
      iout = offset + j;
      iwork_data[j] = idx_data[iout];
      xwork_data[j] = x_data[iout];
    }
    p = 0;
    q = np;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork_data[p] >= xwork_data[q]) {
        idx_data[iout] = iwork_data[p];
        x_data[iout] = xwork_data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx_data[iout] = iwork_data[q];
        x_data[iout] = xwork_data[q];
        if (q + 1 < n_tmp) {
          q++;
        } else {
          q = iout - p;
          for (j = p + 1; j <= np; j++) {
            iout = q + j;
            idx_data[iout] = iwork_data[j - 1];
            x_data[iout] = xwork_data[j - 1];
          }
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

/*
 *
 */
void b_sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_int32_T *iwork;
  emxArray_real_T *xwork;
  double b_xwork[256];
  double x4[4];
  double d;
  double d1;
  double *x_data;
  double *xwork_data;
  int b_iwork[256];
  int idx4[4];
  int b;
  int bLen;
  int bLen2;
  int b_b;
  int b_n;
  int exitg1;
  int i;
  int i1;
  int i4;
  int ib;
  int k;
  int n;
  int nPairs;
  int offset;
  int *idx_data;
  int *iwork_data;
  signed char perm[4];
  x_data = x->data;
  i1 = x->size[0];
  ib = idx->size[0];
  idx->size[0] = i1;
  emxEnsureCapacity_int32_T(idx, ib);
  idx_data = idx->data;
  for (ib = 0; ib < i1; ib++) {
    idx_data[ib] = 0;
  }
  if (x->size[0] != 0) {
    n = x->size[0];
    b_n = x->size[0];
    for (i = 0; i < 4; i++) {
      x4[i] = 0.0;
      idx4[i] = 0;
    }
    emxInit_int32_T(&iwork, 1);
    ib = iwork->size[0];
    iwork->size[0] = i1;
    emxEnsureCapacity_int32_T(iwork, ib);
    iwork_data = iwork->data;
    for (ib = 0; ib < i1; ib++) {
      iwork_data[ib] = 0;
    }
    emxInit_real_T(&xwork, 1);
    i1 = x->size[0];
    ib = xwork->size[0];
    xwork->size[0] = i1;
    emxEnsureCapacity_real_T(xwork, ib);
    xwork_data = xwork->data;
    for (ib = 0; ib < i1; ib++) {
      xwork_data[ib] = 0.0;
    }
    ib = 0;
    for (k = 0; k < b_n; k++) {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = x_data[k];
      if (ib == 4) {
        if (x4[0] <= x4[1]) {
          i1 = 1;
          ib = 2;
        } else {
          i1 = 2;
          ib = 1;
        }
        if (x4[2] <= x4[3]) {
          i = 3;
          i4 = 4;
        } else {
          i = 4;
          i4 = 3;
        }
        d = x4[i1 - 1];
        d1 = x4[i - 1];
        if (d <= d1) {
          d = x4[ib - 1];
          if (d <= d1) {
            perm[0] = (signed char)i1;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)i;
            perm[3] = (signed char)i4;
          } else if (d <= x4[i4 - 1]) {
            perm[0] = (signed char)i1;
            perm[1] = (signed char)i;
            perm[2] = (signed char)ib;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)i1;
            perm[1] = (signed char)i;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)ib;
          }
        } else {
          d1 = x4[i4 - 1];
          if (d <= d1) {
            if (x4[ib - 1] <= d1) {
              perm[0] = (signed char)i;
              perm[1] = (signed char)i1;
              perm[2] = (signed char)ib;
              perm[3] = (signed char)i4;
            } else {
              perm[0] = (signed char)i;
              perm[1] = (signed char)i1;
              perm[2] = (signed char)i4;
              perm[3] = (signed char)ib;
            }
          } else {
            perm[0] = (signed char)i;
            perm[1] = (signed char)i4;
            perm[2] = (signed char)i1;
            perm[3] = (signed char)ib;
          }
        }
        idx_data[k - 3] = idx4[perm[0] - 1];
        idx_data[k - 2] = idx4[perm[1] - 1];
        idx_data[k - 1] = idx4[perm[2] - 1];
        idx_data[k] = idx4[perm[3] - 1];
        x_data[k - 3] = x4[perm[0] - 1];
        x_data[k - 2] = x4[perm[1] - 1];
        x_data[k - 1] = x4[perm[2] - 1];
        x_data[k] = x4[perm[3] - 1];
        ib = 0;
      }
    }
    if (ib > 0) {
      for (i = 0; i < 4; i++) {
        perm[i] = 0;
      }
      if (ib == 1) {
        perm[0] = 1;
      } else if (ib == 2) {
        if (x4[0] <= x4[1]) {
          perm[0] = 1;
          perm[1] = 2;
        } else {
          perm[0] = 2;
          perm[1] = 1;
        }
      } else if (x4[0] <= x4[1]) {
        if (x4[1] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 2;
          perm[2] = 3;
        } else if (x4[0] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 3;
          perm[2] = 2;
        } else {
          perm[0] = 3;
          perm[1] = 1;
          perm[2] = 2;
        }
      } else if (x4[0] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 1;
        perm[2] = 3;
      } else if (x4[1] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 3;
        perm[2] = 1;
      } else {
        perm[0] = 3;
        perm[1] = 2;
        perm[2] = 1;
      }
      for (k = 0; k < ib; k++) {
        i4 = perm[k] - 1;
        i1 = (b_n - ib) + k;
        idx_data[i1] = idx4[i4];
        x_data[i1] = x4[i4];
      }
    }
    i1 = 2;
    if (n > 1) {
      if (n >= 256) {
        b_n = n >> 8;
        for (b = 0; b < b_n; b++) {
          offset = (b << 8) - 1;
          for (b_b = 0; b_b < 6; b_b++) {
            bLen = 1 << (b_b + 2);
            bLen2 = bLen << 1;
            nPairs = 256 >> (b_b + 3);
            for (k = 0; k < nPairs; k++) {
              i4 = (offset + k * bLen2) + 1;
              for (i1 = 0; i1 < bLen2; i1++) {
                ib = i4 + i1;
                b_iwork[i1] = idx_data[ib];
                b_xwork[i1] = x_data[ib];
              }
              i = 0;
              i1 = bLen;
              ib = i4 - 1;
              do {
                exitg1 = 0;
                ib++;
                if (b_xwork[i] <= b_xwork[i1]) {
                  idx_data[ib] = b_iwork[i];
                  x_data[ib] = b_xwork[i];
                  if (i + 1 < bLen) {
                    i++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  idx_data[ib] = b_iwork[i1];
                  x_data[ib] = b_xwork[i1];
                  if (i1 + 1 < bLen2) {
                    i1++;
                  } else {
                    ib -= i;
                    for (i1 = i + 1; i1 <= bLen; i1++) {
                      i4 = ib + i1;
                      idx_data[i4] = b_iwork[i1 - 1];
                      x_data[i4] = b_xwork[i1 - 1];
                    }
                    exitg1 = 1;
                  }
                }
              } while (exitg1 == 0);
            }
          }
        }
        i1 = b_n << 8;
        ib = n - i1;
        if (ib > 0) {
          b_merge_block(idx, x, i1, ib, 2, iwork, xwork);
        }
        i1 = 8;
      }
      b_merge_block(idx, x, 0, n, i1, iwork, xwork);
    }
    emxFree_real_T(&xwork);
    emxFree_int32_T(&iwork);
  }
}

/*
 *
 */
void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int n,
                 int preSortLevel, emxArray_int32_T *iwork,
                 emxArray_real_T *xwork)
{
  int bLen;
  int nPairs;
  int nTail;
  int tailOffset;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }
    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 0; nTail < nPairs; nTail++) {
      merge(idx, x, offset + nTail * tailOffset, bLen, bLen, iwork, xwork);
    }
    bLen = tailOffset;
  }
  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

/*
 *
 */
void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_int32_T *iwork;
  emxArray_real_T *xwork;
  double b_xwork[256];
  double x4[4];
  double d;
  double d1;
  double *x_data;
  double *xwork_data;
  int b_iwork[256];
  int idx4[4];
  int b;
  int bLen;
  int bLen2;
  int b_b;
  int b_n;
  int exitg1;
  int i;
  int i3;
  int i4;
  int ib;
  int j;
  int k;
  int n;
  int nPairs;
  int offset;
  int *idx_data;
  int *iwork_data;
  signed char perm[4];
  x_data = x->data;
  ib = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = x->size[1];
  emxEnsureCapacity_int32_T(idx, ib);
  idx_data = idx->data;
  i = x->size[1];
  for (ib = 0; ib < i; ib++) {
    idx_data[ib] = 0;
  }
  if (x->size[1] != 0) {
    n = x->size[1];
    b_n = x->size[1];
    for (i = 0; i < 4; i++) {
      x4[i] = 0.0;
      idx4[i] = 0;
    }
    emxInit_int32_T(&iwork, 1);
    i = x->size[1];
    ib = iwork->size[0];
    iwork->size[0] = i;
    emxEnsureCapacity_int32_T(iwork, ib);
    iwork_data = iwork->data;
    for (ib = 0; ib < i; ib++) {
      iwork_data[ib] = 0;
    }
    emxInit_real_T(&xwork, 1);
    i = x->size[1];
    ib = xwork->size[0];
    xwork->size[0] = i;
    emxEnsureCapacity_real_T(xwork, ib);
    xwork_data = xwork->data;
    for (ib = 0; ib < i; ib++) {
      xwork_data[ib] = 0.0;
    }
    ib = 0;
    for (k = 0; k < b_n; k++) {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = x_data[k];
      if (ib == 4) {
        if (x4[0] <= x4[1]) {
          i = 1;
          ib = 2;
        } else {
          i = 2;
          ib = 1;
        }
        if (x4[2] <= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }
        d = x4[i - 1];
        d1 = x4[i3 - 1];
        if (d <= d1) {
          d = x4[ib - 1];
          if (d <= d1) {
            perm[0] = (signed char)i;
            perm[1] = (signed char)ib;
            perm[2] = (signed char)i3;
            perm[3] = (signed char)i4;
          } else if (d <= x4[i4 - 1]) {
            perm[0] = (signed char)i;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)ib;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)i;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)ib;
          }
        } else {
          d1 = x4[i4 - 1];
          if (d <= d1) {
            if (x4[ib - 1] <= d1) {
              perm[0] = (signed char)i3;
              perm[1] = (signed char)i;
              perm[2] = (signed char)ib;
              perm[3] = (signed char)i4;
            } else {
              perm[0] = (signed char)i3;
              perm[1] = (signed char)i;
              perm[2] = (signed char)i4;
              perm[3] = (signed char)ib;
            }
          } else {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)i4;
            perm[2] = (signed char)i;
            perm[3] = (signed char)ib;
          }
        }
        idx_data[k - 3] = idx4[perm[0] - 1];
        idx_data[k - 2] = idx4[perm[1] - 1];
        idx_data[k - 1] = idx4[perm[2] - 1];
        idx_data[k] = idx4[perm[3] - 1];
        x_data[k - 3] = x4[perm[0] - 1];
        x_data[k - 2] = x4[perm[1] - 1];
        x_data[k - 1] = x4[perm[2] - 1];
        x_data[k] = x4[perm[3] - 1];
        ib = 0;
      }
    }
    if (ib > 0) {
      for (i = 0; i < 4; i++) {
        perm[i] = 0;
      }
      if (ib == 1) {
        perm[0] = 1;
      } else if (ib == 2) {
        if (x4[0] <= x4[1]) {
          perm[0] = 1;
          perm[1] = 2;
        } else {
          perm[0] = 2;
          perm[1] = 1;
        }
      } else if (x4[0] <= x4[1]) {
        if (x4[1] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 2;
          perm[2] = 3;
        } else if (x4[0] <= x4[2]) {
          perm[0] = 1;
          perm[1] = 3;
          perm[2] = 2;
        } else {
          perm[0] = 3;
          perm[1] = 1;
          perm[2] = 2;
        }
      } else if (x4[0] <= x4[2]) {
        perm[0] = 2;
        perm[1] = 1;
        perm[2] = 3;
      } else if (x4[1] <= x4[2]) {
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
        i = (b_n - ib) + k;
        idx_data[i] = idx4[i3];
        x_data[i] = x4[i3];
      }
    }
    i = 2;
    if (n > 1) {
      if (n >= 256) {
        b_n = n >> 8;
        for (b = 0; b < b_n; b++) {
          offset = (b << 8) - 1;
          for (b_b = 0; b_b < 6; b_b++) {
            bLen = 1 << (b_b + 2);
            bLen2 = bLen << 1;
            nPairs = 256 >> (b_b + 3);
            for (k = 0; k < nPairs; k++) {
              i4 = (offset + k * bLen2) + 1;
              for (j = 0; j < bLen2; j++) {
                ib = i4 + j;
                b_iwork[j] = idx_data[ib];
                b_xwork[j] = x_data[ib];
              }
              i = 0;
              i3 = bLen;
              ib = i4 - 1;
              do {
                exitg1 = 0;
                ib++;
                if (b_xwork[i] <= b_xwork[i3]) {
                  idx_data[ib] = b_iwork[i];
                  x_data[ib] = b_xwork[i];
                  if (i + 1 < bLen) {
                    i++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  idx_data[ib] = b_iwork[i3];
                  x_data[ib] = b_xwork[i3];
                  if (i3 + 1 < bLen2) {
                    i3++;
                  } else {
                    ib -= i;
                    for (j = i + 1; j <= bLen; j++) {
                      i3 = ib + j;
                      idx_data[i3] = b_iwork[j - 1];
                      x_data[i3] = b_xwork[j - 1];
                    }
                    exitg1 = 1;
                  }
                }
              } while (exitg1 == 0);
            }
          }
        }
        i = b_n << 8;
        ib = n - i;
        if (ib > 0) {
          b_merge_block(idx, x, i, ib, 2, iwork, xwork);
        }
        i = 8;
      }
      b_merge_block(idx, x, 0, n, i, iwork, xwork);
    }
    emxFree_real_T(&xwork);
    emxFree_int32_T(&iwork);
  }
}

/* End of code generation (sortIdx.c) */
