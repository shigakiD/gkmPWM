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
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Declarations */
static void b_merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                    int np, int nq, emxArray_int32_T *iwork,
                    emxArray_real_T *xwork);

static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);

static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
                        int n, int preSortLevel, emxArray_int32_T *iwork,
                        emxArray_real_T *xwork);

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
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
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
void b_merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int n,
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
  int b_b;
  int b_n;
  int exitg1;
  int i;
  int i2;
  int i3;
  int i4;
  int ib;
  int idx_tmp;
  int k;
  int n;
  int nBlocks;
  int nNonNaN;
  int *idx_data;
  int *iwork_data;
  signed char perm[4];
  x_data = x->data;
  i = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = x->size[1];
  emxEnsureCapacity_int32_T(idx, i);
  idx_data = idx->data;
  ib = x->size[1];
  for (i = 0; i < ib; i++) {
    idx_data[i] = 0;
  }
  if (x->size[1] != 0) {
    n = x->size[1];
    b_n = x->size[1];
    for (i = 0; i < 4; i++) {
      x4[i] = 0.0;
      idx4[i] = 0;
    }
    emxInit_int32_T(&iwork, 1);
    ib = x->size[1];
    i = iwork->size[0];
    iwork->size[0] = ib;
    emxEnsureCapacity_int32_T(iwork, i);
    iwork_data = iwork->data;
    for (i = 0; i < ib; i++) {
      iwork_data[i] = 0;
    }
    emxInit_real_T(&xwork, 1);
    ib = x->size[1];
    i = xwork->size[0];
    xwork->size[0] = ib;
    emxEnsureCapacity_real_T(xwork, i);
    xwork_data = xwork->data;
    for (i = 0; i < ib; i++) {
      xwork_data[i] = 0.0;
    }
    bLen = 0;
    ib = -1;
    for (k = 0; k < b_n; k++) {
      if (rtIsNaN(x_data[k])) {
        idx_tmp = (b_n - bLen) - 1;
        idx_data[idx_tmp] = k + 1;
        xwork_data[idx_tmp] = x_data[k];
        bLen++;
      } else {
        ib++;
        idx4[ib] = k + 1;
        x4[ib] = x_data[k];
        if (ib + 1 == 4) {
          i = k - bLen;
          if (x4[0] <= x4[1]) {
            ib = 1;
            i2 = 2;
          } else {
            ib = 2;
            i2 = 1;
          }
          if (x4[2] <= x4[3]) {
            i3 = 3;
            i4 = 4;
          } else {
            i3 = 4;
            i4 = 3;
          }
          d = x4[ib - 1];
          d1 = x4[i3 - 1];
          if (d <= d1) {
            d = x4[i2 - 1];
            if (d <= d1) {
              perm[0] = (signed char)ib;
              perm[1] = (signed char)i2;
              perm[2] = (signed char)i3;
              perm[3] = (signed char)i4;
            } else if (d <= x4[i4 - 1]) {
              perm[0] = (signed char)ib;
              perm[1] = (signed char)i3;
              perm[2] = (signed char)i2;
              perm[3] = (signed char)i4;
            } else {
              perm[0] = (signed char)ib;
              perm[1] = (signed char)i3;
              perm[2] = (signed char)i4;
              perm[3] = (signed char)i2;
            }
          } else {
            d1 = x4[i4 - 1];
            if (d <= d1) {
              if (x4[i2 - 1] <= d1) {
                perm[0] = (signed char)i3;
                perm[1] = (signed char)ib;
                perm[2] = (signed char)i2;
                perm[3] = (signed char)i4;
              } else {
                perm[0] = (signed char)i3;
                perm[1] = (signed char)ib;
                perm[2] = (signed char)i4;
                perm[3] = (signed char)i2;
              }
            } else {
              perm[0] = (signed char)i3;
              perm[1] = (signed char)i4;
              perm[2] = (signed char)ib;
              perm[3] = (signed char)i2;
            }
          }
          idx_data[i - 3] = idx4[perm[0] - 1];
          idx_data[i - 2] = idx4[perm[1] - 1];
          idx_data[i - 1] = idx4[perm[2] - 1];
          idx_data[i] = idx4[perm[3] - 1];
          x_data[i - 3] = x4[perm[0] - 1];
          x_data[i - 2] = x4[perm[1] - 1];
          x_data[i - 1] = x4[perm[2] - 1];
          x_data[i] = x4[perm[3] - 1];
          ib = -1;
        }
      }
    }
    i3 = (b_n - bLen) - 1;
    if (ib + 1 > 0) {
      for (i = 0; i < 4; i++) {
        perm[i] = 0;
      }
      if (ib + 1 == 1) {
        perm[0] = 1;
      } else if (ib + 1 == 2) {
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
      for (k = 0; k <= ib; k++) {
        idx_tmp = perm[k] - 1;
        i = (i3 - ib) + k;
        idx_data[i] = idx4[idx_tmp];
        x_data[i] = x4[idx_tmp];
      }
    }
    ib = (bLen >> 1) + 1;
    for (k = 0; k <= ib - 2; k++) {
      i = (i3 + k) + 1;
      i2 = idx_data[i];
      idx_tmp = (b_n - k) - 1;
      idx_data[i] = idx_data[idx_tmp];
      idx_data[idx_tmp] = i2;
      x_data[i] = xwork_data[idx_tmp];
      x_data[idx_tmp] = xwork_data[i];
    }
    if ((bLen & 1) != 0) {
      i = i3 + ib;
      x_data[i] = xwork_data[i];
    }
    nNonNaN = n - bLen;
    ib = 2;
    if (nNonNaN > 1) {
      if (n >= 256) {
        nBlocks = nNonNaN >> 8;
        if (nBlocks > 0) {
          for (b = 0; b < nBlocks; b++) {
            i4 = (b << 8) - 1;
            for (b_b = 0; b_b < 6; b_b++) {
              bLen = 1 << (b_b + 2);
              b_n = bLen << 1;
              n = 256 >> (b_b + 3);
              for (k = 0; k < n; k++) {
                i2 = (i4 + k * b_n) + 1;
                for (i = 0; i < b_n; i++) {
                  ib = i2 + i;
                  b_iwork[i] = idx_data[ib];
                  b_xwork[i] = x_data[ib];
                }
                i3 = 0;
                i = bLen;
                ib = i2 - 1;
                do {
                  exitg1 = 0;
                  ib++;
                  if (b_xwork[i3] <= b_xwork[i]) {
                    idx_data[ib] = b_iwork[i3];
                    x_data[ib] = b_xwork[i3];
                    if (i3 + 1 < bLen) {
                      i3++;
                    } else {
                      exitg1 = 1;
                    }
                  } else {
                    idx_data[ib] = b_iwork[i];
                    x_data[ib] = b_xwork[i];
                    if (i + 1 < b_n) {
                      i++;
                    } else {
                      ib -= i3;
                      for (i = i3 + 1; i <= bLen; i++) {
                        idx_tmp = ib + i;
                        idx_data[idx_tmp] = b_iwork[i - 1];
                        x_data[idx_tmp] = b_xwork[i - 1];
                      }
                      exitg1 = 1;
                    }
                  }
                } while (exitg1 == 0);
              }
            }
          }
          ib = nBlocks << 8;
          i = nNonNaN - ib;
          if (i > 0) {
            merge_block(idx, x, ib, i, 2, iwork, xwork);
          }
          ib = 8;
        }
      }
      merge_block(idx, x, 0, nNonNaN, ib, iwork, xwork);
    }
    emxFree_real_T(&xwork);
    emxFree_int32_T(&iwork);
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
  int b_b;
  int b_n;
  int exitg1;
  int i;
  int i3;
  int i4;
  int ib;
  int idx_tmp;
  int k;
  int n;
  int nBlocks;
  int nNonNaN;
  int quartetOffset;
  int *idx_data;
  int *iwork_data;
  signed char perm[4];
  x_data = x->data;
  quartetOffset = x->size[0];
  i = idx->size[0];
  idx->size[0] = quartetOffset;
  emxEnsureCapacity_int32_T(idx, i);
  idx_data = idx->data;
  for (i = 0; i < quartetOffset; i++) {
    idx_data[i] = 0;
  }
  if (x->size[0] != 0) {
    n = x->size[0];
    b_n = x->size[0];
    for (i = 0; i < 4; i++) {
      x4[i] = 0.0;
      idx4[i] = 0;
    }
    emxInit_int32_T(&iwork, 1);
    i = iwork->size[0];
    iwork->size[0] = quartetOffset;
    emxEnsureCapacity_int32_T(iwork, i);
    iwork_data = iwork->data;
    for (i = 0; i < quartetOffset; i++) {
      iwork_data[i] = 0;
    }
    emxInit_real_T(&xwork, 1);
    ib = x->size[0];
    i = xwork->size[0];
    xwork->size[0] = ib;
    emxEnsureCapacity_real_T(xwork, i);
    xwork_data = xwork->data;
    for (i = 0; i < ib; i++) {
      xwork_data[i] = 0.0;
    }
    bLen = 0;
    ib = -1;
    for (k = 0; k < b_n; k++) {
      if (rtIsNaN(x_data[k])) {
        idx_tmp = (b_n - bLen) - 1;
        idx_data[idx_tmp] = k + 1;
        xwork_data[idx_tmp] = x_data[k];
        bLen++;
      } else {
        ib++;
        idx4[ib] = k + 1;
        x4[ib] = x_data[k];
        if (ib + 1 == 4) {
          quartetOffset = k - bLen;
          if (x4[0] <= x4[1]) {
            ib = 1;
            i = 2;
          } else {
            ib = 2;
            i = 1;
          }
          if (x4[2] <= x4[3]) {
            i3 = 3;
            i4 = 4;
          } else {
            i3 = 4;
            i4 = 3;
          }
          d = x4[ib - 1];
          d1 = x4[i3 - 1];
          if (d <= d1) {
            d = x4[i - 1];
            if (d <= d1) {
              perm[0] = (signed char)ib;
              perm[1] = (signed char)i;
              perm[2] = (signed char)i3;
              perm[3] = (signed char)i4;
            } else if (d <= x4[i4 - 1]) {
              perm[0] = (signed char)ib;
              perm[1] = (signed char)i3;
              perm[2] = (signed char)i;
              perm[3] = (signed char)i4;
            } else {
              perm[0] = (signed char)ib;
              perm[1] = (signed char)i3;
              perm[2] = (signed char)i4;
              perm[3] = (signed char)i;
            }
          } else {
            d1 = x4[i4 - 1];
            if (d <= d1) {
              if (x4[i - 1] <= d1) {
                perm[0] = (signed char)i3;
                perm[1] = (signed char)ib;
                perm[2] = (signed char)i;
                perm[3] = (signed char)i4;
              } else {
                perm[0] = (signed char)i3;
                perm[1] = (signed char)ib;
                perm[2] = (signed char)i4;
                perm[3] = (signed char)i;
              }
            } else {
              perm[0] = (signed char)i3;
              perm[1] = (signed char)i4;
              perm[2] = (signed char)ib;
              perm[3] = (signed char)i;
            }
          }
          idx_data[quartetOffset - 3] = idx4[perm[0] - 1];
          idx_data[quartetOffset - 2] = idx4[perm[1] - 1];
          idx_data[quartetOffset - 1] = idx4[perm[2] - 1];
          idx_data[quartetOffset] = idx4[perm[3] - 1];
          x_data[quartetOffset - 3] = x4[perm[0] - 1];
          x_data[quartetOffset - 2] = x4[perm[1] - 1];
          x_data[quartetOffset - 1] = x4[perm[2] - 1];
          x_data[quartetOffset] = x4[perm[3] - 1];
          ib = -1;
        }
      }
    }
    i3 = (b_n - bLen) - 1;
    if (ib + 1 > 0) {
      for (i = 0; i < 4; i++) {
        perm[i] = 0;
      }
      if (ib + 1 == 1) {
        perm[0] = 1;
      } else if (ib + 1 == 2) {
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
      for (k = 0; k <= ib; k++) {
        idx_tmp = perm[k] - 1;
        quartetOffset = (i3 - ib) + k;
        idx_data[quartetOffset] = idx4[idx_tmp];
        x_data[quartetOffset] = x4[idx_tmp];
      }
    }
    ib = (bLen >> 1) + 1;
    for (k = 0; k <= ib - 2; k++) {
      quartetOffset = (i3 + k) + 1;
      i = idx_data[quartetOffset];
      idx_tmp = (b_n - k) - 1;
      idx_data[quartetOffset] = idx_data[idx_tmp];
      idx_data[idx_tmp] = i;
      x_data[quartetOffset] = xwork_data[idx_tmp];
      x_data[idx_tmp] = xwork_data[quartetOffset];
    }
    if ((bLen & 1) != 0) {
      i = i3 + ib;
      x_data[i] = xwork_data[i];
    }
    nNonNaN = n - bLen;
    ib = 2;
    if (nNonNaN > 1) {
      if (n >= 256) {
        nBlocks = nNonNaN >> 8;
        if (nBlocks > 0) {
          for (b = 0; b < nBlocks; b++) {
            i4 = (b << 8) - 1;
            for (b_b = 0; b_b < 6; b_b++) {
              bLen = 1 << (b_b + 2);
              b_n = bLen << 1;
              n = 256 >> (b_b + 3);
              for (k = 0; k < n; k++) {
                i = (i4 + k * b_n) + 1;
                for (quartetOffset = 0; quartetOffset < b_n; quartetOffset++) {
                  ib = i + quartetOffset;
                  b_iwork[quartetOffset] = idx_data[ib];
                  b_xwork[quartetOffset] = x_data[ib];
                }
                i3 = 0;
                quartetOffset = bLen;
                ib = i - 1;
                do {
                  exitg1 = 0;
                  ib++;
                  if (b_xwork[i3] <= b_xwork[quartetOffset]) {
                    idx_data[ib] = b_iwork[i3];
                    x_data[ib] = b_xwork[i3];
                    if (i3 + 1 < bLen) {
                      i3++;
                    } else {
                      exitg1 = 1;
                    }
                  } else {
                    idx_data[ib] = b_iwork[quartetOffset];
                    x_data[ib] = b_xwork[quartetOffset];
                    if (quartetOffset + 1 < b_n) {
                      quartetOffset++;
                    } else {
                      ib -= i3;
                      for (quartetOffset = i3 + 1; quartetOffset <= bLen;
                           quartetOffset++) {
                        idx_tmp = ib + quartetOffset;
                        idx_data[idx_tmp] = b_iwork[quartetOffset - 1];
                        x_data[idx_tmp] = b_xwork[quartetOffset - 1];
                      }
                      exitg1 = 1;
                    }
                  }
                } while (exitg1 == 0);
              }
            }
          }
          ib = nBlocks << 8;
          quartetOffset = nNonNaN - ib;
          if (quartetOffset > 0) {
            merge_block(idx, x, ib, quartetOffset, 2, iwork, xwork);
          }
          ib = 8;
        }
      }
      merge_block(idx, x, 0, nNonNaN, ib, iwork, xwork);
    }
    emxFree_real_T(&xwork);
    emxFree_int32_T(&iwork);
  }
}

/* End of code generation (sortIdx.c) */
