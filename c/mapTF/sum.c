/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "mapTF_emxutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_sum(const emxArray_real_T *x, emxArray_real_T *y)
{
  const double *x_data;
  double bsum;
  double *y_data;
  unsigned int sz[2];
  int firstBlockLength;
  int hi;
  int ib;
  int k;
  int lastBlockLength;
  int nblocks;
  int npages;
  int xblockoffset;
  int xi;
  int xpageoffset;
  x_data = x->data;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    for (nblocks = 0; nblocks < 2; nblocks++) {
      sz[nblocks] = (unsigned int)x->size[nblocks];
    }
    nblocks = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int)sz[1];
    emxEnsureCapacity_real_T(y, nblocks);
    y_data = y->data;
    firstBlockLength = (int)sz[1];
    for (nblocks = 0; nblocks < firstBlockLength; nblocks++) {
      y_data[nblocks] = 0.0;
    }
  } else {
    npages = x->size[1];
    nblocks = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = x->size[1];
    emxEnsureCapacity_real_T(y, nblocks);
    y_data = y->data;
    if (x->size[0] <= 1024) {
      firstBlockLength = x->size[0];
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = x->size[0] / 1024;
      lastBlockLength = x->size[0] - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    for (xi = 0; xi < npages; xi++) {
      xpageoffset = xi * x->size[0];
      y_data[xi] = x_data[xpageoffset];
      for (k = 2; k <= firstBlockLength; k++) {
        y_data[xi] += x_data[(xpageoffset + k) - 1];
      }
      for (ib = 2; ib <= nblocks; ib++) {
        xblockoffset = xpageoffset + ((ib - 1) << 10);
        bsum = x_data[xblockoffset];
        if (ib == nblocks) {
          hi = lastBlockLength;
        } else {
          hi = 1024;
        }
        for (k = 2; k <= hi; k++) {
          bsum += x_data[(xblockoffset + k) - 1];
        }
        y_data[xi] += bsum;
      }
    }
  }
}

/*
 *
 */
double c_sum(const emxArray_real_T *x)
{
  const double *x_data;
  double bsum;
  double y;
  int firstBlockLength;
  int hi;
  int ib;
  int k;
  int lastBlockLength;
  int nblocks;
  x_data = x->data;
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    if (x->size[1] <= 1024) {
      firstBlockLength = x->size[1];
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = x->size[1] / 1024;
      lastBlockLength = x->size[1] - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    y = x_data[0];
    for (k = 2; k <= firstBlockLength; k++) {
      y += x_data[k - 1];
    }
    for (ib = 2; ib <= nblocks; ib++) {
      firstBlockLength = (ib - 1) << 10;
      bsum = x_data[firstBlockLength];
      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (k = 2; k <= hi; k++) {
        bsum += x_data[(firstBlockLength + k) - 1];
      }
      y += bsum;
    }
  }
  return y;
}

/*
 *
 */
void sum(const emxArray_real_T *x, emxArray_real_T *y)
{
  emxArray_real_T *bsum;
  const double *x_data;
  double *bsum_data;
  double *y_data;
  unsigned int sz[2];
  int bvstride;
  int firstBlockLength;
  int hi;
  int ib;
  int k;
  int lastBlockLength;
  int nblocks;
  int vstride;
  int xj;
  int xoffset;
  x_data = x->data;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    for (hi = 0; hi < 2; hi++) {
      sz[hi] = (unsigned int)x->size[hi];
    }
    hi = y->size[0];
    y->size[0] = (int)sz[0];
    emxEnsureCapacity_real_T(y, hi);
    y_data = y->data;
    firstBlockLength = (int)sz[0];
    for (hi = 0; hi < firstBlockLength; hi++) {
      y_data[hi] = 0.0;
    }
  } else {
    emxInit_real_T(&bsum, 1);
    vstride = x->size[0] - 1;
    bvstride = x->size[0] << 10;
    hi = y->size[0];
    y->size[0] = x->size[0];
    emxEnsureCapacity_real_T(y, hi);
    y_data = y->data;
    hi = bsum->size[0];
    bsum->size[0] = x->size[0];
    emxEnsureCapacity_real_T(bsum, hi);
    bsum_data = bsum->data;
    if (x->size[1] <= 1024) {
      firstBlockLength = x->size[1];
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = x->size[1] / 1024;
      lastBlockLength = x->size[1] - (nblocks << 10);
      if (lastBlockLength > 0) {
        nblocks++;
      } else {
        lastBlockLength = 1024;
      }
    }
    for (xj = 0; xj <= vstride; xj++) {
      y_data[xj] = x_data[xj];
      bsum_data[xj] = 0.0;
    }
    for (k = 2; k <= firstBlockLength; k++) {
      xoffset = (k - 1) * (vstride + 1);
      for (xj = 0; xj <= vstride; xj++) {
        y_data[xj] += x_data[xoffset + xj];
      }
    }
    for (ib = 2; ib <= nblocks; ib++) {
      firstBlockLength = (ib - 1) * bvstride;
      for (xj = 0; xj <= vstride; xj++) {
        bsum_data[xj] = x_data[firstBlockLength + xj];
      }
      if (ib == nblocks) {
        hi = lastBlockLength;
      } else {
        hi = 1024;
      }
      for (k = 2; k <= hi; k++) {
        xoffset = firstBlockLength + (k - 1) * (vstride + 1);
        for (xj = 0; xj <= vstride; xj++) {
          bsum_data[xj] += x_data[xoffset + xj];
        }
      }
      for (xj = 0; xj <= vstride; xj++) {
        y_data[xj] += bsum_data[xj];
      }
    }
    emxFree_real_T(&bsum);
  }
}

/* End of code generation (sum.c) */
