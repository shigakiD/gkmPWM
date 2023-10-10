/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * blockedSummation.c
 *
 * Code generation for function 'blockedSummation'
 *
 */

/* Include files */
#include "blockedSummation.h"
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
double blockedSummation(const emxArray_real_T *x, int vlen)
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
  if ((x->size[0] == 0) || (vlen == 0)) {
    y = 0.0;
  } else {
    if (vlen <= 1024) {
      firstBlockLength = vlen;
      lastBlockLength = 0;
      nblocks = 1;
    } else {
      firstBlockLength = 1024;
      nblocks = vlen / 1024;
      lastBlockLength = vlen - (nblocks << 10);
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
void colMajorFlatIter(const emxArray_real_T *x, int vlen, emxArray_real_T *y)
{
  const double *x_data;
  double bsum;
  double *y_data;
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
  npages = x->size[1];
  firstBlockLength = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = x->size[1];
  emxEnsureCapacity_real_T(y, firstBlockLength);
  y_data = y->data;
  if (vlen <= 1024) {
    firstBlockLength = vlen;
    lastBlockLength = 0;
    nblocks = 1;
  } else {
    firstBlockLength = 1024;
    nblocks = vlen / 1024;
    lastBlockLength = vlen - (nblocks << 10);
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

/* End of code generation (blockedSummation.c) */
