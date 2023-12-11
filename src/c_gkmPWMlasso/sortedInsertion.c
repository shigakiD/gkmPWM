/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sortedInsertion.c
 *
 * Code generation for function 'sortedInsertion'
 *
 */

/* Include files */
#include "sortedInsertion.h"
#include "gkmPWMlasso4_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void sortedInsertion(double x, int ix, emxArray_real_T *b, int *nb, int blen,
                     emxArray_int32_T *idx)
{
  double *b_data;
  int ja;
  int jb;
  int jc;
  int *idx_data;
  idx_data = idx->data;
  b_data = b->data;
  if (*nb == 0) {
    *nb = 1;
    idx_data[0] = ix;
    b_data[0] = x;
  } else if (x > b_data[0]) {
    if (*nb < blen) {
      (*nb)++;
    }
    for (jc = *nb; jc >= 2; jc--) {
      idx_data[jc - 1] = idx_data[jc - 2];
      b_data[jc - 1] = b_data[jc - 2];
    }
    b_data[0] = x;
    idx_data[0] = ix;
  } else if ((*nb > 1) && (x > b_data[*nb - 1])) {
    ja = 1;
    jb = *nb;
    while (ja < jb) {
      jc = ja + ((jb - ja) >> 1);
      if (jc == ja) {
        ja = jb;
      } else if (x > b_data[jc - 1]) {
        jb = jc;
      } else {
        ja = jc;
      }
    }
    if (*nb < blen) {
      (*nb)++;
    }
    jb = ja + 1;
    for (jc = *nb; jc >= jb; jc--) {
      b_data[jc - 1] = b_data[jc - 2];
      idx_data[jc - 1] = idx_data[jc - 2];
    }
    b_data[ja - 1] = x;
    idx_data[ja - 1] = ix;
  } else if (*nb < blen) {
    (*nb)++;
    b_data[*nb - 1] = x;
    idx_data[*nb - 1] = ix;
  }
}

/* End of code generation (sortedInsertion.c) */
