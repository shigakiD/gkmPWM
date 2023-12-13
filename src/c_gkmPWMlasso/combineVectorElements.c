/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * combineVectorElements.c
 *
 * Code generation for function 'combineVectorElements'
 *
 */

/* Include files */
#include "combineVectorElements.h"
#include "gkmPWMlasso_types.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void combineVectorElements(const emxArray_boolean_T *x, int y_data[],
                           int y_size[2])
{
  int i;
  int k;
  int npages;
  int vlen;
  int xpageoffset;
  const bool *x_data;
  x_data = x->data;
  vlen = x->size[0];
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    y_size[0] = 1;
    y_size[1] = x->size[1];
    vlen = x->size[1];
    if (0 <= vlen - 1) {
      memset(&y_data[0], 0, vlen * sizeof(int));
    }
  } else {
    npages = x->size[1];
    y_size[0] = 1;
    y_size[1] = x->size[1];
    for (i = 0; i < npages; i++) {
      xpageoffset = i * x->size[0];
      y_data[i] = x_data[xpageoffset];
      for (k = 2; k <= vlen; k++) {
        y_data[i] += x_data[(xpageoffset + k) - 1];
      }
    }
  }
}

/* End of code generation (combineVectorElements.c) */
