/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * linspace.c
 *
 * Code generation for function 'linspace'
 *
 */

/* Include files */
#include "linspace.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void linspace(double d1, double d2, double n, double y_data[], int y_size[2])
{
  double delta1;
  double delta2;
  int k;
  int y_size_tmp;
  int y_tmp_tmp;
  delta1 = floor(n);
  y_size[0] = 1;
  y_size_tmp = (int)floor(n);
  y_size[1] = (int)delta1;
  y_tmp_tmp = (int)delta1 - 1;
  y_data[y_size_tmp - 1] = d2;
  y_data[0] = d1;
  if (d1 == -d2) {
    delta2 = d2 / ((double)(int)delta1 - 1.0);
    for (k = 2; k <= y_tmp_tmp; k++) {
      y_data[k - 1] = ((double)((k << 1) - (int)delta1) - 1.0) * delta2;
    }
    if (((int)delta1 & 1) == 1) {
      y_data[(int)delta1 >> 1] = 0.0;
    }
  } else if (((d1 < 0.0) != (d2 < 0.0)) &&
             ((fabs(d1) > 8.9884656743115785E+307) ||
              (fabs(d2) > 8.9884656743115785E+307))) {
    delta1 = d1 / ((double)y_size_tmp - 1.0);
    delta2 = d2 / ((double)y_size_tmp - 1.0);
    for (k = 0; k <= y_size_tmp - 3; k++) {
      y_data[k + 1] =
          (d1 + delta2 * ((double)k + 1.0)) - delta1 * ((double)k + 1.0);
    }
  } else {
    delta1 = (d2 - d1) / ((double)y_size_tmp - 1.0);
    for (k = 0; k <= y_size_tmp - 3; k++) {
      y_data[k + 1] = d1 + ((double)k + 1.0) * delta1;
    }
  }
}

/* End of code generation (linspace.c) */
