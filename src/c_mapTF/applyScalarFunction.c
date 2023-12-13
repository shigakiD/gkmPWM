/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * applyScalarFunction.c
 *
 * Code generation for function 'applyScalarFunction'
 *
 */

/* Include files */
#include "applyScalarFunction.h"
#include "mapTF_emxutil.h"
#include "mapTF_rtwutil.h"
#include "mapTF_types.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void applyScalarFunction(const emxArray_real_T *x, emxArray_real_T *z1)
{
  const double *x_data;
  double R;
  double S;
  double absx;
  double s;
  double *z1_data;
  int b;
  int eint;
  int k;
  int ub_loop;
  x_data = x->data;
  ub_loop = z1->size[0];
  z1->size[0] = x->size[0];
  emxEnsureCapacity_real_T(z1, ub_loop);
  z1_data = z1->data;
  ub_loop = x->size[0] - 1;
#pragma omp parallel for num_threads(1 > omp_get_max_threads()                 \
                                         ? omp_get_max_threads()               \
                                         : 1) private(absx, s, S, R, eint, b)

  for (k = 0; k <= ub_loop; k++) {
    /* ========================== COPYRIGHT NOTICE ============================
     */
    /*  The algorithms for calculating ERF(X) and ERFC(X) are derived */
    /*  from FDLIBM, which has the following notice: */
    /*                                                                          */
    /*  Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved. */
    /*                                                                          */
    /*  Developed at SunSoft, a Sun Microsystems, Inc. business. */
    /*  Permission to use, copy, modify, and distribute this */
    /*  software is freely granted, provided that this notice */
    /*  is preserved. */
    /* =============================    END    ================================
     */
    absx = fabs(x_data[k]);
    if (rtIsNaN(x_data[k])) {
      z1_data[k] = x_data[k];
    } else if (rtIsInf(x_data[k])) {
      if (x_data[k] < 0.0) {
        z1_data[k] = -1.0;
      } else {
        z1_data[k] = 1.0;
      }
    } else if (absx < 0.84375) {
      if (absx < 3.7252902984619141E-9) {
        if (absx < 2.8480945388892178E-306) {
          z1_data[k] =
              0.125 * (8.0 * x_data[k] + 1.0270333367641007 * x_data[k]);
        } else {
          z1_data[k] = x_data[k] + 0.12837916709551259 * x_data[k];
        }
      } else {
        s = x_data[k] * x_data[k];
        z1_data[k] =
            x_data[k] +
            x_data[k] * ((s * (s * (s * (s * -2.3763016656650163E-5 +
                                         -0.0057702702964894416) +
                                    -0.02848174957559851) +
                               -0.3250421072470015) +
                          0.12837916709551256) /
                         (s * (s * (s * (s * (s * -3.9602282787753681E-6 +
                                              0.00013249473800432164) +
                                         0.0050813062818757656) +
                                    0.0650222499887673) +
                               0.39791722395915535) +
                          1.0));
      }
    } else if (absx < 1.25) {
      S = (absx - 1.0) *
              ((absx - 1.0) *
                   ((absx - 1.0) *
                        ((absx - 1.0) *
                             ((absx - 1.0) *
                                  ((absx - 1.0) * -0.0021663755948687908 +
                                   0.035478304325618236) +
                              -0.11089469428239668) +
                         0.31834661990116175) +
                    -0.37220787603570132) +
               0.41485611868374833) +
          -0.0023621185607526594;
      s = (absx - 1.0) *
              ((absx - 1.0) *
                   ((absx - 1.0) *
                        ((absx - 1.0) *
                             ((absx - 1.0) *
                                  ((absx - 1.0) * 0.011984499846799107 +
                                   0.013637083912029051) +
                              0.12617121980876164) +
                         0.071828654414196266) +
                    0.540397917702171) +
               0.10642088040084423) +
          1.0;
      if (x_data[k] >= 0.0) {
        z1_data[k] = S / s + 0.84506291151046753;
      } else {
        z1_data[k] = -0.84506291151046753 - S / s;
      }
    } else if (absx > 6.0) {
      if (x_data[k] < 0.0) {
        z1_data[k] = -1.0;
      } else {
        z1_data[k] = 1.0;
      }
    } else {
      s = 1.0 / (absx * absx);
      if (absx < 2.8571434020996094) {
        R = s * (s * (s * (s * (s * (s * (s * -9.8143293441691455 +
                                          -81.2874355063066) +
                                     -184.60509290671104) +
                                -162.39666946257347) +
                           -62.375332450326006) +
                      -10.558626225323291) +
                 -0.69385857270718176) +
            -0.0098649440348471482;
        S = s * (s * (s * (s * (s * (s * (s * (s * -0.0604244152148581 +
                                               6.5702497703192817) +
                                          108.63500554177944) +
                                     429.00814002756783) +
                                645.38727173326788) +
                           434.56587747522923) +
                      137.65775414351904) +
                 19.651271667439257) +
            1.0;
      } else {
        R = s * (s * (s * (s * (s * (s * -483.5191916086514 +
                                     -1025.0951316110772) +
                                -637.56644336838963) +
                           -160.63638485582192) +
                      -17.757954917754752) +
                 -0.799283237680523) +
            -0.0098649429247001;
        S = s * (s * (s * (s * (s * (s * (s * -22.440952446585818 +
                                          474.52854120695537) +
                                     2553.0504064331644) +
                                3199.8582195085955) +
                           1536.729586084437) +
                      325.79251299657392) +
                 30.338060743482458) +
            1.0;
      }
      if (!rtIsNaN(absx)) {
        s = frexp(absx, &eint);
        b = eint;
      } else {
        s = absx;
        b = 0;
      }
      s = floor(s * 2.097152E+6) / 2.097152E+6 * rt_powd_snf(2.0, b);
      s = exp(-s * s - 0.5625) * exp((s - absx) * (s + absx) + R / S) / absx;
      if (x_data[k] < 0.0) {
        z1_data[k] = s - 1.0;
      } else {
        z1_data[k] = 1.0 - s;
      }
    }
  }
}

/* End of code generation (applyScalarFunction.c) */
