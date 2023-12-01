/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xhseqr.c
 *
 * Code generation for function 'xhseqr'
 *
 */

/* Include files */
#include "xhseqr.h"
#include "xdlanv2.h"
#include "xzlarfg.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
int xhseqr(double h[16])
{
  double v[3];
  double aa;
  double ab;
  double ba;
  double bb;
  double rt1r;
  double rt2r;
  double s;
  double s_tmp;
  double tst;
  int L;
  int b_i;
  int b_k;
  int hoffset;
  int i;
  int info;
  int its;
  int j;
  int k;
  int kdefl;
  int m;
  int nr;
  int sum1_tmp;
  bool exitg1;
  bool exitg2;
  bool exitg3;
  bool goto150;
  info = 0;
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
  h[2] = 0.0;
  h[3] = 0.0;
  h[7] = 0.0;
  kdefl = 0;
  i = 3;
  exitg1 = false;
  while ((!exitg1) && (i + 1 >= 1)) {
    L = 1;
    goto150 = false;
    its = 0;
    exitg2 = false;
    while ((!exitg2) && (its < 301)) {
      k = i;
      exitg3 = false;
      while ((!exitg3) && (k + 1 > L)) {
        b_i = k + ((k - 1) << 2);
        ba = fabs(h[b_i]);
        if (ba <= 4.0083367200179456E-292) {
          exitg3 = true;
        } else {
          nr = k + (k << 2);
          bb = fabs(h[nr]);
          tst = fabs(h[b_i - 1]) + bb;
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = fabs(h[(k + ((k - 2) << 2)) - 1]);
            }
            if (k + 2 <= 4) {
              tst += fabs(h[nr + 1]);
            }
          }
          if (ba <= 2.2204460492503131E-16 * tst) {
            tst = fabs(h[nr - 1]);
            if (ba > tst) {
              ab = ba;
              ba = tst;
            } else {
              ab = tst;
            }
            tst = fabs(h[b_i - 1] - h[nr]);
            if (bb > tst) {
              aa = bb;
              bb = tst;
            } else {
              aa = tst;
            }
            s = aa + ab;
            if (ba * (ab / s) <=
                fmax(4.0083367200179456E-292,
                     2.2204460492503131E-16 * (bb * (aa / s)))) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }
      }
      L = k + 1;
      if (k + 1 > 1) {
        h[k + ((k - 1) << 2)] = 0.0;
      }
      if (k + 1 >= i) {
        goto150 = true;
        exitg2 = true;
      } else {
        kdefl++;
        if (kdefl - kdefl / 20 * 20 == 0) {
          s = fabs(h[i + ((i - 1) << 2)]) + fabs(h[(i + ((i - 2) << 2)) - 1]);
          tst = 0.75 * s + h[i + (i << 2)];
          aa = -0.4375 * s;
          ab = s;
          bb = tst;
        } else if (kdefl - kdefl / 10 * 10 == 0) {
          s = fabs(h[(k + (k << 2)) + 1]) + fabs(h[(k + ((k + 1) << 2)) + 2]);
          tst = 0.75 * s + h[k + (k << 2)];
          aa = -0.4375 * s;
          ab = s;
          bb = tst;
        } else {
          nr = i + ((i - 1) << 2);
          tst = h[nr - 1];
          ab = h[nr];
          aa = h[(i + (i << 2)) - 1];
          bb = h[i + (i << 2)];
        }
        s = ((fabs(tst) + fabs(aa)) + fabs(ab)) + fabs(bb);
        if (s == 0.0) {
          rt1r = 0.0;
          ba = 0.0;
          rt2r = 0.0;
          bb = 0.0;
        } else {
          tst /= s;
          ab /= s;
          aa /= s;
          bb /= s;
          ba = (tst + bb) / 2.0;
          tst = (tst - ba) * (bb - ba) - aa * ab;
          ab = sqrt(fabs(tst));
          if (tst >= 0.0) {
            rt1r = ba * s;
            rt2r = rt1r;
            ba = ab * s;
            bb = -ba;
          } else {
            rt1r = ba + ab;
            rt2r = ba - ab;
            if (fabs(rt1r - bb) <= fabs(rt2r - bb)) {
              rt1r *= s;
              rt2r = rt1r;
            } else {
              rt2r *= s;
              rt1r = rt2r;
            }
            ba = 0.0;
            bb = 0.0;
          }
        }
        m = i - 1;
        exitg3 = false;
        while ((!exitg3) && (m >= k + 1)) {
          hoffset = m + ((m - 1) << 2);
          tst = h[hoffset];
          s_tmp = h[hoffset - 1];
          ab = s_tmp - rt2r;
          s = (fabs(ab) + fabs(bb)) + fabs(tst);
          aa = tst / s;
          nr = m + (m << 2);
          v[0] = (aa * h[nr - 1] + (s_tmp - rt1r) * (ab / s)) - ba * (bb / s);
          tst = h[nr];
          v[1] = aa * (((s_tmp + tst) - rt1r) - rt2r);
          v[2] = aa * h[nr + 1];
          s = (fabs(v[0]) + fabs(v[1])) + fabs(v[2]);
          v[0] /= s;
          v[1] /= s;
          v[2] /= s;
          if ((m == k + 1) || (fabs(h[m - 1]) * (fabs(v[1]) + fabs(v[2])) <=
                               2.2204460492503131E-16 * fabs(v[0]) *
                                   ((fabs(h[0]) + fabs(s_tmp)) + fabs(tst)))) {
            exitg3 = true;
          } else {
            m--;
          }
        }
        for (b_k = m; b_k <= i; b_k++) {
          nr = (i - b_k) + 2;
          if (3 <= nr) {
            nr = 3;
          }
          if (b_k > m) {
            hoffset = (b_k + ((b_k - 2) << 2)) - 1;
            for (j = 0; j < nr; j++) {
              v[j] = h[j + hoffset];
            }
          }
          tst = v[0];
          ba = xzlarfg(nr, &tst, v);
          v[0] = tst;
          if (b_k > m) {
            h[(b_k + ((b_k - 2) << 2)) - 1] = tst;
            h[b_k + ((b_k - 2) << 2)] = 0.0;
            if (b_k < i) {
              h[b_k + 1] = 0.0;
            }
          } else if (m > k + 1) {
            h[b_k - 1] *= 1.0 - ba;
          }
          s = v[1];
          ab = ba * v[1];
          if (nr == 3) {
            rt1r = v[2];
            tst = ba * v[2];
            for (j = b_k; j < 5; j++) {
              sum1_tmp = b_k + ((j - 1) << 2);
              aa = (h[sum1_tmp - 1] + s * h[sum1_tmp]) + rt1r * h[sum1_tmp + 1];
              h[sum1_tmp - 1] -= aa * ba;
              h[sum1_tmp] -= aa * ab;
              h[sum1_tmp + 1] -= aa * tst;
            }
            if (b_k + 3 <= i + 1) {
              b_i = b_k + 2;
            } else {
              b_i = i;
            }
            for (j = 0; j <= b_i; j++) {
              sum1_tmp = j + ((b_k - 1) << 2);
              nr = j + (b_k << 2);
              hoffset = j + ((b_k + 1) << 2);
              aa = (h[sum1_tmp] + s * h[nr]) + rt1r * h[hoffset];
              h[sum1_tmp] -= aa * ba;
              h[nr] -= aa * ab;
              h[hoffset] -= aa * tst;
            }
          } else if (nr == 2) {
            for (j = b_k; j < 5; j++) {
              sum1_tmp = b_k + ((j - 1) << 2);
              tst = h[sum1_tmp - 1];
              aa = tst + s * h[sum1_tmp];
              h[sum1_tmp - 1] = tst - aa * ba;
              h[sum1_tmp] -= aa * ab;
            }
            for (j = 0; j <= i; j++) {
              sum1_tmp = j + ((b_k - 1) << 2);
              nr = j + (b_k << 2);
              aa = h[sum1_tmp] + s * h[nr];
              h[sum1_tmp] -= aa * ba;
              h[nr] -= aa * ab;
            }
          }
        }
        its++;
      }
    }
    if (!goto150) {
      info = i + 1;
      exitg1 = true;
    } else {
      if ((L != i + 1) && (L == i)) {
        b_i = i << 2;
        nr = i + b_i;
        s = h[nr - 1];
        b_k = (i - 1) << 2;
        hoffset = i + b_k;
        rt1r = h[hoffset];
        tst = h[nr];
        xdlanv2(&h[(i + ((i - 1) << 2)) - 1], &s, &rt1r, &tst, &ab, &aa, &ba,
                &bb, &s_tmp, &rt2r);
        h[nr - 1] = s;
        h[hoffset] = rt1r;
        h[nr] = tst;
        if (4 > i + 1) {
          nr = 2 - i;
          hoffset = i + ((i + 1) << 2);
          for (k = 0; k <= nr; k++) {
            sum1_tmp = hoffset + (k << 2);
            tst = h[sum1_tmp];
            ab = h[sum1_tmp - 1];
            h[sum1_tmp] = s_tmp * tst - rt2r * ab;
            h[sum1_tmp - 1] = s_tmp * ab + rt2r * tst;
          }
        }
        if (i - 1 >= 1) {
          for (k = 0; k <= i - 2; k++) {
            nr = b_i + k;
            hoffset = b_k + k;
            tst = s_tmp * h[hoffset] + rt2r * h[nr];
            h[nr] = s_tmp * h[nr] - rt2r * h[hoffset];
            h[hoffset] = tst;
          }
        }
      }
      kdefl = 0;
      i = L - 2;
    }
  }
  h[3] = 0.0;
  return info;
}

/* End of code generation (xhseqr.c) */
