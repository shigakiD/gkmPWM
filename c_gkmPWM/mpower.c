/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mpower.c
 *
 * Code generation for function 'mpower'
 *
 */

/* Include files */
#include "mpower.h"
#include "gkmPWM_data.h"
#include "gkmPWM_emxutil.h"
#include "gkmPWM_rtwutil.h"
#include "gkmPWM_types.h"
#include "cblas.h"
#include "lapacke.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void b_mpower(const creal_T a[16], creal_T c[16])
{
  creal_T x[16];
  double ai;
  double ar;
  double bi;
  double br;
  double brm;
  double re;
  double s;
  double smax;
  int b_i;
  int b_tmp;
  int i;
  int i1;
  int j;
  int jA;
  int jBcol;
  int jp1j;
  int k;
  int mmj_tmp;
  int pipk;
  signed char ipiv[4];
  signed char p[4];
  for (i = 0; i < 16; i++) {
    c[i].re = 0.0;
    c[i].im = 0.0;
    x[i] = a[i];
  }
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  for (j = 0; j < 3; j++) {
    mmj_tmp = 2 - j;
    b_tmp = j * 5;
    jp1j = b_tmp + 2;
    pipk = 4 - j;
    jA = 0;
    smax = fabs(x[b_tmp].re) + fabs(x[b_tmp].im);
    for (k = 2; k <= pipk; k++) {
      jBcol = (b_tmp + k) - 1;
      s = fabs(x[jBcol].re) + fabs(x[jBcol].im);
      if (s > smax) {
        jA = k - 1;
        smax = s;
      }
    }
    pipk = b_tmp + jA;
    if ((x[pipk].re != 0.0) || (x[pipk].im != 0.0)) {
      if (jA != 0) {
        pipk = j + jA;
        ipiv[j] = (signed char)(pipk + 1);
        smax = x[j].re;
        s = x[j].im;
        x[j] = x[pipk];
        x[pipk].re = smax;
        x[pipk].im = s;
        smax = x[j + 4].re;
        s = x[j + 4].im;
        x[j + 4] = x[pipk + 4];
        x[pipk + 4].re = smax;
        x[pipk + 4].im = s;
        smax = x[j + 8].re;
        s = x[j + 8].im;
        x[j + 8] = x[pipk + 8];
        x[pipk + 8].re = smax;
        x[pipk + 8].im = s;
        smax = x[j + 12].re;
        s = x[j + 12].im;
        x[j + 12] = x[pipk + 12];
        x[pipk + 12].re = smax;
        x[pipk + 12].im = s;
      }
      i = (b_tmp - j) + 4;
      for (b_i = jp1j; b_i <= i; b_i++) {
        ar = x[b_i - 1].re;
        ai = x[b_i - 1].im;
        br = x[b_tmp].re;
        bi = x[b_tmp].im;
        if (bi == 0.0) {
          if (ai == 0.0) {
            re = ar / br;
            smax = 0.0;
          } else if (ar == 0.0) {
            re = 0.0;
            smax = ai / br;
          } else {
            re = ar / br;
            smax = ai / br;
          }
        } else if (br == 0.0) {
          if (ar == 0.0) {
            re = ai / bi;
            smax = 0.0;
          } else if (ai == 0.0) {
            re = 0.0;
            smax = -(ar / bi);
          } else {
            re = ai / bi;
            smax = -(ar / bi);
          }
        } else {
          brm = fabs(br);
          smax = fabs(bi);
          if (brm > smax) {
            s = bi / br;
            smax = br + s * bi;
            re = (ar + s * ai) / smax;
            smax = (ai - s * ar) / smax;
          } else if (smax == brm) {
            if (br > 0.0) {
              s = 0.5;
            } else {
              s = -0.5;
            }
            if (bi > 0.0) {
              smax = 0.5;
            } else {
              smax = -0.5;
            }
            re = (ar * s + ai * smax) / brm;
            smax = (ai * s - ar * smax) / brm;
          } else {
            s = br / bi;
            smax = bi + s * br;
            re = (s * ar + ai) / smax;
            smax = (s * ai - ar) / smax;
          }
        }
        x[b_i - 1].re = re;
        x[b_i - 1].im = smax;
      }
    }
    jA = b_tmp;
    for (jp1j = 0; jp1j <= mmj_tmp; jp1j++) {
      pipk = (b_tmp + (jp1j << 2)) + 4;
      smax = x[pipk].re;
      s = x[pipk].im;
      if ((smax != 0.0) || (s != 0.0)) {
        i = jA + 6;
        i1 = (jA - j) + 8;
        for (pipk = i; pipk <= i1; pipk++) {
          jBcol = ((b_tmp + pipk) - jA) - 5;
          br = x[jBcol].re;
          ar = x[jBcol].im;
          x[pipk - 1].re += br * -smax - ar * -s;
          x[pipk - 1].im += br * -s + ar * -smax;
        }
      }
      jA += 4;
    }
  }
  p[0] = 1;
  p[1] = 2;
  p[2] = 3;
  p[3] = 4;
  if (ipiv[0] > 1) {
    pipk = p[ipiv[0] - 1];
    p[ipiv[0] - 1] = 1;
    p[0] = (signed char)pipk;
  }
  if (ipiv[1] > 2) {
    pipk = p[ipiv[1] - 1];
    p[ipiv[1] - 1] = p[1];
    p[1] = (signed char)pipk;
  }
  if (ipiv[2] > 3) {
    pipk = p[ipiv[2] - 1];
    p[ipiv[2] - 1] = p[2];
    p[2] = (signed char)pipk;
  }
  for (k = 0; k < 4; k++) {
    i = (p[k] - 1) << 2;
    i1 = k + i;
    c[i1].re = 1.0;
    c[i1].im = 0.0;
    for (j = k + 1; j < 5; j++) {
      jp1j = (j + i) - 1;
      if ((c[jp1j].re != 0.0) || (c[jp1j].im != 0.0)) {
        i1 = j + 1;
        for (b_i = i1; b_i < 5; b_i++) {
          pipk = (b_i + ((j - 1) << 2)) - 1;
          smax = x[pipk].im;
          s = x[pipk].re;
          br = c[jp1j].re;
          pipk = (b_i + i) - 1;
          c[pipk].re -= c[jp1j].re * s - c[jp1j].im * smax;
          c[pipk].im -= br * smax + c[jp1j].im * s;
        }
      }
    }
  }
  for (j = 0; j < 4; j++) {
    jBcol = j << 2;
    for (k = 3; k >= 0; k--) {
      jA = k << 2;
      jp1j = k + jBcol;
      ar = c[jp1j].re;
      ai = c[jp1j].im;
      if ((ar != 0.0) || (ai != 0.0)) {
        pipk = k + jA;
        br = x[pipk].re;
        bi = x[pipk].im;
        if (bi == 0.0) {
          if (ai == 0.0) {
            re = ar / br;
            smax = 0.0;
          } else if (ar == 0.0) {
            re = 0.0;
            smax = ai / br;
          } else {
            re = ar / br;
            smax = ai / br;
          }
        } else if (br == 0.0) {
          if (ar == 0.0) {
            re = ai / bi;
            smax = 0.0;
          } else if (ai == 0.0) {
            re = 0.0;
            smax = -(ar / bi);
          } else {
            re = ai / bi;
            smax = -(ar / bi);
          }
        } else {
          brm = fabs(br);
          smax = fabs(bi);
          if (brm > smax) {
            s = bi / br;
            smax = br + s * bi;
            re = (ar + s * ai) / smax;
            smax = (ai - s * ar) / smax;
          } else if (smax == brm) {
            if (br > 0.0) {
              s = 0.5;
            } else {
              s = -0.5;
            }
            if (bi > 0.0) {
              smax = 0.5;
            } else {
              smax = -0.5;
            }
            re = (ar * s + ai * smax) / brm;
            smax = (ai * s - ar * smax) / brm;
          } else {
            s = br / bi;
            smax = bi + s * br;
            re = (s * ar + ai) / smax;
            smax = (s * ai - ar) / smax;
          }
        }
        c[jp1j].re = re;
        c[jp1j].im = smax;
        for (b_i = 0; b_i < k; b_i++) {
          pipk = b_i + jA;
          smax = x[pipk].im;
          s = x[pipk].re;
          br = c[jp1j].re;
          i = b_i + jBcol;
          c[i].re -= c[jp1j].re * s - c[jp1j].im * smax;
          c[i].im -= br * smax + c[jp1j].im * s;
        }
      }
    }
  }
}

/*
 *
 */
void c_mpower(const double a[16], double c[16])
{
  double x[16];
  double s;
  double smax;
  int b_i;
  int b_tmp;
  int i;
  int i1;
  int j;
  int jA;
  int jp1j;
  int k;
  int kAcol;
  int mmj_tmp;
  signed char ipiv[4];
  signed char p[4];
  for (i = 0; i < 16; i++) {
    c[i] = 0.0;
    x[i] = a[i];
  }
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  for (j = 0; j < 3; j++) {
    mmj_tmp = 2 - j;
    b_tmp = j * 5;
    jp1j = b_tmp + 2;
    jA = 4 - j;
    kAcol = 0;
    smax = fabs(x[b_tmp]);
    for (k = 2; k <= jA; k++) {
      s = fabs(x[(b_tmp + k) - 1]);
      if (s > smax) {
        kAcol = k - 1;
        smax = s;
      }
    }
    if (x[b_tmp + kAcol] != 0.0) {
      if (kAcol != 0) {
        jA = j + kAcol;
        ipiv[j] = (signed char)(jA + 1);
        smax = x[j];
        x[j] = x[jA];
        x[jA] = smax;
        smax = x[j + 4];
        x[j + 4] = x[jA + 4];
        x[jA + 4] = smax;
        smax = x[j + 8];
        x[j + 8] = x[jA + 8];
        x[jA + 8] = smax;
        smax = x[j + 12];
        x[j + 12] = x[jA + 12];
        x[jA + 12] = smax;
      }
      i = (b_tmp - j) + 4;
      for (b_i = jp1j; b_i <= i; b_i++) {
        x[b_i - 1] /= x[b_tmp];
      }
    }
    jA = b_tmp;
    for (kAcol = 0; kAcol <= mmj_tmp; kAcol++) {
      smax = x[(b_tmp + (kAcol << 2)) + 4];
      if (smax != 0.0) {
        i = jA + 6;
        i1 = (jA - j) + 8;
        for (jp1j = i; jp1j <= i1; jp1j++) {
          x[jp1j - 1] += x[((b_tmp + jp1j) - jA) - 5] * -smax;
        }
      }
      jA += 4;
    }
  }
  p[0] = 1;
  p[1] = 2;
  p[2] = 3;
  p[3] = 4;
  if (ipiv[0] > 1) {
    jA = p[ipiv[0] - 1];
    p[ipiv[0] - 1] = 1;
    p[0] = (signed char)jA;
  }
  if (ipiv[1] > 2) {
    jA = p[ipiv[1] - 1];
    p[ipiv[1] - 1] = p[1];
    p[1] = (signed char)jA;
  }
  if (ipiv[2] > 3) {
    jA = p[ipiv[2] - 1];
    p[ipiv[2] - 1] = p[2];
    p[2] = (signed char)jA;
  }
  for (k = 0; k < 4; k++) {
    jp1j = (p[k] - 1) << 2;
    c[k + jp1j] = 1.0;
    for (j = k + 1; j < 5; j++) {
      i = (j + jp1j) - 1;
      if (c[i] != 0.0) {
        i1 = j + 1;
        for (b_i = i1; b_i < 5; b_i++) {
          jA = (b_i + jp1j) - 1;
          c[jA] -= c[i] * x[(b_i + ((j - 1) << 2)) - 1];
        }
      }
    }
  }
  for (j = 0; j < 4; j++) {
    jA = j << 2;
    for (k = 3; k >= 0; k--) {
      kAcol = k << 2;
      i = k + jA;
      smax = c[i];
      if (smax != 0.0) {
        c[i] = smax / x[k + kAcol];
        for (b_i = 0; b_i < k; b_i++) {
          jp1j = b_i + jA;
          c[jp1j] -= c[i] * x[b_i + kAcol];
        }
      }
    }
  }
}

/*
 *
 */
void d_mpower(const creal_T a[9], creal_T c[9])
{
  creal_T x[9];
  double absx11;
  double absx21;
  double absx31;
  double brm;
  double d;
  double im;
  double re;
  double t1_im;
  double t1_re;
  double t2_im;
  double t2_re;
  int itmp;
  int p1;
  int p2;
  int p3;
  memcpy(&x[0], &a[0], 9U * sizeof(creal_T));
  p1 = 0;
  p2 = 3;
  p3 = 6;
  absx11 = rt_hypotd(a[0].re, a[0].im);
  absx21 = rt_hypotd(a[1].re, a[1].im);
  absx31 = rt_hypotd(a[2].re, a[2].im);
  if ((absx21 > absx11) && (absx21 > absx31)) {
    p1 = 3;
    p2 = 0;
    x[0] = a[1];
    x[1] = a[0];
    x[3] = a[4];
    x[4] = a[3];
    x[6] = a[7];
    x[7] = a[6];
  } else if (absx31 > absx11) {
    p1 = 6;
    p3 = 0;
    x[0] = a[2];
    x[2] = a[0];
    x[3] = a[5];
    x[5] = a[3];
    x[6] = a[8];
    x[8] = a[6];
  }
  if (x[0].im == 0.0) {
    if (x[1].im == 0.0) {
      re = x[1].re / x[0].re;
      im = 0.0;
    } else if (x[1].re == 0.0) {
      re = 0.0;
      im = x[1].im / x[0].re;
    } else {
      re = x[1].re / x[0].re;
      im = x[1].im / x[0].re;
    }
  } else if (x[0].re == 0.0) {
    if (x[1].re == 0.0) {
      re = x[1].im / x[0].im;
      im = 0.0;
    } else if (x[1].im == 0.0) {
      re = 0.0;
      im = -(x[1].re / x[0].im);
    } else {
      re = x[1].im / x[0].im;
      im = -(x[1].re / x[0].im);
    }
  } else {
    brm = fabs(x[0].re);
    absx31 = fabs(x[0].im);
    if (brm > absx31) {
      absx31 = x[0].im / x[0].re;
      d = x[0].re + absx31 * x[0].im;
      re = (x[1].re + absx31 * x[1].im) / d;
      im = (x[1].im - absx31 * x[1].re) / d;
    } else if (absx31 == brm) {
      if (x[0].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[0].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      re = (x[1].re * absx31 + x[1].im * d) / brm;
      im = (x[1].im * absx31 - x[1].re * d) / brm;
    } else {
      absx31 = x[0].re / x[0].im;
      d = x[0].im + absx31 * x[0].re;
      re = (absx31 * x[1].re + x[1].im) / d;
      im = (absx31 * x[1].im - x[1].re) / d;
    }
  }
  x[1].re = re;
  x[1].im = im;
  if (x[0].im == 0.0) {
    if (x[2].im == 0.0) {
      absx21 = x[2].re / x[0].re;
      absx11 = 0.0;
    } else if (x[2].re == 0.0) {
      absx21 = 0.0;
      absx11 = x[2].im / x[0].re;
    } else {
      absx21 = x[2].re / x[0].re;
      absx11 = x[2].im / x[0].re;
    }
  } else if (x[0].re == 0.0) {
    if (x[2].re == 0.0) {
      absx21 = x[2].im / x[0].im;
      absx11 = 0.0;
    } else if (x[2].im == 0.0) {
      absx21 = 0.0;
      absx11 = -(x[2].re / x[0].im);
    } else {
      absx21 = x[2].im / x[0].im;
      absx11 = -(x[2].re / x[0].im);
    }
  } else {
    brm = fabs(x[0].re);
    absx31 = fabs(x[0].im);
    if (brm > absx31) {
      absx31 = x[0].im / x[0].re;
      d = x[0].re + absx31 * x[0].im;
      absx21 = (x[2].re + absx31 * x[2].im) / d;
      absx11 = (x[2].im - absx31 * x[2].re) / d;
    } else if (absx31 == brm) {
      if (x[0].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[0].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      absx21 = (x[2].re * absx31 + x[2].im * d) / brm;
      absx11 = (x[2].im * absx31 - x[2].re * d) / brm;
    } else {
      absx31 = x[0].re / x[0].im;
      d = x[0].im + absx31 * x[0].re;
      absx21 = (absx31 * x[2].re + x[2].im) / d;
      absx11 = (absx31 * x[2].im - x[2].re) / d;
    }
  }
  x[2].re = absx21;
  x[2].im = absx11;
  x[4].re -= re * x[3].re - im * x[3].im;
  x[4].im -= re * x[3].im + im * x[3].re;
  x[5].re -= absx21 * x[3].re - absx11 * x[3].im;
  x[5].im -= absx21 * x[3].im + absx11 * x[3].re;
  x[7].re -= re * x[6].re - im * x[6].im;
  x[7].im -= re * x[6].im + im * x[6].re;
  x[8].re -= absx21 * x[6].re - absx11 * x[6].im;
  x[8].im -= absx21 * x[6].im + absx11 * x[6].re;
  if (rt_hypotd(x[5].re, x[5].im) > rt_hypotd(x[4].re, x[4].im)) {
    itmp = p2;
    p2 = p3;
    p3 = itmp;
    x[1] = x[2];
    x[2].re = re;
    x[2].im = im;
    t1_re = x[4].re;
    t1_im = x[4].im;
    x[4] = x[5];
    x[5].re = t1_re;
    x[5].im = t1_im;
    t1_re = x[7].re;
    t1_im = x[7].im;
    x[7] = x[8];
    x[8].re = t1_re;
    x[8].im = t1_im;
  }
  if (x[4].im == 0.0) {
    if (x[5].im == 0.0) {
      re = x[5].re / x[4].re;
      im = 0.0;
    } else if (x[5].re == 0.0) {
      re = 0.0;
      im = x[5].im / x[4].re;
    } else {
      re = x[5].re / x[4].re;
      im = x[5].im / x[4].re;
    }
  } else if (x[4].re == 0.0) {
    if (x[5].re == 0.0) {
      re = x[5].im / x[4].im;
      im = 0.0;
    } else if (x[5].im == 0.0) {
      re = 0.0;
      im = -(x[5].re / x[4].im);
    } else {
      re = x[5].im / x[4].im;
      im = -(x[5].re / x[4].im);
    }
  } else {
    brm = fabs(x[4].re);
    absx31 = fabs(x[4].im);
    if (brm > absx31) {
      absx31 = x[4].im / x[4].re;
      d = x[4].re + absx31 * x[4].im;
      re = (x[5].re + absx31 * x[5].im) / d;
      im = (x[5].im - absx31 * x[5].re) / d;
    } else if (absx31 == brm) {
      if (x[4].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[4].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      re = (x[5].re * absx31 + x[5].im * d) / brm;
      im = (x[5].im * absx31 - x[5].re * d) / brm;
    } else {
      absx31 = x[4].re / x[4].im;
      d = x[4].im + absx31 * x[4].re;
      re = (absx31 * x[5].re + x[5].im) / d;
      im = (absx31 * x[5].im - x[5].re) / d;
    }
  }
  x[8].re -= re * x[7].re - im * x[7].im;
  x[8].im -= re * x[7].im + im * x[7].re;
  absx11 = (x[1].re * re - x[1].im * im) - x[2].re;
  absx21 = (x[1].re * im + x[1].im * re) - x[2].im;
  if (x[8].im == 0.0) {
    if (absx21 == 0.0) {
      t1_re = absx11 / x[8].re;
      t1_im = 0.0;
    } else if (absx11 == 0.0) {
      t1_re = 0.0;
      t1_im = absx21 / x[8].re;
    } else {
      t1_re = absx11 / x[8].re;
      t1_im = absx21 / x[8].re;
    }
  } else if (x[8].re == 0.0) {
    if (absx11 == 0.0) {
      t1_re = absx21 / x[8].im;
      t1_im = 0.0;
    } else if (absx21 == 0.0) {
      t1_re = 0.0;
      t1_im = -(absx11 / x[8].im);
    } else {
      t1_re = absx21 / x[8].im;
      t1_im = -(absx11 / x[8].im);
    }
  } else {
    brm = fabs(x[8].re);
    absx31 = fabs(x[8].im);
    if (brm > absx31) {
      absx31 = x[8].im / x[8].re;
      d = x[8].re + absx31 * x[8].im;
      t1_re = (absx11 + absx31 * absx21) / d;
      t1_im = (absx21 - absx31 * absx11) / d;
    } else if (absx31 == brm) {
      if (x[8].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[8].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      t1_re = (absx11 * absx31 + absx21 * d) / brm;
      t1_im = (absx21 * absx31 - absx11 * d) / brm;
    } else {
      absx31 = x[8].re / x[8].im;
      d = x[8].im + absx31 * x[8].re;
      t1_re = (absx31 * absx11 + absx21) / d;
      t1_im = (absx31 * absx21 - absx11) / d;
    }
  }
  absx11 = -(x[1].re + (x[7].re * t1_re - x[7].im * t1_im));
  absx21 = -(x[1].im + (x[7].re * t1_im + x[7].im * t1_re));
  if (x[4].im == 0.0) {
    if (absx21 == 0.0) {
      t2_re = absx11 / x[4].re;
      t2_im = 0.0;
    } else if (absx11 == 0.0) {
      t2_re = 0.0;
      t2_im = absx21 / x[4].re;
    } else {
      t2_re = absx11 / x[4].re;
      t2_im = absx21 / x[4].re;
    }
  } else if (x[4].re == 0.0) {
    if (absx11 == 0.0) {
      t2_re = absx21 / x[4].im;
      t2_im = 0.0;
    } else if (absx21 == 0.0) {
      t2_re = 0.0;
      t2_im = -(absx11 / x[4].im);
    } else {
      t2_re = absx21 / x[4].im;
      t2_im = -(absx11 / x[4].im);
    }
  } else {
    brm = fabs(x[4].re);
    absx31 = fabs(x[4].im);
    if (brm > absx31) {
      absx31 = x[4].im / x[4].re;
      d = x[4].re + absx31 * x[4].im;
      t2_re = (absx11 + absx31 * absx21) / d;
      t2_im = (absx21 - absx31 * absx11) / d;
    } else if (absx31 == brm) {
      if (x[4].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[4].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      t2_re = (absx11 * absx31 + absx21 * d) / brm;
      t2_im = (absx21 * absx31 - absx11 * d) / brm;
    } else {
      absx31 = x[4].re / x[4].im;
      d = x[4].im + absx31 * x[4].re;
      t2_re = (absx31 * absx11 + absx21) / d;
      t2_im = (absx31 * absx21 - absx11) / d;
    }
  }
  absx11 = (1.0 - (x[3].re * t2_re - x[3].im * t2_im)) -
           (x[6].re * t1_re - x[6].im * t1_im);
  absx21 = (0.0 - (x[3].re * t2_im + x[3].im * t2_re)) -
           (x[6].re * t1_im + x[6].im * t1_re);
  if (x[0].im == 0.0) {
    if (absx21 == 0.0) {
      c[p1].re = absx11 / x[0].re;
      c[p1].im = 0.0;
    } else if (absx11 == 0.0) {
      c[p1].re = 0.0;
      c[p1].im = absx21 / x[0].re;
    } else {
      c[p1].re = absx11 / x[0].re;
      c[p1].im = absx21 / x[0].re;
    }
  } else if (x[0].re == 0.0) {
    if (absx11 == 0.0) {
      c[p1].re = absx21 / x[0].im;
      c[p1].im = 0.0;
    } else if (absx21 == 0.0) {
      c[p1].re = 0.0;
      c[p1].im = -(absx11 / x[0].im);
    } else {
      c[p1].re = absx21 / x[0].im;
      c[p1].im = -(absx11 / x[0].im);
    }
  } else {
    brm = fabs(x[0].re);
    absx31 = fabs(x[0].im);
    if (brm > absx31) {
      absx31 = x[0].im / x[0].re;
      d = x[0].re + absx31 * x[0].im;
      c[p1].re = (absx11 + absx31 * absx21) / d;
      c[p1].im = (absx21 - absx31 * absx11) / d;
    } else if (absx31 == brm) {
      if (x[0].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[0].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      c[p1].re = (absx11 * absx31 + absx21 * d) / brm;
      c[p1].im = (absx21 * absx31 - absx11 * d) / brm;
    } else {
      absx31 = x[0].re / x[0].im;
      d = x[0].im + absx31 * x[0].re;
      c[p1].re = (absx31 * absx11 + absx21) / d;
      c[p1].im = (absx31 * absx21 - absx11) / d;
    }
  }
  c[p1 + 1].re = t2_re;
  c[p1 + 1].im = t2_im;
  c[p1 + 2].re = t1_re;
  c[p1 + 2].im = t1_im;
  if (x[8].im == 0.0) {
    if (-im == 0.0) {
      t1_re = -re / x[8].re;
      t1_im = 0.0;
    } else if (-re == 0.0) {
      t1_re = 0.0;
      t1_im = -im / x[8].re;
    } else {
      t1_re = -re / x[8].re;
      t1_im = -im / x[8].re;
    }
  } else if (x[8].re == 0.0) {
    if (-re == 0.0) {
      t1_re = -im / x[8].im;
      t1_im = 0.0;
    } else if (-im == 0.0) {
      t1_re = 0.0;
      t1_im = -(-re / x[8].im);
    } else {
      t1_re = -im / x[8].im;
      t1_im = -(-re / x[8].im);
    }
  } else {
    brm = fabs(x[8].re);
    absx31 = fabs(x[8].im);
    if (brm > absx31) {
      absx31 = x[8].im / x[8].re;
      d = x[8].re + absx31 * x[8].im;
      t1_re = (-re + absx31 * -im) / d;
      t1_im = (-im - absx31 * -re) / d;
    } else if (absx31 == brm) {
      if (x[8].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[8].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      t1_re = (-re * absx31 + -im * d) / brm;
      t1_im = (-im * absx31 - -re * d) / brm;
    } else {
      absx31 = x[8].re / x[8].im;
      d = x[8].im + absx31 * x[8].re;
      t1_re = (absx31 * -re + -im) / d;
      t1_im = (absx31 * -im - (-re)) / d;
    }
  }
  absx11 = x[7].re * t1_re - x[7].im * t1_im;
  absx21 = x[7].re * t1_im + x[7].im * t1_re;
  if (x[4].im == 0.0) {
    if (0.0 - absx21 == 0.0) {
      t2_re = (1.0 - absx11) / x[4].re;
      t2_im = 0.0;
    } else if (1.0 - absx11 == 0.0) {
      t2_re = 0.0;
      t2_im = (0.0 - absx21) / x[4].re;
    } else {
      t2_re = (1.0 - absx11) / x[4].re;
      t2_im = (0.0 - absx21) / x[4].re;
    }
  } else if (x[4].re == 0.0) {
    if (1.0 - absx11 == 0.0) {
      t2_re = (0.0 - absx21) / x[4].im;
      t2_im = 0.0;
    } else if (0.0 - absx21 == 0.0) {
      t2_re = 0.0;
      t2_im = -((1.0 - absx11) / x[4].im);
    } else {
      t2_re = (0.0 - absx21) / x[4].im;
      t2_im = -((1.0 - absx11) / x[4].im);
    }
  } else {
    brm = fabs(x[4].re);
    absx31 = fabs(x[4].im);
    if (brm > absx31) {
      absx31 = x[4].im / x[4].re;
      d = x[4].re + absx31 * x[4].im;
      t2_re = ((1.0 - absx11) + absx31 * (0.0 - absx21)) / d;
      t2_im = ((0.0 - absx21) - absx31 * (1.0 - absx11)) / d;
    } else if (absx31 == brm) {
      if (x[4].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[4].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      t2_re = ((1.0 - absx11) * absx31 + (0.0 - absx21) * d) / brm;
      t2_im = ((0.0 - absx21) * absx31 - (1.0 - absx11) * d) / brm;
    } else {
      absx31 = x[4].re / x[4].im;
      d = x[4].im + absx31 * x[4].re;
      t2_re = (absx31 * (1.0 - absx11) + (0.0 - absx21)) / d;
      t2_im = (absx31 * (0.0 - absx21) - (1.0 - absx11)) / d;
    }
  }
  absx11 = -((x[3].re * t2_re - x[3].im * t2_im) +
             (x[6].re * t1_re - x[6].im * t1_im));
  absx21 = -((x[3].re * t2_im + x[3].im * t2_re) +
             (x[6].re * t1_im + x[6].im * t1_re));
  if (x[0].im == 0.0) {
    if (absx21 == 0.0) {
      c[p2].re = absx11 / x[0].re;
      c[p2].im = 0.0;
    } else if (absx11 == 0.0) {
      c[p2].re = 0.0;
      c[p2].im = absx21 / x[0].re;
    } else {
      c[p2].re = absx11 / x[0].re;
      c[p2].im = absx21 / x[0].re;
    }
  } else if (x[0].re == 0.0) {
    if (absx11 == 0.0) {
      c[p2].re = absx21 / x[0].im;
      c[p2].im = 0.0;
    } else if (absx21 == 0.0) {
      c[p2].re = 0.0;
      c[p2].im = -(absx11 / x[0].im);
    } else {
      c[p2].re = absx21 / x[0].im;
      c[p2].im = -(absx11 / x[0].im);
    }
  } else {
    brm = fabs(x[0].re);
    absx31 = fabs(x[0].im);
    if (brm > absx31) {
      absx31 = x[0].im / x[0].re;
      d = x[0].re + absx31 * x[0].im;
      c[p2].re = (absx11 + absx31 * absx21) / d;
      c[p2].im = (absx21 - absx31 * absx11) / d;
    } else if (absx31 == brm) {
      if (x[0].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[0].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      c[p2].re = (absx11 * absx31 + absx21 * d) / brm;
      c[p2].im = (absx21 * absx31 - absx11 * d) / brm;
    } else {
      absx31 = x[0].re / x[0].im;
      d = x[0].im + absx31 * x[0].re;
      c[p2].re = (absx31 * absx11 + absx21) / d;
      c[p2].im = (absx31 * absx21 - absx11) / d;
    }
  }
  c[p2 + 1].re = t2_re;
  c[p2 + 1].im = t2_im;
  c[p2 + 2].re = t1_re;
  c[p2 + 2].im = t1_im;
  if (x[8].im == 0.0) {
    t1_re = 1.0 / x[8].re;
    t1_im = 0.0;
  } else if (x[8].re == 0.0) {
    t1_re = 0.0;
    t1_im = -(1.0 / x[8].im);
  } else {
    brm = fabs(x[8].re);
    absx31 = fabs(x[8].im);
    if (brm > absx31) {
      absx31 = x[8].im / x[8].re;
      d = x[8].re + absx31 * x[8].im;
      t1_re = 1.0 / d;
      t1_im = (0.0 - absx31) / d;
    } else if (absx31 == brm) {
      if (x[8].re > 0.0) {
        absx11 = 0.5;
      } else {
        absx11 = -0.5;
      }
      t1_re = absx11 / brm;
      if (x[8].im > 0.0) {
        absx11 = 0.5;
      } else {
        absx11 = -0.5;
      }
      t1_im = (0.0 - absx11) / brm;
    } else {
      absx31 = x[8].re / x[8].im;
      d = x[8].im + absx31 * x[8].re;
      t1_re = absx31 / d;
      t1_im = -1.0 / d;
    }
  }
  absx11 = -x[7].re * t1_re - -x[7].im * t1_im;
  absx21 = -x[7].re * t1_im + -x[7].im * t1_re;
  if (x[4].im == 0.0) {
    if (absx21 == 0.0) {
      t2_re = absx11 / x[4].re;
      t2_im = 0.0;
    } else if (absx11 == 0.0) {
      t2_re = 0.0;
      t2_im = absx21 / x[4].re;
    } else {
      t2_re = absx11 / x[4].re;
      t2_im = absx21 / x[4].re;
    }
  } else if (x[4].re == 0.0) {
    if (absx11 == 0.0) {
      t2_re = absx21 / x[4].im;
      t2_im = 0.0;
    } else if (absx21 == 0.0) {
      t2_re = 0.0;
      t2_im = -(absx11 / x[4].im);
    } else {
      t2_re = absx21 / x[4].im;
      t2_im = -(absx11 / x[4].im);
    }
  } else {
    brm = fabs(x[4].re);
    absx31 = fabs(x[4].im);
    if (brm > absx31) {
      absx31 = x[4].im / x[4].re;
      d = x[4].re + absx31 * x[4].im;
      t2_re = (absx11 + absx31 * absx21) / d;
      t2_im = (absx21 - absx31 * absx11) / d;
    } else if (absx31 == brm) {
      if (x[4].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[4].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      t2_re = (absx11 * absx31 + absx21 * d) / brm;
      t2_im = (absx21 * absx31 - absx11 * d) / brm;
    } else {
      absx31 = x[4].re / x[4].im;
      d = x[4].im + absx31 * x[4].re;
      t2_re = (absx31 * absx11 + absx21) / d;
      t2_im = (absx31 * absx21 - absx11) / d;
    }
  }
  absx11 = -((x[3].re * t2_re - x[3].im * t2_im) +
             (x[6].re * t1_re - x[6].im * t1_im));
  absx21 = -((x[3].re * t2_im + x[3].im * t2_re) +
             (x[6].re * t1_im + x[6].im * t1_re));
  if (x[0].im == 0.0) {
    if (absx21 == 0.0) {
      c[p3].re = absx11 / x[0].re;
      c[p3].im = 0.0;
    } else if (absx11 == 0.0) {
      c[p3].re = 0.0;
      c[p3].im = absx21 / x[0].re;
    } else {
      c[p3].re = absx11 / x[0].re;
      c[p3].im = absx21 / x[0].re;
    }
  } else if (x[0].re == 0.0) {
    if (absx11 == 0.0) {
      c[p3].re = absx21 / x[0].im;
      c[p3].im = 0.0;
    } else if (absx21 == 0.0) {
      c[p3].re = 0.0;
      c[p3].im = -(absx11 / x[0].im);
    } else {
      c[p3].re = absx21 / x[0].im;
      c[p3].im = -(absx11 / x[0].im);
    }
  } else {
    brm = fabs(x[0].re);
    absx31 = fabs(x[0].im);
    if (brm > absx31) {
      absx31 = x[0].im / x[0].re;
      d = x[0].re + absx31 * x[0].im;
      c[p3].re = (absx11 + absx31 * absx21) / d;
      c[p3].im = (absx21 - absx31 * absx11) / d;
    } else if (absx31 == brm) {
      if (x[0].re > 0.0) {
        absx31 = 0.5;
      } else {
        absx31 = -0.5;
      }
      if (x[0].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      c[p3].re = (absx11 * absx31 + absx21 * d) / brm;
      c[p3].im = (absx21 * absx31 - absx11 * d) / brm;
    } else {
      absx31 = x[0].re / x[0].im;
      d = x[0].im + absx31 * x[0].re;
      c[p3].re = (absx31 * absx11 + absx21) / d;
      c[p3].im = (absx31 * absx21 - absx11) / d;
    }
  }
  c[p3 + 1].re = t2_re;
  c[p3 + 1].im = t2_im;
  c[p3 + 2].re = t1_re;
  c[p3 + 2].im = t1_im;
}

/*
 *
 */
void e_mpower(const creal_T a[4], creal_T c[4])
{
  double a_re;
  double b_bim;
  double b_brm;
  double bim;
  double brm;
  double d;
  double r_im;
  double r_re;
  double s;
  double sgnbr;
  double t_im;
  double t_re;
  brm = fabs(a[0].re);
  b_brm = fabs(a[1].re);
  bim = fabs(a[0].im);
  b_bim = fabs(a[1].im);
  if (b_brm + b_bim > brm + bim) {
    if (a[1].im == 0.0) {
      if (a[0].im == 0.0) {
        r_re = a[0].re / a[1].re;
        r_im = 0.0;
      } else if (a[0].re == 0.0) {
        r_re = 0.0;
        r_im = a[0].im / a[1].re;
      } else {
        r_re = a[0].re / a[1].re;
        r_im = a[0].im / a[1].re;
      }
    } else if (a[1].re == 0.0) {
      if (a[0].re == 0.0) {
        r_re = a[0].im / a[1].im;
        r_im = 0.0;
      } else if (a[0].im == 0.0) {
        r_re = 0.0;
        r_im = -(a[0].re / a[1].im);
      } else {
        r_re = a[0].im / a[1].im;
        r_im = -(a[0].re / a[1].im);
      }
    } else if (b_brm > b_bim) {
      s = a[1].im / a[1].re;
      d = a[1].re + s * a[1].im;
      r_re = (a[0].re + s * a[0].im) / d;
      r_im = (a[0].im - s * a[0].re) / d;
    } else if (b_bim == b_brm) {
      if (a[1].re > 0.0) {
        sgnbr = 0.5;
      } else {
        sgnbr = -0.5;
      }
      if (a[1].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      r_re = (a[0].re * sgnbr + a[0].im * d) / b_brm;
      r_im = (a[0].im * sgnbr - a[0].re * d) / b_brm;
    } else {
      s = a[1].re / a[1].im;
      d = a[1].im + s * a[1].re;
      r_re = (s * a[0].re + a[0].im) / d;
      r_im = (s * a[0].im - a[0].re) / d;
    }
    d = (r_re * a[3].re - r_im * a[3].im) - a[2].re;
    sgnbr = (r_re * a[3].im + r_im * a[3].re) - a[2].im;
    if (sgnbr == 0.0) {
      t_re = 1.0 / d;
      t_im = 0.0;
    } else if (d == 0.0) {
      t_re = 0.0;
      t_im = -(1.0 / sgnbr);
    } else {
      brm = fabs(d);
      bim = fabs(sgnbr);
      if (brm > bim) {
        s = sgnbr / d;
        d += s * sgnbr;
        t_re = 1.0 / d;
        t_im = (0.0 - s) / d;
      } else if (bim == brm) {
        if (d > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }
        t_re = d / brm;
        if (sgnbr > 0.0) {
          sgnbr = 0.5;
        } else {
          sgnbr = -0.5;
        }
        t_im = (0.0 - sgnbr) / brm;
      } else {
        s = d / sgnbr;
        d = sgnbr + s * d;
        t_re = s / d;
        t_im = -1.0 / d;
      }
    }
    if (a[1].im == 0.0) {
      if (a[3].im == 0.0) {
        a_re = a[3].re / a[1].re;
        d = 0.0;
      } else if (a[3].re == 0.0) {
        a_re = 0.0;
        d = a[3].im / a[1].re;
      } else {
        a_re = a[3].re / a[1].re;
        d = a[3].im / a[1].re;
      }
    } else if (a[1].re == 0.0) {
      if (a[3].re == 0.0) {
        a_re = a[3].im / a[1].im;
        d = 0.0;
      } else if (a[3].im == 0.0) {
        a_re = 0.0;
        d = -(a[3].re / a[1].im);
      } else {
        a_re = a[3].im / a[1].im;
        d = -(a[3].re / a[1].im);
      }
    } else if (b_brm > b_bim) {
      s = a[1].im / a[1].re;
      d = a[1].re + s * a[1].im;
      a_re = (a[3].re + s * a[3].im) / d;
      d = (a[3].im - s * a[3].re) / d;
    } else if (b_bim == b_brm) {
      if (a[1].re > 0.0) {
        sgnbr = 0.5;
      } else {
        sgnbr = -0.5;
      }
      if (a[1].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      a_re = (a[3].re * sgnbr + a[3].im * d) / b_brm;
      d = (a[3].im * sgnbr - a[3].re * d) / b_brm;
    } else {
      s = a[1].re / a[1].im;
      d = a[1].im + s * a[1].re;
      a_re = (s * a[3].re + a[3].im) / d;
      d = (s * a[3].im - a[3].re) / d;
    }
    c[0].re = a_re * t_re - d * t_im;
    c[0].im = a_re * t_im + d * t_re;
    c[1].re = -t_re;
    c[1].im = -t_im;
    if (a[1].im == 0.0) {
      if (-a[2].im == 0.0) {
        a_re = -a[2].re / a[1].re;
        d = 0.0;
      } else if (-a[2].re == 0.0) {
        a_re = 0.0;
        d = -a[2].im / a[1].re;
      } else {
        a_re = -a[2].re / a[1].re;
        d = -a[2].im / a[1].re;
      }
    } else if (a[1].re == 0.0) {
      if (-a[2].re == 0.0) {
        a_re = -a[2].im / a[1].im;
        d = 0.0;
      } else if (-a[2].im == 0.0) {
        a_re = 0.0;
        d = -(-a[2].re / a[1].im);
      } else {
        a_re = -a[2].im / a[1].im;
        d = -(-a[2].re / a[1].im);
      }
    } else if (b_brm > b_bim) {
      s = a[1].im / a[1].re;
      d = a[1].re + s * a[1].im;
      a_re = (-a[2].re + s * -a[2].im) / d;
      d = (-a[2].im - s * -a[2].re) / d;
    } else if (b_bim == b_brm) {
      if (a[1].re > 0.0) {
        sgnbr = 0.5;
      } else {
        sgnbr = -0.5;
      }
      if (a[1].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      a_re = (-a[2].re * sgnbr + -a[2].im * d) / b_brm;
      d = (-a[2].im * sgnbr - -a[2].re * d) / b_brm;
    } else {
      s = a[1].re / a[1].im;
      d = a[1].im + s * a[1].re;
      a_re = (s * -a[2].re + -a[2].im) / d;
      d = (s * -a[2].im - (-a[2].re)) / d;
    }
    c[2].re = a_re * t_re - d * t_im;
    c[2].im = a_re * t_im + d * t_re;
    c[3].re = r_re * t_re - r_im * t_im;
    c[3].im = r_re * t_im + r_im * t_re;
  } else {
    if (a[0].im == 0.0) {
      if (a[1].im == 0.0) {
        r_re = a[1].re / a[0].re;
        r_im = 0.0;
      } else if (a[1].re == 0.0) {
        r_re = 0.0;
        r_im = a[1].im / a[0].re;
      } else {
        r_re = a[1].re / a[0].re;
        r_im = a[1].im / a[0].re;
      }
    } else if (a[0].re == 0.0) {
      if (a[1].re == 0.0) {
        r_re = a[1].im / a[0].im;
        r_im = 0.0;
      } else if (a[1].im == 0.0) {
        r_re = 0.0;
        r_im = -(a[1].re / a[0].im);
      } else {
        r_re = a[1].im / a[0].im;
        r_im = -(a[1].re / a[0].im);
      }
    } else if (brm > bim) {
      s = a[0].im / a[0].re;
      d = a[0].re + s * a[0].im;
      r_re = (a[1].re + s * a[1].im) / d;
      r_im = (a[1].im - s * a[1].re) / d;
    } else if (bim == brm) {
      if (a[0].re > 0.0) {
        sgnbr = 0.5;
      } else {
        sgnbr = -0.5;
      }
      if (a[0].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      r_re = (a[1].re * sgnbr + a[1].im * d) / brm;
      r_im = (a[1].im * sgnbr - a[1].re * d) / brm;
    } else {
      s = a[0].re / a[0].im;
      d = a[0].im + s * a[0].re;
      r_re = (s * a[1].re + a[1].im) / d;
      r_im = (s * a[1].im - a[1].re) / d;
    }
    d = a[3].re - (r_re * a[2].re - r_im * a[2].im);
    sgnbr = a[3].im - (r_re * a[2].im + r_im * a[2].re);
    if (sgnbr == 0.0) {
      t_re = 1.0 / d;
      t_im = 0.0;
    } else if (d == 0.0) {
      t_re = 0.0;
      t_im = -(1.0 / sgnbr);
    } else {
      b_brm = fabs(d);
      b_bim = fabs(sgnbr);
      if (b_brm > b_bim) {
        s = sgnbr / d;
        d += s * sgnbr;
        t_re = 1.0 / d;
        t_im = (0.0 - s) / d;
      } else if (b_bim == b_brm) {
        if (d > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }
        t_re = d / b_brm;
        if (sgnbr > 0.0) {
          sgnbr = 0.5;
        } else {
          sgnbr = -0.5;
        }
        t_im = (0.0 - sgnbr) / b_brm;
      } else {
        s = d / sgnbr;
        d = sgnbr + s * d;
        t_re = s / d;
        t_im = -1.0 / d;
      }
    }
    if (a[0].im == 0.0) {
      if (a[3].im == 0.0) {
        a_re = a[3].re / a[0].re;
        d = 0.0;
      } else if (a[3].re == 0.0) {
        a_re = 0.0;
        d = a[3].im / a[0].re;
      } else {
        a_re = a[3].re / a[0].re;
        d = a[3].im / a[0].re;
      }
    } else if (a[0].re == 0.0) {
      if (a[3].re == 0.0) {
        a_re = a[3].im / a[0].im;
        d = 0.0;
      } else if (a[3].im == 0.0) {
        a_re = 0.0;
        d = -(a[3].re / a[0].im);
      } else {
        a_re = a[3].im / a[0].im;
        d = -(a[3].re / a[0].im);
      }
    } else if (brm > bim) {
      s = a[0].im / a[0].re;
      d = a[0].re + s * a[0].im;
      a_re = (a[3].re + s * a[3].im) / d;
      d = (a[3].im - s * a[3].re) / d;
    } else if (bim == brm) {
      if (a[0].re > 0.0) {
        sgnbr = 0.5;
      } else {
        sgnbr = -0.5;
      }
      if (a[0].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      a_re = (a[3].re * sgnbr + a[3].im * d) / brm;
      d = (a[3].im * sgnbr - a[3].re * d) / brm;
    } else {
      s = a[0].re / a[0].im;
      d = a[0].im + s * a[0].re;
      a_re = (s * a[3].re + a[3].im) / d;
      d = (s * a[3].im - a[3].re) / d;
    }
    c[0].re = a_re * t_re - d * t_im;
    c[0].im = a_re * t_im + d * t_re;
    c[1].re = -r_re * t_re - -r_im * t_im;
    c[1].im = -r_re * t_im + -r_im * t_re;
    if (a[0].im == 0.0) {
      if (-a[2].im == 0.0) {
        a_re = -a[2].re / a[0].re;
        d = 0.0;
      } else if (-a[2].re == 0.0) {
        a_re = 0.0;
        d = -a[2].im / a[0].re;
      } else {
        a_re = -a[2].re / a[0].re;
        d = -a[2].im / a[0].re;
      }
    } else if (a[0].re == 0.0) {
      if (-a[2].re == 0.0) {
        a_re = -a[2].im / a[0].im;
        d = 0.0;
      } else if (-a[2].im == 0.0) {
        a_re = 0.0;
        d = -(-a[2].re / a[0].im);
      } else {
        a_re = -a[2].im / a[0].im;
        d = -(-a[2].re / a[0].im);
      }
    } else if (brm > bim) {
      s = a[0].im / a[0].re;
      d = a[0].re + s * a[0].im;
      a_re = (-a[2].re + s * -a[2].im) / d;
      d = (-a[2].im - s * -a[2].re) / d;
    } else if (bim == brm) {
      if (a[0].re > 0.0) {
        sgnbr = 0.5;
      } else {
        sgnbr = -0.5;
      }
      if (a[0].im > 0.0) {
        d = 0.5;
      } else {
        d = -0.5;
      }
      a_re = (-a[2].re * sgnbr + -a[2].im * d) / brm;
      d = (-a[2].im * sgnbr - -a[2].re * d) / brm;
    } else {
      s = a[0].re / a[0].im;
      d = a[0].im + s * a[0].re;
      a_re = (s * -a[2].re + -a[2].im) / d;
      d = (s * -a[2].im - (-a[2].re)) / d;
    }
    c[2].re = a_re * t_re - d * t_im;
    c[2].im = a_re * t_im + d * t_re;
    c[3].re = t_re;
    c[3].im = t_im;
  }
}

/*
 *
 */
void mpower(const emxArray_real_T *a, emxArray_real_T *c)
{
  lapack_int info_t;
  lapack_int *ipiv_t_data;
  emxArray_int32_T *ipiv;
  emxArray_int32_T *p;
  emxArray_lapack_int *ipiv_t;
  emxArray_real_T *x;
  const double *a_data;
  double *c_data;
  double *x_data;
  int b_i;
  int b_n;
  int i;
  int k;
  int n;
  int yk;
  int *ipiv_data;
  int *p_data;
  a_data = a->data;
  if ((a->size[0] == 0) || (a->size[1] == 0)) {
    i = c->size[0] * c->size[1];
    c->size[0] = a->size[0];
    c->size[1] = a->size[1];
    emxEnsureCapacity_real_T(c, i);
    c_data = c->data;
    yk = a->size[0] * a->size[1];
    for (i = 0; i < yk; i++) {
      c_data[i] = a_data[i];
    }
  } else {
    n = a->size[0];
    i = c->size[0] * c->size[1];
    c->size[0] = a->size[0];
    c->size[1] = a->size[1];
    emxEnsureCapacity_real_T(c, i);
    c_data = c->data;
    yk = a->size[0] * a->size[1];
    for (i = 0; i < yk; i++) {
      c_data[i] = 0.0;
    }
    emxInit_real_T(&x, 2);
    i = x->size[0] * x->size[1];
    x->size[0] = a->size[0];
    x->size[1] = a->size[1];
    emxEnsureCapacity_real_T(x, i);
    x_data = x->data;
    yk = a->size[0] * a->size[1];
    for (i = 0; i < yk; i++) {
      x_data[i] = a_data[i];
    }
    emxInit_int32_T(&ipiv, 2);
    emxInit_lapack_int(&ipiv_t);
    i = ipiv_t->size[0];
    ipiv_t->size[0] = a->size[0];
    emxEnsureCapacity_lapack_int(ipiv_t, i);
    ipiv_t_data = ipiv_t->data;
    info_t = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, (lapack_int)a->size[0],
                                 (lapack_int)a->size[0], &x_data[0],
                                 (lapack_int)a->size[0], &ipiv_t_data[0]);
    i = ipiv->size[0] * ipiv->size[1];
    ipiv->size[0] = 1;
    ipiv->size[1] = ipiv_t->size[0];
    emxEnsureCapacity_int32_T(ipiv, i);
    ipiv_data = ipiv->data;
    if ((int)info_t < 0) {
      yk = x->size[0];
      b_n = x->size[1];
      i = x->size[0] * x->size[1];
      x->size[0] = yk;
      x->size[1] = b_n;
      emxEnsureCapacity_real_T(x, i);
      x_data = x->data;
      yk *= b_n;
      for (i = 0; i < yk; i++) {
        x_data[i] = 0.0;
      }
      i = ipiv_t->size[0] - 1;
      for (k = 0; k <= i; k++) {
        ipiv_data[k] = k + 1;
      }
    } else {
      i = ipiv_t->size[0] - 1;
      for (k = 0; k <= i; k++) {
        ipiv_data[k] = (int)ipiv_t_data[k];
      }
    }
    emxFree_lapack_int(&ipiv_t);
    emxInit_int32_T(&p, 2);
    b_n = a->size[0];
    i = p->size[0] * p->size[1];
    p->size[0] = 1;
    p->size[1] = a->size[0];
    emxEnsureCapacity_int32_T(p, i);
    p_data = p->data;
    p_data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      p_data[k - 1] = yk;
    }
    i = ipiv->size[1];
    for (k = 0; k < i; k++) {
      b_n = ipiv_data[k];
      if (b_n > k + 1) {
        yk = p_data[b_n - 1];
        p_data[b_n - 1] = p_data[k];
        p_data[k] = yk;
      }
    }
    emxFree_int32_T(&ipiv);
    for (k = 0; k < n; k++) {
      i = p_data[k];
      c_data[k + c->size[0] * (i - 1)] = 1.0;
      for (yk = k + 1; yk <= n; yk++) {
        if (c_data[(yk + c->size[0] * (i - 1)) - 1] != 0.0) {
          b_n = yk + 1;
          for (b_i = b_n; b_i <= n; b_i++) {
            c_data[(b_i + c->size[0] * (i - 1)) - 1] -=
                c_data[(yk + c->size[0] * (i - 1)) - 1] *
                x_data[(b_i + x->size[0] * (yk - 1)) - 1];
          }
        }
      }
    }
    emxFree_int32_T(&p);
    cblas_dtrsm(CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans,
                CblasNonUnit, (blasint)a->size[0], (blasint)a->size[0], 1.0,
                &x_data[0], (blasint)a->size[0], &c_data[0],
                (blasint)a->size[0]);
    emxFree_real_T(&x);
  }
}

/* End of code generation (mpower.c) */
