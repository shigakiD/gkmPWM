/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xzhgeqz.c
 *
 * Code generation for function 'xzhgeqz'
 *
 */

/* Include files */
#include "xzhgeqz.h"
#include "gkmPWM_data.h"
#include "sqrt.h"
#include "xzlartg.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void xzhgeqz(const creal_T A[16], int ilo, int ihi, int *info,
             creal_T alpha1[4], creal_T beta1[4])
{
  creal_T b_A[16];
  creal_T ctemp;
  creal_T stemp;
  creal_T y;
  double ad22_im;
  double ad22_re;
  double anorm;
  double ascale;
  double b_atol;
  double colscale;
  double colssq;
  double eshift_im;
  double eshift_re;
  double scale;
  double ssq;
  double t;
  double x_im;
  int col;
  int exitg1;
  int i;
  int ifirst;
  int iiter;
  int ilast;
  int ilastm;
  int ilastm1;
  int istart;
  int j;
  int jiter;
  int jm1;
  int jp1;
  int y_re_tmp;
  bool b_guard1 = false;
  bool exitg2;
  bool failed;
  bool goto60;
  bool goto70;
  bool goto90;
  bool guard1 = false;
  bool guard2 = false;
  bool guard3 = false;
  memcpy(&b_A[0], &A[0], 16U * sizeof(creal_T));
  *info = 0;
  for (i = 0; i < 4; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }
  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  anorm = 0.0;
  if (ilo <= ihi) {
    scale = 3.3121686421112381E-170;
    ssq = 0.0;
    i = ihi - ilo;
    for (j = 0; j <= i; j++) {
      colscale = 3.3121686421112381E-170;
      colssq = 0.0;
      col = (ilo + j) - 1;
      jm1 = j + 1;
      if (jm1 > i) {
        jm1 = i;
      }
      y_re_tmp = ilo + jm1;
      for (jm1 = ilo; jm1 <= y_re_tmp; jm1++) {
        jp1 = (jm1 + (col << 2)) - 1;
        anorm = fabs(A[jp1].re);
        if (anorm > colscale) {
          t = colscale / anorm;
          colssq = colssq * t * t + 1.0;
          colscale = anorm;
        } else {
          t = anorm / colscale;
          colssq += t * t;
        }
        anorm = fabs(A[jp1].im);
        if (anorm > colscale) {
          t = colscale / anorm;
          colssq = colssq * t * t + 1.0;
          colscale = anorm;
        } else {
          t = anorm / colscale;
          colssq += t * t;
        }
      }
      if (scale >= colscale) {
        t = colscale / scale;
        ssq += t * t * colssq;
      } else {
        t = scale / colscale;
        ssq = colssq + t * t * ssq;
        scale = colscale;
      }
    }
    anorm = scale * sqrt(ssq);
  }
  t = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (t > 2.2250738585072014E-308) {
    b_atol = t;
  }
  t = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    t = anorm;
  }
  ascale = 1.0 / t;
  failed = true;
  y_re_tmp = ihi + 1;
  for (j = y_re_tmp; j < 5; j++) {
    alpha1[j - 1] = A[(j + ((j - 1) << 2)) - 1];
  }
  guard1 = false;
  guard2 = false;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    ilastm = ihi;
    iiter = 0;
    goto60 = false;
    goto70 = false;
    goto90 = false;
    jiter = 0;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1) - 1) {
        b_guard1 = false;
        if (ilast + 1 == ilo) {
          goto60 = true;
          b_guard1 = true;
        } else {
          y_re_tmp = ilastm1 << 2;
          jp1 = ilast + y_re_tmp;
          jm1 = ilast + (ilast << 2);
          y_re_tmp += ilastm1;
          if (fabs(b_A[jp1].re) + fabs(b_A[jp1].im) <=
              fmax(2.2250738585072014E-308,
                   2.2204460492503131E-16 *
                       ((fabs(b_A[jm1].re) + fabs(b_A[jm1].im)) +
                        (fabs(b_A[y_re_tmp].re) + fabs(b_A[y_re_tmp].im))))) {
            b_A[jp1].re = 0.0;
            b_A[jp1].im = 0.0;
            goto60 = true;
            b_guard1 = true;
          } else {
            j = ilastm1;
            guard3 = false;
            exitg2 = false;
            while ((!exitg2) && (j + 1 >= ilo)) {
              if (j + 1 == ilo) {
                guard3 = true;
                exitg2 = true;
              } else {
                y_re_tmp = j + ((j - 1) << 2);
                jp1 = j + (j << 2);
                if (fabs(b_A[y_re_tmp].re) + fabs(b_A[y_re_tmp].im) <=
                    fmax(2.2250738585072014E-308,
                         2.2204460492503131E-16 *
                             ((fabs(b_A[jp1].re) + fabs(b_A[jp1].im)) +
                              (fabs(b_A[y_re_tmp - 1].re) +
                               fabs(b_A[y_re_tmp - 1].im))))) {
                  b_A[y_re_tmp].re = 0.0;
                  b_A[y_re_tmp].im = 0.0;
                  guard3 = true;
                  exitg2 = true;
                } else {
                  j--;
                  guard3 = false;
                }
              }
            }
            if (guard3) {
              ifirst = j + 1;
              goto70 = true;
            }
            if (goto70) {
              b_guard1 = true;
            } else {
              memset(&alpha1[0], 0, 4U * sizeof(creal_T));
              memset(&beta1[0], 0, 4U * sizeof(creal_T));
              *info = 1;
              exitg1 = 1;
            }
          }
        }
        if (b_guard1) {
          if (goto60) {
            goto60 = false;
            alpha1[ilast] = b_A[ilast + (ilast << 2)];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = false;
              guard2 = true;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              ilastm = ilast + 1;
              jiter++;
            }
          } else {
            if (goto70) {
              goto70 = false;
              iiter++;
              if (iiter - iiter / 10 * 10 != 0) {
                jm1 = ilast + (ilast << 2);
                anorm = ascale * b_A[jm1].re;
                t = ascale * b_A[jm1].im;
                if (t == 0.0) {
                  ad22_re = anorm / 0.5;
                  ad22_im = 0.0;
                } else if (anorm == 0.0) {
                  ad22_re = 0.0;
                  ad22_im = t / 0.5;
                } else {
                  ad22_re = anorm / 0.5;
                  ad22_im = t / 0.5;
                }
                jm1 = ilastm1 + (ilast << 2);
                anorm = ascale * b_A[jm1].re;
                t = ascale * b_A[jm1].im;
                if (t == 0.0) {
                  stemp.re = anorm / 0.5;
                  stemp.im = 0.0;
                } else if (anorm == 0.0) {
                  stemp.re = 0.0;
                  stemp.im = t / 0.5;
                } else {
                  stemp.re = anorm / 0.5;
                  stemp.im = t / 0.5;
                }
                b_sqrt(&stemp);
                jm1 = ilast + (ilastm1 << 2);
                anorm = ascale * b_A[jm1].re;
                t = ascale * b_A[jm1].im;
                if (t == 0.0) {
                  y.re = anorm / 0.5;
                  y.im = 0.0;
                } else if (anorm == 0.0) {
                  y.re = 0.0;
                  y.im = t / 0.5;
                } else {
                  y.re = anorm / 0.5;
                  y.im = t / 0.5;
                }
                b_sqrt(&y);
                ctemp.re = stemp.re * y.re - stemp.im * y.im;
                ctemp.im = stemp.re * y.im + stemp.im * y.re;
                if ((ctemp.re != 0.0) || (ctemp.im != 0.0)) {
                  jm1 = ilastm1 + (ilastm1 << 2);
                  anorm = ascale * b_A[jm1].re;
                  t = ascale * b_A[jm1].im;
                  if (t == 0.0) {
                    anorm /= 0.5;
                    t = 0.0;
                  } else if (anorm == 0.0) {
                    anorm = 0.0;
                    t /= 0.5;
                  } else {
                    anorm /= 0.5;
                    t /= 0.5;
                  }
                  colssq = 0.5 * (anorm - ad22_re);
                  x_im = 0.5 * (t - ad22_im);
                  colscale = fabs(colssq) + fabs(x_im);
                  ssq = fmax(fabs(ctemp.re) + fabs(ctemp.im), colscale);
                  if (x_im == 0.0) {
                    stemp.re = colssq / ssq;
                    stemp.im = 0.0;
                  } else if (colssq == 0.0) {
                    stemp.re = 0.0;
                    stemp.im = x_im / ssq;
                  } else {
                    stemp.re = colssq / ssq;
                    stemp.im = x_im / ssq;
                  }
                  if (ctemp.im == 0.0) {
                    y.re = ctemp.re / ssq;
                    y.im = 0.0;
                  } else if (ctemp.re == 0.0) {
                    y.re = 0.0;
                    y.im = ctemp.im / ssq;
                  } else {
                    y.re = ctemp.re / ssq;
                    y.im = ctemp.im / ssq;
                  }
                  anorm = stemp.re * stemp.re - stemp.im * stemp.im;
                  t = stemp.re * stemp.im;
                  scale = y.re * y.im;
                  stemp.re = anorm + (y.re * y.re - y.im * y.im);
                  stemp.im = (t + t) + (scale + scale);
                  b_sqrt(&stemp);
                  y.re = ssq * stemp.re;
                  y.im = ssq * stemp.im;
                  if (colscale > 0.0) {
                    if (x_im == 0.0) {
                      t = colssq / colscale;
                      anorm = 0.0;
                    } else {
                      if (colssq == 0.0) {
                        t = 0.0;
                      } else {
                        t = colssq / colscale;
                      }
                      anorm = x_im / colscale;
                    }
                    if (t * y.re + anorm * y.im < 0.0) {
                      y.re = -y.re;
                      y.im = -y.im;
                    }
                  }
                  scale = colssq + y.re;
                  ssq = x_im + y.im;
                  if (ssq == 0.0) {
                    if (ctemp.im == 0.0) {
                      colssq = ctemp.re / scale;
                      anorm = 0.0;
                    } else if (ctemp.re == 0.0) {
                      colssq = 0.0;
                      anorm = ctemp.im / scale;
                    } else {
                      colssq = ctemp.re / scale;
                      anorm = ctemp.im / scale;
                    }
                  } else if (scale == 0.0) {
                    if (ctemp.re == 0.0) {
                      colssq = ctemp.im / ssq;
                      anorm = 0.0;
                    } else if (ctemp.im == 0.0) {
                      colssq = 0.0;
                      anorm = -(ctemp.re / ssq);
                    } else {
                      colssq = ctemp.im / ssq;
                      anorm = -(ctemp.re / ssq);
                    }
                  } else {
                    colscale = fabs(scale);
                    anorm = fabs(ssq);
                    if (colscale > anorm) {
                      t = ssq / scale;
                      anorm = scale + t * ssq;
                      colssq = (ctemp.re + t * ctemp.im) / anorm;
                      anorm = (ctemp.im - t * ctemp.re) / anorm;
                    } else if (anorm == colscale) {
                      if (scale > 0.0) {
                        t = 0.5;
                      } else {
                        t = -0.5;
                      }
                      if (ssq > 0.0) {
                        anorm = 0.5;
                      } else {
                        anorm = -0.5;
                      }
                      colssq = (ctemp.re * t + ctemp.im * anorm) / colscale;
                      anorm = (ctemp.im * t - ctemp.re * anorm) / colscale;
                    } else {
                      t = scale / ssq;
                      anorm = ssq + t * scale;
                      colssq = (t * ctemp.re + ctemp.im) / anorm;
                      anorm = (t * ctemp.im - ctemp.re) / anorm;
                    }
                  }
                  ad22_re -= ctemp.re * colssq - ctemp.im * anorm;
                  ad22_im -= ctemp.re * anorm + ctemp.im * colssq;
                }
              } else {
                if (iiter - iiter / 20 * 20 == 0) {
                  jm1 = ilast + (ilast << 2);
                  anorm = ascale * b_A[jm1].re;
                  t = ascale * b_A[jm1].im;
                  if (t == 0.0) {
                    anorm /= 0.5;
                    t = 0.0;
                  } else if (anorm == 0.0) {
                    anorm = 0.0;
                    t /= 0.5;
                  } else {
                    anorm /= 0.5;
                    t /= 0.5;
                  }
                  eshift_re += anorm;
                  eshift_im += t;
                } else {
                  jm1 = ilast + (ilastm1 << 2);
                  anorm = ascale * b_A[jm1].re;
                  t = ascale * b_A[jm1].im;
                  if (t == 0.0) {
                    anorm /= 0.5;
                    t = 0.0;
                  } else if (anorm == 0.0) {
                    anorm = 0.0;
                    t /= 0.5;
                  } else {
                    anorm /= 0.5;
                    t /= 0.5;
                  }
                  eshift_re += anorm;
                  eshift_im += t;
                }
                ad22_re = eshift_re;
                ad22_im = eshift_im;
              }
              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = false;
              while ((!exitg2) && (j + 1 > ifirst)) {
                istart = j + 1;
                i = j << 2;
                col = j + i;
                ctemp.re = ascale * b_A[col].re - ad22_re * 0.5;
                ctemp.im = ascale * b_A[col].im - ad22_im * 0.5;
                ssq = fabs(ctemp.re) + fabs(ctemp.im);
                jm1 = jp1 + i;
                colscale = ascale * (fabs(b_A[jm1].re) + fabs(b_A[jm1].im));
                anorm = ssq;
                if (colscale > ssq) {
                  anorm = colscale;
                }
                if ((anorm < 1.0) && (anorm != 0.0)) {
                  ssq /= anorm;
                  colscale /= anorm;
                }
                y_re_tmp = j + ((j - 1) << 2);
                if ((fabs(b_A[y_re_tmp].re) + fabs(b_A[y_re_tmp].im)) *
                        colscale <=
                    ssq * b_atol) {
                  goto90 = true;
                  exitg2 = true;
                } else {
                  jp1 = j;
                  j--;
                }
              }
              if (!goto90) {
                istart = ifirst;
                col = (ifirst + ((ifirst - 1) << 2)) - 1;
                ctemp.re = ascale * b_A[col].re - ad22_re * 0.5;
                ctemp.im = ascale * b_A[col].im - ad22_im * 0.5;
              }
              goto90 = false;
              jm1 = istart + ((istart - 1) << 2);
              stemp.re = ascale * b_A[jm1].re;
              stemp.im = ascale * b_A[jm1].im;
              b_xzlartg(ctemp, stemp, &anorm, &y);
              j = istart;
              jm1 = istart - 2;
              while (j < ilast + 1) {
                if (j > istart) {
                  jp1 = j + (jm1 << 2);
                  xzlartg(b_A[jp1 - 1], b_A[jp1], &anorm, &y,
                          &b_A[(j + (jm1 << 2)) - 1]);
                  b_A[jp1].re = 0.0;
                  b_A[jp1].im = 0.0;
                }
                for (jp1 = j; jp1 <= ilastm; jp1++) {
                  y_re_tmp = j + ((jp1 - 1) << 2);
                  t = b_A[y_re_tmp].im;
                  scale = b_A[y_re_tmp].re;
                  ssq = b_A[y_re_tmp - 1].re;
                  stemp.re = anorm * ssq + (y.re * scale - y.im * t);
                  colscale = b_A[y_re_tmp - 1].im;
                  stemp.im = anorm * colscale + (y.re * t + y.im * scale);
                  b_A[y_re_tmp].re =
                      anorm * scale - (y.re * ssq + y.im * colscale);
                  b_A[y_re_tmp].im = anorm * t - (y.re * colscale - y.im * ssq);
                  b_A[y_re_tmp - 1] = stemp;
                }
                y.re = -y.re;
                y.im = -y.im;
                col = j;
                if (ilast + 1 < j + 2) {
                  col = ilast - 1;
                }
                for (i = ifirst; i <= col + 2; i++) {
                  y_re_tmp = (i + ((j - 1) << 2)) - 1;
                  t = b_A[y_re_tmp].im;
                  scale = b_A[y_re_tmp].re;
                  jm1 = (i + (j << 2)) - 1;
                  ssq = b_A[jm1].re;
                  stemp.re = anorm * ssq + (y.re * scale - y.im * t);
                  colscale = b_A[jm1].im;
                  stemp.im = anorm * colscale + (y.re * t + y.im * scale);
                  b_A[y_re_tmp].re =
                      anorm * scale - (y.re * ssq + y.im * colscale);
                  b_A[y_re_tmp].im = anorm * t - (y.re * colscale - y.im * ssq);
                  b_A[jm1] = stemp;
                }
                jm1 = j - 1;
                j++;
              }
            }
            jiter++;
          }
        }
      } else {
        guard2 = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    guard1 = true;
  }
  if (guard2) {
    if (failed) {
      *info = ilast + 1;
      if (0 <= ilast) {
        memset(&alpha1[0], 0, (ilast + 1) * sizeof(creal_T));
        memset(&beta1[0], 0, (ilast + 1) * sizeof(creal_T));
      }
    } else {
      guard1 = true;
    }
  }
  if (guard1) {
    for (j = 0; j <= ilo - 2; j++) {
      alpha1[j] = b_A[j + (j << 2)];
    }
  }
}

/* End of code generation (xzhgeqz.c) */
