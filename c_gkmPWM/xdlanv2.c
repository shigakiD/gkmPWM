/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * xdlanv2.c
 *
 * Code generation for function 'xdlanv2'
 *
 */

/* Include files */
#include "xdlanv2.h"
#include "gkmPWM_rtwutil.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void xdlanv2(double *a, double *b, double *c, double *d, double *rt1r,
             double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn)
{
  double bcmax;
  double bcmis;
  double p;
  double scale;
  double tau;
  double temp;
  double z;
  int b_c;
  int count;
  if (*c == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (*b == 0.0) {
    *cs = 0.0;
    *sn = 1.0;
    temp = *d;
    *d = *a;
    *a = temp;
    *b = -*c;
    *c = 0.0;
  } else {
    temp = *a - *d;
    if ((temp == 0.0) && ((*b < 0.0) != (*c < 0.0))) {
      *cs = 1.0;
      *sn = 0.0;
    } else {
      p = 0.5 * temp;
      bcmis = fabs(*b);
      scale = fabs(*c);
      bcmax = fmax(bcmis, scale);
      if (*b >= 0.0) {
        count = 1;
      } else {
        count = -1;
      }
      if (*c >= 0.0) {
        b_c = 1;
      } else {
        b_c = -1;
      }
      bcmis = fmin(bcmis, scale) * (double)count * (double)b_c;
      scale = fmax(fabs(p), bcmax);
      z = p / scale * p + bcmax / scale * bcmis;
      if (z >= 8.8817841970012523E-16) {
        *a = sqrt(scale) * sqrt(z);
        if (p < 0.0) {
          *a = -*a;
        }
        z = p + *a;
        *a = *d + z;
        *d -= bcmax / z * bcmis;
        tau = rt_hypotd(*c, z);
        *cs = z / tau;
        *sn = *c / tau;
        *b -= *c;
        *c = 0.0;
      } else {
        bcmis = *b + *c;
        scale = fmax(fabs(temp), fabs(bcmis));
        count = 0;
        while ((scale >= 7.4428285367870146E+137) && (count <= 20)) {
          bcmis *= 1.3435752215134178E-138;
          temp *= 1.3435752215134178E-138;
          scale = fmax(fabs(temp), fabs(bcmis));
          count++;
        }
        while ((scale <= 1.3435752215134178E-138) && (count <= 20)) {
          bcmis *= 7.4428285367870146E+137;
          temp *= 7.4428285367870146E+137;
          scale = fmax(fabs(temp), fabs(bcmis));
          count++;
        }
        tau = rt_hypotd(bcmis, temp);
        *cs = sqrt(0.5 * (fabs(bcmis) / tau + 1.0));
        if (bcmis >= 0.0) {
          count = 1;
        } else {
          count = -1;
        }
        *sn = -(0.5 * temp / (tau * *cs)) * (double)count;
        bcmax = *a * *cs + *b * *sn;
        scale = -*a * *sn + *b * *cs;
        z = *c * *cs + *d * *sn;
        bcmis = -*c * *sn + *d * *cs;
        *b = scale * *cs + bcmis * *sn;
        *c = -bcmax * *sn + z * *cs;
        temp = 0.5 * ((bcmax * *cs + z * *sn) + (-scale * *sn + bcmis * *cs));
        *a = temp;
        *d = temp;
        if (*c != 0.0) {
          if (*b != 0.0) {
            if ((*b < 0.0) == (*c < 0.0)) {
              bcmis = sqrt(fabs(*b));
              scale = sqrt(fabs(*c));
              *a = bcmis * scale;
              if (*c >= 0.0) {
                p = *a;
              } else {
                p = -*a;
              }
              tau = 1.0 / sqrt(fabs(*b + *c));
              *a = temp + p;
              *d = temp - p;
              *b -= *c;
              *c = 0.0;
              bcmax = bcmis * tau;
              bcmis = scale * tau;
              temp = *cs * bcmax - *sn * bcmis;
              *sn = *cs * bcmis + *sn * bcmax;
              *cs = temp;
            }
          } else {
            *b = -*c;
            *c = 0.0;
            temp = *cs;
            *cs = -*sn;
            *sn = temp;
          }
        }
      }
    }
  }
  *rt1r = *a;
  *rt2r = *d;
  if (*c == 0.0) {
    *rt1i = 0.0;
    *rt2i = 0.0;
  } else {
    *rt1i = sqrt(fabs(*b)) * sqrt(fabs(*c));
    *rt2i = -*rt1i;
  }
}

/* End of code generation (xdlanv2.c) */
