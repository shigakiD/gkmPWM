/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eig.c
 *
 * Code generation for function 'eig'
 *
 */

/* Include files */
#include "eig.h"
#include "gkmPWM_data.h"
#include "gkmPWM_rtwutil.h"
#include "schur.h"
#include "xzhgeqz.h"
#include "xzlartg.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 *
 */
void eig(const double A[16], creal_T V[4])
{
  creal_T At[16];
  creal_T beta1[4];
  creal_T atmp;
  double T[16];
  double b_A[16];
  double absxk;
  double anrm;
  double anrmto;
  double brm;
  double cfrom1;
  double cto1;
  double ctoc;
  double re;
  int At_tmp;
  int exitg1;
  int exitg3;
  int i;
  int ihi;
  int ii;
  int ilo;
  int j;
  int jcol;
  int jcolp1;
  int jrow;
  int nzcount;
  bool b_guard1 = false;
  bool exitg2;
  bool exitg4;
  bool guard1 = false;
  bool ilascl;
  bool notdone;
  notdone = true;
  j = 0;
  exitg2 = false;
  while ((!exitg2) && (j < 4)) {
    i = 0;
    do {
      exitg1 = 0;
      if (i <= j) {
        if (A[i + (j << 2)] != A[j + (i << 2)]) {
          notdone = false;
          exitg1 = 1;
        } else {
          i++;
        }
      } else {
        j++;
        exitg1 = 2;
      }
    } while (exitg1 == 0);
    if (exitg1 == 1) {
      exitg2 = true;
    }
  }
  if (notdone) {
    memcpy(&b_A[0], &A[0], 16U * sizeof(double));
    schur(b_A, T);
    for (ii = 0; ii < 4; ii++) {
      V[ii].re = T[ii + (ii << 2)];
      V[ii].im = 0.0;
    }
  } else {
    notdone = true;
    j = 0;
    exitg2 = false;
    while ((!exitg2) && (j < 4)) {
      i = 0;
      do {
        exitg1 = 0;
        if (i <= j) {
          if (A[i + (j << 2)] != -A[j + (i << 2)]) {
            notdone = false;
            exitg1 = 1;
          } else {
            i++;
          }
        } else {
          j++;
          exitg1 = 2;
        }
      } while (exitg1 == 0);
      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
    if (notdone) {
      memcpy(&b_A[0], &A[0], 16U * sizeof(double));
      schur(b_A, T);
      i = 1;
      do {
        exitg1 = 0;
        if (i <= 4) {
          b_guard1 = false;
          if (i != 4) {
            absxk = T[i + ((i - 1) << 2)];
            if (absxk != 0.0) {
              absxk = fabs(absxk);
              V[i - 1].re = 0.0;
              V[i - 1].im = absxk;
              V[i].re = 0.0;
              V[i].im = -absxk;
              i += 2;
            } else {
              b_guard1 = true;
            }
          } else {
            b_guard1 = true;
          }
          if (b_guard1) {
            V[i - 1].re = 0.0;
            V[i - 1].im = 0.0;
            i++;
          }
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    } else {
      anrm = 0.0;
      for (ii = 0; ii < 16; ii++) {
        re = A[ii];
        At[ii].re = re;
        At[ii].im = 0.0;
        absxk = rt_hypotd(re, 0.0);
        if (absxk > anrm) {
          anrm = absxk;
        }
      }
      ilascl = false;
      anrmto = anrm;
      guard1 = false;
      if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
        anrmto = 6.7178761075670888E-139;
        ilascl = true;
        guard1 = true;
      } else if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = true;
        guard1 = true;
      }
      if (guard1) {
        absxk = anrm;
        ctoc = anrmto;
        notdone = true;
        while (notdone) {
          cfrom1 = absxk * 2.0041683600089728E-292;
          cto1 = ctoc / 4.9896007738368E+291;
          if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
            brm = 2.0041683600089728E-292;
            absxk = cfrom1;
          } else if (cto1 > absxk) {
            brm = 4.9896007738368E+291;
            ctoc = cto1;
          } else {
            brm = ctoc / absxk;
            notdone = false;
          }
          for (ii = 0; ii < 16; ii++) {
            At[ii].re *= brm;
            At[ii].im *= brm;
          }
        }
      }
      ilo = 1;
      ihi = 4;
      do {
        exitg3 = 0;
        i = 0;
        j = 0;
        notdone = false;
        ii = ihi;
        exitg2 = false;
        while ((!exitg2) && (ii > 0)) {
          nzcount = 0;
          i = ii;
          j = ihi;
          jcol = 0;
          exitg4 = false;
          while ((!exitg4) && (jcol <= ihi - 1)) {
            At_tmp = (ii + (jcol << 2)) - 1;
            if ((At[At_tmp].re != 0.0) || (At[At_tmp].im != 0.0) ||
                (ii == jcol + 1)) {
              if (nzcount == 0) {
                j = jcol + 1;
                nzcount = 1;
                jcol++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              jcol++;
            }
          }
          if (nzcount < 2) {
            notdone = true;
            exitg2 = true;
          } else {
            ii--;
          }
        }
        if (!notdone) {
          exitg3 = 2;
        } else {
          if (i != ihi) {
            for (ii = 0; ii < 4; ii++) {
              nzcount = ii << 2;
              jcol = (i + nzcount) - 1;
              atmp = At[jcol];
              At_tmp = (ihi + nzcount) - 1;
              At[jcol] = At[At_tmp];
              At[At_tmp] = atmp;
            }
          }
          if (j != ihi) {
            for (ii = 0; ii < ihi; ii++) {
              nzcount = ii + ((j - 1) << 2);
              atmp = At[nzcount];
              At_tmp = ii + ((ihi - 1) << 2);
              At[nzcount] = At[At_tmp];
              At[At_tmp] = atmp;
            }
          }
          ihi--;
          if (ihi == 1) {
            exitg3 = 1;
          }
        }
      } while (exitg3 == 0);
      if (exitg3 != 1) {
        do {
          exitg1 = 0;
          i = 0;
          j = 0;
          notdone = false;
          jcol = ilo;
          exitg2 = false;
          while ((!exitg2) && (jcol <= ihi)) {
            nzcount = 0;
            i = ihi;
            j = jcol;
            ii = ilo;
            exitg4 = false;
            while ((!exitg4) && (ii <= ihi)) {
              At_tmp = (ii + ((jcol - 1) << 2)) - 1;
              if ((At[At_tmp].re != 0.0) || (At[At_tmp].im != 0.0) ||
                  (ii == jcol)) {
                if (nzcount == 0) {
                  i = ii;
                  nzcount = 1;
                  ii++;
                } else {
                  nzcount = 2;
                  exitg4 = true;
                }
              } else {
                ii++;
              }
            }
            if (nzcount < 2) {
              notdone = true;
              exitg2 = true;
            } else {
              jcol++;
            }
          }
          if (!notdone) {
            exitg1 = 1;
          } else {
            if (i != ilo) {
              for (ii = ilo; ii < 5; ii++) {
                nzcount = (ii - 1) << 2;
                jcol = (i + nzcount) - 1;
                atmp = At[jcol];
                At_tmp = (ilo + nzcount) - 1;
                At[jcol] = At[At_tmp];
                At[At_tmp] = atmp;
              }
            }
            if (j != ilo) {
              for (ii = 0; ii < ihi; ii++) {
                nzcount = ii + ((j - 1) << 2);
                atmp = At[nzcount];
                At_tmp = ii + ((ilo - 1) << 2);
                At[nzcount] = At[At_tmp];
                At[At_tmp] = atmp;
              }
            }
            ilo++;
            if (ilo == ihi) {
              exitg1 = 1;
            }
          }
        } while (exitg1 == 0);
      }
      if (ihi >= ilo + 2) {
        for (jcol = ilo - 1; jcol + 1 < ihi - 1; jcol++) {
          jcolp1 = jcol + 2;
          for (jrow = ihi - 1; jrow + 1 > jcol + 2; jrow--) {
            At_tmp = jrow + (jcol << 2);
            xzlartg(At[At_tmp - 1], At[At_tmp], &absxk, &atmp,
                    &At[(jrow + (jcol << 2)) - 1]);
            At[At_tmp].re = 0.0;
            At[At_tmp].im = 0.0;
            for (j = jcolp1; j < 5; j++) {
              nzcount = jrow + ((j - 1) << 2);
              ctoc = At[nzcount].im;
              cfrom1 = At[nzcount].re;
              cto1 = At[nzcount - 1].re;
              brm = At[nzcount - 1].im;
              At[nzcount].re =
                  absxk * cfrom1 - (atmp.re * cto1 + atmp.im * brm);
              At[nzcount].im = absxk * ctoc - (atmp.re * brm - atmp.im * cto1);
              At[nzcount - 1].re =
                  absxk * cto1 + (atmp.re * cfrom1 - atmp.im * ctoc);
              At[nzcount - 1].im =
                  absxk * brm + (atmp.re * ctoc + atmp.im * cfrom1);
            }
            atmp.re = -atmp.re;
            atmp.im = -atmp.im;
            for (i = 1; i <= ihi; i++) {
              nzcount = (i + ((jrow - 1) << 2)) - 1;
              ctoc = At[nzcount].im;
              cfrom1 = At[nzcount].re;
              ii = (i + (jrow << 2)) - 1;
              cto1 = At[ii].re;
              brm = At[ii].im;
              At[nzcount].re =
                  absxk * cfrom1 - (atmp.re * cto1 + atmp.im * brm);
              At[nzcount].im = absxk * ctoc - (atmp.re * brm - atmp.im * cto1);
              At[ii].re = absxk * cto1 + (atmp.re * cfrom1 - atmp.im * ctoc);
              At[ii].im = absxk * brm + (atmp.re * ctoc + atmp.im * cfrom1);
            }
          }
        }
      }
      xzhgeqz(At, ilo, ihi, &ii, V, beta1);
      if ((ii == 0) && ilascl) {
        notdone = true;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((cfrom1 > anrm) && (anrm != 0.0)) {
            brm = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (cto1 > anrmto) {
            brm = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            brm = anrm / anrmto;
            notdone = false;
          }
          for (ii = 0; ii < 4; ii++) {
            V[ii].re *= brm;
            V[ii].im *= brm;
          }
        }
      }
      for (ii = 0; ii < 4; ii++) {
        anrmto = V[ii].re;
        anrm = V[ii].im;
        cfrom1 = beta1[ii].re;
        cto1 = beta1[ii].im;
        if (cto1 == 0.0) {
          if (anrm == 0.0) {
            re = anrmto / cfrom1;
            absxk = 0.0;
          } else if (anrmto == 0.0) {
            re = 0.0;
            absxk = anrm / cfrom1;
          } else {
            re = anrmto / cfrom1;
            absxk = anrm / cfrom1;
          }
        } else if (cfrom1 == 0.0) {
          if (anrmto == 0.0) {
            re = anrm / cto1;
            absxk = 0.0;
          } else if (anrm == 0.0) {
            re = 0.0;
            absxk = -(anrmto / cto1);
          } else {
            re = anrm / cto1;
            absxk = -(anrmto / cto1);
          }
        } else {
          brm = fabs(cfrom1);
          absxk = fabs(cto1);
          if (brm > absxk) {
            ctoc = cto1 / cfrom1;
            absxk = cfrom1 + ctoc * cto1;
            re = (anrmto + ctoc * anrm) / absxk;
            absxk = (anrm - ctoc * anrmto) / absxk;
          } else if (absxk == brm) {
            if (cfrom1 > 0.0) {
              ctoc = 0.5;
            } else {
              ctoc = -0.5;
            }
            if (cto1 > 0.0) {
              absxk = 0.5;
            } else {
              absxk = -0.5;
            }
            re = (anrmto * ctoc + anrm * absxk) / brm;
            absxk = (anrm * ctoc - anrmto * absxk) / brm;
          } else {
            ctoc = cfrom1 / cto1;
            absxk = cto1 + ctoc * cfrom1;
            re = (ctoc * anrmto + anrm) / absxk;
            absxk = (ctoc * anrm - anrmto) / absxk;
          }
        }
        V[ii].re = re;
        V[ii].im = absxk;
      }
    }
  }
}

/* End of code generation (eig.c) */
