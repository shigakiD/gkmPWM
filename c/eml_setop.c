/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eml_setop.c
 *
 * Code generation for function 'eml_setop'
 *
 */

/* Include files */
#include "eml_setop.h"
#include "gkmPWMlasso3_emxutil.h"
#include "gkmPWMlasso3_types.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x);

/* Function Definitions */
/*
 *
 */
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x)
{
  const double *x_data;
  double absx;
  double xk;
  int exponent;
  bool exitg1;
  x_data = x->data;
  xk = x_data[*k - 1];
  exitg1 = false;
  while ((!exitg1) && (*k < x->size[0])) {
    absx = fabs(xk / 2.0);
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = ldexp(1.0, exponent - 53);
    }
    if (fabs(xk - x_data[*k]) < absx) {
      (*k)++;
    } else {
      exitg1 = true;
    }
  }
  return xk;
}

/*
 *
 */
void b_do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
                  emxArray_real_T *c, emxArray_int32_T *ia, int *ib_size)
{
  double absx;
  double ak;
  double bk;
  double *c_data;
  int b_ialast;
  int exponent;
  int iafirst;
  int ialast;
  int iblast;
  int na;
  int nc;
  int nia;
  int *ia_data;
  na = a->size[0];
  iblast = c->size[0];
  c->size[0] = a->size[0];
  emxEnsureCapacity_real_T(c, iblast);
  c_data = c->data;
  iblast = ia->size[0];
  ia->size[0] = a->size[0];
  emxEnsureCapacity_int32_T(ia, iblast);
  ia_data = ia->data;
  *ib_size = 0;
  nc = 0;
  nia = 0;
  iafirst = 0;
  ialast = 1;
  iblast = 1;
  while ((ialast <= na) && (iblast <= b->size[0])) {
    b_ialast = ialast;
    ak = skip_to_last_equal_value(&b_ialast, a);
    ialast = b_ialast;
    bk = skip_to_last_equal_value(&iblast, b);
    absx = fabs(bk / 2.0);
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = ldexp(1.0, exponent - 53);
    }
    if (fabs(bk - ak) < absx) {
      ialast = b_ialast + 1;
      iafirst = b_ialast;
      iblast++;
    } else if (ak < bk) {
      nc++;
      nia++;
      c_data[nc - 1] = ak;
      ia_data[nia - 1] = iafirst + 1;
      ialast = b_ialast + 1;
      iafirst = b_ialast;
    } else {
      iblast++;
    }
  }
  while (ialast <= na) {
    iblast = ialast;
    ak = skip_to_last_equal_value(&iblast, a);
    nc++;
    nia++;
    c_data[nc - 1] = ak;
    ia_data[nia - 1] = iafirst + 1;
    ialast = iblast + 1;
    iafirst = iblast;
  }
  if (a->size[0] > 0) {
    iblast = ia->size[0];
    if (1 > nia) {
      ia->size[0] = 0;
    } else {
      ia->size[0] = nia;
    }
    emxEnsureCapacity_int32_T(ia, iblast);
    iblast = c->size[0];
    if (1 > nc) {
      c->size[0] = 0;
    } else {
      c->size[0] = nc;
    }
    emxEnsureCapacity_real_T(c, iblast);
  }
}

/*
 *
 */
void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
                emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib)
{
  double absx;
  double ak;
  double bk;
  double *c_data;
  int b_ialast;
  int b_iblast;
  int exponent;
  int iafirst;
  int ialast;
  int ibfirst;
  int iblast;
  int nc;
  int ncmax;
  int *ia_data;
  int *ib_data;
  iafirst = a->size[0];
  ncmax = b->size[0];
  if (iafirst <= ncmax) {
    ncmax = iafirst;
  }
  iafirst = c->size[0];
  c->size[0] = ncmax;
  emxEnsureCapacity_real_T(c, iafirst);
  c_data = c->data;
  iafirst = ia->size[0];
  ia->size[0] = ncmax;
  emxEnsureCapacity_int32_T(ia, iafirst);
  ia_data = ia->data;
  iafirst = ib->size[0];
  ib->size[0] = ncmax;
  emxEnsureCapacity_int32_T(ib, iafirst);
  ib_data = ib->data;
  nc = 0;
  iafirst = 0;
  ialast = 1;
  ibfirst = 0;
  iblast = 1;
  while ((ialast <= a->size[0]) && (iblast <= b->size[0])) {
    b_ialast = ialast;
    ak = skip_to_last_equal_value(&b_ialast, a);
    ialast = b_ialast;
    b_iblast = iblast;
    bk = skip_to_last_equal_value(&b_iblast, b);
    iblast = b_iblast;
    absx = fabs(bk / 2.0);
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = ldexp(1.0, exponent - 53);
    }
    if (fabs(bk - ak) < absx) {
      nc++;
      c_data[nc - 1] = ak;
      ia_data[nc - 1] = iafirst + 1;
      ib_data[nc - 1] = ibfirst + 1;
      ialast = b_ialast + 1;
      iafirst = b_ialast;
      iblast = b_iblast + 1;
      ibfirst = b_iblast;
    } else if (ak < bk) {
      ialast = b_ialast + 1;
      iafirst = b_ialast;
    } else {
      iblast = b_iblast + 1;
      ibfirst = b_iblast;
    }
  }
  if (ncmax > 0) {
    iafirst = ia->size[0];
    if (1 > nc) {
      ia->size[0] = 0;
    } else {
      ia->size[0] = nc;
    }
    emxEnsureCapacity_int32_T(ia, iafirst);
    iafirst = ib->size[0];
    if (1 > nc) {
      ib->size[0] = 0;
    } else {
      ib->size[0] = nc;
    }
    emxEnsureCapacity_int32_T(ib, iafirst);
    iafirst = c->size[0];
    if (1 > nc) {
      c->size[0] = 0;
    } else {
      c->size[0] = nc;
    }
    emxEnsureCapacity_real_T(c, iafirst);
  }
}

/* End of code generation (eml_setop.c) */
