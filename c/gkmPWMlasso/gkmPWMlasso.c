/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso.c
 *
 * Code generation for function 'gkmPWMlasso'
 *
 */

/* Include files */
#include "gkmPWMlasso.h"
#include "BGkmer.h"
#include "PWM2kmers.h"
#include "PWM2kmers_norc.h"
#include "blockedSummation.h"
#include "conncomp.h"
#include "corrcoef.h"
#include "diag.h"
#include "eml_setop.h"
#include "extremeKElements.h"
#include "fgetl.h"
#include "fileManager.h"
#include "fileread.h"
#include "find.h"
#include "genIndex.h"
#include "getgkmcounts.h"
#include "getmotif.h"
#include "gkmPWMlasso_data.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_initialize.h"
#include "gkmPWMlasso_rtwutil.h"
#include "gkmPWMlasso_types.h"
#include "graph.h"
#include "lasso_cvmat.h"
#include "minOrMax.h"
#include "mod.h"
#include "movSumProdOrMean.h"
#include "mpower.h"
#include "mtimes.h"
#include "nullAssignment.h"
#include "repmat.h"
#include "sort.h"
#include "sprintf.h"
#include "std.h"
#include "sum.h"
#include "unique.h"
#include "useConstantDim.h"
#include "cblas.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Function Declarations */
static void b_clus_simmat_eig(const emxArray_real_T *simmat,
                              emxArray_cell_wrap_1 *motclus);

static void binary_expand_op(emxArray_real_T *A, int j,
                             const emxArray_real_T *negvec, double b);

static double c_binary_expand_op(
    const emxArray_real_T *B, const emxArray_real_T *cfile2,
    const emxArray_real_T *loc, double varargin_3, emxArray_real_T *weigmat2,
    double expl_temp_data[], int expl_temp_size[2], double b_expl_temp_data[],
    int b_expl_temp_size[2], double a__3_DF_data[], int a__3_DF_size[2],
    double c_expl_temp_data[], int c_expl_temp_size[2]);

static void clus_simmat_eig(const emxArray_real_T *simmat, double b_r,
                            emxArray_cell_wrap_1 *motclus);

static void d_binary_expand_op(emxArray_real_T *cfile2,
                               const emxArray_real_T *negvec, double y,
                               double b);

static void gettopmotifs(const emxArray_real_T *weigvec,
                         const emxArray_real_T *pvec, const emxArray_real_T *E,
                         const emxArray_cell_wrap_1 *motclus,
                         const emxArray_char_T *filename_Value,
                         const emxArray_char_T *memefile, double n, double minL,
                         double minInfo, double C);

static void trim_pwm(emxArray_cell_wrap_0 *p, emxArray_real_T *info,
                     emxArray_real_T *len);

static void u_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *mat,
                               const emxArray_real_T *r1);

/* Function Definitions */
/*
 * function motclus = clus_simmat_eig(simmat,r)
 */
static void b_clus_simmat_eig(const emxArray_real_T *simmat,
                              emxArray_cell_wrap_1 *motclus)
{
  cell_wrap_1 *motclus_data;
  emxArray_boolean_T *b_bin;
  emxArray_boolean_T *b_simmat;
  emxArray_int32_T *b_r;
  emxArray_int32_T *t0_Underlying_Ir;
  emxArray_int32_T *t0_Underlying_Jc;
  emxArray_real_T *bin;
  const double *simmat_data;
  double d;
  double n;
  double *bin_data;
  int b_i;
  int i;
  int k;
  int last;
  int *r1;
  bool *b_simmat_data;
  simmat_data = simmat->data;
  emxInit_boolean_T(&b_simmat, 2);
  /* 'gkmPWMlasso:450' bin = conncomp(graph(simmat>r)); */
  k = b_simmat->size[0] * b_simmat->size[1];
  b_simmat->size[0] = simmat->size[0];
  b_simmat->size[1] = simmat->size[1];
  emxEnsureCapacity_boolean_T(b_simmat, k);
  b_simmat_data = b_simmat->data;
  last = simmat->size[0] * simmat->size[1];
  for (k = 0; k < last; k++) {
    b_simmat_data[k] = (simmat_data[k] > 0.6);
  }
  emxInit_real_T(&bin, 2);
  emxInit_int32_T(&t0_Underlying_Ir, 2);
  emxInit_int32_T(&t0_Underlying_Jc, 1);
  graph_graph(b_simmat, t0_Underlying_Ir, t0_Underlying_Jc);
  graph_conncomp(t0_Underlying_Ir, t0_Underlying_Jc, bin);
  bin_data = bin->data;
  /* 'gkmPWMlasso:451' n = max(bin); */
  last = bin->size[1];
  emxFree_boolean_T(&b_simmat);
  emxFree_int32_T(&t0_Underlying_Jc);
  emxFree_int32_T(&t0_Underlying_Ir);
  if (bin->size[1] <= 2) {
    if (bin->size[1] == 1) {
      n = bin_data[0];
    } else if (bin_data[0] < bin_data[bin->size[1] - 1]) {
      n = bin_data[bin->size[1] - 1];
    } else {
      n = bin_data[0];
    }
  } else {
    n = bin_data[0];
    for (k = 2; k <= last; k++) {
      d = bin_data[k - 1];
      if (n < d) {
        n = d;
      }
    }
  }
  /* 'gkmPWMlasso:452' motclus = cell(n,1); */
  /* 'gkmPWMlasso:453' for i = 1:n */
  k = (int)n;
  i = motclus->size[0];
  motclus->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  emxInit_int32_T(&b_r, 2);
  emxInit_boolean_T(&b_bin, 2);
  for (b_i = 0; b_i < k; b_i++) {
    /* 'gkmPWMlasso:454' motclus{i} = find(bin == i); */
    i = b_bin->size[0] * b_bin->size[1];
    b_bin->size[0] = 1;
    b_bin->size[1] = bin->size[1];
    emxEnsureCapacity_boolean_T(b_bin, i);
    b_simmat_data = b_bin->data;
    last = bin->size[1];
    for (i = 0; i < last; i++) {
      b_simmat_data[i] = (bin_data[i] == (double)b_i + 1.0);
    }
    eml_find(b_bin, b_r);
    r1 = b_r->data;
    i = motclus_data[b_i].f1->size[0] * motclus_data[b_i].f1->size[1];
    motclus_data[b_i].f1->size[0] = 1;
    motclus_data[b_i].f1->size[1] = b_r->size[1];
    emxEnsureCapacity_real_T(motclus_data[b_i].f1, i);
    last = b_r->size[1];
    for (i = 0; i < last; i++) {
      motclus_data[b_i].f1->data[i] = r1[i];
    }
  }
  emxFree_boolean_T(&b_bin);
  emxFree_int32_T(&b_r);
  emxFree_real_T(&bin);
}

static void binary_expand_op(emxArray_real_T *A, int j,
                             const emxArray_real_T *negvec, double b)
{
  emxArray_real_T *b_A;
  const double *negvec_data;
  double *A_data;
  double *b_A_data;
  int c_A;
  int i;
  int stride_0_0;
  int stride_1_0;
  negvec_data = negvec->data;
  A_data = A->data;
  emxInit_real_T(&b_A, 1);
  c_A = A->size[0] - 1;
  i = b_A->size[0];
  if (negvec->size[0] == 1) {
    b_A->size[0] = c_A + 1;
  } else {
    b_A->size[0] = negvec->size[0];
  }
  emxEnsureCapacity_real_T(b_A, i);
  b_A_data = b_A->data;
  stride_0_0 = (c_A + 1 != 1);
  stride_1_0 = (negvec->size[0] != 1);
  if (negvec->size[0] == 1) {
    c_A++;
  } else {
    c_A = negvec->size[0];
  }
  for (i = 0; i < c_A; i++) {
    b_A_data[i] = A_data[i * stride_0_0 + A->size[0] * j] -
                  negvec_data[i * stride_1_0] * b;
  }
  c_A = b_A->size[0];
  for (i = 0; i < c_A; i++) {
    A_data[i + A->size[0] * j] = b_A_data[i];
  }
  emxFree_real_T(&b_A);
}

static double c_binary_expand_op(
    const emxArray_real_T *B, const emxArray_real_T *cfile2,
    const emxArray_real_T *loc, double varargin_3, emxArray_real_T *weigmat2,
    double expl_temp_data[], int expl_temp_size[2], double b_expl_temp_data[],
    int b_expl_temp_size[2], double a__3_DF_data[], int a__3_DF_size[2],
    double c_expl_temp_data[], int c_expl_temp_size[2])
{
  emxArray_real_T *b_B;
  emxArray_real_T *b_cfile2;
  const double *B_data;
  const double *cfile2_data;
  const double *loc_data;
  double expl_temp;
  double *b_cfile2_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  loc_data = loc->data;
  cfile2_data = cfile2->data;
  B_data = B->data;
  emxInit_real_T(&b_cfile2, 1);
  i = b_cfile2->size[0];
  if (loc->size[0] == 1) {
    b_cfile2->size[0] = cfile2->size[0];
  } else {
    b_cfile2->size[0] = loc->size[0];
  }
  emxEnsureCapacity_real_T(b_cfile2, i);
  b_cfile2_data = b_cfile2->data;
  stride_0_0 = (cfile2->size[0] != 1);
  stride_1_0 = (loc->size[0] != 1);
  if (loc->size[0] == 1) {
    loop_ub = cfile2->size[0];
  } else {
    loop_ub = loc->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_cfile2_data[i] = cfile2_data[i * stride_0_0] - loc_data[i * stride_1_0];
  }
  emxInit_real_T(&b_B, 2);
  i = b_B->size[0] * b_B->size[1];
  b_B->size[0] = B->size[0];
  b_B->size[1] = B->size[1];
  emxEnsureCapacity_real_T(b_B, i);
  b_cfile2_data = b_B->data;
  loop_ub = B->size[0] * B->size[1] - 1;
  for (i = 0; i <= loop_ub; i++) {
    b_cfile2_data[i] = B_data[i];
  }
  b_lasso_cvmat(b_B, b_cfile2, varargin_3, weigmat2, expl_temp_data,
                expl_temp_size, b_expl_temp_data, b_expl_temp_size, &expl_temp,
                a__3_DF_data, a__3_DF_size, c_expl_temp_data, c_expl_temp_size);
  emxFree_real_T(&b_B);
  emxFree_real_T(&b_cfile2);
  return expl_temp;
}

/*
 * function motclus = clus_simmat_eig(simmat,r)
 */
static void clus_simmat_eig(const emxArray_real_T *simmat, double b_r,
                            emxArray_cell_wrap_1 *motclus)
{
  cell_wrap_1 *motclus_data;
  emxArray_boolean_T *b_bin;
  emxArray_boolean_T *b_simmat;
  emxArray_int32_T *c_r;
  emxArray_int32_T *t1_Underlying_Ir;
  emxArray_int32_T *t1_Underlying_Jc;
  emxArray_real_T *bin;
  const double *simmat_data;
  double d;
  double n;
  double *bin_data;
  int b_i;
  int i;
  int k;
  int last;
  int *r1;
  bool *b_simmat_data;
  simmat_data = simmat->data;
  emxInit_boolean_T(&b_simmat, 2);
  /* 'gkmPWMlasso:450' bin = conncomp(graph(simmat>r)); */
  k = b_simmat->size[0] * b_simmat->size[1];
  b_simmat->size[0] = simmat->size[0];
  b_simmat->size[1] = simmat->size[1];
  emxEnsureCapacity_boolean_T(b_simmat, k);
  b_simmat_data = b_simmat->data;
  last = simmat->size[0] * simmat->size[1];
  for (k = 0; k < last; k++) {
    b_simmat_data[k] = (simmat_data[k] > b_r);
  }
  emxInit_real_T(&bin, 2);
  emxInit_int32_T(&t1_Underlying_Ir, 2);
  emxInit_int32_T(&t1_Underlying_Jc, 1);
  graph_graph(b_simmat, t1_Underlying_Ir, t1_Underlying_Jc);
  graph_conncomp(t1_Underlying_Ir, t1_Underlying_Jc, bin);
  bin_data = bin->data;
  /* 'gkmPWMlasso:451' n = max(bin); */
  last = bin->size[1];
  emxFree_boolean_T(&b_simmat);
  emxFree_int32_T(&t1_Underlying_Jc);
  emxFree_int32_T(&t1_Underlying_Ir);
  if (bin->size[1] <= 2) {
    if (bin->size[1] == 1) {
      n = bin_data[0];
    } else if (bin_data[0] < bin_data[bin->size[1] - 1]) {
      n = bin_data[bin->size[1] - 1];
    } else {
      n = bin_data[0];
    }
  } else {
    n = bin_data[0];
    for (k = 2; k <= last; k++) {
      d = bin_data[k - 1];
      if (n < d) {
        n = d;
      }
    }
  }
  /* 'gkmPWMlasso:452' motclus = cell(n,1); */
  /* 'gkmPWMlasso:453' for i = 1:n */
  k = (int)n;
  i = motclus->size[0];
  motclus->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  emxInit_int32_T(&c_r, 2);
  emxInit_boolean_T(&b_bin, 2);
  for (b_i = 0; b_i < k; b_i++) {
    /* 'gkmPWMlasso:454' motclus{i} = find(bin == i); */
    i = b_bin->size[0] * b_bin->size[1];
    b_bin->size[0] = 1;
    b_bin->size[1] = bin->size[1];
    emxEnsureCapacity_boolean_T(b_bin, i);
    b_simmat_data = b_bin->data;
    last = bin->size[1];
    for (i = 0; i < last; i++) {
      b_simmat_data[i] = (bin_data[i] == (double)b_i + 1.0);
    }
    eml_find(b_bin, c_r);
    r1 = c_r->data;
    i = motclus_data[b_i].f1->size[0] * motclus_data[b_i].f1->size[1];
    motclus_data[b_i].f1->size[0] = 1;
    motclus_data[b_i].f1->size[1] = c_r->size[1];
    emxEnsureCapacity_real_T(motclus_data[b_i].f1, i);
    last = c_r->size[1];
    for (i = 0; i < last; i++) {
      motclus_data[b_i].f1->data[i] = r1[i];
    }
  }
  emxFree_boolean_T(&b_bin);
  emxFree_int32_T(&c_r);
  emxFree_real_T(&bin);
}

static void d_binary_expand_op(emxArray_real_T *cfile2,
                               const emxArray_real_T *negvec, double y,
                               double b)
{
  emxArray_real_T *b_cfile2;
  const double *negvec_data;
  double *b_cfile2_data;
  double *cfile2_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  negvec_data = negvec->data;
  cfile2_data = cfile2->data;
  emxInit_real_T(&b_cfile2, 1);
  i = b_cfile2->size[0];
  if (negvec->size[0] == 1) {
    b_cfile2->size[0] = cfile2->size[0];
  } else {
    b_cfile2->size[0] = negvec->size[0];
  }
  emxEnsureCapacity_real_T(b_cfile2, i);
  b_cfile2_data = b_cfile2->data;
  stride_0_0 = (cfile2->size[0] != 1);
  stride_1_0 = (negvec->size[0] != 1);
  if (negvec->size[0] == 1) {
    loop_ub = cfile2->size[0];
  } else {
    loop_ub = negvec->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_cfile2_data[i] =
        cfile2_data[i * stride_0_0] - negvec_data[i * stride_1_0] / y * b;
  }
  i = cfile2->size[0];
  cfile2->size[0] = b_cfile2->size[0];
  emxEnsureCapacity_real_T(cfile2, i);
  cfile2_data = cfile2->data;
  loop_ub = b_cfile2->size[0];
  for (i = 0; i < loop_ub; i++) {
    cfile2_data[i] = b_cfile2_data[i];
  }
  emxFree_real_T(&b_cfile2);
}

/*
 * function gettopmotifs(weigvec, pvec, E, motclus, filename, memefile, n, minL
 * , minInfo, C)
 */
static void gettopmotifs(const emxArray_real_T *weigvec,
                         const emxArray_real_T *pvec, const emxArray_real_T *E,
                         const emxArray_cell_wrap_1 *motclus,
                         const emxArray_char_T *filename_Value,
                         const emxArray_char_T *memefile, double n, double minL,
                         double minInfo, double C)
{
  static const char cv[6] = {'\x09', '\x0a', '\x0b', '\x0c', '\x0d', ' '};
  static const char b[5] = {'M', 'O', 'T', 'I', 'F'};
  FILE *b_NULL;
  FILE *c_NULL;
  FILE *filestar;
  const cell_wrap_1 *motclus_data;
  cell_wrap_1 *ordered_motclus_data;
  cell_wrap_8 *names_data;
  emxArray_cell_wrap_1 *ordered_motclus;
  emxArray_cell_wrap_8 *names;
  emxArray_char_T *charStr;
  emxArray_char_T *line;
  emxArray_char_T *varargin_1;
  emxArray_int32_T *iidx;
  emxArray_real_T *a;
  emxArray_real_T *c;
  const double *E_data;
  const double *pvec_data;
  const double *weigvec_data;
  double b_i;
  double *a_data;
  double *c_data;
  int exitg1;
  int exitg2;
  int exitg4;
  int exitg5;
  int exitg6;
  int i;
  int i1;
  int i2;
  int k;
  int l;
  int nbytes;
  int st;
  int *iidx_data;
  char c_a[5];
  const char *filename_Value_data;
  char b_c;
  signed char fileid;
  char *charStr_data;
  char *line_data;
  bool b_a;
  bool exitg7;
  filename_Value_data = filename_Value->data;
  motclus_data = motclus->data;
  E_data = E->data;
  pvec_data = pvec->data;
  weigvec_data = weigvec->data;
  emxInit_real_T(&a, 1);
  /* 'gkmPWMlasso:406' [a, b] = sort(pvec, 'descend'); */
  i = a->size[0];
  a->size[0] = pvec->size[0];
  emxEnsureCapacity_real_T(a, i);
  a_data = a->data;
  nbytes = pvec->size[0];
  for (i = 0; i < nbytes; i++) {
    a_data[i] = pvec_data[i];
  }
  emxInit_real_T(&c, 2);
  emxInit_int32_T(&iidx, 1);
  sort(a, iidx);
  iidx_data = iidx->data;
  a_data = a->data;
  /* 'gkmPWMlasso:407' c = [weigvec(b) a E(b)]; */
  i = c->size[0] * c->size[1];
  c->size[0] = iidx->size[0];
  c->size[1] = 3;
  emxEnsureCapacity_real_T(c, i);
  c_data = c->data;
  nbytes = iidx->size[0];
  for (i = 0; i < nbytes; i++) {
    c_data[i] = weigvec_data[iidx_data[i] - 1];
  }
  nbytes = a->size[0];
  for (i = 0; i < nbytes; i++) {
    c_data[i + c->size[0]] = a_data[i];
  }
  emxFree_real_T(&a);
  nbytes = iidx->size[0];
  for (i = 0; i < nbytes; i++) {
    c_data[i + c->size[0] * 2] = E_data[iidx_data[i] - 1];
  }
  emxInit_cell_wrap_1(&ordered_motclus);
  /*  motclus = motclus(b); */
  /* 'gkmPWMlasso:409' mylen = length(b); */
  /* 'gkmPWMlasso:410' ordered_motclus = cell(mylen, 1); */
  /* 'gkmPWMlasso:411' for idx=1:mylen */
  i = iidx->size[0];
  i1 = ordered_motclus->size[0];
  ordered_motclus->size[0] = iidx->size[0];
  emxEnsureCapacity_cell_wrap_1(ordered_motclus, i1);
  ordered_motclus_data = ordered_motclus->data;
  for (k = 0; k < i; k++) {
    /* 'gkmPWMlasso:412' ordered_motclus{idx} = motclus{b(idx)}; */
    i1 = ordered_motclus_data[k].f1->size[0] *
         ordered_motclus_data[k].f1->size[1];
    ordered_motclus_data[k].f1->size[0] = 1;
    ordered_motclus_data[k].f1->size[1] =
        motclus_data[iidx_data[k] - 1].f1->size[1];
    emxEnsureCapacity_real_T(ordered_motclus_data[k].f1, i1);
    nbytes = motclus_data[iidx_data[k] - 1].f1->size[1];
    for (i1 = 0; i1 < nbytes; i1++) {
      ordered_motclus_data[k].f1->data[i1] =
          motclus_data[iidx_data[k] - 1].f1->data[i1];
    }
  }
  emxFree_int32_T(&iidx);
  emxInit_cell_wrap_8(&names);
  /* 'gkmPWMlasso:414' motclus = ordered_motclus; */
  /* 'gkmPWMlasso:415' [rr,cc] = size(c); */
  /* 'gkmPWMlasso:416' names = cell(n,1); */
  nbytes = (int)n;
  i = names->size[0];
  names->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_8(names, i);
  names_data = names->data;
  for (i = 0; i < nbytes; i++) {
    names_data[i].f1->size[0] = 1;
    names_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso:417' names = coder.nullcopy(names); */
  /* 'gkmPWMlasso:418' i = 0; */
  b_i = 0.0;
  /* 'gkmPWMlasso:419' fid = fopen(memefile, 'r'); */
  fileid = cfopen(memefile, "rb");
  /* 'gkmPWMlasso:420' if fid == -1 */
  if (fileid == -1) {
    /* 'gkmPWMlasso:421' fprintf("ERROR: Cannot open motif database.\n"); */
    printf("ERROR: Cannot open motif database.\n");
    fflush(stdout);
    exit(1);
  }
  /* 'gkmPWMlasso:423' while ~feof(fid) */
  b_NULL = NULL;
  emxInit_char_T(&line, 2);
  do {
    exitg1 = 0;
    getfilestar(fileid, &filestar, &b_a);
    if (filestar == b_NULL) {
      i = 0;
    } else {
      st = feof(filestar);
      i = ((int)st != 0);
    }
    if (i == 0) {
      /* 'gkmPWMlasso:424' line = fgetl(fid); */
      fgetl(fileid, line);
      line_data = line->data;
      /* 'gkmPWMlasso:425' if length(line) >= 5 */
      if (line->size[1] >= 5) {
        /* 'gkmPWMlasso:426' if strcmp(line(1:5), 'MOTIF') */
        for (i = 0; i < 5; i++) {
          c_a[i] = line_data[i];
        }
        nbytes = memcmp(&c_a[0], &b[0], 5);
        if (nbytes == 0) {
          /* 'gkmPWMlasso:427' i = i+1; */
          b_i++;
          /* 'gkmPWMlasso:428' [~,name] = strtok(line); */
          nbytes = line->size[1];
          k = 1;
          do {
            exitg6 = 0;
            if (k <= nbytes) {
              l = 0;
              do {
                exitg5 = 0;
                if (l < 6) {
                  if (line_data[k - 1] == cv[l]) {
                    k++;
                    exitg5 = 1;
                  } else {
                    l++;
                  }
                } else {
                  exitg5 = 2;
                }
              } while (exitg5 == 0);
              if (exitg5 != 1) {
                exitg6 = 1;
              }
            } else {
              exitg6 = 2;
            }
          } while (exitg6 == 0);
          do {
            exitg4 = 0;
            if (k <= nbytes) {
              l = 0;
              do {
                exitg2 = 0;
                if (l < 6) {
                  if (line_data[k - 1] == cv[l]) {
                    exitg2 = 1;
                  } else {
                    l++;
                  }
                } else {
                  k++;
                  exitg2 = 2;
                }
              } while (exitg2 == 0);
              if (exitg2 == 1) {
                exitg4 = 1;
              }
            } else {
              exitg4 = 1;
            }
          } while (exitg4 == 0);
          if (k > line->size[1]) {
            i = -1;
            i1 = 0;
          } else {
            i = k - 2;
            i1 = line->size[1];
          }
          /* 'gkmPWMlasso:429' names{i} = strtrim(name); */
          nbytes = 1;
          exitg7 = false;
          while ((!exitg7) && (nbytes <= (i1 - i) - 1)) {
            b_c = line_data[i + nbytes];
            if ((!bv[(unsigned char)b_c & 127]) || (b_c == '\x00')) {
              exitg7 = true;
            } else {
              nbytes++;
            }
          }
          k = (i1 - i) - 1;
          exitg7 = false;
          while ((!exitg7) && (k > 0)) {
            b_c = line_data[i + k];
            if ((!bv[(unsigned char)b_c & 127]) || (b_c == '\x00')) {
              exitg7 = true;
            } else {
              k--;
            }
          }
          if (nbytes > k) {
            i1 = 0;
            k = 0;
          } else {
            i1 = nbytes - 1;
          }
          i2 = names_data[(int)b_i - 1].f1->size[0] *
               names_data[(int)b_i - 1].f1->size[1];
          names_data[(int)b_i - 1].f1->size[0] = 1;
          nbytes = k - i1;
          names_data[(int)b_i - 1].f1->size[1] = nbytes;
          emxEnsureCapacity_char_T(names_data[(int)b_i - 1].f1, i2);
          for (i2 = 0; i2 < nbytes; i2++) {
            names_data[(int)b_i - 1].f1->data[i2] =
                line_data[((i + i1) + i2) + 1];
          }
        }
      }
    } else {
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  /* 'gkmPWMlasso:433' fclose(fid); */
  cfclose(fileid);
  /* 'gkmPWMlasso:434' fidw = fopen(sprintf('%s_gkmPWMlasso.out', filename),
   * 'w'); */
  i = line->size[0] * line->size[1];
  line->size[0] = 1;
  line->size[1] = filename_Value->size[1] + 1;
  emxEnsureCapacity_char_T(line, i);
  line_data = line->data;
  nbytes = filename_Value->size[1];
  for (i = 0; i < nbytes; i++) {
    line_data[i] = filename_Value_data[i];
  }
  emxInit_char_T(&varargin_1, 2);
  line_data[filename_Value->size[1]] = '\x00';
  i = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = 1;
  varargin_1->size[1] = filename_Value->size[1] + 1;
  emxEnsureCapacity_char_T(varargin_1, i);
  charStr_data = varargin_1->data;
  nbytes = filename_Value->size[1];
  for (i = 0; i < nbytes; i++) {
    charStr_data[i] = filename_Value_data[i];
  }
  emxInit_char_T(&charStr, 2);
  charStr_data[filename_Value->size[1]] = '\x00';
  nbytes = snprintf(NULL, 0, "%s_gkmPWMlasso.out", &charStr_data[0]);
  i = charStr->size[0] * charStr->size[1];
  charStr->size[0] = 1;
  charStr->size[1] = nbytes + 1;
  emxEnsureCapacity_char_T(charStr, i);
  charStr_data = charStr->data;
  snprintf(&charStr_data[0], (size_t)(nbytes + 1), "%s_gkmPWMlasso.out",
           &line_data[0]);
  i = charStr->size[0] * charStr->size[1];
  if (1 > nbytes) {
    charStr->size[1] = 0;
  } else {
    charStr->size[1] = nbytes;
  }
  emxEnsureCapacity_char_T(charStr, i);
  fileid = cfopen(charStr, "wb");
  /* 'gkmPWMlasso:435' if fidw == -1 */
  emxFree_char_T(&varargin_1);
  emxFree_char_T(&charStr);
  if (fileid == -1) {
    /* 'gkmPWMlasso:436' fprintf("ERROR: Cannot create gkmPWMlasso output
     * file.\n"); */
    printf("ERROR: Cannot create gkmPWMlasso output file.\n");
    fflush(stdout);
    exit(1);
  }
  /* 'gkmPWMlasso:438' fprintf(fidw, 'Minimum PWM Length:\t%d\n', int32(minL));
   */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    b_i = rt_roundd(minL);
    if (b_i < 2.147483648E+9) {
      if (b_i >= -2.147483648E+9) {
        i = (int)b_i;
      } else {
        i = MIN_int32_T;
      }
    } else {
      i = MAX_int32_T;
    }
    fprintf(filestar, "Minimum PWM Length:\t%d\n", i);
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso:439' fprintf(fidw, 'Minimum PWM Information:\t%f\n', minInfo);
   */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Minimum PWM Information:\t%f\n", minInfo);
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso:440' fprintf(fidw, 'Correlation with SVM weights:\t%f\n', C);
   */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Correlation with SVM weights:\t%f\n", C);
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso:441' fprintf(fidw, 'Cluster ID\tMotif ID\tMotif
   * Name\tRegression Weight\tZ-score\tImportance\n'); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Cluster ID\tMotif ID\tMotif Name\tRegression "
                      "Weight\tZ-score\tImportance\n");
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso:442' for j = 1:length(E) */
  i = E->size[0];
  for (k = 0; k < i; k++) {
    /* 'gkmPWMlasso:443' for l = 1:length(motclus{j}) */
    i1 = ordered_motclus_data[k].f1->size[1];
    if (0 <= ordered_motclus_data[k].f1->size[1] - 1) {
      c_NULL = NULL;
    }
    for (l = 0; l < i1; l++) {
      /* 'gkmPWMlasso:444' fprintf(fidw, '%d\t%d\t%s\t%0.3f\t%0.3f\t%0.3f\n',
       * int32(j), int32(motclus{j}(l)), names{motclus{j}(l)}, c(j,1), c(j,2),
       * c(j,3)); */
      i2 = line->size[0] * line->size[1];
      line->size[0] = 1;
      line->size[1] =
          names_data[(int)ordered_motclus_data[k].f1->data[l] - 1].f1->size[1] +
          1;
      emxEnsureCapacity_char_T(line, i2);
      line_data = line->data;
      nbytes =
          names_data[(int)ordered_motclus_data[k].f1->data[l] - 1].f1->size[1];
      for (i2 = 0; i2 < nbytes; i2++) {
        line_data[i2] = names_data[(int)ordered_motclus_data[k].f1->data[l] - 1]
                            .f1->data[i2];
      }
      line_data[names_data[(int)ordered_motclus_data[k].f1->data[l] - 1]
                    .f1->size[1]] = '\x00';
      getfilestar(fileid, &filestar, &b_a);
      if (!(filestar == c_NULL)) {
        b_i = rt_roundd(ordered_motclus_data[k].f1->data[l]);
        if (b_i < 2.147483648E+9) {
          if (b_i >= -2.147483648E+9) {
            i2 = (int)b_i;
          } else {
            i2 = MIN_int32_T;
          }
        } else {
          i2 = MAX_int32_T;
        }
        fprintf(filestar, "%d\t%d\t%s\t%0.3f\t%0.3f\t%0.3f\n", k + 1, i2,
                &line_data[0], c_data[k], c_data[k + c->size[0]],
                c_data[k + c->size[0] * 2]);
        if (b_a) {
          fflush(filestar);
        }
      }
    }
  }
  emxFree_char_T(&line);
  emxFree_cell_wrap_8(&names);
  emxFree_cell_wrap_1(&ordered_motclus);
  emxFree_real_T(&c);
  /* 'gkmPWMlasso:447' fclose(fidw); */
  cfclose(fileid);
}

/*
 * function [pp, info, len] = trim_pwm(p,cut)
 */
static void trim_pwm(emxArray_cell_wrap_0 *p, emxArray_real_T *info,
                     emxArray_real_T *len)
{
  cell_wrap_0 *p_data;
  emxArray_real_T *b_mat;
  emxArray_real_T *b_r;
  emxArray_real_T *mat;
  emxArray_real_T *mvec;
  emxArray_real_T *vec;
  double *b_mat_data;
  double *info_data;
  double *len_data;
  double *mat_data;
  double *mvec_data;
  int b_i;
  int c_i;
  int i;
  int idx;
  int j;
  int nrows;
  int nx;
  int nxout;
  p_data = p->data;
  /* 'gkmPWMlasso:458' l = length(p); */
  /* 'gkmPWMlasso:459' info = zeros(l, 1); */
  nx = p->size[0];
  i = info->size[0];
  info->size[0] = nx;
  emxEnsureCapacity_real_T(info, i);
  info_data = info->data;
  /* 'gkmPWMlasso:460' len = zeros(l,1); */
  i = len->size[0];
  len->size[0] = nx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;
  /* 'gkmPWMlasso:461' for i = 1:l */
  i = p->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&vec, 1);
  emxInit_real_T(&mvec, 1);
  emxInit_real_T(&b_r, 2);
  emxInit_real_T(&b_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWMlasso:462' mat = p{i}+(p{i}==0); */
    idx = mat->size[0] * mat->size[1];
    mat->size[0] = p_data[b_i].f1->size[0];
    mat->size[1] = p_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(mat, idx);
    mat_data = mat->data;
    nx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
    for (idx = 0; idx < nx; idx++) {
      mat_data[idx] = p_data[b_i].f1->data[idx] +
                      (double)(p_data[b_i].f1->data[idx] == 0.0);
    }
    /* 'gkmPWMlasso:463' vec = 2+sum(mat.*log(mat)/log(2),2); */
    idx = b_r->size[0] * b_r->size[1];
    b_r->size[0] = mat->size[0];
    b_r->size[1] = mat->size[1];
    emxEnsureCapacity_real_T(b_r, idx);
    mvec_data = b_r->data;
    nx = mat->size[0] * mat->size[1];
    for (idx = 0; idx < nx; idx++) {
      mvec_data[idx] = mat_data[idx];
    }
    nx = mat->size[0] * mat->size[1];
    for (nrows = 0; nrows < nx; nrows++) {
      mvec_data[nrows] = log(mvec_data[nrows]);
    }
    if ((mat->size[0] == b_r->size[0]) && (mat->size[1] == b_r->size[1])) {
      idx = b_mat->size[0] * b_mat->size[1];
      b_mat->size[0] = mat->size[0];
      b_mat->size[1] = mat->size[1];
      emxEnsureCapacity_real_T(b_mat, idx);
      b_mat_data = b_mat->data;
      nx = mat->size[0] * mat->size[1];
      for (idx = 0; idx < nx; idx++) {
        b_mat_data[idx] = mat_data[idx] * mvec_data[idx] / 0.69314718055994529;
      }
      b_sum(b_mat, vec);
      mat_data = vec->data;
    } else {
      u_binary_expand_op(vec, mat, b_r);
      mat_data = vec->data;
    }
    nx = vec->size[0];
    for (idx = 0; idx < nx; idx++) {
      mat_data[idx] += 2.0;
    }
    /* 'gkmPWMlasso:464' mvec = movmean(vec,5); */
    nx = vec->size[0];
    if (2 <= nx) {
      nx = 2;
    }
    vmovfun(vec, vec->size[0], nx, nx, mvec);
    mvec_data = mvec->data;
    /* 'gkmPWMlasso:465' while min(vec(1),mvec(1)) < cut && length(vec) > 1 */
    while ((fmin(mat_data[0], mvec_data[0]) < 0.0) && (vec->size[0] > 1)) {
      /* 'gkmPWMlasso:466' p{i}(1,:) = []; */
      nx = p_data[b_i].f1->size[0] - 2;
      nxout = p_data[b_i].f1->size[1];
      nrows = p_data[b_i].f1->size[0] - 1;
      for (j = 0; j < nxout; j++) {
        for (c_i = 0; c_i < nrows; c_i++) {
          p_data[b_i].f1->data[c_i + p_data[b_i].f1->size[0] * j] =
              p_data[b_i].f1->data[(c_i + p_data[b_i].f1->size[0] * j) + 1];
        }
      }
      if (1 > nrows) {
        nx = 0;
      } else {
        nx++;
      }
      nxout = p_data[b_i].f1->size[1] - 1;
      for (idx = 0; idx <= nxout; idx++) {
        for (nrows = 0; nrows < nx; nrows++) {
          p_data[b_i].f1->data[nrows + nx * idx] =
              p_data[b_i].f1->data[nrows + p_data[b_i].f1->size[0] * idx];
        }
      }
      idx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
      p_data[b_i].f1->size[0] = nx;
      p_data[b_i].f1->size[1] = nxout + 1;
      emxEnsureCapacity_real_T(p_data[b_i].f1, idx);
      /* 'gkmPWMlasso:467' vec(1) = []; */
      nx = vec->size[0];
      nxout = vec->size[0];
      for (nrows = 0; nrows <= nxout - 2; nrows++) {
        mat_data[nrows] = mat_data[nrows + 1];
      }
      idx = vec->size[0];
      vec->size[0] = nx - 1;
      emxEnsureCapacity_real_T(vec, idx);
      mat_data = vec->data;
      /* 'gkmPWMlasso:468' mvec(1)=[]; */
      nx = mvec->size[0];
      nxout = mvec->size[0] - 1;
      for (nrows = 0; nrows < nxout; nrows++) {
        mvec_data[nrows] = mvec_data[nrows + 1];
      }
      idx = mvec->size[0];
      if (1 > nxout) {
        mvec->size[0] = 0;
      } else {
        mvec->size[0] = nx - 1;
      }
      emxEnsureCapacity_real_T(mvec, idx);
      mvec_data = mvec->data;
      /* 'gkmPWMlasso:469' mat(1,:) = []; */
      idx = mat->size[0] * mat->size[1];
      if (1 > mat->size[0] - 1) {
        mat->size[0] = 0;
      } else {
        mat->size[0]--;
      }
      emxEnsureCapacity_real_T(mat, idx);
    }
    /* 'gkmPWMlasso:471' while min(mvec(end),vec(end)) < cut && length(vec) > 1
     */
    while ((fmin(mvec_data[mvec->size[0] - 1], mat_data[vec->size[0] - 1]) <
            0.0) &&
           (vec->size[0] > 1)) {
      /* 'gkmPWMlasso:472' vec(end) = []; */
      idx = vec->size[0];
      nx = vec->size[0];
      nxout = vec->size[0] - 1;
      for (nrows = idx; nrows <= nxout; nrows++) {
        mat_data[nrows - 1] = mat_data[nrows];
      }
      idx = vec->size[0];
      vec->size[0] = nx - 1;
      emxEnsureCapacity_real_T(vec, idx);
      mat_data = vec->data;
      /* 'gkmPWMlasso:473' mvec(end)=[]; */
      idx = mvec->size[0];
      nx = mvec->size[0];
      nxout = mvec->size[0] - 1;
      for (nrows = idx; nrows <= nxout; nrows++) {
        mvec_data[nrows - 1] = mvec_data[nrows];
      }
      idx = mvec->size[0];
      if (1 > nxout) {
        mvec->size[0] = 0;
      } else {
        mvec->size[0] = nx - 1;
      }
      emxEnsureCapacity_real_T(mvec, idx);
      mvec_data = mvec->data;
      /* 'gkmPWMlasso:474' mat(end,:) = []; */
      idx = mat->size[0] * mat->size[1];
      if (1 > mat->size[0] - 1) {
        mat->size[0] = 0;
      } else {
        mat->size[0]--;
      }
      emxEnsureCapacity_real_T(mat, idx);
      /* 'gkmPWMlasso:475' p{i}(end,:) = []; */
      idx = p_data[b_i].f1->size[0];
      nx = p_data[b_i].f1->size[0] - 2;
      nxout = p_data[b_i].f1->size[1];
      nrows = p_data[b_i].f1->size[0] - 1;
      for (j = 0; j < nxout; j++) {
        for (c_i = idx; c_i <= nrows; c_i++) {
          p_data[b_i].f1->data[(c_i + p_data[b_i].f1->size[0] * j) - 1] =
              p_data[b_i].f1->data[c_i + p_data[b_i].f1->size[0] * j];
        }
      }
      if (1 > nrows) {
        nx = 0;
      } else {
        nx++;
      }
      nxout = p_data[b_i].f1->size[1] - 1;
      for (idx = 0; idx <= nxout; idx++) {
        for (nrows = 0; nrows < nx; nrows++) {
          p_data[b_i].f1->data[nrows + nx * idx] =
              p_data[b_i].f1->data[nrows + p_data[b_i].f1->size[0] * idx];
        }
      }
      idx = p_data[b_i].f1->size[0] * p_data[b_i].f1->size[1];
      p_data[b_i].f1->size[0] = nx;
      p_data[b_i].f1->size[1] = nxout + 1;
      emxEnsureCapacity_real_T(p_data[b_i].f1, idx);
    }
    /* 'gkmPWMlasso:477' info(i) = sum(vec); */
    info_data[b_i] = blockedSummation(vec, vec->size[0]);
    /* 'gkmPWMlasso:478' [len(i), ~] = size(mat); */
    len_data[b_i] = mat->size[0];
  }
  emxFree_real_T(&b_mat);
  emxFree_real_T(&b_r);
  emxFree_real_T(&mvec);
  emxFree_real_T(&vec);
  emxFree_real_T(&mat);
  /* 'gkmPWMlasso:480' pp = p; */
}

static void u_binary_expand_op(emxArray_real_T *vec, const emxArray_real_T *mat,
                               const emxArray_real_T *r1)
{
  emxArray_real_T *b_mat;
  const double *b_r;
  const double *mat_data;
  double *b_mat_data;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  b_r = r1->data;
  mat_data = mat->data;
  emxInit_real_T(&b_mat, 2);
  i = b_mat->size[0] * b_mat->size[1];
  if (r1->size[0] == 1) {
    b_mat->size[0] = mat->size[0];
  } else {
    b_mat->size[0] = r1->size[0];
  }
  if (r1->size[1] == 1) {
    b_mat->size[1] = mat->size[1];
  } else {
    b_mat->size[1] = r1->size[1];
  }
  emxEnsureCapacity_real_T(b_mat, i);
  b_mat_data = b_mat->data;
  stride_0_0 = (mat->size[0] != 1);
  stride_0_1 = (mat->size[1] != 1);
  stride_1_0 = (r1->size[0] != 1);
  stride_1_1 = (r1->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  if (r1->size[1] == 1) {
    loop_ub = mat->size[1];
  } else {
    loop_ub = r1->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    if (r1->size[0] == 1) {
      b_loop_ub = mat->size[0];
    } else {
      b_loop_ub = r1->size[0];
    }
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_mat_data[i1 + b_mat->size[0] * i] =
          mat_data[i1 * stride_0_0 + mat->size[0] * aux_0_1] *
          b_r[i1 * stride_1_0 + r1->size[0] * aux_1_1] / 0.69314718055994529;
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  b_sum(b_mat, vec);
  emxFree_real_T(&b_mat);
}

/*
 * function gkmPWMlasso(varargin)
 */
void gkmPWMlasso(const emxArray_char_T *varargin_1,
                 const emxArray_char_T *varargin_2, double varargin_3,
                 double varargin_4, double varargin_5, double varargin_6,
                 double varargin_7, double varargin_8, double varargin_9,
                 double varargin_10, double varargin_11, double varargin_12)
{
  static const char cv[5] = {'M', 'O', 'T', 'I', 'F'};
  cell_wrap_0 *p_data;
  cell_wrap_1 *b_motclus_data;
  cell_wrap_1 *motclus_data;
  cell_wrap_1 *tmp_motclus_data;
  cell_wrap_3 *newFF_data;
  emxArray_boolean_T b_MSE2_data;
  emxArray_boolean_T *b_loc;
  emxArray_cell_wrap_0 *p;
  emxArray_cell_wrap_1 *c_motclus;
  emxArray_cell_wrap_1 *motclus;
  emxArray_cell_wrap_1 *tmp_motclus;
  emxArray_cell_wrap_3 *b_newFF;
  emxArray_cell_wrap_3_1x19 FF;
  emxArray_cell_wrap_3_1x19 b_FF;
  emxArray_cell_wrap_3_1x19 newFF;
  emxArray_cell_wrap_3_20 F;
  emxArray_char_T *text;
  emxArray_int32_T *ib;
  emxArray_int32_T *idx;
  emxArray_int32_T *match_out;
  emxArray_int32_T *matches;
  emxArray_real_T b_B;
  emxArray_real_T *A;
  emxArray_real_T *B;
  emxArray_real_T *BB;
  emxArray_real_T *BY;
  emxArray_real_T *GCmat;
  emxArray_real_T *OLS;
  emxArray_real_T *Pweig;
  emxArray_real_T *b_A;
  emxArray_real_T *b_f;
  emxArray_real_T *b_y;
  emxArray_real_T *c_f;
  emxArray_real_T *cfile2;
  emxArray_real_T *comb;
  emxArray_real_T *comb2;
  emxArray_real_T *diffc;
  emxArray_real_T *f;
  emxArray_real_T *indc;
  emxArray_real_T *indvec;
  emxArray_real_T *loc;
  emxArray_real_T *negvec;
  emxArray_real_T *normvec;
  emxArray_real_T *weigmat;
  emxArray_real_T *xc;
  emxArray_real_T *y;
  double a__3_DF_data[100];
  double b_expl_temp_data[100];
  double c_expl_temp_data[100];
  double d_expl_temp_data[100];
  double e_expl_temp_data[100];
  double expl_temp_data[100];
  double f_expl_temp_data[100];
  double MSE_data[20];
  double MSE2_data[19];
  double mat[16];
  double mat2[16];
  double dv[4];
  double lk_data[2];
  double GCpos1;
  double a;
  double d;
  double lcnum;
  double nfrac;
  double rcnum;
  double *A_data;
  double *BY_data;
  double *B_data;
  double *GCmat_data;
  double *OLS_data;
  double *Pweig_data;
  double *b_f_data;
  double *cfile2_data;
  double *comb_data;
  double *f_data;
  double *indc_data;
  double *indvec_data;
  double *loc_data;
  double *negvec_data;
  double *normvec_data;
  int iidx_data[19];
  int MSE2_size[2];
  int a__3_DF_size[2];
  int b_MSE2_size[2];
  int b_expl_temp_size[2];
  int expl_temp_size[2];
  int lk_size[2];
  int b_i;
  int b_loop_ub;
  int b_motclus;
  int b_nxout;
  int c_loop_ub;
  int d_loop_ub;
  int e_loop_ub;
  int f_loop_ub;
  int g_loop_ub;
  int h_loop_ub;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i_loop_ub;
  int j;
  int j_loop_ub;
  int k;
  int k_loop_ub;
  int l_loop_ub;
  int loop_ub;
  int m_loop_ub;
  int match_idx;
  int n;
  int n_loop_ub;
  int nx;
  int nxout;
  int o_loop_ub;
  int size_tmp_idx_1;
  int *match_out_data;
  int *matches_data;
  char *text_data;
  bool c_MSE2_data[19];
  bool empty_non_axis_sizes;
  bool exitg1;
  bool *b_loc_data;
  if (!isInitialized_gkmPWMlasso) {
    gkmPWMlasso_initialize();
  }
  /*  gkmPWMlasso finds predictive motifs from sequence-based models of
   * regulatory */
  /*      elements using a database of PWMs.   */
  /*   */
  /*      gkmPWMlasso(fileprefix, memefile, motifNum, ...) */
  /*   */
  /*      Works conveniently with the output of the gkmSVM R package or the
   * lsgkm */
  /*      package (https://github.com/Dongwon-Lee/lsgkm).  Can be leveraged to
   */
  /*      extract motifs from other sequence based models as long as you have
   * scores  */
  /*      paired with sequences.  The sequences should be in fasta format with
   * the  */
  /*      suffix *_svseq.fa.  The scores should be in a 2 column tab delimited
   */
  /*      file with the suffix *svalpha.out with the same prefix as the
   * *svseq.fa  */
  /*      file.  The first column containing the labels of each sequence and the
   */
  /*      scores in the second.   */
  /*   */
  /*      Uses a modified version of MATLAB's lasso function.  The changes made
   */
  /*      were to reduce the computation time of working with large matrices. */
  /* 'gkmPWMlasso:19' filename = varargin{1}; */
  /* 'gkmPWMlasso:20' memefile = varargin{2}; */
  /* 'gkmPWMlasso:21' d = varargin{3}; */
  /* 'gkmPWMlasso:22' minL = varargin{4}; */
  /* 'gkmPWMlasso:23' minInfo = varargin{5}; */
  /* 'gkmPWMlasso:24' corrCut = varargin{6}; */
  /* 'gkmPWMlasso:25' l_svm = varargin{7}; */
  /* 'gkmPWMlasso:26' k_svm = varargin{8}; */
  /* 'gkmPWMlasso:27' BG_GC = varargin{9}; */
  /* 'gkmPWMlasso:28' RC = varargin{10}; */
  /* 'gkmPWMlasso:29' nfrac = varargin{11}; */
  /* 'gkmPWMlasso:30' nfracLim = varargin{12}; */
  /* 'gkmPWMlasso:32' lk = 1; */
  lk_size[0] = 1;
  lk_size[1] = 1;
  lk_data[0] = 1.0;
  /* 'gkmPWMlasso:33' if nfrac ~= 1 */
  if (varargin_11 != 1.0) {
    /* 'gkmPWMlasso:34' lk = [l_svm k_svm]; */
    lk_size[0] = 1;
    lk_size[1] = 2;
    lk_data[0] = varargin_7;
    lk_data[1] = varargin_8;
  }
  emxInit_real_T(&comb, 2);
  emxInit_real_T(&comb2, 2);
  emxInit_real_T(&diffc, 1);
  emxInit_real_T(&indc, 1);
  emxInit_real_T(&xc, 2);
  /* 'gkmPWMlasso:38' [comb,comb2,diffc,indc,xc,rcnum] =
   * genIndex(l_svm,k_svm,nfrac); */
  genIndex(varargin_7, varargin_8, varargin_11, comb, comb2, diffc, indc, xc,
           &rcnum);
  /* generate gapped positions, adjusted for reverse complements */
  /* 'gkmPWMlasso:40' if nfracLim && numel(comb)/k_svm*4^k_svm > 5*10^5 */
  if (varargin_12 != 0.0) {
    d = pow(4.0, varargin_8);
    if ((double)(comb->size[0] * comb->size[1]) / varargin_8 * d > 500000.0) {
      /* 'gkmPWMlasso:41' nfrac = round(5*10^7/4^k_svm/numel(comb)*k_svm)/100;
       */
      nfrac = rt_roundd(5.0E+7 / d / (double)(comb->size[0] * comb->size[1]) *
                        varargin_8) /
              100.0;
      /* 'gkmPWMlasso:43' l_svm2 = l_svm; */
      /* 'gkmPWMlasso:44' k_svm2 = k_svm; */
      /* 'gkmPWMlasso:45' lk = ([l_svm k_svm]); */
      lk_size[0] = 1;
      lk_size[1] = 2;
      lk_data[0] = varargin_7;
      lk_data[1] = varargin_8;
      /* 'gkmPWMlasso:46' [comb,comb2,diffc,indc,xc,rcnum] =
       * genIndex(l_svm,k_svm,nfrac); */
      genIndex(varargin_7, varargin_8, nfrac, comb, comb2, diffc, indc, xc,
               &rcnum);
      printf("WARNING: Using %d gapped kmers\n", (int)((double)(comb->size[0] * comb->size[1]) / varargin_8 * pow(4.0, varargin_8)));
      fflush(stdout); 
    } else {
      /* 'gkmPWMlasso:47' else */
      /* 'gkmPWMlasso:48' l_svm2 = l_svm; */
      /* 'gkmPWMlasso:49' k_svm2 = k_svm; */
    }
  } else {
    /* 'gkmPWMlasso:47' else */
    /* 'gkmPWMlasso:48' l_svm2 = l_svm; */
    /* 'gkmPWMlasso:49' k_svm2 = k_svm; */
  }
  emxInit_real_T(&cfile2, 1);
  /* 'gkmPWMlasso:52' fprintf('Counting gapped kmers\n'); */
  printf("Counting gapped kmers\n");
  fflush(stdout);
  /* 'gkmPWMlasso:54' [cfile, GCpos1, GCneg1,mat,mat2] = getgkmcounts(filename,
   * l_svm, k_svm, lk, RC, comb,rcnum); */
  getgkmcounts(varargin_1, varargin_7, varargin_8, lk_data, lk_size,
               varargin_10, comb, rcnum, cfile2, &GCpos1, &nfrac, mat, mat2);
  cfile2_data = cfile2->data;
  /* 'gkmPWMlasso:55' if BG_GC == 1 */
  if (varargin_9 == 1.0) {
    /* 'gkmPWMlasso:56' mat = (mat+mat2)/2; */
    for (i = 0; i < 16; i++) {
      mat[i] = (mat[i] + mat2[i]) / 2.0;
    }
    /* 'gkmPWMlasso:57' GCpos1 = (GCpos1+GCneg1)/2; */
    GCpos1 = (GCpos1 + nfrac) / 2.0;
    /* 'gkmPWMlasso:58' GCneg1 = GCpos1; */
    nfrac = GCpos1;
  }
  emxInit_real_T(&negvec, 1);
  emxInit_char_T(&text, 2);
  /* 'gkmPWMlasso:60' fprintf('Generating negative set\n'); */
  printf("Generating negative set\n");
  fflush(stdout);
  /* 'gkmPWMlasso:61' negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC); */
  BGkmer(mat, nfrac, comb, rcnum, varargin_7, varargin_8, varargin_10, negvec);
  negvec_data = negvec->data;
  /* 'gkmPWMlasso:63' fprintf('Filtering motifs\n'); */
  printf("Filtering motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso:64' num = length(strfind(fileread(memefile),'MOTIF')); */
  fileread(varargin_2, text);
  text_data = text->data;
  emxInit_int32_T(&match_out, 2);
  if (text->size[1] == 0) {
    size_tmp_idx_1 = 0;
  } else {
    emxInit_int32_T(&matches, 2);
    i = matches->size[0] * matches->size[1];
    matches->size[0] = 1;
    matches->size[1] = text->size[1];
    emxEnsureCapacity_int32_T(matches, i);
    matches_data = matches->data;
    match_idx = 0;
    i = text->size[1];
    for (b_i = 0; b_i <= i - 5; b_i++) {
      j = 1;
      while ((j <= 5) && (text_data[(b_i + j) - 1] == cv[j - 1])) {
        j++;
      }
      if (j > 5) {
        matches_data[match_idx] = b_i + 1;
        match_idx++;
      }
    }
    i = match_out->size[0] * match_out->size[1];
    match_out->size[0] = 1;
    match_out->size[1] = match_idx;
    emxEnsureCapacity_int32_T(match_out, i);
    match_out_data = match_out->data;
    for (b_i = 0; b_i < match_idx; b_i++) {
      match_out_data[b_i] = matches_data[b_i];
    }
    emxFree_int32_T(&matches);
    size_tmp_idx_1 = match_out->size[1];
  }
  /* 'gkmPWMlasso:65' p = getmotif(memefile,1:num); */
  emxInit_real_T(&f, 2);
  if (size_tmp_idx_1 < 1) {
    f->size[0] = 1;
    f->size[1] = 0;
  } else {
    i = f->size[0] * f->size[1];
    f->size[0] = 1;
    f->size[1] = size_tmp_idx_1;
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = size_tmp_idx_1 - 1;
    for (i = 0; i <= loop_ub; i++) {
      f_data[i] = (double)i + 1.0;
    }
  }
  emxInit_cell_wrap_0(&p);
  getmotif(varargin_2, f, p);
  p_data = p->data;
  /* 'gkmPWMlasso:66' for i = 1:num */
  for (b_i = 0; b_i < size_tmp_idx_1; b_i++) {
    /* 'gkmPWMlasso:67' [r c] = size(p{i}); */
    /* 'gkmPWMlasso:68' for j = 1:r */
    i = p_data[b_i].f1->size[0];
    for (j = 0; j < i; j++) {
      /* 'gkmPWMlasso:69' a = sum(p{i}(j,:)); */
      loop_ub = p_data[b_i].f1->size[1];
      i1 = f->size[0] * f->size[1];
      f->size[0] = 1;
      f->size[1] = p_data[b_i].f1->size[1];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = p_data[b_i].f1->data[j + p_data[b_i].f1->size[0] * i1];
      }
      a = sum(f);
      /* 'gkmPWMlasso:70' if abs(a-1)>0 */
      if (fabs(a - 1.0) > 0.0) {
        /* 'gkmPWMlasso:71' [b1 loc] = max(p{i}(j,:)); */
        loop_ub = p_data[b_i].f1->size[1];
        i1 = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = p_data[b_i].f1->size[1];
        emxEnsureCapacity_real_T(f, i1);
        f_data = f->data;
        for (i1 = 0; i1 < loop_ub; i1++) {
          f_data[i1] = p_data[b_i].f1->data[j + p_data[b_i].f1->size[0] * i1];
        }
        b_maximum(f, &nfrac, &match_idx);
        /* 'gkmPWMlasso:72' p{i}(j,loc) = b1-a+1; */
        p_data[b_i].f1->data[j + p_data[b_i].f1->size[0] * (match_idx - 1)] =
            (nfrac - a) + 1.0;
      }
    }
  }
  emxInit_real_T(&b_f, 1);
  emxInit_real_T(&loc, 1);
  /* 'gkmPWMlasso:76' [p, info, lenvec] = trim_pwm(p,0.0); */
  trim_pwm(p, loc, b_f);
  b_f_data = b_f->data;
  loc_data = loc->data;
  p_data = p->data;
  /* 'gkmPWMlasso:77' indvec =
   * intersect(find(info./lenvec>=minInfo),find(lenvec>=minL)); */
  emxInit_int32_T(&idx, 1);
  emxInit_boolean_T(&b_loc, 1);
  if (loc->size[0] == b_f->size[0]) {
    i = b_loc->size[0];
    b_loc->size[0] = loc->size[0];
    emxEnsureCapacity_boolean_T(b_loc, i);
    b_loc_data = b_loc->data;
    loop_ub = loc->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_loc_data[i] = (loc_data[i] / b_f_data[i] >= varargin_5);
    }
    b_eml_find(b_loc, idx);
    matches_data = idx->data;
  } else {
    e_binary_expand_op(idx, loc, b_f, varargin_5);
    matches_data = idx->data;
  }
  i = b_loc->size[0];
  b_loc->size[0] = b_f->size[0];
  emxEnsureCapacity_boolean_T(b_loc, i);
  b_loc_data = b_loc->data;
  loop_ub = b_f->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_loc_data[i] = (b_f_data[i] >= varargin_4);
  }
  emxInit_real_T(&BY, 1);
  emxInit_int32_T(&ib, 1);
  b_eml_find(b_loc, ib);
  match_out_data = ib->data;
  i = BY->size[0];
  BY->size[0] = idx->size[0];
  emxEnsureCapacity_real_T(BY, i);
  BY_data = BY->data;
  loop_ub = idx->size[0];
  for (i = 0; i < loop_ub; i++) {
    BY_data[i] = matches_data[i];
  }
  emxInit_real_T(&A, 1);
  i = A->size[0];
  A->size[0] = ib->size[0];
  emxEnsureCapacity_real_T(A, i);
  f_data = A->data;
  loop_ub = ib->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = match_out_data[i];
  }
  emxInit_real_T(&indvec, 1);
  emxInit_real_T(&b_A, 2);
  do_vectors(BY, A, indvec, idx, ib);
  indvec_data = indvec->data;
  /* 'gkmPWMlasso:78' n = length(indvec); */
  /* 'gkmPWMlasso:79' fprintf('Mapping PWMs to gkm space\n'); */
  printf("Mapping PWMs to gkm space\n");
  fflush(stdout);
  /* 'gkmPWMlasso:80' lcnum = numel(comb)/k_svm; */
  lcnum = (double)(comb->size[0] * comb->size[1]) / varargin_8;
  /* 'gkmPWMlasso:81' A=zeros(lcnum*4^k_svm,n); */
  i = (int)(lcnum * pow(4.0, varargin_8));
  i1 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = i;
  b_A->size[1] = indvec->size[0];
  emxEnsureCapacity_real_T(b_A, i1);
  A_data = b_A->data;
  loop_ub = i * indvec->size[0];
  emxFree_int32_T(&ib);
  for (i1 = 0; i1 < loop_ub; i1++) {
    A_data[i1] = 0.0;
  }
  emxInit_real_T(&GCmat, 2);
  emxInit_real_T(&normvec, 1);
  /* 'gkmPWMlasso:82' GCmat = repmat([0.5-GCpos1/2 GCpos1/2 GCpos1/2
   * 0.5-GCpos1/2],l_svm-1,1); */
  d = 0.5 - GCpos1 / 2.0;
  dv[0] = d;
  dv[1] = GCpos1 / 2.0;
  dv[2] = GCpos1 / 2.0;
  dv[3] = d;
  b_repmat(dv, varargin_7 - 1.0, GCmat);
  GCmat_data = GCmat->data;
  /* 'gkmPWMlasso:83' per = 10; */
  a = 10.0;
  /* 'gkmPWMlasso:84' normvec = zeros(n,1); */
  /* 'gkmPWMlasso:85' for j = 1:n */
  i1 = indvec->size[0];
  i2 = normvec->size[0];
  normvec->size[0] = indvec->size[0];
  emxEnsureCapacity_real_T(normvec, i2);
  normvec_data = normvec->data;
  for (j = 0; j < i1; j++) {
    /* 'gkmPWMlasso:86' if mod(j, floor(n/10))==0 */
    if (b_mod((double)j + 1.0, floor((double)indvec->size[0] / 10.0)) == 0.0) {
      /* 'gkmPWMlasso:87' fprintf('%d...', int32(per)); */
      if (a < 2.147483648E+9) {
        i2 = (int)a;
      } else {
        i2 = MAX_int32_T;
      }
      printf("%d...", i2);
      fflush(stdout);
      /* 'gkmPWMlasso:88' per = per+10; */
      a += 10.0;
    }
    /* 'gkmPWMlasso:90' loc = zeros(l_svm*2-2+lenvec(indvec(j)), 1); */
    i2 = (int)indvec_data[j] - 1;
    d = b_f_data[i2];
    loop_ub = (int)((varargin_7 * 2.0 - 2.0) + d);
    b_i = loc->size[0];
    loc->size[0] = loop_ub;
    emxEnsureCapacity_real_T(loc, b_i);
    loc_data = loc->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      loc_data[b_i] = 0.0;
    }
    /* 'gkmPWMlasso:91' loc(l_svm:lenvec(indvec(j))+l_svm-1) = 1; */
    nfrac = (d + varargin_7) - 1.0;
    if (varargin_7 > nfrac) {
      b_i = -1;
      i3 = 0;
    } else {
      b_i = (int)varargin_7 - 2;
      i3 = (int)nfrac;
    }
    loop_ub = (i3 - b_i) - 1;
    for (i3 = 0; i3 < loop_ub; i3++) {
      loc_data[(b_i + i3) + 1] = 1.0;
    }
    /* 'gkmPWMlasso:92' if RC */
    if (varargin_10 != 0.0) {
      /* 'gkmPWMlasso:93' A(:,j) =
       * PWM2kmers([GCmat;p{indvec(j)};GCmat],mat,comb2,diffc,indc,loc,xc,l_svm,k_svm,rcnum);
       */
      if (GCmat->size[0] != 0) {
        loop_ub = 4;
      } else if ((p_data[(int)indvec_data[j] - 1].f1->size[0] != 0) &&
                 (p_data[(int)indvec_data[j] - 1].f1->size[1] != 0)) {
        loop_ub = p_data[(int)indvec_data[j] - 1].f1->size[1];
      } else {
        loop_ub = 4;
      }
      empty_non_axis_sizes = (loop_ub == 0);
      if (empty_non_axis_sizes || (GCmat->size[0] != 0)) {
        b_loop_ub = GCmat->size[0];
      } else {
        b_loop_ub = 0;
      }
      if (empty_non_axis_sizes ||
          ((p_data[(int)indvec_data[j] - 1].f1->size[0] != 0) &&
           (p_data[(int)indvec_data[j] - 1].f1->size[1] != 0))) {
        match_idx = p_data[(int)indvec_data[j] - 1].f1->size[0];
      } else {
        match_idx = 0;
      }
      if (empty_non_axis_sizes || (GCmat->size[0] != 0)) {
        nx = GCmat->size[0];
      } else {
        nx = 0;
      }
      b_i = comb->size[0] * comb->size[1];
      comb->size[0] = (b_loop_ub + match_idx) + nx;
      comb->size[1] = loop_ub;
      emxEnsureCapacity_real_T(comb, b_i);
      comb_data = comb->data;
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          comb_data[i3 + comb->size[0] * b_i] =
              GCmat_data[i3 + b_loop_ub * b_i];
        }
      }
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < match_idx; i3++) {
          comb_data[(i3 + b_loop_ub) + comb->size[0] * b_i] =
              p_data[i2].f1->data[i3 + match_idx * b_i];
        }
      }
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (b_i = 0; b_i < nx; b_i++) {
          comb_data[((b_i + b_loop_ub) + match_idx) + comb->size[0] * i2] =
              GCmat_data[b_i + nx * i2];
        }
      }
      PWM2kmers(comb, mat, comb2, diffc, indc, loc, xc, varargin_7, varargin_8,
                rcnum, A);
      f_data = A->data;
      loop_ub = A->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        A_data[i2 + b_A->size[0] * j] = f_data[i2];
      }
    } else {
      /* 'gkmPWMlasso:94' else */
      /* 'gkmPWMlasso:95' A(:,j) =
       * PWM2kmers_norc([GCmat;p{indvec(j)};GCmat],mat,comb2,diffc,indc,loc,xc,l_svm,k_svm,rcnum);
       */
      if (GCmat->size[0] != 0) {
        loop_ub = 4;
      } else if ((p_data[(int)indvec_data[j] - 1].f1->size[0] != 0) &&
                 (p_data[(int)indvec_data[j] - 1].f1->size[1] != 0)) {
        loop_ub = p_data[(int)indvec_data[j] - 1].f1->size[1];
      } else {
        loop_ub = 4;
      }
      empty_non_axis_sizes = (loop_ub == 0);
      if (empty_non_axis_sizes || (GCmat->size[0] != 0)) {
        b_loop_ub = GCmat->size[0];
      } else {
        b_loop_ub = 0;
      }
      if (empty_non_axis_sizes ||
          ((p_data[(int)indvec_data[j] - 1].f1->size[0] != 0) &&
           (p_data[(int)indvec_data[j] - 1].f1->size[1] != 0))) {
        match_idx = p_data[(int)indvec_data[j] - 1].f1->size[0];
      } else {
        match_idx = 0;
      }
      if (empty_non_axis_sizes || (GCmat->size[0] != 0)) {
        nx = GCmat->size[0];
      } else {
        nx = 0;
      }
      b_i = comb->size[0] * comb->size[1];
      comb->size[0] = (b_loop_ub + match_idx) + nx;
      comb->size[1] = loop_ub;
      emxEnsureCapacity_real_T(comb, b_i);
      comb_data = comb->data;
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          comb_data[i3 + comb->size[0] * b_i] =
              GCmat_data[i3 + b_loop_ub * b_i];
        }
      }
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < match_idx; i3++) {
          comb_data[(i3 + b_loop_ub) + comb->size[0] * b_i] =
              p_data[i2].f1->data[i3 + match_idx * b_i];
        }
      }
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (b_i = 0; b_i < nx; b_i++) {
          comb_data[((b_i + b_loop_ub) + match_idx) + comb->size[0] * i2] =
              GCmat_data[b_i + nx * i2];
        }
      }
      PWM2kmers_norc(comb, mat, comb2, diffc, indc, loc, xc, varargin_7,
                     varargin_8, A);
      f_data = A->data;
      loop_ub = A->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        A_data[i2 + b_A->size[0] * j] = f_data[i2];
      }
    }
    /* 'gkmPWMlasso:97' A(:,j) = A(:,j) - negvec*(l_svm-1+lenvec(indvec(j))); */
    GCpos1 = (varargin_7 - 1.0) + d;
    if (b_A->size[0] == negvec->size[0]) {
      match_idx = b_A->size[0] - 1;
      i2 = A->size[0];
      A->size[0] = b_A->size[0];
      emxEnsureCapacity_real_T(A, i2);
      f_data = A->data;
      for (i2 = 0; i2 <= match_idx; i2++) {
        f_data[i2] = A_data[i2 + b_A->size[0] * j] - negvec_data[i2] * GCpos1;
      }
      loop_ub = A->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        A_data[i2 + b_A->size[0] * j] = f_data[i2];
      }
    } else {
      binary_expand_op(b_A, j, negvec, GCpos1);
      A_data = b_A->data;
    }
    /* 'gkmPWMlasso:98' normvec(j) = sqrt(A(:,j)'*A(:,j)); */
    loop_ub = b_A->size[0];
    i2 = loc->size[0];
    loc->size[0] = b_A->size[0];
    emxEnsureCapacity_real_T(loc, i2);
    loc_data = loc->data;
    i2 = BY->size[0];
    BY->size[0] = b_A->size[0];
    emxEnsureCapacity_real_T(BY, i2);
    BY_data = BY->data;
    for (i2 = 0; i2 < loop_ub; i2++) {
      d = A_data[i2 + b_A->size[0] * j];
      loc_data[i2] = d;
      BY_data[i2] = d;
    }
    if (b_A->size[0] < 1) {
      GCpos1 = 0.0;
    } else {
      GCpos1 = cblas_ddot((blasint)b_A->size[0], &loc_data[0], (blasint)1,
                          &BY_data[0], (blasint)1);
    }
    normvec_data[j] = sqrt(GCpos1);
  }
  emxFree_real_T(&GCmat);
  emxFree_cell_wrap_0(&p);
  /* 'gkmPWMlasso:100' fprintf('\n'); */
  printf("\n");
  fflush(stdout);
  /* 'gkmPWMlasso:101' fprintf('Clustering motifs\n'); */
  printf("Clustering motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso:102' simmat = diag(normvec.^-1)*(A'*A)*diag(normvec.^-1); */
  b_mtimes(b_A, b_A, comb2);
  i1 = loc->size[0];
  loc->size[0] = normvec->size[0];
  emxEnsureCapacity_real_T(loc, i1);
  loc_data = loc->data;
  loop_ub = normvec->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    nfrac = normvec_data[i1];
    loc_data[i1] = pow(nfrac, -1.0);
  }
  emxInit_cell_wrap_1(&motclus);
  emxInit_real_T(&y, 2);
  emxInit_real_T(&b_y, 2);
  diag(loc, comb);
  c_mtimes(comb, comb2, y);
  c_mtimes(y, comb, b_y);
  /* 'gkmPWMlasso:103' motclus = clus_simmat_eig(simmat,corrCut); */
  clus_simmat_eig(b_y, varargin_6, motclus);
  motclus_data = motclus->data;
  /* 'gkmPWMlasso:104' fprintf('Number of motif clusters: %d\n',
   * int32(length(motclus))); */
  printf("Number of motif clusters: %d\n", motclus->size[0]);
  fflush(stdout);
  /* 'gkmPWMlasso:106' fprintf('Selecting Motifs\n'); */
  printf("Selecting Motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso:107' cfile2 = cfile-negvec/sum(negvec)*sum(cfile); */
  GCpos1 = blockedSummation(cfile2, cfile2->size[0]);
  nfrac = blockedSummation(negvec, negvec->size[0]);
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
  if (cfile2->size[0] == negvec->size[0]) {
    loop_ub = cfile2->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      cfile2_data[i1] -= negvec_data[i1] / nfrac * GCpos1;
    }
  } else {
    d_binary_expand_op(cfile2, negvec, nfrac, GCpos1);
    cfile2_data = cfile2->data;
  }
  /* 'gkmPWMlasso:108' cfile2 = cfile2/std(cfile2); */
  d = b_std(cfile2);
  loop_ub = cfile2->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    cfile2_data[i1] /= d;
  }
  emxInit_real_T(&B, 2);
  /* 'gkmPWMlasso:109' B = zeros((4^k_svm)*lcnum, length(motclus)); */
  i1 = B->size[0] * B->size[1];
  B->size[0] = i;
  B->size[1] = motclus->size[0];
  emxEnsureCapacity_real_T(B, i1);
  B_data = B->data;
  loop_ub = i * motclus->size[0];
  for (i = 0; i < loop_ub; i++) {
    B_data[i] = 0.0;
  }
  /* 'gkmPWMlasso:110' corrvec = zeros(n,1); */
  /* 'gkmPWMlasso:111' zvec = zeros(n,1); */
  /* 'gkmPWMlasso:112' Z = zeros(length(motclus),1); */
  i = negvec->size[0];
  negvec->size[0] = motclus->size[0];
  emxEnsureCapacity_real_T(negvec, i);
  negvec_data = negvec->data;
  loop_ub = motclus->size[0];
  for (i = 0; i < loop_ub; i++) {
    negvec_data[i] = 0.0;
  }
  /* 'gkmPWMlasso:114' for i = 1:n */
  i = indvec->size[0];
  i1 = normvec->size[0];
  normvec->size[0] = indvec->size[0];
  emxEnsureCapacity_real_T(normvec, i1);
  normvec_data = normvec->data;
  i1 = diffc->size[0];
  diffc->size[0] = indvec->size[0];
  emxEnsureCapacity_real_T(diffc, i1);
  GCmat_data = diffc->data;
  if (0 <= indvec->size[0] - 1) {
    nfrac = lcnum * 10.0;
    if (nfrac <= b_A->size[0]) {
      k = (int)nfrac;
    } else {
      k = b_A->size[0];
    }
    n = b_A->size[0];
    c_loop_ub = b_A->size[0];
    if (1.0 > lcnum) {
      d_loop_ub = 0;
    } else {
      d_loop_ub = (int)lcnum;
    }
  }
  for (b_i = 0; b_i < i; b_i++) {
    /*  [~,I] = sort(A(:,i),'descend'); */
    /*  zvec(i) = mean(cfile2(I(1:lcnum))); */
    /* 'gkmPWMlasso:117' [~,I2] = maxk(A(:,i), lcnum*10); */
    i1 = A->size[0];
    A->size[0] = n;
    emxEnsureCapacity_real_T(A, i1);
    f_data = A->data;
    for (i1 = 0; i1 < c_loop_ub; i1++) {
      f_data[i1] = A_data[i1 + b_A->size[0] * b_i];
    }
    exkib(A, k, idx, BY);
    matches_data = idx->data;
    loop_ub = idx->size[0];
    i1 = loc->size[0];
    loc->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(loc, i1);
    loc_data = loc->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      loc_data[i1] = matches_data[i1];
    }
    /* 'gkmPWMlasso:118' zvec(i) = mean(cfile2(I2(1:lcnum))); */
    i1 = BY->size[0];
    BY->size[0] = d_loop_ub;
    emxEnsureCapacity_real_T(BY, i1);
    BY_data = BY->data;
    for (i1 = 0; i1 < d_loop_ub; i1++) {
      BY_data[i1] = cfile2_data[(int)loc_data[i1] - 1];
    }
    d = blockedSummation(BY, d_loop_ub) / (double)d_loop_ub;
    normvec_data[b_i] = d;
    /* 'gkmPWMlasso:119' if zvec(i) < 0 */
    if (d < 0.0) {
      /*  Alternative to corr */
      /*  corrvec(i) = -1*corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10))); */
      /* 'gkmPWMlasso:122' correlation_matrix = -1*corrcoef(A(I2,i),
       * cfile2(I2)); */
      /* 'gkmPWMlasso:123' corrvec(i) = correlation_matrix(1,2); */
      i1 = A->size[0];
      A->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(A, i1);
      f_data = A->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = A_data[((int)loc_data[i1] + b_A->size[0] * b_i) - 1];
      }
      i1 = BY->size[0];
      BY->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(BY, i1);
      BY_data = BY->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        BY_data[i1] = cfile2_data[(int)loc_data[i1] - 1];
      }
      corrcoef(A, BY, dv);
      GCmat_data[b_i] = -dv[2];
    } else {
      /* 'gkmPWMlasso:124' else */
      /*  Alternative to corr */
      /*  corrvec(i) = corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10))); */
      /* 'gkmPWMlasso:127' correlation_matrix = corrcoef(A(I2,i), cfile2(I2));
       */
      /* 'gkmPWMlasso:128' corrvec(i) = correlation_matrix(1,2); */
      i1 = A->size[0];
      A->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(A, i1);
      f_data = A->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = A_data[((int)loc_data[i1] + b_A->size[0] * b_i) - 1];
      }
      i1 = BY->size[0];
      BY->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(BY, i1);
      BY_data = BY->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        BY_data[i1] = cfile2_data[(int)loc_data[i1] - 1];
      }
      corrcoef(A, BY, dv);
      GCmat_data[b_i] = dv[2];
    }
  }
  emxFree_real_T(&A);
  /* 'gkmPWMlasso:131' for i = 1:length(motclus) */
  i = motclus->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWMlasso:132' [a,b] = sort(zvec(motclus{i}),'descend'); */
    i1 = loc->size[0];
    loc->size[0] = motclus_data[b_i].f1->size[1];
    emxEnsureCapacity_real_T(loc, i1);
    loc_data = loc->data;
    loop_ub = motclus_data[b_i].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      loc_data[i1] = normvec_data[(int)motclus_data[b_i].f1->data[i1] - 1];
    }
    sort(loc, idx);
    matches_data = idx->data;
    loc_data = loc->data;
    i1 = BY->size[0];
    BY->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(BY, i1);
    BY_data = BY->data;
    loop_ub = idx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      BY_data[i1] = matches_data[i1];
    }
    /* 'gkmPWMlasso:133' f = find(a == a(1)); */
    i1 = b_loc->size[0];
    b_loc->size[0] = loc->size[0];
    emxEnsureCapacity_boolean_T(b_loc, i1);
    b_loc_data = b_loc->data;
    loop_ub = loc->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_loc_data[i1] = (loc_data[i1] == loc_data[0]);
    }
    b_eml_find(b_loc, idx);
    matches_data = idx->data;
    i1 = b_f->size[0];
    b_f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(b_f, i1);
    b_f_data = b_f->data;
    loop_ub = idx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_f_data[i1] = matches_data[i1];
    }
    /* 'gkmPWMlasso:134' if length(f) > 1 */
    if (b_f->size[0] > 1) {
      /* 'gkmPWMlasso:135' [~,bb] = sort(abs(corrvec(motclus{i}(b(f)))),
       * 'descend'); */
      i1 = f->size[0] * f->size[1];
      f->size[0] = 1;
      f->size[1] = b_f->size[0];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      loop_ub = b_f->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] =
            motclus_data[b_i].f1->data[(int)BY_data[(int)b_f_data[i1] - 1] - 1];
      }
      nx = b_f->size[0];
      i1 = indc->size[0];
      indc->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(indc, i1);
      indc_data = indc->data;
      for (k = 0; k < nx; k++) {
        indc_data[k] = fabs(GCmat_data[(int)f_data[k] - 1]);
      }
      sort(indc, idx);
      matches_data = idx->data;
      /* 'gkmPWMlasso:136' b(1:length(f)) = b(bb); */
      i1 = f->size[0] * f->size[1];
      f->size[0] = 1;
      f->size[1] = b_f->size[0];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      loop_ub = b_f->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = BY_data[matches_data[i1] - 1];
      }
      loop_ub = f->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        BY_data[i1] = f_data[i1];
      }
    }
    /* 'gkmPWMlasso:138' B(:,i) = A(:,motclus{i}(b(1))); */
    nx = (int)motclus_data[b_i].f1->data[(int)BY_data[0] - 1];
    loop_ub = b_A->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      B_data[i1 + B->size[0] * b_i] = A_data[i1 + b_A->size[0] * (nx - 1)];
    }
    /* 'gkmPWMlasso:139' motclus{i} = motclus{i}(b); */
    i1 = f->size[0] * f->size[1];
    f->size[0] = 1;
    f->size[1] = BY->size[0];
    emxEnsureCapacity_real_T(f, i1);
    f_data = f->data;
    loop_ub = BY->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[i1] = motclus_data[b_i].f1->data[(int)BY_data[i1] - 1];
    }
    i1 = motclus_data[b_i].f1->size[0] * motclus_data[b_i].f1->size[1];
    motclus_data[b_i].f1->size[0] = 1;
    motclus_data[b_i].f1->size[1] = f->size[1];
    emxEnsureCapacity_real_T(motclus_data[b_i].f1, i1);
    loop_ub = f->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      motclus_data[b_i].f1->data[i1] = f_data[i1];
    }
    /* 'gkmPWMlasso:140' Z(i) = zvec(motclus{i}(1)); */
    negvec_data[b_i] = normvec_data[(int)motclus_data[b_i].f1->data[0] - 1];
  }
  emxFree_real_T(&b_A);
  /* 'gkmPWMlasso:142' f = find(abs(Z)>1); */
  nx = negvec->size[0];
  i = indc->size[0];
  indc->size[0] = negvec->size[0];
  emxEnsureCapacity_real_T(indc, i);
  indc_data = indc->data;
  for (k = 0; k < nx; k++) {
    indc_data[k] = fabs(negvec_data[k]);
  }
  i = b_loc->size[0];
  b_loc->size[0] = indc->size[0];
  emxEnsureCapacity_boolean_T(b_loc, i);
  b_loc_data = b_loc->data;
  loop_ub = indc->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_loc_data[i] = (indc_data[i] > 1.0);
  }
  b_eml_find(b_loc, idx);
  matches_data = idx->data;
  i = b_f->size[0];
  b_f->size[0] = idx->size[0];
  emxEnsureCapacity_real_T(b_f, i);
  b_f_data = b_f->data;
  loop_ub = idx->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_f_data[i] = matches_data[i];
  }
  /* 'gkmPWMlasso:143' B = B(:,f); */
  nx = B->size[0] - 1;
  i = xc->size[0] * xc->size[1];
  xc->size[0] = B->size[0];
  xc->size[1] = b_f->size[0];
  emxEnsureCapacity_real_T(xc, i);
  A_data = xc->data;
  loop_ub = b_f->size[0];
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 <= nx; i1++) {
      A_data[i1 + xc->size[0] * i] =
          B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
    }
  }
  i = B->size[0] * B->size[1];
  B->size[0] = xc->size[0];
  B->size[1] = xc->size[1];
  emxEnsureCapacity_real_T(B, i);
  B_data = B->data;
  loop_ub = xc->size[0] * xc->size[1];
  for (i = 0; i < loop_ub; i++) {
    B_data[i] = A_data[i];
  }
  /* 'gkmPWMlasso:144' B = B/std(B(:))'; */
  nx = B->size[0] * B->size[1];
  b_B = *B;
  b_motclus = nx;
  b_B.size = &b_motclus;
  b_B.numDimensions = 1;
  d = b_std(&b_B);
  loop_ub = B->size[0] * B->size[1];
  for (i = 0; i < loop_ub; i++) {
    B_data[i] /= d;
  }
  emxInit_cell_wrap_1(&tmp_motclus);
  /*  motclus = motclus(f); */
  /* 'gkmPWMlasso:146' tmp_motclus = cell(length(f),1); */
  match_idx = b_f->size[0];
  i = tmp_motclus->size[0];
  tmp_motclus->size[0] = b_f->size[0];
  emxEnsureCapacity_cell_wrap_1(tmp_motclus, i);
  tmp_motclus_data = tmp_motclus->data;
  for (i = 0; i < match_idx; i++) {
    tmp_motclus_data[i].f1->size[0] = 1;
    tmp_motclus_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso:147' tmp_motclus = coder.nullcopy(tmp_motclus); */
  /* 'gkmPWMlasso:148' for idx = 1:length(f) */
  i = b_f->size[0];
  for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
    /* 'gkmPWMlasso:149' tmp_motclus{idx} = motclus{f(idx)}; */
    i1 = tmp_motclus_data[b_loop_ub].f1->size[0] *
         tmp_motclus_data[b_loop_ub].f1->size[1];
    tmp_motclus_data[b_loop_ub].f1->size[0] = 1;
    tmp_motclus_data[b_loop_ub].f1->size[1] =
        motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->size[1];
    emxEnsureCapacity_real_T(tmp_motclus_data[b_loop_ub].f1, i1);
    loop_ub = motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      tmp_motclus_data[b_loop_ub].f1->data[i1] =
          motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->data[i1];
    }
  }
  emxInit_cell_wrap_1(&c_motclus);
  /* 'gkmPWMlasso:151' motclus = cell(length(f),1); */
  match_idx = b_f->size[0];
  i = c_motclus->size[0];
  c_motclus->size[0] = b_f->size[0];
  emxEnsureCapacity_cell_wrap_1(c_motclus, i);
  b_motclus_data = c_motclus->data;
  for (i = 0; i < match_idx; i++) {
    b_motclus_data[i].f1->size[0] = 1;
    b_motclus_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso:152' motclus = coder.nullcopy(motclus); */
  i = motclus->size[0];
  motclus->size[0] = c_motclus->size[0];
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  /* 'gkmPWMlasso:153' for idx = 1:length(motclus) */
  i = c_motclus->size[0];
  for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
    /* 'gkmPWMlasso:154' motclus{idx} = tmp_motclus{idx}; */
    i1 = motclus_data[b_loop_ub].f1->size[0] *
         motclus_data[b_loop_ub].f1->size[1];
    motclus_data[b_loop_ub].f1->size[0] = 1;
    motclus_data[b_loop_ub].f1->size[1] =
        tmp_motclus_data[b_loop_ub].f1->size[1];
    emxEnsureCapacity_real_T(motclus_data[b_loop_ub].f1, i1);
    loop_ub = tmp_motclus_data[b_loop_ub].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      motclus_data[b_loop_ub].f1->data[i1] =
          tmp_motclus_data[b_loop_ub].f1->data[i1];
    }
  }
  /* 'gkmPWMlasso:156' Z = Z(f); */
  i = BY->size[0];
  BY->size[0] = b_f->size[0];
  emxEnsureCapacity_real_T(BY, i);
  BY_data = BY->data;
  loop_ub = b_f->size[0];
  for (i = 0; i < loop_ub; i++) {
    BY_data[i] = negvec_data[(int)b_f_data[i] - 1];
  }
  i = negvec->size[0];
  negvec->size[0] = BY->size[0];
  emxEnsureCapacity_real_T(negvec, i);
  negvec_data = negvec->data;
  loop_ub = BY->size[0];
  for (i = 0; i < loop_ub; i++) {
    negvec_data[i] = BY_data[i];
  }
  /*  clear A loc mat GMmat */
  /* 'gkmPWMlasso:159' if d == 0 */
  emxInit_real_T(&OLS, 1);
  emxInit_real_T(&BB, 2);
  if (varargin_3 == 0.0) {
    /* 'gkmPWMlasso:160' fprintf('Running LASSO\n'); */
    printf("Running LASSO\n");
    fflush(stdout);
    /*  [weigmat, FitInfo] = lasso_cvmat(B, cfile2,'DFmax',
     * length(Z),'Standardize', false, 'NumLambda', 20); */
    /* 'gkmPWMlasso:162' [weigmat, FitInfo] = lasso_cvmat(B, cfile2, length(Z),
     * false, 20); */
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = B->size[1];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    loop_ub = B->size[0] * B->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      comb_data[i] = B_data[i];
    }
    i = diffc->size[0];
    diffc->size[0] = cfile2->size[0];
    emxEnsureCapacity_real_T(diffc, i);
    GCmat_data = diffc->data;
    loop_ub = cfile2->size[0] - 1;
    for (i = 0; i <= loop_ub; i++) {
      GCmat_data[i] = cfile2_data[i];
    }
    emxInit_real_T(&weigmat, 2);
    lasso_cvmat(comb, diffc, negvec->size[0], weigmat, expl_temp_data, lk_size,
                b_expl_temp_data, expl_temp_size, &nfrac, a__3_DF_data,
                a__3_DF_size, c_expl_temp_data, b_expl_temp_size);
    GCmat_data = weigmat->data;
    /* 'gkmPWMlasso:163' MSE = zeros(length(FitInfo.DF),1); */
    if (a__3_DF_size[1] == 0) {
      n = 0;
    } else {
      n = a__3_DF_size[1];
    }
    if (0 <= n - 1) {
      memset(&MSE_data[0], 0, n * sizeof(double));
    }
    /* 'gkmPWMlasso:164' cnorm = cfile2'*cfile2; */
    if (cfile2->size[0] < 1) {
      nfrac = 0.0;
    } else {
      nfrac = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                         &cfile2_data[0], (blasint)1);
    }
    emxInit_cell_wrap_3_20(&F);
    /* 'gkmPWMlasso:165' F = cell(numel(FitInfo.DF),1); */
    match_idx = a__3_DF_size[1];
    i = F.size[0];
    F.size[0] = a__3_DF_size[1];
    emxEnsureCapacity_cell_wrap_3(F.data, a__3_DF_size[1], i);
    for (i = 0; i < match_idx; i++) {
      F.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso:166' F = coder.nullcopy(F); */
    /* 'gkmPWMlasso:168' f = find(weigmat(:,1)~=0); */
    loop_ub = weigmat->size[0];
    i = b_loc->size[0];
    b_loc->size[0] = weigmat->size[0];
    emxEnsureCapacity_boolean_T(b_loc, i);
    b_loc_data = b_loc->data;
    for (i = 0; i < loop_ub; i++) {
      b_loc_data[i] = (GCmat_data[i] != 0.0);
    }
    b_eml_find(b_loc, idx);
    matches_data = idx->data;
    i = b_f->size[0];
    b_f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(b_f, i);
    b_f_data = b_f->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_f_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso:169' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    loop_ub = B->size[0];
    i = xc->size[0] * xc->size[1];
    xc->size[0] = B->size[0];
    xc->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(xc, i);
    A_data = xc->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + xc->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    b_mtimes(xc, comb, comb2);
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (b_f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)b_f->size[0], (blasint)1, (blasint)B->size[0], 1.0,
                  &comb_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)b_f->size[0]);
    }
    mpower(comb2, xc);
    A_data = xc->data;
    if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
      i = loc->size[0];
      loc->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(loc, i);
    } else {
      i = loc->size[0];
      loc->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                  &A_data[0], (blasint)xc->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &loc_data[0],
                  (blasint)xc->size[0]);
    }
    i = OLS->size[0];
    OLS->size[0] = loc->size[0];
    emxEnsureCapacity_real_T(OLS, i);
    /* 'gkmPWMlasso:171' for i = 1:length(FitInfo.DF) */
    if (a__3_DF_size[1] == 0) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = a__3_DF_size[1];
    }
    if (0 <= b_loop_ub - 1) {
      i4 = weigmat->size[0];
      e_loop_ub = weigmat->size[0];
      f_loop_ub = B->size[0];
      g_loop_ub = B->size[0];
      h_loop_ub = B->size[0];
      i5 = B->size[0];
      i_loop_ub = B->size[0];
      i6 = B->size[0];
    }
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      /* 'gkmPWMlasso:172' f = find(weigmat(:,i)~=0); */
      i = b_loc->size[0];
      b_loc->size[0] = i4;
      emxEnsureCapacity_boolean_T(b_loc, i);
      b_loc_data = b_loc->data;
      for (i = 0; i < e_loop_ub; i++) {
        b_loc_data[i] = (GCmat_data[i + weigmat->size[0] * b_i] != 0.0);
      }
      b_eml_find(b_loc, idx);
      matches_data = idx->data;
      i = b_f->size[0];
      b_f->size[0] = idx->size[0];
      emxEnsureCapacity_real_T(b_f, i);
      b_f_data = b_f->data;
      loop_ub = idx->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_f_data[i] = matches_data[i];
      }
      /* 'gkmPWMlasso:173' F{i} = f; */
      i = F.data[b_i].f1->size[0];
      F.data[b_i].f1->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(F.data[b_i].f1, i);
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        F.data[b_i].f1->data[i] = b_f_data[i];
      }
      /* 'gkmPWMlasso:174' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
      i = xc->size[0] * xc->size[1];
      xc->size[0] = f_loop_ub;
      xc->size[1] = b_f->size[0];
      emxEnsureCapacity_real_T(xc, i);
      A_data = xc->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < f_loop_ub; i1++) {
          A_data[i1 + xc->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
        }
      }
      i = comb->size[0] * comb->size[1];
      comb->size[0] = g_loop_ub;
      comb->size[1] = b_f->size[0];
      emxEnsureCapacity_real_T(comb, i);
      comb_data = comb->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < g_loop_ub; i1++) {
          comb_data[i1 + comb->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
        }
      }
      b_mtimes(xc, comb, comb2);
      i = comb->size[0] * comb->size[1];
      comb->size[0] = h_loop_ub;
      comb->size[1] = b_f->size[0];
      emxEnsureCapacity_real_T(comb, i);
      comb_data = comb->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < h_loop_ub; i1++) {
          comb_data[i1 + comb->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
        }
      }
      if ((i5 == 0) || (b_f->size[0] == 0) || (cfile2->size[0] == 0)) {
        i = indc->size[0];
        indc->size[0] = b_f->size[0];
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        loop_ub = b_f->size[0];
        for (i = 0; i < loop_ub; i++) {
          indc_data[i] = 0.0;
        }
      } else {
        i = indc->size[0];
        indc->size[0] = b_f->size[0];
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    (blasint)b_f->size[0], (blasint)1, (blasint)B->size[0], 1.0,
                    &comb_data[0], (blasint)B->size[0], &cfile2_data[0],
                    (blasint)cfile2->size[0], 0.0, &indc_data[0],
                    (blasint)b_f->size[0]);
      }
      mpower(comb2, xc);
      A_data = xc->data;
      if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
        i = OLS->size[0];
        OLS->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(OLS, i);
        OLS_data = OLS->data;
        loop_ub = xc->size[0];
        for (i = 0; i < loop_ub; i++) {
          OLS_data[i] = 0.0;
        }
      } else {
        i = OLS->size[0];
        OLS->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(OLS, i);
        OLS_data = OLS->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                    &A_data[0], (blasint)xc->size[0], &indc_data[0],
                    (blasint)indc->size[0], 0.0, &OLS_data[0],
                    (blasint)xc->size[0]);
      }
      /* 'gkmPWMlasso:175' res = B(:,f)*OLS; */
      i = comb->size[0] * comb->size[1];
      comb->size[0] = i_loop_ub;
      comb->size[1] = b_f->size[0];
      emxEnsureCapacity_real_T(comb, i);
      comb_data = comb->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < i_loop_ub; i1++) {
          comb_data[i1 + comb->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
        }
      }
      if ((i6 == 0) || (b_f->size[0] == 0) || (OLS->size[0] == 0)) {
        loop_ub = B->size[0];
        i = loc->size[0];
        loc->size[0] = B->size[0];
        emxEnsureCapacity_real_T(loc, i);
        loc_data = loc->data;
        for (i = 0; i < loop_ub; i++) {
          loc_data[i] = 0.0;
        }
      } else {
        i = loc->size[0];
        loc->size[0] = B->size[0];
        emxEnsureCapacity_real_T(loc, i);
        loc_data = loc->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)B->size[0], (blasint)1, (blasint)b_f->size[0], 1.0,
                    &comb_data[0], (blasint)B->size[0], &OLS_data[0],
                    (blasint)OLS->size[0], 0.0, &loc_data[0],
                    (blasint)B->size[0]);
      }
      /* 'gkmPWMlasso:176' if sum(abs(res)>0) */
      nx = loc->size[0];
      i = indc->size[0];
      indc->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      for (k = 0; k < nx; k++) {
        indc_data[k] = fabs(loc_data[k]);
      }
      i = b_loc->size[0];
      b_loc->size[0] = indc->size[0];
      emxEnsureCapacity_boolean_T(b_loc, i);
      b_loc_data = b_loc->data;
      loop_ub = indc->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_loc_data[i] = (indc_data[i] > 0.0);
      }
      nx = b_loc->size[0];
      if (b_loc->size[0] == 0) {
        match_idx = 0;
      } else {
        match_idx = b_loc_data[0];
        for (k = 2; k <= nx; k++) {
          match_idx += b_loc_data[k - 1];
        }
      }
      if (match_idx != 0) {
        /* 'gkmPWMlasso:177' MSE(i) = (cfile2'*res)^2/(res'*res)/cnorm; */
        if (cfile2->size[0] < 1) {
          a = 0.0;
        } else {
          a = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                         &loc_data[0], (blasint)1);
        }
        if (loc->size[0] < 1) {
          GCpos1 = 0.0;
        } else {
          GCpos1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                              &loc_data[0], (blasint)1);
        }
        MSE_data[b_i] = a * a / GCpos1 / nfrac;
      }
    }
    emxFree_real_T(&weigmat);
    emxInit_cell_wrap_3_1x19(&FF);
    /* 'gkmPWMlasso:180' FF = cell(1,length(MSE)-1); */
    match_idx = n - 1;
    i = FF.size[0] * FF.size[1];
    FF.size[0] = 1;
    FF.size[1] = n - 1;
    emxEnsureCapacity_cell_wrap_31(FF.data, FF.size, i);
    for (i = 0; i < match_idx; i++) {
      FF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso:181' FF = coder.nullcopy(FF); */
    /* 'gkmPWMlasso:182' MSE2 = zeros(1,length(MSE)-1); */
    MSE2_size[0] = 1;
    if (0 <= match_idx - 1) {
      memset(&MSE2_data[0], 0, match_idx * sizeof(double));
    }
    /* 'gkmPWMlasso:183' count = 1; */
    nx = 0;
    /* 'gkmPWMlasso:184' for i = 1:length(MSE)-1 */
    for (b_i = 0; b_i <= n - 2; b_i++) {
      /* 'gkmPWMlasso:185' if numel(setdiff(F{i}, F{i+1})) > 0 */
      b_do_vectors(F.data[b_i].f1, F.data[b_i + 1].f1, loc, idx, &match_idx);
      if (loc->size[0] > 0) {
        /* 'gkmPWMlasso:186' FF{count} = setdiff(F{i}, F{i+1}); */
        b_do_vectors(F.data[b_i].f1, F.data[b_i + 1].f1, FF.data[nx].f1, idx,
                     &match_idx);
        /* 'gkmPWMlasso:187' MSE2(count) = (MSE(i)-MSE(i+1)); */
        MSE2_data[nx] = MSE_data[b_i] - MSE_data[b_i + 1];
        /* 'gkmPWMlasso:188' count = count +1; */
        nx++;
      }
    }
    emxFree_cell_wrap_3_20(&F);
    emxInit_cell_wrap_3_1x19(&newFF);
    /* 'gkmPWMlasso:191' MSE2 = MSE2(1:count-1); */
    if (1 > nx) {
      MSE2_size[1] = 0;
    } else {
      MSE2_size[1] = nx;
    }
    /* 'gkmPWMlasso:192' newFF = cell(1,count-1); */
    i = newFF.size[0] * newFF.size[1];
    newFF.size[0] = 1;
    newFF.size[1] = nx;
    emxEnsureCapacity_cell_wrap_31(newFF.data, newFF.size, i);
    for (i = 0; i < nx; i++) {
      newFF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso:193' newFF = coder.nullcopy(newFF); */
    /* 'gkmPWMlasso:194' for idx = 1:count-1 */
    for (b_loop_ub = 0; b_loop_ub < nx; b_loop_ub++) {
      /* 'gkmPWMlasso:195' newFF{idx} = FF{idx}; */
      loop_ub = FF.data[b_loop_ub].f1->size[0];
      i = newFF.data[b_loop_ub].f1->size[0];
      newFF.data[b_loop_ub].f1->size[0] = FF.data[b_loop_ub].f1->size[0];
      emxEnsureCapacity_real_T(newFF.data[b_loop_ub].f1, i);
      for (i = 0; i < loop_ub; i++) {
        newFF.data[b_loop_ub].f1->data[i] = FF.data[b_loop_ub].f1->data[i];
      }
    }
    emxInit_cell_wrap_3_1x19(&b_FF);
    /* 'gkmPWMlasso:197' FF = cell(1,count-1); */
    i = b_FF.size[0] * b_FF.size[1];
    b_FF.size[0] = 1;
    b_FF.size[1] = nx;
    emxEnsureCapacity_cell_wrap_31(b_FF.data, b_FF.size, i);
    for (i = 0; i < nx; i++) {
      b_FF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso:198' FF = coder.nullcopy(FF); */
    i = FF.size[0] * FF.size[1];
    FF.size[0] = 1;
    FF.size[1] = b_FF.size[1];
    emxEnsureCapacity_cell_wrap_31(FF.data, FF.size, i);
    /* 'gkmPWMlasso:199' for idx = 1:count-1 */
    emxFree_cell_wrap_3_1x19(&b_FF);
    for (b_loop_ub = 0; b_loop_ub < nx; b_loop_ub++) {
      /* 'gkmPWMlasso:200' FF{idx} = newFF{idx}; */
      loop_ub = newFF.data[b_loop_ub].f1->size[0];
      i = FF.data[b_loop_ub].f1->size[0];
      FF.data[b_loop_ub].f1->size[0] = newFF.data[b_loop_ub].f1->size[0];
      emxEnsureCapacity_real_T(FF.data[b_loop_ub].f1, i);
      for (i = 0; i < loop_ub; i++) {
        FF.data[b_loop_ub].f1->data[i] = newFF.data[b_loop_ub].f1->data[i];
      }
    }
    emxFree_cell_wrap_3_1x19(&newFF);
    /* 'gkmPWMlasso:203' res = B*((B.'*B)^-1*(B.'*cfile2)); */
    b_mtimes(B, B, comb2);
    if ((B->size[0] == 0) || (B->size[1] == 0) || (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = B->size[1];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = B->size[1];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = B->size[1];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)B->size[1],
                  (blasint)1, (blasint)B->size[0], 1.0, &B_data[0],
                  (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)B->size[1]);
    }
    mpower(comb2, xc);
    A_data = xc->data;
    if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = xc->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                  &A_data[0], (blasint)xc->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &BY_data[0],
                  (blasint)xc->size[0]);
    }
    if ((B->size[0] == 0) || (B->size[1] == 0) || (BY->size[0] == 0)) {
      i = loc->size[0];
      loc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = B->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = 0.0;
      }
    } else {
      i = loc->size[0];
      loc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)B->size[0], (blasint)1, (blasint)B->size[1], 1.0,
                  &B_data[0], (blasint)B->size[0], &BY_data[0],
                  (blasint)BY->size[0], 0.0, &loc_data[0], (blasint)B->size[0]);
    }
    /* 'gkmPWMlasso:204' csm = (cfile2'*res)^2/(res'*res)/cnorm; */
    if (cfile2->size[0] < 1) {
      a = 0.0;
    } else {
      a = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                     &loc_data[0], (blasint)1);
    }
    if (loc->size[0] < 1) {
      GCpos1 = 0.0;
    } else {
      GCpos1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                          &loc_data[0], (blasint)1);
    }
    nfrac = a * a / GCpos1 / nfrac;
    /* 'gkmPWMlasso:205' [a,b] = sort(MSE2,'descend'); */
    e_sort(MSE2_data, MSE2_size, iidx_data, lk_size);
    /* 'gkmPWMlasso:206' cs = cumsum(a); */
    useConstantDim(MSE2_data, MSE2_size);
    /* 'gkmPWMlasso:207' cs = cs/csm; */
    loop_ub = MSE2_size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      MSE2_data[i] /= nfrac;
    }
    /* 'gkmPWMlasso:208' f = find(cs>0.9); */
    b_MSE2_size[0] = 1;
    b_MSE2_size[1] = MSE2_size[1];
    loop_ub = MSE2_size[1];
    for (i = 0; i < loop_ub; i++) {
      c_MSE2_data[i] = (MSE2_data[i] > 0.9);
    }
    b_MSE2_data.data = &c_MSE2_data[0];
    b_MSE2_data.size = &b_MSE2_size[0];
    b_MSE2_data.allocatedSize = 19;
    b_MSE2_data.numDimensions = 2;
    b_MSE2_data.canFreeData = false;
    eml_find(&b_MSE2_data, match_out);
    match_out_data = match_out->data;
    i = f->size[0] * f->size[1];
    f->size[0] = 1;
    f->size[1] = match_out->size[1];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = match_out->size[1];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = match_out_data[i];
    }
    /* 'gkmPWMlasso:210' if isempty(f) */
    emxInit_real_T(&c_f, 2);
    if (f->size[1] == 0) {
      /* 'gkmPWMlasso:211' f = 1:length(OLS); */
      if (OLS->size[0] < 1) {
        f->size[0] = 1;
        f->size[1] = 0;
      } else {
        i = f->size[0] * f->size[1];
        f->size[0] = 1;
        f->size[1] = OLS->size[0];
        emxEnsureCapacity_real_T(f, i);
        f_data = f->data;
        loop_ub = OLS->size[0] - 1;
        for (i = 0; i <= loop_ub; i++) {
          f_data[i] = (double)i + 1.0;
        }
      }
      i = c_f->size[0] * c_f->size[1];
      c_f->size[0] = 1;
      c_f->size[1] = f->size[1];
      emxEnsureCapacity_real_T(c_f, i);
      b_f_data = c_f->data;
      loop_ub = f->size[1];
      for (i = 0; i < loop_ub; i++) {
        b_f_data[i] = f_data[i];
      }
    } else {
      emxInit_cell_wrap_3(&b_newFF, 2);
      /* 'gkmPWMlasso:212' else */
      /*  FF = FF(b(1:f(1))); */
      /* 'gkmPWMlasso:214' endIdx = f(1); */
      /* 'gkmPWMlasso:215' newFF = cell(1, endIdx); */
      match_idx = (int)f_data[0];
      i = b_newFF->size[0] * b_newFF->size[1];
      b_newFF->size[0] = 1;
      b_newFF->size[1] = match_idx;
      emxEnsureCapacity_cell_wrap_32(b_newFF, i);
      newFF_data = b_newFF->data;
      for (i = 0; i < match_idx; i++) {
        newFF_data[i].f1->size[0] = 0;
      }
      /* 'gkmPWMlasso:216' newFF = coder.nullcopy(newFF); */
      /* 'gkmPWMlasso:217' for idx = 1:endIdx */
      for (b_loop_ub = 0; b_loop_ub < match_idx; b_loop_ub++) {
        /* 'gkmPWMlasso:218' newFF{idx} = FF{b(idx)}; */
        i = newFF_data[b_loop_ub].f1->size[0];
        i1 = iidx_data[b_loop_ub];
        newFF_data[b_loop_ub].f1->size[0] = FF.data[i1 - 1].f1->size[0];
        emxEnsureCapacity_real_T(newFF_data[b_loop_ub].f1, i);
        loop_ub = FF.data[i1 - 1].f1->size[0];
        for (i = 0; i < loop_ub; i++) {
          newFF_data[b_loop_ub].f1->data[i] =
              FF.data[iidx_data[b_loop_ub] - 1].f1->data[i];
        }
      }
      /* 'gkmPWMlasso:220' FF = newFF; */
      /* 'gkmPWMlasso:222' f = []; */
      b_f->size[0] = 0;
      /* 'gkmPWMlasso:223' for i = 1:length(FF) */
      i = b_newFF->size[1];
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso:224' f = [f;FF{i}]; */
        i1 = b_f->size[0];
        loop_ub = newFF_data[b_i].f1->size[0];
        i2 = b_f->size[0];
        b_f->size[0] += newFF_data[b_i].f1->size[0];
        emxEnsureCapacity_real_T(b_f, i2);
        b_f_data = b_f->data;
        for (i2 = 0; i2 < loop_ub; i2++) {
          b_f_data[i1 + i2] = newFF_data[b_i].f1->data[i2];
        }
      }
      emxFree_cell_wrap_3(&b_newFF);
      /* 'gkmPWMlasso:226' f=unique(f); */
      unique_vector(b_f, loc);
      loc_data = loc->data;
      i = c_f->size[0] * c_f->size[1];
      c_f->size[0] = loc->size[0];
      c_f->size[1] = 1;
      emxEnsureCapacity_real_T(c_f, i);
      b_f_data = c_f->data;
      loop_ub = loc->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_f_data[i] = loc_data[i];
      }
    }
    emxFree_cell_wrap_3_1x19(&FF);
    /* f = find(weigmat(:,1)~=0); */
    /* 'gkmPWMlasso:231' F = length(f); */
    /* 'gkmPWMlasso:232' fprintf('Selecting top motifs\n'); */
    printf("Selecting top motifs\n");
    fflush(stdout);
    /* 'gkmPWMlasso:233' ind = true; */
    empty_non_axis_sizes = true;
    /* 'gkmPWMlasso:234' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    match_idx = c_f->size[0] * c_f->size[1];
    nx = c_f->size[0] * c_f->size[1];
    loop_ub = B->size[0];
    i = xc->size[0] * xc->size[1];
    xc->size[0] = B->size[0];
    xc->size[1] = match_idx;
    emxEnsureCapacity_real_T(xc, i);
    A_data = xc->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + xc->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = nx;
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    b_mtimes(xc, comb, comb2);
    match_idx = c_f->size[0] * c_f->size[1];
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = match_idx;
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (c_f->size[0] * c_f->size[1] == 0) ||
        (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = c_f->size[0] * c_f->size[1];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = c_f->size[0] * c_f->size[1];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = c_f->size[0] * c_f->size[1];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)(c_f->size[0] * c_f->size[1]), (blasint)1,
                  (blasint)B->size[0], 1.0, &comb_data[0], (blasint)B->size[0],
                  &cfile2_data[0], (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)(c_f->size[0] * c_f->size[1]));
    }
    mpower(comb2, xc);
    A_data = xc->data;
    if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
      i = OLS->size[0];
      OLS->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      loop_ub = xc->size[0];
      for (i = 0; i < loop_ub; i++) {
        OLS_data[i] = 0.0;
      }
    } else {
      i = OLS->size[0];
      OLS->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                  &A_data[0], (blasint)xc->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &OLS_data[0],
                  (blasint)xc->size[0]);
    }
    emxInit_real_T(&Pweig, 2);
    /* 'gkmPWMlasso:235' Pweig = Z(f); */
    i = Pweig->size[0] * Pweig->size[1];
    Pweig->size[0] = c_f->size[0];
    Pweig->size[1] = c_f->size[1];
    emxEnsureCapacity_real_T(Pweig, i);
    Pweig_data = Pweig->data;
    loop_ub = c_f->size[0] * c_f->size[1];
    for (i = 0; i < loop_ub; i++) {
      Pweig_data[i] = negvec_data[(int)b_f_data[i] - 1];
    }
    /* 'gkmPWMlasso:236' while ind */
    while (empty_non_axis_sizes) {
      /* 'gkmPWMlasso:237' ff = []; */
      loc->size[0] = 0;
      /* 'gkmPWMlasso:238' for i = 1:length(f) */
      if ((c_f->size[0] == 0) || (c_f->size[1] == 0)) {
        nx = 0;
      } else {
        match_idx = c_f->size[0];
        nx = c_f->size[1];
        if (match_idx >= nx) {
          nx = match_idx;
        }
      }
      for (b_i = 0; b_i < nx; b_i++) {
        /* 'gkmPWMlasso:239' if sign(OLS(i)) ~= sign(Pweig(i)) */
        a = OLS_data[b_i];
        if (OLS_data[b_i] < 0.0) {
          a = -1.0;
        } else if (OLS_data[b_i] > 0.0) {
          a = 1.0;
        }
        nfrac = Pweig_data[b_i];
        if (Pweig_data[b_i] < 0.0) {
          nfrac = -1.0;
        } else if (Pweig_data[b_i] > 0.0) {
          nfrac = 1.0;
        }
        if (a != nfrac) {
          /* 'gkmPWMlasso:240' ff = [ff;i]; */
          i = loc->size[0];
          i1 = loc->size[0];
          loc->size[0]++;
          emxEnsureCapacity_real_T(loc, i1);
          loc_data = loc->data;
          loc_data[i] = (double)b_i + 1.0;
        }
      }
      /* 'gkmPWMlasso:243' if length(ff) > 0 */
      if (loc->size[0] > 0) {
        /* 'gkmPWMlasso:244' f(ff) = []; */
        i = idx->size[0];
        idx->size[0] = loc->size[0];
        emxEnsureCapacity_int32_T(idx, i);
        matches_data = idx->data;
        loop_ub = loc->size[0];
        for (i = 0; i < loop_ub; i++) {
          matches_data[i] = (int)loc_data[i];
        }
        c_nullAssignment(c_f, idx);
        b_f_data = c_f->data;
        /* 'gkmPWMlasso:245' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
        match_idx = c_f->size[0] * c_f->size[1];
        nx = c_f->size[0] * c_f->size[1];
        loop_ub = B->size[0];
        i = xc->size[0] * xc->size[1];
        xc->size[0] = B->size[0];
        xc->size[1] = match_idx;
        emxEnsureCapacity_real_T(xc, i);
        A_data = xc->data;
        for (i = 0; i < match_idx; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            A_data[i1 + xc->size[0] * i] =
                B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
          }
        }
        loop_ub = B->size[0];
        i = comb->size[0] * comb->size[1];
        comb->size[0] = B->size[0];
        comb->size[1] = nx;
        emxEnsureCapacity_real_T(comb, i);
        comb_data = comb->data;
        for (i = 0; i < nx; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            comb_data[i1 + comb->size[0] * i] =
                B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
          }
        }
        b_mtimes(xc, comb, comb2);
        match_idx = c_f->size[0] * c_f->size[1];
        loop_ub = B->size[0];
        i = comb->size[0] * comb->size[1];
        comb->size[0] = B->size[0];
        comb->size[1] = match_idx;
        emxEnsureCapacity_real_T(comb, i);
        comb_data = comb->data;
        for (i = 0; i < match_idx; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            comb_data[i1 + comb->size[0] * i] =
                B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
          }
        }
        if ((B->size[0] == 0) || (c_f->size[0] * c_f->size[1] == 0) ||
            (cfile2->size[0] == 0)) {
          i = indc->size[0];
          indc->size[0] = c_f->size[0] * c_f->size[1];
          emxEnsureCapacity_real_T(indc, i);
          indc_data = indc->data;
          loop_ub = c_f->size[0] * c_f->size[1];
          for (i = 0; i < loop_ub; i++) {
            indc_data[i] = 0.0;
          }
        } else {
          i = indc->size[0];
          indc->size[0] = c_f->size[0] * c_f->size[1];
          emxEnsureCapacity_real_T(indc, i);
          indc_data = indc->data;
          cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                      (blasint)(c_f->size[0] * c_f->size[1]), (blasint)1,
                      (blasint)B->size[0], 1.0, &comb_data[0],
                      (blasint)B->size[0], &cfile2_data[0],
                      (blasint)cfile2->size[0], 0.0, &indc_data[0],
                      (blasint)(c_f->size[0] * c_f->size[1]));
        }
        mpower(comb2, xc);
        A_data = xc->data;
        if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
          i = OLS->size[0];
          OLS->size[0] = xc->size[0];
          emxEnsureCapacity_real_T(OLS, i);
          OLS_data = OLS->data;
          loop_ub = xc->size[0];
          for (i = 0; i < loop_ub; i++) {
            OLS_data[i] = 0.0;
          }
        } else {
          i = OLS->size[0];
          OLS->size[0] = xc->size[0];
          emxEnsureCapacity_real_T(OLS, i);
          OLS_data = OLS->data;
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                      (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1],
                      1.0, &A_data[0], (blasint)xc->size[0], &indc_data[0],
                      (blasint)indc->size[0], 0.0, &OLS_data[0],
                      (blasint)xc->size[0]);
        }
        /* 'gkmPWMlasso:246' Pweig = Z(f); */
        i = Pweig->size[0] * Pweig->size[1];
        Pweig->size[0] = c_f->size[0];
        Pweig->size[1] = c_f->size[1];
        emxEnsureCapacity_real_T(Pweig, i);
        Pweig_data = Pweig->data;
        loop_ub = c_f->size[0] * c_f->size[1];
        for (i = 0; i < loop_ub; i++) {
          Pweig_data[i] = negvec_data[(int)b_f_data[i] - 1];
        }
      } else {
        /* 'gkmPWMlasso:247' else */
        /* 'gkmPWMlasso:248' ind = false; */
        empty_non_axis_sizes = false;
      }
    }
    /* 'gkmPWMlasso:251' BB = B(:,f); */
    match_idx = c_f->size[0] * c_f->size[1];
    loop_ub = B->size[0];
    i = BB->size[0] * BB->size[1];
    BB->size[0] = B->size[0];
    BB->size[1] = match_idx;
    emxEnsureCapacity_real_T(BB, i);
    GCmat_data = BB->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + BB->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    /* 'gkmPWMlasso:252' BX = B(:,f)'*B(:,f); */
    match_idx = c_f->size[0] * c_f->size[1];
    nx = c_f->size[0] * c_f->size[1];
    loop_ub = B->size[0];
    i = xc->size[0] * xc->size[1];
    xc->size[0] = B->size[0];
    xc->size[1] = match_idx;
    emxEnsureCapacity_real_T(xc, i);
    A_data = xc->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + xc->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = nx;
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    b_mtimes(xc, comb, comb2);
    f_data = comb2->data;
    /* 'gkmPWMlasso:253' BY = B(:,f)'*cfile2; */
    match_idx = c_f->size[0] * c_f->size[1];
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = match_idx;
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (c_f->size[0] * c_f->size[1] == 0) ||
        (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = c_f->size[0] * c_f->size[1];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = c_f->size[0] * c_f->size[1];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = c_f->size[0] * c_f->size[1];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)(c_f->size[0] * c_f->size[1]), (blasint)1,
                  (blasint)B->size[0], 1.0, &comb_data[0], (blasint)B->size[0],
                  &cfile2_data[0], (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)(c_f->size[0] * c_f->size[1]));
    }
    /* 'gkmPWMlasso:254' E = zeros(length(f),1); */
    if ((c_f->size[0] == 0) || (c_f->size[1] == 0)) {
      nx = 0;
    } else {
      match_idx = c_f->size[0];
      nx = c_f->size[1];
      if (match_idx >= nx) {
        nx = match_idx;
      }
    }
    i = normvec->size[0];
    normvec->size[0] = nx;
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    for (i = 0; i < nx; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso:255' for i = 1:length(f) */
    if ((c_f->size[0] == 0) || (c_f->size[1] == 0)) {
      nx = 0;
    } else {
      match_idx = c_f->size[0];
      nx = c_f->size[1];
      if (match_idx >= nx) {
        nx = match_idx;
      }
      j_loop_ub = comb2->size[0] * comb2->size[1];
      k_loop_ub = BY->size[0];
      nxout = BY->size[0] - 1;
      l_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < nx; b_i++) {
      /* 'gkmPWMlasso:256' B = BB; */
      /* 'gkmPWMlasso:257' BBX = BX; */
      /* 'gkmPWMlasso:258' BBY = BY; */
      /* 'gkmPWMlasso:259' B(:,i) = []; */
      /* 'gkmPWMlasso:260' BBX(:,i) = []; */
      /* 'gkmPWMlasso:261' BBX(i,:) = []; */
      /* 'gkmPWMlasso:262' BBY(i) = []; */
      /* 'gkmPWMlasso:263' res = cfile2-B*(BBX^-1*BBY); */
      i = comb->size[0] * comb->size[1];
      comb->size[0] = comb2->size[0];
      comb->size[1] = comb2->size[1];
      emxEnsureCapacity_real_T(comb, i);
      comb_data = comb->data;
      for (i = 0; i < j_loop_ub; i++) {
        comb_data[i] = f_data[i];
      }
      d_nullAssignment(comb, b_i + 1);
      b_nullAssignment(comb, b_i + 1);
      mpower(comb, xc);
      A_data = xc->data;
      b_loop_ub = b_i + 1;
      i = loc->size[0];
      loc->size[0] = BY->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      for (i = 0; i < k_loop_ub; i++) {
        loc_data[i] = BY_data[i];
      }
      for (k = b_loop_ub; k <= nxout; k++) {
        loc_data[k - 1] = loc_data[k];
      }
      i = loc->size[0];
      if (1 > BY->size[0] - 1) {
        loc->size[0] = 0;
      } else {
        loc->size[0] = BY->size[0] - 1;
      }
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      if ((xc->size[0] == 0) || (xc->size[1] == 0) || (loc->size[0] == 0)) {
        i = indc->size[0];
        indc->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        loop_ub = xc->size[0];
        for (i = 0; i < loop_ub; i++) {
          indc_data[i] = 0.0;
        }
      } else {
        i = indc->size[0];
        indc->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                    &A_data[0], (blasint)xc->size[0], &loc_data[0],
                    (blasint)loc->size[0], 0.0, &indc_data[0],
                    (blasint)xc->size[0]);
      }
      i = comb->size[0] * comb->size[1];
      comb->size[0] = BB->size[0];
      comb->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb, i);
      comb_data = comb->data;
      for (i = 0; i < l_loop_ub; i++) {
        comb_data[i] = GCmat_data[i];
      }
      d_nullAssignment(comb, b_i + 1);
      comb_data = comb->data;
      if ((comb->size[0] == 0) || (comb->size[1] == 0) ||
          (indc->size[0] == 0)) {
        i = loc->size[0];
        loc->size[0] = comb->size[0];
        emxEnsureCapacity_real_T(loc, i);
        loc_data = loc->data;
        loop_ub = comb->size[0];
        for (i = 0; i < loop_ub; i++) {
          loc_data[i] = 0.0;
        }
      } else {
        i = loc->size[0];
        loc->size[0] = comb->size[0];
        emxEnsureCapacity_real_T(loc, i);
        loc_data = loc->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)comb->size[0], (blasint)1, (blasint)comb->size[1],
                    1.0, &comb_data[0], (blasint)comb->size[0], &indc_data[0],
                    (blasint)indc->size[0], 0.0, &loc_data[0],
                    (blasint)comb->size[0]);
      }
      loop_ub = cfile2->size[0];
      if (cfile2->size[0] == loc->size[0]) {
        i = loc->size[0];
        loc->size[0] = cfile2->size[0];
        emxEnsureCapacity_real_T(loc, i);
        loc_data = loc->data;
        for (i = 0; i < loop_ub; i++) {
          loc_data[i] = cfile2_data[i] - loc_data[i];
        }
      } else {
        minus(loc, cfile2);
        loc_data = loc->data;
      }
      /* 'gkmPWMlasso:264' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        GCpos1 = 0.0;
      } else {
        GCpos1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                            &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(GCpos1);
    }
    /* 'gkmPWMlasso:266' res = cfile2-BB*OLS; */
    if ((B->size[0] == 0) || (c_f->size[0] * c_f->size[1] == 0) ||
        (OLS->size[0] == 0)) {
      loop_ub = B->size[0];
      i = loc->size[0];
      loc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = 0.0;
      }
    } else {
      i = loc->size[0];
      loc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)B->size[0], (blasint)1,
                  (blasint)(c_f->size[0] * c_f->size[1]), 1.0, &GCmat_data[0],
                  (blasint)B->size[0], &OLS_data[0], (blasint)OLS->size[0], 0.0,
                  &loc_data[0], (blasint)B->size[0]);
    }
    if (cfile2->size[0] == loc->size[0]) {
      i = loc->size[0];
      loc->size[0] = cfile2->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = cfile2->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = cfile2_data[i] - loc_data[i];
      }
    } else {
      minus(loc, cfile2);
      loc_data = loc->data;
    }
    /* 'gkmPWMlasso:267' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      GCpos1 = 0.0;
    } else {
      GCpos1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                          &loc_data[0], (blasint)1);
    }
    nfrac = sqrt(GCpos1);
    /* 'gkmPWMlasso:268' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - nfrac) / nfrac;
    }
    /*  motclus = motclus(f); */
    /* 'gkmPWMlasso:270' mylen = length(f); */
    if ((c_f->size[0] == 0) || (c_f->size[1] == 0)) {
      nx = 0;
    } else {
      match_idx = c_f->size[0];
      nx = c_f->size[1];
      if (match_idx >= nx) {
        nx = match_idx;
      }
    }
    /* 'gkmPWMlasso:271' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso:272' for idx=1:mylen */
    i = c_motclus->size[0];
    c_motclus->size[0] = nx;
    emxEnsureCapacity_cell_wrap_1(c_motclus, i);
    b_motclus_data = c_motclus->data;
    for (b_loop_ub = 0; b_loop_ub < nx; b_loop_ub++) {
      /* 'gkmPWMlasso:273' newMotclus{idx} = motclus{f(idx)}; */
      i = b_motclus_data[b_loop_ub].f1->size[0] *
          b_motclus_data[b_loop_ub].f1->size[1];
      b_motclus_data[b_loop_ub].f1->size[0] = 1;
      b_motclus_data[b_loop_ub].f1->size[1] =
          motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->size[1];
      emxEnsureCapacity_real_T(b_motclus_data[b_loop_ub].f1, i);
      loop_ub = motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->size[1];
      for (i = 0; i < loop_ub; i++) {
        b_motclus_data[b_loop_ub].f1->data[i] =
            motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->data[i];
      }
    }
    emxFree_real_T(&c_f);
    /* 'gkmPWMlasso:275' motclus = newMotclus; */
    /* 'gkmPWMlasso:277' f = find(E/max(E) >= 0.01); */
    nfrac = maximum(normvec);
    i = b_loc->size[0];
    b_loc->size[0] = normvec->size[0];
    emxEnsureCapacity_boolean_T(b_loc, i);
    b_loc_data = b_loc->data;
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_loc_data[i] = (normvec_data[i] / nfrac >= 0.01);
    }
    b_eml_find(b_loc, idx);
    matches_data = idx->data;
    i = b_f->size[0];
    b_f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(b_f, i);
    b_f_data = b_f->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_f_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso:278' B = BB; */
    i = B->size[0] * B->size[1];
    B->size[0] = BB->size[0];
    B->size[1] = BB->size[1];
    emxEnsureCapacity_real_T(B, i);
    B_data = B->data;
    loop_ub = BB->size[0] * BB->size[1];
    for (i = 0; i < loop_ub; i++) {
      B_data[i] = GCmat_data[i];
    }
    /* 'gkmPWMlasso:279' BB = B(:,f); */
    match_idx = BB->size[0] - 1;
    i = comb->size[0] * comb->size[1];
    comb->size[0] = BB->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    loop_ub = b_f->size[0];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 <= match_idx; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            GCmat_data[i1 + BB->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    i = BB->size[0] * BB->size[1];
    BB->size[0] = comb->size[0];
    BB->size[1] = comb->size[1];
    emxEnsureCapacity_real_T(BB, i);
    GCmat_data = BB->data;
    loop_ub = comb->size[0] * comb->size[1];
    for (i = 0; i < loop_ub; i++) {
      GCmat_data[i] = comb_data[i];
    }
    /* 'gkmPWMlasso:280' BX = B(:,f)'*B(:,f); */
    loop_ub = B->size[0];
    i = xc->size[0] * xc->size[1];
    xc->size[0] = B->size[0];
    xc->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(xc, i);
    A_data = xc->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + xc->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    b_mtimes(xc, comb, comb2);
    f_data = comb2->data;
    /* 'gkmPWMlasso:281' BY = B(:,f)'*cfile2; */
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (b_f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)b_f->size[0], (blasint)1, (blasint)B->size[0], 1.0,
                  &comb_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)b_f->size[0]);
    }
    /* 'gkmPWMlasso:282' E = zeros(length(f),1); */
    i = normvec->size[0];
    normvec->size[0] = b_f->size[0];
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    loop_ub = b_f->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso:283' for i = 1:length(f) */
    i = b_f->size[0];
    if (0 <= b_f->size[0] - 1) {
      m_loop_ub = comb2->size[0] * comb2->size[1];
      n_loop_ub = BY->size[0];
      b_nxout = BY->size[0] - 1;
      o_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso:284' B = BB; */
      /* 'gkmPWMlasso:285' BBX = BX; */
      /* 'gkmPWMlasso:286' BBY = BY; */
      /* 'gkmPWMlasso:287' B(:,i) = []; */
      /* 'gkmPWMlasso:288' BBX(:,i) = []; */
      /* 'gkmPWMlasso:289' BBX(i,:) = []; */
      /* 'gkmPWMlasso:290' BBY(i) = []; */
      /* 'gkmPWMlasso:291' res = cfile2-B*(BBX^-1*BBY); */
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = comb2->size[0];
      comb->size[1] = comb2->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      comb_data = comb->data;
      for (i1 = 0; i1 < m_loop_ub; i1++) {
        comb_data[i1] = f_data[i1];
      }
      d_nullAssignment(comb, b_i + 1);
      b_nullAssignment(comb, b_i + 1);
      mpower(comb, xc);
      A_data = xc->data;
      b_loop_ub = b_i + 1;
      i1 = loc->size[0];
      loc->size[0] = BY->size[0];
      emxEnsureCapacity_real_T(loc, i1);
      loc_data = loc->data;
      for (i1 = 0; i1 < n_loop_ub; i1++) {
        loc_data[i1] = BY_data[i1];
      }
      for (k = b_loop_ub; k <= b_nxout; k++) {
        loc_data[k - 1] = loc_data[k];
      }
      i1 = loc->size[0];
      if (1 > BY->size[0] - 1) {
        loc->size[0] = 0;
      } else {
        loc->size[0] = BY->size[0] - 1;
      }
      emxEnsureCapacity_real_T(loc, i1);
      loc_data = loc->data;
      if ((xc->size[0] == 0) || (xc->size[1] == 0) || (loc->size[0] == 0)) {
        i1 = indc->size[0];
        indc->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(indc, i1);
        indc_data = indc->data;
        loop_ub = xc->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          indc_data[i1] = 0.0;
        }
      } else {
        i1 = indc->size[0];
        indc->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(indc, i1);
        indc_data = indc->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                    &A_data[0], (blasint)xc->size[0], &loc_data[0],
                    (blasint)loc->size[0], 0.0, &indc_data[0],
                    (blasint)xc->size[0]);
      }
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = BB->size[0];
      comb->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      comb_data = comb->data;
      for (i1 = 0; i1 < o_loop_ub; i1++) {
        comb_data[i1] = GCmat_data[i1];
      }
      d_nullAssignment(comb, b_i + 1);
      comb_data = comb->data;
      if ((comb->size[0] == 0) || (comb->size[1] == 0) ||
          (indc->size[0] == 0)) {
        i1 = loc->size[0];
        loc->size[0] = comb->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        loop_ub = comb->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          loc_data[i1] = 0.0;
        }
      } else {
        i1 = loc->size[0];
        loc->size[0] = comb->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)comb->size[0], (blasint)1, (blasint)comb->size[1],
                    1.0, &comb_data[0], (blasint)comb->size[0], &indc_data[0],
                    (blasint)indc->size[0], 0.0, &loc_data[0],
                    (blasint)comb->size[0]);
      }
      loop_ub = cfile2->size[0];
      if (cfile2->size[0] == loc->size[0]) {
        i1 = loc->size[0];
        loc->size[0] = cfile2->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        for (i1 = 0; i1 < loop_ub; i1++) {
          loc_data[i1] = cfile2_data[i1] - loc_data[i1];
        }
      } else {
        minus(loc, cfile2);
        loc_data = loc->data;
      }
      /* 'gkmPWMlasso:292' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        GCpos1 = 0.0;
      } else {
        GCpos1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                            &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(GCpos1);
    }
    /* 'gkmPWMlasso:294' OLS = (BB.'*BB)^-1*(BB.'*cfile2); */
    b_mtimes(BB, BB, comb2);
    if ((BB->size[0] == 0) || (BB->size[1] == 0) || (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = BB->size[1];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = BB->size[1];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = BB->size[1];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)BB->size[1],
                  (blasint)1, (blasint)BB->size[0], 1.0, &GCmat_data[0],
                  (blasint)BB->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)BB->size[1]);
    }
    mpower(comb2, xc);
    A_data = xc->data;
    if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
      i = OLS->size[0];
      OLS->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      loop_ub = xc->size[0];
      for (i = 0; i < loop_ub; i++) {
        OLS_data[i] = 0.0;
      }
    } else {
      i = OLS->size[0];
      OLS->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                  &A_data[0], (blasint)xc->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &OLS_data[0],
                  (blasint)xc->size[0]);
    }
    /* 'gkmPWMlasso:295' Pweig = Pweig(f); */
    /* 'gkmPWMlasso:296' res = cfile2-BB*OLS; */
    if ((BB->size[0] == 0) || (BB->size[1] == 0) || (OLS->size[0] == 0)) {
      i = loc->size[0];
      loc->size[0] = BB->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = BB->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = 0.0;
      }
    } else {
      i = loc->size[0];
      loc->size[0] = BB->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)BB->size[0], (blasint)1, (blasint)BB->size[1], 1.0,
                  &GCmat_data[0], (blasint)BB->size[0], &OLS_data[0],
                  (blasint)OLS->size[0], 0.0, &loc_data[0],
                  (blasint)BB->size[0]);
    }
    if (cfile2->size[0] == loc->size[0]) {
      i = loc->size[0];
      loc->size[0] = cfile2->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = cfile2->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = cfile2_data[i] - loc_data[i];
      }
    } else {
      minus(loc, cfile2);
      loc_data = loc->data;
    }
    /* 'gkmPWMlasso:297' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      GCpos1 = 0.0;
    } else {
      GCpos1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                          &loc_data[0], (blasint)1);
    }
    nfrac = sqrt(GCpos1);
    /* 'gkmPWMlasso:298' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - nfrac) / nfrac;
    }
    /* 'gkmPWMlasso:299' for i = 1:length(motclus) */
    i = c_motclus->size[0];
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso:300' motclus{i} = indvec(motclus{i})'; */
      i1 = b_motclus_data[b_i].f1->size[0] * b_motclus_data[b_i].f1->size[1];
      b_motclus_data[b_i].f1->size[0] = 1;
      emxEnsureCapacity_real_T(b_motclus_data[b_i].f1, i1);
      loop_ub = b_motclus_data[b_i].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_motclus_data[b_i].f1->data[i1] =
            indvec_data[(int)b_motclus_data[b_i].f1->data[i1] - 1];
      }
    }
    /*  motclus = motclus(f); */
    /* 'gkmPWMlasso:303' mylen = length(f); */
    /* 'gkmPWMlasso:304' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso:305' for idx=1:mylen */
    i = b_f->size[0];
    i1 = tmp_motclus->size[0];
    tmp_motclus->size[0] = b_f->size[0];
    emxEnsureCapacity_cell_wrap_1(tmp_motclus, i1);
    tmp_motclus_data = tmp_motclus->data;
    for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
      /* 'gkmPWMlasso:306' newMotclus{idx} = motclus{f(idx)}; */
      i1 = tmp_motclus_data[b_loop_ub].f1->size[0] *
           tmp_motclus_data[b_loop_ub].f1->size[1];
      tmp_motclus_data[b_loop_ub].f1->size[0] = 1;
      tmp_motclus_data[b_loop_ub].f1->size[1] =
          b_motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->size[1];
      emxEnsureCapacity_real_T(tmp_motclus_data[b_loop_ub].f1, i1);
      loop_ub = b_motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        tmp_motclus_data[b_loop_ub].f1->data[i1] =
            b_motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->data[i1];
      }
    }
    /* 'gkmPWMlasso:308' motclus = newMotclus; */
    /* 'gkmPWMlasso:310' correlation = corrcoef(cfile2, BB*OLS); */
    if ((BB->size[0] == 0) || (BB->size[1] == 0) || (OLS->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = BB->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = BB->size[0];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = BB->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)BB->size[0], (blasint)1, (blasint)BB->size[1], 1.0,
                  &GCmat_data[0], (blasint)BB->size[0], &OLS_data[0],
                  (blasint)OLS->size[0], 0.0, &indc_data[0],
                  (blasint)BB->size[0]);
    }
    /* 'gkmPWMlasso:311' gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus,
     * sprintf("%s_%d_%d", filename, int32(l_svm2), int32(k_svm2)),
     * memefile,num,minL, minInfo, correlation(1,2)); */
    nfrac = maximum(OLS);
    GCpos1 = maximum(normvec);
    d = rt_roundd(varargin_7);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        i = (int)d;
      } else {
        i = MIN_int32_T;
      }
    } else {
      i = MAX_int32_T;
    }
    d = rt_roundd(varargin_8);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        i1 = (int)d;
      } else {
        i1 = MIN_int32_T;
      }
    } else {
      i1 = MAX_int32_T;
    }
    b_sprintf(varargin_1, i, i1, text);
    loop_ub = OLS->size[0];
    for (i = 0; i < loop_ub; i++) {
      OLS_data[i] /= nfrac;
    }
    i = BY->size[0];
    BY->size[0] = b_f->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    loop_ub = b_f->size[0];
    for (i = 0; i < loop_ub; i++) {
      BY_data[i] = Pweig_data[(int)b_f_data[i] - 1];
    }
    emxFree_real_T(&Pweig);
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] /= GCpos1;
    }
    corrcoef(cfile2, indc, dv);
    gettopmotifs(OLS, BY, normvec, tmp_motclus, text, varargin_2,
                 size_tmp_idx_1, varargin_4, varargin_5, dv[2]);
  } else {
    /* 'gkmPWMlasso:313' else */
    /* 'gkmPWMlasso:315' fprintf('Running LASSO (1)\n'); */
    printf("Running LASSO (1)\n");
    fflush(stdout);
    /*  weigmat = lasso_cvmat(B, cfile2,'DFmax', d,'Standardize', false); */
    /* 'gkmPWMlasso:317' [weigmat, ~] = lasso_cvmat(B, cfile2, d, false, 100);
     */
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = B->size[1];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    loop_ub = B->size[0] * B->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      comb_data[i] = B_data[i];
    }
    i = diffc->size[0];
    diffc->size[0] = cfile2->size[0];
    emxEnsureCapacity_real_T(diffc, i);
    GCmat_data = diffc->data;
    loop_ub = cfile2->size[0] - 1;
    for (i = 0; i <= loop_ub; i++) {
      GCmat_data[i] = cfile2_data[i];
    }
    emxInit_real_T(&weigmat, 2);
    b_lasso_cvmat(comb, diffc, varargin_3, weigmat, expl_temp_data, lk_size,
                  b_expl_temp_data, expl_temp_size, &nfrac, a__3_DF_data,
                  a__3_DF_size, c_expl_temp_data, b_expl_temp_size);
    GCmat_data = weigmat->data;
    /* 'gkmPWMlasso:319' f = find(weigmat(:,1)~=0); */
    loop_ub = weigmat->size[0];
    i = b_loc->size[0];
    b_loc->size[0] = weigmat->size[0];
    emxEnsureCapacity_boolean_T(b_loc, i);
    b_loc_data = b_loc->data;
    for (i = 0; i < loop_ub; i++) {
      b_loc_data[i] = (GCmat_data[i] != 0.0);
    }
    b_eml_find(b_loc, idx);
    matches_data = idx->data;
    i = b_f->size[0];
    b_f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(b_f, i);
    b_f_data = b_f->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_f_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso:320' cfile3 = cfile2 -
     * B(:,f)*((B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2)); */
    loop_ub = B->size[0];
    i = xc->size[0] * xc->size[1];
    xc->size[0] = B->size[0];
    xc->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(xc, i);
    A_data = xc->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + xc->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    b_mtimes(xc, comb, comb2);
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (b_f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)b_f->size[0], (blasint)1, (blasint)B->size[0], 1.0,
                  &comb_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)b_f->size[0]);
    }
    mpower(comb2, xc);
    A_data = xc->data;
    if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = xc->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                  &A_data[0], (blasint)xc->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &BY_data[0],
                  (blasint)xc->size[0]);
    }
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (b_f->size[0] == 0) || (BY->size[0] == 0)) {
      loop_ub = B->size[0];
      i = loc->size[0];
      loc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = 0.0;
      }
    } else {
      i = loc->size[0];
      loc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)B->size[0], (blasint)1, (blasint)b_f->size[0], 1.0,
                  &comb_data[0], (blasint)B->size[0], &BY_data[0],
                  (blasint)BY->size[0], 0.0, &loc_data[0], (blasint)B->size[0]);
    }
    /* 'gkmPWMlasso:322' fprintf('Running LASSO (2)\n'); */
    printf("Running LASSO (2)\n");
    fflush(stdout);
    /*  weigmat2 = lasso_cvmat(B, cfile3,'DFmax', d,'Standardize', false); */
    /* 'gkmPWMlasso:324' [weigmat2, ~] = lasso_cvmat(B, cfile3, d, false, 100);
     */
    emxInit_real_T(&c_f, 2);
    if (cfile2->size[0] == loc->size[0]) {
      i = loc->size[0];
      loc->size[0] = cfile2->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = cfile2->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = cfile2_data[i] - loc_data[i];
      }
      i = comb->size[0] * comb->size[1];
      comb->size[0] = B->size[0];
      comb->size[1] = B->size[1];
      emxEnsureCapacity_real_T(comb, i);
      comb_data = comb->data;
      loop_ub = B->size[0] * B->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        comb_data[i] = B_data[i];
      }
      b_lasso_cvmat(comb, loc, varargin_3, c_f, d_expl_temp_data, lk_size,
                    e_expl_temp_data, expl_temp_size, &nfrac, a__3_DF_data,
                    a__3_DF_size, f_expl_temp_data, b_expl_temp_size);
      f_data = c_f->data;
      loc_data = loc->data;
    } else {
      c_binary_expand_op(B, cfile2, loc, varargin_3, c_f, d_expl_temp_data,
                         lk_size, e_expl_temp_data, expl_temp_size,
                         a__3_DF_data, a__3_DF_size, f_expl_temp_data,
                         b_expl_temp_size);
      f_data = c_f->data;
    }
    /* 'gkmPWMlasso:326' weigmat = abs(weigmat(:,1)) + abs(weigmat2(:,1)); */
    /* 'gkmPWMlasso:327' f = find(weigmat~=0); */
    nx = weigmat->size[0] - 1;
    i = indc->size[0];
    indc->size[0] = weigmat->size[0];
    emxEnsureCapacity_real_T(indc, i);
    indc_data = indc->data;
    for (k = 0; k <= nx; k++) {
      indc_data[k] = fabs(GCmat_data[k]);
    }
    emxFree_real_T(&weigmat);
    nx = c_f->size[0] - 1;
    i = BY->size[0];
    BY->size[0] = c_f->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    for (k = 0; k <= nx; k++) {
      BY_data[k] = fabs(f_data[k]);
    }
    emxFree_real_T(&c_f);
    if (indc->size[0] == BY->size[0]) {
      i = b_loc->size[0];
      b_loc->size[0] = indc->size[0];
      emxEnsureCapacity_boolean_T(b_loc, i);
      b_loc_data = b_loc->data;
      loop_ub = indc->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_loc_data[i] = (indc_data[i] + BY_data[i] != 0.0);
      }
      b_eml_find(b_loc, idx);
      matches_data = idx->data;
    } else {
      b_binary_expand_op(idx, indc, BY);
      matches_data = idx->data;
    }
    i = b_f->size[0];
    b_f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(b_f, i);
    b_f_data = b_f->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_f_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso:328' motclus2 = clus_simmat_eig(corrcoef(B(:,f)),0.6); */
    loop_ub = B->size[0];
    i = xc->size[0] * xc->size[1];
    xc->size[0] = B->size[0];
    xc->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(xc, i);
    A_data = xc->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + xc->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    b_corrcoef(xc, comb);
    b_clus_simmat_eig(comb, tmp_motclus);
    tmp_motclus_data = tmp_motclus->data;
    /* 'gkmPWMlasso:329' if length(length(motclus2)) ~= length(f) */
    if (1 != b_f->size[0]) {
      /* 'gkmPWMlasso:330' f2 = f; */
      i = loc->size[0];
      loc->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = b_f_data[i];
      }
      /* 'gkmPWMlasso:331' f = zeros(length(motclus2),1); */
      i = b_f->size[0];
      b_f->size[0] = tmp_motclus->size[0];
      emxEnsureCapacity_real_T(b_f, i);
      b_f_data = b_f->data;
      loop_ub = tmp_motclus->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_f_data[i] = 0.0;
      }
      /* 'gkmPWMlasso:332' for i = 1:length(motclus2) */
      i = tmp_motclus->size[0];
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso:333' if length(motclus2{i}) > 1 */
        if (tmp_motclus_data[b_i].f1->size[1] > 1) {
          /* 'gkmPWMlasso:334' [~,f3] = max(abs(Z(f2(motclus2{i})))); */
          nx = tmp_motclus_data[b_i].f1->size[1];
          i1 = BY->size[0];
          BY->size[0] = tmp_motclus_data[b_i].f1->size[1];
          emxEnsureCapacity_real_T(BY, i1);
          BY_data = BY->data;
          for (k = 0; k < nx; k++) {
            BY_data[k] = fabs(
                negvec_data
                    [(int)loc_data[(int)tmp_motclus_data[b_i].f1->data[k] - 1] -
                     1]);
          }
          d_maximum(BY, &nfrac, &match_idx);
          /* 'gkmPWMlasso:335' f(i) = f2(motclus2{i}(f3)); */
          b_f_data[b_i] =
              loc_data[(int)tmp_motclus_data[b_i].f1->data[match_idx - 1] - 1];
        } else {
          /* 'gkmPWMlasso:336' else */
          /* 'gkmPWMlasso:337' f(i) = f2(motclus2{i}(1)); */
          b_f_data[b_i] = loc_data[(int)tmp_motclus_data[b_i].f1->data[0] - 1];
        }
      }
    }
    /* 'gkmPWMlasso:341' [a,b] = sort(abs(Z(f)), 'descend'); */
    nx = b_f->size[0];
    i = indc->size[0];
    indc->size[0] = b_f->size[0];
    emxEnsureCapacity_real_T(indc, i);
    indc_data = indc->data;
    for (k = 0; k < nx; k++) {
      indc_data[k] = fabs(negvec_data[(int)b_f_data[k] - 1]);
    }
    sort(indc, idx);
    matches_data = idx->data;
    /* 'gkmPWMlasso:342' f = f(b); */
    i = BY->size[0];
    BY->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      BY_data[i] = b_f_data[matches_data[i] - 1];
    }
    i = b_f->size[0];
    b_f->size[0] = BY->size[0];
    emxEnsureCapacity_real_T(b_f, i);
    b_f_data = b_f->data;
    loop_ub = BY->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_f_data[i] = BY_data[i];
    }
    /* 'gkmPWMlasso:343' F = length(f); */
    nfrac = b_f->size[0];
    /* 'gkmPWMlasso:344' fprintf('Selecting top motifs\n') */
    printf("Selecting top motifs\n");
    fflush(stdout);
    /* 'gkmPWMlasso:345' ind = true; */
    /* 'gkmPWMlasso:346' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    loop_ub = B->size[0];
    i = xc->size[0] * xc->size[1];
    xc->size[0] = B->size[0];
    xc->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(xc, i);
    A_data = xc->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + xc->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    b_mtimes(xc, comb, comb2);
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (b_f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)b_f->size[0], (blasint)1, (blasint)B->size[0], 1.0,
                  &comb_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)b_f->size[0]);
    }
    mpower(comb2, xc);
    A_data = xc->data;
    if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
      i = OLS->size[0];
      OLS->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      loop_ub = xc->size[0];
      for (i = 0; i < loop_ub; i++) {
        OLS_data[i] = 0.0;
      }
    } else {
      i = OLS->size[0];
      OLS->size[0] = xc->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                  &A_data[0], (blasint)xc->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &OLS_data[0],
                  (blasint)xc->size[0]);
    }
    /* 'gkmPWMlasso:347' Pweig = Z(f); */
    i = diffc->size[0];
    diffc->size[0] = b_f->size[0];
    emxEnsureCapacity_real_T(diffc, i);
    GCmat_data = diffc->data;
    loop_ub = b_f->size[0];
    for (i = 0; i < loop_ub; i++) {
      GCmat_data[i] = negvec_data[(int)b_f_data[i] - 1];
    }
    /* 'gkmPWMlasso:348' while ind && F >= d */
    exitg1 = false;
    while ((!exitg1) && (nfrac >= varargin_3)) {
      /* 'gkmPWMlasso:349' OLS =
       * (B(:,f(1:d)).'*B(:,f(1:d)))^-1*(B(:,f(1:d)).'*cfile2); */
      if (1.0 > varargin_3) {
        loop_ub = 0;
        b_loop_ub = 0;
        match_idx = 0;
      } else {
        loop_ub = (int)varargin_3;
        b_loop_ub = (int)varargin_3;
        match_idx = (int)varargin_3;
      }
      nx = B->size[0];
      i = xc->size[0] * xc->size[1];
      xc->size[0] = B->size[0];
      xc->size[1] = loop_ub;
      emxEnsureCapacity_real_T(xc, i);
      A_data = xc->data;
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < nx; i1++) {
          A_data[i1 + xc->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
        }
      }
      loop_ub = B->size[0];
      i = comb->size[0] * comb->size[1];
      comb->size[0] = B->size[0];
      comb->size[1] = b_loop_ub;
      emxEnsureCapacity_real_T(comb, i);
      comb_data = comb->data;
      for (i = 0; i < b_loop_ub; i++) {
        for (i1 = 0; i1 < loop_ub; i1++) {
          comb_data[i1 + comb->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
        }
      }
      b_mtimes(xc, comb, comb2);
      loop_ub = B->size[0];
      i = comb->size[0] * comb->size[1];
      comb->size[0] = B->size[0];
      comb->size[1] = match_idx;
      emxEnsureCapacity_real_T(comb, i);
      comb_data = comb->data;
      for (i = 0; i < match_idx; i++) {
        for (i1 = 0; i1 < loop_ub; i1++) {
          comb_data[i1 + comb->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
        }
      }
      if ((B->size[0] == 0) || (match_idx == 0) || (cfile2->size[0] == 0)) {
        i = indc->size[0];
        indc->size[0] = match_idx;
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        for (i = 0; i < match_idx; i++) {
          indc_data[i] = 0.0;
        }
      } else {
        i = indc->size[0];
        indc->size[0] = match_idx;
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)match_idx,
                    (blasint)1, (blasint)B->size[0], 1.0, &comb_data[0],
                    (blasint)B->size[0], &cfile2_data[0],
                    (blasint)cfile2->size[0], 0.0, &indc_data[0],
                    (blasint)match_idx);
      }
      mpower(comb2, xc);
      A_data = xc->data;
      if ((xc->size[0] == 0) || (xc->size[1] == 0) || (indc->size[0] == 0)) {
        i = OLS->size[0];
        OLS->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(OLS, i);
        OLS_data = OLS->data;
        loop_ub = xc->size[0];
        for (i = 0; i < loop_ub; i++) {
          OLS_data[i] = 0.0;
        }
      } else {
        i = OLS->size[0];
        OLS->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(OLS, i);
        OLS_data = OLS->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                    &A_data[0], (blasint)xc->size[0], &indc_data[0],
                    (blasint)indc->size[0], 0.0, &OLS_data[0],
                    (blasint)xc->size[0]);
      }
      /* 'gkmPWMlasso:350' Pweig = Z(f(1:d)); */
      if (1.0 > varargin_3) {
        loop_ub = 0;
      } else {
        loop_ub = (int)varargin_3;
      }
      i = diffc->size[0];
      diffc->size[0] = loop_ub;
      emxEnsureCapacity_real_T(diffc, i);
      GCmat_data = diffc->data;
      for (i = 0; i < loop_ub; i++) {
        GCmat_data[i] = negvec_data[(int)b_f_data[i] - 1];
      }
      /* 'gkmPWMlasso:351' ff = []; */
      loc->size[0] = 0;
      /* 'gkmPWMlasso:352' for i = 1:d */
      i = (int)varargin_3;
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso:353' if sign(OLS(i)) ~= sign(Pweig(i)) */
        a = OLS_data[b_i];
        if (OLS_data[b_i] < 0.0) {
          a = -1.0;
        } else if (OLS_data[b_i] > 0.0) {
          a = 1.0;
        }
        GCpos1 = negvec_data[(int)b_f_data[b_i] - 1];
        if (GCpos1 < 0.0) {
          GCpos1 = -1.0;
        } else if (GCpos1 > 0.0) {
          GCpos1 = 1.0;
        }
        if (a != GCpos1) {
          /* 'gkmPWMlasso:354' ff = [ff;i]; */
          i1 = loc->size[0];
          i2 = loc->size[0];
          loc->size[0]++;
          emxEnsureCapacity_real_T(loc, i2);
          loc_data = loc->data;
          loc_data[i1] = (double)b_i + 1.0;
        }
      }
      /* 'gkmPWMlasso:357' if length(ff) > 0 && F - length(ff) >= d */
      if ((loc->size[0] > 0) && (nfrac - (double)loc->size[0] >= varargin_3)) {
        /* 'gkmPWMlasso:358' f(ff) = []; */
        i = idx->size[0];
        idx->size[0] = loc->size[0];
        emxEnsureCapacity_int32_T(idx, i);
        matches_data = idx->data;
        loop_ub = loc->size[0];
        for (i = 0; i < loop_ub; i++) {
          matches_data[i] = (int)loc_data[i];
        }
        e_nullAssignment(b_f, idx);
        b_f_data = b_f->data;
        /* 'gkmPWMlasso:359' F = F - length(ff); */
        nfrac -= (double)loc->size[0];
      } else {
        /* 'gkmPWMlasso:360' else */
        /* 'gkmPWMlasso:361' ind = false; */
        /* 'gkmPWMlasso:362' f = f(1:d); */
        i = b_f->size[0];
        if (1.0 > varargin_3) {
          b_f->size[0] = 0;
        } else {
          b_f->size[0] = (int)varargin_3;
        }
        emxEnsureCapacity_real_T(b_f, i);
        b_f_data = b_f->data;
        exitg1 = true;
      }
    }
    /*  if F < d */
    /*      OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    /*      Pweig = Z(f); */
    /*  end */
    /* 'gkmPWMlasso:369' BB = B(:,f); */
    loop_ub = B->size[0];
    i = BB->size[0] * BB->size[1];
    BB->size[0] = B->size[0];
    BB->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(BB, i);
    GCmat_data = BB->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + BB->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    /* 'gkmPWMlasso:370' BX = B(:,f)'*B(:,f); */
    loop_ub = B->size[0];
    i = xc->size[0] * xc->size[1];
    xc->size[0] = B->size[0];
    xc->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(xc, i);
    A_data = xc->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        A_data[i1 + xc->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    b_mtimes(xc, comb, comb2);
    f_data = comb2->data;
    /* 'gkmPWMlasso:371' BY = B(:,f)'*cfile2; */
    loop_ub = B->size[0];
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = b_f->size[0];
    emxEnsureCapacity_real_T(comb, i);
    comb_data = comb->data;
    b_loop_ub = b_f->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        comb_data[i1 + comb->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)b_f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (b_f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = b_f->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = b_f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)b_f->size[0], (blasint)1, (blasint)B->size[0], 1.0,
                  &comb_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)b_f->size[0]);
    }
    /* 'gkmPWMlasso:372' E = zeros(length(f),1); */
    i = normvec->size[0];
    normvec->size[0] = b_f->size[0];
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    loop_ub = b_f->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso:373' for i = 1:length(f) */
    i = b_f->size[0];
    if (0 <= b_f->size[0] - 1) {
      e_loop_ub = comb2->size[0] * comb2->size[1];
      f_loop_ub = BY->size[0];
      nxout = BY->size[0] - 1;
      g_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso:374' B = BB; */
      /* 'gkmPWMlasso:375' BBX = BX; */
      /* 'gkmPWMlasso:376' BBY = BY; */
      /* 'gkmPWMlasso:377' B(:,i) = []; */
      /* 'gkmPWMlasso:378' BBX(:,i) = []; */
      /* 'gkmPWMlasso:379' BBX(i,:) = []; */
      /* 'gkmPWMlasso:380' BBY(i) = []; */
      /* 'gkmPWMlasso:381' res = cfile2-B*(BBX^-1*BBY); */
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = comb2->size[0];
      comb->size[1] = comb2->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      comb_data = comb->data;
      for (i1 = 0; i1 < e_loop_ub; i1++) {
        comb_data[i1] = f_data[i1];
      }
      d_nullAssignment(comb, b_i + 1);
      b_nullAssignment(comb, b_i + 1);
      mpower(comb, xc);
      A_data = xc->data;
      b_loop_ub = b_i + 1;
      i1 = loc->size[0];
      loc->size[0] = BY->size[0];
      emxEnsureCapacity_real_T(loc, i1);
      loc_data = loc->data;
      for (i1 = 0; i1 < f_loop_ub; i1++) {
        loc_data[i1] = BY_data[i1];
      }
      for (k = b_loop_ub; k <= nxout; k++) {
        loc_data[k - 1] = loc_data[k];
      }
      i1 = loc->size[0];
      if (1 > BY->size[0] - 1) {
        loc->size[0] = 0;
      } else {
        loc->size[0] = BY->size[0] - 1;
      }
      emxEnsureCapacity_real_T(loc, i1);
      loc_data = loc->data;
      if ((xc->size[0] == 0) || (xc->size[1] == 0) || (loc->size[0] == 0)) {
        i1 = indc->size[0];
        indc->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(indc, i1);
        indc_data = indc->data;
        loop_ub = xc->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          indc_data[i1] = 0.0;
        }
      } else {
        i1 = indc->size[0];
        indc->size[0] = xc->size[0];
        emxEnsureCapacity_real_T(indc, i1);
        indc_data = indc->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)xc->size[0], (blasint)1, (blasint)xc->size[1], 1.0,
                    &A_data[0], (blasint)xc->size[0], &loc_data[0],
                    (blasint)loc->size[0], 0.0, &indc_data[0],
                    (blasint)xc->size[0]);
      }
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = BB->size[0];
      comb->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      comb_data = comb->data;
      for (i1 = 0; i1 < g_loop_ub; i1++) {
        comb_data[i1] = GCmat_data[i1];
      }
      d_nullAssignment(comb, b_i + 1);
      comb_data = comb->data;
      if ((comb->size[0] == 0) || (comb->size[1] == 0) ||
          (indc->size[0] == 0)) {
        i1 = loc->size[0];
        loc->size[0] = comb->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        loop_ub = comb->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          loc_data[i1] = 0.0;
        }
      } else {
        i1 = loc->size[0];
        loc->size[0] = comb->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    (blasint)comb->size[0], (blasint)1, (blasint)comb->size[1],
                    1.0, &comb_data[0], (blasint)comb->size[0], &indc_data[0],
                    (blasint)indc->size[0], 0.0, &loc_data[0],
                    (blasint)comb->size[0]);
      }
      loop_ub = cfile2->size[0];
      if (cfile2->size[0] == loc->size[0]) {
        i1 = loc->size[0];
        loc->size[0] = cfile2->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        for (i1 = 0; i1 < loop_ub; i1++) {
          loc_data[i1] = cfile2_data[i1] - loc_data[i1];
        }
      } else {
        minus(loc, cfile2);
        loc_data = loc->data;
      }
      /* 'gkmPWMlasso:382' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        GCpos1 = 0.0;
      } else {
        GCpos1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                            &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(GCpos1);
    }
    /* 'gkmPWMlasso:384' res = cfile2-BB*OLS; */
    if ((B->size[0] == 0) || (b_f->size[0] == 0) || (OLS->size[0] == 0)) {
      loop_ub = B->size[0];
      i = loc->size[0];
      loc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = 0.0;
      }
    } else {
      i = loc->size[0];
      loc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)B->size[0], (blasint)1, (blasint)b_f->size[0], 1.0,
                  &GCmat_data[0], (blasint)B->size[0], &OLS_data[0],
                  (blasint)OLS->size[0], 0.0, &loc_data[0],
                  (blasint)B->size[0]);
    }
    if (cfile2->size[0] == loc->size[0]) {
      i = loc->size[0];
      loc->size[0] = cfile2->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = cfile2->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = cfile2_data[i] - loc_data[i];
      }
    } else {
      minus(loc, cfile2);
      loc_data = loc->data;
    }
    /* 'gkmPWMlasso:385' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      GCpos1 = 0.0;
    } else {
      GCpos1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                          &loc_data[0], (blasint)1);
    }
    nfrac = sqrt(GCpos1);
    /* 'gkmPWMlasso:386' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - nfrac) / nfrac;
    }
    /* 'gkmPWMlasso:387' for i = 1:length(motclus) */
    i = motclus->size[0];
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso:388' motclus{i} = indvec(motclus{i})'; */
      i1 = motclus_data[b_i].f1->size[0] * motclus_data[b_i].f1->size[1];
      motclus_data[b_i].f1->size[0] = 1;
      emxEnsureCapacity_real_T(motclus_data[b_i].f1, i1);
      loop_ub = motclus_data[b_i].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        motclus_data[b_i].f1->data[i1] =
            indvec_data[(int)motclus_data[b_i].f1->data[i1] - 1];
      }
    }
    /*  motclus = motclus(f); */
    /* 'gkmPWMlasso:391' mylen = length(f); */
    /* 'gkmPWMlasso:392' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso:393' for idx=1:mylen */
    i = b_f->size[0];
    i1 = c_motclus->size[0];
    c_motclus->size[0] = b_f->size[0];
    emxEnsureCapacity_cell_wrap_1(c_motclus, i1);
    b_motclus_data = c_motclus->data;
    for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
      /* 'gkmPWMlasso:394' newMotclus{idx} = motclus{f(idx)}; */
      i1 = b_motclus_data[b_loop_ub].f1->size[0] *
           b_motclus_data[b_loop_ub].f1->size[1];
      b_motclus_data[b_loop_ub].f1->size[0] = 1;
      b_motclus_data[b_loop_ub].f1->size[1] =
          motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->size[1];
      emxEnsureCapacity_real_T(b_motclus_data[b_loop_ub].f1, i1);
      loop_ub = motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_motclus_data[b_loop_ub].f1->data[i1] =
            motclus_data[(int)b_f_data[b_loop_ub] - 1].f1->data[i1];
      }
    }
    /* 'gkmPWMlasso:396' motclus = newMotclus; */
    /* 'gkmPWMlasso:398' correlation = corrcoef(cfile2, BB*OLS); */
    if ((B->size[0] == 0) || (b_f->size[0] == 0) || (OLS->size[0] == 0)) {
      loop_ub = B->size[0];
      i = indc->size[0];
      indc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = B->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)B->size[0], (blasint)1, (blasint)b_f->size[0], 1.0,
                  &GCmat_data[0], (blasint)B->size[0], &OLS_data[0],
                  (blasint)OLS->size[0], 0.0, &indc_data[0],
                  (blasint)B->size[0]);
    }
    /* 'gkmPWMlasso:399' gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus,
     * sprintf("%s_%d_%d_%d", filename, int32(l_svm2), int32(k_svm2), int32(d)),
     * memefile,num,minL, minInfo, correlation(1,2)); */
    nfrac = maximum(OLS);
    GCpos1 = maximum(normvec);
    d = rt_roundd(varargin_7);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        i = (int)d;
      } else {
        i = MIN_int32_T;
      }
    } else {
      i = MAX_int32_T;
    }
    d = rt_roundd(varargin_8);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        i1 = (int)d;
      } else {
        i1 = MIN_int32_T;
      }
    } else {
      i1 = MAX_int32_T;
    }
    d = rt_roundd(varargin_3);
    if (d < 2.147483648E+9) {
      if (d >= -2.147483648E+9) {
        i2 = (int)d;
      } else {
        i2 = MIN_int32_T;
      }
    } else {
      i2 = MAX_int32_T;
    }
    c_sprintf(varargin_1, i, i1, i2, text);
    loop_ub = OLS->size[0];
    for (i = 0; i < loop_ub; i++) {
      OLS_data[i] /= nfrac;
    }
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] /= GCpos1;
    }
    corrcoef(cfile2, indc, dv);
    gettopmotifs(OLS, diffc, normvec, c_motclus, text, varargin_2,
                 size_tmp_idx_1, varargin_4, varargin_5, dv[2]);
  }
  emxFree_boolean_T(&b_loc);
  emxFree_int32_T(&idx);
  emxFree_int32_T(&match_out);
  emxFree_char_T(&text);
  emxFree_real_T(&f);
  emxFree_cell_wrap_1(&c_motclus);
  emxFree_real_T(&loc);
  emxFree_real_T(&BY);
  emxFree_real_T(&BB);
  emxFree_real_T(&OLS);
  emxFree_cell_wrap_1(&tmp_motclus);
  emxFree_real_T(&b_f);
  emxFree_real_T(&B);
  emxFree_real_T(&cfile2);
  emxFree_cell_wrap_1(&motclus);
  emxFree_real_T(&normvec);
  emxFree_real_T(&indvec);
  emxFree_real_T(&negvec);
  emxFree_real_T(&xc);
  emxFree_real_T(&indc);
  emxFree_real_T(&diffc);
  emxFree_real_T(&comb2);
  emxFree_real_T(&comb);
  /* 'gkmPWMlasso:403' fprintf('Done\n'); */
  printf("Done\n");
  fflush(stdout);
}

void minus(emxArray_real_T *loc, const emxArray_real_T *cfile2)
{
  emxArray_real_T *b_cfile2;
  const double *cfile2_data;
  double *b_cfile2_data;
  double *loc_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  cfile2_data = cfile2->data;
  loc_data = loc->data;
  emxInit_real_T(&b_cfile2, 1);
  i = b_cfile2->size[0];
  if (loc->size[0] == 1) {
    b_cfile2->size[0] = cfile2->size[0];
  } else {
    b_cfile2->size[0] = loc->size[0];
  }
  emxEnsureCapacity_real_T(b_cfile2, i);
  b_cfile2_data = b_cfile2->data;
  stride_0_0 = (cfile2->size[0] != 1);
  stride_1_0 = (loc->size[0] != 1);
  if (loc->size[0] == 1) {
    loop_ub = cfile2->size[0];
  } else {
    loop_ub = loc->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_cfile2_data[i] = cfile2_data[i * stride_0_0] - loc_data[i * stride_1_0];
  }
  i = loc->size[0];
  loc->size[0] = b_cfile2->size[0];
  emxEnsureCapacity_real_T(loc, i);
  loc_data = loc->data;
  loop_ub = b_cfile2->size[0];
  for (i = 0; i < loop_ub; i++) {
    loc_data[i] = b_cfile2_data[i];
  }
  emxFree_real_T(&b_cfile2);
}

/* End of code generation (gkmPWMlasso.c) */
