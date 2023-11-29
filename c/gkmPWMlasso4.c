/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso4.c
 *
 * Code generation for function 'gkmPWMlasso4'
 *
 */

/* Include files */
#include "gkmPWMlasso4.h"
#include "BGkmer.h"
#include "PWM2kmers.h"
#include "PWM2kmers_norc.h"
#include "blockedSummation.h"
#include "conncomp.h"
#include "corrcoef.h"
#include "eml_setop.h"
#include "extremeKElements.h"
#include "fgetl.h"
#include "fileManager.h"
#include "fileread.h"
#include "find.h"
#include "genIndex.h"
#include "getgkmcounts.h"
#include "getmotif.h"
#include "gkmPWMlasso4_data.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_initialize.h"
#include "gkmPWMlasso4_rtwutil.h"
#include "gkmPWMlasso4_types.h"
#include "graph.h"
#include "lasso_cvmat.h"
#include "minOrMax.h"
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
#include "lapacke.h"
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
    const emxArray_real_T *loc, double varargin_10, emxArray_real_T *weigmat2,
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
  /* 'gkmPWMlasso4:462' bin = conncomp(graph(simmat>r)); */
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
  /* 'gkmPWMlasso4:463' n = max(bin); */
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
  /* 'gkmPWMlasso4:464' motclus = cell(n,1); */
  /* 'gkmPWMlasso4:465' for i = 1:n */
  k = (int)n;
  i = motclus->size[0];
  motclus->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  emxInit_int32_T(&b_r, 2);
  emxInit_boolean_T(&b_bin, 2);
  for (b_i = 0; b_i < k; b_i++) {
    /* 'gkmPWMlasso4:466' motclus{i} = find(bin == i); */
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
    const emxArray_real_T *loc, double varargin_10, emxArray_real_T *weigmat2,
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
  b_lasso_cvmat(b_B, b_cfile2, varargin_10, weigmat2, expl_temp_data,
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
  /* 'gkmPWMlasso4:462' bin = conncomp(graph(simmat>r)); */
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
  /* 'gkmPWMlasso4:463' n = max(bin); */
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
  /* 'gkmPWMlasso4:464' motclus = cell(n,1); */
  /* 'gkmPWMlasso4:465' for i = 1:n */
  k = (int)n;
  i = motclus->size[0];
  motclus->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  emxInit_int32_T(&c_r, 2);
  emxInit_boolean_T(&b_bin, 2);
  for (b_i = 0; b_i < k; b_i++) {
    /* 'gkmPWMlasso4:466' motclus{i} = find(bin == i); */
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
  /* 'gkmPWMlasso4:424' [a, b] = sort(pvec, 'descend'); */
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
  /* 'gkmPWMlasso4:425' c = [weigvec(b) a E(b)]; */
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
  /* 'gkmPWMlasso4:427' mylen = length(b); */
  /* 'gkmPWMlasso4:428' ordered_motclus = cell(mylen, 1); */
  /* 'gkmPWMlasso4:429' for idx=1:mylen */
  i = iidx->size[0];
  i1 = ordered_motclus->size[0];
  ordered_motclus->size[0] = iidx->size[0];
  emxEnsureCapacity_cell_wrap_1(ordered_motclus, i1);
  ordered_motclus_data = ordered_motclus->data;
  for (k = 0; k < i; k++) {
    /* 'gkmPWMlasso4:430' ordered_motclus{idx} = motclus{b(idx)}; */
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
  /* 'gkmPWMlasso4:432' motclus = ordered_motclus; */
  /* 'gkmPWMlasso4:433' [rr,cc] = size(c); */
  /* 'gkmPWMlasso4:434' names = cell(n,1); */
  nbytes = (int)n;
  i = names->size[0];
  names->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_8(names, i);
  names_data = names->data;
  for (i = 0; i < nbytes; i++) {
    names_data[i].f1->size[0] = 1;
    names_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso4:435' names = coder.nullcopy(names); */
  /* 'gkmPWMlasso4:436' i = 0; */
  b_i = 0.0;
  /* 'gkmPWMlasso4:437' fid = fopen(memefile, 'r'); */
  fileid = cfopen(memefile, "rb");
  /* 'gkmPWMlasso4:438' while ~feof(fid) */
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
      /* 'gkmPWMlasso4:439' line = fgetl(fid); */
      fgetl(fileid, line);
      line_data = line->data;
      /* 'gkmPWMlasso4:440' if length(line) >= 5 */
      if (line->size[1] >= 5) {
        /* 'gkmPWMlasso4:441' if strcmp(line(1:5), 'MOTIF') */
        for (i = 0; i < 5; i++) {
          c_a[i] = line_data[i];
        }
        nbytes = memcmp(&c_a[0], &b[0], 5);
        if (nbytes == 0) {
          /* 'gkmPWMlasso4:442' i = i+1; */
          b_i++;
          /* 'gkmPWMlasso4:443' [~,name] = strtok(line); */
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
          /* 'gkmPWMlasso4:444' names{i} = strtrim(name); */
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
  /* 'gkmPWMlasso4:448' fclose(fid); */
  cfclose(fileid);
  /* 'gkmPWMlasso4:449' fidw = fopen(sprintf('%s_gkmPWMlasso4.out', filename),
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
  nbytes = snprintf(NULL, 0, "%s_gkmPWMlasso4.out", &charStr_data[0]);
  i = charStr->size[0] * charStr->size[1];
  charStr->size[0] = 1;
  charStr->size[1] = nbytes + 1;
  emxEnsureCapacity_char_T(charStr, i);
  charStr_data = charStr->data;
  snprintf(&charStr_data[0], (size_t)(nbytes + 1), "%s_gkmPWMlasso4.out",
           &line_data[0]);
  i = charStr->size[0] * charStr->size[1];
  if (1 > nbytes) {
    charStr->size[1] = 0;
  } else {
    charStr->size[1] = nbytes;
  }
  emxEnsureCapacity_char_T(charStr, i);
  fileid = cfopen(charStr, "wb");
  /* 'gkmPWMlasso4:450' fprintf(fidw, 'Minimum PWM Length:\t%d\n', int32(minL));
   */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  emxFree_char_T(&varargin_1);
  emxFree_char_T(&charStr);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Minimum PWM Length:\t%d\n", (int)rt_roundd(minL));
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso4:451' fprintf(fidw, 'Minimum PWM Information:\t%f\n',
   * minInfo); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Minimum PWM Information:\t%f\n", minInfo);
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso4:452' fprintf(fidw, 'Correlation with SVM weights:\t%f\n', C);
   */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Correlation with SVM weights:\t%f\n", C);
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso4:453' fprintf(fidw, 'Cluster ID\tMotif ID\tMotif
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
  /* 'gkmPWMlasso4:454' for j = 1:length(E) */
  i = E->size[0];
  for (k = 0; k < i; k++) {
    /* 'gkmPWMlasso4:455' for l = 1:length(motclus{j}) */
    i1 = ordered_motclus_data[k].f1->size[1];
    if (0 <= ordered_motclus_data[k].f1->size[1] - 1) {
      c_NULL = NULL;
    }
    for (l = 0; l < i1; l++) {
      /* 'gkmPWMlasso4:456' fprintf(fidw, '%d\t%d\t%s\t%0.3f\t%0.3f\t%0.3f\n',
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
        fprintf(filestar, "%d\t%d\t%s\t%0.3f\t%0.3f\t%0.3f\n", k + 1,
                (int)rt_roundd(ordered_motclus_data[k].f1->data[l]),
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
  /* 'gkmPWMlasso4:459' fclose(fidw); */
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
  /* 'gkmPWMlasso4:470' l = length(p); */
  /* 'gkmPWMlasso4:471' info = zeros(l, 1); */
  nx = p->size[0];
  i = info->size[0];
  info->size[0] = nx;
  emxEnsureCapacity_real_T(info, i);
  info_data = info->data;
  /* 'gkmPWMlasso4:472' len = zeros(l,1); */
  i = len->size[0];
  len->size[0] = nx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;
  /* 'gkmPWMlasso4:473' for i = 1:l */
  i = p->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&vec, 1);
  emxInit_real_T(&mvec, 1);
  emxInit_real_T(&b_r, 2);
  emxInit_real_T(&b_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWMlasso4:474' mat = p{i}+(p{i}==0); */
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
    /* 'gkmPWMlasso4:475' vec = 2+sum(mat.*log(mat)/log(2),2); */
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
    /* 'gkmPWMlasso4:476' mvec = movmean(vec,5); */
    nx = vec->size[0];
    if (2 <= nx) {
      nx = 2;
    }
    vmovfun(vec, vec->size[0], nx, nx, mvec);
    mvec_data = mvec->data;
    /* 'gkmPWMlasso4:477' while min(vec(1),mvec(1)) < cut && length(vec) > 1 */
    while ((fmin(mat_data[0], mvec_data[0]) < 0.0) && (vec->size[0] > 1)) {
      /* 'gkmPWMlasso4:478' p{i}(1,:) = []; */
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
      /* 'gkmPWMlasso4:479' vec(1) = []; */
      nx = vec->size[0];
      nxout = vec->size[0];
      for (nrows = 0; nrows <= nxout - 2; nrows++) {
        mat_data[nrows] = mat_data[nrows + 1];
      }
      idx = vec->size[0];
      vec->size[0] = nx - 1;
      emxEnsureCapacity_real_T(vec, idx);
      mat_data = vec->data;
      /* 'gkmPWMlasso4:480' mvec(1)=[]; */
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
      /* 'gkmPWMlasso4:481' mat(1,:) = []; */
      idx = mat->size[0] * mat->size[1];
      if (1 > mat->size[0] - 1) {
        mat->size[0] = 0;
      } else {
        mat->size[0]--;
      }
      emxEnsureCapacity_real_T(mat, idx);
    }
    /* 'gkmPWMlasso4:483' while min(mvec(end),vec(end)) < cut && length(vec) > 1
     */
    while ((fmin(mvec_data[mvec->size[0] - 1], mat_data[vec->size[0] - 1]) <
            0.0) &&
           (vec->size[0] > 1)) {
      /* 'gkmPWMlasso4:484' vec(end) = []; */
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
      /* 'gkmPWMlasso4:485' mvec(end)=[]; */
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
      /* 'gkmPWMlasso4:486' mat(end,:) = []; */
      idx = mat->size[0] * mat->size[1];
      if (1 > mat->size[0] - 1) {
        mat->size[0] = 0;
      } else {
        mat->size[0]--;
      }
      emxEnsureCapacity_real_T(mat, idx);
      /* 'gkmPWMlasso4:487' p{i}(end,:) = []; */
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
    /* 'gkmPWMlasso4:489' info(i) = sum(vec); */
    info_data[b_i] = blockedSummation(vec, vec->size[0]);
    /* 'gkmPWMlasso4:490' [len(i), ~] = size(mat); */
    len_data[b_i] = mat->size[0];
  }
  emxFree_real_T(&b_mat);
  emxFree_real_T(&b_r);
  emxFree_real_T(&mvec);
  emxFree_real_T(&vec);
  emxFree_real_T(&mat);
  /* 'gkmPWMlasso4:492' pp = p; */
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
void gkmPWMlasso4(const emxArray_char_T *varargin_1,
                  const emxArray_char_T *varargin_2, double varargin_3,
                  double varargin_4, double varargin_5, double varargin_6,
                  double varargin_7, double varargin_8, bool varargin_9,
                  double varargin_10)
{
  static const char cv[5] = {'M', 'O', 'T', 'I', 'F'};
  cell_wrap_0 *p_data;
  cell_wrap_1 *b_motclus_data;
  cell_wrap_1 *motclus_data;
  cell_wrap_1 *tmp_motclus_data;
  emxArray_boolean_T b_MSE2_data;
  emxArray_boolean_T *b_loc;
  emxArray_cell_wrap_0 *p;
  emxArray_cell_wrap_1 *c_motclus;
  emxArray_cell_wrap_1 *motclus;
  emxArray_cell_wrap_1 *tmp_motclus;
  emxArray_cell_wrap_3_1x19 FF;
  emxArray_cell_wrap_3_1x19 b_FF;
  emxArray_cell_wrap_3_1x19 b_newFF;
  emxArray_cell_wrap_3_1x19 newFF;
  emxArray_cell_wrap_3_20 F;
  emxArray_char_T *text;
  emxArray_int32_T *ib;
  emxArray_int32_T *idx;
  emxArray_int32_T *match_out;
  emxArray_int32_T *matches;
  emxArray_real_T b_B;
  emxArray_real_T *A;
  emxArray_real_T *AA;
  emxArray_real_T *B;
  emxArray_real_T *BB;
  emxArray_real_T *BY;
  emxArray_real_T *GCmat;
  emxArray_real_T *Z;
  emxArray_real_T *b_GCmat;
  emxArray_real_T *b_xc;
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
  emxArray_real_T *weigmat2;
  emxArray_real_T *xc;
  emxArray_real_T *y;
  double comb2_data[2508];
  double comb_data[2508];
  double diffc_data[418];
  double indc_data[418];
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
  double GCneg1;
  double GCpos1;
  double a;
  double b1;
  double k_svm;
  double l_svm;
  double rcnum;
  double *AA_data;
  double *BY_data;
  double *B_data;
  double *Z_data;
  double *b_comb_data;
  double *b_diffc_data;
  double *b_indc_data;
  double *cfile2_data;
  double *f_data;
  double *indvec_data;
  double *loc_data;
  double *negvec_data;
  double *normvec_data;
  double *xc_data;
  int b_f_data[19];
  int iidx_data[19];
  int MSE2_size[2];
  int a__3_DF_size[2];
  int b_MSE2_size[2];
  int comb2_size[2];
  int comb_size[2];
  int expl_temp_size[2];
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
  int i7;
  int i_loop_ub;
  int indc_size;
  int j;
  int j_loop_ub;
  int k;
  int k_loop_ub;
  int l_loop_ub;
  int lk;
  int loop_ub;
  int m_loop_ub;
  int n_loop_ub;
  int nx;
  int nxout;
  int size_tmp_idx_1;
  int u1;
  int *match_out_data;
  int *matches_data;
  char *text_data;
  bool c_MSE2_data[19];
  bool empty_non_axis_sizes;
  bool exitg1;
  bool *b_loc_data;
  if (!isInitialized_gkmPWMlasso4) {
    gkmPWMlasso4_initialize();
  }
  emxInit_real_T(&comb, 2);
  emxInit_real_T(&comb2, 2);
  emxInit_real_T(&diffc, 1);
  emxInit_real_T(&indc, 1);
  emxInit_real_T(&xc, 2);
  /* 'gkmPWMlasso4:3' if nargin < 2 */
  /* 'gkmPWMlasso4:7' filename = varargin{1}; */
  /* 'gkmPWMlasso4:8' memefile = varargin{2}; */
  /* 'gkmPWMlasso4:9' minL = varargin{3}; */
  /* 'gkmPWMlasso4:10' minInfo = varargin{4}; */
  /* 'gkmPWMlasso4:11' corrCut = varargin{5}; */
  /* 'gkmPWMlasso4:12' l_svm = varargin{6}; */
  l_svm = varargin_6;
  /* 'gkmPWMlasso4:13' k_svm = varargin{7}; */
  k_svm = varargin_7;
  /* 'gkmPWMlasso4:14' BG_GC = varargin{8}; */
  /* 'gkmPWMlasso4:15' RC = varargin{9}; */
  /* 'gkmPWMlasso4:16' d = varargin{10}; */
  /*  minL = 8; */
  /*  minInfo = 0.5; */
  /*  corrCut = 0.85; */
  /*  l_svm = 10; */
  /*  k_svm = 6; */
  /*  BG_GC = 0; */
  /*  RC = true; */
  /*  if nargin == 3 */
  /*      d = varargin{3}; */
  /*  else */
  /*      d = 0; */
  /*  end */
  /*  if nargin > 2 */
  /*      f = find(strcmp('MinLength', varargin)); */
  /*      if ~isempty(f); */
  /*          minL = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('MinInfo', varargin)); */
  /*      if ~isempty(f); */
  /*          minInfo = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('CorrCutoff', varargin)); */
  /*      if ~isempty(f); */
  /*          corrCut = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('l', varargin)); */
  /*      if ~isempty(f); */
  /*          l_svm = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('k', varargin)); */
  /*      if ~isempty(f); */
  /*          k_svm = varargin{f+1}; */
  /*      end */
  /*      f = find(strcmp('Mode', varargin)); */
  /*      if ~isempty(f) && strcmp('Compare',varargin{f+1}); */
  /*          BG_GC = 1; */
  /*      end */
  /*      f = find(strcmp('RC', varargin)); */
  /*      if ~isempty(f); */
  /*          RC = varargin{f+1}; */
  /*      end */
  /*  end */
  /* 'gkmPWMlasso4:62' [comb,comb2,diffc,indc,xc,rcnum] = genIndex(l_svm,k_svm);
   */
  genIndex(varargin_6, varargin_7, comb, comb2, diffc, indc, xc, &rcnum);
  /* generate gapped positions, adjusted for reverse complements */
  /* 'gkmPWMlasso4:64' if length(comb)*4^k_svm > 10^6 */
  if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
    u1 = 0;
  } else {
    nx = comb->size[0];
    u1 = comb->size[1];
    if (nx >= u1) {
      u1 = nx;
    }
  }
  if ((double)u1 * pow(4.0, varargin_7) > 1.0E+6) {
    emxInit_real_T(&b_xc, 2);
    /* error([num2str(length(comb)*4^k_svm) ' exceeds the maximum number of
     * gapped kmers allowed (10^6)']) */
    /* 'gkmPWMlasso4:66' fprintf('l = 10, k = 6\n'); */
    printf("l = 10, k = 6\n");
    fflush(stdout);
    /* 'gkmPWMlasso4:67' lk = 0; */
    lk = 0;
    /* 'gkmPWMlasso4:68' l_svm2 = l_svm; */
    /* 'gkmPWMlasso4:69' k_svm2 = k_svm; */
    /* 'gkmPWMlasso4:70' l_svm = 10; */
    l_svm = 10.0;
    /* 'gkmPWMlasso4:71' k_svm = 6; */
    k_svm = 6.0;
    /* 'gkmPWMlasso4:72' [comb,comb2,diffc,indc,xc,rcnum] =
     * genIndex(l_svm,k_svm); */
    b_genIndex(comb_data, comb_size, comb2_data, comb2_size, diffc_data, &nx,
               indc_data, &indc_size, b_xc, &rcnum);
    AA_data = b_xc->data;
    i = comb->size[0] * comb->size[1];
    comb->size[0] = comb_size[0];
    comb->size[1] = 6;
    emxEnsureCapacity_real_T(comb, i);
    b_comb_data = comb->data;
    loop_ub = comb_size[0] * 6;
    for (i = 0; i < loop_ub; i++) {
      b_comb_data[i] = comb_data[i];
    }
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = comb2_size[0];
    comb2->size[1] = 6;
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    loop_ub = comb2_size[0] * 6;
    for (i = 0; i < loop_ub; i++) {
      b_comb_data[i] = comb2_data[i];
    }
    i = diffc->size[0];
    diffc->size[0] = nx;
    emxEnsureCapacity_real_T(diffc, i);
    b_diffc_data = diffc->data;
    for (i = 0; i < nx; i++) {
      b_diffc_data[i] = diffc_data[i];
    }
    i = indc->size[0];
    indc->size[0] = indc_size;
    emxEnsureCapacity_real_T(indc, i);
    b_indc_data = indc->data;
    for (i = 0; i < indc_size; i++) {
      b_indc_data[i] = indc_data[i];
    }
    i = xc->size[0] * xc->size[1];
    xc->size[0] = b_xc->size[0];
    xc->size[1] = b_xc->size[1];
    emxEnsureCapacity_real_T(xc, i);
    xc_data = xc->data;
    loop_ub = b_xc->size[0] * b_xc->size[1];
    for (i = 0; i < loop_ub; i++) {
      xc_data[i] = AA_data[i];
    }
    emxFree_real_T(&b_xc);
  } else {
    /* 'gkmPWMlasso4:73' else */
    /* 'gkmPWMlasso4:74' lk = 1; */
    lk = 1;
    /* 'gkmPWMlasso4:75' l_svm2 = l_svm; */
    /* 'gkmPWMlasso4:76' k_svm2 = k_svm; */
  }
  emxInit_real_T(&cfile2, 1);
  /* 'gkmPWMlasso4:79' fprintf('Counting gapped kmers\n'); */
  printf("Counting gapped kmers\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:81' [cfile, GCpos1, GCneg1,mat,mat2] = getgkmcounts(filename,
   * l_svm2, k_svm2, lk, RC); */
  getgkmcounts(varargin_1, varargin_6, varargin_7, lk, varargin_9, cfile2,
               &GCpos1, &GCneg1, mat, mat2);
  cfile2_data = cfile2->data;
  /* 'gkmPWMlasso4:82' if BG_GC == 1 */
  if (varargin_8 == 1.0) {
    /* 'gkmPWMlasso4:83' mat = (mat+mat2)/2; */
    for (i = 0; i < 16; i++) {
      mat[i] = (mat[i] + mat2[i]) / 2.0;
    }
    /* 'gkmPWMlasso4:84' GCpos1 = (GCpos1+GCneg1)/2; */
    GCpos1 = (GCpos1 + GCneg1) / 2.0;
    /* 'gkmPWMlasso4:85' GCneg1 = GCpos1; */
    GCneg1 = GCpos1;
  }
  emxInit_real_T(&negvec, 1);
  emxInit_char_T(&text, 2);
  /* 'gkmPWMlasso4:87' fprintf('Generating negative set\n'); */
  printf("Generating negative set\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:88' negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC);
   */
  BGkmer(mat, GCneg1, comb, rcnum, l_svm, k_svm, varargin_9, negvec);
  negvec_data = negvec->data;
  /* 'gkmPWMlasso4:90' fprintf('Filtering motifs\n'); */
  printf("Filtering motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:91' num = length(strfind(fileread(memefile),'MOTIF')); */
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
    lk = 0;
    i = text->size[1];
    for (b_i = 0; b_i <= i - 5; b_i++) {
      j = 1;
      while ((j <= 5) && (text_data[(b_i + j) - 1] == cv[j - 1])) {
        j++;
      }
      if (j > 5) {
        matches_data[lk] = b_i + 1;
        lk++;
      }
    }
    i = match_out->size[0] * match_out->size[1];
    match_out->size[0] = 1;
    match_out->size[1] = lk;
    emxEnsureCapacity_int32_T(match_out, i);
    match_out_data = match_out->data;
    for (b_i = 0; b_i < lk; b_i++) {
      match_out_data[b_i] = matches_data[b_i];
    }
    emxFree_int32_T(&matches);
    size_tmp_idx_1 = match_out->size[1];
  }
  /* 'gkmPWMlasso4:92' p = getmotif(memefile,1:num); */
  emxInit_real_T(&y, 2);
  if (size_tmp_idx_1 < 1) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = size_tmp_idx_1;
    emxEnsureCapacity_real_T(y, i);
    AA_data = y->data;
    loop_ub = size_tmp_idx_1 - 1;
    for (i = 0; i <= loop_ub; i++) {
      AA_data[i] = (double)i + 1.0;
    }
  }
  emxInit_cell_wrap_0(&p);
  getmotif(varargin_2, y, p);
  p_data = p->data;
  /* 'gkmPWMlasso4:93' for i = 1:num */
  for (b_i = 0; b_i < size_tmp_idx_1; b_i++) {
    /* 'gkmPWMlasso4:94' [r c] = size(p{i}); */
    /* 'gkmPWMlasso4:95' for j = 1:r */
    i = p_data[b_i].f1->size[0];
    for (j = 0; j < i; j++) {
      /* 'gkmPWMlasso4:96' a = sum(p{i}(j,:)); */
      loop_ub = p_data[b_i].f1->size[1];
      i1 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = p_data[b_i].f1->size[1];
      emxEnsureCapacity_real_T(y, i1);
      AA_data = y->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        AA_data[i1] = p_data[b_i].f1->data[j + p_data[b_i].f1->size[0] * i1];
      }
      a = sum(y);
      /* 'gkmPWMlasso4:97' if abs(a-1)>0 */
      if (fabs(a - 1.0) > 0.0) {
        /* 'gkmPWMlasso4:98' [b1 loc] = max(p{i}(j,:)); */
        loop_ub = p_data[b_i].f1->size[1];
        i1 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = p_data[b_i].f1->size[1];
        emxEnsureCapacity_real_T(y, i1);
        AA_data = y->data;
        for (i1 = 0; i1 < loop_ub; i1++) {
          AA_data[i1] = p_data[b_i].f1->data[j + p_data[b_i].f1->size[0] * i1];
        }
        b_maximum(y, &b1, &nx);
        /* 'gkmPWMlasso4:99' p{i}(j,loc) = b1-a+1; */
        p_data[b_i].f1->data[j + p_data[b_i].f1->size[0] * (nx - 1)] =
            (b1 - a) + 1.0;
      }
    }
  }
  emxInit_real_T(&Z, 1);
  emxInit_real_T(&loc, 1);
  /* 'gkmPWMlasso4:103' [p, info, lenvec] = trim_pwm(p,0.0); */
  trim_pwm(p, loc, Z);
  Z_data = Z->data;
  loc_data = loc->data;
  p_data = p->data;
  /* 'gkmPWMlasso4:104' indvec =
   * intersect(find(info./lenvec>=minInfo),find(lenvec>=minL)); */
  emxInit_int32_T(&idx, 1);
  emxInit_boolean_T(&b_loc, 1);
  if (loc->size[0] == Z->size[0]) {
    i = b_loc->size[0];
    b_loc->size[0] = loc->size[0];
    emxEnsureCapacity_boolean_T(b_loc, i);
    b_loc_data = b_loc->data;
    loop_ub = loc->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_loc_data[i] = (loc_data[i] / Z_data[i] >= varargin_4);
    }
    b_eml_find(b_loc, idx);
    matches_data = idx->data;
  } else {
    e_binary_expand_op(idx, loc, Z, varargin_4);
    matches_data = idx->data;
  }
  i = b_loc->size[0];
  b_loc->size[0] = Z->size[0];
  emxEnsureCapacity_boolean_T(b_loc, i);
  b_loc_data = b_loc->data;
  loop_ub = Z->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_loc_data[i] = (Z_data[i] >= varargin_3);
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
  emxInit_real_T(&f, 1);
  i = f->size[0];
  f->size[0] = ib->size[0];
  emxEnsureCapacity_real_T(f, i);
  f_data = f->data;
  loop_ub = ib->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = match_out_data[i];
  }
  emxInit_real_T(&indvec, 1);
  do_vectors(BY, f, indvec, idx, ib);
  indvec_data = indvec->data;
  /* 'gkmPWMlasso4:105' n = length(indvec); */
  /* 'gkmPWMlasso4:106' fprintf('Mapping PWMs to gkm space\n'); */
  printf("Mapping PWMs to gkm space\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:107' lcnum = length(comb); */
  emxFree_int32_T(&ib);
  if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
    u1 = 0;
  } else {
    nx = comb->size[0];
    u1 = comb->size[1];
    if (nx >= u1) {
      u1 = nx;
    }
  }
  emxInit_real_T(&A, 2);
  /* 'gkmPWMlasso4:108' A=zeros(lcnum*4^k_svm,n); */
  i = (int)((double)u1 * pow(4.0, k_svm));
  i1 = A->size[0] * A->size[1];
  A->size[0] = i;
  A->size[1] = indvec->size[0];
  emxEnsureCapacity_real_T(A, i1);
  xc_data = A->data;
  loop_ub = i * indvec->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    xc_data[i1] = 0.0;
  }
  emxInit_real_T(&AA, 2);
  /* 'gkmPWMlasso4:109' AA=zeros(lcnum*4^k_svm,n); */
  i1 = AA->size[0] * AA->size[1];
  AA->size[0] = i;
  AA->size[1] = indvec->size[0];
  emxEnsureCapacity_real_T(AA, i1);
  AA_data = AA->data;
  loop_ub = i * indvec->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    AA_data[i1] = 0.0;
  }
  emxInit_real_T(&GCmat, 2);
  emxInit_real_T(&normvec, 1);
  /* 'gkmPWMlasso4:110' GCmat = repmat([0.5-GCpos1/2 GCpos1/2 GCpos1/2
   * 0.5-GCpos1/2],l_svm-1,1); */
  b1 = 0.5 - GCpos1 / 2.0;
  dv[0] = b1;
  dv[1] = GCpos1 / 2.0;
  dv[2] = GCpos1 / 2.0;
  dv[3] = b1;
  b_repmat(dv, l_svm - 1.0, GCmat);
  b_comb_data = GCmat->data;
  /* 'gkmPWMlasso4:111' per = 10; */
  a = 10.0;
  /* 'gkmPWMlasso4:112' normvec = zeros(n,1); */
  /* 'gkmPWMlasso4:113' for j = 1:n */
  i1 = indvec->size[0];
  emxInit_real_T(&b_GCmat, 2);
  for (j = 0; j < i1; j++) {
    /* 'gkmPWMlasso4:114' if mod(j, floor(n/10))==0 */
    indc_size = (int)floor((double)indvec->size[0] / 10.0);
    nx = j + 1;
    if (indc_size != 0) {
      nx = (int)fmod((double)j + 1.0, indc_size);
    }
    if (nx == 0) {
      /* 'gkmPWMlasso4:115' fprintf('%d...', int32(per)); */
      printf("%d...", (int)a);
      fflush(stdout);
      /* 'gkmPWMlasso4:116' per = per+10; */
      a += 10.0;
    }
    /* 'gkmPWMlasso4:118' loc = zeros(l_svm*2-2+lenvec(indvec(j)), 1); */
    i2 = (int)indvec_data[j] - 1;
    b1 = Z_data[i2];
    loop_ub = (int)((l_svm * 2.0 - 2.0) + b1);
    b_i = loc->size[0];
    loc->size[0] = loop_ub;
    emxEnsureCapacity_real_T(loc, b_i);
    loc_data = loc->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      loc_data[b_i] = 0.0;
    }
    /* 'gkmPWMlasso4:119' loc(l_svm:lenvec(indvec(j))+l_svm-1) = 1; */
    GCneg1 = (b1 + l_svm) - 1.0;
    if (l_svm > GCneg1) {
      b_i = -1;
      i3 = 0;
    } else {
      b_i = (int)l_svm - 2;
      i3 = (int)GCneg1;
    }
    loop_ub = (i3 - b_i) - 1;
    for (i3 = 0; i3 < loop_ub; i3++) {
      loc_data[(b_i + i3) + 1] = 1.0;
    }
    /* 'gkmPWMlasso4:120' if RC */
    if (varargin_9) {
      /* 'gkmPWMlasso4:121' A(:,j) =
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
        nx = GCmat->size[0];
      } else {
        nx = 0;
      }
      if (empty_non_axis_sizes ||
          ((p_data[(int)indvec_data[j] - 1].f1->size[0] != 0) &&
           (p_data[(int)indvec_data[j] - 1].f1->size[1] != 0))) {
        lk = p_data[(int)indvec_data[j] - 1].f1->size[0];
      } else {
        lk = 0;
      }
      if (empty_non_axis_sizes || (GCmat->size[0] != 0)) {
        indc_size = GCmat->size[0];
      } else {
        indc_size = 0;
      }
      b_i = b_GCmat->size[0] * b_GCmat->size[1];
      b_GCmat->size[0] = (nx + lk) + indc_size;
      b_GCmat->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_GCmat, b_i);
      normvec_data = b_GCmat->data;
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < nx; i3++) {
          normvec_data[i3 + b_GCmat->size[0] * b_i] =
              b_comb_data[i3 + nx * b_i];
        }
      }
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < lk; i3++) {
          normvec_data[(i3 + nx) + b_GCmat->size[0] * b_i] =
              p_data[i2].f1->data[i3 + lk * b_i];
        }
      }
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (b_i = 0; b_i < indc_size; b_i++) {
          normvec_data[((b_i + nx) + lk) + b_GCmat->size[0] * i2] =
              b_comb_data[b_i + indc_size * i2];
        }
      }
      PWM2kmers(b_GCmat, mat, comb2, diffc, indc, loc, xc, l_svm, k_svm, rcnum,
                f);
      f_data = f->data;
      loop_ub = f->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        xc_data[i2 + A->size[0] * j] = f_data[i2];
      }
    } else {
      /* 'gkmPWMlasso4:122' else */
      /* 'gkmPWMlasso4:123' A(:,j) =
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
        nx = GCmat->size[0];
      } else {
        nx = 0;
      }
      if (empty_non_axis_sizes ||
          ((p_data[(int)indvec_data[j] - 1].f1->size[0] != 0) &&
           (p_data[(int)indvec_data[j] - 1].f1->size[1] != 0))) {
        lk = p_data[(int)indvec_data[j] - 1].f1->size[0];
      } else {
        lk = 0;
      }
      if (empty_non_axis_sizes || (GCmat->size[0] != 0)) {
        indc_size = GCmat->size[0];
      } else {
        indc_size = 0;
      }
      b_i = b_GCmat->size[0] * b_GCmat->size[1];
      b_GCmat->size[0] = (nx + lk) + indc_size;
      b_GCmat->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_GCmat, b_i);
      normvec_data = b_GCmat->data;
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < nx; i3++) {
          normvec_data[i3 + b_GCmat->size[0] * b_i] =
              b_comb_data[i3 + nx * b_i];
        }
      }
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < lk; i3++) {
          normvec_data[(i3 + nx) + b_GCmat->size[0] * b_i] =
              p_data[i2].f1->data[i3 + lk * b_i];
        }
      }
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (b_i = 0; b_i < indc_size; b_i++) {
          normvec_data[((b_i + nx) + lk) + b_GCmat->size[0] * i2] =
              b_comb_data[b_i + indc_size * i2];
        }
      }
      PWM2kmers_norc(b_GCmat, mat, comb2, diffc, indc, loc, xc, l_svm, k_svm,
                     f);
      f_data = f->data;
      loop_ub = f->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        xc_data[i2 + A->size[0] * j] = f_data[i2];
      }
    }
    /* 'gkmPWMlasso4:125' A(:,j) = A(:,j) - negvec*(l_svm-1+lenvec(indvec(j)));
     */
    GCpos1 = (l_svm - 1.0) + b1;
    if (A->size[0] == negvec->size[0]) {
      nx = A->size[0] - 1;
      i2 = f->size[0];
      f->size[0] = A->size[0];
      emxEnsureCapacity_real_T(f, i2);
      f_data = f->data;
      for (i2 = 0; i2 <= nx; i2++) {
        f_data[i2] = xc_data[i2 + A->size[0] * j] - negvec_data[i2] * GCpos1;
      }
      loop_ub = f->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        xc_data[i2 + A->size[0] * j] = f_data[i2];
      }
    } else {
      binary_expand_op(A, j, negvec, GCpos1);
      xc_data = A->data;
    }
    /* 'gkmPWMlasso4:126' normvec(j) = (A(:,j)'*A(:,j))^0.5; */
    loop_ub = A->size[0];
    i2 = loc->size[0];
    loc->size[0] = A->size[0];
    emxEnsureCapacity_real_T(loc, i2);
    loc_data = loc->data;
    i2 = BY->size[0];
    BY->size[0] = A->size[0];
    emxEnsureCapacity_real_T(BY, i2);
    BY_data = BY->data;
    for (i2 = 0; i2 < loop_ub; i2++) {
      b1 = xc_data[i2 + A->size[0] * j];
      loc_data[i2] = b1;
      BY_data[i2] = b1;
    }
    loop_ub = A->size[0];
    if (A->size[0] < 1) {
      b1 = 0.0;
    } else {
      b1 = cblas_ddot((blasint)A->size[0], &loc_data[0], (blasint)1,
                      &BY_data[0], (blasint)1);
    }
    GCpos1 = sqrt(b1);
    /* 'gkmPWMlasso4:127' AA(:,j) = A(:,j)/normvec(j); */
    for (i2 = 0; i2 < loop_ub; i2++) {
      AA_data[i2 + AA->size[0] * j] = xc_data[i2 + A->size[0] * j] / GCpos1;
    }
  }
  emxFree_real_T(&GCmat);
  emxFree_cell_wrap_0(&p);
  emxInit_cell_wrap_1(&motclus);
  /* 'gkmPWMlasso4:129' fprintf('\n'); */
  printf("\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:130' fprintf('Clustering motifs\n'); */
  printf("Clustering motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:131' simmat = AA'*AA; */
  b_mtimes(AA, AA, comb);
  /*  clear AA */
  /* 'gkmPWMlasso4:133' motclus = clus_simmat_eig(simmat,corrCut); */
  clus_simmat_eig(comb, varargin_5, motclus);
  motclus_data = motclus->data;
  /* 'gkmPWMlasso4:134' fprintf('Number of motif clusters: %d\n',
   * int32(length(motclus))); */
  printf("Number of motif clusters: %d\n", motclus->size[0]);
  fflush(stdout);
  /* 'gkmPWMlasso4:136' fprintf('Selecting Motifs\n'); */
  printf("Selecting Motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:137' cfile2 = cfile-negvec/sum(negvec)*sum(cfile); */
  GCpos1 = blockedSummation(cfile2, cfile2->size[0]);
  GCneg1 = blockedSummation(negvec, negvec->size[0]);
  emxFree_real_T(&AA);
  if (cfile2->size[0] == negvec->size[0]) {
    loop_ub = cfile2->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      cfile2_data[i1] -= negvec_data[i1] / GCneg1 * GCpos1;
    }
  } else {
    d_binary_expand_op(cfile2, negvec, GCneg1, GCpos1);
    cfile2_data = cfile2->data;
  }
  /* 'gkmPWMlasso4:138' cfile2 = cfile2/std(cfile2); */
  b1 = b_std(cfile2);
  loop_ub = cfile2->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    cfile2_data[i1] /= b1;
  }
  emxInit_real_T(&B, 2);
  /* 'gkmPWMlasso4:139' B = zeros((4^k_svm)*lcnum, length(motclus)); */
  i1 = B->size[0] * B->size[1];
  B->size[0] = i;
  B->size[1] = motclus->size[0];
  emxEnsureCapacity_real_T(B, i1);
  B_data = B->data;
  loop_ub = i * motclus->size[0];
  for (i = 0; i < loop_ub; i++) {
    B_data[i] = 0.0;
  }
  /* 'gkmPWMlasso4:140' corrvec = zeros(n,1); */
  /* 'gkmPWMlasso4:141' zvec = zeros(n,1); */
  /* 'gkmPWMlasso4:142' Z = zeros(length(motclus),1); */
  i = Z->size[0];
  Z->size[0] = motclus->size[0];
  emxEnsureCapacity_real_T(Z, i);
  Z_data = Z->data;
  loop_ub = motclus->size[0];
  for (i = 0; i < loop_ub; i++) {
    Z_data[i] = 0.0;
  }
  /* 'gkmPWMlasso4:143' for i = 1:n */
  i = indvec->size[0];
  i1 = normvec->size[0];
  normvec->size[0] = indvec->size[0];
  emxEnsureCapacity_real_T(normvec, i1);
  normvec_data = normvec->data;
  i1 = diffc->size[0];
  diffc->size[0] = indvec->size[0];
  emxEnsureCapacity_real_T(diffc, i1);
  b_diffc_data = diffc->data;
  if (0 <= indvec->size[0] - 1) {
    GCneg1 = (double)u1 * 10.0;
    if (GCneg1 <= A->size[0]) {
      k = (int)GCneg1;
    } else {
      k = A->size[0];
    }
    i4 = A->size[0];
    b_loop_ub = A->size[0];
    if (1 > u1) {
      c_loop_ub = 0;
    } else {
      c_loop_ub = u1;
    }
  }
  for (b_i = 0; b_i < i; b_i++) {
    /*  [~,I] = sort(A(:,i),'descend'); */
    /*  zvec(i) = mean(cfile2(I(1:lcnum))); */
    /* 'gkmPWMlasso4:146' [~,I2] = maxk(A(:,i), lcnum*10); */
    i1 = f->size[0];
    f->size[0] = i4;
    emxEnsureCapacity_real_T(f, i1);
    f_data = f->data;
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      f_data[i1] = xc_data[i1 + A->size[0] * b_i];
    }
    exkib(f, k, idx, BY);
    matches_data = idx->data;
    loop_ub = idx->size[0];
    i1 = loc->size[0];
    loc->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(loc, i1);
    loc_data = loc->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      loc_data[i1] = matches_data[i1];
    }
    /* 'gkmPWMlasso4:147' zvec(i) = mean(cfile2(I2(1:lcnum))); */
    i1 = BY->size[0];
    BY->size[0] = c_loop_ub;
    emxEnsureCapacity_real_T(BY, i1);
    BY_data = BY->data;
    for (i1 = 0; i1 < c_loop_ub; i1++) {
      BY_data[i1] = cfile2_data[(int)loc_data[i1] - 1];
    }
    b1 = blockedSummation(BY, c_loop_ub) / (double)c_loop_ub;
    normvec_data[b_i] = b1;
    /* 'gkmPWMlasso4:148' if zvec(i) < 0 */
    if (b1 < 0.0) {
      /*  Alternative to corr */
      /*  corrvec(i) = -1*corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10))); */
      /* 'gkmPWMlasso4:151' correlation_matrix = -1*corrcoef(A(I2,i),
       * cfile2(I2)); */
      /* 'gkmPWMlasso4:152' corrvec(i) = correlation_matrix(1,2); */
      i1 = f->size[0];
      f->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = xc_data[((int)loc_data[i1] + A->size[0] * b_i) - 1];
      }
      i1 = BY->size[0];
      BY->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(BY, i1);
      BY_data = BY->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        BY_data[i1] = cfile2_data[(int)loc_data[i1] - 1];
      }
      corrcoef(f, BY, dv);
      b_diffc_data[b_i] = -dv[2];
    } else {
      /* 'gkmPWMlasso4:153' else */
      /*  Alternative to corr */
      /*  corrvec(i) = corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10))); */
      /* 'gkmPWMlasso4:156' correlation_matrix = corrcoef(A(I2,i), cfile2(I2));
       */
      /* 'gkmPWMlasso4:157' corrvec(i) = correlation_matrix(1,2); */
      i1 = f->size[0];
      f->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = xc_data[((int)loc_data[i1] + A->size[0] * b_i) - 1];
      }
      i1 = BY->size[0];
      BY->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(BY, i1);
      BY_data = BY->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        BY_data[i1] = cfile2_data[(int)loc_data[i1] - 1];
      }
      corrcoef(f, BY, dv);
      b_diffc_data[b_i] = dv[2];
    }
  }
  /* 'gkmPWMlasso4:160' for i = 1:length(motclus) */
  i = motclus->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWMlasso4:161' [a,b] = sort(zvec(motclus{i}),'descend'); */
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
    /* 'gkmPWMlasso4:162' f = find(a == a(1)); */
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
    i1 = f->size[0];
    f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(f, i1);
    f_data = f->data;
    loop_ub = idx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      f_data[i1] = matches_data[i1];
    }
    /* 'gkmPWMlasso4:163' if length(f) > 1 */
    if (f->size[0] > 1) {
      /* 'gkmPWMlasso4:164' [~,bb] = sort(abs(corrvec(motclus{i}(b(f)))),
       * 'descend'); */
      i1 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = f->size[0];
      emxEnsureCapacity_real_T(y, i1);
      AA_data = y->data;
      loop_ub = f->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        AA_data[i1] =
            motclus_data[b_i].f1->data[(int)BY_data[(int)f_data[i1] - 1] - 1];
      }
      nx = f->size[0];
      i1 = negvec->size[0];
      negvec->size[0] = f->size[0];
      emxEnsureCapacity_real_T(negvec, i1);
      negvec_data = negvec->data;
      for (k = 0; k < nx; k++) {
        negvec_data[k] = fabs(b_diffc_data[(int)AA_data[k] - 1]);
      }
      sort(negvec, idx);
      matches_data = idx->data;
      /* 'gkmPWMlasso4:165' b(1:length(f)) = b(bb); */
      i1 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = f->size[0];
      emxEnsureCapacity_real_T(y, i1);
      AA_data = y->data;
      loop_ub = f->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        AA_data[i1] = BY_data[matches_data[i1] - 1];
      }
      loop_ub = y->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        BY_data[i1] = AA_data[i1];
      }
    }
    /* 'gkmPWMlasso4:167' B(:,i) = A(:,motclus{i}(b(1))); */
    nx = (int)motclus_data[b_i].f1->data[(int)BY_data[0] - 1];
    loop_ub = A->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      B_data[i1 + B->size[0] * b_i] = xc_data[i1 + A->size[0] * (nx - 1)];
    }
    /* 'gkmPWMlasso4:168' motclus{i} = motclus{i}(b); */
    i1 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = BY->size[0];
    emxEnsureCapacity_real_T(y, i1);
    AA_data = y->data;
    loop_ub = BY->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      AA_data[i1] = motclus_data[b_i].f1->data[(int)BY_data[i1] - 1];
    }
    i1 = motclus_data[b_i].f1->size[0] * motclus_data[b_i].f1->size[1];
    motclus_data[b_i].f1->size[0] = 1;
    motclus_data[b_i].f1->size[1] = y->size[1];
    emxEnsureCapacity_real_T(motclus_data[b_i].f1, i1);
    loop_ub = y->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      motclus_data[b_i].f1->data[i1] = AA_data[i1];
    }
    /* 'gkmPWMlasso4:169' Z(i) = zvec(motclus{i}(1)); */
    Z_data[b_i] = normvec_data[(int)motclus_data[b_i].f1->data[0] - 1];
  }
  emxFree_real_T(&y);
  emxFree_real_T(&A);
  /* 'gkmPWMlasso4:171' f = find(abs(Z)>1); */
  nx = Z->size[0];
  i = negvec->size[0];
  negvec->size[0] = Z->size[0];
  emxEnsureCapacity_real_T(negvec, i);
  negvec_data = negvec->data;
  for (k = 0; k < nx; k++) {
    negvec_data[k] = fabs(Z_data[k]);
  }
  i = b_loc->size[0];
  b_loc->size[0] = negvec->size[0];
  emxEnsureCapacity_boolean_T(b_loc, i);
  b_loc_data = b_loc->data;
  loop_ub = negvec->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_loc_data[i] = (negvec_data[i] > 1.0);
  }
  b_eml_find(b_loc, idx);
  matches_data = idx->data;
  i = f->size[0];
  f->size[0] = idx->size[0];
  emxEnsureCapacity_real_T(f, i);
  f_data = f->data;
  loop_ub = idx->size[0];
  for (i = 0; i < loop_ub; i++) {
    f_data[i] = matches_data[i];
  }
  /* 'gkmPWMlasso4:172' B = B(:,f); */
  nx = B->size[0] - 1;
  i = b_GCmat->size[0] * b_GCmat->size[1];
  b_GCmat->size[0] = B->size[0];
  b_GCmat->size[1] = f->size[0];
  emxEnsureCapacity_real_T(b_GCmat, i);
  normvec_data = b_GCmat->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 <= nx; i1++) {
      normvec_data[i1 + b_GCmat->size[0] * i] =
          B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
    }
  }
  i = B->size[0] * B->size[1];
  B->size[0] = b_GCmat->size[0];
  B->size[1] = b_GCmat->size[1];
  emxEnsureCapacity_real_T(B, i);
  B_data = B->data;
  loop_ub = b_GCmat->size[0] * b_GCmat->size[1];
  for (i = 0; i < loop_ub; i++) {
    B_data[i] = normvec_data[i];
  }
  /* 'gkmPWMlasso4:173' B = B/std(B(:))'; */
  nx = B->size[0] * B->size[1];
  b_B = *B;
  b_motclus = nx;
  b_B.size = &b_motclus;
  b_B.numDimensions = 1;
  b1 = b_std(&b_B);
  loop_ub = B->size[0] * B->size[1];
  for (i = 0; i < loop_ub; i++) {
    B_data[i] /= b1;
  }
  emxInit_cell_wrap_1(&tmp_motclus);
  /*  motclus = motclus(f); */
  /* 'gkmPWMlasso4:175' tmp_motclus = cell(length(f),1); */
  nx = f->size[0];
  i = tmp_motclus->size[0];
  tmp_motclus->size[0] = f->size[0];
  emxEnsureCapacity_cell_wrap_1(tmp_motclus, i);
  tmp_motclus_data = tmp_motclus->data;
  for (i = 0; i < nx; i++) {
    tmp_motclus_data[i].f1->size[0] = 1;
    tmp_motclus_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso4:176' tmp_motclus = coder.nullcopy(tmp_motclus); */
  /* 'gkmPWMlasso4:177' for idx = 1:length(f) */
  i = f->size[0];
  for (nx = 0; nx < i; nx++) {
    /* 'gkmPWMlasso4:178' tmp_motclus{idx} = motclus{f(idx)}; */
    i1 = tmp_motclus_data[nx].f1->size[0] * tmp_motclus_data[nx].f1->size[1];
    tmp_motclus_data[nx].f1->size[0] = 1;
    tmp_motclus_data[nx].f1->size[1] =
        motclus_data[(int)f_data[nx] - 1].f1->size[1];
    emxEnsureCapacity_real_T(tmp_motclus_data[nx].f1, i1);
    loop_ub = motclus_data[(int)f_data[nx] - 1].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      tmp_motclus_data[nx].f1->data[i1] =
          motclus_data[(int)f_data[nx] - 1].f1->data[i1];
    }
  }
  emxInit_cell_wrap_1(&c_motclus);
  /* 'gkmPWMlasso4:180' motclus = cell(length(f),1); */
  nx = f->size[0];
  i = c_motclus->size[0];
  c_motclus->size[0] = f->size[0];
  emxEnsureCapacity_cell_wrap_1(c_motclus, i);
  b_motclus_data = c_motclus->data;
  for (i = 0; i < nx; i++) {
    b_motclus_data[i].f1->size[0] = 1;
    b_motclus_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso4:181' motclus = coder.nullcopy(motclus); */
  i = motclus->size[0];
  motclus->size[0] = c_motclus->size[0];
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  /* 'gkmPWMlasso4:182' for idx = 1:length(motclus) */
  i = c_motclus->size[0];
  for (nx = 0; nx < i; nx++) {
    /* 'gkmPWMlasso4:183' motclus{idx} = tmp_motclus{idx}; */
    i1 = motclus_data[nx].f1->size[0] * motclus_data[nx].f1->size[1];
    motclus_data[nx].f1->size[0] = 1;
    motclus_data[nx].f1->size[1] = tmp_motclus_data[nx].f1->size[1];
    emxEnsureCapacity_real_T(motclus_data[nx].f1, i1);
    loop_ub = tmp_motclus_data[nx].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      motclus_data[nx].f1->data[i1] = tmp_motclus_data[nx].f1->data[i1];
    }
  }
  /* 'gkmPWMlasso4:185' Z = Z(f); */
  i = BY->size[0];
  BY->size[0] = f->size[0];
  emxEnsureCapacity_real_T(BY, i);
  BY_data = BY->data;
  loop_ub = f->size[0];
  for (i = 0; i < loop_ub; i++) {
    BY_data[i] = Z_data[(int)f_data[i] - 1];
  }
  i = Z->size[0];
  Z->size[0] = BY->size[0];
  emxEnsureCapacity_real_T(Z, i);
  Z_data = Z->data;
  loop_ub = BY->size[0];
  for (i = 0; i < loop_ub; i++) {
    Z_data[i] = BY_data[i];
  }
  /*  clear A AA loc mat GMmat */
  /* 'gkmPWMlasso4:188' if d == 0 */
  emxInit_real_T(&BB, 2);
  if (varargin_10 == 0.0) {
    /* 'gkmPWMlasso4:189' fprintf('Running LASSO\n'); */
    printf("Running LASSO\n");
    fflush(stdout);
    /*  [weigmat, FitInfo] = lasso_cvmat(B, cfile2,'DFmax',
     * length(Z),'Standardize', false, 'NumLambda', 20); */
    /* 'gkmPWMlasso4:191' [weigmat, FitInfo] = lasso_cvmat(B, cfile2, length(Z),
     * false, 20); */
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = B->size[1];
    emxEnsureCapacity_real_T(comb, i);
    b_comb_data = comb->data;
    loop_ub = B->size[0] * B->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      b_comb_data[i] = B_data[i];
    }
    i = diffc->size[0];
    diffc->size[0] = cfile2->size[0];
    emxEnsureCapacity_real_T(diffc, i);
    b_diffc_data = diffc->data;
    loop_ub = cfile2->size[0] - 1;
    for (i = 0; i <= loop_ub; i++) {
      b_diffc_data[i] = cfile2_data[i];
    }
    emxInit_real_T(&weigmat, 2);
    lasso_cvmat(comb, diffc, Z->size[0], weigmat, expl_temp_data, comb_size,
                b_expl_temp_data, comb2_size, &GCneg1, a__3_DF_data,
                a__3_DF_size, c_expl_temp_data, expl_temp_size);
    xc_data = weigmat->data;
    /* 'gkmPWMlasso4:192' MSE = zeros(length(FitInfo.DF),1); */
    if (a__3_DF_size[1] == 0) {
      u1 = 0;
    } else {
      u1 = a__3_DF_size[1];
    }
    if (0 <= u1 - 1) {
      memset(&MSE_data[0], 0, u1 * sizeof(double));
    }
    /* 'gkmPWMlasso4:193' cnorm = cfile2'*cfile2; */
    if (cfile2->size[0] < 1) {
      GCneg1 = 0.0;
    } else {
      GCneg1 = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                          &cfile2_data[0], (blasint)1);
    }
    emxInit_cell_wrap_3_20(&F);
    /* 'gkmPWMlasso4:194' F = cell(numel(FitInfo.DF),1); */
    nx = a__3_DF_size[1];
    i = F.size[0];
    F.size[0] = a__3_DF_size[1];
    emxEnsureCapacity_cell_wrap_3(F.data, a__3_DF_size[1], i);
    for (i = 0; i < nx; i++) {
      F.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:195' F = coder.nullcopy(F); */
    /* 'gkmPWMlasso4:196' for i = 1:length(FitInfo.DF) */
    if (a__3_DF_size[1] == 0) {
      lk = 0;
    } else {
      lk = a__3_DF_size[1];
    }
    if (0 <= lk - 1) {
      i5 = weigmat->size[0];
      d_loop_ub = weigmat->size[0];
      e_loop_ub = B->size[0];
      f_loop_ub = B->size[0];
      g_loop_ub = B->size[0];
      i6 = B->size[0];
      h_loop_ub = B->size[0];
      i7 = B->size[0];
    }
    for (b_i = 0; b_i < lk; b_i++) {
      /* 'gkmPWMlasso4:197' f = find(weigmat(:,i)~=0); */
      i = b_loc->size[0];
      b_loc->size[0] = i5;
      emxEnsureCapacity_boolean_T(b_loc, i);
      b_loc_data = b_loc->data;
      for (i = 0; i < d_loop_ub; i++) {
        b_loc_data[i] = (xc_data[i + weigmat->size[0] * b_i] != 0.0);
      }
      b_eml_find(b_loc, idx);
      matches_data = idx->data;
      i = f->size[0];
      f->size[0] = idx->size[0];
      emxEnsureCapacity_real_T(f, i);
      f_data = f->data;
      loop_ub = idx->size[0];
      for (i = 0; i < loop_ub; i++) {
        f_data[i] = matches_data[i];
      }
      /* 'gkmPWMlasso4:198' F{i} = f; */
      i = F.data[b_i].f1->size[0];
      F.data[b_i].f1->size[0] = f->size[0];
      emxEnsureCapacity_real_T(F.data[b_i].f1, i);
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        F.data[b_i].f1->data[i] = f_data[i];
      }
      /* 'gkmPWMlasso4:199' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
      i = b_GCmat->size[0] * b_GCmat->size[1];
      b_GCmat->size[0] = e_loop_ub;
      b_GCmat->size[1] = f->size[0];
      emxEnsureCapacity_real_T(b_GCmat, i);
      normvec_data = b_GCmat->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < e_loop_ub; i1++) {
          normvec_data[i1 + b_GCmat->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
        }
      }
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = f_loop_ub;
      comb2->size[1] = f->size[0];
      emxEnsureCapacity_real_T(comb2, i);
      b_comb_data = comb2->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < f_loop_ub; i1++) {
          b_comb_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
        }
      }
      b_mtimes(b_GCmat, comb2, comb);
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = g_loop_ub;
      comb2->size[1] = f->size[0];
      emxEnsureCapacity_real_T(comb2, i);
      b_comb_data = comb2->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < g_loop_ub; i1++) {
          b_comb_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
        }
      }
      if ((i6 == 0) || (f->size[0] == 0) || (cfile2->size[0] == 0)) {
        i = negvec->size[0];
        negvec->size[0] = f->size[0];
        emxEnsureCapacity_real_T(negvec, i);
        negvec_data = negvec->data;
        loop_ub = f->size[0];
        for (i = 0; i < loop_ub; i++) {
          negvec_data[i] = 0.0;
        }
      } else {
        i = negvec->size[0];
        negvec->size[0] = f->size[0];
        emxEnsureCapacity_real_T(negvec, i);
        negvec_data = negvec->data;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    (blasint)f->size[0], (blasint)1, (blasint)B->size[0], 1.0,
                    &b_comb_data[0], (blasint)B->size[0], &cfile2_data[0],
                    (blasint)cfile2->size[0], 0.0, &negvec_data[0],
                    (blasint)f->size[0]);
      }
      mpower(comb, comb2);
      b_comb_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (negvec->size[0] == 0)) {
        i = diffc->size[0];
        diffc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(diffc, i);
        b_diffc_data = diffc->data;
        loop_ub = comb2->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_diffc_data[i] = 0.0;
        }
      } else {
        i = diffc->size[0];
        diffc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(diffc, i);
        b_diffc_data = diffc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &b_comb_data[0],
            (blasint)comb2->size[0], &negvec_data[0], (blasint)negvec->size[0],
            0.0, &b_diffc_data[0], (blasint)comb2->size[0]);
      }
      /* 'gkmPWMlasso4:200' res = B(:,f)*OLS; */
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = h_loop_ub;
      comb2->size[1] = f->size[0];
      emxEnsureCapacity_real_T(comb2, i);
      b_comb_data = comb2->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < h_loop_ub; i1++) {
          b_comb_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
        }
      }
      if ((i7 == 0) || (f->size[0] == 0) || (diffc->size[0] == 0)) {
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
                    (blasint)B->size[0], (blasint)1, (blasint)f->size[0], 1.0,
                    &b_comb_data[0], (blasint)B->size[0], &b_diffc_data[0],
                    (blasint)diffc->size[0], 0.0, &loc_data[0],
                    (blasint)B->size[0]);
      }
      /* 'gkmPWMlasso4:201' if sum(abs(res)>0) */
      nx = loc->size[0];
      i = negvec->size[0];
      negvec->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      for (k = 0; k < nx; k++) {
        negvec_data[k] = fabs(loc_data[k]);
      }
      i = b_loc->size[0];
      b_loc->size[0] = negvec->size[0];
      emxEnsureCapacity_boolean_T(b_loc, i);
      b_loc_data = b_loc->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_loc_data[i] = (negvec_data[i] > 0.0);
      }
      nx = b_loc->size[0];
      if (b_loc->size[0] == 0) {
        indc_size = 0;
      } else {
        indc_size = b_loc_data[0];
        for (k = 2; k <= nx; k++) {
          indc_size += b_loc_data[k - 1];
        }
      }
      if (indc_size != 0) {
        /* 'gkmPWMlasso4:202' MSE(i) = (cfile2'*res)^2/(res'*res)/cnorm; */
        if (cfile2->size[0] < 1) {
          a = 0.0;
        } else {
          a = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                         &loc_data[0], (blasint)1);
        }
        if (loc->size[0] < 1) {
          b1 = 0.0;
        } else {
          b1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                          &loc_data[0], (blasint)1);
        }
        MSE_data[b_i] = a * a / b1 / GCneg1;
      }
    }
    emxFree_real_T(&weigmat);
    emxInit_cell_wrap_3_1x19(&FF);
    /* 'gkmPWMlasso4:205' FF = cell(1,length(MSE)-1); */
    nx = u1 - 1;
    i = FF.size[0] * FF.size[1];
    FF.size[0] = 1;
    FF.size[1] = u1 - 1;
    emxEnsureCapacity_cell_wrap_31(FF.data, FF.size, i);
    for (i = 0; i < nx; i++) {
      FF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:206' FF = coder.nullcopy(FF); */
    /* 'gkmPWMlasso4:207' MSE2 = zeros(1,length(MSE)-1); */
    MSE2_size[0] = 1;
    if (0 <= nx - 1) {
      memset(&MSE2_data[0], 0, nx * sizeof(double));
    }
    /* 'gkmPWMlasso4:208' count = 1; */
    lk = 0;
    /* 'gkmPWMlasso4:209' for i = 1:length(MSE)-1 */
    for (b_i = 0; b_i <= u1 - 2; b_i++) {
      /* 'gkmPWMlasso4:210' if numel(setdiff(F{i}, F{i+1})) > 0 */
      b_do_vectors(F.data[b_i].f1, F.data[b_i + 1].f1, loc, idx, &nx);
      if (loc->size[0] > 0) {
        /* 'gkmPWMlasso4:211' FF{count} = setdiff(F{i}, F{i+1}); */
        b_do_vectors(F.data[b_i].f1, F.data[b_i + 1].f1, FF.data[lk].f1, idx,
                     &nx);
        /* 'gkmPWMlasso4:212' MSE2(count) = (MSE(i)-MSE(i+1)); */
        MSE2_data[lk] = MSE_data[b_i] - MSE_data[b_i + 1];
        /* 'gkmPWMlasso4:213' count = count +1; */
        lk++;
      }
    }
    emxFree_cell_wrap_3_20(&F);
    emxInit_cell_wrap_3_1x19(&newFF);
    /* 'gkmPWMlasso4:216' MSE2 = MSE2(1:count-1); */
    if (1 > lk) {
      MSE2_size[1] = 0;
    } else {
      MSE2_size[1] = lk;
    }
    /* 'gkmPWMlasso4:217' newFF = cell(1,count-1); */
    i = newFF.size[0] * newFF.size[1];
    newFF.size[0] = 1;
    newFF.size[1] = lk;
    emxEnsureCapacity_cell_wrap_31(newFF.data, newFF.size, i);
    for (i = 0; i < lk; i++) {
      newFF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:218' newFF = coder.nullcopy(newFF); */
    /* 'gkmPWMlasso4:219' for idx = 1:count-1 */
    for (nx = 0; nx < lk; nx++) {
      /* 'gkmPWMlasso4:220' newFF{idx} = FF{idx}; */
      loop_ub = FF.data[nx].f1->size[0];
      i = newFF.data[nx].f1->size[0];
      newFF.data[nx].f1->size[0] = FF.data[nx].f1->size[0];
      emxEnsureCapacity_real_T(newFF.data[nx].f1, i);
      for (i = 0; i < loop_ub; i++) {
        newFF.data[nx].f1->data[i] = FF.data[nx].f1->data[i];
      }
    }
    emxInit_cell_wrap_3_1x19(&b_FF);
    /* 'gkmPWMlasso4:222' FF = cell(1,count-1); */
    i = b_FF.size[0] * b_FF.size[1];
    b_FF.size[0] = 1;
    b_FF.size[1] = lk;
    emxEnsureCapacity_cell_wrap_31(b_FF.data, b_FF.size, i);
    for (i = 0; i < lk; i++) {
      b_FF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:223' FF = coder.nullcopy(FF); */
    i = FF.size[0] * FF.size[1];
    FF.size[0] = 1;
    FF.size[1] = b_FF.size[1];
    emxEnsureCapacity_cell_wrap_31(FF.data, FF.size, i);
    /* 'gkmPWMlasso4:224' for idx = 1:count-1 */
    emxFree_cell_wrap_3_1x19(&b_FF);
    for (nx = 0; nx < lk; nx++) {
      /* 'gkmPWMlasso4:225' FF{idx} = newFF{idx}; */
      loop_ub = newFF.data[nx].f1->size[0];
      i = FF.data[nx].f1->size[0];
      FF.data[nx].f1->size[0] = newFF.data[nx].f1->size[0];
      emxEnsureCapacity_real_T(FF.data[nx].f1, i);
      for (i = 0; i < loop_ub; i++) {
        FF.data[nx].f1->data[i] = newFF.data[nx].f1->data[i];
      }
    }
    /* 'gkmPWMlasso4:228' res = B*((B.'*B)^-1*(B.'*cfile2)); */
    b_mtimes(B, B, comb);
    if ((B->size[0] == 0) || (B->size[1] == 0) || (cfile2->size[0] == 0)) {
      i = negvec->size[0];
      negvec->size[0] = B->size[1];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      loop_ub = B->size[1];
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = 0.0;
      }
    } else {
      i = negvec->size[0];
      negvec->size[0] = B->size[1];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)B->size[1],
                  (blasint)1, (blasint)B->size[0], 1.0, &B_data[0],
                  (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &negvec_data[0],
                  (blasint)B->size[1]);
    }
    mpower(comb, comb2);
    b_comb_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (negvec->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = comb2->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &b_comb_data[0], (blasint)comb2->size[0],
                  &negvec_data[0], (blasint)negvec->size[0], 0.0, &BY_data[0],
                  (blasint)comb2->size[0]);
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
    /* 'gkmPWMlasso4:229' csm = (cfile2'*res)^2/(res'*res)/cnorm; */
    if (cfile2->size[0] < 1) {
      a = 0.0;
    } else {
      a = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                     &loc_data[0], (blasint)1);
    }
    if (loc->size[0] < 1) {
      b1 = 0.0;
    } else {
      b1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                      &loc_data[0], (blasint)1);
    }
    GCneg1 = a * a / b1 / GCneg1;
    /* 'gkmPWMlasso4:230' [a,b] = sort(MSE2,'descend'); */
    c_sort(MSE2_data, MSE2_size, iidx_data, comb_size);
    /* 'gkmPWMlasso4:231' cs = cumsum(a); */
    useConstantDim(MSE2_data, MSE2_size);
    /* 'gkmPWMlasso4:232' cs = cs/csm; */
    loop_ub = MSE2_size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      MSE2_data[i] /= GCneg1;
    }
    /* 'gkmPWMlasso4:233' f = find(cs>0.9); */
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
    loop_ub = match_out->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_f_data[i] = match_out_data[i];
    }
    emxInit_cell_wrap_3_1x19(&b_newFF);
    /*  FF = FF(b(1:f(1))); */
    /* 'gkmPWMlasso4:235' endIdx = f(1); */
    /* 'gkmPWMlasso4:236' newFF = cell(1, endIdx); */
    nx = b_f_data[0];
    i = b_newFF.size[0] * b_newFF.size[1];
    b_newFF.size[0] = 1;
    b_newFF.size[1] = b_f_data[0];
    emxEnsureCapacity_cell_wrap_31(b_newFF.data, b_newFF.size, i);
    for (i = 0; i < nx; i++) {
      b_newFF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:237' newFF = coder.nullcopy(newFF); */
    i = newFF.size[0] * newFF.size[1];
    newFF.size[0] = 1;
    newFF.size[1] = b_newFF.size[1];
    emxEnsureCapacity_cell_wrap_31(newFF.data, newFF.size, i);
    /* 'gkmPWMlasso4:238' for idx = 1:endIdx */
    i = b_f_data[0];
    emxFree_cell_wrap_3_1x19(&b_newFF);
    for (nx = 0; nx < i; nx++) {
      /* 'gkmPWMlasso4:239' newFF{idx} = FF{b(idx)}; */
      i1 = newFF.data[nx].f1->size[0];
      i2 = iidx_data[nx];
      newFF.data[nx].f1->size[0] = FF.data[i2 - 1].f1->size[0];
      emxEnsureCapacity_real_T(newFF.data[nx].f1, i1);
      loop_ub = FF.data[i2 - 1].f1->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        newFF.data[nx].f1->data[i1] = FF.data[iidx_data[nx] - 1].f1->data[i1];
      }
    }
    emxFree_cell_wrap_3_1x19(&FF);
    /* 'gkmPWMlasso4:241' FF = newFF; */
    /* 'gkmPWMlasso4:243' f = []; */
    f->size[0] = 0;
    /* 'gkmPWMlasso4:244' for i = 1:length(FF) */
    i = newFF.size[1];
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:245' f = [f;FF{i}]; */
      i1 = f->size[0];
      loop_ub = newFF.data[b_i].f1->size[0];
      i2 = f->size[0];
      f->size[0] += newFF.data[b_i].f1->size[0];
      emxEnsureCapacity_real_T(f, i2);
      f_data = f->data;
      for (i2 = 0; i2 < loop_ub; i2++) {
        f_data[i1 + i2] = newFF.data[b_i].f1->data[i2];
      }
    }
    emxFree_cell_wrap_3_1x19(&newFF);
    /* 'gkmPWMlasso4:247' f=unique(f); */
    i = diffc->size[0];
    diffc->size[0] = f->size[0];
    emxEnsureCapacity_real_T(diffc, i);
    b_diffc_data = diffc->data;
    loop_ub = f->size[0] - 1;
    for (i = 0; i <= loop_ub; i++) {
      b_diffc_data[i] = f_data[i];
    }
    unique_vector(diffc, f);
    f_data = f->data;
    /* f = find(weigmat(:,1)~=0); */
    /* 'gkmPWMlasso4:249' F = length(f); */
    /* 'gkmPWMlasso4:250' fprintf('Selecting top motifs\n'); */
    printf("Selecting top motifs\n");
    fflush(stdout);
    /* 'gkmPWMlasso4:251' ind = true; */
    empty_non_axis_sizes = true;
    /* 'gkmPWMlasso4:252' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = f->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, comb);
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = negvec->size[0];
      negvec->size[0] = f->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = 0.0;
      }
    } else {
      i = negvec->size[0];
      negvec->size[0] = f->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)f->size[0],
                  (blasint)1, (blasint)B->size[0], 1.0, &b_comb_data[0],
                  (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &negvec_data[0],
                  (blasint)f->size[0]);
    }
    mpower(comb, comb2);
    b_comb_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (negvec->size[0] == 0)) {
      i = diffc->size[0];
      diffc->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(diffc, i);
      b_diffc_data = diffc->data;
      loop_ub = comb2->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_diffc_data[i] = 0.0;
      }
    } else {
      i = diffc->size[0];
      diffc->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(diffc, i);
      b_diffc_data = diffc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &b_comb_data[0], (blasint)comb2->size[0],
                  &negvec_data[0], (blasint)negvec->size[0], 0.0,
                  &b_diffc_data[0], (blasint)comb2->size[0]);
    }
    /* 'gkmPWMlasso4:253' Pweig = Z(f); */
    i = indc->size[0];
    indc->size[0] = f->size[0];
    emxEnsureCapacity_real_T(indc, i);
    b_indc_data = indc->data;
    loop_ub = f->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_indc_data[i] = Z_data[(int)f_data[i] - 1];
    }
    /* 'gkmPWMlasso4:254' while ind */
    while (empty_non_axis_sizes) {
      /* 'gkmPWMlasso4:255' ff = []; */
      loc->size[0] = 0;
      /* 'gkmPWMlasso4:256' for i = 1:length(f) */
      i = f->size[0];
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso4:257' if sign(OLS(i)) ~= sign(Pweig(i)) */
        b1 = b_diffc_data[b_i];
        if (b_diffc_data[b_i] < 0.0) {
          b1 = -1.0;
        } else if (b_diffc_data[b_i] > 0.0) {
          b1 = 1.0;
        }
        GCneg1 = b_indc_data[b_i];
        if (b_indc_data[b_i] < 0.0) {
          GCneg1 = -1.0;
        } else if (b_indc_data[b_i] > 0.0) {
          GCneg1 = 1.0;
        }
        if (b1 != GCneg1) {
          /* 'gkmPWMlasso4:258' ff = [ff;i]; */
          i1 = loc->size[0];
          i2 = loc->size[0];
          loc->size[0]++;
          emxEnsureCapacity_real_T(loc, i2);
          loc_data = loc->data;
          loc_data[i1] = (double)b_i + 1.0;
        }
      }
      /* 'gkmPWMlasso4:261' if length(ff) > 0 */
      if (loc->size[0] > 0) {
        /* 'gkmPWMlasso4:262' f(ff) = []; */
        i = idx->size[0];
        idx->size[0] = loc->size[0];
        emxEnsureCapacity_int32_T(idx, i);
        matches_data = idx->data;
        loop_ub = loc->size[0];
        for (i = 0; i < loop_ub; i++) {
          matches_data[i] = (int)loc_data[i];
        }
        b_nullAssignment(f, idx);
        f_data = f->data;
        /* 'gkmPWMlasso4:263' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
        loop_ub = B->size[0];
        i = b_GCmat->size[0] * b_GCmat->size[1];
        b_GCmat->size[0] = B->size[0];
        b_GCmat->size[1] = f->size[0];
        emxEnsureCapacity_real_T(b_GCmat, i);
        normvec_data = b_GCmat->data;
        nx = f->size[0];
        i = comb2->size[0] * comb2->size[1];
        comb2->size[0] = B->size[0];
        comb2->size[1] = f->size[0];
        emxEnsureCapacity_real_T(comb2, i);
        b_comb_data = comb2->data;
        for (i = 0; i < nx; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            GCneg1 = B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
            normvec_data[i1 + b_GCmat->size[0] * i] = GCneg1;
            b_comb_data[i1 + comb2->size[0] * i] = GCneg1;
          }
        }
        b_mtimes(b_GCmat, comb2, comb);
        loop_ub = B->size[0];
        i = comb2->size[0] * comb2->size[1];
        comb2->size[0] = B->size[0];
        comb2->size[1] = f->size[0];
        emxEnsureCapacity_real_T(comb2, i);
        b_comb_data = comb2->data;
        nx = f->size[0];
        for (i = 0; i < nx; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            b_comb_data[i1 + comb2->size[0] * i] =
                B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
          }
        }
        if ((B->size[0] == 0) || (f->size[0] == 0) || (cfile2->size[0] == 0)) {
          i = negvec->size[0];
          negvec->size[0] = f->size[0];
          emxEnsureCapacity_real_T(negvec, i);
          negvec_data = negvec->data;
          loop_ub = f->size[0];
          for (i = 0; i < loop_ub; i++) {
            negvec_data[i] = 0.0;
          }
        } else {
          i = negvec->size[0];
          negvec->size[0] = f->size[0];
          emxEnsureCapacity_real_T(negvec, i);
          negvec_data = negvec->data;
          cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                      (blasint)f->size[0], (blasint)1, (blasint)B->size[0], 1.0,
                      &b_comb_data[0], (blasint)B->size[0], &cfile2_data[0],
                      (blasint)cfile2->size[0], 0.0, &negvec_data[0],
                      (blasint)f->size[0]);
        }
        mpower(comb, comb2);
        b_comb_data = comb2->data;
        if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
            (negvec->size[0] == 0)) {
          i = diffc->size[0];
          diffc->size[0] = comb2->size[0];
          emxEnsureCapacity_real_T(diffc, i);
          b_diffc_data = diffc->data;
          loop_ub = comb2->size[0];
          for (i = 0; i < loop_ub; i++) {
            b_diffc_data[i] = 0.0;
          }
        } else {
          i = diffc->size[0];
          diffc->size[0] = comb2->size[0];
          emxEnsureCapacity_real_T(diffc, i);
          b_diffc_data = diffc->data;
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                      (blasint)comb2->size[0], (blasint)1,
                      (blasint)comb2->size[1], 1.0, &b_comb_data[0],
                      (blasint)comb2->size[0], &negvec_data[0],
                      (blasint)negvec->size[0], 0.0, &b_diffc_data[0],
                      (blasint)comb2->size[0]);
        }
        /* 'gkmPWMlasso4:264' Pweig = Z(f); */
        i = indc->size[0];
        indc->size[0] = f->size[0];
        emxEnsureCapacity_real_T(indc, i);
        b_indc_data = indc->data;
        loop_ub = f->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_indc_data[i] = Z_data[(int)f_data[i] - 1];
        }
      } else {
        /* 'gkmPWMlasso4:265' else */
        /* 'gkmPWMlasso4:266' ind = false; */
        empty_non_axis_sizes = false;
      }
    }
    /* 'gkmPWMlasso4:269' BB = B(:,f); */
    loop_ub = B->size[0];
    i = BB->size[0] * BB->size[1];
    BB->size[0] = B->size[0];
    BB->size[1] = f->size[0];
    emxEnsureCapacity_real_T(BB, i);
    AA_data = BB->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        AA_data[i1 + BB->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    /* 'gkmPWMlasso4:270' BX = B(:,f)'*B(:,f); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = f->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, xc);
    xc_data = xc->data;
    /* 'gkmPWMlasso4:271' BY = B(:,f)'*cfile2; */
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)f->size[0],
                  (blasint)1, (blasint)B->size[0], 1.0, &b_comb_data[0],
                  (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)f->size[0]);
    }
    /* 'gkmPWMlasso4:272' E = zeros(length(f),1); */
    i = normvec->size[0];
    normvec->size[0] = f->size[0];
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    loop_ub = f->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso4:273' for i = 1:length(f) */
    i = f->size[0];
    if (0 <= f->size[0] - 1) {
      i_loop_ub = xc->size[0] * xc->size[1];
      j_loop_ub = BY->size[0];
      nxout = BY->size[0] - 1;
      k_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:274' B = BB; */
      /* 'gkmPWMlasso4:275' BBX = BX; */
      /* 'gkmPWMlasso4:276' BBY = BY; */
      /* 'gkmPWMlasso4:277' B(:,i) = []; */
      /* 'gkmPWMlasso4:278' BBX(:,i) = []; */
      /* 'gkmPWMlasso4:279' BBX(i,:) = []; */
      /* 'gkmPWMlasso4:280' BBY(i) = []; */
      /* 'gkmPWMlasso4:281' res = cfile2-B*(BBX^-1*BBY); */
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = xc->size[0];
      comb->size[1] = xc->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      b_comb_data = comb->data;
      for (i1 = 0; i1 < i_loop_ub; i1++) {
        b_comb_data[i1] = xc_data[i1];
      }
      c_nullAssignment(comb, b_i + 1);
      nullAssignment(comb, b_i + 1);
      mpower(comb, comb2);
      b_comb_data = comb2->data;
      nx = b_i + 1;
      i1 = loc->size[0];
      loc->size[0] = BY->size[0];
      emxEnsureCapacity_real_T(loc, i1);
      loc_data = loc->data;
      for (i1 = 0; i1 < j_loop_ub; i1++) {
        loc_data[i1] = BY_data[i1];
      }
      for (k = nx; k <= nxout; k++) {
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
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (loc->size[0] == 0)) {
        i1 = negvec->size[0];
        negvec->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(negvec, i1);
        negvec_data = negvec->data;
        loop_ub = comb2->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          negvec_data[i1] = 0.0;
        }
      } else {
        i1 = negvec->size[0];
        negvec->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(negvec, i1);
        negvec_data = negvec->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &b_comb_data[0],
            (blasint)comb2->size[0], &loc_data[0], (blasint)loc->size[0], 0.0,
            &negvec_data[0], (blasint)comb2->size[0]);
      }
      i1 = comb2->size[0] * comb2->size[1];
      comb2->size[0] = BB->size[0];
      comb2->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb2, i1);
      b_comb_data = comb2->data;
      for (i1 = 0; i1 < k_loop_ub; i1++) {
        b_comb_data[i1] = AA_data[i1];
      }
      c_nullAssignment(comb2, b_i + 1);
      b_comb_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (negvec->size[0] == 0)) {
        i1 = loc->size[0];
        loc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        loop_ub = comb2->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          loc_data[i1] = 0.0;
        }
      } else {
        i1 = loc->size[0];
        loc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &b_comb_data[0],
            (blasint)comb2->size[0], &negvec_data[0], (blasint)negvec->size[0],
            0.0, &loc_data[0], (blasint)comb2->size[0]);
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
      /* 'gkmPWMlasso4:282' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        b1 = 0.0;
      } else {
        b1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                        &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(b1);
    }
    /* 'gkmPWMlasso4:284' res = cfile2-BB*OLS; */
    if ((B->size[0] == 0) || (f->size[0] == 0) || (diffc->size[0] == 0)) {
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
                  (blasint)B->size[0], (blasint)1, (blasint)f->size[0], 1.0,
                  &AA_data[0], (blasint)B->size[0], &b_diffc_data[0],
                  (blasint)diffc->size[0], 0.0, &loc_data[0],
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
    /* 'gkmPWMlasso4:285' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      b1 = 0.0;
    } else {
      b1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                      &loc_data[0], (blasint)1);
    }
    GCneg1 = sqrt(b1);
    /* 'gkmPWMlasso4:286' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - GCneg1) / GCneg1;
    }
    /*  motclus = motclus(f); */
    /* 'gkmPWMlasso4:288' mylen = length(f); */
    /* 'gkmPWMlasso4:289' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso4:290' for idx=1:mylen */
    i = f->size[0];
    i1 = c_motclus->size[0];
    c_motclus->size[0] = f->size[0];
    emxEnsureCapacity_cell_wrap_1(c_motclus, i1);
    b_motclus_data = c_motclus->data;
    for (nx = 0; nx < i; nx++) {
      /* 'gkmPWMlasso4:291' newMotclus{idx} = motclus{f(idx)}; */
      i1 = b_motclus_data[nx].f1->size[0] * b_motclus_data[nx].f1->size[1];
      b_motclus_data[nx].f1->size[0] = 1;
      b_motclus_data[nx].f1->size[1] =
          motclus_data[(int)f_data[nx] - 1].f1->size[1];
      emxEnsureCapacity_real_T(b_motclus_data[nx].f1, i1);
      loop_ub = motclus_data[(int)f_data[nx] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_motclus_data[nx].f1->data[i1] =
            motclus_data[(int)f_data[nx] - 1].f1->data[i1];
      }
    }
    /* 'gkmPWMlasso4:293' motclus = newMotclus; */
    /* 'gkmPWMlasso4:295' f = find(E/max(E) >= 0.01); */
    GCneg1 = maximum(normvec);
    i = b_loc->size[0];
    b_loc->size[0] = normvec->size[0];
    emxEnsureCapacity_boolean_T(b_loc, i);
    b_loc_data = b_loc->data;
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_loc_data[i] = (normvec_data[i] / GCneg1 >= 0.01);
    }
    b_eml_find(b_loc, idx);
    matches_data = idx->data;
    i = f->size[0];
    f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso4:296' B = BB; */
    i = B->size[0] * B->size[1];
    B->size[0] = BB->size[0];
    B->size[1] = BB->size[1];
    emxEnsureCapacity_real_T(B, i);
    B_data = B->data;
    loop_ub = BB->size[0] * BB->size[1];
    for (i = 0; i < loop_ub; i++) {
      B_data[i] = AA_data[i];
    }
    /* 'gkmPWMlasso4:297' BB = B(:,f); */
    nx = BB->size[0] - 1;
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = BB->size[0];
    b_GCmat->size[1] = f->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    loop_ub = f->size[0];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 <= nx; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            AA_data[i1 + BB->size[0] * ((int)f_data[i] - 1)];
      }
    }
    i = BB->size[0] * BB->size[1];
    BB->size[0] = b_GCmat->size[0];
    BB->size[1] = b_GCmat->size[1];
    emxEnsureCapacity_real_T(BB, i);
    AA_data = BB->data;
    loop_ub = b_GCmat->size[0] * b_GCmat->size[1];
    for (i = 0; i < loop_ub; i++) {
      AA_data[i] = normvec_data[i];
    }
    /* 'gkmPWMlasso4:298' BX = B(:,f)'*B(:,f); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = f->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, xc);
    xc_data = xc->data;
    /* 'gkmPWMlasso4:299' BY = B(:,f)'*cfile2; */
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)f->size[0],
                  (blasint)1, (blasint)B->size[0], 1.0, &b_comb_data[0],
                  (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)f->size[0]);
    }
    /* 'gkmPWMlasso4:300' E = zeros(length(f),1); */
    i = normvec->size[0];
    normvec->size[0] = f->size[0];
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    loop_ub = f->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso4:301' for i = 1:length(f) */
    i = f->size[0];
    if (0 <= f->size[0] - 1) {
      l_loop_ub = xc->size[0] * xc->size[1];
      m_loop_ub = BY->size[0];
      b_nxout = BY->size[0] - 1;
      n_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:302' B = BB; */
      /* 'gkmPWMlasso4:303' BBX = BX; */
      /* 'gkmPWMlasso4:304' BBY = BY; */
      /* 'gkmPWMlasso4:305' B(:,i) = []; */
      /* 'gkmPWMlasso4:306' BBX(:,i) = []; */
      /* 'gkmPWMlasso4:307' BBX(i,:) = []; */
      /* 'gkmPWMlasso4:308' BBY(i) = []; */
      /* 'gkmPWMlasso4:309' res = cfile2-B*(BBX^-1*BBY); */
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = xc->size[0];
      comb->size[1] = xc->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      b_comb_data = comb->data;
      for (i1 = 0; i1 < l_loop_ub; i1++) {
        b_comb_data[i1] = xc_data[i1];
      }
      c_nullAssignment(comb, b_i + 1);
      nullAssignment(comb, b_i + 1);
      mpower(comb, comb2);
      b_comb_data = comb2->data;
      nx = b_i + 1;
      i1 = loc->size[0];
      loc->size[0] = BY->size[0];
      emxEnsureCapacity_real_T(loc, i1);
      loc_data = loc->data;
      for (i1 = 0; i1 < m_loop_ub; i1++) {
        loc_data[i1] = BY_data[i1];
      }
      for (k = nx; k <= b_nxout; k++) {
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
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (loc->size[0] == 0)) {
        i1 = negvec->size[0];
        negvec->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(negvec, i1);
        negvec_data = negvec->data;
        loop_ub = comb2->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          negvec_data[i1] = 0.0;
        }
      } else {
        i1 = negvec->size[0];
        negvec->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(negvec, i1);
        negvec_data = negvec->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &b_comb_data[0],
            (blasint)comb2->size[0], &loc_data[0], (blasint)loc->size[0], 0.0,
            &negvec_data[0], (blasint)comb2->size[0]);
      }
      i1 = comb2->size[0] * comb2->size[1];
      comb2->size[0] = BB->size[0];
      comb2->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb2, i1);
      b_comb_data = comb2->data;
      for (i1 = 0; i1 < n_loop_ub; i1++) {
        b_comb_data[i1] = AA_data[i1];
      }
      c_nullAssignment(comb2, b_i + 1);
      b_comb_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (negvec->size[0] == 0)) {
        i1 = loc->size[0];
        loc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        loop_ub = comb2->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          loc_data[i1] = 0.0;
        }
      } else {
        i1 = loc->size[0];
        loc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &b_comb_data[0],
            (blasint)comb2->size[0], &negvec_data[0], (blasint)negvec->size[0],
            0.0, &loc_data[0], (blasint)comb2->size[0]);
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
      /* 'gkmPWMlasso4:310' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        b1 = 0.0;
      } else {
        b1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                        &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(b1);
    }
    /* 'gkmPWMlasso4:312' OLS = (BB.'*BB)^-1*(BB.'*cfile2); */
    b_mtimes(BB, BB, comb);
    if ((BB->size[0] == 0) || (BB->size[1] == 0) || (cfile2->size[0] == 0)) {
      i = negvec->size[0];
      negvec->size[0] = BB->size[1];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      loop_ub = BB->size[1];
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = 0.0;
      }
    } else {
      i = negvec->size[0];
      negvec->size[0] = BB->size[1];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)BB->size[1],
                  (blasint)1, (blasint)BB->size[0], 1.0, &AA_data[0],
                  (blasint)BB->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &negvec_data[0],
                  (blasint)BB->size[1]);
    }
    mpower(comb, comb2);
    b_comb_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (negvec->size[0] == 0)) {
      i = diffc->size[0];
      diffc->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(diffc, i);
      b_diffc_data = diffc->data;
      loop_ub = comb2->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_diffc_data[i] = 0.0;
      }
    } else {
      i = diffc->size[0];
      diffc->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(diffc, i);
      b_diffc_data = diffc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &b_comb_data[0], (blasint)comb2->size[0],
                  &negvec_data[0], (blasint)negvec->size[0], 0.0,
                  &b_diffc_data[0], (blasint)comb2->size[0]);
    }
    /* 'gkmPWMlasso4:313' Pweig = Pweig(f); */
    i = BY->size[0];
    BY->size[0] = f->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    loop_ub = f->size[0];
    for (i = 0; i < loop_ub; i++) {
      BY_data[i] = b_indc_data[(int)f_data[i] - 1];
    }
    i = indc->size[0];
    indc->size[0] = BY->size[0];
    emxEnsureCapacity_real_T(indc, i);
    b_indc_data = indc->data;
    loop_ub = BY->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_indc_data[i] = BY_data[i];
    }
    /* 'gkmPWMlasso4:314' res = cfile2-BB*OLS; */
    if ((BB->size[0] == 0) || (BB->size[1] == 0) || (diffc->size[0] == 0)) {
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
                  &AA_data[0], (blasint)BB->size[0], &b_diffc_data[0],
                  (blasint)diffc->size[0], 0.0, &loc_data[0],
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
    /* 'gkmPWMlasso4:315' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      b1 = 0.0;
    } else {
      b1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                      &loc_data[0], (blasint)1);
    }
    GCneg1 = sqrt(b1);
    /* 'gkmPWMlasso4:316' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - GCneg1) / GCneg1;
    }
    /* 'gkmPWMlasso4:317' for i = 1:length(motclus) */
    i = c_motclus->size[0];
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:318' motclus{i} = indvec(motclus{i})'; */
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
    /* 'gkmPWMlasso4:321' mylen = length(f); */
    /* 'gkmPWMlasso4:322' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso4:323' for idx=1:mylen */
    i = f->size[0];
    i1 = tmp_motclus->size[0];
    tmp_motclus->size[0] = f->size[0];
    emxEnsureCapacity_cell_wrap_1(tmp_motclus, i1);
    tmp_motclus_data = tmp_motclus->data;
    for (nx = 0; nx < i; nx++) {
      /* 'gkmPWMlasso4:324' newMotclus{idx} = motclus{f(idx)}; */
      i1 = tmp_motclus_data[nx].f1->size[0] * tmp_motclus_data[nx].f1->size[1];
      tmp_motclus_data[nx].f1->size[0] = 1;
      tmp_motclus_data[nx].f1->size[1] =
          b_motclus_data[(int)f_data[nx] - 1].f1->size[1];
      emxEnsureCapacity_real_T(tmp_motclus_data[nx].f1, i1);
      loop_ub = b_motclus_data[(int)f_data[nx] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        tmp_motclus_data[nx].f1->data[i1] =
            b_motclus_data[(int)f_data[nx] - 1].f1->data[i1];
      }
    }
    /* 'gkmPWMlasso4:326' motclus = newMotclus; */
    /* 'gkmPWMlasso4:328' correlation = corrcoef(cfile2, BB*OLS); */
    if ((BB->size[0] == 0) || (BB->size[1] == 0) || (diffc->size[0] == 0)) {
      i = negvec->size[0];
      negvec->size[0] = BB->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      loop_ub = BB->size[0];
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = 0.0;
      }
    } else {
      i = negvec->size[0];
      negvec->size[0] = BB->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)BB->size[0], (blasint)1, (blasint)BB->size[1], 1.0,
                  &AA_data[0], (blasint)BB->size[0], &b_diffc_data[0],
                  (blasint)diffc->size[0], 0.0, &negvec_data[0],
                  (blasint)BB->size[0]);
    }
    /* 'gkmPWMlasso4:329' gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus,
     * sprintf("%s_%d_%d", filename, int32(l_svm2),
     * int32(k_svm2)),memefile,num,minL, minInfo, correlation(1,2)); */
    GCneg1 = maximum(diffc);
    GCpos1 = maximum(normvec);
    b_sprintf(varargin_1, (int)rt_roundd(varargin_6),
              (int)rt_roundd(varargin_7), text);
    loop_ub = diffc->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_diffc_data[i] /= GCneg1;
    }
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] /= GCpos1;
    }
    corrcoef(cfile2, negvec, dv);
    gettopmotifs(diffc, indc, normvec, tmp_motclus, text, varargin_2,
                 size_tmp_idx_1, varargin_3, varargin_4, dv[2]);
  } else {
    /* 'gkmPWMlasso4:331' else */
    /* 'gkmPWMlasso4:333' fprintf('Running LASSO (1)\n'); */
    printf("Running LASSO (1)\n");
    fflush(stdout);
    /*  weigmat = lasso_cvmat(B, cfile2,'DFmax', d,'Standardize', false); */
    /* 'gkmPWMlasso4:335' [weigmat, ~] = lasso_cvmat(B, cfile2, d, false, 100);
     */
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = B->size[1];
    emxEnsureCapacity_real_T(comb, i);
    b_comb_data = comb->data;
    loop_ub = B->size[0] * B->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      b_comb_data[i] = B_data[i];
    }
    i = diffc->size[0];
    diffc->size[0] = cfile2->size[0];
    emxEnsureCapacity_real_T(diffc, i);
    b_diffc_data = diffc->data;
    loop_ub = cfile2->size[0] - 1;
    for (i = 0; i <= loop_ub; i++) {
      b_diffc_data[i] = cfile2_data[i];
    }
    emxInit_real_T(&weigmat, 2);
    b_lasso_cvmat(comb, diffc, varargin_10, weigmat, expl_temp_data, comb_size,
                  b_expl_temp_data, comb2_size, &GCneg1, a__3_DF_data,
                  a__3_DF_size, c_expl_temp_data, expl_temp_size);
    xc_data = weigmat->data;
    /* 'gkmPWMlasso4:337' f = find(weigmat(:,1)~=0); */
    loop_ub = weigmat->size[0];
    i = b_loc->size[0];
    b_loc->size[0] = weigmat->size[0];
    emxEnsureCapacity_boolean_T(b_loc, i);
    b_loc_data = b_loc->data;
    for (i = 0; i < loop_ub; i++) {
      b_loc_data[i] = (xc_data[i] != 0.0);
    }
    b_eml_find(b_loc, idx);
    matches_data = idx->data;
    i = f->size[0];
    f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso4:338' cfile3 = cfile2 -
     * B(:,f)*((B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2)); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = f->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, comb);
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = negvec->size[0];
      negvec->size[0] = f->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = 0.0;
      }
    } else {
      i = negvec->size[0];
      negvec->size[0] = f->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)f->size[0],
                  (blasint)1, (blasint)B->size[0], 1.0, &b_comb_data[0],
                  (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &negvec_data[0],
                  (blasint)f->size[0]);
    }
    mpower(comb, comb2);
    b_comb_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (negvec->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = comb2->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &b_comb_data[0], (blasint)comb2->size[0],
                  &negvec_data[0], (blasint)negvec->size[0], 0.0, &BY_data[0],
                  (blasint)comb2->size[0]);
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (f->size[0] == 0) || (BY->size[0] == 0)) {
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
                  (blasint)B->size[0], (blasint)1, (blasint)f->size[0], 1.0,
                  &b_comb_data[0], (blasint)B->size[0], &BY_data[0],
                  (blasint)BY->size[0], 0.0, &loc_data[0], (blasint)B->size[0]);
    }
    /* 'gkmPWMlasso4:340' fprintf('Running LASSO (2)\n'); */
    printf("Running LASSO (2)\n");
    fflush(stdout);
    /*  weigmat2 = lasso_cvmat(B, cfile3,'DFmax', d,'Standardize', false); */
    /* 'gkmPWMlasso4:342' [weigmat2, ~] = lasso_cvmat(B, cfile3, d, false, 100);
     */
    emxInit_real_T(&weigmat2, 2);
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
      b_comb_data = comb->data;
      loop_ub = B->size[0] * B->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        b_comb_data[i] = B_data[i];
      }
      b_lasso_cvmat(comb, loc, varargin_10, weigmat2, d_expl_temp_data,
                    comb_size, e_expl_temp_data, comb2_size, &GCneg1,
                    a__3_DF_data, a__3_DF_size, f_expl_temp_data,
                    expl_temp_size);
      AA_data = weigmat2->data;
      loc_data = loc->data;
    } else {
      c_binary_expand_op(B, cfile2, loc, varargin_10, weigmat2,
                         d_expl_temp_data, comb_size, e_expl_temp_data,
                         comb2_size, a__3_DF_data, a__3_DF_size,
                         f_expl_temp_data, expl_temp_size);
      AA_data = weigmat2->data;
    }
    /* 'gkmPWMlasso4:344' weigmat = abs(weigmat(:,1)) + abs(weigmat2(:,1)); */
    /* 'gkmPWMlasso4:345' f = find(weigmat~=0); */
    nx = weigmat->size[0] - 1;
    i = negvec->size[0];
    negvec->size[0] = weigmat->size[0];
    emxEnsureCapacity_real_T(negvec, i);
    negvec_data = negvec->data;
    for (k = 0; k <= nx; k++) {
      negvec_data[k] = fabs(xc_data[k]);
    }
    emxFree_real_T(&weigmat);
    nx = weigmat2->size[0] - 1;
    i = BY->size[0];
    BY->size[0] = weigmat2->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    for (k = 0; k <= nx; k++) {
      BY_data[k] = fabs(AA_data[k]);
    }
    emxFree_real_T(&weigmat2);
    if (negvec->size[0] == BY->size[0]) {
      i = b_loc->size[0];
      b_loc->size[0] = negvec->size[0];
      emxEnsureCapacity_boolean_T(b_loc, i);
      b_loc_data = b_loc->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_loc_data[i] = (negvec_data[i] + BY_data[i] != 0.0);
      }
      b_eml_find(b_loc, idx);
      matches_data = idx->data;
    } else {
      b_binary_expand_op(idx, negvec, BY);
      matches_data = idx->data;
    }
    i = f->size[0];
    f->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso4:346' motclus2 = clus_simmat_eig(corrcoef(B(:,f)),0.6); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = f->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    b_corrcoef(b_GCmat, comb);
    b_clus_simmat_eig(comb, tmp_motclus);
    tmp_motclus_data = tmp_motclus->data;
    /* 'gkmPWMlasso4:347' if length(length(motclus2)) ~= length(f) */
    if (1 != f->size[0]) {
      /* 'gkmPWMlasso4:348' f2 = f; */
      i = loc->size[0];
      loc->size[0] = f->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = f_data[i];
      }
      /* 'gkmPWMlasso4:349' f = zeros(length(motclus2),1); */
      i = f->size[0];
      f->size[0] = tmp_motclus->size[0];
      emxEnsureCapacity_real_T(f, i);
      f_data = f->data;
      loop_ub = tmp_motclus->size[0];
      for (i = 0; i < loop_ub; i++) {
        f_data[i] = 0.0;
      }
      /* 'gkmPWMlasso4:350' for i = 1:length(motclus2) */
      i = tmp_motclus->size[0];
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso4:351' if length(motclus2{i}) > 1 */
        if (tmp_motclus_data[b_i].f1->size[1] > 1) {
          /* 'gkmPWMlasso4:352' [~,f3] = max(abs(Z(f2(motclus2{i})))); */
          nx = tmp_motclus_data[b_i].f1->size[1];
          i1 = BY->size[0];
          BY->size[0] = tmp_motclus_data[b_i].f1->size[1];
          emxEnsureCapacity_real_T(BY, i1);
          BY_data = BY->data;
          for (k = 0; k < nx; k++) {
            BY_data[k] = fabs(
                Z_data[(int)loc_data[(int)tmp_motclus_data[b_i].f1->data[k] -
                                     1] -
                       1]);
          }
          d_maximum(BY, &GCneg1, &nx);
          /* 'gkmPWMlasso4:353' f(i) = f2(motclus2{i}(f3)); */
          f_data[b_i] =
              loc_data[(int)tmp_motclus_data[b_i].f1->data[nx - 1] - 1];
        } else {
          /* 'gkmPWMlasso4:354' else */
          /* 'gkmPWMlasso4:355' f(i) = f2(motclus2{i}(1)); */
          f_data[b_i] = loc_data[(int)tmp_motclus_data[b_i].f1->data[0] - 1];
        }
      }
    }
    /* 'gkmPWMlasso4:359' [a,b] = sort(abs(Z(f)), 'descend'); */
    nx = f->size[0];
    i = negvec->size[0];
    negvec->size[0] = f->size[0];
    emxEnsureCapacity_real_T(negvec, i);
    negvec_data = negvec->data;
    for (k = 0; k < nx; k++) {
      negvec_data[k] = fabs(Z_data[(int)f_data[k] - 1]);
    }
    sort(negvec, idx);
    matches_data = idx->data;
    /* 'gkmPWMlasso4:360' f = f(b); */
    i = BY->size[0];
    BY->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      BY_data[i] = f_data[matches_data[i] - 1];
    }
    i = f->size[0];
    f->size[0] = BY->size[0];
    emxEnsureCapacity_real_T(f, i);
    f_data = f->data;
    loop_ub = BY->size[0];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = BY_data[i];
    }
    /* 'gkmPWMlasso4:361' F = length(f); */
    GCneg1 = f->size[0];
    /* 'gkmPWMlasso4:362' fprintf('Selecting top motifs\n') */
    printf("Selecting top motifs\n");
    fflush(stdout);
    /* 'gkmPWMlasso4:363' ind = true; */
    /* 'gkmPWMlasso4:364' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = f->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, comb);
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = negvec->size[0];
      negvec->size[0] = f->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = 0.0;
      }
    } else {
      i = negvec->size[0];
      negvec->size[0] = f->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)f->size[0],
                  (blasint)1, (blasint)B->size[0], 1.0, &b_comb_data[0],
                  (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &negvec_data[0],
                  (blasint)f->size[0]);
    }
    mpower(comb, comb2);
    b_comb_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (negvec->size[0] == 0)) {
      i = diffc->size[0];
      diffc->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(diffc, i);
      b_diffc_data = diffc->data;
      loop_ub = comb2->size[0];
      for (i = 0; i < loop_ub; i++) {
        b_diffc_data[i] = 0.0;
      }
    } else {
      i = diffc->size[0];
      diffc->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(diffc, i);
      b_diffc_data = diffc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &b_comb_data[0], (blasint)comb2->size[0],
                  &negvec_data[0], (blasint)negvec->size[0], 0.0,
                  &b_diffc_data[0], (blasint)comb2->size[0]);
    }
    /* 'gkmPWMlasso4:365' Pweig = Z(f); */
    i = indc->size[0];
    indc->size[0] = f->size[0];
    emxEnsureCapacity_real_T(indc, i);
    b_indc_data = indc->data;
    loop_ub = f->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_indc_data[i] = Z_data[(int)f_data[i] - 1];
    }
    /* 'gkmPWMlasso4:366' while ind && F >= d */
    exitg1 = false;
    while ((!exitg1) && (GCneg1 >= varargin_10)) {
      /* 'gkmPWMlasso4:367' OLS =
       * (B(:,f(1:d)).'*B(:,f(1:d)))^-1*(B(:,f(1:d)).'*cfile2); */
      if (1.0 > varargin_10) {
        loop_ub = 0;
        nx = 0;
        lk = 0;
      } else {
        loop_ub = (int)varargin_10;
        nx = (int)varargin_10;
        lk = (int)varargin_10;
      }
      indc_size = B->size[0];
      i = b_GCmat->size[0] * b_GCmat->size[1];
      b_GCmat->size[0] = B->size[0];
      b_GCmat->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_GCmat, i);
      normvec_data = b_GCmat->data;
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < indc_size; i1++) {
          normvec_data[i1 + b_GCmat->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
        }
      }
      loop_ub = B->size[0];
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = B->size[0];
      comb2->size[1] = nx;
      emxEnsureCapacity_real_T(comb2, i);
      b_comb_data = comb2->data;
      for (i = 0; i < nx; i++) {
        for (i1 = 0; i1 < loop_ub; i1++) {
          b_comb_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
        }
      }
      b_mtimes(b_GCmat, comb2, comb);
      loop_ub = B->size[0];
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = B->size[0];
      comb2->size[1] = lk;
      emxEnsureCapacity_real_T(comb2, i);
      b_comb_data = comb2->data;
      for (i = 0; i < lk; i++) {
        for (i1 = 0; i1 < loop_ub; i1++) {
          b_comb_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
        }
      }
      if ((B->size[0] == 0) || (lk == 0) || (cfile2->size[0] == 0)) {
        i = negvec->size[0];
        negvec->size[0] = lk;
        emxEnsureCapacity_real_T(negvec, i);
        negvec_data = negvec->data;
        for (i = 0; i < lk; i++) {
          negvec_data[i] = 0.0;
        }
      } else {
        i = negvec->size[0];
        negvec->size[0] = lk;
        emxEnsureCapacity_real_T(negvec, i);
        negvec_data = negvec->data;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)lk,
                    (blasint)1, (blasint)B->size[0], 1.0, &b_comb_data[0],
                    (blasint)B->size[0], &cfile2_data[0],
                    (blasint)cfile2->size[0], 0.0, &negvec_data[0],
                    (blasint)lk);
      }
      mpower(comb, comb2);
      b_comb_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (negvec->size[0] == 0)) {
        i = diffc->size[0];
        diffc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(diffc, i);
        b_diffc_data = diffc->data;
        loop_ub = comb2->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_diffc_data[i] = 0.0;
        }
      } else {
        i = diffc->size[0];
        diffc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(diffc, i);
        b_diffc_data = diffc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &b_comb_data[0],
            (blasint)comb2->size[0], &negvec_data[0], (blasint)negvec->size[0],
            0.0, &b_diffc_data[0], (blasint)comb2->size[0]);
      }
      /* 'gkmPWMlasso4:368' Pweig = Z(f(1:d)); */
      if (1.0 > varargin_10) {
        loop_ub = 0;
      } else {
        loop_ub = (int)varargin_10;
      }
      i = indc->size[0];
      indc->size[0] = loop_ub;
      emxEnsureCapacity_real_T(indc, i);
      b_indc_data = indc->data;
      for (i = 0; i < loop_ub; i++) {
        b_indc_data[i] = Z_data[(int)f_data[i] - 1];
      }
      /* 'gkmPWMlasso4:369' ff = []; */
      loc->size[0] = 0;
      /* 'gkmPWMlasso4:370' for i = 1:d */
      i = (int)varargin_10;
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso4:371' if sign(OLS(i)) ~= sign(Pweig(i)) */
        b1 = b_diffc_data[b_i];
        if (b_diffc_data[b_i] < 0.0) {
          b1 = -1.0;
        } else if (b_diffc_data[b_i] > 0.0) {
          b1 = 1.0;
        }
        GCpos1 = Z_data[(int)f_data[b_i] - 1];
        if (GCpos1 < 0.0) {
          GCpos1 = -1.0;
        } else if (GCpos1 > 0.0) {
          GCpos1 = 1.0;
        }
        if (b1 != GCpos1) {
          /* 'gkmPWMlasso4:372' ff = [ff;i]; */
          i1 = loc->size[0];
          i2 = loc->size[0];
          loc->size[0]++;
          emxEnsureCapacity_real_T(loc, i2);
          loc_data = loc->data;
          loc_data[i1] = (double)b_i + 1.0;
        }
      }
      /* 'gkmPWMlasso4:375' if length(ff) > 0 && F - length(ff) >= d */
      if ((loc->size[0] > 0) &&
          (GCneg1 - (double)loc->size[0] >= varargin_10)) {
        /* 'gkmPWMlasso4:376' f(ff) = []; */
        i = idx->size[0];
        idx->size[0] = loc->size[0];
        emxEnsureCapacity_int32_T(idx, i);
        matches_data = idx->data;
        loop_ub = loc->size[0];
        for (i = 0; i < loop_ub; i++) {
          matches_data[i] = (int)loc_data[i];
        }
        b_nullAssignment(f, idx);
        f_data = f->data;
        /* 'gkmPWMlasso4:377' F = F - length(ff); */
        GCneg1 -= (double)loc->size[0];
      } else {
        /* 'gkmPWMlasso4:378' else */
        /* 'gkmPWMlasso4:379' ind = false; */
        /* 'gkmPWMlasso4:380' f = f(1:d); */
        i = f->size[0];
        if (1.0 > varargin_10) {
          f->size[0] = 0;
        } else {
          f->size[0] = (int)varargin_10;
        }
        emxEnsureCapacity_real_T(f, i);
        f_data = f->data;
        exitg1 = true;
      }
    }
    /*  if F < d */
    /*      OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    /*      Pweig = Z(f); */
    /*  end */
    /* 'gkmPWMlasso4:387' BB = B(:,f); */
    loop_ub = B->size[0];
    i = BB->size[0] * BB->size[1];
    BB->size[0] = B->size[0];
    BB->size[1] = f->size[0];
    emxEnsureCapacity_real_T(BB, i);
    AA_data = BB->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        AA_data[i1 + BB->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    /* 'gkmPWMlasso4:388' BX = B(:,f)'*B(:,f); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = f->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, xc);
    xc_data = xc->data;
    /* 'gkmPWMlasso4:389' BY = B(:,f)'*cfile2; */
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = f->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    b_comb_data = comb2->data;
    nx = f->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_comb_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)f_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (f->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = f->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = f->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)f->size[0],
                  (blasint)1, (blasint)B->size[0], 1.0, &b_comb_data[0],
                  (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)f->size[0]);
    }
    /* 'gkmPWMlasso4:390' E = zeros(length(f),1); */
    i = normvec->size[0];
    normvec->size[0] = f->size[0];
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    loop_ub = f->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso4:391' for i = 1:length(f) */
    i = f->size[0];
    if (0 <= f->size[0] - 1) {
      d_loop_ub = xc->size[0] * xc->size[1];
      e_loop_ub = BY->size[0];
      nxout = BY->size[0] - 1;
      f_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:392' B = BB; */
      /* 'gkmPWMlasso4:393' BBX = BX; */
      /* 'gkmPWMlasso4:394' BBY = BY; */
      /* 'gkmPWMlasso4:395' B(:,i) = []; */
      /* 'gkmPWMlasso4:396' BBX(:,i) = []; */
      /* 'gkmPWMlasso4:397' BBX(i,:) = []; */
      /* 'gkmPWMlasso4:398' BBY(i) = []; */
      /* 'gkmPWMlasso4:399' res = cfile2-B*(BBX^-1*BBY); */
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = xc->size[0];
      comb->size[1] = xc->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      b_comb_data = comb->data;
      for (i1 = 0; i1 < d_loop_ub; i1++) {
        b_comb_data[i1] = xc_data[i1];
      }
      c_nullAssignment(comb, b_i + 1);
      nullAssignment(comb, b_i + 1);
      mpower(comb, comb2);
      b_comb_data = comb2->data;
      nx = b_i + 1;
      i1 = loc->size[0];
      loc->size[0] = BY->size[0];
      emxEnsureCapacity_real_T(loc, i1);
      loc_data = loc->data;
      for (i1 = 0; i1 < e_loop_ub; i1++) {
        loc_data[i1] = BY_data[i1];
      }
      for (k = nx; k <= nxout; k++) {
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
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (loc->size[0] == 0)) {
        i1 = negvec->size[0];
        negvec->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(negvec, i1);
        negvec_data = negvec->data;
        loop_ub = comb2->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          negvec_data[i1] = 0.0;
        }
      } else {
        i1 = negvec->size[0];
        negvec->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(negvec, i1);
        negvec_data = negvec->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &b_comb_data[0],
            (blasint)comb2->size[0], &loc_data[0], (blasint)loc->size[0], 0.0,
            &negvec_data[0], (blasint)comb2->size[0]);
      }
      i1 = comb2->size[0] * comb2->size[1];
      comb2->size[0] = BB->size[0];
      comb2->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb2, i1);
      b_comb_data = comb2->data;
      for (i1 = 0; i1 < f_loop_ub; i1++) {
        b_comb_data[i1] = AA_data[i1];
      }
      c_nullAssignment(comb2, b_i + 1);
      b_comb_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (negvec->size[0] == 0)) {
        i1 = loc->size[0];
        loc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        loop_ub = comb2->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          loc_data[i1] = 0.0;
        }
      } else {
        i1 = loc->size[0];
        loc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(loc, i1);
        loc_data = loc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &b_comb_data[0],
            (blasint)comb2->size[0], &negvec_data[0], (blasint)negvec->size[0],
            0.0, &loc_data[0], (blasint)comb2->size[0]);
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
      /* 'gkmPWMlasso4:400' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        b1 = 0.0;
      } else {
        b1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                        &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(b1);
    }
    /* 'gkmPWMlasso4:402' res = cfile2-BB*OLS; */
    if ((B->size[0] == 0) || (f->size[0] == 0) || (diffc->size[0] == 0)) {
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
                  (blasint)B->size[0], (blasint)1, (blasint)f->size[0], 1.0,
                  &AA_data[0], (blasint)B->size[0], &b_diffc_data[0],
                  (blasint)diffc->size[0], 0.0, &loc_data[0],
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
    /* 'gkmPWMlasso4:403' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      b1 = 0.0;
    } else {
      b1 = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                      &loc_data[0], (blasint)1);
    }
    GCneg1 = sqrt(b1);
    /* 'gkmPWMlasso4:404' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - GCneg1) / GCneg1;
    }
    /* 'gkmPWMlasso4:405' for i = 1:length(motclus) */
    i = motclus->size[0];
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:406' motclus{i} = indvec(motclus{i})'; */
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
    /* 'gkmPWMlasso4:409' mylen = length(f); */
    /* 'gkmPWMlasso4:410' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso4:411' for idx=1:mylen */
    i = f->size[0];
    i1 = c_motclus->size[0];
    c_motclus->size[0] = f->size[0];
    emxEnsureCapacity_cell_wrap_1(c_motclus, i1);
    b_motclus_data = c_motclus->data;
    for (nx = 0; nx < i; nx++) {
      /* 'gkmPWMlasso4:412' newMotclus{idx} = motclus{f(idx)}; */
      i1 = b_motclus_data[nx].f1->size[0] * b_motclus_data[nx].f1->size[1];
      b_motclus_data[nx].f1->size[0] = 1;
      b_motclus_data[nx].f1->size[1] =
          motclus_data[(int)f_data[nx] - 1].f1->size[1];
      emxEnsureCapacity_real_T(b_motclus_data[nx].f1, i1);
      loop_ub = motclus_data[(int)f_data[nx] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_motclus_data[nx].f1->data[i1] =
            motclus_data[(int)f_data[nx] - 1].f1->data[i1];
      }
    }
    /* 'gkmPWMlasso4:414' motclus = newMotclus; */
    /* 'gkmPWMlasso4:416' correlation = corrcoef(cfile2, BB*OLS); */
    if ((B->size[0] == 0) || (f->size[0] == 0) || (diffc->size[0] == 0)) {
      loop_ub = B->size[0];
      i = negvec->size[0];
      negvec->size[0] = B->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = 0.0;
      }
    } else {
      i = negvec->size[0];
      negvec->size[0] = B->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)B->size[0], (blasint)1, (blasint)f->size[0], 1.0,
                  &AA_data[0], (blasint)B->size[0], &b_diffc_data[0],
                  (blasint)diffc->size[0], 0.0, &negvec_data[0],
                  (blasint)B->size[0]);
    }
    /* 'gkmPWMlasso4:417' gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus,
     * sprintf("%s_%d_%d", filename, int32(l_svm2),
     * int32(k_svm2)),memefile,num,minL, minInfo, correlation(1,2)); */
    GCneg1 = maximum(diffc);
    GCpos1 = maximum(normvec);
    b_sprintf(varargin_1, (int)rt_roundd(varargin_6),
              (int)rt_roundd(varargin_7), text);
    loop_ub = diffc->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_diffc_data[i] /= GCneg1;
    }
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] /= GCpos1;
    }
    corrcoef(cfile2, negvec, dv);
    gettopmotifs(diffc, indc, normvec, c_motclus, text, varargin_2,
                 size_tmp_idx_1, varargin_3, varargin_4, dv[2]);
  }
  emxFree_real_T(&b_GCmat);
  emxFree_boolean_T(&b_loc);
  emxFree_int32_T(&idx);
  emxFree_int32_T(&match_out);
  emxFree_char_T(&text);
  emxFree_cell_wrap_1(&c_motclus);
  emxFree_real_T(&loc);
  emxFree_real_T(&BY);
  emxFree_real_T(&BB);
  emxFree_cell_wrap_1(&tmp_motclus);
  emxFree_real_T(&f);
  emxFree_real_T(&Z);
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
  /* 'gkmPWMlasso4:421' fprintf('Done\n'); */
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

/* End of code generation (gkmPWMlasso4.c) */
