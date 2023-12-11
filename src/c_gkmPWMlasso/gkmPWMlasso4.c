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
  /* 'gkmPWMlasso4:481' bin = conncomp(graph(simmat>r)); */
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
  /* 'gkmPWMlasso4:482' n = max(bin); */
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
  /* 'gkmPWMlasso4:483' motclus = cell(n,1); */
  /* 'gkmPWMlasso4:484' for i = 1:n */
  k = (int)n;
  i = motclus->size[0];
  motclus->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  emxInit_int32_T(&b_r, 2);
  emxInit_boolean_T(&b_bin, 2);
  for (b_i = 0; b_i < k; b_i++) {
    /* 'gkmPWMlasso4:485' motclus{i} = find(bin == i); */
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
  /* 'gkmPWMlasso4:481' bin = conncomp(graph(simmat>r)); */
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
  /* 'gkmPWMlasso4:482' n = max(bin); */
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
  /* 'gkmPWMlasso4:483' motclus = cell(n,1); */
  /* 'gkmPWMlasso4:484' for i = 1:n */
  k = (int)n;
  i = motclus->size[0];
  motclus->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  emxInit_int32_T(&c_r, 2);
  emxInit_boolean_T(&b_bin, 2);
  for (b_i = 0; b_i < k; b_i++) {
    /* 'gkmPWMlasso4:485' motclus{i} = find(bin == i); */
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
  /* 'gkmPWMlasso4:437' [a, b] = sort(pvec, 'descend'); */
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
  /* 'gkmPWMlasso4:438' c = [weigvec(b) a E(b)]; */
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
  /* 'gkmPWMlasso4:440' mylen = length(b); */
  /* 'gkmPWMlasso4:441' ordered_motclus = cell(mylen, 1); */
  /* 'gkmPWMlasso4:442' for idx=1:mylen */
  i = iidx->size[0];
  i1 = ordered_motclus->size[0];
  ordered_motclus->size[0] = iidx->size[0];
  emxEnsureCapacity_cell_wrap_1(ordered_motclus, i1);
  ordered_motclus_data = ordered_motclus->data;
  for (k = 0; k < i; k++) {
    /* 'gkmPWMlasso4:443' ordered_motclus{idx} = motclus{b(idx)}; */
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
  /* 'gkmPWMlasso4:445' motclus = ordered_motclus; */
  /* 'gkmPWMlasso4:446' [rr,cc] = size(c); */
  /* 'gkmPWMlasso4:447' names = cell(n,1); */
  nbytes = (int)n;
  i = names->size[0];
  names->size[0] = (int)n;
  emxEnsureCapacity_cell_wrap_8(names, i);
  names_data = names->data;
  for (i = 0; i < nbytes; i++) {
    names_data[i].f1->size[0] = 1;
    names_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso4:448' names = coder.nullcopy(names); */
  /* 'gkmPWMlasso4:449' i = 0; */
  b_i = 0.0;
  /* 'gkmPWMlasso4:450' fid = fopen(memefile, 'r'); */
  fileid = cfopen(memefile, "rb");
  /* 'gkmPWMlasso4:451' if fid == -1 */
  if (fileid == -1) {
    /* 'gkmPWMlasso4:452' fprintf("ERROR: Cannot open motif database.\n"); */
    printf("ERROR: Cannot open motif database.\n");
    fflush(stdout);
  }
  /* 'gkmPWMlasso4:454' while ~feof(fid) */
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
      /* 'gkmPWMlasso4:455' line = fgetl(fid); */
      fgetl(fileid, line);
      line_data = line->data;
      /* 'gkmPWMlasso4:456' if length(line) >= 5 */
      if (line->size[1] >= 5) {
        /* 'gkmPWMlasso4:457' if strcmp(line(1:5), 'MOTIF') */
        for (i = 0; i < 5; i++) {
          c_a[i] = line_data[i];
        }
        nbytes = memcmp(&c_a[0], &b[0], 5);
        if (nbytes == 0) {
          /* 'gkmPWMlasso4:458' i = i+1; */
          b_i++;
          /* 'gkmPWMlasso4:459' [~,name] = strtok(line); */
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
          /* 'gkmPWMlasso4:460' names{i} = strtrim(name); */
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
  /* 'gkmPWMlasso4:464' fclose(fid); */
  cfclose(fileid);
  /* 'gkmPWMlasso4:465' fidw = fopen(sprintf('%s_gkmPWMlasso4.out', filename),
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
  /* 'gkmPWMlasso4:466' if fidw == -1 */
  emxFree_char_T(&varargin_1);
  emxFree_char_T(&charStr);
  if (fileid == -1) {
    /* 'gkmPWMlasso4:467' fprintf("ERROR: Cannot create gkmPWMlasso output
     * file.\n"); */
    printf("ERROR: Cannot create gkmPWMlasso output file.\n");
    fflush(stdout);
  }
  /* 'gkmPWMlasso4:469' fprintf(fidw, 'Minimum PWM Length:\t%d\n', int32(minL));
   */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Minimum PWM Length:\t%d\n", (int)rt_roundd(minL));
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso4:470' fprintf(fidw, 'Minimum PWM Information:\t%f\n',
   * minInfo); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Minimum PWM Information:\t%f\n", minInfo);
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso4:471' fprintf(fidw, 'Correlation with SVM weights:\t%f\n', C);
   */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &b_a);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "Correlation with SVM weights:\t%f\n", C);
    if (b_a) {
      fflush(filestar);
    }
  }
  /* 'gkmPWMlasso4:472' fprintf(fidw, 'Cluster ID\tMotif ID\tMotif
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
  /* 'gkmPWMlasso4:473' for j = 1:length(E) */
  i = E->size[0];
  for (k = 0; k < i; k++) {
    /* 'gkmPWMlasso4:474' for l = 1:length(motclus{j}) */
    i1 = ordered_motclus_data[k].f1->size[1];
    if (0 <= ordered_motclus_data[k].f1->size[1] - 1) {
      c_NULL = NULL;
    }
    for (l = 0; l < i1; l++) {
      /* 'gkmPWMlasso4:475' fprintf(fidw, '%d\t%d\t%s\t%0.3f\t%0.3f\t%0.3f\n',
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
  /* 'gkmPWMlasso4:478' fclose(fidw); */
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
  /* 'gkmPWMlasso4:489' l = length(p); */
  /* 'gkmPWMlasso4:490' info = zeros(l, 1); */
  nx = p->size[0];
  i = info->size[0];
  info->size[0] = nx;
  emxEnsureCapacity_real_T(info, i);
  info_data = info->data;
  /* 'gkmPWMlasso4:491' len = zeros(l,1); */
  i = len->size[0];
  len->size[0] = nx;
  emxEnsureCapacity_real_T(len, i);
  len_data = len->data;
  /* 'gkmPWMlasso4:492' for i = 1:l */
  i = p->size[0];
  emxInit_real_T(&mat, 2);
  emxInit_real_T(&vec, 1);
  emxInit_real_T(&mvec, 1);
  emxInit_real_T(&b_r, 2);
  emxInit_real_T(&b_mat, 2);
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWMlasso4:493' mat = p{i}+(p{i}==0); */
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
    /* 'gkmPWMlasso4:494' vec = 2+sum(mat.*log(mat)/log(2),2); */
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
    /* 'gkmPWMlasso4:495' mvec = movmean(vec,5); */
    nx = vec->size[0];
    if (2 <= nx) {
      nx = 2;
    }
    vmovfun(vec, vec->size[0], nx, nx, mvec);
    mvec_data = mvec->data;
    /* 'gkmPWMlasso4:496' while min(vec(1),mvec(1)) < cut && length(vec) > 1 */
    while ((fmin(mat_data[0], mvec_data[0]) < 0.0) && (vec->size[0] > 1)) {
      /* 'gkmPWMlasso4:497' p{i}(1,:) = []; */
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
      /* 'gkmPWMlasso4:498' vec(1) = []; */
      nx = vec->size[0];
      nxout = vec->size[0];
      for (nrows = 0; nrows <= nxout - 2; nrows++) {
        mat_data[nrows] = mat_data[nrows + 1];
      }
      idx = vec->size[0];
      vec->size[0] = nx - 1;
      emxEnsureCapacity_real_T(vec, idx);
      mat_data = vec->data;
      /* 'gkmPWMlasso4:499' mvec(1)=[]; */
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
      /* 'gkmPWMlasso4:500' mat(1,:) = []; */
      idx = mat->size[0] * mat->size[1];
      if (1 > mat->size[0] - 1) {
        mat->size[0] = 0;
      } else {
        mat->size[0]--;
      }
      emxEnsureCapacity_real_T(mat, idx);
    }
    /* 'gkmPWMlasso4:502' while min(mvec(end),vec(end)) < cut && length(vec) > 1
     */
    while ((fmin(mvec_data[mvec->size[0] - 1], mat_data[vec->size[0] - 1]) <
            0.0) &&
           (vec->size[0] > 1)) {
      /* 'gkmPWMlasso4:503' vec(end) = []; */
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
      /* 'gkmPWMlasso4:504' mvec(end)=[]; */
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
      /* 'gkmPWMlasso4:505' mat(end,:) = []; */
      idx = mat->size[0] * mat->size[1];
      if (1 > mat->size[0] - 1) {
        mat->size[0] = 0;
      } else {
        mat->size[0]--;
      }
      emxEnsureCapacity_real_T(mat, idx);
      /* 'gkmPWMlasso4:506' p{i}(end,:) = []; */
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
    /* 'gkmPWMlasso4:508' info(i) = sum(vec); */
    info_data[b_i] = blockedSummation(vec, vec->size[0]);
    /* 'gkmPWMlasso4:509' [len(i), ~] = size(mat); */
    len_data[b_i] = mat->size[0];
  }
  emxFree_real_T(&b_mat);
  emxFree_real_T(&b_r);
  emxFree_real_T(&mvec);
  emxFree_real_T(&vec);
  emxFree_real_T(&mat);
  /* 'gkmPWMlasso4:511' pp = p; */
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
                  double varargin_10, double varargin_11)
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
  emxArray_real_T *AA;
  emxArray_real_T *B;
  emxArray_real_T *BB;
  emxArray_real_T *BY;
  emxArray_real_T *GCmat;
  emxArray_real_T *OLS;
  emxArray_real_T *Pweig;
  emxArray_real_T *Z;
  emxArray_real_T *b_A;
  emxArray_real_T *b_GCmat;
  emxArray_real_T *b_f;
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
  double c;
  double nfrac;
  double rcnum;
  double *A_data;
  double *BY_data;
  double *B_data;
  double *GCmat_data;
  double *OLS_data;
  double *Pweig_data;
  double *Z_data;
  double *cfile2_data;
  double *f_data;
  double *indc_data;
  double *indvec_data;
  double *loc_data;
  double *negvec_data;
  double *normvec_data;
  double *xc_data;
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
  int i7;
  int i_loop_ub;
  int j;
  int j_loop_ub;
  int k;
  int k_loop_ub;
  int l_loop_ub;
  int loop_ub;
  int m_loop_ub;
  int match_idx;
  int n_loop_ub;
  int nx;
  int nxout;
  int o_loop_ub;
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
  /*  if nargin < 2 */
  /*      error('Need at least 3 inputs') */
  /*  end */
  /* 'gkmPWMlasso4:7' filename = varargin{1}; */
  /* 'gkmPWMlasso4:8' memefile = varargin{2}; */
  /* 'gkmPWMlasso4:9' minL = varargin{3}; */
  /* 'gkmPWMlasso4:10' minInfo = varargin{4}; */
  /* 'gkmPWMlasso4:11' corrCut = varargin{5}; */
  /* 'gkmPWMlasso4:12' l_svm = varargin{6}; */
  /* 'gkmPWMlasso4:13' k_svm = varargin{7}; */
  /* 'gkmPWMlasso4:14' BG_GC = varargin{8}; */
  /* 'gkmPWMlasso4:15' RC = varargin{9}; */
  /* 'gkmPWMlasso4:16' d = varargin{10}; */
  /* 'gkmPWMlasso4:17' nfrac = varargin{11}; */
  /* 'gkmPWMlasso4:18' lk = 1; */
  lk_size[0] = 1;
  lk_size[1] = 1;
  lk_data[0] = 1.0;
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
  /* 'gkmPWMlasso4:63' if nfrac ~= 1 */
  if (varargin_11 != 1.0) {
    /* 'gkmPWMlasso4:64' lk = [l_svm k_svm]; */
    lk_size[0] = 1;
    lk_size[1] = 2;
    lk_data[0] = varargin_6;
    lk_data[1] = varargin_7;
  }
  emxInit_real_T(&comb, 2);
  emxInit_real_T(&comb2, 2);
  emxInit_real_T(&diffc, 1);
  emxInit_real_T(&indc, 1);
  emxInit_real_T(&xc, 2);
  /* 'gkmPWMlasso4:67' [comb,comb2,diffc,indc,xc,rcnum] =
   * genIndex(l_svm,k_svm,nfrac); */
  genIndex(varargin_6, varargin_7, varargin_11, comb, comb2, diffc, indc, xc,
           &rcnum);
  /* generate gapped positions, adjusted for reverse complements */
  /* 'gkmPWMlasso4:69' if length(comb)*4^k_svm > 6*10^5 */
  if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
    u1 = 0;
  } else {
    nx = comb->size[0];
    u1 = comb->size[1];
    if (nx >= u1) {
      u1 = nx;
    }
  }
  c = pow(4.0, varargin_7);
  if ((double)u1 * c > 600000.0) {
    /* 'gkmPWMlasso4:70' nfrac = round(5*10^7/4^k_svm/length(comb))/100; */
    if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
      u1 = 0;
    } else {
      nx = comb->size[0];
      u1 = comb->size[1];
      if (nx >= u1) {
        u1 = nx;
      }
    }
    nfrac = rt_roundd(5.0E+7 / c / (double)u1) / 100.0;
    /* 'gkmPWMlasso4:71' fprintf('WARNING: Combination of (l,k) yields too many
     * gapped kmers.  Using %f of the total gapped kmers', nfrac); */
    printf("WARNING: Combination of (l,k) yields too many gapped kmers.  Using "
           "%f of the total gapped kmers",
           nfrac);
    fflush(stdout);
    /* 'gkmPWMlasso4:72' l_svm2 = l_svm; */
    /* 'gkmPWMlasso4:73' k_svm2 = k_svm; */
    /* 'gkmPWMlasso4:74' lk = ([l_svm k_svm]); */
    lk_size[0] = 1;
    lk_size[1] = 2;
    lk_data[0] = varargin_6;
    lk_data[1] = varargin_7;
    /* 'gkmPWMlasso4:75' [comb,comb2,diffc,indc,xc,rcnum] =
     * genIndex(l_svm,k_svm,nfrac); */
    genIndex(varargin_6, varargin_7, nfrac, comb, comb2, diffc, indc, xc,
             &rcnum);
  } else {
    /* 'gkmPWMlasso4:76' else */
    /* 'gkmPWMlasso4:77' l_svm2 = l_svm; */
    /* 'gkmPWMlasso4:78' k_svm2 = k_svm; */
  }
  emxInit_real_T(&cfile2, 1);
  /* 'gkmPWMlasso4:81' fprintf('Counting gapped kmers\n'); */
  printf("Counting gapped kmers\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:83' [cfile, GCpos1, GCneg1,mat,mat2] = getgkmcounts(filename,
   * l_svm, k_svm, lk, RC, comb,rcnum); */
  getgkmcounts(varargin_1, varargin_6, varargin_7, lk_data, lk_size, varargin_9,
               comb, rcnum, cfile2, &GCpos1, &nfrac, mat, mat2);
  cfile2_data = cfile2->data;
  /* 'gkmPWMlasso4:84' if BG_GC == 1 */
  if (varargin_8 == 1.0) {
    /* 'gkmPWMlasso4:85' mat = (mat+mat2)/2; */
    for (i = 0; i < 16; i++) {
      mat[i] = (mat[i] + mat2[i]) / 2.0;
    }
    /* 'gkmPWMlasso4:86' GCpos1 = (GCpos1+GCneg1)/2; */
    GCpos1 = (GCpos1 + nfrac) / 2.0;
    /* 'gkmPWMlasso4:87' GCneg1 = GCpos1; */
    nfrac = GCpos1;
  }
  emxInit_real_T(&negvec, 1);
  emxInit_char_T(&text, 2);
  /* 'gkmPWMlasso4:89' fprintf('Generating negative set\n'); */
  printf("Generating negative set\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:90' negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC);
   */
  BGkmer(mat, nfrac, comb, rcnum, varargin_6, varargin_7, varargin_9, negvec);
  negvec_data = negvec->data;
  /* 'gkmPWMlasso4:92' fprintf('Filtering motifs\n'); */
  printf("Filtering motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:93' num = length(strfind(fileread(memefile),'MOTIF')); */
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
  /* 'gkmPWMlasso4:94' p = getmotif(memefile,1:num); */
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
  /* 'gkmPWMlasso4:95' for i = 1:num */
  for (b_i = 0; b_i < size_tmp_idx_1; b_i++) {
    /* 'gkmPWMlasso4:96' [r c] = size(p{i}); */
    /* 'gkmPWMlasso4:97' for j = 1:r */
    i = p_data[b_i].f1->size[0];
    for (j = 0; j < i; j++) {
      /* 'gkmPWMlasso4:98' a = sum(p{i}(j,:)); */
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
      /* 'gkmPWMlasso4:99' if abs(a-1)>0 */
      if (fabs(a - 1.0) > 0.0) {
        /* 'gkmPWMlasso4:100' [b1 loc] = max(p{i}(j,:)); */
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
        /* 'gkmPWMlasso4:101' p{i}(j,loc) = b1-a+1; */
        p_data[b_i].f1->data[j + p_data[b_i].f1->size[0] * (match_idx - 1)] =
            (nfrac - a) + 1.0;
      }
    }
  }
  emxInit_real_T(&Z, 1);
  emxInit_real_T(&loc, 1);
  /* 'gkmPWMlasso4:105' [p, info, lenvec] = trim_pwm(p,0.0); */
  trim_pwm(p, loc, Z);
  Z_data = Z->data;
  loc_data = loc->data;
  p_data = p->data;
  /* 'gkmPWMlasso4:106' indvec =
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
  emxInit_real_T(&A, 1);
  i = A->size[0];
  A->size[0] = ib->size[0];
  emxEnsureCapacity_real_T(A, i);
  xc_data = A->data;
  loop_ub = ib->size[0];
  for (i = 0; i < loop_ub; i++) {
    xc_data[i] = match_out_data[i];
  }
  emxInit_real_T(&indvec, 1);
  do_vectors(BY, A, indvec, idx, ib);
  indvec_data = indvec->data;
  /* 'gkmPWMlasso4:107' n = length(indvec); */
  /* 'gkmPWMlasso4:108' fprintf('Mapping PWMs to gkm space\n'); */
  printf("Mapping PWMs to gkm space\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:109' lcnum = length(comb); */
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
  emxInit_real_T(&b_A, 2);
  /* 'gkmPWMlasso4:110' A=zeros(lcnum*4^k_svm,n); */
  i = (int)((double)u1 * c);
  i1 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = i;
  b_A->size[1] = indvec->size[0];
  emxEnsureCapacity_real_T(b_A, i1);
  A_data = b_A->data;
  loop_ub = i * indvec->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    A_data[i1] = 0.0;
  }
  emxInit_real_T(&AA, 2);
  /* 'gkmPWMlasso4:111' AA=zeros(lcnum*4^k_svm,n); */
  i1 = AA->size[0] * AA->size[1];
  AA->size[0] = i;
  AA->size[1] = indvec->size[0];
  emxEnsureCapacity_real_T(AA, i1);
  f_data = AA->data;
  loop_ub = i * indvec->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    f_data[i1] = 0.0;
  }
  emxInit_real_T(&GCmat, 2);
  emxInit_real_T(&normvec, 1);
  /* 'gkmPWMlasso4:112' GCmat = repmat([0.5-GCpos1/2 GCpos1/2 GCpos1/2
   * 0.5-GCpos1/2],l_svm-1,1); */
  c = 0.5 - GCpos1 / 2.0;
  dv[0] = c;
  dv[1] = GCpos1 / 2.0;
  dv[2] = GCpos1 / 2.0;
  dv[3] = c;
  b_repmat(dv, varargin_6 - 1.0, GCmat);
  GCmat_data = GCmat->data;
  /* 'gkmPWMlasso4:113' per = 10; */
  GCpos1 = 10.0;
  /* 'gkmPWMlasso4:114' normvec = zeros(n,1); */
  /* 'gkmPWMlasso4:115' for j = 1:n */
  i1 = indvec->size[0];
  emxInit_real_T(&b_GCmat, 2);
  for (j = 0; j < i1; j++) {
    /* 'gkmPWMlasso4:116' if mod(j, floor(n/10))==0 */
    if (b_mod((double)j + 1.0, floor((double)indvec->size[0] / 10.0)) == 0.0) {
      /* 'gkmPWMlasso4:117' fprintf('%d...', int32(per)); */
      printf("%d...", (int)GCpos1);
      fflush(stdout);
      /* 'gkmPWMlasso4:118' per = per+10; */
      GCpos1 += 10.0;
    }
    /* 'gkmPWMlasso4:120' loc = zeros(l_svm*2-2+lenvec(indvec(j)), 1); */
    i2 = (int)indvec_data[j] - 1;
    c = Z_data[i2];
    loop_ub = (int)((varargin_6 * 2.0 - 2.0) + c);
    b_i = loc->size[0];
    loc->size[0] = loop_ub;
    emxEnsureCapacity_real_T(loc, b_i);
    loc_data = loc->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      loc_data[b_i] = 0.0;
    }
    /* 'gkmPWMlasso4:121' loc(l_svm:lenvec(indvec(j))+l_svm-1) = 1; */
    nfrac = (c + varargin_6) - 1.0;
    if (varargin_6 > nfrac) {
      b_i = -1;
      i3 = 0;
    } else {
      b_i = (int)varargin_6 - 2;
      i3 = (int)nfrac;
    }
    loop_ub = (i3 - b_i) - 1;
    for (i3 = 0; i3 < loop_ub; i3++) {
      loc_data[(b_i + i3) + 1] = 1.0;
    }
    /* 'gkmPWMlasso4:122' if RC */
    if (varargin_9) {
      /* 'gkmPWMlasso4:123' A(:,j) =
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
      b_i = b_GCmat->size[0] * b_GCmat->size[1];
      b_GCmat->size[0] = (b_loop_ub + match_idx) + nx;
      b_GCmat->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_GCmat, b_i);
      normvec_data = b_GCmat->data;
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          normvec_data[i3 + b_GCmat->size[0] * b_i] =
              GCmat_data[i3 + b_loop_ub * b_i];
        }
      }
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < match_idx; i3++) {
          normvec_data[(i3 + b_loop_ub) + b_GCmat->size[0] * b_i] =
              p_data[i2].f1->data[i3 + match_idx * b_i];
        }
      }
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (b_i = 0; b_i < nx; b_i++) {
          normvec_data[((b_i + b_loop_ub) + match_idx) +
                       b_GCmat->size[0] * i2] = GCmat_data[b_i + nx * i2];
        }
      }
      PWM2kmers(b_GCmat, mat, comb2, diffc, indc, loc, xc, varargin_6,
                varargin_7, rcnum, A);
      xc_data = A->data;
      loop_ub = A->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        A_data[i2 + b_A->size[0] * j] = xc_data[i2];
      }
    } else {
      /* 'gkmPWMlasso4:124' else */
      /* 'gkmPWMlasso4:125' A(:,j) =
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
      b_i = b_GCmat->size[0] * b_GCmat->size[1];
      b_GCmat->size[0] = (b_loop_ub + match_idx) + nx;
      b_GCmat->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_GCmat, b_i);
      normvec_data = b_GCmat->data;
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          normvec_data[i3 + b_GCmat->size[0] * b_i] =
              GCmat_data[i3 + b_loop_ub * b_i];
        }
      }
      for (b_i = 0; b_i < loop_ub; b_i++) {
        for (i3 = 0; i3 < match_idx; i3++) {
          normvec_data[(i3 + b_loop_ub) + b_GCmat->size[0] * b_i] =
              p_data[i2].f1->data[i3 + match_idx * b_i];
        }
      }
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (b_i = 0; b_i < nx; b_i++) {
          normvec_data[((b_i + b_loop_ub) + match_idx) +
                       b_GCmat->size[0] * i2] = GCmat_data[b_i + nx * i2];
        }
      }
      PWM2kmers_norc(b_GCmat, mat, comb2, diffc, indc, loc, xc, varargin_6,
                     varargin_7, A);
      xc_data = A->data;
      loop_ub = A->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        A_data[i2 + b_A->size[0] * j] = xc_data[i2];
      }
    }
    /* 'gkmPWMlasso4:127' A(:,j) = A(:,j) - negvec*(l_svm-1+lenvec(indvec(j)));
     */
    c += varargin_6 - 1.0;
    if (b_A->size[0] == negvec->size[0]) {
      match_idx = b_A->size[0] - 1;
      i2 = A->size[0];
      A->size[0] = b_A->size[0];
      emxEnsureCapacity_real_T(A, i2);
      xc_data = A->data;
      for (i2 = 0; i2 <= match_idx; i2++) {
        xc_data[i2] = A_data[i2 + b_A->size[0] * j] - negvec_data[i2] * c;
      }
      loop_ub = A->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        A_data[i2 + b_A->size[0] * j] = xc_data[i2];
      }
    } else {
      binary_expand_op(b_A, j, negvec, c);
      A_data = b_A->data;
    }
    /* 'gkmPWMlasso4:128' normvec(j) = (A(:,j)'*A(:,j))^0.5; */
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
      c = A_data[i2 + b_A->size[0] * j];
      loc_data[i2] = c;
      BY_data[i2] = c;
    }
    loop_ub = b_A->size[0];
    if (b_A->size[0] < 1) {
      c = 0.0;
    } else {
      c = cblas_ddot((blasint)b_A->size[0], &loc_data[0], (blasint)1,
                     &BY_data[0], (blasint)1);
    }
    c = sqrt(c);
    /* 'gkmPWMlasso4:129' AA(:,j) = A(:,j)/normvec(j); */
    for (i2 = 0; i2 < loop_ub; i2++) {
      f_data[i2 + AA->size[0] * j] = A_data[i2 + b_A->size[0] * j] / c;
    }
  }
  emxFree_real_T(&GCmat);
  emxFree_cell_wrap_0(&p);
  emxInit_cell_wrap_1(&motclus);
  /* 'gkmPWMlasso4:131' fprintf('\n'); */
  printf("\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:132' fprintf('Clustering motifs\n'); */
  printf("Clustering motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:133' simmat = AA'*AA; */
  b_mtimes(AA, AA, comb);
  /*  clear AA */
  /* 'gkmPWMlasso4:135' motclus = clus_simmat_eig(simmat,corrCut); */
  clus_simmat_eig(comb, varargin_5, motclus);
  motclus_data = motclus->data;
  /* 'gkmPWMlasso4:136' fprintf('Number of motif clusters: %d\n',
   * int32(length(motclus))); */
  printf("Number of motif clusters: %d\n", motclus->size[0]);
  fflush(stdout);
  /* 'gkmPWMlasso4:138' fprintf('Selecting Motifs\n'); */
  printf("Selecting Motifs\n");
  fflush(stdout);
  /* 'gkmPWMlasso4:139' cfile2 = cfile-negvec/sum(negvec)*sum(cfile); */
  c = blockedSummation(cfile2, cfile2->size[0]);
  nfrac = blockedSummation(negvec, negvec->size[0]);
  emxFree_real_T(&AA);
  if (cfile2->size[0] == negvec->size[0]) {
    loop_ub = cfile2->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      cfile2_data[i1] -= negvec_data[i1] / nfrac * c;
    }
  } else {
    d_binary_expand_op(cfile2, negvec, nfrac, c);
    cfile2_data = cfile2->data;
  }
  /* 'gkmPWMlasso4:140' cfile2 = cfile2/std(cfile2); */
  c = b_std(cfile2);
  loop_ub = cfile2->size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    cfile2_data[i1] /= c;
  }
  emxInit_real_T(&B, 2);
  /* 'gkmPWMlasso4:141' B = zeros((4^k_svm)*lcnum, length(motclus)); */
  i1 = B->size[0] * B->size[1];
  B->size[0] = i;
  B->size[1] = motclus->size[0];
  emxEnsureCapacity_real_T(B, i1);
  B_data = B->data;
  loop_ub = i * motclus->size[0];
  for (i = 0; i < loop_ub; i++) {
    B_data[i] = 0.0;
  }
  /* 'gkmPWMlasso4:142' corrvec = zeros(n,1); */
  /* 'gkmPWMlasso4:143' zvec = zeros(n,1); */
  /* 'gkmPWMlasso4:144' Z = zeros(length(motclus),1); */
  i = Z->size[0];
  Z->size[0] = motclus->size[0];
  emxEnsureCapacity_real_T(Z, i);
  Z_data = Z->data;
  loop_ub = motclus->size[0];
  for (i = 0; i < loop_ub; i++) {
    Z_data[i] = 0.0;
  }
  /* 'gkmPWMlasso4:145' for i = 1:n */
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
    nfrac = (double)u1 * 10.0;
    if (nfrac <= b_A->size[0]) {
      k = (int)nfrac;
    } else {
      k = b_A->size[0];
    }
    i4 = b_A->size[0];
    c_loop_ub = b_A->size[0];
    if (1 > u1) {
      d_loop_ub = 0;
    } else {
      d_loop_ub = u1;
    }
  }
  for (b_i = 0; b_i < i; b_i++) {
    /*  [~,I] = sort(A(:,i),'descend'); */
    /*  zvec(i) = mean(cfile2(I(1:lcnum))); */
    /* 'gkmPWMlasso4:148' [~,I2] = maxk(A(:,i), lcnum*10); */
    i1 = A->size[0];
    A->size[0] = i4;
    emxEnsureCapacity_real_T(A, i1);
    xc_data = A->data;
    for (i1 = 0; i1 < c_loop_ub; i1++) {
      xc_data[i1] = A_data[i1 + b_A->size[0] * b_i];
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
    /* 'gkmPWMlasso4:149' zvec(i) = mean(cfile2(I2(1:lcnum))); */
    i1 = BY->size[0];
    BY->size[0] = d_loop_ub;
    emxEnsureCapacity_real_T(BY, i1);
    BY_data = BY->data;
    for (i1 = 0; i1 < d_loop_ub; i1++) {
      BY_data[i1] = cfile2_data[(int)loc_data[i1] - 1];
    }
    c = blockedSummation(BY, d_loop_ub) / (double)d_loop_ub;
    normvec_data[b_i] = c;
    /* 'gkmPWMlasso4:150' if zvec(i) < 0 */
    if (c < 0.0) {
      /*  Alternative to corr */
      /*  corrvec(i) = -1*corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10))); */
      /* 'gkmPWMlasso4:153' correlation_matrix = -1*corrcoef(A(I2,i),
       * cfile2(I2)); */
      /* 'gkmPWMlasso4:154' corrvec(i) = correlation_matrix(1,2); */
      i1 = A->size[0];
      A->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(A, i1);
      xc_data = A->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        xc_data[i1] = A_data[((int)loc_data[i1] + b_A->size[0] * b_i) - 1];
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
      /* 'gkmPWMlasso4:155' else */
      /*  Alternative to corr */
      /*  corrvec(i) = corr(A(I(1:lcnum*10),i), cfile2(I(1:lcnum*10))); */
      /* 'gkmPWMlasso4:158' correlation_matrix = corrcoef(A(I2,i), cfile2(I2));
       */
      /* 'gkmPWMlasso4:159' corrvec(i) = correlation_matrix(1,2); */
      i1 = A->size[0];
      A->size[0] = loc->size[0];
      emxEnsureCapacity_real_T(A, i1);
      xc_data = A->data;
      loop_ub = loc->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        xc_data[i1] = A_data[((int)loc_data[i1] + b_A->size[0] * b_i) - 1];
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
  /* 'gkmPWMlasso4:162' for i = 1:length(motclus) */
  i = motclus->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    /* 'gkmPWMlasso4:163' [a,b] = sort(zvec(motclus{i}),'descend'); */
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
    /* 'gkmPWMlasso4:164' f = find(a == a(1)); */
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
    i1 = negvec->size[0];
    negvec->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(negvec, i1);
    negvec_data = negvec->data;
    loop_ub = idx->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      negvec_data[i1] = matches_data[i1];
    }
    /* 'gkmPWMlasso4:165' if length(f) > 1 */
    if (negvec->size[0] > 1) {
      /* 'gkmPWMlasso4:166' [~,bb] = sort(abs(corrvec(motclus{i}(b(f)))),
       * 'descend'); */
      i1 = f->size[0] * f->size[1];
      f->size[0] = 1;
      f->size[1] = negvec->size[0];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      loop_ub = negvec->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = motclus_data[b_i]
                         .f1->data[(int)BY_data[(int)negvec_data[i1] - 1] - 1];
      }
      nx = negvec->size[0];
      i1 = indc->size[0];
      indc->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(indc, i1);
      indc_data = indc->data;
      for (k = 0; k < nx; k++) {
        indc_data[k] = fabs(GCmat_data[(int)f_data[k] - 1]);
      }
      sort(indc, idx);
      matches_data = idx->data;
      /* 'gkmPWMlasso4:167' b(1:length(f)) = b(bb); */
      i1 = f->size[0] * f->size[1];
      f->size[0] = 1;
      f->size[1] = negvec->size[0];
      emxEnsureCapacity_real_T(f, i1);
      f_data = f->data;
      loop_ub = negvec->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1] = BY_data[matches_data[i1] - 1];
      }
      loop_ub = f->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        BY_data[i1] = f_data[i1];
      }
    }
    /* 'gkmPWMlasso4:169' B(:,i) = A(:,motclus{i}(b(1))); */
    nx = (int)motclus_data[b_i].f1->data[(int)BY_data[0] - 1];
    loop_ub = b_A->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      B_data[i1 + B->size[0] * b_i] = A_data[i1 + b_A->size[0] * (nx - 1)];
    }
    /* 'gkmPWMlasso4:170' motclus{i} = motclus{i}(b); */
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
    /* 'gkmPWMlasso4:171' Z(i) = zvec(motclus{i}(1)); */
    Z_data[b_i] = normvec_data[(int)motclus_data[b_i].f1->data[0] - 1];
  }
  emxFree_real_T(&b_A);
  /* 'gkmPWMlasso4:173' f = find(abs(Z)>1); */
  nx = Z->size[0];
  i = indc->size[0];
  indc->size[0] = Z->size[0];
  emxEnsureCapacity_real_T(indc, i);
  indc_data = indc->data;
  for (k = 0; k < nx; k++) {
    indc_data[k] = fabs(Z_data[k]);
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
  i = negvec->size[0];
  negvec->size[0] = idx->size[0];
  emxEnsureCapacity_real_T(negvec, i);
  negvec_data = negvec->data;
  loop_ub = idx->size[0];
  for (i = 0; i < loop_ub; i++) {
    negvec_data[i] = matches_data[i];
  }
  /* 'gkmPWMlasso4:174' B = B(:,f); */
  nx = B->size[0] - 1;
  i = b_GCmat->size[0] * b_GCmat->size[1];
  b_GCmat->size[0] = B->size[0];
  b_GCmat->size[1] = negvec->size[0];
  emxEnsureCapacity_real_T(b_GCmat, i);
  normvec_data = b_GCmat->data;
  loop_ub = negvec->size[0];
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 <= nx; i1++) {
      normvec_data[i1 + b_GCmat->size[0] * i] =
          B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
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
  /* 'gkmPWMlasso4:175' B = B/std(B(:))'; */
  nx = B->size[0] * B->size[1];
  b_B = *B;
  b_motclus = nx;
  b_B.size = &b_motclus;
  b_B.numDimensions = 1;
  c = b_std(&b_B);
  loop_ub = B->size[0] * B->size[1];
  for (i = 0; i < loop_ub; i++) {
    B_data[i] /= c;
  }
  emxInit_cell_wrap_1(&tmp_motclus);
  /*  motclus = motclus(f); */
  /* 'gkmPWMlasso4:177' tmp_motclus = cell(length(f),1); */
  match_idx = negvec->size[0];
  i = tmp_motclus->size[0];
  tmp_motclus->size[0] = negvec->size[0];
  emxEnsureCapacity_cell_wrap_1(tmp_motclus, i);
  tmp_motclus_data = tmp_motclus->data;
  for (i = 0; i < match_idx; i++) {
    tmp_motclus_data[i].f1->size[0] = 1;
    tmp_motclus_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso4:178' tmp_motclus = coder.nullcopy(tmp_motclus); */
  /* 'gkmPWMlasso4:179' for idx = 1:length(f) */
  i = negvec->size[0];
  for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
    /* 'gkmPWMlasso4:180' tmp_motclus{idx} = motclus{f(idx)}; */
    i1 = tmp_motclus_data[b_loop_ub].f1->size[0] *
         tmp_motclus_data[b_loop_ub].f1->size[1];
    tmp_motclus_data[b_loop_ub].f1->size[0] = 1;
    tmp_motclus_data[b_loop_ub].f1->size[1] =
        motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->size[1];
    emxEnsureCapacity_real_T(tmp_motclus_data[b_loop_ub].f1, i1);
    loop_ub = motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      tmp_motclus_data[b_loop_ub].f1->data[i1] =
          motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->data[i1];
    }
  }
  emxInit_cell_wrap_1(&c_motclus);
  /* 'gkmPWMlasso4:182' motclus = cell(length(f),1); */
  match_idx = negvec->size[0];
  i = c_motclus->size[0];
  c_motclus->size[0] = negvec->size[0];
  emxEnsureCapacity_cell_wrap_1(c_motclus, i);
  b_motclus_data = c_motclus->data;
  for (i = 0; i < match_idx; i++) {
    b_motclus_data[i].f1->size[0] = 1;
    b_motclus_data[i].f1->size[1] = 0;
  }
  /* 'gkmPWMlasso4:183' motclus = coder.nullcopy(motclus); */
  i = motclus->size[0];
  motclus->size[0] = c_motclus->size[0];
  emxEnsureCapacity_cell_wrap_1(motclus, i);
  motclus_data = motclus->data;
  /* 'gkmPWMlasso4:184' for idx = 1:length(motclus) */
  i = c_motclus->size[0];
  for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
    /* 'gkmPWMlasso4:185' motclus{idx} = tmp_motclus{idx}; */
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
  /* 'gkmPWMlasso4:187' Z = Z(f); */
  i = BY->size[0];
  BY->size[0] = negvec->size[0];
  emxEnsureCapacity_real_T(BY, i);
  BY_data = BY->data;
  loop_ub = negvec->size[0];
  for (i = 0; i < loop_ub; i++) {
    BY_data[i] = Z_data[(int)negvec_data[i] - 1];
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
  /* 'gkmPWMlasso4:190' if d == 0 */
  emxInit_real_T(&OLS, 1);
  emxInit_real_T(&BB, 2);
  if (varargin_10 == 0.0) {
    /* 'gkmPWMlasso4:191' fprintf('Running LASSO\n'); */
    printf("Running LASSO\n");
    fflush(stdout);
    /*  [weigmat, FitInfo] = lasso_cvmat(B, cfile2,'DFmax',
     * length(Z),'Standardize', false, 'NumLambda', 20); */
    /* 'gkmPWMlasso4:193' [weigmat, FitInfo] = lasso_cvmat(B, cfile2, length(Z),
     * false, 20); */
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = B->size[1];
    emxEnsureCapacity_real_T(comb, i);
    GCmat_data = comb->data;
    loop_ub = B->size[0] * B->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      GCmat_data[i] = B_data[i];
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
    lasso_cvmat(comb, diffc, Z->size[0], weigmat, expl_temp_data, lk_size,
                b_expl_temp_data, expl_temp_size, &nfrac, a__3_DF_data,
                a__3_DF_size, c_expl_temp_data, b_expl_temp_size);
    xc_data = weigmat->data;
    /* 'gkmPWMlasso4:194' MSE = zeros(length(FitInfo.DF),1); */
    if (a__3_DF_size[1] == 0) {
      u1 = 0;
    } else {
      u1 = a__3_DF_size[1];
    }
    if (0 <= u1 - 1) {
      memset(&MSE_data[0], 0, u1 * sizeof(double));
    }
    /* 'gkmPWMlasso4:195' cnorm = cfile2'*cfile2; */
    if (cfile2->size[0] < 1) {
      nfrac = 0.0;
    } else {
      nfrac = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                         &cfile2_data[0], (blasint)1);
    }
    emxInit_cell_wrap_3_20(&F);
    /* 'gkmPWMlasso4:196' F = cell(numel(FitInfo.DF),1); */
    match_idx = a__3_DF_size[1];
    i = F.size[0];
    F.size[0] = a__3_DF_size[1];
    emxEnsureCapacity_cell_wrap_3(F.data, a__3_DF_size[1], i);
    for (i = 0; i < match_idx; i++) {
      F.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:197' F = coder.nullcopy(F); */
    /* 'gkmPWMlasso4:199' f = find(weigmat(:,1)~=0); */
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
    i = negvec->size[0];
    negvec->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(negvec, i);
    negvec_data = negvec->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      negvec_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso4:200' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, comb);
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (negvec->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)negvec->size[0], (blasint)1, (blasint)B->size[0],
                  1.0, &GCmat_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)negvec->size[0]);
    }
    mpower(comb, comb2);
    GCmat_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (indc->size[0] == 0)) {
      i = loc->size[0];
      loc->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(loc, i);
    } else {
      i = loc->size[0];
      loc->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &GCmat_data[0], (blasint)comb2->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &loc_data[0],
                  (blasint)comb2->size[0]);
    }
    i = OLS->size[0];
    OLS->size[0] = loc->size[0];
    emxEnsureCapacity_real_T(OLS, i);
    /* 'gkmPWMlasso4:202' for i = 1:length(FitInfo.DF) */
    if (a__3_DF_size[1] == 0) {
      b_loop_ub = 0;
    } else {
      b_loop_ub = a__3_DF_size[1];
    }
    if (0 <= b_loop_ub - 1) {
      i5 = weigmat->size[0];
      e_loop_ub = weigmat->size[0];
      f_loop_ub = B->size[0];
      g_loop_ub = B->size[0];
      h_loop_ub = B->size[0];
      i6 = B->size[0];
      i_loop_ub = B->size[0];
      i7 = B->size[0];
    }
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      /* 'gkmPWMlasso4:203' f = find(weigmat(:,i)~=0); */
      i = b_loc->size[0];
      b_loc->size[0] = i5;
      emxEnsureCapacity_boolean_T(b_loc, i);
      b_loc_data = b_loc->data;
      for (i = 0; i < e_loop_ub; i++) {
        b_loc_data[i] = (xc_data[i + weigmat->size[0] * b_i] != 0.0);
      }
      b_eml_find(b_loc, idx);
      matches_data = idx->data;
      i = negvec->size[0];
      negvec->size[0] = idx->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      loop_ub = idx->size[0];
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = matches_data[i];
      }
      /* 'gkmPWMlasso4:204' F{i} = f; */
      i = F.data[b_i].f1->size[0];
      F.data[b_i].f1->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(F.data[b_i].f1, i);
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        F.data[b_i].f1->data[i] = negvec_data[i];
      }
      /* 'gkmPWMlasso4:205' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
      i = b_GCmat->size[0] * b_GCmat->size[1];
      b_GCmat->size[0] = f_loop_ub;
      b_GCmat->size[1] = negvec->size[0];
      emxEnsureCapacity_real_T(b_GCmat, i);
      normvec_data = b_GCmat->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < f_loop_ub; i1++) {
          normvec_data[i1 + b_GCmat->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
        }
      }
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = g_loop_ub;
      comb2->size[1] = negvec->size[0];
      emxEnsureCapacity_real_T(comb2, i);
      GCmat_data = comb2->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < g_loop_ub; i1++) {
          GCmat_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
        }
      }
      b_mtimes(b_GCmat, comb2, comb);
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = h_loop_ub;
      comb2->size[1] = negvec->size[0];
      emxEnsureCapacity_real_T(comb2, i);
      GCmat_data = comb2->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < h_loop_ub; i1++) {
          GCmat_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
        }
      }
      if ((i6 == 0) || (negvec->size[0] == 0) || (cfile2->size[0] == 0)) {
        i = indc->size[0];
        indc->size[0] = negvec->size[0];
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        loop_ub = negvec->size[0];
        for (i = 0; i < loop_ub; i++) {
          indc_data[i] = 0.0;
        }
      } else {
        i = indc->size[0];
        indc->size[0] = negvec->size[0];
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                    (blasint)negvec->size[0], (blasint)1, (blasint)B->size[0],
                    1.0, &GCmat_data[0], (blasint)B->size[0], &cfile2_data[0],
                    (blasint)cfile2->size[0], 0.0, &indc_data[0],
                    (blasint)negvec->size[0]);
      }
      mpower(comb, comb2);
      GCmat_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (indc->size[0] == 0)) {
        i = OLS->size[0];
        OLS->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(OLS, i);
        OLS_data = OLS->data;
        loop_ub = comb2->size[0];
        for (i = 0; i < loop_ub; i++) {
          OLS_data[i] = 0.0;
        }
      } else {
        i = OLS->size[0];
        OLS->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(OLS, i);
        OLS_data = OLS->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &GCmat_data[0],
            (blasint)comb2->size[0], &indc_data[0], (blasint)indc->size[0], 0.0,
            &OLS_data[0], (blasint)comb2->size[0]);
      }
      /* 'gkmPWMlasso4:206' res = B(:,f)*OLS; */
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = i_loop_ub;
      comb2->size[1] = negvec->size[0];
      emxEnsureCapacity_real_T(comb2, i);
      GCmat_data = comb2->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < i_loop_ub; i1++) {
          GCmat_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
        }
      }
      if ((i7 == 0) || (negvec->size[0] == 0) || (OLS->size[0] == 0)) {
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
                    (blasint)B->size[0], (blasint)1, (blasint)negvec->size[0],
                    1.0, &GCmat_data[0], (blasint)B->size[0], &OLS_data[0],
                    (blasint)OLS->size[0], 0.0, &loc_data[0],
                    (blasint)B->size[0]);
      }
      /* 'gkmPWMlasso4:207' if sum(abs(res)>0) */
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
        /* 'gkmPWMlasso4:208' MSE(i) = (cfile2'*res)^2/(res'*res)/cnorm; */
        if (cfile2->size[0] < 1) {
          a = 0.0;
        } else {
          a = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                         &loc_data[0], (blasint)1);
        }
        if (loc->size[0] < 1) {
          c = 0.0;
        } else {
          c = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                         &loc_data[0], (blasint)1);
        }
        MSE_data[b_i] = a * a / c / nfrac;
      }
    }
    emxFree_real_T(&weigmat);
    emxInit_cell_wrap_3_1x19(&FF);
    /* 'gkmPWMlasso4:211' FF = cell(1,length(MSE)-1); */
    match_idx = u1 - 1;
    i = FF.size[0] * FF.size[1];
    FF.size[0] = 1;
    FF.size[1] = u1 - 1;
    emxEnsureCapacity_cell_wrap_31(FF.data, FF.size, i);
    for (i = 0; i < match_idx; i++) {
      FF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:212' FF = coder.nullcopy(FF); */
    /* 'gkmPWMlasso4:213' MSE2 = zeros(1,length(MSE)-1); */
    MSE2_size[0] = 1;
    if (0 <= match_idx - 1) {
      memset(&MSE2_data[0], 0, match_idx * sizeof(double));
    }
    /* 'gkmPWMlasso4:214' count = 1; */
    nx = 0;
    /* 'gkmPWMlasso4:215' for i = 1:length(MSE)-1 */
    for (b_i = 0; b_i <= u1 - 2; b_i++) {
      /* 'gkmPWMlasso4:216' if numel(setdiff(F{i}, F{i+1})) > 0 */
      b_do_vectors(F.data[b_i].f1, F.data[b_i + 1].f1, loc, idx, &match_idx);
      if (loc->size[0] > 0) {
        /* 'gkmPWMlasso4:217' FF{count} = setdiff(F{i}, F{i+1}); */
        b_do_vectors(F.data[b_i].f1, F.data[b_i + 1].f1, FF.data[nx].f1, idx,
                     &match_idx);
        /* 'gkmPWMlasso4:218' MSE2(count) = (MSE(i)-MSE(i+1)); */
        MSE2_data[nx] = MSE_data[b_i] - MSE_data[b_i + 1];
        /* 'gkmPWMlasso4:219' count = count +1; */
        nx++;
      }
    }
    emxFree_cell_wrap_3_20(&F);
    emxInit_cell_wrap_3_1x19(&newFF);
    /* 'gkmPWMlasso4:222' MSE2 = MSE2(1:count-1); */
    if (1 > nx) {
      MSE2_size[1] = 0;
    } else {
      MSE2_size[1] = nx;
    }
    /* 'gkmPWMlasso4:223' newFF = cell(1,count-1); */
    i = newFF.size[0] * newFF.size[1];
    newFF.size[0] = 1;
    newFF.size[1] = nx;
    emxEnsureCapacity_cell_wrap_31(newFF.data, newFF.size, i);
    for (i = 0; i < nx; i++) {
      newFF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:224' newFF = coder.nullcopy(newFF); */
    /* 'gkmPWMlasso4:225' for idx = 1:count-1 */
    for (b_loop_ub = 0; b_loop_ub < nx; b_loop_ub++) {
      /* 'gkmPWMlasso4:226' newFF{idx} = FF{idx}; */
      loop_ub = FF.data[b_loop_ub].f1->size[0];
      i = newFF.data[b_loop_ub].f1->size[0];
      newFF.data[b_loop_ub].f1->size[0] = FF.data[b_loop_ub].f1->size[0];
      emxEnsureCapacity_real_T(newFF.data[b_loop_ub].f1, i);
      for (i = 0; i < loop_ub; i++) {
        newFF.data[b_loop_ub].f1->data[i] = FF.data[b_loop_ub].f1->data[i];
      }
    }
    emxInit_cell_wrap_3_1x19(&b_FF);
    /* 'gkmPWMlasso4:228' FF = cell(1,count-1); */
    i = b_FF.size[0] * b_FF.size[1];
    b_FF.size[0] = 1;
    b_FF.size[1] = nx;
    emxEnsureCapacity_cell_wrap_31(b_FF.data, b_FF.size, i);
    for (i = 0; i < nx; i++) {
      b_FF.data[i].f1->size[0] = 0;
    }
    /* 'gkmPWMlasso4:229' FF = coder.nullcopy(FF); */
    i = FF.size[0] * FF.size[1];
    FF.size[0] = 1;
    FF.size[1] = b_FF.size[1];
    emxEnsureCapacity_cell_wrap_31(FF.data, FF.size, i);
    /* 'gkmPWMlasso4:230' for idx = 1:count-1 */
    emxFree_cell_wrap_3_1x19(&b_FF);
    for (b_loop_ub = 0; b_loop_ub < nx; b_loop_ub++) {
      /* 'gkmPWMlasso4:231' FF{idx} = newFF{idx}; */
      loop_ub = newFF.data[b_loop_ub].f1->size[0];
      i = FF.data[b_loop_ub].f1->size[0];
      FF.data[b_loop_ub].f1->size[0] = newFF.data[b_loop_ub].f1->size[0];
      emxEnsureCapacity_real_T(FF.data[b_loop_ub].f1, i);
      for (i = 0; i < loop_ub; i++) {
        FF.data[b_loop_ub].f1->data[i] = newFF.data[b_loop_ub].f1->data[i];
      }
    }
    emxFree_cell_wrap_3_1x19(&newFF);
    /* 'gkmPWMlasso4:234' res = B*((B.'*B)^-1*(B.'*cfile2)); */
    b_mtimes(B, B, comb);
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
    mpower(comb, comb2);
    GCmat_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (indc->size[0] == 0)) {
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
                  1.0, &GCmat_data[0], (blasint)comb2->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &BY_data[0],
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
    /* 'gkmPWMlasso4:235' csm = (cfile2'*res)^2/(res'*res)/cnorm; */
    if (cfile2->size[0] < 1) {
      a = 0.0;
    } else {
      a = cblas_ddot((blasint)cfile2->size[0], &cfile2_data[0], (blasint)1,
                     &loc_data[0], (blasint)1);
    }
    if (loc->size[0] < 1) {
      c = 0.0;
    } else {
      c = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                     &loc_data[0], (blasint)1);
    }
    nfrac = a * a / c / nfrac;
    /* 'gkmPWMlasso4:236' [a,b] = sort(MSE2,'descend'); */
    e_sort(MSE2_data, MSE2_size, iidx_data, lk_size);
    /* 'gkmPWMlasso4:237' cs = cumsum(a); */
    useConstantDim(MSE2_data, MSE2_size);
    /* 'gkmPWMlasso4:238' cs = cs/csm; */
    loop_ub = MSE2_size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      MSE2_data[i] /= nfrac;
    }
    /* 'gkmPWMlasso4:239' f = find(cs>0.9); */
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
    /* 'gkmPWMlasso4:241' if isempty(f) */
    emxInit_real_T(&b_f, 2);
    if (f->size[1] == 0) {
      /* 'gkmPWMlasso4:242' f = 1:length(OLS); */
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
      i = b_f->size[0] * b_f->size[1];
      b_f->size[0] = 1;
      b_f->size[1] = f->size[1];
      emxEnsureCapacity_real_T(b_f, i);
      A_data = b_f->data;
      loop_ub = f->size[1];
      for (i = 0; i < loop_ub; i++) {
        A_data[i] = f_data[i];
      }
    } else {
      emxInit_cell_wrap_3(&b_newFF, 2);
      /* 'gkmPWMlasso4:243' else */
      /*  FF = FF(b(1:f(1))); */
      /* 'gkmPWMlasso4:245' endIdx = f(1); */
      /* 'gkmPWMlasso4:246' newFF = cell(1, endIdx); */
      match_idx = (int)f_data[0];
      i = b_newFF->size[0] * b_newFF->size[1];
      b_newFF->size[0] = 1;
      b_newFF->size[1] = match_idx;
      emxEnsureCapacity_cell_wrap_32(b_newFF, i);
      newFF_data = b_newFF->data;
      for (i = 0; i < match_idx; i++) {
        newFF_data[i].f1->size[0] = 0;
      }
      /* 'gkmPWMlasso4:247' newFF = coder.nullcopy(newFF); */
      /* 'gkmPWMlasso4:248' for idx = 1:endIdx */
      for (b_loop_ub = 0; b_loop_ub < match_idx; b_loop_ub++) {
        /* 'gkmPWMlasso4:249' newFF{idx} = FF{b(idx)}; */
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
      /* 'gkmPWMlasso4:251' FF = newFF; */
      /* 'gkmPWMlasso4:253' f = []; */
      negvec->size[0] = 0;
      /* 'gkmPWMlasso4:254' for i = 1:length(FF) */
      i = b_newFF->size[1];
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso4:255' f = [f;FF{i}]; */
        i1 = negvec->size[0];
        loop_ub = newFF_data[b_i].f1->size[0];
        i2 = negvec->size[0];
        negvec->size[0] += newFF_data[b_i].f1->size[0];
        emxEnsureCapacity_real_T(negvec, i2);
        negvec_data = negvec->data;
        for (i2 = 0; i2 < loop_ub; i2++) {
          negvec_data[i1 + i2] = newFF_data[b_i].f1->data[i2];
        }
      }
      emxFree_cell_wrap_3(&b_newFF);
      /* 'gkmPWMlasso4:257' f=unique(f); */
      unique_vector(negvec, loc);
      loc_data = loc->data;
      i = b_f->size[0] * b_f->size[1];
      b_f->size[0] = loc->size[0];
      b_f->size[1] = 1;
      emxEnsureCapacity_real_T(b_f, i);
      A_data = b_f->data;
      loop_ub = loc->size[0];
      for (i = 0; i < loop_ub; i++) {
        A_data[i] = loc_data[i];
      }
    }
    emxFree_cell_wrap_3_1x19(&FF);
    /* f = find(weigmat(:,1)~=0); */
    /* 'gkmPWMlasso4:262' F = length(f); */
    /* 'gkmPWMlasso4:263' fprintf('Selecting top motifs\n'); */
    printf("Selecting top motifs\n");
    fflush(stdout);
    /* 'gkmPWMlasso4:264' ind = true; */
    empty_non_axis_sizes = true;
    /* 'gkmPWMlasso4:265' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    match_idx = b_f->size[0] * b_f->size[1];
    nx = b_f->size[0] * b_f->size[1];
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = match_idx;
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = nx;
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, comb);
    match_idx = b_f->size[0] * b_f->size[1];
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = match_idx;
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (b_f->size[0] * b_f->size[1] == 0) ||
        (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = b_f->size[0] * b_f->size[1];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = b_f->size[0] * b_f->size[1];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = b_f->size[0] * b_f->size[1];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)(b_f->size[0] * b_f->size[1]), (blasint)1,
                  (blasint)B->size[0], 1.0, &GCmat_data[0], (blasint)B->size[0],
                  &cfile2_data[0], (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)(b_f->size[0] * b_f->size[1]));
    }
    mpower(comb, comb2);
    GCmat_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (indc->size[0] == 0)) {
      i = OLS->size[0];
      OLS->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      loop_ub = comb2->size[0];
      for (i = 0; i < loop_ub; i++) {
        OLS_data[i] = 0.0;
      }
    } else {
      i = OLS->size[0];
      OLS->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &GCmat_data[0], (blasint)comb2->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &OLS_data[0],
                  (blasint)comb2->size[0]);
    }
    emxInit_real_T(&Pweig, 2);
    /* 'gkmPWMlasso4:266' Pweig = Z(f); */
    i = Pweig->size[0] * Pweig->size[1];
    Pweig->size[0] = b_f->size[0];
    Pweig->size[1] = b_f->size[1];
    emxEnsureCapacity_real_T(Pweig, i);
    Pweig_data = Pweig->data;
    loop_ub = b_f->size[0] * b_f->size[1];
    for (i = 0; i < loop_ub; i++) {
      Pweig_data[i] = Z_data[(int)A_data[i] - 1];
    }
    /* 'gkmPWMlasso4:267' while ind */
    while (empty_non_axis_sizes) {
      /* 'gkmPWMlasso4:268' ff = []; */
      loc->size[0] = 0;
      /* 'gkmPWMlasso4:269' for i = 1:length(f) */
      if ((b_f->size[0] == 0) || (b_f->size[1] == 0)) {
        u1 = 0;
      } else {
        nx = b_f->size[0];
        u1 = b_f->size[1];
        if (nx >= u1) {
          u1 = nx;
        }
      }
      for (b_i = 0; b_i < u1; b_i++) {
        /* 'gkmPWMlasso4:270' if sign(OLS(i)) ~= sign(Pweig(i)) */
        GCpos1 = OLS_data[b_i];
        if (OLS_data[b_i] < 0.0) {
          GCpos1 = -1.0;
        } else if (OLS_data[b_i] > 0.0) {
          GCpos1 = 1.0;
        }
        nfrac = Pweig_data[b_i];
        if (Pweig_data[b_i] < 0.0) {
          nfrac = -1.0;
        } else if (Pweig_data[b_i] > 0.0) {
          nfrac = 1.0;
        }
        if (GCpos1 != nfrac) {
          /* 'gkmPWMlasso4:271' ff = [ff;i]; */
          i = loc->size[0];
          i1 = loc->size[0];
          loc->size[0]++;
          emxEnsureCapacity_real_T(loc, i1);
          loc_data = loc->data;
          loc_data[i] = (double)b_i + 1.0;
        }
      }
      /* 'gkmPWMlasso4:274' if length(ff) > 0 */
      if (loc->size[0] > 0) {
        /* 'gkmPWMlasso4:275' f(ff) = []; */
        i = idx->size[0];
        idx->size[0] = loc->size[0];
        emxEnsureCapacity_int32_T(idx, i);
        matches_data = idx->data;
        loop_ub = loc->size[0];
        for (i = 0; i < loop_ub; i++) {
          matches_data[i] = (int)loc_data[i];
        }
        c_nullAssignment(b_f, idx);
        A_data = b_f->data;
        /* 'gkmPWMlasso4:276' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
        match_idx = b_f->size[0] * b_f->size[1];
        nx = b_f->size[0] * b_f->size[1];
        loop_ub = B->size[0];
        i = b_GCmat->size[0] * b_GCmat->size[1];
        b_GCmat->size[0] = B->size[0];
        b_GCmat->size[1] = match_idx;
        emxEnsureCapacity_real_T(b_GCmat, i);
        normvec_data = b_GCmat->data;
        for (i = 0; i < match_idx; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            normvec_data[i1 + b_GCmat->size[0] * i] =
                B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
          }
        }
        loop_ub = B->size[0];
        i = comb2->size[0] * comb2->size[1];
        comb2->size[0] = B->size[0];
        comb2->size[1] = nx;
        emxEnsureCapacity_real_T(comb2, i);
        GCmat_data = comb2->data;
        for (i = 0; i < nx; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            GCmat_data[i1 + comb2->size[0] * i] =
                B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
          }
        }
        b_mtimes(b_GCmat, comb2, comb);
        match_idx = b_f->size[0] * b_f->size[1];
        loop_ub = B->size[0];
        i = comb2->size[0] * comb2->size[1];
        comb2->size[0] = B->size[0];
        comb2->size[1] = match_idx;
        emxEnsureCapacity_real_T(comb2, i);
        GCmat_data = comb2->data;
        for (i = 0; i < match_idx; i++) {
          for (i1 = 0; i1 < loop_ub; i1++) {
            GCmat_data[i1 + comb2->size[0] * i] =
                B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
          }
        }
        if ((B->size[0] == 0) || (b_f->size[0] * b_f->size[1] == 0) ||
            (cfile2->size[0] == 0)) {
          i = indc->size[0];
          indc->size[0] = b_f->size[0] * b_f->size[1];
          emxEnsureCapacity_real_T(indc, i);
          indc_data = indc->data;
          loop_ub = b_f->size[0] * b_f->size[1];
          for (i = 0; i < loop_ub; i++) {
            indc_data[i] = 0.0;
          }
        } else {
          i = indc->size[0];
          indc->size[0] = b_f->size[0] * b_f->size[1];
          emxEnsureCapacity_real_T(indc, i);
          indc_data = indc->data;
          cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                      (blasint)(b_f->size[0] * b_f->size[1]), (blasint)1,
                      (blasint)B->size[0], 1.0, &GCmat_data[0],
                      (blasint)B->size[0], &cfile2_data[0],
                      (blasint)cfile2->size[0], 0.0, &indc_data[0],
                      (blasint)(b_f->size[0] * b_f->size[1]));
        }
        mpower(comb, comb2);
        GCmat_data = comb2->data;
        if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
            (indc->size[0] == 0)) {
          i = OLS->size[0];
          OLS->size[0] = comb2->size[0];
          emxEnsureCapacity_real_T(OLS, i);
          OLS_data = OLS->data;
          loop_ub = comb2->size[0];
          for (i = 0; i < loop_ub; i++) {
            OLS_data[i] = 0.0;
          }
        } else {
          i = OLS->size[0];
          OLS->size[0] = comb2->size[0];
          emxEnsureCapacity_real_T(OLS, i);
          OLS_data = OLS->data;
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                      (blasint)comb2->size[0], (blasint)1,
                      (blasint)comb2->size[1], 1.0, &GCmat_data[0],
                      (blasint)comb2->size[0], &indc_data[0],
                      (blasint)indc->size[0], 0.0, &OLS_data[0],
                      (blasint)comb2->size[0]);
        }
        /* 'gkmPWMlasso4:277' Pweig = Z(f); */
        i = Pweig->size[0] * Pweig->size[1];
        Pweig->size[0] = b_f->size[0];
        Pweig->size[1] = b_f->size[1];
        emxEnsureCapacity_real_T(Pweig, i);
        Pweig_data = Pweig->data;
        loop_ub = b_f->size[0] * b_f->size[1];
        for (i = 0; i < loop_ub; i++) {
          Pweig_data[i] = Z_data[(int)A_data[i] - 1];
        }
      } else {
        /* 'gkmPWMlasso4:278' else */
        /* 'gkmPWMlasso4:279' ind = false; */
        empty_non_axis_sizes = false;
      }
    }
    /* 'gkmPWMlasso4:282' BB = B(:,f); */
    match_idx = b_f->size[0] * b_f->size[1];
    loop_ub = B->size[0];
    i = BB->size[0] * BB->size[1];
    BB->size[0] = B->size[0];
    BB->size[1] = match_idx;
    emxEnsureCapacity_real_T(BB, i);
    f_data = BB->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1 + BB->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
      }
    }
    /* 'gkmPWMlasso4:283' BX = B(:,f)'*B(:,f); */
    match_idx = b_f->size[0] * b_f->size[1];
    nx = b_f->size[0] * b_f->size[1];
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = match_idx;
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = nx;
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, xc);
    xc_data = xc->data;
    /* 'gkmPWMlasso4:284' BY = B(:,f)'*cfile2; */
    match_idx = b_f->size[0] * b_f->size[1];
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = match_idx;
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    for (i = 0; i < match_idx; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)A_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (b_f->size[0] * b_f->size[1] == 0) ||
        (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = b_f->size[0] * b_f->size[1];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = b_f->size[0] * b_f->size[1];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = b_f->size[0] * b_f->size[1];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)(b_f->size[0] * b_f->size[1]), (blasint)1,
                  (blasint)B->size[0], 1.0, &GCmat_data[0], (blasint)B->size[0],
                  &cfile2_data[0], (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)(b_f->size[0] * b_f->size[1]));
    }
    /* 'gkmPWMlasso4:285' E = zeros(length(f),1); */
    if ((b_f->size[0] == 0) || (b_f->size[1] == 0)) {
      u1 = 0;
    } else {
      nx = b_f->size[0];
      u1 = b_f->size[1];
      if (nx >= u1) {
        u1 = nx;
      }
    }
    i = normvec->size[0];
    normvec->size[0] = u1;
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    for (i = 0; i < u1; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso4:286' for i = 1:length(f) */
    if ((b_f->size[0] == 0) || (b_f->size[1] == 0)) {
      u1 = 0;
    } else {
      nx = b_f->size[0];
      u1 = b_f->size[1];
      if (nx >= u1) {
        u1 = nx;
      }
      j_loop_ub = xc->size[0] * xc->size[1];
      k_loop_ub = BY->size[0];
      nxout = BY->size[0] - 1;
      l_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < u1; b_i++) {
      /* 'gkmPWMlasso4:287' B = BB; */
      /* 'gkmPWMlasso4:288' BBX = BX; */
      /* 'gkmPWMlasso4:289' BBY = BY; */
      /* 'gkmPWMlasso4:290' B(:,i) = []; */
      /* 'gkmPWMlasso4:291' BBX(:,i) = []; */
      /* 'gkmPWMlasso4:292' BBX(i,:) = []; */
      /* 'gkmPWMlasso4:293' BBY(i) = []; */
      /* 'gkmPWMlasso4:294' res = cfile2-B*(BBX^-1*BBY); */
      i = comb->size[0] * comb->size[1];
      comb->size[0] = xc->size[0];
      comb->size[1] = xc->size[1];
      emxEnsureCapacity_real_T(comb, i);
      GCmat_data = comb->data;
      for (i = 0; i < j_loop_ub; i++) {
        GCmat_data[i] = xc_data[i];
      }
      d_nullAssignment(comb, b_i + 1);
      b_nullAssignment(comb, b_i + 1);
      mpower(comb, comb2);
      GCmat_data = comb2->data;
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
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (loc->size[0] == 0)) {
        i = indc->size[0];
        indc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        loop_ub = comb2->size[0];
        for (i = 0; i < loop_ub; i++) {
          indc_data[i] = 0.0;
        }
      } else {
        i = indc->size[0];
        indc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(indc, i);
        indc_data = indc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &GCmat_data[0],
            (blasint)comb2->size[0], &loc_data[0], (blasint)loc->size[0], 0.0,
            &indc_data[0], (blasint)comb2->size[0]);
      }
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = BB->size[0];
      comb2->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb2, i);
      GCmat_data = comb2->data;
      for (i = 0; i < l_loop_ub; i++) {
        GCmat_data[i] = f_data[i];
      }
      d_nullAssignment(comb2, b_i + 1);
      GCmat_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (indc->size[0] == 0)) {
        i = loc->size[0];
        loc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(loc, i);
        loc_data = loc->data;
        loop_ub = comb2->size[0];
        for (i = 0; i < loop_ub; i++) {
          loc_data[i] = 0.0;
        }
      } else {
        i = loc->size[0];
        loc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(loc, i);
        loc_data = loc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &GCmat_data[0],
            (blasint)comb2->size[0], &indc_data[0], (blasint)indc->size[0], 0.0,
            &loc_data[0], (blasint)comb2->size[0]);
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
      /* 'gkmPWMlasso4:295' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        c = 0.0;
      } else {
        c = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                       &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(c);
    }
    /* 'gkmPWMlasso4:297' res = cfile2-BB*OLS; */
    if ((B->size[0] == 0) || (b_f->size[0] * b_f->size[1] == 0) ||
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
                  (blasint)(b_f->size[0] * b_f->size[1]), 1.0, &f_data[0],
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
    /* 'gkmPWMlasso4:298' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      c = 0.0;
    } else {
      c = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                     &loc_data[0], (blasint)1);
    }
    nfrac = sqrt(c);
    /* 'gkmPWMlasso4:299' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - nfrac) / nfrac;
    }
    /*  motclus = motclus(f); */
    /* 'gkmPWMlasso4:301' mylen = length(f); */
    if ((b_f->size[0] == 0) || (b_f->size[1] == 0)) {
      u1 = 0;
    } else {
      nx = b_f->size[0];
      u1 = b_f->size[1];
      if (nx >= u1) {
        u1 = nx;
      }
    }
    /* 'gkmPWMlasso4:302' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso4:303' for idx=1:mylen */
    i = c_motclus->size[0];
    c_motclus->size[0] = u1;
    emxEnsureCapacity_cell_wrap_1(c_motclus, i);
    b_motclus_data = c_motclus->data;
    for (b_loop_ub = 0; b_loop_ub < u1; b_loop_ub++) {
      /* 'gkmPWMlasso4:304' newMotclus{idx} = motclus{f(idx)}; */
      i = b_motclus_data[b_loop_ub].f1->size[0] *
          b_motclus_data[b_loop_ub].f1->size[1];
      b_motclus_data[b_loop_ub].f1->size[0] = 1;
      b_motclus_data[b_loop_ub].f1->size[1] =
          motclus_data[(int)A_data[b_loop_ub] - 1].f1->size[1];
      emxEnsureCapacity_real_T(b_motclus_data[b_loop_ub].f1, i);
      loop_ub = motclus_data[(int)A_data[b_loop_ub] - 1].f1->size[1];
      for (i = 0; i < loop_ub; i++) {
        b_motclus_data[b_loop_ub].f1->data[i] =
            motclus_data[(int)A_data[b_loop_ub] - 1].f1->data[i];
      }
    }
    emxFree_real_T(&b_f);
    /* 'gkmPWMlasso4:306' motclus = newMotclus; */
    /* 'gkmPWMlasso4:308' f = find(E/max(E) >= 0.01); */
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
    i = negvec->size[0];
    negvec->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(negvec, i);
    negvec_data = negvec->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      negvec_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso4:309' B = BB; */
    i = B->size[0] * B->size[1];
    B->size[0] = BB->size[0];
    B->size[1] = BB->size[1];
    emxEnsureCapacity_real_T(B, i);
    B_data = B->data;
    loop_ub = BB->size[0] * BB->size[1];
    for (i = 0; i < loop_ub; i++) {
      B_data[i] = f_data[i];
    }
    /* 'gkmPWMlasso4:310' BB = B(:,f); */
    match_idx = BB->size[0] - 1;
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = BB->size[0];
    b_GCmat->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    loop_ub = negvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      for (i1 = 0; i1 <= match_idx; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            f_data[i1 + BB->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    i = BB->size[0] * BB->size[1];
    BB->size[0] = b_GCmat->size[0];
    BB->size[1] = b_GCmat->size[1];
    emxEnsureCapacity_real_T(BB, i);
    f_data = BB->data;
    loop_ub = b_GCmat->size[0] * b_GCmat->size[1];
    for (i = 0; i < loop_ub; i++) {
      f_data[i] = normvec_data[i];
    }
    /* 'gkmPWMlasso4:311' BX = B(:,f)'*B(:,f); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, xc);
    xc_data = xc->data;
    /* 'gkmPWMlasso4:312' BY = B(:,f)'*cfile2; */
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (negvec->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)negvec->size[0], (blasint)1, (blasint)B->size[0],
                  1.0, &GCmat_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)negvec->size[0]);
    }
    /* 'gkmPWMlasso4:313' E = zeros(length(f),1); */
    i = normvec->size[0];
    normvec->size[0] = negvec->size[0];
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    loop_ub = negvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso4:314' for i = 1:length(f) */
    i = negvec->size[0];
    if (0 <= negvec->size[0] - 1) {
      m_loop_ub = xc->size[0] * xc->size[1];
      n_loop_ub = BY->size[0];
      b_nxout = BY->size[0] - 1;
      o_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:315' B = BB; */
      /* 'gkmPWMlasso4:316' BBX = BX; */
      /* 'gkmPWMlasso4:317' BBY = BY; */
      /* 'gkmPWMlasso4:318' B(:,i) = []; */
      /* 'gkmPWMlasso4:319' BBX(:,i) = []; */
      /* 'gkmPWMlasso4:320' BBX(i,:) = []; */
      /* 'gkmPWMlasso4:321' BBY(i) = []; */
      /* 'gkmPWMlasso4:322' res = cfile2-B*(BBX^-1*BBY); */
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = xc->size[0];
      comb->size[1] = xc->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      GCmat_data = comb->data;
      for (i1 = 0; i1 < m_loop_ub; i1++) {
        GCmat_data[i1] = xc_data[i1];
      }
      d_nullAssignment(comb, b_i + 1);
      b_nullAssignment(comb, b_i + 1);
      mpower(comb, comb2);
      GCmat_data = comb2->data;
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
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (loc->size[0] == 0)) {
        i1 = indc->size[0];
        indc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(indc, i1);
        indc_data = indc->data;
        loop_ub = comb2->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          indc_data[i1] = 0.0;
        }
      } else {
        i1 = indc->size[0];
        indc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(indc, i1);
        indc_data = indc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &GCmat_data[0],
            (blasint)comb2->size[0], &loc_data[0], (blasint)loc->size[0], 0.0,
            &indc_data[0], (blasint)comb2->size[0]);
      }
      i1 = comb2->size[0] * comb2->size[1];
      comb2->size[0] = BB->size[0];
      comb2->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb2, i1);
      GCmat_data = comb2->data;
      for (i1 = 0; i1 < o_loop_ub; i1++) {
        GCmat_data[i1] = f_data[i1];
      }
      d_nullAssignment(comb2, b_i + 1);
      GCmat_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (indc->size[0] == 0)) {
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
            (blasint)1, (blasint)comb2->size[1], 1.0, &GCmat_data[0],
            (blasint)comb2->size[0], &indc_data[0], (blasint)indc->size[0], 0.0,
            &loc_data[0], (blasint)comb2->size[0]);
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
      /* 'gkmPWMlasso4:323' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        c = 0.0;
      } else {
        c = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                       &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(c);
    }
    /* 'gkmPWMlasso4:325' OLS = (BB.'*BB)^-1*(BB.'*cfile2); */
    b_mtimes(BB, BB, comb);
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
                  (blasint)1, (blasint)BB->size[0], 1.0, &f_data[0],
                  (blasint)BB->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)BB->size[1]);
    }
    mpower(comb, comb2);
    GCmat_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (indc->size[0] == 0)) {
      i = OLS->size[0];
      OLS->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      loop_ub = comb2->size[0];
      for (i = 0; i < loop_ub; i++) {
        OLS_data[i] = 0.0;
      }
    } else {
      i = OLS->size[0];
      OLS->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &GCmat_data[0], (blasint)comb2->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &OLS_data[0],
                  (blasint)comb2->size[0]);
    }
    /* 'gkmPWMlasso4:326' Pweig = Pweig(f); */
    /* 'gkmPWMlasso4:327' res = cfile2-BB*OLS; */
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
                  &f_data[0], (blasint)BB->size[0], &OLS_data[0],
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
    /* 'gkmPWMlasso4:328' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      c = 0.0;
    } else {
      c = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                     &loc_data[0], (blasint)1);
    }
    nfrac = sqrt(c);
    /* 'gkmPWMlasso4:329' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - nfrac) / nfrac;
    }
    /* 'gkmPWMlasso4:330' for i = 1:length(motclus) */
    i = c_motclus->size[0];
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:331' motclus{i} = indvec(motclus{i})'; */
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
    /* 'gkmPWMlasso4:334' mylen = length(f); */
    /* 'gkmPWMlasso4:335' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso4:336' for idx=1:mylen */
    i = negvec->size[0];
    i1 = tmp_motclus->size[0];
    tmp_motclus->size[0] = negvec->size[0];
    emxEnsureCapacity_cell_wrap_1(tmp_motclus, i1);
    tmp_motclus_data = tmp_motclus->data;
    for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
      /* 'gkmPWMlasso4:337' newMotclus{idx} = motclus{f(idx)}; */
      i1 = tmp_motclus_data[b_loop_ub].f1->size[0] *
           tmp_motclus_data[b_loop_ub].f1->size[1];
      tmp_motclus_data[b_loop_ub].f1->size[0] = 1;
      tmp_motclus_data[b_loop_ub].f1->size[1] =
          b_motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->size[1];
      emxEnsureCapacity_real_T(tmp_motclus_data[b_loop_ub].f1, i1);
      loop_ub = b_motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        tmp_motclus_data[b_loop_ub].f1->data[i1] =
            b_motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->data[i1];
      }
    }
    /* 'gkmPWMlasso4:339' motclus = newMotclus; */
    /* 'gkmPWMlasso4:341' correlation = corrcoef(cfile2, BB*OLS); */
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
                  &f_data[0], (blasint)BB->size[0], &OLS_data[0],
                  (blasint)OLS->size[0], 0.0, &indc_data[0],
                  (blasint)BB->size[0]);
    }
    /* 'gkmPWMlasso4:342' gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus,
     * sprintf("%s_%d_%d_%d", filename, int32(l_svm2), int32(k_svm2), int32(d)),
     * memefile,num,minL, minInfo, correlation(1,2)); */
    nfrac = maximum(OLS);
    c = maximum(normvec);
    c_sprintf(varargin_1, (int)rt_roundd(varargin_6),
              (int)rt_roundd(varargin_7), text);
    loop_ub = OLS->size[0];
    for (i = 0; i < loop_ub; i++) {
      OLS_data[i] /= nfrac;
    }
    i = BY->size[0];
    BY->size[0] = negvec->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    loop_ub = negvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      BY_data[i] = Pweig_data[(int)negvec_data[i] - 1];
    }
    emxFree_real_T(&Pweig);
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] /= c;
    }
    corrcoef(cfile2, indc, dv);
    gettopmotifs(OLS, BY, normvec, tmp_motclus, text, varargin_2,
                 size_tmp_idx_1, varargin_3, varargin_4, dv[2]);
  } else {
    /* 'gkmPWMlasso4:344' else */
    /* 'gkmPWMlasso4:346' fprintf('Running LASSO (1)\n'); */
    printf("Running LASSO (1)\n");
    fflush(stdout);
    /*  weigmat = lasso_cvmat(B, cfile2,'DFmax', d,'Standardize', false); */
    /* 'gkmPWMlasso4:348' [weigmat, ~] = lasso_cvmat(B, cfile2, d, false, 100);
     */
    i = comb->size[0] * comb->size[1];
    comb->size[0] = B->size[0];
    comb->size[1] = B->size[1];
    emxEnsureCapacity_real_T(comb, i);
    GCmat_data = comb->data;
    loop_ub = B->size[0] * B->size[1] - 1;
    for (i = 0; i <= loop_ub; i++) {
      GCmat_data[i] = B_data[i];
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
    b_lasso_cvmat(comb, diffc, varargin_10, weigmat, expl_temp_data, lk_size,
                  b_expl_temp_data, expl_temp_size, &nfrac, a__3_DF_data,
                  a__3_DF_size, c_expl_temp_data, b_expl_temp_size);
    xc_data = weigmat->data;
    /* 'gkmPWMlasso4:350' f = find(weigmat(:,1)~=0); */
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
    i = negvec->size[0];
    negvec->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(negvec, i);
    negvec_data = negvec->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      negvec_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso4:351' cfile3 = cfile2 -
     * B(:,f)*((B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2)); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, comb);
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (negvec->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)negvec->size[0], (blasint)1, (blasint)B->size[0],
                  1.0, &GCmat_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)negvec->size[0]);
    }
    mpower(comb, comb2);
    GCmat_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (indc->size[0] == 0)) {
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
                  1.0, &GCmat_data[0], (blasint)comb2->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &BY_data[0],
                  (blasint)comb2->size[0]);
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (negvec->size[0] == 0) || (BY->size[0] == 0)) {
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
                  (blasint)B->size[0], (blasint)1, (blasint)negvec->size[0],
                  1.0, &GCmat_data[0], (blasint)B->size[0], &BY_data[0],
                  (blasint)BY->size[0], 0.0, &loc_data[0], (blasint)B->size[0]);
    }
    /* 'gkmPWMlasso4:353' fprintf('Running LASSO (2)\n'); */
    printf("Running LASSO (2)\n");
    fflush(stdout);
    /*  weigmat2 = lasso_cvmat(B, cfile3,'DFmax', d,'Standardize', false); */
    /* 'gkmPWMlasso4:355' [weigmat2, ~] = lasso_cvmat(B, cfile3, d, false, 100);
     */
    emxInit_real_T(&b_f, 2);
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
      GCmat_data = comb->data;
      loop_ub = B->size[0] * B->size[1] - 1;
      for (i = 0; i <= loop_ub; i++) {
        GCmat_data[i] = B_data[i];
      }
      b_lasso_cvmat(comb, loc, varargin_10, b_f, d_expl_temp_data, lk_size,
                    e_expl_temp_data, expl_temp_size, &nfrac, a__3_DF_data,
                    a__3_DF_size, f_expl_temp_data, b_expl_temp_size);
      f_data = b_f->data;
      loc_data = loc->data;
    } else {
      c_binary_expand_op(B, cfile2, loc, varargin_10, b_f, d_expl_temp_data,
                         lk_size, e_expl_temp_data, expl_temp_size,
                         a__3_DF_data, a__3_DF_size, f_expl_temp_data,
                         b_expl_temp_size);
      f_data = b_f->data;
    }
    /* 'gkmPWMlasso4:357' weigmat = abs(weigmat(:,1)) + abs(weigmat2(:,1)); */
    /* 'gkmPWMlasso4:358' f = find(weigmat~=0); */
    nx = weigmat->size[0] - 1;
    i = indc->size[0];
    indc->size[0] = weigmat->size[0];
    emxEnsureCapacity_real_T(indc, i);
    indc_data = indc->data;
    for (k = 0; k <= nx; k++) {
      indc_data[k] = fabs(xc_data[k]);
    }
    emxFree_real_T(&weigmat);
    nx = b_f->size[0] - 1;
    i = BY->size[0];
    BY->size[0] = b_f->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    for (k = 0; k <= nx; k++) {
      BY_data[k] = fabs(f_data[k]);
    }
    emxFree_real_T(&b_f);
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
    i = negvec->size[0];
    negvec->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(negvec, i);
    negvec_data = negvec->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      negvec_data[i] = matches_data[i];
    }
    /* 'gkmPWMlasso4:359' motclus2 = clus_simmat_eig(corrcoef(B(:,f)),0.6); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    b_corrcoef(b_GCmat, comb);
    b_clus_simmat_eig(comb, tmp_motclus);
    tmp_motclus_data = tmp_motclus->data;
    /* 'gkmPWMlasso4:360' if length(length(motclus2)) ~= length(f) */
    if (1 != negvec->size[0]) {
      /* 'gkmPWMlasso4:361' f2 = f; */
      i = loc->size[0];
      loc->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(loc, i);
      loc_data = loc->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        loc_data[i] = negvec_data[i];
      }
      /* 'gkmPWMlasso4:362' f = zeros(length(motclus2),1); */
      i = negvec->size[0];
      negvec->size[0] = tmp_motclus->size[0];
      emxEnsureCapacity_real_T(negvec, i);
      negvec_data = negvec->data;
      loop_ub = tmp_motclus->size[0];
      for (i = 0; i < loop_ub; i++) {
        negvec_data[i] = 0.0;
      }
      /* 'gkmPWMlasso4:363' for i = 1:length(motclus2) */
      i = tmp_motclus->size[0];
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso4:364' if length(motclus2{i}) > 1 */
        if (tmp_motclus_data[b_i].f1->size[1] > 1) {
          /* 'gkmPWMlasso4:365' [~,f3] = max(abs(Z(f2(motclus2{i})))); */
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
          d_maximum(BY, &nfrac, &match_idx);
          /* 'gkmPWMlasso4:366' f(i) = f2(motclus2{i}(f3)); */
          negvec_data[b_i] =
              loc_data[(int)tmp_motclus_data[b_i].f1->data[match_idx - 1] - 1];
        } else {
          /* 'gkmPWMlasso4:367' else */
          /* 'gkmPWMlasso4:368' f(i) = f2(motclus2{i}(1)); */
          negvec_data[b_i] =
              loc_data[(int)tmp_motclus_data[b_i].f1->data[0] - 1];
        }
      }
    }
    /* 'gkmPWMlasso4:372' [a,b] = sort(abs(Z(f)), 'descend'); */
    nx = negvec->size[0];
    i = indc->size[0];
    indc->size[0] = negvec->size[0];
    emxEnsureCapacity_real_T(indc, i);
    indc_data = indc->data;
    for (k = 0; k < nx; k++) {
      indc_data[k] = fabs(Z_data[(int)negvec_data[k] - 1]);
    }
    sort(indc, idx);
    matches_data = idx->data;
    /* 'gkmPWMlasso4:373' f = f(b); */
    i = BY->size[0];
    BY->size[0] = idx->size[0];
    emxEnsureCapacity_real_T(BY, i);
    BY_data = BY->data;
    loop_ub = idx->size[0];
    for (i = 0; i < loop_ub; i++) {
      BY_data[i] = negvec_data[matches_data[i] - 1];
    }
    i = negvec->size[0];
    negvec->size[0] = BY->size[0];
    emxEnsureCapacity_real_T(negvec, i);
    negvec_data = negvec->data;
    loop_ub = BY->size[0];
    for (i = 0; i < loop_ub; i++) {
      negvec_data[i] = BY_data[i];
    }
    /* 'gkmPWMlasso4:374' F = length(f); */
    nfrac = negvec->size[0];
    /* 'gkmPWMlasso4:375' fprintf('Selecting top motifs\n') */
    printf("Selecting top motifs\n");
    fflush(stdout);
    /* 'gkmPWMlasso4:376' ind = true; */
    /* 'gkmPWMlasso4:377' OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, comb);
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (negvec->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = indc->size[0];
      indc->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        indc_data[i] = 0.0;
      }
    } else {
      i = indc->size[0];
      indc->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(indc, i);
      indc_data = indc->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)negvec->size[0], (blasint)1, (blasint)B->size[0],
                  1.0, &GCmat_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &indc_data[0],
                  (blasint)negvec->size[0]);
    }
    mpower(comb, comb2);
    GCmat_data = comb2->data;
    if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
        (indc->size[0] == 0)) {
      i = OLS->size[0];
      OLS->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      loop_ub = comb2->size[0];
      for (i = 0; i < loop_ub; i++) {
        OLS_data[i] = 0.0;
      }
    } else {
      i = OLS->size[0];
      OLS->size[0] = comb2->size[0];
      emxEnsureCapacity_real_T(OLS, i);
      OLS_data = OLS->data;
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  (blasint)comb2->size[0], (blasint)1, (blasint)comb2->size[1],
                  1.0, &GCmat_data[0], (blasint)comb2->size[0], &indc_data[0],
                  (blasint)indc->size[0], 0.0, &OLS_data[0],
                  (blasint)comb2->size[0]);
    }
    /* 'gkmPWMlasso4:378' Pweig = Z(f); */
    i = diffc->size[0];
    diffc->size[0] = negvec->size[0];
    emxEnsureCapacity_real_T(diffc, i);
    GCmat_data = diffc->data;
    loop_ub = negvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      GCmat_data[i] = Z_data[(int)negvec_data[i] - 1];
    }
    /* 'gkmPWMlasso4:379' while ind && F >= d */
    exitg1 = false;
    while ((!exitg1) && (nfrac >= varargin_10)) {
      /* 'gkmPWMlasso4:380' OLS =
       * (B(:,f(1:d)).'*B(:,f(1:d)))^-1*(B(:,f(1:d)).'*cfile2); */
      if (1.0 > varargin_10) {
        loop_ub = 0;
        b_loop_ub = 0;
        match_idx = 0;
      } else {
        loop_ub = (int)varargin_10;
        b_loop_ub = (int)varargin_10;
        match_idx = (int)varargin_10;
      }
      nx = B->size[0];
      i = b_GCmat->size[0] * b_GCmat->size[1];
      b_GCmat->size[0] = B->size[0];
      b_GCmat->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_GCmat, i);
      normvec_data = b_GCmat->data;
      for (i = 0; i < loop_ub; i++) {
        for (i1 = 0; i1 < nx; i1++) {
          normvec_data[i1 + b_GCmat->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
        }
      }
      loop_ub = B->size[0];
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = B->size[0];
      comb2->size[1] = b_loop_ub;
      emxEnsureCapacity_real_T(comb2, i);
      GCmat_data = comb2->data;
      for (i = 0; i < b_loop_ub; i++) {
        for (i1 = 0; i1 < loop_ub; i1++) {
          GCmat_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
        }
      }
      b_mtimes(b_GCmat, comb2, comb);
      loop_ub = B->size[0];
      i = comb2->size[0] * comb2->size[1];
      comb2->size[0] = B->size[0];
      comb2->size[1] = match_idx;
      emxEnsureCapacity_real_T(comb2, i);
      GCmat_data = comb2->data;
      for (i = 0; i < match_idx; i++) {
        for (i1 = 0; i1 < loop_ub; i1++) {
          GCmat_data[i1 + comb2->size[0] * i] =
              B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
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
                    (blasint)1, (blasint)B->size[0], 1.0, &GCmat_data[0],
                    (blasint)B->size[0], &cfile2_data[0],
                    (blasint)cfile2->size[0], 0.0, &indc_data[0],
                    (blasint)match_idx);
      }
      mpower(comb, comb2);
      GCmat_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (indc->size[0] == 0)) {
        i = OLS->size[0];
        OLS->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(OLS, i);
        OLS_data = OLS->data;
        loop_ub = comb2->size[0];
        for (i = 0; i < loop_ub; i++) {
          OLS_data[i] = 0.0;
        }
      } else {
        i = OLS->size[0];
        OLS->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(OLS, i);
        OLS_data = OLS->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &GCmat_data[0],
            (blasint)comb2->size[0], &indc_data[0], (blasint)indc->size[0], 0.0,
            &OLS_data[0], (blasint)comb2->size[0]);
      }
      /* 'gkmPWMlasso4:381' Pweig = Z(f(1:d)); */
      if (1.0 > varargin_10) {
        loop_ub = 0;
      } else {
        loop_ub = (int)varargin_10;
      }
      i = diffc->size[0];
      diffc->size[0] = loop_ub;
      emxEnsureCapacity_real_T(diffc, i);
      GCmat_data = diffc->data;
      for (i = 0; i < loop_ub; i++) {
        GCmat_data[i] = Z_data[(int)negvec_data[i] - 1];
      }
      /* 'gkmPWMlasso4:382' ff = []; */
      loc->size[0] = 0;
      /* 'gkmPWMlasso4:383' for i = 1:d */
      i = (int)varargin_10;
      for (b_i = 0; b_i < i; b_i++) {
        /* 'gkmPWMlasso4:384' if sign(OLS(i)) ~= sign(Pweig(i)) */
        GCpos1 = OLS_data[b_i];
        if (OLS_data[b_i] < 0.0) {
          GCpos1 = -1.0;
        } else if (OLS_data[b_i] > 0.0) {
          GCpos1 = 1.0;
        }
        c = Z_data[(int)negvec_data[b_i] - 1];
        if (c < 0.0) {
          c = -1.0;
        } else if (c > 0.0) {
          c = 1.0;
        }
        if (GCpos1 != c) {
          /* 'gkmPWMlasso4:385' ff = [ff;i]; */
          i1 = loc->size[0];
          i2 = loc->size[0];
          loc->size[0]++;
          emxEnsureCapacity_real_T(loc, i2);
          loc_data = loc->data;
          loc_data[i1] = (double)b_i + 1.0;
        }
      }
      /* 'gkmPWMlasso4:388' if length(ff) > 0 && F - length(ff) >= d */
      if ((loc->size[0] > 0) && (nfrac - (double)loc->size[0] >= varargin_10)) {
        /* 'gkmPWMlasso4:389' f(ff) = []; */
        i = idx->size[0];
        idx->size[0] = loc->size[0];
        emxEnsureCapacity_int32_T(idx, i);
        matches_data = idx->data;
        loop_ub = loc->size[0];
        for (i = 0; i < loop_ub; i++) {
          matches_data[i] = (int)loc_data[i];
        }
        e_nullAssignment(negvec, idx);
        negvec_data = negvec->data;
        /* 'gkmPWMlasso4:390' F = F - length(ff); */
        nfrac -= (double)loc->size[0];
      } else {
        /* 'gkmPWMlasso4:391' else */
        /* 'gkmPWMlasso4:392' ind = false; */
        /* 'gkmPWMlasso4:393' f = f(1:d); */
        i = negvec->size[0];
        if (1.0 > varargin_10) {
          negvec->size[0] = 0;
        } else {
          negvec->size[0] = (int)varargin_10;
        }
        emxEnsureCapacity_real_T(negvec, i);
        negvec_data = negvec->data;
        exitg1 = true;
      }
    }
    /*  if F < d */
    /*      OLS = (B(:,f).'*B(:,f))^-1*(B(:,f).'*cfile2); */
    /*      Pweig = Z(f); */
    /*  end */
    /* 'gkmPWMlasso4:400' BB = B(:,f); */
    loop_ub = B->size[0];
    i = BB->size[0] * BB->size[1];
    BB->size[0] = B->size[0];
    BB->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(BB, i);
    f_data = BB->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        f_data[i1 + BB->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    /* 'gkmPWMlasso4:401' BX = B(:,f)'*B(:,f); */
    loop_ub = B->size[0];
    i = b_GCmat->size[0] * b_GCmat->size[1];
    b_GCmat->size[0] = B->size[0];
    b_GCmat->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(b_GCmat, i);
    normvec_data = b_GCmat->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        normvec_data[i1 + b_GCmat->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    b_mtimes(b_GCmat, comb2, xc);
    xc_data = xc->data;
    /* 'gkmPWMlasso4:402' BY = B(:,f)'*cfile2; */
    loop_ub = B->size[0];
    i = comb2->size[0] * comb2->size[1];
    comb2->size[0] = B->size[0];
    comb2->size[1] = negvec->size[0];
    emxEnsureCapacity_real_T(comb2, i);
    GCmat_data = comb2->data;
    b_loop_ub = negvec->size[0];
    for (i = 0; i < b_loop_ub; i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        GCmat_data[i1 + comb2->size[0] * i] =
            B_data[i1 + B->size[0] * ((int)negvec_data[i] - 1)];
      }
    }
    if ((B->size[0] == 0) || (negvec->size[0] == 0) || (cfile2->size[0] == 0)) {
      i = BY->size[0];
      BY->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      loop_ub = negvec->size[0];
      for (i = 0; i < loop_ub; i++) {
        BY_data[i] = 0.0;
      }
    } else {
      i = BY->size[0];
      BY->size[0] = negvec->size[0];
      emxEnsureCapacity_real_T(BY, i);
      BY_data = BY->data;
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                  (blasint)negvec->size[0], (blasint)1, (blasint)B->size[0],
                  1.0, &GCmat_data[0], (blasint)B->size[0], &cfile2_data[0],
                  (blasint)cfile2->size[0], 0.0, &BY_data[0],
                  (blasint)negvec->size[0]);
    }
    /* 'gkmPWMlasso4:403' E = zeros(length(f),1); */
    i = normvec->size[0];
    normvec->size[0] = negvec->size[0];
    emxEnsureCapacity_real_T(normvec, i);
    normvec_data = normvec->data;
    loop_ub = negvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = 0.0;
    }
    /* 'gkmPWMlasso4:404' for i = 1:length(f) */
    i = negvec->size[0];
    if (0 <= negvec->size[0] - 1) {
      e_loop_ub = xc->size[0] * xc->size[1];
      f_loop_ub = BY->size[0];
      nxout = BY->size[0] - 1;
      g_loop_ub = BB->size[0] * BB->size[1];
    }
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:405' B = BB; */
      /* 'gkmPWMlasso4:406' BBX = BX; */
      /* 'gkmPWMlasso4:407' BBY = BY; */
      /* 'gkmPWMlasso4:408' B(:,i) = []; */
      /* 'gkmPWMlasso4:409' BBX(:,i) = []; */
      /* 'gkmPWMlasso4:410' BBX(i,:) = []; */
      /* 'gkmPWMlasso4:411' BBY(i) = []; */
      /* 'gkmPWMlasso4:412' res = cfile2-B*(BBX^-1*BBY); */
      i1 = comb->size[0] * comb->size[1];
      comb->size[0] = xc->size[0];
      comb->size[1] = xc->size[1];
      emxEnsureCapacity_real_T(comb, i1);
      GCmat_data = comb->data;
      for (i1 = 0; i1 < e_loop_ub; i1++) {
        GCmat_data[i1] = xc_data[i1];
      }
      d_nullAssignment(comb, b_i + 1);
      b_nullAssignment(comb, b_i + 1);
      mpower(comb, comb2);
      GCmat_data = comb2->data;
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
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (loc->size[0] == 0)) {
        i1 = indc->size[0];
        indc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(indc, i1);
        indc_data = indc->data;
        loop_ub = comb2->size[0];
        for (i1 = 0; i1 < loop_ub; i1++) {
          indc_data[i1] = 0.0;
        }
      } else {
        i1 = indc->size[0];
        indc->size[0] = comb2->size[0];
        emxEnsureCapacity_real_T(indc, i1);
        indc_data = indc->data;
        cblas_dgemm(
            CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)comb2->size[0],
            (blasint)1, (blasint)comb2->size[1], 1.0, &GCmat_data[0],
            (blasint)comb2->size[0], &loc_data[0], (blasint)loc->size[0], 0.0,
            &indc_data[0], (blasint)comb2->size[0]);
      }
      i1 = comb2->size[0] * comb2->size[1];
      comb2->size[0] = BB->size[0];
      comb2->size[1] = BB->size[1];
      emxEnsureCapacity_real_T(comb2, i1);
      GCmat_data = comb2->data;
      for (i1 = 0; i1 < g_loop_ub; i1++) {
        GCmat_data[i1] = f_data[i1];
      }
      d_nullAssignment(comb2, b_i + 1);
      GCmat_data = comb2->data;
      if ((comb2->size[0] == 0) || (comb2->size[1] == 0) ||
          (indc->size[0] == 0)) {
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
            (blasint)1, (blasint)comb2->size[1], 1.0, &GCmat_data[0],
            (blasint)comb2->size[0], &indc_data[0], (blasint)indc->size[0], 0.0,
            &loc_data[0], (blasint)comb2->size[0]);
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
      /* 'gkmPWMlasso4:413' E(i) = sqrt(res'*res); */
      if (loc->size[0] < 1) {
        c = 0.0;
      } else {
        c = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                       &loc_data[0], (blasint)1);
      }
      normvec_data[b_i] = sqrt(c);
    }
    /* 'gkmPWMlasso4:415' res = cfile2-BB*OLS; */
    if ((B->size[0] == 0) || (negvec->size[0] == 0) || (OLS->size[0] == 0)) {
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
                  (blasint)B->size[0], (blasint)1, (blasint)negvec->size[0],
                  1.0, &f_data[0], (blasint)B->size[0], &OLS_data[0],
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
    /* 'gkmPWMlasso4:416' EE = sqrt(res'*res); */
    if (loc->size[0] < 1) {
      c = 0.0;
    } else {
      c = cblas_ddot((blasint)loc->size[0], &loc_data[0], (blasint)1,
                     &loc_data[0], (blasint)1);
    }
    nfrac = sqrt(c);
    /* 'gkmPWMlasso4:417' E = (E-EE)/EE; */
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] = (normvec_data[i] - nfrac) / nfrac;
    }
    /* 'gkmPWMlasso4:418' for i = 1:length(motclus) */
    i = motclus->size[0];
    for (b_i = 0; b_i < i; b_i++) {
      /* 'gkmPWMlasso4:419' motclus{i} = indvec(motclus{i})'; */
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
    /* 'gkmPWMlasso4:422' mylen = length(f); */
    /* 'gkmPWMlasso4:423' newMotclus = cell(mylen, 1); */
    /* 'gkmPWMlasso4:424' for idx=1:mylen */
    i = negvec->size[0];
    i1 = c_motclus->size[0];
    c_motclus->size[0] = negvec->size[0];
    emxEnsureCapacity_cell_wrap_1(c_motclus, i1);
    b_motclus_data = c_motclus->data;
    for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
      /* 'gkmPWMlasso4:425' newMotclus{idx} = motclus{f(idx)}; */
      i1 = b_motclus_data[b_loop_ub].f1->size[0] *
           b_motclus_data[b_loop_ub].f1->size[1];
      b_motclus_data[b_loop_ub].f1->size[0] = 1;
      b_motclus_data[b_loop_ub].f1->size[1] =
          motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->size[1];
      emxEnsureCapacity_real_T(b_motclus_data[b_loop_ub].f1, i1);
      loop_ub = motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_motclus_data[b_loop_ub].f1->data[i1] =
            motclus_data[(int)negvec_data[b_loop_ub] - 1].f1->data[i1];
      }
    }
    /* 'gkmPWMlasso4:427' motclus = newMotclus; */
    /* 'gkmPWMlasso4:429' correlation = corrcoef(cfile2, BB*OLS); */
    if ((B->size[0] == 0) || (negvec->size[0] == 0) || (OLS->size[0] == 0)) {
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
                  (blasint)B->size[0], (blasint)1, (blasint)negvec->size[0],
                  1.0, &f_data[0], (blasint)B->size[0], &OLS_data[0],
                  (blasint)OLS->size[0], 0.0, &indc_data[0],
                  (blasint)B->size[0]);
    }
    /* 'gkmPWMlasso4:430' gettopmotifs(OLS/max(OLS), Pweig, E/max(E), motclus,
     * sprintf("%s_%d_%d", filename, int32(l_svm2),
     * int32(k_svm2)),memefile,num,minL, minInfo, correlation(1,2)); */
    nfrac = maximum(OLS);
    c = maximum(normvec);
    b_sprintf(varargin_1, (int)rt_roundd(varargin_6),
              (int)rt_roundd(varargin_7), text);
    loop_ub = OLS->size[0];
    for (i = 0; i < loop_ub; i++) {
      OLS_data[i] /= nfrac;
    }
    loop_ub = normvec->size[0];
    for (i = 0; i < loop_ub; i++) {
      normvec_data[i] /= c;
    }
    corrcoef(cfile2, indc, dv);
    gettopmotifs(OLS, diffc, normvec, c_motclus, text, varargin_2,
                 size_tmp_idx_1, varargin_3, varargin_4, dv[2]);
  }
  emxFree_real_T(&b_GCmat);
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
  /* 'gkmPWMlasso4:434' fprintf('Done\n'); */
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
