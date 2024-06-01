/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getgkmweights.c
 *
 * Code generation for function 'getgkmweights'
 *
 */

/* Include files */
#include "getgkmweights.h"
#include "blockedSummation.h"
#include "diff.h"
#include "fileManager.h"
#include "find.h"
#include "genIndex.h"
#include "getgkmcounts.h"
#include "getgkmweights_data.h"
#include "getgkmweights_emxutil.h"
#include "getgkmweights_initialize.h"
#include "getgkmweights_rtwutil.h"
#include "getgkmweights_types.h"
#include "mod.h"
#include "std.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* Function Declarations */
static void binary_expand_op(emxArray_real_T *cfile,
                             const emxArray_real_T *negvec, double y, double b);

/* Function Definitions */
static void binary_expand_op(emxArray_real_T *cfile,
                             const emxArray_real_T *negvec, double y, double b)
{
  emxArray_real_T *b_cfile;
  const double *negvec_data;
  double *b_cfile_data;
  double *cfile_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  negvec_data = negvec->data;
  cfile_data = cfile->data;
  emxInit_real_T(&b_cfile, 1);
  i = b_cfile->size[0];
  if (negvec->size[0] == 1) {
    b_cfile->size[0] = cfile->size[0];
  } else {
    b_cfile->size[0] = negvec->size[0];
  }
  emxEnsureCapacity_real_T(b_cfile, i);
  b_cfile_data = b_cfile->data;
  stride_0_0 = (cfile->size[0] != 1);
  stride_1_0 = (negvec->size[0] != 1);
  if (negvec->size[0] == 1) {
    loop_ub = cfile->size[0];
  } else {
    loop_ub = negvec->size[0];
  }
  for (i = 0; i < loop_ub; i++) {
    b_cfile_data[i] =
        cfile_data[i * stride_0_0] - negvec_data[i * stride_1_0] / y * b;
  }
  i = cfile->size[0];
  cfile->size[0] = b_cfile->size[0];
  emxEnsureCapacity_real_T(cfile, i);
  cfile_data = cfile->data;
  loop_ub = b_cfile->size[0];
  for (i = 0; i < loop_ub; i++) {
    cfile_data[i] = b_cfile_data[i];
  }
  emxFree_real_T(&b_cfile);
}

/*
 * function getgkmweights(varargin)
 */
void getgkmweights(const emxArray_char_T *varargin_1, double varargin_2,
                   double varargin_3, double varargin_4, double varargin_5,
                   double varargin_6)
{
  FILE *b_NULL;
  FILE *filestar;
  cell_wrap_10 *p_data;
  emxArray_boolean_T *b_dc2;
  emxArray_cell_wrap_10 *p;
  emxArray_char_T *S;
  emxArray_char_T *b_varargin_1;
  emxArray_char_T *charStr;
  emxArray_char_T *filename;
  emxArray_char_T *kmer;
  emxArray_char_T *s;
  emxArray_int32_T *r;
  emxArray_int32_T *r1;
  emxArray_real_T *b_comb;
  emxArray_real_T *c2;
  emxArray_real_T *cfile;
  emxArray_real_T *comb;
  emxArray_real_T *dc;
  emxArray_real_T *dc2;
  emxArray_real_T *negvec;
  emxArray_real_T *seqvec;
  emxArray_real_T *seqvec2;
  double mat[16];
  double mat2[16];
  double startvec[4];
  double lk_data[2];
  double GCmat_idx_0_tmp;
  double GCmat_idx_1;
  double GCmat_idx_2;
  double GCneg1;
  double c_tmp;
  double nfrac;
  double rcnum;
  double tmp;
  double validatedHoleFilling_idx_0;
  double xtmp;
  double *b_comb_data;
  double *c2_data;
  double *cfile_data;
  double *comb_data;
  double *dc2_data;
  double *dc_data;
  double *negvec_data;
  double *seqvec2_data;
  double *seqvec_data;
  int lk_size[2];
  unsigned int a;
  int b_i;
  int b_j1;
  int b_loop_ub;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int j2;
  int loop_ub;
  int m;
  int nd2;
  int *r2;
  int *r3;
  const char *varargin_1_data;
  signed char fileid;
  char *S_data;
  char *b_varargin_1_data;
  char *kmer_data;
  char *s_data;
  bool autoflush;
  bool *b_dc2_data;
  if (!isInitialized_getgkmweights) {
    getgkmweights_initialize();
  }
  varargin_1_data = varargin_1->data;
  /*  getgkmweights saves the gapped kmers weights for a given (l,k) */
  /*  */
  /*      getgkmweights(fileprefix, l, k ...) */
  /*  */
  /*      Works conveniently with the output of the gkmSVM R package or the
   * lsgkm */
  /*      package (https://github.com/Dongwon-Lee/lsgkm).  Can be leveraged to
   */
  /*      extract motifs from other sequence based models as long as you have
   * scores */
  /*      paired with sequences.  The sequences should be in fasta format with
   * the */
  /*      suffix *_svseq.fa.  The scores should be in a 2 column tab delimited
   */
  /*      file with the suffix *svalpha.out with the same prefix as the
   * *svseq.fa */
  /*      file.  The first column containing the labels of each sequence and the
   */
  /*      scores in the second. */
  /* 'getgkmweights:15' fn = varargin{1}; */
  /* 'getgkmweights:16' l_svm = varargin{2}; */
  /* 'getgkmweights:17' k_svm = varargin{3}; */
  /* 'getgkmweights:18' RC = varargin{4}; */
  /* 'getgkmweights:19' nfrac = varargin{5}; */
  /* 'getgkmweights:20' nfracLim = varargin{6}; */
  /* 'getgkmweights:21' lk = 1; */
  lk_size[0] = 1;
  lk_size[1] = 1;
  lk_data[0] = 1.0;
  /* 'getgkmweights:23' if nfrac ~= 1 */
  if (varargin_5 != 1.0) {
    /* 'getgkmweights:24' lk = [l_svm k_svm]; */
    lk_size[0] = 1;
    lk_size[1] = 2;
    lk_data[0] = varargin_2;
    lk_data[1] = varargin_3;
  }
  emxInit_real_T(&comb, 2);
  emxInit_real_T(&cfile, 1);
  emxInit_real_T(&negvec, 1);
  emxInit_real_T(&c2, 2);
  emxInit_real_T(&seqvec2, 2);
  /* 'getgkmweights:27' [comb,comb2,diffc,indc,xc,rcnum] =
   * genIndex(l_svm,k_svm,nfrac); */
  genIndex(varargin_2, varargin_3, varargin_5, comb, c2, negvec, cfile, seqvec2,
           &rcnum);
  comb_data = comb->data;
  /* 'getgkmweights:28' if nfracLim && length(comb)*4^k_svm > 5*10^5 */
  if (varargin_6 != 0.0) {
    if ((comb->size[0] == 0) || (comb->size[1] == 0)) {
      nd2 = 0;
    } else {
      m = comb->size[0];
      nd2 = comb->size[1];
      if (m >= nd2) {
        nd2 = m;
      }
    }
    nfrac = pow(4.0, varargin_3);
    if ((double)nd2 * nfrac > 500000.0) {
      /* 'getgkmweights:29' nfrac = round(5*10^7/4^k_svm/numel(comb)*k_svm)/100;
       */
      nfrac = rt_roundd(5.0E+7 / nfrac /
                        (double)(comb->size[0] * comb->size[1]) * varargin_3) /
              100.0;
    
      // /* 'getgkmweights:30' fprintf('Combination of (l,k) yields too many gapped
      //  * kmers.  Using %f of the total gapped kmers\n', nfrac); */
      // printf("Combination of (l,k) yields too many gapped kmers.  Using %f of "
      //        "the total gapped kmers\n",
      //        nfrac);
      // fflush(stdout);
    
      /* 'getgkmweights:31' lk = ([l_svm k_svm]); */
      lk_size[0] = 1;
      lk_size[1] = 2;
      lk_data[0] = varargin_2;
      lk_data[1] = varargin_3;
      /* 'getgkmweights:32' [comb,comb2,diffc,indc,xc,rcnum] =
       * genIndex(l_svm,k_svm,nfrac); */
      genIndex(varargin_2, varargin_3, nfrac, comb, c2, negvec, cfile, seqvec2,
               &rcnum);
      comb_data = comb->data;
        
      printf("WARNING: Using %d gapped kmers\n", (int)((double)(comb->size[0] * comb->size[1]) / varargin_3 * pow(4.0, varargin_3)));
      fflush(stdout); 
    }
  }
  /* 'getgkmweights:34' fprintf('Counting gapped kmers\n'); */
  printf("Counting gapped kmers\n");
  fflush(stdout);
  /* 'getgkmweights:35' [cfile, GCpos1, GCneg1,mat,mat2] = getgkmcounts(fn,
   * l_svm, k_svm, lk, RC, comb,rcnum); */
  getgkmcounts(varargin_1, varargin_2, varargin_3, lk_data, lk_size, varargin_4,
               comb, rcnum, cfile, &validatedHoleFilling_idx_0, &GCneg1, mat,
               mat2);
  cfile_data = cfile->data;
  /* 'getgkmweights:36' negvec = BGkmer(mat, GCneg1,comb,rcnum,l_svm,k_svm,RC);
   */
  /* 'BGkmer:2' len = numel(c)/k; */
  nfrac = (double)(comb->size[0] * comb->size[1]) / varargin_3;
  /* 'BGkmer:3' alen = len-rcnum; */
  /* 'BGkmer:4' negweights = zeros(len*4^k,1); */
  c_tmp = pow(4.0, varargin_3);
  m = (int)(nfrac * c_tmp);
  i = negvec->size[0];
  negvec->size[0] = m;
  emxEnsureCapacity_real_T(negvec, i);
  negvec_data = negvec->data;
  for (i = 0; i < m; i++) {
    negvec_data[i] = 0.0;
  }
  emxInit_real_T(&dc2, 2);
  emxInit_cell_wrap_10(&p);
  /* 'BGkmer:5' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  GCmat_idx_0_tmp = 0.5 - GCneg1 / 2.0;
  GCmat_idx_1 = GCneg1 / 2.0;
  GCmat_idx_2 = GCneg1 / 2.0;
  /* 'BGkmer:6' c2 = 0; */
  i = c2->size[0] * c2->size[1];
  c2->size[0] = 1;
  c2->size[1] = 1;
  emxEnsureCapacity_real_T(c2, i);
  c2_data = c2->data;
  c2_data[0] = 0.0;
  /*  To appease Coder */
  /* 'BGkmer:7' seqvec2 = 0; */
  i = seqvec2->size[0] * seqvec2->size[1];
  seqvec2->size[0] = 1;
  seqvec2->size[1] = 1;
  emxEnsureCapacity_real_T(seqvec2, i);
  seqvec2_data = seqvec2->data;
  seqvec2_data[0] = 0.0;
  /*  To appease Coder */
  /* 'BGkmer:8' dc2 = 0; */
  i = dc2->size[0] * dc2->size[1];
  dc2->size[0] = 1;
  dc2->size[1] = 1;
  emxEnsureCapacity_real_T(dc2, i);
  dc2_data = dc2->data;
  dc2_data[0] = 0.0;
  /*  To appease Coder */
  /* 'BGkmer:9' tmp = 0; */
  tmp = 0.0;
  /*  To appease Coder */
  /* 'BGkmer:10' p = cell(l,1); */
  i = (int)varargin_2;
  i1 = p->size[0];
  p->size[0] = (int)varargin_2;
  emxEnsureCapacity_cell_wrap_10(p, i1);
  p_data = p->data;
  /* 'BGkmer:11' p = coder.nullcopy(p); */
  /* 'BGkmer:12' p{1} = eye(4); */
  for (i1 = 0; i1 < 16; i1++) {
    p_data[0].f1[i1] = 0.0;
  }
  p_data[0].f1[0] = 1.0;
  p_data[0].f1[5] = 1.0;
  p_data[0].f1[10] = 1.0;
  p_data[0].f1[15] = 1.0;
  /* 'BGkmer:13' for i = 1:l-1 */
  i1 = (int)(varargin_2 - 1.0);
  for (b_i = 0; b_i < i1; b_i++) {
    /* 'BGkmer:14' p{i+1} = p{i}*mat; */
    for (i2 = 0; i2 < 4; i2++) {
      for (i3 = 0; i3 < 4; i3++) {
        i4 = i3 << 2;
        mat2[i2 + i4] = ((p_data[b_i].f1[i2] * mat[i4] +
                          p_data[b_i].f1[i2 + 4] * mat[i4 + 1]) +
                         p_data[b_i].f1[i2 + 8] * mat[i4 + 2]) +
                        p_data[b_i].f1[i2 + 12] * mat[i4 + 3];
      }
    }
    for (i2 = 0; i2 < 16; i2++) {
      p_data[(int)(((double)b_i + 1.0) + 1.0) - 1].f1[i2] = mat2[i2];
    }
  }
  emxInit_real_T(&seqvec, 2);
  /* 'BGkmer:16' seqvec = zeros(4^k, k); */
  i1 = seqvec->size[0] * seqvec->size[1];
  i2 = (int)c_tmp;
  seqvec->size[0] = (int)c_tmp;
  i3 = (int)varargin_3;
  seqvec->size[1] = (int)varargin_3;
  emxEnsureCapacity_real_T(seqvec, i1);
  seqvec_data = seqvec->data;
  m = (int)c_tmp * (int)varargin_3;
  for (i1 = 0; i1 < m; i1++) {
    seqvec_data[i1] = 0.0;
  }
  /* 'BGkmer:17' for i = 1:k */
  for (b_i = 0; b_i < i3; b_i++) {
    /* 'BGkmer:18' for j = 1:4^k */
    for (nd2 = 0; nd2 < i2; nd2++) {
      /* 'BGkmer:19' seqvec(j,i) = mod(floor((j-1)/4^(i-1)), 4)+1; */
      seqvec_data[nd2 + seqvec->size[0] * b_i] =
          b_mod(floor((((double)nd2 + 1.0) - 1.0) /
                      pow(4.0, ((double)b_i + 1.0) - 1.0)),
                4.0) +
          1.0;
    }
  }
  /* 'BGkmer:22' if RC */
  if (varargin_4 != 0.0) {
    /* 'BGkmer:23' seqvec2 = 5-fliplr(seqvec); */
    i1 = seqvec2->size[0] * seqvec2->size[1];
    seqvec2->size[0] = seqvec->size[0];
    seqvec2->size[1] = seqvec->size[1];
    emxEnsureCapacity_real_T(seqvec2, i1);
    seqvec2_data = seqvec2->data;
    m = seqvec->size[0] * seqvec->size[1];
    for (i1 = 0; i1 < m; i1++) {
      seqvec2_data[i1] = seqvec_data[i1];
    }
    m = seqvec->size[0];
    nd2 = seqvec->size[1] >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (seqvec->size[1] - b_j1) - 1;
      for (b_i = 0; b_i < m; b_i++) {
        xtmp = seqvec2_data[b_i + seqvec2->size[0] * b_j1];
        seqvec2_data[b_i + seqvec2->size[0] * b_j1] =
            seqvec2_data[b_i + seqvec2->size[0] * j2];
        seqvec2_data[b_i + seqvec2->size[0] * j2] = xtmp;
      }
    }
    m = seqvec2->size[0] * seqvec2->size[1];
    for (i1 = 0; i1 < m; i1++) {
      seqvec2_data[i1] = 5.0 - seqvec2_data[i1];
    }
    /* 'BGkmer:24' c2 = l+1-fliplr(c); */
    i1 = c2->size[0] * c2->size[1];
    c2->size[0] = comb->size[0];
    c2->size[1] = comb->size[1];
    emxEnsureCapacity_real_T(c2, i1);
    c2_data = c2->data;
    m = comb->size[0] * comb->size[1];
    for (i1 = 0; i1 < m; i1++) {
      c2_data[i1] = comb_data[i1];
    }
    m = comb->size[0];
    nd2 = comb->size[1] >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (comb->size[1] - b_j1) - 1;
      for (b_i = 0; b_i < m; b_i++) {
        xtmp = c2_data[b_i + c2->size[0] * b_j1];
        c2_data[b_i + c2->size[0] * b_j1] = c2_data[b_i + c2->size[0] * j2];
        c2_data[b_i + c2->size[0] * j2] = xtmp;
      }
    }
    m = c2->size[0] * c2->size[1];
    for (i1 = 0; i1 < m; i1++) {
      c2_data[i1] = (varargin_2 + 1.0) - c2_data[i1];
    }
  }
  /* 'BGkmer:26' a = 1; */
  a = 1U;
  /* 'BGkmer:27' for i = 1:len */
  i1 = (int)nfrac;
  emxInit_real_T(&dc, 2);
  dc_data = dc->data;
  emxInit_real_T(&b_comb, 2);
  for (b_i = 0; b_i < i1; b_i++) {
    /* 'BGkmer:28' dc = diff(c(i,:)); */
    m = comb->size[1];
    i4 = b_comb->size[0] * b_comb->size[1];
    b_comb->size[0] = 1;
    b_comb->size[1] = comb->size[1];
    emxEnsureCapacity_real_T(b_comb, i4);
    b_comb_data = b_comb->data;
    for (i4 = 0; i4 < m; i4++) {
      b_comb_data[i4] = comb_data[b_i + comb->size[0] * i4];
    }
    diff(b_comb, dc);
    dc_data = dc->data;
    /* 'BGkmer:29' startvec = GCmat*p{c(i,1)}; */
    for (i4 = 0; i4 < 16; i4++) {
      mat2[i4] = p_data[(int)comb_data[b_i] - 1].f1[i4];
    }
    for (i4 = 0; i4 < 4; i4++) {
      b_j1 = i4 << 2;
      startvec[i4] =
          ((GCmat_idx_0_tmp * mat2[b_j1] + GCmat_idx_1 * mat2[b_j1 + 1]) +
           GCmat_idx_2 * mat2[b_j1 + 2]) +
          GCmat_idx_0_tmp * mat2[b_j1 + 3];
    }
    /* 'BGkmer:30' if RC */
    if (varargin_4 != 0.0) {
      /* 'BGkmer:31' dc2 = diff(c2(i,:)); */
      m = c2->size[1];
      i4 = b_comb->size[0] * b_comb->size[1];
      b_comb->size[0] = 1;
      b_comb->size[1] = c2->size[1];
      emxEnsureCapacity_real_T(b_comb, i4);
      b_comb_data = b_comb->data;
      for (i4 = 0; i4 < m; i4++) {
        b_comb_data[i4] = c2_data[b_i + c2->size[0] * i4];
      }
      diff(b_comb, dc2);
      dc2_data = dc2->data;
      /* 'BGkmer:32' startvec2 = GCmat*p{c2(i,1)}; */
    }
    /* 'BGkmer:34' for ii = 1:4^k */
    i4 = (int)pow(4.0, varargin_3);
    for (nd2 = 0; nd2 < i4; nd2++) {
      /* 'BGkmer:35' negweights(a) = startvec(seqvec(ii,1)); */
      m = (int)(a + nd2) - 1;
      negvec_data[m] = startvec[(int)seqvec_data[nd2] - 1];
      /* 'BGkmer:36' if RC */
      if (varargin_4 != 0.0) {
        /* 'BGkmer:37' tmp = startvec(seqvec2(ii,1)); */
        tmp = startvec[(int)seqvec2_data[nd2] - 1];
      }
      /* 'BGkmer:39' for iii = 1:k-1 */
      b_j1 = (int)(varargin_3 - 1.0);
      for (j2 = 0; j2 < b_j1; j2++) {
        /* 'BGkmer:40' matt = p{dc(iii)+1}; */
        /* 'BGkmer:41' negweights(a) = negweights(a)*matt(seqvec(ii,iii),
         * seqvec(ii,iii+1)); */
        negvec_data[m] *=
            p_data[(int)(dc_data[j2] + 1.0) - 1]
                .f1[((int)seqvec_data[nd2 + seqvec->size[0] * j2] +
                     (((int)seqvec_data[nd2 + seqvec->size[0] * (j2 + 1)] - 1)
                      << 2)) -
                    1];
        /* 'BGkmer:42' if RC */
        if (varargin_4 != 0.0) {
          /* 'BGkmer:43' matt2 = p{dc2(iii)+1}; */
          /* 'BGkmer:44' tmp = tmp*matt2(seqvec2(ii,iii), seqvec2(ii,iii+1)); */
          tmp *=
              p_data[(int)(dc2_data[j2] + 1.0) - 1].f1
                  [((int)seqvec2_data[nd2 + seqvec2->size[0] * j2] +
                    (((int)seqvec2_data[nd2 + seqvec2->size[0] * (j2 + 1)] - 1)
                     << 2)) -
                   1];
        }
      }
      /* 'BGkmer:47' if RC */
      if (varargin_4 != 0.0) {
        /* 'BGkmer:48' negweights(a) = (negweights(a)+tmp); */
        negvec_data[m] += tmp;
      }
      /* 'BGkmer:50' a = a+1; */
    }
    a += i4;
  }
  emxFree_cell_wrap_10(&p);
  emxFree_real_T(&seqvec);
  emxFree_real_T(&seqvec2);
  emxFree_real_T(&c2);
  /* 'BGkmer:53' if rcnum > 0 && RC */
  if ((rcnum > 0.0) && (varargin_4 != 0.0)) {
    /* 'BGkmer:54' negweights(4^k*alen+1:end) =
     * negweights(4^k*alen+1:end)/sqrt(2); */
    nfrac = c_tmp * (nfrac - rcnum) + 1.0;
    if (nfrac > negvec->size[0]) {
      i1 = 1;
      i4 = -1;
      b_j1 = 0;
    } else {
      i1 = (int)nfrac;
      i4 = (int)nfrac - 2;
      b_j1 = negvec->size[0];
    }
    nd2 = (b_j1 - i4) - 1;
    b_j1 = b_comb->size[0] * b_comb->size[1];
    b_comb->size[0] = 1;
    b_comb->size[1] = nd2;
    emxEnsureCapacity_real_T(b_comb, b_j1);
    b_comb_data = b_comb->data;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      b_comb_data[b_j1] = negvec_data[(i1 + b_j1) - 1] / 1.4142135623730951;
    }
    m = b_comb->size[1];
    for (i1 = 0; i1 < m; i1++) {
      negvec_data[(i4 + i1) + 1] = b_comb_data[i1];
    }
  }
  emxFree_real_T(&b_comb);
  /* 'getgkmweights:37' cfile = cfile-negvec/sum(negvec)*sum(cfile); */
  nfrac = blockedSummation(negvec, negvec->size[0]);
  GCmat_idx_0_tmp = blockedSummation(cfile, cfile->size[0]);
  if (cfile->size[0] == negvec->size[0]) {
    m = cfile->size[0];
    for (i1 = 0; i1 < m; i1++) {
      cfile_data[i1] -= negvec_data[i1] / nfrac * GCmat_idx_0_tmp;
    }
  } else {
    binary_expand_op(cfile, negvec, nfrac, GCmat_idx_0_tmp);
    cfile_data = cfile->data;
  }
  emxFree_real_T(&negvec);
  /* 'getgkmweights:38' cfile = cfile/std(cfile); */
  nfrac = b_std(cfile);
  m = cfile->size[0];
  for (i1 = 0; i1 < m; i1++) {
    cfile_data[i1] /= nfrac;
  }
  emxInit_char_T(&s, 2);
  s_data = s->data;
  emxInit_char_T(&S, 2);
  /* 'getgkmweights:39' s = ''; */
  s->size[0] = 1;
  s->size[1] = 0;
  /* 'getgkmweights:40' for i = 1:k_svm */
  for (b_i = 0; b_i < i3; b_i++) {
    /* 'getgkmweights:41' s = [s 'A']; */
    i1 = s->size[1];
    i4 = s->size[0] * s->size[1];
    s->size[1]++;
    emxEnsureCapacity_char_T(s, i4);
    s_data = s->data;
    s_data[i1] = 'A';
  }
  emxInit_char_T(&kmer, 2);
  /* 'getgkmweights:43' kmer = repmat(s,4^k_svm,1); */
  i1 = kmer->size[0] * kmer->size[1];
  kmer->size[0] = (int)c_tmp;
  kmer->size[1] = s->size[1];
  emxEnsureCapacity_char_T(kmer, i1);
  kmer_data = kmer->data;
  nd2 = s->size[1];
  for (m = 0; m < nd2; m++) {
    j2 = m * (int)c_tmp;
    for (b_j1 = 0; b_j1 < i2; b_j1++) {
      kmer_data[j2 + b_j1] = 'A';
    }
  }
  /* 'getgkmweights:44' vec = fliplr(0:4^k_svm-1); */
  if (c_tmp - 1.0 < 0.0) {
    dc->size[0] = 1;
    dc->size[1] = 0;
  } else {
    i1 = dc->size[0] * dc->size[1];
    dc->size[0] = 1;
    m = (int)floor(c_tmp - 1.0);
    dc->size[1] = m + 1;
    emxEnsureCapacity_real_T(dc, i1);
    dc_data = dc->data;
    for (i1 = 0; i1 <= m; i1++) {
      dc_data[i1] = i1;
    }
  }
  m = dc->size[1] - 1;
  nd2 = dc->size[1] >> 1;
  for (b_j1 = 0; b_j1 < nd2; b_j1++) {
    j2 = m - b_j1;
    xtmp = dc_data[b_j1];
    dc_data[b_j1] = dc_data[j2];
    dc_data[j2] = xtmp;
  }
  /* 'getgkmweights:45' for i = 1:k_svm */
  emxInit_int32_T(&r, 1);
  emxInit_int32_T(&r1, 2);
  emxInit_boolean_T(&b_dc2, 2);
  for (b_i = 0; b_i < i3; b_i++) {
    /* 'getgkmweights:46' vec2 = mod(floor(vec/4^(k_svm-i)),4); */
    nfrac = pow(4.0, varargin_3 - ((double)b_i + 1.0));
    i1 = dc2->size[0] * dc2->size[1];
    dc2->size[0] = 1;
    dc2->size[1] = dc->size[1];
    emxEnsureCapacity_real_T(dc2, i1);
    dc2_data = dc2->data;
    m = dc->size[1];
    for (i1 = 0; i1 < m; i1++) {
      dc2_data[i1] = dc_data[i1] / nfrac;
    }
    nd2 = dc2->size[1];
    for (m = 0; m < nd2; m++) {
      dc2_data[m] = floor(dc2_data[m]);
    }
    i1 = dc2->size[0] * dc2->size[1];
    dc2->size[0] = 1;
    emxEnsureCapacity_real_T(dc2, i1);
    dc2_data = dc2->data;
    m = dc2->size[1] - 1;
    for (i1 = 0; i1 <= m; i1++) {
      nfrac = dc2_data[i1];
      if (nfrac == 0.0) {
        nfrac = 0.0;
      } else {
        nfrac = fmod(nfrac, 4.0);
        if (nfrac == 0.0) {
          nfrac = 0.0;
        }
      }
      dc2_data[i1] = nfrac;
    }
    /* 'getgkmweights:47' f = find(vec2==1); */
    /* 'getgkmweights:48' kmer(f,i) = 'C'; */
    i1 = b_dc2->size[0] * b_dc2->size[1];
    b_dc2->size[0] = 1;
    b_dc2->size[1] = dc2->size[1];
    emxEnsureCapacity_boolean_T(b_dc2, i1);
    b_dc2_data = b_dc2->data;
    m = dc2->size[1];
    for (i1 = 0; i1 < m; i1++) {
      b_dc2_data[i1] = (dc2_data[i1] == 1.0);
    }
    eml_find(b_dc2, r1);
    r2 = r1->data;
    i1 = r->size[0];
    r->size[0] = r1->size[1];
    emxEnsureCapacity_int32_T(r, i1);
    r3 = r->data;
    m = r1->size[1];
    for (i1 = 0; i1 < m; i1++) {
      r3[i1] = r2[i1];
    }
    m = r->size[0];
    for (i1 = 0; i1 < m; i1++) {
      kmer_data[(r3[i1] + kmer->size[0] * b_i) - 1] = 'C';
    }
    /* 'getgkmweights:49' f = find(vec2==2); */
    /* 'getgkmweights:50' kmer(f,i) = 'G'; */
    i1 = b_dc2->size[0] * b_dc2->size[1];
    b_dc2->size[0] = 1;
    b_dc2->size[1] = dc2->size[1];
    emxEnsureCapacity_boolean_T(b_dc2, i1);
    b_dc2_data = b_dc2->data;
    m = dc2->size[1];
    for (i1 = 0; i1 < m; i1++) {
      b_dc2_data[i1] = (dc2_data[i1] == 2.0);
    }
    eml_find(b_dc2, r1);
    r2 = r1->data;
    i1 = r->size[0];
    r->size[0] = r1->size[1];
    emxEnsureCapacity_int32_T(r, i1);
    r3 = r->data;
    m = r1->size[1];
    for (i1 = 0; i1 < m; i1++) {
      r3[i1] = r2[i1];
    }
    m = r->size[0];
    for (i1 = 0; i1 < m; i1++) {
      kmer_data[(r3[i1] + kmer->size[0] * b_i) - 1] = 'G';
    }
    /* 'getgkmweights:51' f = find(vec2==3); */
    /* 'getgkmweights:52' kmer(f,i) = 'T'; */
    i1 = b_dc2->size[0] * b_dc2->size[1];
    b_dc2->size[0] = 1;
    b_dc2->size[1] = dc2->size[1];
    emxEnsureCapacity_boolean_T(b_dc2, i1);
    b_dc2_data = b_dc2->data;
    m = dc2->size[1];
    for (i1 = 0; i1 < m; i1++) {
      b_dc2_data[i1] = (dc2_data[i1] == 3.0);
    }
    eml_find(b_dc2, r1);
    r2 = r1->data;
    i1 = r->size[0];
    r->size[0] = r1->size[1];
    emxEnsureCapacity_int32_T(r, i1);
    r3 = r->data;
    m = r1->size[1];
    for (i1 = 0; i1 < m; i1++) {
      r3[i1] = r2[i1];
    }
    m = r->size[0];
    for (i1 = 0; i1 < m; i1++) {
      kmer_data[(r3[i1] + kmer->size[0] * b_i) - 1] = 'T';
    }
  }
  emxFree_boolean_T(&b_dc2);
  emxFree_int32_T(&r1);
  emxFree_real_T(&dc);
  emxFree_real_T(&dc2);
  emxFree_int32_T(&r);
  /* 'getgkmweights:54' s = ''; */
  s->size[0] = 1;
  s->size[1] = 0;
  /* 'getgkmweights:55' for i = 1:l_svm */
  for (b_i = 0; b_i < i; b_i++) {
    /* 'getgkmweights:56' s = [s '-']; */
    i1 = s->size[1];
    i2 = s->size[0] * s->size[1];
    s->size[1]++;
    emxEnsureCapacity_char_T(s, i2);
    s_data = s->data;
    s_data[i1] = '-';
  }
  emxInit_char_T(&b_varargin_1, 2);
  /* 'getgkmweights:58' a = 1; */
  a = 1U;
  /* 'getgkmweights:59' fid = fopen(sprintf('%s_%d_%d_gkmweights.out', fn,
   * int32(l_svm), int32(k_svm)), 'w'); */
  nd2 = (int)rt_roundd(varargin_2);
  j2 = (int)rt_roundd(varargin_3);
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  b_varargin_1_data = b_varargin_1->data;
  m = varargin_1->size[1];
  for (i = 0; i < m; i++) {
    b_varargin_1_data[i] = varargin_1_data[i];
  }
  b_varargin_1_data[varargin_1->size[1]] = '\x00';
  i = S->size[0] * S->size[1];
  S->size[0] = 1;
  S->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(S, i);
  S_data = S->data;
  m = varargin_1->size[1];
  for (i = 0; i < m; i++) {
    S_data[i] = varargin_1_data[i];
  }
  emxInit_char_T(&charStr, 2);
  S_data[varargin_1->size[1]] = '\x00';
  m = snprintf(NULL, 0, "%s_%d_%d_gkmweights.out", &S_data[0], nd2, j2);
  i = charStr->size[0] * charStr->size[1];
  charStr->size[0] = 1;
  charStr->size[1] = m + 1;
  emxEnsureCapacity_char_T(charStr, i);
  S_data = charStr->data;
  snprintf(&S_data[0], (size_t)(m + 1), "%s_%d_%d_gkmweights.out",
           &b_varargin_1_data[0], nd2, j2);
  i = charStr->size[0] * charStr->size[1];
  if (1 > m) {
    charStr->size[1] = 0;
  } else {
    charStr->size[1] = m;
  }
  emxEnsureCapacity_char_T(charStr, i);
  fileid = cfopen(charStr, "wb");
  /* 'getgkmweights:60' for i = 1:numel(comb)/k_svm */
  i = (int)((double)(comb->size[0] * comb->size[1]) / varargin_3);
  emxFree_char_T(&charStr);
  if (0 <= i - 1) {
    i5 = (int)pow(4.0, varargin_3);
    if (0 <= i5 - 1) {
      loop_ub = s->size[1];
      b_loop_ub = kmer->size[1];
      b_NULL = NULL;
    }
  }
  for (b_i = 0; b_i < i; b_i++) {
    /* 'getgkmweights:61' for j = 1:4^k_svm */
    for (nd2 = 0; nd2 < i5; nd2++) {
      /* 'getgkmweights:62' S = s; */
      i1 = S->size[0] * S->size[1];
      S->size[0] = 1;
      S->size[1] = s->size[1];
      emxEnsureCapacity_char_T(S, i1);
      S_data = S->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        S_data[i1] = s_data[i1];
      }
      /* 'getgkmweights:63' S(comb(i,:)) = kmer(j,:); */
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        S_data[(int)comb_data[b_i + comb->size[0] * i1] - 1] =
            kmer_data[nd2 + kmer->size[0] * i1];
      }
      /* 'getgkmweights:64' fprintf(fid, '%s\t%0.5f\n', S, cfile(a)); */
      i1 = b_varargin_1->size[0] * b_varargin_1->size[1];
      b_varargin_1->size[0] = 1;
      b_varargin_1->size[1] = S->size[1] + 1;
      emxEnsureCapacity_char_T(b_varargin_1, i1);
      b_varargin_1_data = b_varargin_1->data;
      m = S->size[1];
      for (i1 = 0; i1 < m; i1++) {
        b_varargin_1_data[i1] = S_data[i1];
      }
      b_varargin_1_data[S->size[1]] = '\x00';
      getfilestar(fileid, &filestar, &autoflush);
      if (!(filestar == b_NULL)) {
        fprintf(filestar, "%s\t%0.5f\n", &b_varargin_1_data[0],
                cfile_data[(int)(a + nd2) - 1]);
        if (autoflush) {
          fflush(filestar);
        }
      }
      /* 'getgkmweights:65' a = a+1; */
    }
    a += i5;
  }
  emxFree_char_T(&kmer);
  emxFree_char_T(&s);
  emxFree_real_T(&cfile);
  emxFree_real_T(&comb);
  /* 'getgkmweights:68' fclose(fid); */
  cfclose(fileid);
  /*  dlmwrite([fn '_negmat.out'], [mat; GCpos1 GCneg1 0 0], 'precision',10); */
  /* 'getgkmweights:70' fid = fopen(sprintf('%s_negmat.out', fn), 'w'); */
  i = b_varargin_1->size[0] * b_varargin_1->size[1];
  b_varargin_1->size[0] = 1;
  b_varargin_1->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(b_varargin_1, i);
  b_varargin_1_data = b_varargin_1->data;
  m = varargin_1->size[1];
  for (i = 0; i < m; i++) {
    b_varargin_1_data[i] = varargin_1_data[i];
  }
  b_varargin_1_data[varargin_1->size[1]] = '\x00';
  i = S->size[0] * S->size[1];
  S->size[0] = 1;
  S->size[1] = varargin_1->size[1] + 1;
  emxEnsureCapacity_char_T(S, i);
  S_data = S->data;
  m = varargin_1->size[1];
  for (i = 0; i < m; i++) {
    S_data[i] = varargin_1_data[i];
  }
  emxInit_char_T(&filename, 2);
  S_data[varargin_1->size[1]] = '\x00';
  m = snprintf(NULL, 0, "%s_negmat.out", &S_data[0]);
  i = filename->size[0] * filename->size[1];
  filename->size[0] = 1;
  filename->size[1] = m + 1;
  emxEnsureCapacity_char_T(filename, i);
  S_data = filename->data;
  snprintf(&S_data[0], (size_t)(m + 1), "%s_negmat.out", &b_varargin_1_data[0]);
  i = filename->size[0] * filename->size[1];
  if (1 > m) {
    filename->size[1] = 0;
  } else {
    filename->size[1] = m;
  }
  emxEnsureCapacity_char_T(filename, i);
  fileid = cfopen(filename, "wb");
  /* 'getgkmweights:71' matsize = size(mat); */
  /* 'getgkmweights:72' for i=1:matsize(1) */
  b_NULL = NULL;
  /* 'getgkmweights:73' fprintf(fid, '%.10f,%.10f,%.10f,%.10f\n', mat(i,1),
   * mat(i,2), mat(i,3), mat(i,4)); */
  getfilestar(fileid, &filestar, &autoflush);
  emxFree_char_T(&b_varargin_1);
  emxFree_char_T(&filename);
  emxFree_char_T(&S);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%.10f,%.10f,%.10f,%.10f\n", mat[0], mat[4], mat[8],
            mat[12]);
    if (autoflush) {
      fflush(filestar);
    }
  }
  /* 'getgkmweights:73' fprintf(fid, '%.10f,%.10f,%.10f,%.10f\n', mat(i,1),
   * mat(i,2), mat(i,3), mat(i,4)); */
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%.10f,%.10f,%.10f,%.10f\n", mat[1], mat[5], mat[9],
            mat[13]);
    if (autoflush) {
      fflush(filestar);
    }
  }
  /* 'getgkmweights:73' fprintf(fid, '%.10f,%.10f,%.10f,%.10f\n', mat(i,1),
   * mat(i,2), mat(i,3), mat(i,4)); */
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%.10f,%.10f,%.10f,%.10f\n", mat[2], mat[6], mat[10],
            mat[14]);
    if (autoflush) {
      fflush(filestar);
    }
  }
  /* 'getgkmweights:73' fprintf(fid, '%.10f,%.10f,%.10f,%.10f\n', mat(i,1),
   * mat(i,2), mat(i,3), mat(i,4)); */
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%.10f,%.10f,%.10f,%.10f\n", mat[3], mat[7], mat[11],
            mat[15]);
    if (autoflush) {
      fflush(filestar);
    }
  }
  /* 'getgkmweights:75' fprintf(fid, '%.10f,%.10f,%.10f,%.10f\n', GCpos1,
   * GCneg1, 0.0, 0.0); */
  b_NULL = NULL;
  getfilestar(fileid, &filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    fprintf(filestar, "%.10f,%.10f,%.10f,%.10f\n", validatedHoleFilling_idx_0,
            GCneg1, 0.0, 0.0);
    if (autoflush) {
      fflush(filestar);
    }
  }
}

/* End of code generation (getgkmweights.c) */
