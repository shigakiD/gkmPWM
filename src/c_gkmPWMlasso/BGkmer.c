/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * BGkmer.c
 *
 * Code generation for function 'BGkmer'
 *
 */

/* Include files */
#include "BGkmer.h"
#include "diff.h"
#include "gkmPWMlasso4_emxutil.h"
#include "gkmPWMlasso4_types.h"
#include "mod.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
/*
 * function [negweights] = BGkmer(mat,GC,c,rcnum,l,k,RC)
 */
void BGkmer(const double mat[16], double GC, const emxArray_real_T *c,
            double rcnum, double l, double k, bool RC,
            emxArray_real_T *negweights)
{
  cell_wrap_12 *p_data;
  emxArray_cell_wrap_12 *p;
  emxArray_real_T *b_c;
  emxArray_real_T *c2;
  emxArray_real_T *dc;
  emxArray_real_T *dc2;
  emxArray_real_T *seqvec;
  emxArray_real_T *seqvec2;
  double b_p[16];
  double startvec[4];
  const double *c_data;
  double GCmat_idx_0_tmp;
  double GCmat_idx_1;
  double GCmat_idx_2;
  double c_tmp;
  double tmp;
  double xtmp;
  double *b_c_data;
  double *c2_data;
  double *dc2_data;
  double *dc_data;
  double *negweights_data;
  double *seqvec2_data;
  double *seqvec_data;
  unsigned int a;
  int b_i;
  int b_j1;
  int i;
  int j2;
  int m;
  int nd2;
  int u1;
  c_data = c->data;
  /* 'BGkmer:2' len = length(c); */
  if ((c->size[0] == 0) || (c->size[1] == 0)) {
    u1 = 0;
  } else {
    m = c->size[0];
    u1 = c->size[1];
    if (m >= u1) {
      u1 = m;
    }
  }
  /* 'BGkmer:3' alen = len-rcnum; */
  /* 'BGkmer:4' negweights = zeros(len*4^k,1); */
  c_tmp = pow(4.0, k);
  m = (int)((double)u1 * c_tmp);
  b_j1 = negweights->size[0];
  negweights->size[0] = m;
  emxEnsureCapacity_real_T(negweights, b_j1);
  negweights_data = negweights->data;
  for (b_j1 = 0; b_j1 < m; b_j1++) {
    negweights_data[b_j1] = 0.0;
  }
  emxInit_real_T(&c2, 2);
  emxInit_real_T(&seqvec2, 2);
  emxInit_real_T(&dc2, 2);
  emxInit_cell_wrap_12(&p);
  /* 'BGkmer:5' GCmat = [0.5-GC/2 GC/2 GC/2 0.5-GC/2]; */
  GCmat_idx_0_tmp = 0.5 - GC / 2.0;
  GCmat_idx_1 = GC / 2.0;
  GCmat_idx_2 = GC / 2.0;
  /* 'BGkmer:6' c2 = 0; */
  b_j1 = c2->size[0] * c2->size[1];
  c2->size[0] = 1;
  c2->size[1] = 1;
  emxEnsureCapacity_real_T(c2, b_j1);
  c2_data = c2->data;
  c2_data[0] = 0.0;
  /*  To appease Coder */
  /* 'BGkmer:7' seqvec2 = 0; */
  b_j1 = seqvec2->size[0] * seqvec2->size[1];
  seqvec2->size[0] = 1;
  seqvec2->size[1] = 1;
  emxEnsureCapacity_real_T(seqvec2, b_j1);
  seqvec2_data = seqvec2->data;
  seqvec2_data[0] = 0.0;
  /*  To appease Coder */
  /* 'BGkmer:8' dc2 = 0; */
  b_j1 = dc2->size[0] * dc2->size[1];
  dc2->size[0] = 1;
  dc2->size[1] = 1;
  emxEnsureCapacity_real_T(dc2, b_j1);
  dc2_data = dc2->data;
  dc2_data[0] = 0.0;
  /*  To appease Coder */
  /* 'BGkmer:9' tmp = 0; */
  tmp = 0.0;
  /*  To appease Coder */
  /* 'BGkmer:10' p = cell(l,1); */
  b_j1 = p->size[0];
  p->size[0] = (int)l;
  emxEnsureCapacity_cell_wrap_12(p, b_j1);
  p_data = p->data;
  /* 'BGkmer:11' p = coder.nullcopy(p); */
  /* 'BGkmer:12' p{1} = eye(4); */
  for (b_j1 = 0; b_j1 < 16; b_j1++) {
    p_data[0].f1[b_j1] = 0.0;
  }
  p_data[0].f1[0] = 1.0;
  p_data[0].f1[5] = 1.0;
  p_data[0].f1[10] = 1.0;
  p_data[0].f1[15] = 1.0;
  /* 'BGkmer:13' for i = 1:l-1 */
  b_j1 = (int)(l - 1.0);
  for (i = 0; i < b_j1; i++) {
    /* 'BGkmer:14' p{i+1} = p{i}*mat; */
    for (b_i = 0; b_i < 4; b_i++) {
      for (nd2 = 0; nd2 < 4; nd2++) {
        m = nd2 << 2;
        b_p[b_i + m] =
            ((p_data[i].f1[b_i] * mat[m] + p_data[i].f1[b_i + 4] * mat[m + 1]) +
             p_data[i].f1[b_i + 8] * mat[m + 2]) +
            p_data[i].f1[b_i + 12] * mat[m + 3];
      }
    }
    for (b_i = 0; b_i < 16; b_i++) {
      p_data[(int)(((double)i + 1.0) + 1.0) - 1].f1[b_i] = b_p[b_i];
    }
  }
  emxInit_real_T(&seqvec, 2);
  /* 'BGkmer:16' seqvec = zeros(4^k, k); */
  b_j1 = seqvec->size[0] * seqvec->size[1];
  b_i = (int)c_tmp;
  seqvec->size[0] = (int)c_tmp;
  nd2 = (int)k;
  seqvec->size[1] = (int)k;
  emxEnsureCapacity_real_T(seqvec, b_j1);
  seqvec_data = seqvec->data;
  m = (int)c_tmp * (int)k;
  for (b_j1 = 0; b_j1 < m; b_j1++) {
    seqvec_data[b_j1] = 0.0;
  }
  /* 'BGkmer:17' for i = 1:k */
  for (i = 0; i < nd2; i++) {
    /* 'BGkmer:18' for j = 1:4^k */
    for (m = 0; m < b_i; m++) {
      /* 'BGkmer:19' seqvec(j,i) = mod(floor((j-1)/4^(i-1)), 4)+1; */
      seqvec_data[m + seqvec->size[0] * i] =
          b_mod(floor((((double)m + 1.0) - 1.0) /
                      pow(4.0, ((double)i + 1.0) - 1.0)),
                4.0) +
          1.0;
    }
  }
  /* 'BGkmer:22' if RC */
  if (RC) {
    /* 'BGkmer:23' seqvec2 = 5-fliplr(seqvec); */
    b_j1 = seqvec2->size[0] * seqvec2->size[1];
    seqvec2->size[0] = seqvec->size[0];
    seqvec2->size[1] = seqvec->size[1];
    emxEnsureCapacity_real_T(seqvec2, b_j1);
    seqvec2_data = seqvec2->data;
    m = seqvec->size[0] * seqvec->size[1];
    for (b_j1 = 0; b_j1 < m; b_j1++) {
      seqvec2_data[b_j1] = seqvec_data[b_j1];
    }
    m = seqvec->size[0];
    nd2 = seqvec->size[1] >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (seqvec->size[1] - b_j1) - 1;
      for (i = 0; i < m; i++) {
        xtmp = seqvec2_data[i + seqvec2->size[0] * b_j1];
        seqvec2_data[i + seqvec2->size[0] * b_j1] =
            seqvec2_data[i + seqvec2->size[0] * j2];
        seqvec2_data[i + seqvec2->size[0] * j2] = xtmp;
      }
    }
    m = seqvec2->size[0] * seqvec2->size[1];
    for (b_j1 = 0; b_j1 < m; b_j1++) {
      seqvec2_data[b_j1] = 5.0 - seqvec2_data[b_j1];
    }
    /* 'BGkmer:24' c2 = l+1-fliplr(c); */
    b_j1 = c2->size[0] * c2->size[1];
    c2->size[0] = c->size[0];
    c2->size[1] = c->size[1];
    emxEnsureCapacity_real_T(c2, b_j1);
    c2_data = c2->data;
    m = c->size[0] * c->size[1];
    for (b_j1 = 0; b_j1 < m; b_j1++) {
      c2_data[b_j1] = c_data[b_j1];
    }
    m = c->size[0];
    nd2 = c->size[1] >> 1;
    for (b_j1 = 0; b_j1 < nd2; b_j1++) {
      j2 = (c->size[1] - b_j1) - 1;
      for (i = 0; i < m; i++) {
        xtmp = c2_data[i + c2->size[0] * b_j1];
        c2_data[i + c2->size[0] * b_j1] = c2_data[i + c2->size[0] * j2];
        c2_data[i + c2->size[0] * j2] = xtmp;
      }
    }
    m = c2->size[0] * c2->size[1];
    for (b_j1 = 0; b_j1 < m; b_j1++) {
      c2_data[b_j1] = (l + 1.0) - c2_data[b_j1];
    }
  }
  /* 'BGkmer:26' a = 1; */
  a = 1U;
  /* 'BGkmer:27' for i = 1:len */
  emxInit_real_T(&dc, 2);
  emxInit_real_T(&b_c, 2);
  for (i = 0; i < u1; i++) {
    /* 'BGkmer:28' dc = diff(c(i,:)); */
    m = c->size[1];
    b_j1 = b_c->size[0] * b_c->size[1];
    b_c->size[0] = 1;
    b_c->size[1] = c->size[1];
    emxEnsureCapacity_real_T(b_c, b_j1);
    b_c_data = b_c->data;
    for (b_j1 = 0; b_j1 < m; b_j1++) {
      b_c_data[b_j1] = c_data[i + c->size[0] * b_j1];
    }
    diff(b_c, dc);
    dc_data = dc->data;
    /* 'BGkmer:29' startvec = GCmat*p{c(i,1)}; */
    for (b_j1 = 0; b_j1 < 16; b_j1++) {
      b_p[b_j1] = p_data[(int)c_data[i] - 1].f1[b_j1];
    }
    for (b_j1 = 0; b_j1 < 4; b_j1++) {
      b_i = b_j1 << 2;
      startvec[b_j1] =
          ((GCmat_idx_0_tmp * b_p[b_i] + GCmat_idx_1 * b_p[b_i + 1]) +
           GCmat_idx_2 * b_p[b_i + 2]) +
          GCmat_idx_0_tmp * b_p[b_i + 3];
    }
    /* 'BGkmer:30' if RC */
    if (RC) {
      /* 'BGkmer:31' dc2 = diff(c2(i,:)); */
      m = c2->size[1];
      b_j1 = b_c->size[0] * b_c->size[1];
      b_c->size[0] = 1;
      b_c->size[1] = c2->size[1];
      emxEnsureCapacity_real_T(b_c, b_j1);
      b_c_data = b_c->data;
      for (b_j1 = 0; b_j1 < m; b_j1++) {
        b_c_data[b_j1] = c2_data[i + c2->size[0] * b_j1];
      }
      diff(b_c, dc2);
      dc2_data = dc2->data;
      /* 'BGkmer:32' startvec2 = GCmat*p{c2(i,1)}; */
    }
    /* 'BGkmer:34' for ii = 1:4^k */
    b_j1 = (int)pow(4.0, k);
    for (nd2 = 0; nd2 < b_j1; nd2++) {
      /* 'BGkmer:35' negweights(a) = startvec(seqvec(ii,1)); */
      m = (int)(a + nd2) - 1;
      negweights_data[m] = startvec[(int)seqvec_data[nd2] - 1];
      /* 'BGkmer:36' if RC */
      if (RC) {
        /* 'BGkmer:37' tmp = startvec(seqvec2(ii,1)); */
        tmp = startvec[(int)seqvec2_data[nd2] - 1];
      }
      /* 'BGkmer:39' for iii = 1:k-1 */
      b_i = (int)(k - 1.0);
      for (j2 = 0; j2 < b_i; j2++) {
        /* 'BGkmer:40' matt = p{dc(iii)+1}; */
        /* 'BGkmer:41' negweights(a) = negweights(a)*matt(seqvec(ii,iii),
         * seqvec(ii,iii+1)); */
        negweights_data[m] *=
            p_data[(int)(dc_data[j2] + 1.0) - 1]
                .f1[((int)seqvec_data[nd2 + seqvec->size[0] * j2] +
                     (((int)seqvec_data[nd2 + seqvec->size[0] * (j2 + 1)] - 1)
                      << 2)) -
                    1];
        /* 'BGkmer:42' if RC */
        if (RC) {
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
      if (RC) {
        /* 'BGkmer:48' negweights(a) = (negweights(a)+tmp); */
        negweights_data[m] += tmp;
      }
      /* 'BGkmer:50' a = a+1; */
    }
    a += b_j1;
  }
  emxFree_cell_wrap_12(&p);
  emxFree_real_T(&dc);
  emxFree_real_T(&seqvec);
  emxFree_real_T(&dc2);
  emxFree_real_T(&seqvec2);
  emxFree_real_T(&c2);
  /* 'BGkmer:53' if rcnum > 0 && RC */
  if ((rcnum > 0.0) && RC) {
    /* 'BGkmer:54' negweights(4^k*alen+1:end) =
     * negweights(4^k*alen+1:end)/sqrt(2); */
    c_tmp = c_tmp * ((double)u1 - rcnum) + 1.0;
    if (c_tmp > negweights->size[0]) {
      b_j1 = 1;
      b_i = -1;
      nd2 = 0;
    } else {
      b_j1 = (int)c_tmp;
      b_i = (int)c_tmp - 2;
      nd2 = negweights->size[0];
    }
    m = (nd2 - b_i) - 1;
    nd2 = b_c->size[0] * b_c->size[1];
    b_c->size[0] = 1;
    b_c->size[1] = m;
    emxEnsureCapacity_real_T(b_c, nd2);
    b_c_data = b_c->data;
    for (nd2 = 0; nd2 < m; nd2++) {
      b_c_data[nd2] = negweights_data[(b_j1 + nd2) - 1] / 1.4142135623730951;
    }
    m = b_c->size[1];
    for (b_j1 = 0; b_j1 < m; b_j1++) {
      negweights_data[(b_i + b_j1) + 1] = b_c_data[b_j1];
    }
  }
  emxFree_real_T(&b_c);
}

/* End of code generation (BGkmer.c) */
