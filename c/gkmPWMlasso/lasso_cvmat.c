/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * lasso_cvmat.c
 *
 * Code generation for function 'lasso_cvmat'
 *
 */

/* Include files */
#include "lasso_cvmat.h"
#include "any.h"
#include "blockedSummation.h"
#include "bsxfun.h"
#include "combineVectorElements.h"
#include "div.h"
#include "gkmPWMlasso.h"
#include "gkmPWMlasso_emxutil.h"
#include "gkmPWMlasso_types.h"
#include "isequal.h"
#include "linspace.h"
#include "mean.h"
#include "minOrMax.h"
#include "mtimes.h"
#include "cblas.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/* Function Declarations */
static void cdescentCycle(const emxArray_real_T *XX0,
                          const emxArray_real_T *XY0, emxArray_real_T *b,
                          emxArray_boolean_T *active, double totalweight,
                          const emxArray_real_T *shrinkFactor,
                          double threshold);

static void computeLambdaMax(const emxArray_real_T *X, const emxArray_real_T *Y,
                             double *lambdaMax, double *nullMSE);

static void lassoFit(const emxArray_real_T *X, const emxArray_real_T *Y,
                     const double lambda_data[], int lambda_size[2],
                     double dfmax, double lambdaMax,
                     emxArray_boolean_T *ever_active, double nullMSE,
                     emxArray_real_T *B, double Intercept_data[],
                     int Intercept_size[2], emxArray_real_T *mspe);

static void o_binary_expand_op(const emxArray_real_T *X,
                               const emxArray_real_T *Y,
                               const double lambda_data[], int lambda_size[2],
                               double dfMax, double lambdaMax,
                               const emxArray_real_T *ex,
                               const emxArray_real_T *b_ex, double nullMSE,
                               emxArray_real_T *B, double x_data[],
                               int x_size[2], emxArray_real_T *mse);

static void or
    (emxArray_boolean_T * potentially_active, const emxArray_boolean_T *active);

static void p_binary_expand_op(emxArray_boolean_T *okrows,
                               const unsigned int outsize[2],
                               const unsigned int b_outsize[2]);

static void w_binary_expand_op(emxArray_boolean_T *ever_active,
                               const emxArray_real_T *A,
                               const emxArray_real_T *y);

/* Function Definitions */
/*
 * function [b,active] = cdescentCycle(XX0, wX0, XY0, ...
 *     b, active, totalweight, shrinkFactor, threshold)
 */
static void cdescentCycle(const emxArray_real_T *XX0,
                          const emxArray_real_T *XY0, emxArray_real_T *b,
                          emxArray_boolean_T *active, double totalweight,
                          const emxArray_real_T *shrinkFactor, double threshold)
{
  emxArray_real_T *a;
  const double *XX0_data;
  const double *XY0_data;
  const double *shrinkFactor_data;
  double bj;
  double x;
  double *a_data;
  double *b_data;
  int i;
  int i1;
  int j;
  int loop_ub;
  bool *active_data;
  shrinkFactor_data = shrinkFactor->data;
  active_data = active->data;
  b_data = b->data;
  XY0_data = XY0->data;
  XX0_data = XX0->data;
  /* -lassoFit */
  /*  =================================================== */
  /*                  cdescentCycle()  */
  /*  =================================================== */
  /*  */
  /* 'lasso_cvmat:686' [m,n] = size(XX0); */
  /*  for j=find(active) */
  /* 'lasso_cvmat:689' for j = 1:n */
  i = XX0->size[1];
  emxInit_real_T(&a, 2);
  for (j = 0; j < i; j++) {
    /* 'lasso_cvmat:690' if active(j) == 0 */
    if (active_data[j]) {
      /* 'lasso_cvmat:693' bjold = b(j); */
      /*  Regress j-th partial residuals on j-th predictor */
      /* 'lasso_cvmat:695' bj = (XY0(j)-XX0(j,:)*b+XX0(j,j)*b(j)) / totalweight;
       */
      loop_ub = XX0->size[1];
      i1 = a->size[0] * a->size[1];
      a->size[0] = 1;
      a->size[1] = XX0->size[1];
      emxEnsureCapacity_real_T(a, i1);
      a_data = a->data;
      for (i1 = 0; i1 < loop_ub; i1++) {
        a_data[i1] = XX0_data[j + XX0->size[0] * i1];
      }
      if (XX0->size[1] < 1) {
        bj = 0.0;
      } else {
        bj = cblas_ddot((blasint)XX0->size[1], &a_data[0], (blasint)1,
                        &b_data[0], (blasint)1);
      }
      bj = ((XY0_data[j] - bj) + XX0_data[j + XX0->size[0] * j] * b_data[j]) /
           totalweight;
      /*  Soft thresholding */
      /* 'lasso_cvmat:697' b(j) = sign(bj) .* max((abs(bj) - threshold), 0) ./
       * shrinkFactor(j); */
      x = bj;
      if (bj < 0.0) {
        x = -1.0;
      } else if (bj > 0.0) {
        x = 1.0;
      }
      b_data[j] = x * fmax(fabs(bj) - threshold, 0.0) / shrinkFactor_data[j];
      /* 'lasso_cvmat:698' if b(j) == 0 */
      if (b_data[j] == 0.0) {
        /* 'lasso_cvmat:699' active(j) = 0; */
        active_data[j] = false;
      }
    }
  }
  emxFree_real_T(&a);
}

/*
 * function [lambdaMax, nullMSE] = computeLambdaMax(X, Y, weights, alpha,
 * standardize)
 */
static void computeLambdaMax(const emxArray_real_T *X, const emxArray_real_T *Y,
                             double *lambdaMax, double *nullMSE)
{
  emxArray_real_T *A;
  emxArray_real_T *Y0;
  emxArray_real_T *b_y;
  emxArray_real_T *c_y;
  emxArray_real_T *dotp;
  emxArray_real_T *y;
  const double *Y_data;
  double d;
  double muY;
  double *Y0_data;
  double *dotp_data;
  double *y_data;
  int k;
  int nx;
  Y_data = Y->data;
  /* -thresholdScreen */
  /*  =================================================== */
  /*                  computeLambdaMaX()  */
  /*  =================================================== */
  /*  */
  /*  lambdaMax is the penalty term (lambda) beyond which coefficients */
  /*  are guaranteed to be all zero. */
  /*  */
  /*  nullMse is the mse of the fit using just a constant term. */
  /*  It is provided in this function as a convenience, because it needs  */
  /*  to be calculated in the same context as lambdaMax whenever */
  /*  lambdaMax is calculated. */
  /* 'lasso_cvmat:732' if ~isempty(weights) */
  /* 'lasso_cvmat:737' else */
  /* 'lasso_cvmat:738' observationWeights = false; */
  /* 'lasso_cvmat:741' [N,~] = size(X); */
  /*  If we were asked to standardize the predictors, do so here because */
  /*  the calculation of lambdaMax needs the predictors as we will use */
  /*  them to perform fits. */
  /* 'lasso_cvmat:747' if standardize */
  /* 'lasso_cvmat:764' else */
  /* 'lasso_cvmat:765' if ~observationWeights */
  /*  Center */
  /* 'lasso_cvmat:767' muX = mean(X,1); */
  emxInit_real_T(&y, 2);
  y_data = y->data;
  if (X->size[1] == 0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    colMajorFlatIter(X, X->size[0], y);
    y_data = y->data;
  }
  emxInit_real_T(&Y0, 1);
  /* 'lasso_cvmat:768' X0 = bsxfun(@minus,X,muX); */
  /*  If using observation weights, make a weighted copy of the  */
  /*  predictor matrix, for use in weighted dot products. */
  /* 'lasso_cvmat:779' if observationWeights */
  /* 'lasso_cvmat:783' if ~observationWeights */
  /* 'lasso_cvmat:784' muY = mean(Y); */
  muY = blockedSummation(Y, Y->size[0]) / (double)Y->size[0];
  /*  Y0 = bsxfun(@minus,Y,muY); */
  /* 'lasso_cvmat:789' Y0 = Y - muY; */
  k = Y0->size[0];
  Y0->size[0] = Y->size[0];
  emxEnsureCapacity_real_T(Y0, k);
  Y0_data = Y0->data;
  nx = Y->size[0];
  for (k = 0; k < nx; k++) {
    Y0_data[k] = Y_data[k] - muY;
  }
  emxInit_real_T(&b_y, 1);
  emxInit_real_T(&c_y, 2);
  /*  Calculate max lambda that permits non-zero coefficients */
  /*  */
  /* 'lasso_cvmat:793' if ~observationWeights */
  /* 'lasso_cvmat:794' dotp = abs(X0' * Y0); */
  k = c_y->size[0] * c_y->size[1];
  c_y->size[0] = 1;
  c_y->size[1] = y->size[1];
  emxEnsureCapacity_real_T(c_y, k);
  dotp_data = c_y->data;
  nx = y->size[1];
  for (k = 0; k < nx; k++) {
    dotp_data[k] = y_data[k] / (double)X->size[0];
  }
  emxFree_real_T(&y);
  emxInit_real_T(&A, 2);
  bsxfun(X, c_y, A);
  dotp_data = A->data;
  emxFree_real_T(&c_y);
  if ((A->size[0] == 0) || (A->size[1] == 0) || (Y0->size[0] == 0)) {
    k = b_y->size[0];
    b_y->size[0] = A->size[1];
    emxEnsureCapacity_real_T(b_y, k);
    y_data = b_y->data;
    nx = A->size[1];
    for (k = 0; k < nx; k++) {
      y_data[k] = 0.0;
    }
  } else {
    k = b_y->size[0];
    b_y->size[0] = A->size[1];
    emxEnsureCapacity_real_T(b_y, k);
    y_data = b_y->data;
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)A->size[1],
                (blasint)1, (blasint)A->size[0], 1.0, &dotp_data[0],
                (blasint)A->size[0], &Y0_data[0], (blasint)Y0->size[0], 0.0,
                &y_data[0], (blasint)A->size[1]);
  }
  emxFree_real_T(&A);
  emxInit_real_T(&dotp, 1);
  nx = b_y->size[0];
  k = dotp->size[0];
  dotp->size[0] = b_y->size[0];
  emxEnsureCapacity_real_T(dotp, k);
  dotp_data = dotp->data;
  for (k = 0; k < nx; k++) {
    dotp_data[k] = fabs(y_data[k]);
  }
  /* 'lasso_cvmat:795' lambdaMax = max(dotp) / (N*alpha); */
  nx = dotp->size[0];
  if (dotp->size[0] <= 2) {
    if (dotp->size[0] == 1) {
      muY = dotp_data[0];
    } else if (dotp_data[0] < dotp_data[dotp->size[0] - 1]) {
      muY = dotp_data[dotp->size[0] - 1];
    } else {
      muY = dotp_data[0];
    }
  } else {
    muY = dotp_data[0];
    for (k = 2; k <= nx; k++) {
      d = dotp_data[k - 1];
      if (muY < d) {
        muY = d;
      }
    }
  }
  emxFree_real_T(&dotp);
  *lambdaMax = muY / (double)X->size[0];
  /* 'lasso_cvmat:801' if ~observationWeights */
  /* 'lasso_cvmat:802' nullMSE = mean(Y0.^2); */
  k = b_y->size[0];
  b_y->size[0] = Y0->size[0];
  emxEnsureCapacity_real_T(b_y, k);
  y_data = b_y->data;
  nx = Y0->size[0];
  for (k = 0; k < nx; k++) {
    muY = Y0_data[k];
    y_data[k] = pow(muY, 2.0);
  }
  emxFree_real_T(&Y0);
  *nullMSE = blockedSummation(b_y, b_y->size[0]) / (double)b_y->size[0];
  emxFree_real_T(&b_y);
}

/*
 * function [B,Intercept,lambda,mspe] = ...
 *     lassoFit(X,Y,weights,lambda,alpha,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE,maxIter)
 */
static void lassoFit(const emxArray_real_T *X, const emxArray_real_T *Y,
                     const double lambda_data[], int lambda_size[2],
                     double dfmax, double lambdaMax,
                     emxArray_boolean_T *ever_active, double nullMSE,
                     emxArray_real_T *B, double Intercept_data[],
                     int Intercept_size[2], emxArray_real_T *mspe)
{
  emxArray_boolean_T *active;
  emxArray_boolean_T *potentially_active;
  emxArray_int32_T *r1;
  emxArray_int32_T *r2;
  emxArray_real_T *A;
  emxArray_real_T *Intercept;
  emxArray_real_T *X0;
  emxArray_real_T *XX0;
  emxArray_real_T *XY0;
  emxArray_real_T *Y0;
  emxArray_real_T *b;
  emxArray_real_T *b_B;
  emxArray_real_T *b_b;
  emxArray_real_T *b_r;
  emxArray_real_T *b_y;
  emxArray_real_T *bold;
  emxArray_real_T *c_y;
  emxArray_real_T *fit;
  emxArray_real_T *muX;
  emxArray_real_T *shrinkFactor;
  emxArray_real_T *y;
  const double *X_data;
  const double *Y_data;
  double a;
  double absx;
  double lam;
  double muY;
  double *B_data;
  double *X0_data;
  double *XY0_data;
  double *Y0_data;
  double *b_b_data;
  double *b_data;
  double *bold_data;
  double *muX_data;
  double *y_data;
  int N;
  int acoef;
  int b_i;
  int end;
  int i;
  int i1;
  int k;
  int loop_ub;
  int sizes_idx_1;
  int trueCount;
  int u1;
  int *r3;
  signed char input_sizes_idx_0;
  signed char input_sizes_idx_1;
  bool empty_non_axis_sizes;
  bool exitg1;
  bool exitg2;
  bool guard1 = false;
  bool *active_data;
  bool *ever_active_data;
  bool *potentially_active_data;
  Y_data = Y->data;
  X_data = X->data;
  /*  =================================================== */
  /*                  lassoFit()  */
  /*  =================================================== */
  /*  */
  /*  ------------------------------------------------------ */
  /*  Perform model fit for each lambda and the given alpha */
  /*  ------------------------------------------------------ */
  /* 'lasso_cvmat:486' [N,P] = size(X); */
  N = X->size[0];
  /* 'lasso_cvmat:487' nLambda = length(lambda); */
  if ((lambda_size[0] == 0) || (lambda_size[1] == 0)) {
    u1 = 0;
  } else {
    acoef = lambda_size[0];
    u1 = lambda_size[1];
    if (acoef >= u1) {
      u1 = acoef;
    }
  }
  emxInit_real_T(&y, 2);
  emxInit_real_T(&A, 2);
  /*  If X has any constant columns, we want to exclude them from the */
  /*  coordinate descent calculations.  The corresponding coefficients */
  /*  will be returned as zero. */
  /*  constantPredictors = (range(X)==0); */
  /* 'lasso_cvmat:493' constantPredictors = ((max(X)-min(X))==0); */
  /* 'lasso_cvmat:494' ever_active = ever_active & ~constantPredictors; */
  d_maximum(X, A);
  bold_data = A->data;
  minimum(X, y);
  y_data = y->data;
  if (A->size[1] == 1) {
    acoef = y->size[1];
  } else {
    acoef = A->size[1];
  }
  if ((A->size[1] == y->size[1]) && (ever_active->size[1] == acoef)) {
    loop_ub = ever_active->size[1] - 1;
    i = ever_active->size[0] * ever_active->size[1];
    ever_active->size[0] = 1;
    emxEnsureCapacity_boolean_T(ever_active, i);
    ever_active_data = ever_active->data;
    for (i = 0; i <= loop_ub; i++) {
      ever_active_data[i] =
          (ever_active_data[i] && (bold_data[i] - y_data[i] != 0.0));
    }
  } else {
    w_binary_expand_op(ever_active, A, y);
    ever_active_data = ever_active->data;
  }
  emxInit_real_T(&Y0, 1);
  /*  === standardization and weights === */
  /*  */
  /* 'lasso_cvmat:498' observationWeights = ~isempty(weights); */
  /* 'lasso_cvmat:499' if ~isempty(weights) */
  /* 'lasso_cvmat:506' if ~observationWeights */
  /* 'lasso_cvmat:507' muY = mean(Y); */
  muY = blockedSummation(Y, Y->size[0]) / (double)Y->size[0];
  /* 'lasso_cvmat:511' Y0 = bsxfun(@minus,Y,muY); */
  i = Y0->size[0];
  Y0->size[0] = Y->size[0];
  emxEnsureCapacity_real_T(Y0, i);
  Y0_data = Y0->data;
  if (Y->size[0] != 0) {
    acoef = (Y->size[0] != 1);
    i = Y->size[0] - 1;
    for (k = 0; k <= i; k++) {
      Y0_data[k] = Y_data[acoef * k] - muY;
    }
  }
  /* 'lasso_cvmat:513' if standardize */
  /* 'lasso_cvmat:528' else */
  /* 'lasso_cvmat:529' if ~observationWeights */
  /*  Center */
  /* 'lasso_cvmat:531' muX = mean(X,1); */
  emxInit_real_T(&muX, 2);
  if (X->size[1] == 0) {
    muX->size[1] = 0;
  } else {
    colMajorFlatIter(X, X->size[0], muX);
  }
  i = muX->size[0] * muX->size[1];
  muX->size[0] = 1;
  emxEnsureCapacity_real_T(muX, i);
  muX_data = muX->data;
  loop_ub = muX->size[1] - 1;
  for (i = 0; i <= loop_ub; i++) {
    muX_data[i] /= (double)X->size[0];
  }
  emxInit_real_T(&X0, 2);
  emxInit_real_T(&b, 1);
  /* 'lasso_cvmat:532' X0 = bsxfun(@minus,X,muX); */
  bsxfun(X, muX, X0);
  X0_data = X0->data;
  /* 'lasso_cvmat:533' sigmaX = 1; */
  /*  If using observation weights, make a weighted copy of the  */
  /*  predictor matrix, to save time in the weighted partial regressions. */
  /* 'lasso_cvmat:544' if observationWeights */
  /* 'lasso_cvmat:547' else */
  /* 'lasso_cvmat:548' wX0 = X0; */
  /* 'lasso_cvmat:549' totalweight = N; */
  /*  b will be the current coefficient estimate, iteratively updated. */
  /*  Because we retain b from one value of lambda to the next, */
  /*  we get a de facto warm start. */
  /* 'lasso_cvmat:555' b = zeros(P,1); */
  i = b->size[0];
  b->size[0] = X->size[1];
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  loop_ub = X->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_data[i] = 0.0;
  }
  /*  Preallocate the returned matrix of coefficients, B, and the intercepts. */
  /* 'lasso_cvmat:558' B = zeros(P,nLambda); */
  i = B->size[0] * B->size[1];
  B->size[0] = X->size[1];
  B->size[1] = u1;
  emxEnsureCapacity_real_T(B, i);
  B_data = B->data;
  loop_ub = X->size[1] * u1;
  for (i = 0; i < loop_ub; i++) {
    B_data[i] = 0.0;
  }
  emxInit_boolean_T(&active, 2);
  /* 'lasso_cvmat:560' active = false(1,P); */
  i = active->size[0] * active->size[1];
  active->size[0] = 1;
  active->size[1] = X->size[1];
  emxEnsureCapacity_boolean_T(active, i);
  active_data = active->data;
  loop_ub = X->size[1];
  for (i = 0; i < loop_ub; i++) {
    active_data[i] = false;
  }
  /* 'lasso_cvmat:562' XY0 = X0'*Y0; */
  emxInit_real_T(&XY0, 1);
  if ((X0->size[0] == 0) || (X0->size[1] == 0) || (Y0->size[0] == 0)) {
    i = XY0->size[0];
    XY0->size[0] = X0->size[1];
    emxEnsureCapacity_real_T(XY0, i);
    XY0_data = XY0->data;
    loop_ub = X0->size[1];
    for (i = 0; i < loop_ub; i++) {
      XY0_data[i] = 0.0;
    }
  } else {
    i = XY0->size[0];
    XY0->size[0] = X0->size[1];
    emxEnsureCapacity_real_T(XY0, i);
    XY0_data = XY0->data;
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)X0->size[1],
                (blasint)1, (blasint)X0->size[0], 1.0, &X0_data[0],
                (blasint)X0->size[0], &Y0_data[0], (blasint)Y0->size[0], 0.0,
                &XY0_data[0], (blasint)X0->size[1]);
  }
  emxInit_real_T(&XX0, 2);
  /* 'lasso_cvmat:563' XX0 = X0'*X0; */
  b_mtimes(X0, X0, XX0);
  /* 'lasso_cvmat:564' for i = 1:nLambda */
  b_i = 0;
  emxInit_real_T(&shrinkFactor, 2);
  emxInit_real_T(&bold, 1);
  emxInit_boolean_T(&potentially_active, 2);
  emxInit_real_T(&fit, 1);
  emxInit_real_T(&b_b, 1);
  emxInit_real_T(&b_r, 2);
  emxInit_int32_T(&r1, 2);
  emxInit_real_T(&b_y, 2);
  exitg1 = false;
  while ((!exitg1) && (b_i <= u1 - 1)) {
    /* 'lasso_cvmat:566' lam = lambda(i); */
    lam = lambda_data[b_i];
    /* 'lasso_cvmat:567' if lam >= lambdaMax */
    if (lambda_data[b_i] >= lambdaMax) {
      b_i++;
    } else {
      /* 'lasso_cvmat:570' threshold = lam * alpha; */
      /*  Denominator in coordinate descent update */
      /* 'lasso_cvmat:573' if standardize */
      /* 'lasso_cvmat:579' else */
      /* 'lasso_cvmat:580' if observationWeights */
      /* 'lasso_cvmat:582' else */
      /* 'lasso_cvmat:583' shrinkFactor = (1/N) * ones(1,N)*(X0.^2) + lam*(1 -
       * alpha); */
      a = 1.0 / (double)N;
      i = b_y->size[0] * b_y->size[1];
      b_y->size[0] = X0->size[0];
      b_y->size[1] = X0->size[1];
      emxEnsureCapacity_real_T(b_y, i);
      y_data = b_y->data;
      loop_ub = X0->size[0] * X0->size[1];
      for (i = 0; i < loop_ub; i++) {
        absx = X0_data[i];
        y_data[i] = pow(absx, 2.0);
      }
      i = A->size[0] * A->size[1];
      A->size[0] = 1;
      A->size[1] = N;
      emxEnsureCapacity_real_T(A, i);
      bold_data = A->data;
      for (i = 0; i < N; i++) {
        bold_data[i] = a;
      }
      if ((N == 0) || (b_y->size[0] == 0) || (b_y->size[1] == 0)) {
        i = shrinkFactor->size[0] * shrinkFactor->size[1];
        shrinkFactor->size[0] = 1;
        shrinkFactor->size[1] = b_y->size[1];
        emxEnsureCapacity_real_T(shrinkFactor, i);
        XY0_data = shrinkFactor->data;
        loop_ub = b_y->size[1];
        for (i = 0; i < loop_ub; i++) {
          XY0_data[i] = 0.0;
        }
      } else {
        i = shrinkFactor->size[0] * shrinkFactor->size[1];
        shrinkFactor->size[0] = 1;
        shrinkFactor->size[1] = b_y->size[1];
        emxEnsureCapacity_real_T(shrinkFactor, i);
        XY0_data = shrinkFactor->data;
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)1,
                    (blasint)b_y->size[1], (blasint)N, 1.0, &bold_data[0],
                    (blasint)1, &y_data[0], (blasint)b_y->size[0], 0.0,
                    &XY0_data[0], (blasint)1);
      }
      /*  Iterative coordinate descent until converged */
      /* 'lasso_cvmat:588' for numIter = 1:maxIter */
      sizes_idx_1 = 0;
      exitg2 = false;
      while ((!exitg2) && (sizes_idx_1 < 100000)) {
        /* 'lasso_cvmat:590' bold = b; */
        i = bold->size[0];
        bold->size[0] = b->size[0];
        emxEnsureCapacity_real_T(bold, i);
        bold_data = bold->data;
        loop_ub = b->size[0];
        for (i = 0; i < loop_ub; i++) {
          bold_data[i] = b_data[i];
        }
        /* 'lasso_cvmat:591' [b,active] = cdescentCycle(XX0,wX0,XY0, ... */
        /* 'lasso_cvmat:592' b,active,totalweight,shrinkFactor,threshold); */
        i = b_b->size[0];
        b_b->size[0] = b->size[0];
        emxEnsureCapacity_real_T(b_b, i);
        b_b_data = b_b->data;
        loop_ub = b->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_b_data[i] = b_data[i];
        }
        cdescentCycle(XX0, XY0, b_b, active, N, shrinkFactor, lam);
        active_data = active->data;
        b_b_data = b_b->data;
        i = b->size[0];
        b->size[0] = b_b->size[0];
        emxEnsureCapacity_real_T(b, i);
        b_data = b->data;
        loop_ub = b_b->size[0];
        for (i = 0; i < loop_ub; i++) {
          b_data[i] = b_b_data[i];
        }
        /* 'lasso_cvmat:593' if norm( (b-bold) ./ (1.0 + abs(bold)), Inf ) <
         * reltol */
        acoef = bold->size[0];
        i = fit->size[0];
        fit->size[0] = bold->size[0];
        emxEnsureCapacity_real_T(fit, i);
        XY0_data = fit->data;
        for (k = 0; k < acoef; k++) {
          XY0_data[k] = fabs(bold_data[k]);
        }
        if (b_b->size[0] == 1) {
          acoef = bold->size[0];
        } else {
          acoef = b_b->size[0];
        }
        if ((b_b->size[0] == bold->size[0]) && (acoef == fit->size[0])) {
          i = bold->size[0];
          bold->size[0] = b_b->size[0];
          emxEnsureCapacity_real_T(bold, i);
          bold_data = bold->data;
          loop_ub = b_b->size[0];
          for (i = 0; i < loop_ub; i++) {
            bold_data[i] = (b_b_data[i] - bold_data[i]) / (XY0_data[i] + 1.0);
          }
        } else {
          v_binary_expand_op(bold, b_b, fit);
          bold_data = bold->data;
        }
        if (bold->size[0] == 0) {
          a = 0.0;
        } else {
          a = 0.0;
          i = bold->size[0];
          for (k = 0; k < i; k++) {
            absx = fabs(bold_data[k]);
            if (absx > a) {
              a = absx;
            }
          }
        }
        guard1 = false;
        if (a < 0.0001) {
          /*  Cycling over the active set converged. */
          /*  Do one full pass through the predictors. */
          /*  If there is no predictor added to the active set, we're done. */
          /*  Otherwise, resume the coordinate descent iterations. */
          /* 'lasso_cvmat:598' bold = b; */
          /* 'lasso_cvmat:599' potentially_active =
           * thresholdScreen(X0,wX0,Y0,b,ever_active,threshold); */
          /* -cdescentCycle */
          /*  =================================================== */
          /*                  thresholdScreen()  */
          /*  =================================================== */
          /* 'lasso_cvmat:712' r = Y0 - X0(:,active)*b(active); */
          k = ever_active->size[1] - 1;
          trueCount = 0;
          for (loop_ub = 0; loop_ub <= k; loop_ub++) {
            if (ever_active_data[loop_ub]) {
              trueCount++;
            }
          }
          i = b_r->size[0] * b_r->size[1];
          b_r->size[0] = 1;
          b_r->size[1] = trueCount;
          emxEnsureCapacity_real_T(b_r, i);
          XY0_data = b_r->data;
          acoef = 0;
          end = ever_active->size[1] - 1;
          trueCount = 0;
          for (loop_ub = 0; loop_ub <= k; loop_ub++) {
            if (ever_active_data[loop_ub]) {
              XY0_data[acoef] = b_b_data[loop_ub];
              acoef++;
              trueCount++;
            }
          }
          i = r1->size[0] * r1->size[1];
          r1->size[0] = 1;
          r1->size[1] = trueCount;
          emxEnsureCapacity_int32_T(r1, i);
          r3 = r1->data;
          acoef = 0;
          for (loop_ub = 0; loop_ub <= end; loop_ub++) {
            if (ever_active_data[loop_ub]) {
              r3[acoef] = loop_ub + 1;
              acoef++;
            }
          }
          loop_ub = X0->size[0];
          i = b_y->size[0] * b_y->size[1];
          b_y->size[0] = X0->size[0];
          b_y->size[1] = r1->size[1];
          emxEnsureCapacity_real_T(b_y, i);
          y_data = b_y->data;
          trueCount = r1->size[1];
          for (i = 0; i < trueCount; i++) {
            for (i1 = 0; i1 < loop_ub; i1++) {
              y_data[i1 + b_y->size[0] * i] =
                  X0_data[i1 + X0->size[0] * (r3[i] - 1)];
            }
          }
          if ((r1->size[1] == 0) || (b_r->size[1] == 0)) {
            loop_ub = X0->size[0];
            i = bold->size[0];
            bold->size[0] = X0->size[0];
            emxEnsureCapacity_real_T(bold, i);
            bold_data = bold->data;
            for (i = 0; i < loop_ub; i++) {
              bold_data[i] = 0.0;
            }
          } else {
            i = bold->size[0];
            bold->size[0] = X0->size[0];
            emxEnsureCapacity_real_T(bold, i);
            bold_data = bold->data;
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                        (blasint)X0->size[0], (blasint)1, (blasint)r1->size[1],
                        1.0, &y_data[0], (blasint)X0->size[0], &XY0_data[0],
                        (blasint)1, 0.0, &bold_data[0], (blasint)X0->size[0]);
          }
          if (Y0->size[0] == bold->size[0]) {
            i = bold->size[0];
            bold->size[0] = Y0->size[0];
            emxEnsureCapacity_real_T(bold, i);
            bold_data = bold->data;
            loop_ub = Y0->size[0];
            for (i = 0; i < loop_ub; i++) {
              bold_data[i] = Y0_data[i] - bold_data[i];
            }
          } else {
            minus(bold, Y0);
            bold_data = bold->data;
          }
          /*  We don't need the (b.*wX2)' term that one might expect, because it
           */
          /*  is zero for the inactive predictors. */
          /* 'lasso_cvmat:715' potentially_active = abs(r' *wX0) > threshold; */
          if ((bold->size[0] == 0) || (X0->size[1] == 0)) {
            i = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = X0->size[1];
            emxEnsureCapacity_real_T(y, i);
            y_data = y->data;
            loop_ub = X0->size[1];
            for (i = 0; i < loop_ub; i++) {
              y_data[i] = 0.0;
            }
          } else {
            i = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = X0->size[1];
            emxEnsureCapacity_real_T(y, i);
            y_data = y->data;
            cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, (blasint)1,
                        (blasint)X0->size[1], (blasint)bold->size[0], 1.0,
                        &bold_data[0], (blasint)bold->size[0], &X0_data[0],
                        (blasint)X0->size[0], 0.0, &y_data[0], (blasint)1);
          }
          acoef = y->size[1];
          i = A->size[0] * A->size[1];
          A->size[0] = 1;
          A->size[1] = y->size[1];
          emxEnsureCapacity_real_T(A, i);
          bold_data = A->data;
          for (k = 0; k < acoef; k++) {
            bold_data[k] = fabs(y_data[k]);
          }
          i = potentially_active->size[0] * potentially_active->size[1];
          potentially_active->size[0] = 1;
          potentially_active->size[1] = A->size[1];
          emxEnsureCapacity_boolean_T(potentially_active, i);
          potentially_active_data = potentially_active->data;
          loop_ub = A->size[1];
          for (i = 0; i < loop_ub; i++) {
            potentially_active_data[i] = (bold_data[i] > lam);
          }
          /* 'lasso_cvmat:600' if any(potentially_active) */
          if (any(potentially_active)) {
            /* 'lasso_cvmat:601' new_active = active | potentially_active; */
            /* 'lasso_cvmat:602' [b,new_active] = cdescentCycle(XX0,wX0,XY0, ...
             */
            /* 'lasso_cvmat:603'
             * b,new_active,totalweight,shrinkFactor,threshold); */
            if (active->size[1] == potentially_active->size[1]) {
              loop_ub = active->size[1] - 1;
              i = potentially_active->size[0] * potentially_active->size[1];
              potentially_active->size[0] = 1;
              potentially_active->size[1] = active->size[1];
              emxEnsureCapacity_boolean_T(potentially_active, i);
              potentially_active_data = potentially_active->data;
              for (i = 0; i <= loop_ub; i++) {
                potentially_active_data[i] =
                    (active_data[i] || potentially_active_data[i]);
              }
            } else {
              or (potentially_active, active);
            }
            cdescentCycle(XX0, XY0, b, potentially_active, N, shrinkFactor,
                          lam);
            potentially_active_data = potentially_active->data;
            b_data = b->data;
          } else {
            /* 'lasso_cvmat:604' else */
            /* 'lasso_cvmat:605' new_active = active; */
            i = potentially_active->size[0] * potentially_active->size[1];
            potentially_active->size[0] = 1;
            potentially_active->size[1] = active->size[1];
            emxEnsureCapacity_boolean_T(potentially_active, i);
            potentially_active_data = potentially_active->data;
            loop_ub = active->size[1];
            for (i = 0; i < loop_ub; i++) {
              potentially_active_data[i] = active_data[i];
            }
          }
          /* 'lasso_cvmat:608' if isequal(new_active, active) */
          if (isequal(potentially_active, active)) {
            exitg2 = true;
          } else {
            /* 'lasso_cvmat:610' else */
            /* 'lasso_cvmat:611' active = new_active; */
            i = active->size[0] * active->size[1];
            active->size[0] = 1;
            active->size[1] = potentially_active->size[1];
            emxEnsureCapacity_boolean_T(active, i);
            active_data = active->data;
            loop_ub = potentially_active->size[1];
            for (i = 0; i < loop_ub; i++) {
              active_data[i] = potentially_active_data[i];
            }
            /* 'lasso_cvmat:614' if norm( (b-bold) ./ (1.0 + abs(bold)), Inf ) <
             * reltol */
            acoef = b_b->size[0];
            i = fit->size[0];
            fit->size[0] = b_b->size[0];
            emxEnsureCapacity_real_T(fit, i);
            XY0_data = fit->data;
            for (k = 0; k < acoef; k++) {
              XY0_data[k] = fabs(b_b_data[k]);
            }
            if (b->size[0] == 1) {
              acoef = b_b->size[0];
            } else {
              acoef = b->size[0];
            }
            if ((b->size[0] == b_b->size[0]) && (acoef == fit->size[0])) {
              i = b_b->size[0];
              b_b->size[0] = b->size[0];
              emxEnsureCapacity_real_T(b_b, i);
              b_b_data = b_b->data;
              loop_ub = b->size[0];
              for (i = 0; i < loop_ub; i++) {
                b_b_data[i] = (b_data[i] - b_b_data[i]) / (XY0_data[i] + 1.0);
              }
            } else {
              v_binary_expand_op(b_b, b, fit);
              b_b_data = b_b->data;
            }
            if (b_b->size[0] == 0) {
              a = 0.0;
            } else {
              a = 0.0;
              i = b_b->size[0];
              for (k = 0; k < i; k++) {
                absx = fabs(b_b_data[k]);
                if (absx > a) {
                  a = absx;
                }
              }
            }
            if (a < 0.0001) {
              exitg2 = true;
            } else {
              guard1 = true;
            }
          }
        } else {
          guard1 = true;
        }
        if (guard1) {
          /* 'lasso_cvmat:619' if numIter == maxIter */
          if (sizes_idx_1 + 1 == 100000) {
            /* 'lasso_cvmat:620' fprintf('stats:lasso:MaxIterReached %f\n',
             * lam); */
            printf("stats:lasso:MaxIterReached %f\n", lam);
            fflush(stdout);
          }
          sizes_idx_1++;
        }
      }
      /* 'lasso_cvmat:624' B(:,i) = b; */
      loop_ub = b->size[0];
      for (i = 0; i < loop_ub; i++) {
        B_data[i + B->size[0] * b_i] = b_data[i];
      }
      /*  Halt if maximum model size ('DFmax') has been met or exceeded. */
      /* 'lasso_cvmat:627' if sum(active) > dfmax */
      acoef = active->size[1];
      if (active->size[1] == 0) {
        trueCount = 0;
      } else {
        trueCount = active_data[0];
        for (k = 2; k <= acoef; k++) {
          trueCount += active_data[k - 1];
        }
      }
      if (trueCount > dfmax) {
        /*  truncate B and lambda output arguments */
        /* 'lasso_cvmat:629' lambda = lambda(1:(i-1)); */
        lambda_size[0] = 1;
        if (1 > b_i) {
          lambda_size[1] = 0;
        } else {
          lambda_size[1] = b_i;
        }
        /* 'lasso_cvmat:630' B = B(:,1:(i-1)); */
        if (1 > b_i) {
          loop_ub = 0;
        } else {
          loop_ub = b_i;
        }
        acoef = B->size[0] - 1;
        trueCount = B->size[0];
        for (i = 0; i < loop_ub; i++) {
          for (i1 = 0; i1 < trueCount; i1++) {
            B_data[i1 + (acoef + 1) * i] = B_data[i1 + B->size[0] * i];
          }
        }
        i = B->size[0] * B->size[1];
        B->size[0] = acoef + 1;
        B->size[1] = loop_ub;
        emxEnsureCapacity_real_T(B, i);
        B_data = B->data;
        exitg1 = true;
      } else {
        /*  Halt if we have exceeded a threshold on the percent of */
        /*  residual variance left unexplained. */
        /* 'lasso_cvmat:636' if ~userSuppliedLambda */
        /*  Calculate mse of the current fit */
        /* 'lasso_cvmat:638' bsig = b ./ sigmaX'; */
        /* 'lasso_cvmat:639' fit = [ones(size(X,1),1) X] * [(muY-muX*bsig);
         * bsig]; */
        if (X->size[0] != 0) {
          end = X->size[0];
        } else {
          end = 0;
        }
        empty_non_axis_sizes = (end == 0);
        if (empty_non_axis_sizes || (X->size[0] != 0)) {
          input_sizes_idx_1 = 1;
        } else {
          input_sizes_idx_1 = 0;
        }
        if (empty_non_axis_sizes || ((X->size[0] != 0) && (X->size[1] != 0))) {
          sizes_idx_1 = X->size[1];
        } else {
          sizes_idx_1 = 0;
        }
        if (muX->size[1] < 1) {
          a = 0.0;
        } else {
          a = cblas_ddot((blasint)muX->size[1], &muX_data[0], (blasint)1,
                         &b_data[0], (blasint)1);
        }
        i = b_y->size[0] * b_y->size[1];
        b_y->size[0] = end;
        b_y->size[1] = input_sizes_idx_1 + sizes_idx_1;
        emxEnsureCapacity_real_T(b_y, i);
        y_data = b_y->data;
        loop_ub = input_sizes_idx_1;
        if (0 <= loop_ub - 1) {
          for (i = 0; i < end; i++) {
            y_data[i] = 1.0;
          }
        }
        for (i = 0; i < sizes_idx_1; i++) {
          for (i1 = 0; i1 < end; i1++) {
            y_data[i1 + b_y->size[0] * (i + input_sizes_idx_1)] =
                X_data[i1 + end * i];
          }
        }
        i = bold->size[0];
        bold->size[0] = b->size[0] + 1;
        emxEnsureCapacity_real_T(bold, i);
        bold_data = bold->data;
        bold_data[0] = muY - a;
        loop_ub = b->size[0];
        for (i = 0; i < loop_ub; i++) {
          bold_data[i + 1] = b_data[i];
        }
        if ((b_y->size[0] == 0) || (b_y->size[1] == 0)) {
          i = fit->size[0];
          fit->size[0] = b_y->size[0];
          emxEnsureCapacity_real_T(fit, i);
          XY0_data = fit->data;
          loop_ub = b_y->size[0];
          for (i = 0; i < loop_ub; i++) {
            XY0_data[i] = 0.0;
          }
        } else {
          i = fit->size[0];
          fit->size[0] = b_y->size[0];
          emxEnsureCapacity_real_T(fit, i);
          XY0_data = fit->data;
          cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                      (blasint)b_y->size[0], (blasint)1, (blasint)b_y->size[1],
                      1.0, &y_data[0], (blasint)b_y->size[0], &bold_data[0],
                      (blasint)bold->size[0], 0.0, &XY0_data[0],
                      (blasint)b_y->size[0]);
        }
        /* 'lasso_cvmat:640' residuals = bsxfun(@minus, Y, fit); */
        /* 'lasso_cvmat:641' if ~observationWeights */
        /* 'lasso_cvmat:642' mspe = mean(residuals.^2); */
        /* 'lasso_cvmat:647' if mspe < 1.0e-3 * nullMSE */
        b_bsxfun(Y, fit, bold);
        bold_data = bold->data;
        i = fit->size[0];
        fit->size[0] = bold->size[0];
        emxEnsureCapacity_real_T(fit, i);
        XY0_data = fit->data;
        loop_ub = bold->size[0];
        for (i = 0; i < loop_ub; i++) {
          absx = bold_data[i];
          XY0_data[i] = pow(absx, 2.0);
        }
        if (blockedSummation(fit, fit->size[0]) / (double)fit->size[0] <
            0.001 * nullMSE) {
          /* 'lasso_cvmat:648' lambda = lambda(1:i); */
          lambda_size[0] = 1;
          lambda_size[1] = b_i + 1;
          /* 'lasso_cvmat:649' B = B(:,1:i); */
          acoef = B->size[0] - 1;
          loop_ub = B->size[0];
          for (i = 0; i <= b_i; i++) {
            for (i1 = 0; i1 < loop_ub; i1++) {
              B_data[i1 + (acoef + 1) * i] = B_data[i1 + B->size[0] * i];
            }
          }
          i = B->size[0] * B->size[1];
          B->size[0] = acoef + 1;
          B->size[1] = b_i + 1;
          emxEnsureCapacity_real_T(B, i);
          B_data = B->data;
          exitg1 = true;
        } else {
          b_i++;
        }
      }
    }
  }
  emxFree_real_T(&A);
  emxFree_int32_T(&r1);
  emxFree_real_T(&b_r);
  emxFree_real_T(&b_b);
  emxFree_real_T(&fit);
  emxFree_boolean_T(&potentially_active);
  emxFree_real_T(&bold);
  emxFree_real_T(&shrinkFactor);
  emxFree_real_T(&XX0);
  emxFree_real_T(&XY0);
  emxFree_boolean_T(&active);
  emxFree_real_T(&b);
  emxFree_real_T(&Y0);
  emxInit_real_T(&b_B, 2);
  /*  of lambda sequence */
  /*  ------------------------------------------ */
  /*  Unwind the centering and scaling (if any) */
  /*  ------------------------------------------ */
  /* 'lasso_cvmat:660' B = bsxfun(@rdivide, B, sigmaX'); */
  i = b_B->size[0] * b_B->size[1];
  b_B->size[0] = B->size[0];
  b_B->size[1] = B->size[1];
  emxEnsureCapacity_real_T(b_B, i);
  XY0_data = b_B->data;
  loop_ub = B->size[0] * B->size[1] - 1;
  for (i = 0; i <= loop_ub; i++) {
    XY0_data[i] = B_data[i];
  }
  c_bsxfun(b_B, B);
  B_data = B->data;
  /* 'lasso_cvmat:661' B(~ever_active,:) = 0; */
  k = ever_active->size[1] - 1;
  trueCount = 0;
  for (b_i = 0; b_i <= k; b_i++) {
    if (!ever_active_data[b_i]) {
      trueCount++;
    }
  }
  emxInit_int32_T(&r2, 2);
  i = r2->size[0] * r2->size[1];
  r2->size[0] = 1;
  r2->size[1] = trueCount;
  emxEnsureCapacity_int32_T(r2, i);
  r3 = r2->data;
  acoef = 0;
  for (b_i = 0; b_i <= k; b_i++) {
    if (!ever_active_data[b_i]) {
      r3[acoef] = b_i + 1;
      acoef++;
    }
  }
  loop_ub = B->size[1];
  for (i = 0; i < loop_ub; i++) {
    trueCount = r2->size[1];
    for (i1 = 0; i1 < trueCount; i1++) {
      B_data[(r3[i1] + B->size[0] * i) - 1] = 0.0;
    }
  }
  emxFree_int32_T(&r2);
  /* 'lasso_cvmat:662' Intercept = muY-muX*B; */
  if ((muX->size[1] == 0) || (B->size[0] == 0) || (B->size[1] == 0)) {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = B->size[1];
    emxEnsureCapacity_real_T(y, i);
    y_data = y->data;
    loop_ub = B->size[1];
    for (i = 0; i < loop_ub; i++) {
      y_data[i] = 0.0;
    }
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = B->size[1];
    emxEnsureCapacity_real_T(y, i);
    y_data = y->data;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, (blasint)1,
                (blasint)B->size[1], (blasint)muX->size[1], 1.0, &muX_data[0],
                (blasint)1, &B_data[0], (blasint)B->size[0], 0.0, &y_data[0],
                (blasint)1);
  }
  emxFree_real_T(&muX);
  Intercept_size[0] = 1;
  Intercept_size[1] = y->size[1];
  loop_ub = y->size[1];
  for (i = 0; i < loop_ub; i++) {
    Intercept_data[i] = muY - y_data[i];
  }
  emxFree_real_T(&y);
  /*  ------------------------------------------ */
  /*  Calculate Mean Prediction Squared Error */
  /*  ------------------------------------------ */
  /* 'lasso_cvmat:667' BwithI = [Intercept; B]; */
  if (Intercept_size[1] != 0) {
    end = Intercept_size[1];
  } else if ((B->size[0] != 0) && (B->size[1] != 0)) {
    end = B->size[1];
  } else {
    end = 0;
    if (B->size[1] > 0) {
      end = B->size[1];
    }
  }
  empty_non_axis_sizes = (end == 0);
  if (empty_non_axis_sizes || (Intercept_size[1] != 0)) {
    input_sizes_idx_0 = 1;
  } else {
    input_sizes_idx_0 = 0;
  }
  if (empty_non_axis_sizes || ((B->size[0] != 0) && (B->size[1] != 0))) {
    acoef = B->size[0];
  } else {
    acoef = 0;
  }
  /* 'lasso_cvmat:668' fits = [ones(size(X,1),1) X]*BwithI; */
  if (X->size[0] != 0) {
    trueCount = X->size[0];
  } else {
    trueCount = 0;
  }
  empty_non_axis_sizes = (trueCount == 0);
  if (empty_non_axis_sizes || (X->size[0] != 0)) {
    input_sizes_idx_1 = 1;
  } else {
    input_sizes_idx_1 = 0;
  }
  if (empty_non_axis_sizes || ((X->size[0] != 0) && (X->size[1] != 0))) {
    sizes_idx_1 = X->size[1];
  } else {
    sizes_idx_1 = 0;
  }
  k = input_sizes_idx_0;
  i = X0->size[0] * X0->size[1];
  X0->size[0] = trueCount;
  X0->size[1] = input_sizes_idx_1 + sizes_idx_1;
  emxEnsureCapacity_real_T(X0, i);
  X0_data = X0->data;
  loop_ub = input_sizes_idx_1;
  if (0 <= loop_ub - 1) {
    for (i = 0; i < trueCount; i++) {
      X0_data[i] = 1.0;
    }
  }
  for (i = 0; i < sizes_idx_1; i++) {
    for (i1 = 0; i1 < trueCount; i1++) {
      X0_data[i1 + X0->size[0] * (i + input_sizes_idx_1)] =
          X_data[i1 + trueCount * i];
    }
  }
  emxInit_real_T(&Intercept, 2);
  i = Intercept->size[0] * Intercept->size[1];
  Intercept->size[0] = input_sizes_idx_0 + acoef;
  Intercept->size[1] = end;
  emxEnsureCapacity_real_T(Intercept, i);
  XY0_data = Intercept->data;
  for (i = 0; i < end; i++) {
    for (i1 = 0; i1 < k; i1++) {
      XY0_data[Intercept->size[0] * i] = Intercept_data[input_sizes_idx_0 * i];
    }
  }
  for (i = 0; i < end; i++) {
    for (i1 = 0; i1 < acoef; i1++) {
      XY0_data[(i1 + input_sizes_idx_0) + Intercept->size[0] * i] =
          B_data[i1 + acoef * i];
    }
  }
  emxInit_real_T(&c_y, 2);
  c_mtimes(X0, Intercept, c_y);
  /* 'lasso_cvmat:669' residuals = bsxfun(@minus, Y, fits); */
  /* 'lasso_cvmat:670' if ~observationWeights */
  /* 'lasso_cvmat:671' mspe = mean(residuals.^2); */
  d_bsxfun(Y, c_y, b_B);
  XY0_data = b_B->data;
  i = b_y->size[0] * b_y->size[1];
  b_y->size[0] = b_B->size[0];
  b_y->size[1] = b_B->size[1];
  emxEnsureCapacity_real_T(b_y, i);
  y_data = b_y->data;
  loop_ub = b_B->size[0] * b_B->size[1];
  emxFree_real_T(&Intercept);
  emxFree_real_T(&c_y);
  emxFree_real_T(&X0);
  for (i = 0; i < loop_ub; i++) {
    absx = XY0_data[i];
    y_data[i] = pow(absx, 2.0);
  }
  emxFree_real_T(&b_B);
  mean(b_y, mspe);
  emxFree_real_T(&b_y);
}

static void o_binary_expand_op(const emxArray_real_T *X,
                               const emxArray_real_T *Y,
                               const double lambda_data[], int lambda_size[2],
                               double dfMax, double lambdaMax,
                               const emxArray_real_T *ex,
                               const emxArray_real_T *b_ex, double nullMSE,
                               emxArray_real_T *B, double x_data[],
                               int x_size[2], emxArray_real_T *mse)
{
  emxArray_boolean_T *c_ex;
  const double *b_ex_data;
  const double *ex_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  bool *c_ex_data;
  ex_data = b_ex->data;
  b_ex_data = ex->data;
  emxInit_boolean_T(&c_ex, 2);
  i = c_ex->size[0] * c_ex->size[1];
  c_ex->size[0] = 1;
  if (b_ex->size[1] == 1) {
    c_ex->size[1] = ex->size[1];
  } else {
    c_ex->size[1] = b_ex->size[1];
  }
  emxEnsureCapacity_boolean_T(c_ex, i);
  c_ex_data = c_ex->data;
  stride_0_1 = (ex->size[1] != 1);
  stride_1_1 = (b_ex->size[1] != 1);
  if (b_ex->size[1] == 1) {
    loop_ub = ex->size[1];
  } else {
    loop_ub = b_ex->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    c_ex_data[i] = (b_ex_data[i * stride_0_1] - ex_data[i * stride_1_1] != 0.0);
  }
  lassoFit(X, Y, lambda_data, lambda_size, dfMax, lambdaMax, c_ex, nullMSE, B,
           x_data, x_size, mse);
  emxFree_boolean_T(&c_ex);
}

static void or
    (emxArray_boolean_T * potentially_active, const emxArray_boolean_T *active)
{
  emxArray_boolean_T *b_active;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  const bool *active_data;
  bool *b_active_data;
  bool *potentially_active_data;
  active_data = active->data;
  potentially_active_data = potentially_active->data;
  emxInit_boolean_T(&b_active, 2);
  i = b_active->size[0] * b_active->size[1];
  b_active->size[0] = 1;
  if (potentially_active->size[1] == 1) {
    b_active->size[1] = active->size[1];
  } else {
    b_active->size[1] = potentially_active->size[1];
  }
  emxEnsureCapacity_boolean_T(b_active, i);
  b_active_data = b_active->data;
  stride_0_1 = (active->size[1] != 1);
  stride_1_1 = (potentially_active->size[1] != 1);
  if (potentially_active->size[1] == 1) {
    loop_ub = active->size[1];
  } else {
    loop_ub = potentially_active->size[1];
  }
  for (i = 0; i < loop_ub; i++) {
    b_active_data[i] = (active_data[i * stride_0_1] ||
                        potentially_active_data[i * stride_1_1]);
  }
  i = potentially_active->size[0] * potentially_active->size[1];
  potentially_active->size[0] = 1;
  potentially_active->size[1] = b_active->size[1];
  emxEnsureCapacity_boolean_T(potentially_active, i);
  potentially_active_data = potentially_active->data;
  loop_ub = b_active->size[1];
  for (i = 0; i < loop_ub; i++) {
    potentially_active_data[i] = b_active_data[i];
  }
  emxFree_boolean_T(&b_active);
}

static void p_binary_expand_op(emxArray_boolean_T *okrows,
                               const unsigned int outsize[2],
                               const unsigned int b_outsize[2])
{
  int b_outsize_idx_0;
  int i;
  int outsize_idx_0;
  bool *okrows_data;
  outsize_idx_0 = (int)outsize[0];
  b_outsize_idx_0 = (int)b_outsize[0];
  i = okrows->size[0];
  if (b_outsize_idx_0 == 1) {
    okrows->size[0] = outsize_idx_0;
  } else {
    okrows->size[0] = b_outsize_idx_0;
  }
  emxEnsureCapacity_boolean_T(okrows, i);
  okrows_data = okrows->data;
  if (b_outsize_idx_0 != 1) {
    outsize_idx_0 = b_outsize_idx_0;
  }
  for (i = 0; i < outsize_idx_0; i++) {
    okrows_data[i] = true;
  }
}

static void w_binary_expand_op(emxArray_boolean_T *ever_active,
                               const emxArray_real_T *A,
                               const emxArray_real_T *y)
{
  emxArray_boolean_T *b_ever_active;
  const double *A_data;
  const double *y_data;
  int b_y;
  int i;
  int stride_0_1;
  int stride_1_1;
  int stride_2_1;
  bool *b_ever_active_data;
  bool *ever_active_data;
  y_data = y->data;
  A_data = A->data;
  ever_active_data = ever_active->data;
  emxInit_boolean_T(&b_ever_active, 2);
  i = b_ever_active->size[0] * b_ever_active->size[1];
  b_ever_active->size[0] = 1;
  if (y->size[1] == 1) {
    b_y = A->size[1];
  } else {
    b_y = y->size[1];
  }
  if (b_y == 1) {
    b_ever_active->size[1] = ever_active->size[1];
  } else if (y->size[1] == 1) {
    b_ever_active->size[1] = A->size[1];
  } else {
    b_ever_active->size[1] = y->size[1];
  }
  emxEnsureCapacity_boolean_T(b_ever_active, i);
  b_ever_active_data = b_ever_active->data;
  stride_0_1 = (ever_active->size[1] != 1);
  stride_1_1 = (A->size[1] != 1);
  stride_2_1 = (y->size[1] != 1);
  if (y->size[1] == 1) {
    b_y = A->size[1];
  } else {
    b_y = y->size[1];
  }
  if (b_y == 1) {
    b_y = ever_active->size[1];
  } else if (y->size[1] == 1) {
    b_y = A->size[1];
  } else {
    b_y = y->size[1];
  }
  for (i = 0; i < b_y; i++) {
    b_ever_active_data[i] =
        (ever_active_data[i * stride_0_1] &&
         (A_data[i * stride_1_1] - y_data[i * stride_2_1] != 0.0));
  }
  i = ever_active->size[0] * ever_active->size[1];
  ever_active->size[0] = 1;
  ever_active->size[1] = b_ever_active->size[1];
  emxEnsureCapacity_boolean_T(ever_active, i);
  ever_active_data = ever_active->data;
  b_y = b_ever_active->size[1];
  for (i = 0; i < b_y; i++) {
    ever_active_data[i] = b_ever_active_data[i];
  }
  emxFree_boolean_T(&b_ever_active);
}

/*
 * function [B,stats] = lasso_mat(X, Y, dfMax, Std, numLambda)
 */
void b_lasso_cvmat(emxArray_real_T *X, emxArray_real_T *Y, double dfMax,
                   emxArray_real_T *B, double stats_Intercept_data[],
                   int stats_Intercept_size[2], double stats_Lambda_data[],
                   int stats_Lambda_size[2], double *stats_Alpha,
                   double stats_DF_data[], int stats_DF_size[2],
                   double stats_MSE_data[], int stats_MSE_size[2])
{
  emxArray_boolean_T *c_B;
  emxArray_boolean_T *c_ex;
  emxArray_boolean_T *okrows;
  emxArray_int32_T *b_r;
  emxArray_real_T *b_B;
  emxArray_real_T *b_X;
  emxArray_real_T *b_ex;
  emxArray_real_T *ex;
  emxArray_real_T *mse;
  double lambda_data[100];
  double x_data[100];
  double lambdaMax;
  double lambdaMin;
  double nullMSE;
  double *X_data;
  double *b_X_data;
  double *ex_data;
  int b_df_data[100];
  int df_data[100];
  int reverseIndices_data[100];
  unsigned int b_outsize[2];
  int lambda_size[2];
  unsigned int outsize[2];
  int x_size[2];
  int i;
  int loop_ub;
  int n;
  int nx;
  int *r1;
  bool *okrows_data;
  X_data = X->data;
  /* LASSO Perform lasso or elastic net regularization for linear regression. */
  /*    [B,STATS] = lasso(X,Y,...) Performs L1-constrained linear least   */
  /*    squares fits (lasso) or L1- and L2-constrained fits (elastic net) */
  /*    relating the predictors in X to the responses in Y. The default is a */
  /*    lasso fit, or constraint on the L1-norm of the coefficients B. */
  /*  */
  /*    Positional parameters: */
  /*  */
  /*      X                A numeric matrix (dimension, say, NxP) */
  /*      Y                A numeric vector of length N */
  /*     */
  /*    Optional input parameters:   */
  /*  */
  /*      'Weights'        Observation weights.  Must be a vector of
   * non-negative */
  /*                       values, of the same length as columns of X.  At least
   */
  /*                       two values must be positive. (default ones(N,1) or */
  /*                       equivalently (1/N)*ones(N,1)). */
  /*      'Alpha'          Elastic net mixing value, or the relative balance */
  /*                       between L2 and L1 penalty (default 1, range (0,1]).
   */
  /*                       Alpha=1 ==> lasso, otherwise elastic net. */
  /*                       Alpha near zero ==> nearly ridge regression. */
  /*      'NumLambda'      The number of lambda values to use, if the parameter
   */
  /*                       'Lambda' is not supplied (default 100).  Ignored */
  /*                       if 'Lambda' is supplied.  LASSO may return fewer */
  /*                       fits than specified by 'NumLambda' if the residual */
  /*                       error of the fits drops below a threshold percentage
   */
  /*                       of the variance of Y. */
  /*      'LambdaRatio'    Ratio between the minimum value and maximum value of
   */
  /*                       lambda to generate, if the  parameter "Lambda" is not
   */
  /*                       supplied.  Legal range is [0,1). Default is 0.0001.
   */
  /*                       If 'LambdaRatio' is zero, LASSO will generate its */
  /*                       default sequence of lambda values but replace the */
  /*                       smallest value in this sequence with the value zero.
   */
  /*                       'LambdaRatio' is ignored if 'Lambda' is supplied. */
  /*      'Lambda'         Lambda values. Will be returned in return argument */
  /*                       STATS in ascending order. The default is to have
   * LASSO */
  /*                       generate a sequence of lambda values, based on
   * 'NumLambda' */
  /*                       and 'LambdaRatio'. LASSO will generate a sequence,
   * based */
  /*                       on the values in X and Y, such that the largest
   * LAMBDA                  */
  /*                       value is just sufficient to produce all zero
   * coefficients B. */
  /*                       You may supply a vector of real, non-negative values
   * of  */
  /*                       lambda for LASSO to use, in place of its default
   * sequence. */
  /*                       If you supply a value for 'Lambda', 'NumLambda' and
   */
  /*                       'LambdaRatio' are ignored. */
  /*      'DFmax'          Maximum number of non-zero coefficients in the model.
   */
  /*                       Can be useful with large numbers of predictors. */
  /*                       Results only for lambda values that satisfy this */
  /*                       degree of sparseness will be returned. Default is */
  /*                       to not limit the number of non-zero coefficients. */
  /*      'Standardize'    Whether to scale X prior to fitting the model */
  /*                       sequence. This affects whether the regularization is
   */
  /*                       applied to the coefficients on the standardized */
  /*                       scale or the original scale. The results are always
   */
  /*                       presented on the original data scale. Default is */
  /*                       TRUE, do scale X. */
  /*                       Note: X and Y are always centered. */
  /*      'RelTol'         Convergence threshold for coordinate descent
   * algorithm. */
  /*                       The coordinate descent iterations will terminate */
  /*                       when the relative change in the size of the */
  /*                       estimated coefficients B drops below this threshold.
   */
  /*                       Default: 1e-4. Legal range is (0,1). */
  /*      'MaxIter'        Maximum number of iterations allowed.  Default is
   * 1e5. */
  /*      'PredictorNames' A string/cell array of names for the predictor
   * variables, */
  /*                       in the order in which they appear in X.  */
  /*                       Default: {} */
  /*      'Options'        A structure that contains options specifying whether
   * to */
  /*                       conduct cross-validation evaluations in parallel, and
   */
  /*                       options specifying how to use random numbers when
   * computing */
  /*                       cross validation partitions. This argument can be
   * created */
  /*                       by a call to STATSET. CROSSVAL uses the following
   * fields: */
  /*                         'UseParallel' */
  /*                         'UseSubstreams' */
  /*                         'Streams' */
  /*                       For information on these fields see PARALLELSTATS. */
  /*                       NOTE: If supplied, 'Streams' must be of length one.
   */
  /*     */
  /*    Return values: */
  /*      B                The fitted coefficients for each model.  */
  /*                       B will have dimension PxL, where  */
  /*                       P = size(X,2) is the number of predictors, and */
  /*                       L = length(lambda). */
  /*      STATS            STATS is a struct that contains information about the
   */
  /*                       sequence of model fits corresponding to the columns
   */
  /*                       of B. STATS contains the following fields: */
  /*  */
  /*        'Intercept'    The intercept term for each model. Dimension 1xL. */
  /*        'Lambda'       The sequence of lambda penalties used, in ascending
   * order.  */
  /*                       Dimension 1xL. */
  /*        'Alpha'        The elastic net mixing value that was used. */
  /*        'DF'           The number of nonzero coefficients in B for each */
  /*                       value of lambda. Dimension 1xL. */
  /*        'MSE'          The mean squared error of the fitted model for each
   */
  /*                       value of lambda. If cross-validation was performed,
   */
  /*                       the values for 'MSE' represent Mean Prediction */
  /*                       Squared Error for each value of lambda, as calculated
   */
  /*                       by cross-validation. Otherwise, 'MSE' is the mean */
  /*                       sum of squared residuals obtained from the model */
  /*                       with B and STATS.Intercept. */
  /*  */
  /*      If cross-validation was performed, STATS also includes the following
   */
  /*      fields: */
  /*  */
  /*        'SE'           The standard error of MSE for each lambda, as */
  /*                       calculated during cross-validation. Dimension 1xL. */
  /*        'LambdaMinMSE' The lambda value with minimum MSE. Scalar. */
  /*        'Lambda1SE'    The largest lambda such that MSE is within  */
  /*                       one standard error of the minimum. Scalar. */
  /*        'IndexMinMSE'  The index of Lambda with value LambdaMinMSE. */
  /*        'Index1SE'     The index of Lambda with value Lambda1SE. */
  /*  */
  /*      Examples: */
  /*  */
  /*         % (1) Run the lasso on data obtained from the 1985 Auto Imports
   * Database  */
  /*         % of the UCI repository.   */
  /*         %
   * http://archive.ics.uci.edu/ml/machine-learning-databases/autos/imports-85.names
   */
  /*         load imports-85; */
  /*         Description */
  /*  */
  /*         % Extract Price as the response variable and extract
   * non-categorical */
  /*         % variables related to auto construction and performance */
  /*         % */
  /*         X = X(~any(isnan(X(:,1:16)),2),:); */
  /*         Y = X(:,16); */
  /*         Y = log(Y); */
  /*         X = X(:,3:15); */
  /*         predictorNames = {'wheel-base' 'length' 'width' 'height' ... */
  /*             'curb-weight' 'engine-size' 'bore' 'stroke' 'compression-ratio'
   * ... */
  /*             'horsepower' 'peak-rpm' 'city-mpg' 'highway-mpg'}; */
  /*  */
  /*         % Compute the default sequence of lasso fits. */
  /*         [B,S] = lasso(X,Y,'CV',10,'PredictorNames',predictorNames); */
  /*  */
  /*         % Display a trace plot of the lasso fits. */
  /*         axTrace = lassoPlot(B,S); */
  /*         % Display the sequence of cross-validated predictive MSEs. */
  /*         axCV = lassoPlot(B,S,'PlotType','CV'); */
  /*         % Look at the kind of fit information returned by lasso. */
  /*         S */
  /*  */
  /*         % What variables are in the model corresponding to minimum  */
  /*         % cross-validated MSE, and in the sparsest model within one  */
  /*         % standard error of that minimum. */
  /*         minMSEModel = S.PredictorNames(B(:,S.IndexMinMSE)~=0) */
  /*         sparseModel = S.PredictorNames(B(:,S.Index1SE)~=0) */
  /*  */
  /*         % Fit the sparse model and examine residuals. */
  /*         Xplus = [ones(size(X,1),1) X]; */
  /*         fitSparse = Xplus * [S.Intercept(S.Index1SE); B(:,S.Index1SE)]; */
  /*         corr(fitSparse,Y-fitSparse) */
  /*         figure */
  /*         plot(fitSparse,Y-fitSparse,'o') */
  /*  */
  /*         % Consider a slightly richer model. A model with 6 variables may be
   * a  */
  /*         % reasonable alternative.  Find the index for a corresponding fit.
   */
  /*         df6index = min(find(S.DF==6)); */
  /*         fitDF6 = Xplus * [S.Intercept(df6index); B(:,df6index)]; */
  /*         corr(fitDF6,Y-fitDF6) */
  /*         plot(fitDF6,Y-fitDF6,'o')          */
  /*           */
  /*         % (2) Run lasso on some random data with 250 predictors */
  /*         % */
  /*         n = 1000; p = 250; */
  /*         X = randn(n,p); */
  /*         beta = randn(p,1); beta0 = randn; */
  /*         Y = beta0 + X*beta + randn(n,1); */
  /*         lambda = 0:.01:.5; */
  /*         [B,S] = lasso(X,Y,'Lambda',lambda); */
  /*         lassoPlot(B,S); */
  /*  */
  /*         % compare against OLS */
  /*         % */
  /*         figure */
  /*         bls = [ones(size(X,1),1) X] \ Y; */
  /*         plot(bls,[S.Intercept; B],'.'); */
  /*  */
  /*         % Run the same lasso fit but restricting the number of */
  /*         % non-zero coefficients in the fitted model. */
  /*         % */
  /*         [B2,S2] = lasso(X,Y,'Lambda',lambda,'DFmax',12); */
  /*  */
  /*    See also lassoPlot, ridge, parallelstats. */
  /*    References:  */
  /*    [1] Tibshirani, R. (1996) Regression shrinkage and selection */
  /*        via the lasso. Journal of the Royal Statistical Society, */
  /*        Series B, Vol 58, No. 1, pp. 267-288. */
  /*    [2] Zou, H. and T. Hastie. (2005) Regularization and variable */
  /*        selection via the elastic net. Journal of the Royal Statistical */
  /*        Society, Series B, Vol. 67, No. 2, pp. 301-320. */
  /*    [3] Friedman, J., R. Tibshirani, and T. Hastie. (2010) Regularization */
  /*        paths for generalized linear models via coordinate descent. */
  /*        Journal of Statistical Software, Vol 33, No. 1, */
  /*        http://www.jstatsoft.org/v33/i01. */
  /*    [4] Hastie, T., R. Tibshirani, and J. Friedman. (2008) The Elements */
  /*        of Statistical Learning, 2nd edition, Springer, New York. */
  /*    Copyright 2011-2016 The MathWorks, Inc. */
  /*  This screen (okrows) selects all the predictions and response we can use,
   */
  /*  but there may be further reductions for zero observation weights. */
  /*  Defer further checking of X,Y until weights are pre-screened. */
  /* 'lasso_cvmat:203' okrows = all(isfinite(X),2) & all(isfinite(Y),2); */
  outsize[0] = (unsigned int)X->size[0];
  outsize[1] = (unsigned int)X->size[1];
  b_outsize[0] = (unsigned int)Y->size[0];
  b_outsize[1] = 1U;
  emxInit_boolean_T(&okrows, 1);
  if (X->size[0] == Y->size[0]) {
    i = okrows->size[0];
    okrows->size[0] = X->size[0];
    emxEnsureCapacity_boolean_T(okrows, i);
    okrows_data = okrows->data;
    loop_ub = X->size[0];
    for (i = 0; i < loop_ub; i++) {
      okrows_data[i] = true;
    }
  } else {
    p_binary_expand_op(okrows, outsize, b_outsize);
  }
  /*  -------------------------------------------------------------------- */
  /*  Parse and process the optional parameters */
  /*  -------------------------------------------------------------------- */
  /* 'lasso_cvmat:209' LRdefault = 0.0001; */
  /* 'lasso_cvmat:210' coder.varsize('B'); */
  /* 'lasso_cvmat:211' coder.varsize('lambda'); */
  /*  **parseArgs** is not supported */
  /*  pnames = { 'weights' 'alpha' 'numlambda' 'lambdaratio' 'lambda' ... */
  /*      'dfmax' 'standardize' 'reltol' 'cv' 'mcreps' 'maxiter'... */
  /*      'predictornames' 'options' }; */
  /*  dflts  = { []        1       100       LRdefault     []      ... */
  /*       []      true          1e-4    'resubstitution'  1  1e5 ... */
  /*       {}               []}; */
  /*  [weights, alpha, nLambda, lambdaRatio, lambda, ... */
  /*      dfmax, standardize, reltol, cvp, mcreps, maxIter,predictorNames,
   * ParOptions] ... */
  /*       = internal.stats.parseArgs(pnames, dflts, varargin{:}); */
  /* 'lasso_cvmat:224' weights = []; */
  /* 'lasso_cvmat:225' alpha = 1; */
  /* 'lasso_cvmat:226' nLambda = numLambda; */
  /* 'lasso_cvmat:227' lambdaRatio = LRdefault; */
  /* 'lasso_cvmat:228' lambda = []; */
  /* 'lasso_cvmat:229' dfmax = dfMax; */
  /* 'lasso_cvmat:230' standardize = Std; */
  /* 'lasso_cvmat:231' reltol = 1e-4; */
  /* 'lasso_cvmat:232' cvp = 'resubstitution'; */
  /* 'lasso_cvmat:233' mcreps = 1; */
  /* 'lasso_cvmat:234' maxIter = 1e5; */
  /* 'lasso_cvmat:235' predictorNames = {}; */
  /* 'lasso_cvmat:236' ParOptions = []; */
  /*  Remove observations with NaNs and Infs in the predictor or response */
  /*  or with zero observation weight. */
  /* 'lasso_cvmat:240' X = X(okrows,:); */
  nx = okrows->size[0] - 1;
  loop_ub = 0;
  for (i = 0; i <= nx; i++) {
    loop_ub++;
  }
  emxInit_int32_T(&b_r, 1);
  i = b_r->size[0];
  b_r->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(b_r, i);
  r1 = b_r->data;
  for (i = 0; i <= nx; i++) {
    r1[i] = i + 1;
  }
  emxInit_real_T(&b_X, 2);
  nx = X->size[1] - 1;
  i = b_X->size[0] * b_X->size[1];
  b_X->size[0] = b_r->size[0];
  b_X->size[1] = X->size[1];
  emxEnsureCapacity_real_T(b_X, i);
  b_X_data = b_X->data;
  for (i = 0; i <= nx; i++) {
    loop_ub = b_r->size[0];
    for (n = 0; n < loop_ub; n++) {
      b_X_data[n + b_X->size[0] * i] = X_data[(r1[n] + X->size[0] * i) - 1];
    }
  }
  emxFree_int32_T(&b_r);
  i = X->size[0] * X->size[1];
  X->size[0] = b_X->size[0];
  X->size[1] = b_X->size[1];
  emxEnsureCapacity_real_T(X, i);
  X_data = X->data;
  loop_ub = b_X->size[0] * b_X->size[1];
  for (i = 0; i < loop_ub; i++) {
    X_data[i] = b_X_data[i];
  }
  emxFree_real_T(&b_X);
  /* 'lasso_cvmat:241' Y = Y(okrows); */
  nx = okrows->size[0] - 1;
  loop_ub = 0;
  emxFree_boolean_T(&okrows);
  for (i = 0; i <= nx; i++) {
    loop_ub++;
  }
  i = Y->size[0];
  Y->size[0] = loop_ub;
  emxEnsureCapacity_real_T(Y, i);
  /*  We need at least two observations after stripping NaNs and Infs. */
  /* 'lasso_cvmat:244' if size(X,1) < 2 */
  /*  If X has any constant columns, we want to exclude them from the */
  /*  coordinate descent calculations.  The corresponding coefficients */
  /*  will be returned as zero. */
  /*  constantPredictors = (range(X)==0); */
  /* 'lasso_cvmat:252' constantPredictors = ((max(X)-min(X))==0); */
  /* 'lasso_cvmat:253' ever_active = ~constantPredictors; */
  /* 'lasso_cvmat:255' [~,P] = size(X); */
  /*  === 'Standardize' option === */
  /*  Require a logical value. */
  /*  if ~isscalar(standardize) || (~islogical(standardize) && standardize~=0 &&
   * standardize~=1) */
  /*      error('stats:lasso:error', 'InvalidStandardize'); */
  /*  end */
  /*  === 'Lambda' sequence === */
  /*  lambdaMax is the penalty term (lambda) beyond which coefficients */
  /*  are guaranteed to be all zero.  If the command line does not provide */
  /*  a lambda sequence, we use lambdaMax in constructing the default  */
  /*  lambda sequence.  We always skip computation with lambda > lambdaMax */
  /*  because we know a priori that the computed coefficients will be zero. */
  /*  */
  /*  nullMSE is the mse of the fit using just a constant term. */
  /*  It is used to terminate the (ever-less penalized) fits when it becomes */
  /*  clear that we are overfitting. */
  /* 'lasso_cvmat:276' [lambdaMax, nullMSE] = computeLambdaMax(X, Y, weights,
   * alpha, standardize); */
  computeLambdaMax(X, Y, &lambdaMax, &nullMSE);
  /*  Used with nullMSE (calculated below) to terminate */
  /*  (ever-less penalized) fits when overfitting is detected. */
  /* 'lasso_cvmat:280' userSuppliedLambda = true; */
  /* 'lasso_cvmat:282' if isempty(lambda) */
  /*  Used with nullMSE (calculated below) to terminate  */
  /*  (ever-less penalized) fits when overfitting is detected. */
  /* 'lasso_cvmat:286' userSuppliedLambda = false; */
  /*  Sanity-check of 'NumLambda', should be positive integer. */
  /* 'lasso_cvmat:289' if ~isreal(nLambda) || ~isfinite(nLambda) || nLambda < 1
   */
  /* 'lasso_cvmat:291' else */
  /* 'lasso_cvmat:292' nLambda = floor(nLambda); */
  /*  Sanity-checking of LambdaRatio, should be in [0,1). */
  /* 'lasso_cvmat:296' if ~isreal(lambdaRatio) || lambdaRatio <0 || lambdaRatio
   * >= 1 */
  /* 'lasso_cvmat:300' if nLambda==1 */
  /* 'lasso_cvmat:302' else */
  /*  Fill in a number "nLambda" of smaller values, on a log scale. */
  /* 'lasso_cvmat:304' if lambdaRatio==0 */
  /* 'lasso_cvmat:307' else */
  /* 'lasso_cvmat:308' addZeroLambda = false; */
  /* 'lasso_cvmat:310' lambdaMin = lambdaMax * lambdaRatio; */
  lambdaMin = lambdaMax * 0.0001;
  /* 'lasso_cvmat:311' loghi = log(lambdaMax); */
  /* 'lasso_cvmat:312' loglo = log(lambdaMin); */
  /* 'lasso_cvmat:313' lambda = exp(linspace(loghi,loglo,nLambda)); */
  linspace(log(lambdaMax), log(lambdaMin), 100.0, x_data, x_size);
  nx = x_size[1];
  for (loop_ub = 0; loop_ub < nx; loop_ub++) {
    x_data[loop_ub] = exp(x_data[loop_ub]);
  }
  lambda_size[0] = 1;
  lambda_size[1] = x_size[1];
  loop_ub = x_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&lambda_data[0], &x_data[0], loop_ub * sizeof(double));
  }
  emxInit_real_T(&ex, 2);
  /* 'lasso_cvmat:314' if addZeroLambda */
  /* 'lasso_cvmat:316' else */
  /* 'lasso_cvmat:317' lambda(end) = lambdaMin; */
  lambda_data[99] = lambdaMin;
  /*  === 'RelTol' parameter === */
  /*  */
  /* 'lasso_cvmat:327' if ~isscalar(reltol) || ~isreal(reltol) ||
   * ~isfinite(reltol) || reltol <= 0 || reltol >= 1 */
  /*  === 'Mcreps' parameter === */
  /*  */
  /* 'lasso_cvmat:333' if ~isscalar(mcreps) || ~isreal(mcreps) ||
   * ~isfinite(mcreps) || mcreps < 1 */
  /* 'lasso_cvmat:336' mcreps = fix(mcreps); */
  /*  === 'MaxIter' parameter === */
  /* 'lasso_cvmat:339' validateattributes(maxIter, {'numeric'},... */
  /* 'lasso_cvmat:340'     {'scalar','positive','finite','integer'},... */
  /* 'lasso_cvmat:341'     mfilename,'''MaxIter'' parameter'); */
  /*  === 'CV' parameter === */
  /*  */
  /*  Cross Validation is Not Used */
  /*  if isnumeric(cvp) && isscalar(cvp) && (cvp==round(cvp)) && (0<cvp) */
  /*      % cvp is a kfold value. Create a cvpartition to pass to crossval.  */
  /*      if (cvp > size(X,1)) */
  /*          error('stats:lasso:error', 'InvalidCVforX'); */
  /*      end */
  /*      cvp = cvpartition(size(X,1),'Kfold',cvp); */
  /*  elseif isa(cvp,'cvpartition') */
  /*      if strcmpi(cvp.Type,'resubstitution') */
  /*          cvp = 'resubstitution'; */
  /*      elseif strcmpi(cvp.Type,'leaveout') */
  /*          error('stats:lasso:error', 'InvalidCVtype'); */
  /*      elseif strcmpi(cvp.Type,'holdout') && mcreps<=1 */
  /*          error('stats:lasso:error', 'InvalidMCReps'); */
  /*      end */
  /*  elseif strncmpi(cvp,'resubstitution',length(cvp)) */
  /*      % This may have been set as the default, or may have been */
  /*      % provided at the command line.  In case it's the latter, we */
  /*      % expand abbreviations. */
  /*      cvp = 'resubstitution'; */
  /*  else */
  /*      error('stats:lasso:error', 'InvalidCVtype'); */
  /*  end */
  /*  if strcmp(cvp,'resubstitution') && mcreps ~= 1 */
  /*      error('stats:lasso:error', 'InvalidMCReps'); */
  /*  end */
  /* 'lasso_cvmat:373' if isa(cvp,'cvpartition') */
  /*  === 'PredictorNames' parameter === */
  /*  */
  /*  If PredictorNames is not supplied or is supplied as empty, we just  */
  /*  leave it that way. Otherwise, confirm that it is a cell array of strings.
   */
  /*  */
  /* 'lasso_cvmat:387' if ~isempty(predictorNames) */
  /*  === 'Options' parameter === */
  /*  The 'Options' parameter is passed to crossval for handling. */
  /*  crossval will do sanity checking. */
  /*  -------------------- */
  /*  Lasso model fits */
  /*  -------------------- */
  /*  The struct 'stats' will comprise the second return argument. */
  /*  Put place holders for ever-present fields to secure the order */
  /*  we want in the struct. */
  /* 'lasso_cvmat:406' stats = struct(); */
  /* 'lasso_cvmat:407' stats.Intercept      = []; */
  /* 'lasso_cvmat:408' stats.Lambda         = []; */
  /* 'lasso_cvmat:409' stats.Alpha          = alpha; */
  /* 'lasso_cvmat:410' stats.DF             = []; */
  /* 'lasso_cvmat:411' stats.MSE            = []; */
  /* 'lasso_cvmat:412' stats.PredictorNames = predictorNames; */
  /* 'lasso_cvmat:414' [B,Intercept,lambda,mse] = ... */
  /* 'lasso_cvmat:415'     lassoFit(X,Y, ... */
  /* 'lasso_cvmat:416'
   * weights,lambda,alpha,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE,maxIter);
   */
  nx = X->size[0];
  n = X->size[1];
  i = ex->size[0] * ex->size[1];
  ex->size[0] = 1;
  ex->size[1] = X->size[1];
  emxEnsureCapacity_real_T(ex, i);
  b_X_data = ex->data;
  if (X->size[1] >= 1) {
    for (loop_ub = 0; loop_ub < n; loop_ub++) {
      b_X_data[loop_ub] = X_data[X->size[0] * loop_ub];
      for (i = 2; i <= nx; i++) {
        lambdaMin = X_data[(i + X->size[0] * loop_ub) - 1];
        if (b_X_data[loop_ub] < lambdaMin) {
          b_X_data[loop_ub] = lambdaMin;
        }
      }
    }
  }
  emxInit_real_T(&b_ex, 2);
  nx = X->size[0];
  n = X->size[1];
  i = b_ex->size[0] * b_ex->size[1];
  b_ex->size[0] = 1;
  b_ex->size[1] = X->size[1];
  emxEnsureCapacity_real_T(b_ex, i);
  ex_data = b_ex->data;
  if (X->size[1] >= 1) {
    for (loop_ub = 0; loop_ub < n; loop_ub++) {
      ex_data[loop_ub] = X_data[X->size[0] * loop_ub];
      for (i = 2; i <= nx; i++) {
        lambdaMin = X_data[(i + X->size[0] * loop_ub) - 1];
        if (ex_data[loop_ub] > lambdaMin) {
          ex_data[loop_ub] = lambdaMin;
        }
      }
    }
  }
  emxInit_real_T(&b_B, 2);
  emxInit_real_T(&mse, 2);
  if (ex->size[1] == b_ex->size[1]) {
    emxInit_boolean_T(&c_ex, 2);
    i = c_ex->size[0] * c_ex->size[1];
    c_ex->size[0] = 1;
    c_ex->size[1] = ex->size[1];
    emxEnsureCapacity_boolean_T(c_ex, i);
    okrows_data = c_ex->data;
    loop_ub = ex->size[1];
    for (i = 0; i < loop_ub; i++) {
      okrows_data[i] = (b_X_data[i] - ex_data[i] != 0.0);
    }
    lassoFit(X, Y, lambda_data, lambda_size, dfMax, lambdaMax, c_ex, nullMSE,
             b_B, x_data, x_size, mse);
    b_X_data = mse->data;
    ex_data = b_B->data;
    emxFree_boolean_T(&c_ex);
  } else {
    o_binary_expand_op(X, Y, lambda_data, lambda_size, dfMax, lambdaMax, ex,
                       b_ex, nullMSE, b_B, x_data, x_size, mse);
    b_X_data = mse->data;
    ex_data = b_B->data;
  }
  emxFree_real_T(&b_ex);
  emxFree_real_T(&ex);
  emxInit_boolean_T(&c_B, 2);
  /*  Store the number of non-zero coefficients for each lambda. */
  /* 'lasso_cvmat:419' df = sum(B~=0,1); */
  i = c_B->size[0] * c_B->size[1];
  c_B->size[0] = b_B->size[0];
  c_B->size[1] = b_B->size[1];
  emxEnsureCapacity_boolean_T(c_B, i);
  okrows_data = c_B->data;
  loop_ub = b_B->size[0] * b_B->size[1];
  for (i = 0; i < loop_ub; i++) {
    okrows_data[i] = (ex_data[i] != 0.0);
  }
  combineVectorElements(c_B, df_data, x_size);
  /*  ------------------------------------------ */
  /*  Order results by ascending lambda */
  /*  ------------------------------------------ */
  /* 'lasso_cvmat:425' nLambda = length(lambda); */
  emxFree_boolean_T(&c_B);
  if ((lambda_size[0] == 0) || (lambda_size[1] == 0)) {
    n = 0;
  } else {
    n = lambda_size[1];
  }
  /* 'lasso_cvmat:426' reverseIndices = nLambda:-1:1; */
  if (n < 1) {
    nx = 0;
  } else {
    nx = n;
    loop_ub = n - 1;
    for (i = 0; i <= loop_ub; i++) {
      reverseIndices_data[i] = n - i;
    }
  }
  /* 'lasso_cvmat:427' lambda = lambda(reverseIndices); */
  /* 'lasso_cvmat:428' lambda = reshape(lambda,1,nLambda); */
  stats_Lambda_size[0] = 1;
  stats_Lambda_size[1] = n;
  for (i = 0; i < n; i++) {
    stats_Lambda_data[i] = lambda_data[reverseIndices_data[i] - 1];
  }
  /* 'lasso_cvmat:429' B = B(:,reverseIndices); */
  loop_ub = b_B->size[0];
  i = B->size[0] * B->size[1];
  B->size[0] = b_B->size[0];
  B->size[1] = nx;
  emxEnsureCapacity_real_T(B, i);
  X_data = B->data;
  for (i = 0; i < nx; i++) {
    for (n = 0; n < loop_ub; n++) {
      X_data[n + B->size[0] * i] =
          ex_data[n + b_B->size[0] * (reverseIndices_data[i] - 1)];
    }
  }
  emxFree_real_T(&b_B);
  /* 'lasso_cvmat:430' Intercept = Intercept(reverseIndices); */
  /* 'lasso_cvmat:431' df = df(reverseIndices); */
  for (i = 0; i < nx; i++) {
    b_df_data[i] = df_data[reverseIndices_data[i] - 1];
  }
  if (0 <= nx - 1) {
    memcpy(&df_data[0], &b_df_data[0], nx * sizeof(int));
  }
  /* 'lasso_cvmat:432' mse = mse(reverseIndices); */
  /*  if ~isequal(cvp,'resubstitution') */
  /*      stats.SE          = stats.SE(reverseIndices); */
  /*      stats.IndexMinMSE = nLambda - stats.IndexMinMSE + 1; */
  /*      stats.Index1SE    = nLambda - stats.Index1SE + 1; */
  /*  end */
  /* 'lasso_cvmat:439' stats.Intercept = Intercept; */
  stats_Intercept_size[0] = 1;
  stats_Intercept_size[1] = nx;
  for (i = 0; i < nx; i++) {
    stats_Intercept_data[i] = x_data[reverseIndices_data[i] - 1];
  }
  /* 'lasso_cvmat:440' stats.Lambda = lambda; */
  /* 'lasso_cvmat:441' stats.DF = df; */
  stats_DF_size[0] = 1;
  stats_DF_size[1] = nx;
  for (i = 0; i < nx; i++) {
    stats_DF_data[i] = df_data[i];
  }
  /* 'lasso_cvmat:442' stats.MSE = mse; */
  stats_MSE_size[0] = 1;
  stats_MSE_size[1] = nx;
  for (i = 0; i < nx; i++) {
    stats_MSE_data[i] = b_X_data[reverseIndices_data[i] - 1];
  }
  emxFree_real_T(&mse);
  *stats_Alpha = 1.0;
}

/*
 * function [B,stats] = lasso_mat(X, Y, dfMax, Std, numLambda)
 */
void lasso_cvmat(emxArray_real_T *X, emxArray_real_T *Y, double dfMax,
                 emxArray_real_T *B, double stats_Intercept_data[],
                 int stats_Intercept_size[2], double stats_Lambda_data[],
                 int stats_Lambda_size[2], double *stats_Alpha,
                 double stats_DF_data[], int stats_DF_size[2],
                 double stats_MSE_data[], int stats_MSE_size[2])
{
  emxArray_boolean_T *c_B;
  emxArray_boolean_T *c_ex;
  emxArray_boolean_T *okrows;
  emxArray_int32_T *b_r;
  emxArray_real_T *b_B;
  emxArray_real_T *b_X;
  emxArray_real_T *b_ex;
  emxArray_real_T *ex;
  emxArray_real_T *mse;
  double b_lambda_data[100];
  double x_data[100];
  double lambda_data[20];
  double lambdaMax;
  double lambdaMin;
  double nullMSE;
  double *X_data;
  double *b_X_data;
  double *ex_data;
  int b_y_data[100];
  int y_data[100];
  unsigned int b_outsize[2];
  int lambda_size[2];
  unsigned int outsize[2];
  int x_size[2];
  int i;
  int lambda_size_idx_1;
  int loop_ub;
  int n;
  int nx;
  int *r1;
  bool *okrows_data;
  X_data = X->data;
  /* LASSO Perform lasso or elastic net regularization for linear regression. */
  /*    [B,STATS] = lasso(X,Y,...) Performs L1-constrained linear least   */
  /*    squares fits (lasso) or L1- and L2-constrained fits (elastic net) */
  /*    relating the predictors in X to the responses in Y. The default is a */
  /*    lasso fit, or constraint on the L1-norm of the coefficients B. */
  /*  */
  /*    Positional parameters: */
  /*  */
  /*      X                A numeric matrix (dimension, say, NxP) */
  /*      Y                A numeric vector of length N */
  /*     */
  /*    Optional input parameters:   */
  /*  */
  /*      'Weights'        Observation weights.  Must be a vector of
   * non-negative */
  /*                       values, of the same length as columns of X.  At least
   */
  /*                       two values must be positive. (default ones(N,1) or */
  /*                       equivalently (1/N)*ones(N,1)). */
  /*      'Alpha'          Elastic net mixing value, or the relative balance */
  /*                       between L2 and L1 penalty (default 1, range (0,1]).
   */
  /*                       Alpha=1 ==> lasso, otherwise elastic net. */
  /*                       Alpha near zero ==> nearly ridge regression. */
  /*      'NumLambda'      The number of lambda values to use, if the parameter
   */
  /*                       'Lambda' is not supplied (default 100).  Ignored */
  /*                       if 'Lambda' is supplied.  LASSO may return fewer */
  /*                       fits than specified by 'NumLambda' if the residual */
  /*                       error of the fits drops below a threshold percentage
   */
  /*                       of the variance of Y. */
  /*      'LambdaRatio'    Ratio between the minimum value and maximum value of
   */
  /*                       lambda to generate, if the  parameter "Lambda" is not
   */
  /*                       supplied.  Legal range is [0,1). Default is 0.0001.
   */
  /*                       If 'LambdaRatio' is zero, LASSO will generate its */
  /*                       default sequence of lambda values but replace the */
  /*                       smallest value in this sequence with the value zero.
   */
  /*                       'LambdaRatio' is ignored if 'Lambda' is supplied. */
  /*      'Lambda'         Lambda values. Will be returned in return argument */
  /*                       STATS in ascending order. The default is to have
   * LASSO */
  /*                       generate a sequence of lambda values, based on
   * 'NumLambda' */
  /*                       and 'LambdaRatio'. LASSO will generate a sequence,
   * based */
  /*                       on the values in X and Y, such that the largest
   * LAMBDA                  */
  /*                       value is just sufficient to produce all zero
   * coefficients B. */
  /*                       You may supply a vector of real, non-negative values
   * of  */
  /*                       lambda for LASSO to use, in place of its default
   * sequence. */
  /*                       If you supply a value for 'Lambda', 'NumLambda' and
   */
  /*                       'LambdaRatio' are ignored. */
  /*      'DFmax'          Maximum number of non-zero coefficients in the model.
   */
  /*                       Can be useful with large numbers of predictors. */
  /*                       Results only for lambda values that satisfy this */
  /*                       degree of sparseness will be returned. Default is */
  /*                       to not limit the number of non-zero coefficients. */
  /*      'Standardize'    Whether to scale X prior to fitting the model */
  /*                       sequence. This affects whether the regularization is
   */
  /*                       applied to the coefficients on the standardized */
  /*                       scale or the original scale. The results are always
   */
  /*                       presented on the original data scale. Default is */
  /*                       TRUE, do scale X. */
  /*                       Note: X and Y are always centered. */
  /*      'RelTol'         Convergence threshold for coordinate descent
   * algorithm. */
  /*                       The coordinate descent iterations will terminate */
  /*                       when the relative change in the size of the */
  /*                       estimated coefficients B drops below this threshold.
   */
  /*                       Default: 1e-4. Legal range is (0,1). */
  /*      'MaxIter'        Maximum number of iterations allowed.  Default is
   * 1e5. */
  /*      'PredictorNames' A string/cell array of names for the predictor
   * variables, */
  /*                       in the order in which they appear in X.  */
  /*                       Default: {} */
  /*      'Options'        A structure that contains options specifying whether
   * to */
  /*                       conduct cross-validation evaluations in parallel, and
   */
  /*                       options specifying how to use random numbers when
   * computing */
  /*                       cross validation partitions. This argument can be
   * created */
  /*                       by a call to STATSET. CROSSVAL uses the following
   * fields: */
  /*                         'UseParallel' */
  /*                         'UseSubstreams' */
  /*                         'Streams' */
  /*                       For information on these fields see PARALLELSTATS. */
  /*                       NOTE: If supplied, 'Streams' must be of length one.
   */
  /*     */
  /*    Return values: */
  /*      B                The fitted coefficients for each model.  */
  /*                       B will have dimension PxL, where  */
  /*                       P = size(X,2) is the number of predictors, and */
  /*                       L = length(lambda). */
  /*      STATS            STATS is a struct that contains information about the
   */
  /*                       sequence of model fits corresponding to the columns
   */
  /*                       of B. STATS contains the following fields: */
  /*  */
  /*        'Intercept'    The intercept term for each model. Dimension 1xL. */
  /*        'Lambda'       The sequence of lambda penalties used, in ascending
   * order.  */
  /*                       Dimension 1xL. */
  /*        'Alpha'        The elastic net mixing value that was used. */
  /*        'DF'           The number of nonzero coefficients in B for each */
  /*                       value of lambda. Dimension 1xL. */
  /*        'MSE'          The mean squared error of the fitted model for each
   */
  /*                       value of lambda. If cross-validation was performed,
   */
  /*                       the values for 'MSE' represent Mean Prediction */
  /*                       Squared Error for each value of lambda, as calculated
   */
  /*                       by cross-validation. Otherwise, 'MSE' is the mean */
  /*                       sum of squared residuals obtained from the model */
  /*                       with B and STATS.Intercept. */
  /*  */
  /*      If cross-validation was performed, STATS also includes the following
   */
  /*      fields: */
  /*  */
  /*        'SE'           The standard error of MSE for each lambda, as */
  /*                       calculated during cross-validation. Dimension 1xL. */
  /*        'LambdaMinMSE' The lambda value with minimum MSE. Scalar. */
  /*        'Lambda1SE'    The largest lambda such that MSE is within  */
  /*                       one standard error of the minimum. Scalar. */
  /*        'IndexMinMSE'  The index of Lambda with value LambdaMinMSE. */
  /*        'Index1SE'     The index of Lambda with value Lambda1SE. */
  /*  */
  /*      Examples: */
  /*  */
  /*         % (1) Run the lasso on data obtained from the 1985 Auto Imports
   * Database  */
  /*         % of the UCI repository.   */
  /*         %
   * http://archive.ics.uci.edu/ml/machine-learning-databases/autos/imports-85.names
   */
  /*         load imports-85; */
  /*         Description */
  /*  */
  /*         % Extract Price as the response variable and extract
   * non-categorical */
  /*         % variables related to auto construction and performance */
  /*         % */
  /*         X = X(~any(isnan(X(:,1:16)),2),:); */
  /*         Y = X(:,16); */
  /*         Y = log(Y); */
  /*         X = X(:,3:15); */
  /*         predictorNames = {'wheel-base' 'length' 'width' 'height' ... */
  /*             'curb-weight' 'engine-size' 'bore' 'stroke' 'compression-ratio'
   * ... */
  /*             'horsepower' 'peak-rpm' 'city-mpg' 'highway-mpg'}; */
  /*  */
  /*         % Compute the default sequence of lasso fits. */
  /*         [B,S] = lasso(X,Y,'CV',10,'PredictorNames',predictorNames); */
  /*  */
  /*         % Display a trace plot of the lasso fits. */
  /*         axTrace = lassoPlot(B,S); */
  /*         % Display the sequence of cross-validated predictive MSEs. */
  /*         axCV = lassoPlot(B,S,'PlotType','CV'); */
  /*         % Look at the kind of fit information returned by lasso. */
  /*         S */
  /*  */
  /*         % What variables are in the model corresponding to minimum  */
  /*         % cross-validated MSE, and in the sparsest model within one  */
  /*         % standard error of that minimum. */
  /*         minMSEModel = S.PredictorNames(B(:,S.IndexMinMSE)~=0) */
  /*         sparseModel = S.PredictorNames(B(:,S.Index1SE)~=0) */
  /*  */
  /*         % Fit the sparse model and examine residuals. */
  /*         Xplus = [ones(size(X,1),1) X]; */
  /*         fitSparse = Xplus * [S.Intercept(S.Index1SE); B(:,S.Index1SE)]; */
  /*         corr(fitSparse,Y-fitSparse) */
  /*         figure */
  /*         plot(fitSparse,Y-fitSparse,'o') */
  /*  */
  /*         % Consider a slightly richer model. A model with 6 variables may be
   * a  */
  /*         % reasonable alternative.  Find the index for a corresponding fit.
   */
  /*         df6index = min(find(S.DF==6)); */
  /*         fitDF6 = Xplus * [S.Intercept(df6index); B(:,df6index)]; */
  /*         corr(fitDF6,Y-fitDF6) */
  /*         plot(fitDF6,Y-fitDF6,'o')          */
  /*           */
  /*         % (2) Run lasso on some random data with 250 predictors */
  /*         % */
  /*         n = 1000; p = 250; */
  /*         X = randn(n,p); */
  /*         beta = randn(p,1); beta0 = randn; */
  /*         Y = beta0 + X*beta + randn(n,1); */
  /*         lambda = 0:.01:.5; */
  /*         [B,S] = lasso(X,Y,'Lambda',lambda); */
  /*         lassoPlot(B,S); */
  /*  */
  /*         % compare against OLS */
  /*         % */
  /*         figure */
  /*         bls = [ones(size(X,1),1) X] \ Y; */
  /*         plot(bls,[S.Intercept; B],'.'); */
  /*  */
  /*         % Run the same lasso fit but restricting the number of */
  /*         % non-zero coefficients in the fitted model. */
  /*         % */
  /*         [B2,S2] = lasso(X,Y,'Lambda',lambda,'DFmax',12); */
  /*  */
  /*    See also lassoPlot, ridge, parallelstats. */
  /*    References:  */
  /*    [1] Tibshirani, R. (1996) Regression shrinkage and selection */
  /*        via the lasso. Journal of the Royal Statistical Society, */
  /*        Series B, Vol 58, No. 1, pp. 267-288. */
  /*    [2] Zou, H. and T. Hastie. (2005) Regularization and variable */
  /*        selection via the elastic net. Journal of the Royal Statistical */
  /*        Society, Series B, Vol. 67, No. 2, pp. 301-320. */
  /*    [3] Friedman, J., R. Tibshirani, and T. Hastie. (2010) Regularization */
  /*        paths for generalized linear models via coordinate descent. */
  /*        Journal of Statistical Software, Vol 33, No. 1, */
  /*        http://www.jstatsoft.org/v33/i01. */
  /*    [4] Hastie, T., R. Tibshirani, and J. Friedman. (2008) The Elements */
  /*        of Statistical Learning, 2nd edition, Springer, New York. */
  /*    Copyright 2011-2016 The MathWorks, Inc. */
  /*  This screen (okrows) selects all the predictions and response we can use,
   */
  /*  but there may be further reductions for zero observation weights. */
  /*  Defer further checking of X,Y until weights are pre-screened. */
  /* 'lasso_cvmat:203' okrows = all(isfinite(X),2) & all(isfinite(Y),2); */
  outsize[0] = (unsigned int)X->size[0];
  outsize[1] = (unsigned int)X->size[1];
  b_outsize[0] = (unsigned int)Y->size[0];
  b_outsize[1] = 1U;
  emxInit_boolean_T(&okrows, 1);
  if (X->size[0] == Y->size[0]) {
    i = okrows->size[0];
    okrows->size[0] = X->size[0];
    emxEnsureCapacity_boolean_T(okrows, i);
    okrows_data = okrows->data;
    loop_ub = X->size[0];
    for (i = 0; i < loop_ub; i++) {
      okrows_data[i] = true;
    }
  } else {
    p_binary_expand_op(okrows, outsize, b_outsize);
  }
  /*  -------------------------------------------------------------------- */
  /*  Parse and process the optional parameters */
  /*  -------------------------------------------------------------------- */
  /* 'lasso_cvmat:209' LRdefault = 0.0001; */
  /* 'lasso_cvmat:210' coder.varsize('B'); */
  /* 'lasso_cvmat:211' coder.varsize('lambda'); */
  /*  **parseArgs** is not supported */
  /*  pnames = { 'weights' 'alpha' 'numlambda' 'lambdaratio' 'lambda' ... */
  /*      'dfmax' 'standardize' 'reltol' 'cv' 'mcreps' 'maxiter'... */
  /*      'predictornames' 'options' }; */
  /*  dflts  = { []        1       100       LRdefault     []      ... */
  /*       []      true          1e-4    'resubstitution'  1  1e5 ... */
  /*       {}               []}; */
  /*  [weights, alpha, nLambda, lambdaRatio, lambda, ... */
  /*      dfmax, standardize, reltol, cvp, mcreps, maxIter,predictorNames,
   * ParOptions] ... */
  /*       = internal.stats.parseArgs(pnames, dflts, varargin{:}); */
  /* 'lasso_cvmat:224' weights = []; */
  /* 'lasso_cvmat:225' alpha = 1; */
  /* 'lasso_cvmat:226' nLambda = numLambda; */
  /* 'lasso_cvmat:227' lambdaRatio = LRdefault; */
  /* 'lasso_cvmat:228' lambda = []; */
  /* 'lasso_cvmat:229' dfmax = dfMax; */
  /* 'lasso_cvmat:230' standardize = Std; */
  /* 'lasso_cvmat:231' reltol = 1e-4; */
  /* 'lasso_cvmat:232' cvp = 'resubstitution'; */
  /* 'lasso_cvmat:233' mcreps = 1; */
  /* 'lasso_cvmat:234' maxIter = 1e5; */
  /* 'lasso_cvmat:235' predictorNames = {}; */
  /* 'lasso_cvmat:236' ParOptions = []; */
  /*  Remove observations with NaNs and Infs in the predictor or response */
  /*  or with zero observation weight. */
  /* 'lasso_cvmat:240' X = X(okrows,:); */
  nx = okrows->size[0] - 1;
  loop_ub = 0;
  for (i = 0; i <= nx; i++) {
    loop_ub++;
  }
  emxInit_int32_T(&b_r, 1);
  i = b_r->size[0];
  b_r->size[0] = loop_ub;
  emxEnsureCapacity_int32_T(b_r, i);
  r1 = b_r->data;
  for (i = 0; i <= nx; i++) {
    r1[i] = i + 1;
  }
  emxInit_real_T(&b_X, 2);
  nx = X->size[1] - 1;
  i = b_X->size[0] * b_X->size[1];
  b_X->size[0] = b_r->size[0];
  b_X->size[1] = X->size[1];
  emxEnsureCapacity_real_T(b_X, i);
  b_X_data = b_X->data;
  for (i = 0; i <= nx; i++) {
    loop_ub = b_r->size[0];
    for (lambda_size_idx_1 = 0; lambda_size_idx_1 < loop_ub;
         lambda_size_idx_1++) {
      b_X_data[lambda_size_idx_1 + b_X->size[0] * i] =
          X_data[(r1[lambda_size_idx_1] + X->size[0] * i) - 1];
    }
  }
  emxFree_int32_T(&b_r);
  i = X->size[0] * X->size[1];
  X->size[0] = b_X->size[0];
  X->size[1] = b_X->size[1];
  emxEnsureCapacity_real_T(X, i);
  X_data = X->data;
  loop_ub = b_X->size[0] * b_X->size[1];
  for (i = 0; i < loop_ub; i++) {
    X_data[i] = b_X_data[i];
  }
  emxFree_real_T(&b_X);
  /* 'lasso_cvmat:241' Y = Y(okrows); */
  nx = okrows->size[0] - 1;
  loop_ub = 0;
  emxFree_boolean_T(&okrows);
  for (i = 0; i <= nx; i++) {
    loop_ub++;
  }
  i = Y->size[0];
  Y->size[0] = loop_ub;
  emxEnsureCapacity_real_T(Y, i);
  /*  We need at least two observations after stripping NaNs and Infs. */
  /* 'lasso_cvmat:244' if size(X,1) < 2 */
  /*  If X has any constant columns, we want to exclude them from the */
  /*  coordinate descent calculations.  The corresponding coefficients */
  /*  will be returned as zero. */
  /*  constantPredictors = (range(X)==0); */
  /* 'lasso_cvmat:252' constantPredictors = ((max(X)-min(X))==0); */
  /* 'lasso_cvmat:253' ever_active = ~constantPredictors; */
  /* 'lasso_cvmat:255' [~,P] = size(X); */
  /*  === 'Standardize' option === */
  /*  Require a logical value. */
  /*  if ~isscalar(standardize) || (~islogical(standardize) && standardize~=0 &&
   * standardize~=1) */
  /*      error('stats:lasso:error', 'InvalidStandardize'); */
  /*  end */
  /*  === 'Lambda' sequence === */
  /*  lambdaMax is the penalty term (lambda) beyond which coefficients */
  /*  are guaranteed to be all zero.  If the command line does not provide */
  /*  a lambda sequence, we use lambdaMax in constructing the default  */
  /*  lambda sequence.  We always skip computation with lambda > lambdaMax */
  /*  because we know a priori that the computed coefficients will be zero. */
  /*  */
  /*  nullMSE is the mse of the fit using just a constant term. */
  /*  It is used to terminate the (ever-less penalized) fits when it becomes */
  /*  clear that we are overfitting. */
  /* 'lasso_cvmat:276' [lambdaMax, nullMSE] = computeLambdaMax(X, Y, weights,
   * alpha, standardize); */
  computeLambdaMax(X, Y, &lambdaMax, &nullMSE);
  /*  Used with nullMSE (calculated below) to terminate */
  /*  (ever-less penalized) fits when overfitting is detected. */
  /* 'lasso_cvmat:280' userSuppliedLambda = true; */
  /* 'lasso_cvmat:282' if isempty(lambda) */
  /*  Used with nullMSE (calculated below) to terminate  */
  /*  (ever-less penalized) fits when overfitting is detected. */
  /* 'lasso_cvmat:286' userSuppliedLambda = false; */
  /*  Sanity-check of 'NumLambda', should be positive integer. */
  /* 'lasso_cvmat:289' if ~isreal(nLambda) || ~isfinite(nLambda) || nLambda < 1
   */
  /* 'lasso_cvmat:291' else */
  /* 'lasso_cvmat:292' nLambda = floor(nLambda); */
  /*  Sanity-checking of LambdaRatio, should be in [0,1). */
  /* 'lasso_cvmat:296' if ~isreal(lambdaRatio) || lambdaRatio <0 || lambdaRatio
   * >= 1 */
  /* 'lasso_cvmat:300' if nLambda==1 */
  /* 'lasso_cvmat:302' else */
  /*  Fill in a number "nLambda" of smaller values, on a log scale. */
  /* 'lasso_cvmat:304' if lambdaRatio==0 */
  /* 'lasso_cvmat:307' else */
  /* 'lasso_cvmat:308' addZeroLambda = false; */
  /* 'lasso_cvmat:310' lambdaMin = lambdaMax * lambdaRatio; */
  lambdaMin = lambdaMax * 0.0001;
  /* 'lasso_cvmat:311' loghi = log(lambdaMax); */
  /* 'lasso_cvmat:312' loglo = log(lambdaMin); */
  /* 'lasso_cvmat:313' lambda = exp(linspace(loghi,loglo,nLambda)); */
  linspace(log(lambdaMax), log(lambdaMin), 20.0, x_data, x_size);
  nx = x_size[1];
  for (loop_ub = 0; loop_ub < nx; loop_ub++) {
    x_data[loop_ub] = exp(x_data[loop_ub]);
  }
  lambda_size_idx_1 = x_size[1];
  loop_ub = x_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&lambda_data[0], &x_data[0], loop_ub * sizeof(double));
  }
  emxInit_real_T(&ex, 2);
  /* 'lasso_cvmat:314' if addZeroLambda */
  /* 'lasso_cvmat:316' else */
  /* 'lasso_cvmat:317' lambda(end) = lambdaMin; */
  lambda_data[19] = lambdaMin;
  /*  === 'RelTol' parameter === */
  /*  */
  /* 'lasso_cvmat:327' if ~isscalar(reltol) || ~isreal(reltol) ||
   * ~isfinite(reltol) || reltol <= 0 || reltol >= 1 */
  /*  === 'Mcreps' parameter === */
  /*  */
  /* 'lasso_cvmat:333' if ~isscalar(mcreps) || ~isreal(mcreps) ||
   * ~isfinite(mcreps) || mcreps < 1 */
  /* 'lasso_cvmat:336' mcreps = fix(mcreps); */
  /*  === 'MaxIter' parameter === */
  /* 'lasso_cvmat:339' validateattributes(maxIter, {'numeric'},... */
  /* 'lasso_cvmat:340'     {'scalar','positive','finite','integer'},... */
  /* 'lasso_cvmat:341'     mfilename,'''MaxIter'' parameter'); */
  /*  === 'CV' parameter === */
  /*  */
  /*  Cross Validation is Not Used */
  /*  if isnumeric(cvp) && isscalar(cvp) && (cvp==round(cvp)) && (0<cvp) */
  /*      % cvp is a kfold value. Create a cvpartition to pass to crossval.  */
  /*      if (cvp > size(X,1)) */
  /*          error('stats:lasso:error', 'InvalidCVforX'); */
  /*      end */
  /*      cvp = cvpartition(size(X,1),'Kfold',cvp); */
  /*  elseif isa(cvp,'cvpartition') */
  /*      if strcmpi(cvp.Type,'resubstitution') */
  /*          cvp = 'resubstitution'; */
  /*      elseif strcmpi(cvp.Type,'leaveout') */
  /*          error('stats:lasso:error', 'InvalidCVtype'); */
  /*      elseif strcmpi(cvp.Type,'holdout') && mcreps<=1 */
  /*          error('stats:lasso:error', 'InvalidMCReps'); */
  /*      end */
  /*  elseif strncmpi(cvp,'resubstitution',length(cvp)) */
  /*      % This may have been set as the default, or may have been */
  /*      % provided at the command line.  In case it's the latter, we */
  /*      % expand abbreviations. */
  /*      cvp = 'resubstitution'; */
  /*  else */
  /*      error('stats:lasso:error', 'InvalidCVtype'); */
  /*  end */
  /*  if strcmp(cvp,'resubstitution') && mcreps ~= 1 */
  /*      error('stats:lasso:error', 'InvalidMCReps'); */
  /*  end */
  /* 'lasso_cvmat:373' if isa(cvp,'cvpartition') */
  /*  === 'PredictorNames' parameter === */
  /*  */
  /*  If PredictorNames is not supplied or is supplied as empty, we just  */
  /*  leave it that way. Otherwise, confirm that it is a cell array of strings.
   */
  /*  */
  /* 'lasso_cvmat:387' if ~isempty(predictorNames) */
  /*  === 'Options' parameter === */
  /*  The 'Options' parameter is passed to crossval for handling. */
  /*  crossval will do sanity checking. */
  /*  -------------------- */
  /*  Lasso model fits */
  /*  -------------------- */
  /*  The struct 'stats' will comprise the second return argument. */
  /*  Put place holders for ever-present fields to secure the order */
  /*  we want in the struct. */
  /* 'lasso_cvmat:406' stats = struct(); */
  /* 'lasso_cvmat:407' stats.Intercept      = []; */
  /* 'lasso_cvmat:408' stats.Lambda         = []; */
  /* 'lasso_cvmat:409' stats.Alpha          = alpha; */
  /* 'lasso_cvmat:410' stats.DF             = []; */
  /* 'lasso_cvmat:411' stats.MSE            = []; */
  /* 'lasso_cvmat:412' stats.PredictorNames = predictorNames; */
  /* 'lasso_cvmat:414' [B,Intercept,lambda,mse] = ... */
  /* 'lasso_cvmat:415'     lassoFit(X,Y, ... */
  /* 'lasso_cvmat:416'
   * weights,lambda,alpha,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE,maxIter);
   */
  nx = X->size[0];
  n = X->size[1];
  i = ex->size[0] * ex->size[1];
  ex->size[0] = 1;
  ex->size[1] = X->size[1];
  emxEnsureCapacity_real_T(ex, i);
  b_X_data = ex->data;
  if (X->size[1] >= 1) {
    for (loop_ub = 0; loop_ub < n; loop_ub++) {
      b_X_data[loop_ub] = X_data[X->size[0] * loop_ub];
      for (i = 2; i <= nx; i++) {
        lambdaMin = X_data[(i + X->size[0] * loop_ub) - 1];
        if (b_X_data[loop_ub] < lambdaMin) {
          b_X_data[loop_ub] = lambdaMin;
        }
      }
    }
  }
  emxInit_real_T(&b_ex, 2);
  nx = X->size[0];
  n = X->size[1];
  i = b_ex->size[0] * b_ex->size[1];
  b_ex->size[0] = 1;
  b_ex->size[1] = X->size[1];
  emxEnsureCapacity_real_T(b_ex, i);
  ex_data = b_ex->data;
  if (X->size[1] >= 1) {
    for (loop_ub = 0; loop_ub < n; loop_ub++) {
      ex_data[loop_ub] = X_data[X->size[0] * loop_ub];
      for (i = 2; i <= nx; i++) {
        lambdaMin = X_data[(i + X->size[0] * loop_ub) - 1];
        if (ex_data[loop_ub] > lambdaMin) {
          ex_data[loop_ub] = lambdaMin;
        }
      }
    }
  }
  lambda_size[0] = 1;
  lambda_size[1] = x_size[1];
  if (0 <= lambda_size_idx_1 - 1) {
    memcpy(&b_lambda_data[0], &lambda_data[0],
           lambda_size_idx_1 * sizeof(double));
  }
  emxInit_real_T(&mse, 2);
  emxInit_real_T(&b_B, 2);
  if (ex->size[1] == b_ex->size[1]) {
    emxInit_boolean_T(&c_ex, 2);
    i = c_ex->size[0] * c_ex->size[1];
    c_ex->size[0] = 1;
    c_ex->size[1] = ex->size[1];
    emxEnsureCapacity_boolean_T(c_ex, i);
    okrows_data = c_ex->data;
    loop_ub = ex->size[1];
    for (i = 0; i < loop_ub; i++) {
      okrows_data[i] = (b_X_data[i] - ex_data[i] != 0.0);
    }
    lassoFit(X, Y, b_lambda_data, lambda_size, dfMax, lambdaMax, c_ex, nullMSE,
             b_B, x_data, x_size, mse);
    b_X_data = mse->data;
    ex_data = b_B->data;
    emxFree_boolean_T(&c_ex);
  } else {
    o_binary_expand_op(X, Y, b_lambda_data, lambda_size, dfMax, lambdaMax, ex,
                       b_ex, nullMSE, b_B, x_data, x_size, mse);
    b_X_data = mse->data;
    ex_data = b_B->data;
  }
  emxFree_real_T(&b_ex);
  emxFree_real_T(&ex);
  emxInit_boolean_T(&c_B, 2);
  /*  Store the number of non-zero coefficients for each lambda. */
  /* 'lasso_cvmat:419' df = sum(B~=0,1); */
  i = c_B->size[0] * c_B->size[1];
  c_B->size[0] = b_B->size[0];
  c_B->size[1] = b_B->size[1];
  emxEnsureCapacity_boolean_T(c_B, i);
  okrows_data = c_B->data;
  loop_ub = b_B->size[0] * b_B->size[1];
  for (i = 0; i < loop_ub; i++) {
    okrows_data[i] = (ex_data[i] != 0.0);
  }
  combineVectorElements(c_B, y_data, x_size);
  /*  ------------------------------------------ */
  /*  Order results by ascending lambda */
  /*  ------------------------------------------ */
  /* 'lasso_cvmat:425' nLambda = length(lambda); */
  emxFree_boolean_T(&c_B);
  if ((lambda_size[0] == 0) || (lambda_size[1] == 0)) {
    n = 0;
  } else {
    n = lambda_size[1];
  }
  /* 'lasso_cvmat:426' reverseIndices = nLambda:-1:1; */
  if (n < 1) {
    nx = 0;
  } else {
    nx = n;
    loop_ub = n - 1;
    for (i = 0; i <= loop_ub; i++) {
      b_y_data[i] = n - i;
    }
  }
  /* 'lasso_cvmat:427' lambda = lambda(reverseIndices); */
  /* 'lasso_cvmat:428' lambda = reshape(lambda,1,nLambda); */
  /* 'lasso_cvmat:429' B = B(:,reverseIndices); */
  loop_ub = b_B->size[0];
  i = B->size[0] * B->size[1];
  B->size[0] = b_B->size[0];
  B->size[1] = nx;
  emxEnsureCapacity_real_T(B, i);
  X_data = B->data;
  for (i = 0; i < nx; i++) {
    for (lambda_size_idx_1 = 0; lambda_size_idx_1 < loop_ub;
         lambda_size_idx_1++) {
      X_data[lambda_size_idx_1 + B->size[0] * i] =
          ex_data[lambda_size_idx_1 + b_B->size[0] * (b_y_data[i] - 1)];
    }
  }
  emxFree_real_T(&b_B);
  /* 'lasso_cvmat:430' Intercept = Intercept(reverseIndices); */
  /* 'lasso_cvmat:431' df = df(reverseIndices); */
  /* 'lasso_cvmat:432' mse = mse(reverseIndices); */
  /*  if ~isequal(cvp,'resubstitution') */
  /*      stats.SE          = stats.SE(reverseIndices); */
  /*      stats.IndexMinMSE = nLambda - stats.IndexMinMSE + 1; */
  /*      stats.Index1SE    = nLambda - stats.Index1SE + 1; */
  /*  end */
  /* 'lasso_cvmat:439' stats.Intercept = Intercept; */
  stats_Intercept_size[0] = 1;
  stats_Intercept_size[1] = nx;
  for (i = 0; i < nx; i++) {
    stats_Intercept_data[i] = x_data[b_y_data[i] - 1];
  }
  /* 'lasso_cvmat:440' stats.Lambda = lambda; */
  stats_Lambda_size[0] = 1;
  stats_Lambda_size[1] = n;
  for (i = 0; i < n; i++) {
    stats_Lambda_data[i] = b_lambda_data[b_y_data[i] - 1];
  }
  /* 'lasso_cvmat:441' stats.DF = df; */
  stats_DF_size[0] = 1;
  stats_DF_size[1] = nx;
  for (i = 0; i < nx; i++) {
    stats_DF_data[i] = y_data[b_y_data[i] - 1];
  }
  /* 'lasso_cvmat:442' stats.MSE = mse; */
  stats_MSE_size[0] = 1;
  stats_MSE_size[1] = nx;
  for (i = 0; i < nx; i++) {
    stats_MSE_data[i] = b_X_data[b_y_data[i] - 1];
  }
  emxFree_real_T(&mse);
  *stats_Alpha = 1.0;
}

/* End of code generation (lasso_cvmat.c) */
