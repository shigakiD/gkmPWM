/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rand.c
 *
 * Code generation for function 'rand'
 *
 */

/* Include files */
#include "rand.h"
#include "eml_rand_mt19937ar.h"
#include "gkmPWM_data.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
double b_rand(void)
{
  return eml_rand_mt19937ar(state);
}

/*
 *
 */
void c_rand(double r[4])
{
  int k;
  for (k = 0; k < 4; k++) {
    r[k] = eml_rand_mt19937ar(state);
  }
}

/* End of code generation (rand.c) */
