/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * eml_rand_mt19937ar_stateful.c
 *
 * Code generation for function 'eml_rand_mt19937ar_stateful'
 *
 */

/* Include files */
#include "eml_rand_mt19937ar_stateful.h"
#include "gkmPWMlasso_data.h"
#include <string.h>

/* Function Definitions */
/*
 *
 */
void c_eml_rand_mt19937ar_stateful_i(void)
{
  unsigned int b_r;
  int mti;
  memset(&state[0], 0, 625U * sizeof(unsigned int));
  b_r = 5489U;
  state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    b_r = ((b_r ^ b_r >> 30U) * 1812433253U + mti) + 1U;
    state[mti + 1] = b_r;
  }
  state[624] = 624U;
}

/* End of code generation (eml_rand_mt19937ar_stateful.c) */
