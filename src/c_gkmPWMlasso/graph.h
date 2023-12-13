/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * graph.h
 *
 * Code generation for function 'graph'
 *
 */

#ifndef GRAPH_H
#define GRAPH_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include "omp.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void graph_graph(const emxArray_boolean_T *varargin_1,
                 emxArray_int32_T *G_Underlying_Ir,
                 emxArray_int32_T *G_Underlying_Jc);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (graph.h) */
