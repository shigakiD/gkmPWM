/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso_emxutil.h
 *
 * Code generation for function 'gkmPWMlasso_emxutil'
 *
 */

#ifndef GKMPWMLASSO_EMXUTIL_H
#define GKMPWMLASSO_EMXUTIL_H

/* Include files */
#include "gkmPWMlasso_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void c_emxFreeStruct_coder_internal_(coder_internal_list *pStruct);

extern void c_emxInitStruct_coder_internal_(coder_internal_list *pStruct);

extern void emxCopyStruct_cell_wrap_10(cell_wrap_10 *dst,
                                       const cell_wrap_10 *src);

extern void emxCopy_char_T(emxArray_char_T **dst, emxArray_char_T *const *src);

extern void emxEnsureCapacity_boolean_T(emxArray_boolean_T *emxArray,
                                        int oldNumel);

extern void emxEnsureCapacity_cell_wrap_0(emxArray_cell_wrap_0 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_1(emxArray_cell_wrap_1 *emxArray,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_10(emxArray_cell_wrap_10 *emxArray,
                                           int oldNumel);

extern void emxEnsureCapacity_cell_wrap_11(emxArray_cell_wrap_11 *emxArray,
                                           int oldNumel);

extern void emxEnsureCapacity_cell_wrap_14(emxArray_cell_wrap_14 *emxArray,
                                           int oldNumel);

extern void emxEnsureCapacity_cell_wrap_3(cell_wrap_3 data[20], int size,
                                          int oldNumel);

extern void emxEnsureCapacity_cell_wrap_31(cell_wrap_3 data[19],
                                           const int size[2], int oldNumel);

extern void emxEnsureCapacity_cell_wrap_32(emxArray_cell_wrap_3 *emxArray,
                                           int oldNumel);

extern void emxEnsureCapacity_char_T(emxArray_char_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);

extern void emxEnsureCapacity_struct_T(emxArray_struct_T *emxArray,
                                       int oldNumel);

extern void emxEnsureCapacity_uint32_T(emxArray_uint32_T *emxArray,
                                       int oldNumel);

extern void emxExpand_cell_wrap_0(emxArray_cell_wrap_0 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_1(emxArray_cell_wrap_1 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_10(emxArray_cell_wrap_10 *emxArray,
                                   int fromIndex, int toIndex);

extern void emxExpand_cell_wrap_11(emxArray_cell_wrap_11 *emxArray,
                                   int fromIndex, int toIndex);

extern void emxExpand_cell_wrap_3(emxArray_cell_wrap_3 *emxArray, int fromIndex,
                                  int toIndex);

extern void emxExpand_cell_wrap_3_1x19(cell_wrap_3 data[19], int fromIndex,
                                       int toIndex);

extern void emxExpand_cell_wrap_3_20(cell_wrap_3 data[20], int fromIndex,
                                     int toIndex);

extern void emxFreeMatrix_cell_wrap_0(cell_wrap_0 pMatrix[8]);

extern void emxFreeMatrix_cell_wrap_10(cell_wrap_10 *pMatrix);

extern void emxFreeMatrix_cell_wrap_101(cell_wrap_10 pMatrix[3]);

extern void emxFreeStruct_cell_wrap_0(cell_wrap_0 *pStruct);

extern void emxFreeStruct_cell_wrap_1(cell_wrap_1 *pStruct);

extern void emxFreeStruct_cell_wrap_10(cell_wrap_10 *pStruct);

extern void emxFreeStruct_cell_wrap_11(cell_wrap_11 *pStruct);

extern void emxFreeStruct_cell_wrap_3(cell_wrap_3 *pStruct);

extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);

extern void emxFree_cell_wrap_0(emxArray_cell_wrap_0 **pEmxArray);

extern void emxFree_cell_wrap_1(emxArray_cell_wrap_1 **pEmxArray);

extern void emxFree_cell_wrap_10(emxArray_cell_wrap_10 **pEmxArray);

extern void emxFree_cell_wrap_11(emxArray_cell_wrap_11 **pEmxArray);

extern void emxFree_cell_wrap_14(emxArray_cell_wrap_14 **pEmxArray);

extern void emxFree_cell_wrap_3(emxArray_cell_wrap_3 **pEmxArray);

extern void emxFree_cell_wrap_3_1x19(emxArray_cell_wrap_3_1x19 *pEmxArray);

extern void emxFree_cell_wrap_3_20(emxArray_cell_wrap_3_20 *pEmxArray);

extern void emxFree_char_T(emxArray_char_T **pEmxArray);

extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);

extern void emxFree_real_T(emxArray_real_T **pEmxArray);

extern void emxFree_struct_T(emxArray_struct_T **pEmxArray);

extern void emxFree_uint32_T(emxArray_uint32_T **pEmxArray);

extern void emxInitMatrix_cell_wrap_0(cell_wrap_0 pMatrix[8]);

extern void emxInitMatrix_cell_wrap_10(cell_wrap_10 *pMatrix);

extern void emxInitMatrix_cell_wrap_101(cell_wrap_10 pMatrix[3]);

extern void emxInitStruct_cell_wrap_0(cell_wrap_0 *pStruct);

extern void emxInitStruct_cell_wrap_1(cell_wrap_1 *pStruct);

extern void emxInitStruct_cell_wrap_10(cell_wrap_10 *pStruct);

extern void emxInitStruct_cell_wrap_11(cell_wrap_11 *pStruct);

extern void emxInitStruct_cell_wrap_3(cell_wrap_3 *pStruct);

extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray,
                              int numDimensions);

extern void emxInit_cell_wrap_0(emxArray_cell_wrap_0 **pEmxArray);

extern void emxInit_cell_wrap_1(emxArray_cell_wrap_1 **pEmxArray);

extern void emxInit_cell_wrap_10(emxArray_cell_wrap_10 **pEmxArray);

extern void emxInit_cell_wrap_11(emxArray_cell_wrap_11 **pEmxArray);

extern void emxInit_cell_wrap_14(emxArray_cell_wrap_14 **pEmxArray);

extern void emxInit_cell_wrap_3(emxArray_cell_wrap_3 **pEmxArray,
                                int numDimensions);

extern void emxInit_cell_wrap_3_1x19(emxArray_cell_wrap_3_1x19 *pEmxArray);

extern void emxInit_cell_wrap_3_20(emxArray_cell_wrap_3_20 *pEmxArray);

extern void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);

extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);

extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

extern void emxInit_struct_T(emxArray_struct_T **pEmxArray, int numDimensions);

extern void emxInit_uint32_T(emxArray_uint32_T **pEmxArray);

extern void emxTrim_cell_wrap_0(emxArray_cell_wrap_0 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_1(emxArray_cell_wrap_1 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_10(emxArray_cell_wrap_10 *emxArray, int fromIndex,
                                 int toIndex);

extern void emxTrim_cell_wrap_11(emxArray_cell_wrap_11 *emxArray, int fromIndex,
                                 int toIndex);

extern void emxTrim_cell_wrap_3(emxArray_cell_wrap_3 *emxArray, int fromIndex,
                                int toIndex);

extern void emxTrim_cell_wrap_3_1x19(cell_wrap_3 data[19], int fromIndex,
                                     int toIndex);

extern void emxTrim_cell_wrap_3_20(cell_wrap_3 data[20], int fromIndex,
                                   int toIndex);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (gkmPWMlasso_emxutil.h) */
