/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWMlasso_types.h
 *
 * Code generation for function 'gkmPWMlasso'
 *
 */

#ifndef GKMPWMLASSO_TYPES_H
#define GKMPWMLASSO_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T
typedef struct {
  int addr;
  int next;
  int prev;
} struct_T;
#endif /* typedef_struct_T */

#ifndef typedef_cell_wrap_14
#define typedef_cell_wrap_14
typedef struct {
  double f1[16];
} cell_wrap_14;
#endif /* typedef_cell_wrap_14 */

#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T
struct emxArray_char_T {
  char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};
#endif /* struct_emxArray_char_T */
#ifndef typedef_emxArray_char_T
#define typedef_emxArray_char_T
typedef struct emxArray_char_T emxArray_char_T;
#endif /* typedef_emxArray_char_T */

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
struct emxArray_real_T {
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};
#endif /* struct_emxArray_real_T */
#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T
typedef struct emxArray_real_T emxArray_real_T;
#endif /* typedef_emxArray_real_T */

#ifndef typedef_cell_wrap_0
#define typedef_cell_wrap_0
typedef struct {
  emxArray_real_T *f1;
} cell_wrap_0;
#endif /* typedef_cell_wrap_0 */

#ifndef typedef_emxArray_cell_wrap_0
#define typedef_emxArray_cell_wrap_0
typedef struct {
  cell_wrap_0 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_0;
#endif /* typedef_emxArray_cell_wrap_0 */

#ifndef typedef_cell_wrap_1
#define typedef_cell_wrap_1
typedef struct {
  emxArray_real_T *f1;
} cell_wrap_1;
#endif /* typedef_cell_wrap_1 */

#ifndef typedef_emxArray_cell_wrap_1
#define typedef_emxArray_cell_wrap_1
typedef struct {
  cell_wrap_1 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_1;
#endif /* typedef_emxArray_cell_wrap_1 */

#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3
typedef struct {
  emxArray_real_T *f1;
} cell_wrap_3;
#endif /* typedef_cell_wrap_3 */

#ifndef typedef_emxArray_cell_wrap_3_1x19
#define typedef_emxArray_cell_wrap_3_1x19
typedef struct {
  cell_wrap_3 data[19];
  int size[2];
} emxArray_cell_wrap_3_1x19;
#endif /* typedef_emxArray_cell_wrap_3_1x19 */

#ifndef typedef_emxArray_cell_wrap_3_20
#define typedef_emxArray_cell_wrap_3_20
typedef struct {
  cell_wrap_3 data[20];
  int size[1];
} emxArray_cell_wrap_3_20;
#endif /* typedef_emxArray_cell_wrap_3_20 */

#ifndef typedef_emxArray_cell_wrap_3
#define typedef_emxArray_cell_wrap_3
typedef struct {
  cell_wrap_3 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_3;
#endif /* typedef_emxArray_cell_wrap_3 */

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T
struct emxArray_int32_T {
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};
#endif /* struct_emxArray_int32_T */
#ifndef typedef_emxArray_int32_T
#define typedef_emxArray_int32_T
typedef struct emxArray_int32_T emxArray_int32_T;
#endif /* typedef_emxArray_int32_T */

#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T
struct emxArray_boolean_T {
  bool *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};
#endif /* struct_emxArray_boolean_T */
#ifndef typedef_emxArray_boolean_T
#define typedef_emxArray_boolean_T
typedef struct emxArray_boolean_T emxArray_boolean_T;
#endif /* typedef_emxArray_boolean_T */

#ifndef struct_emxArray_uint32_T
#define struct_emxArray_uint32_T
struct emxArray_uint32_T {
  unsigned int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};
#endif /* struct_emxArray_uint32_T */
#ifndef typedef_emxArray_uint32_T
#define typedef_emxArray_uint32_T
typedef struct emxArray_uint32_T emxArray_uint32_T;
#endif /* typedef_emxArray_uint32_T */

#ifndef typedef_cell_wrap_10
#define typedef_cell_wrap_10
typedef struct {
  emxArray_char_T *f1;
} cell_wrap_10;
#endif /* typedef_cell_wrap_10 */

#ifndef typedef_cell_wrap_11
#define typedef_cell_wrap_11
typedef struct {
  cell_wrap_10 f1[1];
} cell_wrap_11;
#endif /* typedef_cell_wrap_11 */

#ifndef typedef_emxArray_cell_wrap_11
#define typedef_emxArray_cell_wrap_11
typedef struct {
  cell_wrap_11 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_11;
#endif /* typedef_emxArray_cell_wrap_11 */

#ifndef typedef_emxArray_cell_wrap_14
#define typedef_emxArray_cell_wrap_14
typedef struct {
  cell_wrap_14 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_14;
#endif /* typedef_emxArray_cell_wrap_14 */

#ifndef typedef_emxArray_struct_T
#define typedef_emxArray_struct_T
typedef struct {
  struct_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_struct_T;
#endif /* typedef_emxArray_struct_T */

#ifndef typedef_coder_internal_list
#define typedef_coder_internal_list
typedef struct {
  emxArray_struct_T *nodePool;
  emxArray_int32_T *valuePool;
  int unusedAddr;
  int frontAddr;
  int backAddr;
  int len;
} coder_internal_list;
#endif /* typedef_coder_internal_list */

#ifndef typedef_emxArray_cell_wrap_10
#define typedef_emxArray_cell_wrap_10
typedef struct {
  cell_wrap_10 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_10;
#endif /* typedef_emxArray_cell_wrap_10 */

#endif
/* End of code generation (gkmPWMlasso_types.h) */
