/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * gkmPWM_types.h
 *
 * Code generation for function 'gkmPWM'
 *
 */

#ifndef GKMPWM_TYPES_H
#define GKMPWM_TYPES_H

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_cell_wrap_13
#define typedef_cell_wrap_13
typedef struct {
  double f1[16];
} cell_wrap_13;
#endif /* typedef_cell_wrap_13 */

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
  emxArray_char_T *f1;
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

#ifndef typedef_cell_wrap_2
#define typedef_cell_wrap_2
typedef struct {
  emxArray_real_T *f1;
} cell_wrap_2;
#endif /* typedef_cell_wrap_2 */

#ifndef typedef_emxArray_cell_wrap_2
#define typedef_emxArray_cell_wrap_2
typedef struct {
  cell_wrap_2 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_2;
#endif /* typedef_emxArray_cell_wrap_2 */

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

#ifndef typedef_emxArray_cell_wrap_13
#define typedef_emxArray_cell_wrap_13
typedef struct {
  cell_wrap_13 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_13;
#endif /* typedef_emxArray_cell_wrap_13 */

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

#ifndef typedef_cell_wrap_4
#define typedef_cell_wrap_4
typedef struct {
  emxArray_real_T *f1;
} cell_wrap_4;
#endif /* typedef_cell_wrap_4 */

#ifndef typedef_cell_wrap_10
#define typedef_cell_wrap_10
typedef struct {
  cell_wrap_0 f1[1];
} cell_wrap_10;
#endif /* typedef_cell_wrap_10 */

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

#ifndef struct_emxArray_int8_T
#define struct_emxArray_int8_T
struct emxArray_int8_T {
  signed char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
};
#endif /* struct_emxArray_int8_T */
#ifndef typedef_emxArray_int8_T
#define typedef_emxArray_int8_T
typedef struct emxArray_int8_T emxArray_int8_T;
#endif /* typedef_emxArray_int8_T */

#ifndef typedef_cell_wrap_14
#define typedef_cell_wrap_14
typedef struct {
  emxArray_real_T *f1;
} cell_wrap_14;
#endif /* typedef_cell_wrap_14 */

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

#ifndef typedef_emxArray_cell_wrap_4
#define typedef_emxArray_cell_wrap_4
typedef struct {
  cell_wrap_4 *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  bool canFreeData;
} emxArray_cell_wrap_4;
#endif /* typedef_emxArray_cell_wrap_4 */

#endif
/* End of code generation (gkmPWM_types.h) */
