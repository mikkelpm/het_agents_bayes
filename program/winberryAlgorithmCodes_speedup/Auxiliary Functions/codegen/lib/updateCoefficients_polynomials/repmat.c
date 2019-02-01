/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "updateCoefficients_polynomials.h"
#include "repmat.h"
#include "updateCoefficients_polynomials_emxutil.h"

/* Function Definitions */
void repmat(const double a_data[], const int a_size[2], const double varargin_1
            [2], emxArray_real_T *b)
{
  int ntilerows;
  boolean_T p;
  int jcol;
  int ibmat;
  int itilerow;
  ntilerows = b->size[0] * b->size[1];
  b->size[0] = (int)varargin_1[0];
  b->size[1] = a_size[1];
  emxEnsureCapacity_real_T(b, ntilerows);
  if (!(a_size[1] == 0)) {
    if ((int)varargin_1[0] == 0) {
      p = true;
    } else if (a_size[1] == 0) {
      p = true;
    } else {
      p = false;
    }

    if (!p) {
      ntilerows = (int)varargin_1[0];
      for (jcol = 0; jcol + 1 <= a_size[1]; jcol++) {
        ibmat = jcol * ntilerows;
        for (itilerow = 1; itilerow <= ntilerows; itilerow++) {
          b->data[(ibmat + itilerow) - 1] = a_data[jcol];
        }
      }
    }
  }
}

/* End of code generation (repmat.c) */
