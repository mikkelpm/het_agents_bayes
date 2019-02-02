/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * exp.c
 *
 * Code generation for function 'exp'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "updateCoefficients_polynomials.h"
#include "exp.h"

/* Function Definitions */
void b_exp(double x[50])
{
  int k;
  for (k = 0; k < 50; k++) {
    x[k] = exp(x[k]);
  }
}

/* End of code generation (exp.c) */
