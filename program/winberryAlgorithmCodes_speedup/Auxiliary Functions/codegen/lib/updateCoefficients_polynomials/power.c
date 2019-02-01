/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * power.c
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "updateCoefficients_polynomials.h"
#include "power.h"
#include "updateCoefficients_polynomials_rtwutil.h"

/* Function Definitions */
void power(const double a[50], double b, double y[50])
{
  int k;
  for (k = 0; k < 50; k++) {
    y[k] = rt_powd_snf(a[k], b);
  }
}

/* End of code generation (power.c) */
