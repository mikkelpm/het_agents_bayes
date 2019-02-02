/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "updateCoefficients_polynomials.h"
#include "main.h"
#include "updateCoefficients_polynomials_terminate.h"
#include "updateCoefficients_polynomials_emxAPI.h"
#include "updateCoefficients_polynomials_initialize.h"

/* Function Declarations */
static void argInit_25x1_real_T(double result[25]);
static void argInit_25x25_real_T(double result[625]);
static void argInit_2x25_real_T(double result[50]);
static void argInit_2x2_real_T(double result[4]);
static void argInit_2x50_real_T(double result[100]);
static double argInit_real_T(void);
static void main_updateCoefficients_polynomials(void);

/* Function Definitions */
static void argInit_25x1_real_T(double result[25])
{
  int idx0;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 25; idx0++) {
    /* Set the value of the array element.
       Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

static void argInit_25x25_real_T(double result[625])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 25; idx0++) {
    for (idx1 = 0; idx1 < 25; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + 25 * idx1] = argInit_real_T();
    }
  }
}

static void argInit_2x25_real_T(double result[50])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    for (idx1 = 0; idx1 < 25; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + (idx1 << 1)] = argInit_real_T();
    }
  }
}

static void argInit_2x2_real_T(double result[4])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    for (idx1 = 0; idx1 < 2; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + (idx1 << 1)] = argInit_real_T();
    }
  }
}

static void argInit_2x50_real_T(double result[100])
{
  int idx0;
  int idx1;

  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 2; idx0++) {
    for (idx1 = 0; idx1 < 50; idx1++) {
      /* Set the value of the array element.
         Change this value to the value that the application requires. */
      result[idx0 + (idx1 << 1)] = argInit_real_T();
    }
  }
}

static double argInit_real_T(void)
{
  return 0.0;
}

static void main_updateCoefficients_polynomials(void)
{
  emxArray_real_T *mCoefficientsNew;
  double mCoefficients[50];
  double bbeta;
  double ssigma;
  double aaBar;
  double mmu;
  double ttau;
  double mEpsilonTransition[4];
  double dv0[50];
  double dv1[50];
  double dv2[625];
  double dv3[25];
  double dv4[100];
  emxInitArray_real_T(&mCoefficientsNew, 2);

  /* Initialize function 'updateCoefficients_polynomials' input arguments. */
  /* Initialize function input argument 'mCoefficients'. */
  argInit_2x25_real_T(mCoefficients);
  bbeta = argInit_real_T();
  ssigma = argInit_real_T();
  aaBar = argInit_real_T();
  mmu = argInit_real_T();
  ttau = argInit_real_T();

  /* Initialize function input argument 'mEpsilonTransition'. */
  argInit_2x2_real_T(mEpsilonTransition);

  /* Initialize function input argument 'mEpsilonGrid'. */
  /* Initialize function input argument 'mAssetsGrid'. */
  /* Initialize function input argument 'vAssetsPoly'. */
  /* Initialize function input argument 'vAssetsPolySquared'. */
  /* Initialize function input argument 'mEpsilonPrimeGrid'. */
  /* Call the entry-point 'updateCoefficients_polynomials'. */
  argInit_2x25_real_T(dv0);
  argInit_2x25_real_T(dv1);
  argInit_25x25_real_T(dv2);
  argInit_25x1_real_T(dv3);
  argInit_2x50_real_T(dv4);
  updateCoefficients_polynomials(mCoefficients, bbeta, ssigma, aaBar, mmu, ttau,
    mEpsilonTransition, argInit_real_T(), argInit_real_T(), argInit_real_T(),
    argInit_real_T(), argInit_real_T(), dv0, dv1, dv2, dv3, dv4, argInit_real_T(),
    argInit_real_T(), mCoefficientsNew);
  emxDestroyArray_real_T(mCoefficientsNew);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  /* Initialize the application.
     You do not need to do this more than one time. */
  updateCoefficients_polynomials_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  main_updateCoefficients_polynomials();

  /* Terminate the application.
     You do not need to do this more than one time. */
  updateCoefficients_polynomials_terminate();
  return 0;
}

/* End of code generation (main.c) */
