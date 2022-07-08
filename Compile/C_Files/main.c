/*
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
#include "main.h"
#include "bem_calc.h"
#include "bem_calc_emxutil.h"
#include "bem_calc_initialize.h"
#include "bem_calc_emxAPI.h"
#include "bem_calc_terminate.h"
#include "bem_calc_types.h"
#include "rt_nonfinite.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
static emxArray_char_T *argInit_1xUnbounded_char_T(char *argv);

static char argInit_char_T(void);

static void main_bem_calc(char *argv);

/* Function Definitions */
static emxArray_char_T *argInit_1xUnbounded_char_T(char *path)
{
  emxArray_char_T *result;
  /* Set the size of the array.
Change this size to the value that the application requires. */
  result = emxCreate_char_T(1, 100);
  result->size[1U];
  result->data = path;
  return result;
}

static char argInit_char_T(void)
{
  return '?';
}

static void main_bem_calc(char *path)
{
  cell_0 output_details;
  emxArray_char_T *file;
  emxArray_real_T *output;
  struct0_T BEM;
  emxInitArray_real_T(&output, 2);
  emxInit_struct0_T(&BEM);
  /* Initialize function 'bem_calc' input arguments. */
  /* Initialize function input argument 'file'. */
  file = argInit_1xUnbounded_char_T(path);
  /* Call the entry-point 'bem_calc'. */
  bem_calc(file, output, &output_details, &BEM);
  emxDestroy_struct0_T(BEM);
  emxDestroyArray_real_T(output);
  emxDestroyArray_char_T(file);
}

int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* The initialize function is being called automatically from your entry-point
   * function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_bem_calc(argv[1]);
  /* Terminate the application.
You do not need to do this more than one time. */
  bem_calc_terminate();
  return 0;
}

/* End of code generation (main.c) */
