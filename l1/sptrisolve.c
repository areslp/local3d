/***************************************************************
 * Sparse Triangular Solve for MATLAB.
 * Copyright (C) Haim Avron, July 2006.
 **************************************************************/

#include "mex.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "matrix.h"

/* Actual mex function */
void mexFunction(int nargout, mxArray *argout[], int nargin, const mxArray *argin[])
{
  const mxArray *matlab_A;
  const mxArray *matlab_B;
  mwSize An, nrhs;
  mwIndex *colptr, *rowind;
  bool lower_flag, notranspose_flag;
  double *values, *B, *X, Aij, Ajj, Aii;
  mwIndex i, j, r, ii, jj, ip, jp;
  char tmp[255];

  /* Make sure A is sparse */
  matlab_A = argin[0];
  if (!mxIsSparse(matlab_A))
    {
      mexPrintf("input matrix A must be sparse\n");
      return;
    }
  An = mxGetN(matlab_A);
  colptr = mxGetJc(matlab_A);
  rowind = mxGetIr(matlab_A);
  values = mxGetPr(matlab_A);

  /* B side vector */
  matlab_B = argin[1];
  nrhs = mxGetN(matlab_B);
  B = mxGetPr(matlab_B);

  /* Lower or upper? Transpose or not? */
  if (nargin < 3)
    lower_flag = true;
  else {
    mxGetString(argin[2], tmp, 255);
    if (strcmp(tmp, "upper") == 0)
      lower_flag = false;
    else
      lower_flag = true;
  }
  if (nargin < 4)
    notranspose_flag = true;
  else {
    mxGetString(argin[3], tmp, 255);
    if (strcmp(tmp, "transpose") == 0 || strcmp(tmp, "transp") == 0)
      notranspose_flag = false;
    else
      notranspose_flag = true;
  }

  /* Allocate space for right hand side */
  argout[0] = mxCreateDoubleMatrix(An, nrhs, mxREAL);
  X = mxGetPr(argout[0]);

  /* Copy B into X, so we modify it only */
  memcpy(X, B, nrhs * An * sizeof(double));

  /* Lower no transpose */
  if (lower_flag &&  notranspose_flag) {
    for(r = 0; r < nrhs; r++)
      for (j = 0; j < An; j++) {

	/* Lower triangular - diagonal element first */
	ip = colptr[j];
	i = rowind[ip];
	assert(i == j);
	Ajj = values[ip];

	X[j + r * An] /=  Ajj;

	for (ip = colptr[j] + 1; ip < colptr[j+1]; ip++) {
	  i = rowind[ip];
	  Aij = values[ip];

	  X[i + r * An] -= X[j + r * An] * Aij;
	}
      }
  }

  /* Lower transpose */
  if (lower_flag && !notranspose_flag)
    {
      for(r = 0; r < nrhs; r++) {
	for (ii = An; ii > 0; ii--) {
	  i = ii - 1;
	  for (jp = colptr[i] + 1; jp < colptr[i + 1]; jp++) {
	    j = rowind[jp];
	    Aij = values[jp];
	    X[i + r * An] -= X[j + r * An] * Aij;
	  }


	  jp = colptr[i];
	  Aii = values[jp];
	  X[i + r * An] /= Aii;
	}
      }
    }

  /* Upper no transpose */
  if (!lower_flag &&  notranspose_flag){
    for(r = 0; r < nrhs; r++)
      for (jj = An; jj > 0; jj--) {
	/* Upper triangular - diagonal element last */
	j = jj  - 1;
	ip = colptr[j + 1] - 1;
	i = rowind[ip];
	assert(i == j);
	Ajj = values[ip];
	X[j + r * An] /=  Ajj;
	for (ip = colptr[j]; ip < colptr[j+1] - 1; ip++) {
	  i = rowind[ip];
	  Aij = values[ip];

	  X[i + r * An] -= X[j + r * An] * Aij;
	}
      }
  }

  /* Upper transpose */
  if (!lower_flag && !notranspose_flag) {
    for(r = 0; r < nrhs; r++)
      for (i = 0; i < An; i++) {
	for (jp = colptr[i]; jp < colptr[i + 1] - 1; jp++) {
	  j = rowind[jp];
	  Aij = values[jp];
	  X[i + r * An] -= X[j + r * An] * Aij;
	}


	jp = colptr[i + 1] - 1;
	Aii = values[jp];
	X[i + r * An] /= Aii;
      }
  }
}
