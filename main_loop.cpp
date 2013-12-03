#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <omp.h>
#include <iostream>

using namespace std;

// for iteration, vertex, normals, index, lambda
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// declare variables
	double *vertex, *normals, *index;
	double *normals_new;
    double lambda;
	int num;

	//associate inputs
	vertex = mxGetPr(prhs[0]);
	normals = mxGetPr(prhs[1]);
	index = mxGetPr(prhs[2]);
	lambda = mxGetScalar(prhs[3]);

	// figure out dimensions
	num = mxGetN(prhs[0]);

	// output
	plhs[0] = mxCreateDoubleMatrix(num, 3, mxREAL);
	normals_new = mxGetPr(plhs[0]);

    omp_set_num_threads(4);

	#pragma omp parallel for
	for (int i=0;i<num;i++)
	{
        // int nThreads = omp_get_num_threads();
        // printf("nThreads is %d\n",nThreads);

        int num_out, num_in;

        // step1: [X,mapping,idx]=genrealdata_batch(i,index,vertex,normals); 
        num_out = 3;
        num_in = 4;
        mxArray *output_array1[3], *input_array1[4];
        input_array1[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        double *p = mxGetPr(input_array1[0]);
        *p = i+1;
        input_array1[1] = const_cast<mxArray*>(prhs[2]);
        input_array1[2] = const_cast<mxArray*>(prhs[0]);
        input_array1[3] = const_cast<mxArray*>(prhs[1]);
        mexCallMATLAB(num_out, output_array1, num_in, input_array1, 
            "genrealdata_batch");

        // // step2: [Z,E]=ladmp_lrr_fast(X,lambda,rho,DEBUG);
        // num_out = 2;
        // num_in = 2;
        // mxArray *output_array2[3], *input_array2[2];
        // input_array2[0] = output_array1[0];
        // input_array2[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        // p = mxGetPr(input_array2[1]);
        // *p = lambda;
        // mexCallMATLAB(num_out, output_array2, num_in, input_array2, 
            // "ladmp_lrr_fast");

        // // step3: idxs=cut(X,Z,idx,mapping,false);
        // num_out = 1;
        // num_in = 5;
        // mxArray *output_array3[1], *input_array3[5];
        // input_array3[0] = output_array1[0];
        // input_array3[1] = output_array2[0];
        // input_array3[2] = output_array1[2];
        // input_array3[3] = output_array1[1];
        // input_array3[4] = mxCreateDoubleMatrix(1, 1, mxREAL);
        // p = mxGetPr(input_array3[4]);
        // *p = 0;
        // mexCallMATLAB(num_out, output_array3, num_in, input_array3, 
            // "cut");

        // // step4: normals_new(i,:)=lsqnormest2(X,idxs);
        // num_out = 1;
        // num_in = 2;
        // mxArray *output_array4[1], *input_array4[2];
        // input_array4[0] = const_cast<mxArray*>(prhs[0]);
        // input_array4[1] = output_array3[0];
        // mexCallMATLAB(num_out, output_array4, num_in, input_array4, 
            // "lsqnormest2");

        // // ith normal
        // double* ni;
        // ni = mxGetPr(output_array4[0]);
        // // ith row, 1 column
        // int start=1+i*3;
        // normals_new[start] = ni[0];
        // normals_new[start+1] = ni[1];
        // normals_new[start+2] = ni[2];

        // mxDestroyArray(input_array2[1]);
        // mxDestroyArray(input_array3[4]);

        // mxDestroyArray(output_array1[0]);
        // mxDestroyArray(output_array1[1]);
        // mxDestroyArray(output_array1[2]);
        // mxDestroyArray(output_array2[0]);
        // mxDestroyArray(output_array2[1]);
        // mxDestroyArray(output_array3[0]);
        // mxDestroyArray(output_array4[0]);
	}

}
