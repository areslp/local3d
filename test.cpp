#include <omp.h>
#include "mex.h"
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
{ 
    int nthreads, tid;
    // omp_set_num_threads(4);
    // printf("Number of processors: %d\n", omp_get_num_procs());
// #pragma omp parallel private(nthreads, tid) 
    {
        tid = omp_get_thread_num();
        // printf("Hello World from thread = %d\n", omp_get_thread_num() );
        int num_out = 0;
        int num_in = 2;
        mxArray *output_array[1], *input_array[2];
        input_array[0] = const_cast<mxArray *>(prhs[0]);
        input_array[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
        double *p = mxGetPr(input_array[1]);
        // *p = omp_get_thread_num()+1;
        *p = 3;
        mexCallMATLAB(num_out, output_array, num_in, input_array, 
            "echorow");
    }
}

