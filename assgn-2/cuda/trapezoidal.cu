/*
Title- Numerical integration of sin function using Trapezoidal method (parallel execution)
Author- Harshit Agrawal
Description- This integrates the sin function using Trapezoidal method in [0,PI] parallely using CUDA library.
To compile using openmpi run "$ nvcc trapezoidal.c -o trapezoidal.out"
To run the executable run "$ ./montecarlo.out"
Change the value of N, THREADS_PER_BLOCK according to the experimental setup.
*/


#include<stdio.h>
#include<math.h>
#include<time.h>
#define N 1000000
#define THREADS_PER_BLOCK (100)
#define DX M_PI/N


//Device function to compute the area of differential element using trapezoidal approximation
__global__ void area( float *dev_sum) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    dev_sum[index] = ((sinpif(((float)index)/N)+sinpif(((float)(index+1))/N))/2)*DX;
}



int main( void ) {
    float *sum, *dev_sum;
    double size = N * sizeof( float );
    cudaMalloc( (void**)&dev_sum, size );
    sum = (float*)malloc( size );
    clock_t begin = clock();
    area<<< N/THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>(dev_sum);
    
    cudaMemcpy( sum, dev_sum, size, cudaMemcpyDeviceToHost );
    float net_sum = 0;
    for (int i=0; i<N; i++){
        net_sum += sum[i];
    }
    clock_t end = clock();
    double timex = (double)(end-begin)/CLOCKS_PER_SEC;
    printf("Total time of execution: %f \n", timex);
    printf("sum= %f \n", net_sum);
    free( sum );
    cudaFree( dev_sum );
    return 0;
}
