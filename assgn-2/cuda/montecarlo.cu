/*
Title- Numerical integration of sin function using Monte Carlo method (parallel execution)
Author- Harshit Agrawal
Description- This integrates the sin function using Monte Carlo method in [0,PI] parallely using CUDA library.
To compile using CUDA run "$ nvcc montecarlo.cu -o montecarlo.out"
To run the executable run "$ ./montecarlo.out"
Change the value of N, THREADS_PER_BLOCK according to the experimental setup.
*/
#include<stdio.h>
#include<math.h>
#include<curand.h>
#include<curand_kernel.h>
#include<time.h>
#define N 1000000
#define THREADS_PER_BLOCK (100)


//This function generates the random numbers
__device__ float generate( curandState* globalState, int ind )
{
    curandState localState = globalState[ind];
    float RANDOM = curand_uniform( &localState );
    globalState[ind] = localState;
    return RANDOM;
}

//This function sets up the random number generator
__global__ void setup_kernel ( curandState * state, unsigned long seed )
{
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init ( seed, id, 0, &state[id] );        
}

//This function executes the main monte carlo method
__global__ void area( int *dev_sum, curandState* globalState) {
    __shared__ int success[THREADS_PER_BLOCK];
    float r1 = generate(globalState, threadIdx.x + 2 * blockIdx.x * blockDim.x)*M_PI, r2 = generate(globalState, threadIdx.x + 2*blockIdx.x * blockDim.x+blockDim.x);
    float sinr1 = sinf(r1);
    if (sinr1>=r2){
        success[threadIdx.x] = 1;
    } else{
        success[threadIdx.x] = 0;
    }
    __syncthreads();
    int suc_num = 0;
    if (0==threadIdx.x){
        for (int i=0; i<THREADS_PER_BLOCK; i++){
            if (success[i]){
                ++suc_num;
            }
        }
        dev_sum[blockIdx.x] = suc_num;
    }
}    



int main( void ) {
    int *sum, *dev_sum;
    float integral = 0;
    curandState* devStates;
    cudaMalloc ( &devStates, N*sizeof( curandState ) );
    double size = (N/THREADS_PER_BLOCK)*sizeof(int);
    clock_t begin = clock();

    setup_kernel <<< 2*N/THREADS_PER_BLOCK, THREADS_PER_BLOCK>>> ( devStates,unsigned(time(NULL)) );

    cudaMalloc( (void**)&dev_sum, size );
    sum = (int*)malloc( size );

    area<<< N/THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>(dev_sum, devStates);

    cudaMemcpy( sum, dev_sum, size, cudaMemcpyDeviceToHost );
    float net_sum = 0;
    for (int i=0; i<(N/THREADS_PER_BLOCK); i++){
        net_sum += sum[i];
    }

    integral = ((float)net_sum)*M_PI/N;
    clock_t end = clock();
    double timex = (double)(end-begin)/CLOCKS_PER_SEC;
    printf("Total time of execution: %f \n", timex);
    printf("integral= %f \n", integral);
    free( sum );
    cudaFree( dev_sum ), cudaFree(devStates);
    return 0;
}
