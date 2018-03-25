/*
Title- Numerical integration of sin function using Monte Carlo method (parallel execution)
Author- Harshit Agrawal
Description- This integrates the sin function using Monte Carlo method in [0,PI] parallely using openmpi library.
To compile using openmpi run "$ mpicc montecarlo.c -o montecarlo.out -lm"
To run the executable run (with 4 mpi threads here) "$ mpirun -np 4 montecarlo.out"
Change the value of DIVS according to the experimental setup.
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"mpi.h"
#include<time.h>
#define TRIALS 10000000 //change this as per need


int main(int argc, char** argv){
    int my_PE_num, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (my_PE_num == 0){ //Master Thread
        clock_t begin = clock();
        unsigned long recv_sum;
        MPI_Status status;
        unsigned long sum2 = 0;
        for (int i=0; i<num_procs-1; i++){
            MPI_Recv(&recv_sum, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); //recieving slave thread sums
            sum2 = sum2 + recv_sum; //adding slave sums to a common sum
        }
        clock_t end = clock();
        double timex = (double)(end-begin)/CLOCKS_PER_SEC;
        printf("Total time of execution: %f \n", timex);
        printf("Sum is %f \n", (double)sum2*M_PI/TRIALS);
    } else{ //Slave Threads
        int count, remainder = TRIALS % (num_procs-1);
        unsigned long success = 0, denom = RAND_MAX;
        double r1=0, r2=0, sinr1=0;
        srand(time(NULL));
        //division of work load
        if (my_PE_num<remainder+1){
            count = TRIALS/(num_procs-1);
        } else{
            count = TRIALS/(num_procs-1)-1;
        }

        for (int i=0; i<count; i++){
            r1 = rand()*M_PI/denom;
            r2 = rand()*1.0/denom;
            sinr1 = sin(r1);
            if(sinr1>=r2){
                success+=1;
            }
        }
        MPI_Send(&success, 1, MPI_LONG, 0, my_PE_num, MPI_COMM_WORLD); //sending to master thread

    }

    MPI_Finalize();
}