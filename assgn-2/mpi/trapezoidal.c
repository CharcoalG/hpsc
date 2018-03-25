/*
Title- Numerical integration of sin function using Trapezoidal method (parallel execution)
Author- Harshit Agrawal
Description- This integrates the sin function using Trapezoidal method in [0,PI] parallely using openmpi library.
To compile using openmpi run "$ mpicc trapezoidal.c -o trapezoidal.out -lm"
To run the executable run (with 4 mpi threads here) "$ mpirun -np 4 trapezoidal.out"
Change the value of DIVS according to the experimental setup.
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"mpi.h"
#include<time.h>
#define DIVS 10000000 //change this as per need


int main(int argc, char** argv){
    double dx = M_PI/DIVS;
    int my_PE_num, num_procs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_PE_num);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (my_PE_num == 0){ //Master Thread
        clock_t begin = clock();
        double recv_sum;
        MPI_Status status;
        double sum2 = 0;
        for (int i=0; i<num_procs-1; i++){
            MPI_Recv(&recv_sum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status); //recieving slave thread sums
            sum2 = sum2 + recv_sum; //adding slave sums to a common sum
        }
        clock_t end = clock();
        double timex = (double)(end-begin)/CLOCKS_PER_SEC;
        printf("Total time of execution: %f \n", timex);
        printf("Sum is %f \n", sum2);
    } else{ //Slave Threads
        int count = DIVS/(num_procs-1), remainder = DIVS % (num_procs-1);
        int start, stop;
        double sum1 = 0, prev_value = 0;
        //division of work load
        if (my_PE_num<remainder+1){
            start = (my_PE_num-1)*(count+1);
            stop = start+count;
        } else{
            start = (my_PE_num-1)*count+remainder;
            stop = start+(count-1);
        }

        //calculation of sum
        prev_value = sin((start-1)*dx);
        for (int i = start; i<stop; i++){
            double curr_value = sin(i*dx);
            sum1 = sum1 + (curr_value+prev_value)*dx/2;
            prev_value = curr_value;
        }
        MPI_Send(&sum1, 1, MPI_DOUBLE, 0, my_PE_num, MPI_COMM_WORLD); //sending to master thread

    }

    MPI_Finalize();
}