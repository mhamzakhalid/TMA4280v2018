#include <mpi.h>
#include <stdio.h>
#include <math.h>

const double pi = 3.14159265358979323846;

int main(int argc, char **argv)
{
    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status stat;

    if( size & (size - 1) ){

        if(rank == 0){
            printf("Please enter a power of 2 for the number of processors.\n");
        }
        MPI_Finalize();
        return 0;

    }

    int n = 999999;
    double local_sum = 0;
    double time_start;

    if(rank == 0){

        time_start = MPI_Wtime();

    }

    MPI_Bcast(&n,1,MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);


    for(int i = rank+1; i<=n; i += size){

            local_sum += 1/(pow(i,2));

    }

    double sigma = local_sum;
    double sigmaq;

    for(int d = 0; d <=log2(size) - 1; d++){
        int q = rank^((int)pow(2,d));
        MPI_Send(&sigma, 1, MPI_DOUBLE, q, 1, MPI_COMM_WORLD);
        MPI_Recv(&sigmaq, 1, MPI_DOUBLE, q, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        sigma += sigmaq;
    }


    double global_sum = sqrt(6*sigma);
    // Uncomment line below to see result on each process
    // printf("Process %d has the solution as: %.17g\n", rank, global_sum);

    MPI_Barrier(MPI_COMM_WORLD);

    if(rank == 0){

        double duration = MPI_Wtime() - time_start;
        printf("Recursive doubling calculated the sum on each process in: %.17g seconds.", duration);

    }

    MPI_Finalize();

    return 0;

}
