#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>

int main(int argc, char **argv)
{
  int rank, size;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status stat;
if (size &(size-1)){
 	{ if (rank==0){printf("Enter number of power of two processors\n");}
	MPI_Finalize();
	return 0;
	}
  int N=8; //Number of cells
  int n=N/size;
  double ver[N];
  int offset=rank*n,i;
  for (i=offset+1;i<=offset+n;i++)
	{
	printf("Process %d responsible for vertices ver[%d]",rank,i);
	}
MPI_Finalize();
	return 0;
}

