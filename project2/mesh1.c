#include<stdio.h>
#include<math.h>
#include<mpi.h>
#include<stdlib.h>

int main(int argc, char **argv)
{
  int rank, size;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Status stat;
  /*
//Assert for power of 2 processors
  if (size & (size - 1))
  {
    if (rank == 0)
    {
      printf("Please enter a power of 2 for the number of processors.\n");
    }
    MPI_Finalize();
    return 0;
  }
*/
 int N=8; //Number of cells
  int n=N/size;
  double ver[N];
  int offset=rank*n,i;
if(rank==0){
  for (i=offset;i<=offset+n;i++)
	{
	ver[i]=i;
	printf("Process %d responsible for vertices ver[%d]=%f\n",rank,i,ver[i]);
	}
  
  }
if(rank!=0){
  for (int j=offset+1;j<=offset+n;j++)
	{
	ver[j]=j;
	printf("Process %d responsible for vertices ver[%d]=%f\n",rank,j,ver[j]);
	}
  
  }
  MPI_Finalize();
  return 0;
}
