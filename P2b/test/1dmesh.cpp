#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <mpi.h>

using namespace std;


int main(int argc, char** argv) {

    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n;

    if(rank == 0){
        n=8;
    }

    MPI_Bcast(&n,1,MPI_INT, 0, MPI_COMM_WORLD);

    int EPP = n/size; // Elements per processor
    double v_local[EPP][2];
    int v[n][2];
    double v_global[n]={};
    double meshStep = 1./n;
    int globalv[n][2]={};

    for(int i=rank*EPP;i<rank*EPP+EPP;i++){
      v[i][0]=i;
      v[i][1]=i+1; 
     fflush(stdout); }   
    
    MPI_Gatherv(&v,n*2, MPI_INT,
                 &globalv,n*2,MPI_INT, MPI_COMM_WORLD);
/*    
   if(rank==0){
    for(int i=rank*EPP;i<rank*EPP+EPP;i++){
        for(int j=0;j<2;j++){cout<<"v["<<i<<"]["<<j<<"]="<<v[i][j]<<" ";}
					cout<<endl;   }
     }	
*/
   if(rank==0){
	for(int i=0;i<n;i++){
        for(int j=0;j<2;j++){cout<<"v["<<i<<"]["<<j<<"]="<<globalv[i][j]<<" ";}
					cout<<endl;   }
     }			

    
    MPI_Finalize();

}
