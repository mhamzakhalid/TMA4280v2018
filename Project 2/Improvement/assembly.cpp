#include <stdio.h>
#include <cmath>
#include <vector>
#include <mpi.h>
#include <iomanip>

using namespace std;

int main(int argc, char** argv) {

    int rank, size;

    MPI_Request request, request2,request3,request4;
    MPI_Status status;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n, i, j;

    if(rank == 0){

       n=8;

    }

    MPI_Bcast(&n,1,MPI_INT, 0, MPI_COMM_WORLD);

    int localn = n/size;

    int localTopology[localn][2] = {};

    double localGeometry[localn+1] = {};

    int ownership[localn+1][2] = {};

    int range[2] = {};

    int globalMap[localn+1] = {};

    for (i = 0; i < localn; i++){

        localTopology[i][0] = i;
        localTopology[i][1] = i+1;

        localGeometry[i] = (1./size)*rank + i*1./n;
        globalMap[i] = i + localn*rank;

        if(rank != 0 && i == 0){
            ownership[i][0] = 0;
            ownership[i][1] = rank - 1;
        }else if(rank != (size-1) && i==localn){

            ownership[i][0] = 2;
            ownership[i][1] = rank+1;

        }else{

            ownership[i][0] = 1;
            ownership[i][1] = rank;

        }

  

    localGeometry[i] = (1./size)*rank + i*1./n;
    globalMap[i] = i + localn*rank;
    ownership[i][0] = 1;
  }

    range[0] = globalMap[0];
    range[1] = localn+1;

    double assemA[2][2] = {};
    assemA[0][0] = 1./3.;
    assemA[0][1] = 1./6.;
    assemA[1][0] = 1./6.;
    assemA[1][1] = 1./3.;

    double assemM[2][2] = {};
    assemM[0][0] = 1.;
    assemM[0][1] = -1.;
    assemM[1][0] = -1.;
    assemM[1][1] = 1.;

    double localA[localn+(rank==0)][n+1]={};
    double localM[localn+(rank==0)][n+1]={};
    double tempSendA, tempRecvA, tempSendM, tempRecvM;
    double jacobian;
    double globalv[n+1]={};
    double globalA[n+1][n+1]={};
    double globalM[n+1][n+1]={};
    for (i = 0; i <localn; i++){

        jacobian = abs(localGeometry[localTopology[i][1]] - localGeometry[localTopology[i][0]]);


            for (j = 0; j <localn; j++){

                if(ownership[i][0] == 1){

                    localA[i - (rank!=0)][globalMap[j]] += jacobian*assemA[i][j];
                    localM[i - (rank!=0)][globalMap[j]] += jacobian*assemM[i][j];

                }else if(ownership[i][0] == 0){

                    tempSendA = jacobian*assemA[i][j];
                    MPI_Isend(&tempSendA, 1, MPI_DOUBLE, ownership[i][1], 0, MPI_COMM_WORLD,&request);
                    tempSendM = jacobian*assemM[i][j];
                    MPI_Isend(&tempSendM, 1, MPI_DOUBLE, ownership[i][1], 1, MPI_COMM_WORLD,&request2);
			MPI_Wait(&request, &status);
                       MPI_Wait(&request2, &status);

                }else if(ownership[i][0] == 2){

                    MPI_Irecv(&tempRecvA, 1, MPI_DOUBLE, ownership[i][1], 0, MPI_COMM_WORLD,&request3);
                    localA[i - (rank!=0)][globalMap[j]] += jacobian*assemA[i][j] + tempRecvA;
                    MPI_Irecv(&tempRecvM, 1, MPI_DOUBLE, ownership[i][1], 1, MPI_COMM_WORLD,&request4);
                    localM[i - (rank!=0)][globalMap[j]] += jacobian*assemM[i][j] + tempRecvM;
			 MPI_Wait(&request3, &status);
    			MPI_Wait(&request4, &status);
                
            }
        }
    }
   MPI_Reduce(&localA, &globalA, (n+1)*(n+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&localM, &globalM, (n+1)*(n+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//************************************For Printing Local Topology uncomment this*************************************//
/*   if(rank==0){
		for(int i=0;i<localn;i++)
 	cout <<" localTopology["<<i<<"][0] ="<< localTopology[i][0] <<"  "<< "localTopology["<<i<<"][1] ="<< localTopology[i][1]<<endl; 
} 
*/
if (rank == 0){
        cout<<endl<<"***************************** Mass matrix *******************"<<endl;
        for(int i = 0; i <n; i++){

            for(int j = 0; j<n; j++){

               cout <<globalA[i][j] << "     ";

          }

          cout << endl;

        }
    cout<<endl<<"***************************** Stiffness matrix *******************"<<endl;
        for(int i = 0; i <n; i++){

            for(int j = 0; j<n; j++){

               cout <<globalM[i][j] << "     ";

          }

          cout << endl;

        }
    }








































    MPI_Finalize();

}
