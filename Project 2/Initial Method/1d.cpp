#include <stdio.h>
#include <cmath>
#include <vector>
#include <mpi.h>
#include <iomanip>
#define pi 3.141592653589793238462643383279502
using namespace std;

class Gaussian{
    public:

        double GQ1d(double P[1][2]){
            //compute jacobian
            double v[4][2]={};
           v[0][0] = 0.33981;
           v[0][1] = 0.652145;
           v[1][0] = -0.333981;
           v[1][1] = 0.652145;
           v[2][0] = 0.861136;
           v[2][1] = 0.347855;
           v[3][0] = -0.861136;
           v[3][1] = 0.347855;

            double answer= 0.;
            double z[4] = {};
            for (int i=0 ;i<4;i++){
            z[i] = (P[0][1] - P[0][0])/2.*v[i][0] + (P[0][1] + P[0][1])/2.;
            answer += (P[0][1] - P[0][0])/2.*v[i][1]*(4*pow(pi,2))*cos(pi*z[i])*sin(pi*z[i]);
            }
            return answer;
        }

};


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

    vector<double> localPositions(n+1,0);
    vector<vector<int> > localTopology(n,vector<int>(2,0));

    double meshStep = 1./n;

    int EPP = n/size; // Elements per processor

    if(rank == 0)
        for(int i = 1; i <= (rank+1)*EPP + 1; i++)
            localPositions[i-1] = (i-1)*meshStep;
    else
        for(int i = 1 + rank*EPP; i <= (rank+1)*EPP + 1; i++)
            localPositions[i-1] = (i-1)*meshStep;

    for(int i = rank*EPP + 1; i <= (rank+1)*EPP; i++){

        localTopology[i-1][0] = i;
        localTopology[i-1][1] = i+1;

    }


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

    double A[n+1][n+1] = {};
    double M[n+1][n+1] = {};
    double F[n+1] = {};

    double globalA[n+1][n+1];
    double globalM[n+1][n+1];
    double globalF[n+1];
    double jacobian;
    double f;

    for(int i = rank*EPP; i < (rank+1)*EPP; i++){

        jacobian = abs(localPositions[localTopology[i][1]-1] - localPositions[localTopology[i][0]-1]);

        A[localTopology[i][0]-1][localTopology[i][0]-1] += jacobian * assemA[0][0];
        A[localTopology[i][0]-1][localTopology[i][1]-1] += jacobian * assemA[0][1];
        A[localTopology[i][1]-1][localTopology[i][0]-1] += jacobian * assemA[1][0];
        A[localTopology[i][1]-1][localTopology[i][1]-1] += jacobian * assemA[1][1];

        M[localTopology[i][0]-1][localTopology[i][0]-1] += jacobian * assemM[0][0];
        M[localTopology[i][0]-1][localTopology[i][1]-1] += jacobian * assemM[0][1];
        M[localTopology[i][1]-1][localTopology[i][0]-1] += jacobian * assemM[1][0];
        M[localTopology[i][1]-1][localTopology[i][1]-1] += jacobian * assemM[1][1];
        Gaussian GQ;

        double y[1][2]={};
        y[0][1]=localPositions[localTopology[i][1]-1];
        y[0][0]=localPositions[localTopology[i][0]-1];
        f =GQ.GQ1d(y);
        F[localTopology[i][0]-1] += (2*pow(pi,2)+1)*cos(pi*localPositions[localTopology[i][0]-1]) + f;
        F[localTopology[i][1]-1] += (2*pow(pi,2)+1)*cos(pi*localPositions[localTopology[i][1]-1]) + f;
    }

    MPI_Reduce(&A, &globalA, (n+1)*(n+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&M, &globalM, (n+1)*(n+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&F, &globalF, (n+1), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   if (rank == 0){
	cout <<endl<< "****************** Mass matrix is************"<<endl;
        for(int i = 0; i <= n; i++){

            for(int j = 0; j<= n; j++){

               cout <<globalM[i][j] << "       ";

          }

          cout << endl;

        }
    
       cout <<endl<< "****************** Stiffness matrix is************"<<endl;
        for(int i = 0; i <= n; i++){

            for(int j = 0; j<= n; j++){

               cout <<globalM[i][j] << "       ";

          }

          cout << endl;

        }
    

      cout <<endl<< "****************** Right Hand side ************"<<endl;
     for(int i =0; i <=n; i++){

                cout << globalF[i] << " " << endl;
            }

      }



    MPI_Finalize();

}
