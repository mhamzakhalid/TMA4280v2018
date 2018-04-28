#include <stdio.h>
#include <cmath>
#include <vector>
#include <mpi.h>
#include <iomanip>
#include <iostream>
using namespace std;
#define pi 3.141592653589793238462643383279502

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
struct addScalar{

    double val;
    addScalar(double val) : val(val){};

    void operator()(double &elem) const{

        elem += val;

    }

};
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

class Gaussian{
    public:

        double GQ2d(double P[3][2]){

            double J[2][2];
            J[0][0] = P[0][0] - P[2][0];
            J[1][0] = P[0][1] - P[2][1];
            J[0][1] = P[1][0] - P[2][0];
            J[1][1] = P[1][1] - P[2][1];

            double det;
            det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

            double K[2][2];
            K[0][0] = + J[1][1] / det;
            K[1][0] = - J[1][0] / det;
            K[0][1] = - J[0][1] / det;
            K[1][1] = + J[0][0] / det;

            double lambda[3][3];
            double rho=1./3.;
            //Making objects from Classes


            double K1;
            K1=abs(det/2.);
            double X[3][2];

            X[0][0]=lambda[0][0]*P[0][0] + lambda[0][1]*P[1][0] + lambda[0][2]*P[2][0];
            X[0][1]=lambda[0][0]*P[0][1] + lambda[0][1]*P[1][1] + lambda[0][2]*P[2][1];
            X[1][0]=lambda[1][0]*P[0][0] + lambda[1][1]*P[1][0] + lambda[1][2]*P[2][0];
            X[1][1]=lambda[1][0]*P[0][1] + lambda[1][1]*P[1][1] + lambda[1][2]*P[2][1];
            X[2][0]=lambda[2][0]*P[0][0] + lambda[2][1]*P[1][0] + lambda[2][2]*P[2][0];
            X[2][1]=lambda[2][0]*P[0][1] + lambda[2][1]*P[1][1] + lambda[2][2]*P[2][1];
            double solution;
            double solution1;
            double solution2;
            double solution3;

            solution1=K1*rho*2*(pow(pi,2)+1)*cos(pi*X[0][0])*sin(pi*X[0][1]);
            solution2=K1*rho*2*(pow(pi,2)+1)*cos(pi*X[1][0])*sin(pi*X[1][1]);
            solution3=K1*rho*2*(pow(pi,2)+1)*cos(pi*X[2][0])*sin(pi*X[2][1]);
            solution=solution1+solution2+solution3;
            return solution;
        }

};
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN FUNCTION STARTS HERE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%///


int main(int argc, char** argv) {

    int rank, size;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int n;

    if(rank == 0){

        printf("Please enter the number of square cells in each direction: ");
        fflush(stdout);
        cin >> n;

    }

    MPI_Bcast(&n,1,MPI_INT, 0, MPI_COMM_WORLD);

    int SMS = (int)n/sqrt(size); // SUB MESH SIZE
    int NOSC = pow(n,2); // NUMBER OF SQUARE CELLS
    double meshStep = 1./n;

    vector<vector<double> > localPositions(pow(n+1,2), vector<double>(2,0));
    vector<vector<int> > localSquareCells(NOSC, vector<int>(4,0));
    vector<vector<int> > localTriangleCells(NOSC*2, vector<int>(3,0));

    int squareNum;
    int pointNum;

    int PIED = sqrt(size); // PROCESSORS IN EACH DIRECTION

    int PSSV = SMS*n*floor(rank/PIED) + (rank % PIED)*SMS; // PROCESSOR SQUARE SHIFT VALUE
    int PPSV = SMS*(n+1)*floor(rank/PIED) + (rank % PIED)*SMS; // PROCESSOR POINT SHIFT VALUE

    for (int i = 0; i < pow(SMS,2); i++) {

        squareNum = n*floor(i/SMS) + (i % SMS) + 1 + PSSV;
        pointNum = (n+1)*floor(i/SMS) + (i % SMS) + 1 + PPSV;

        localSquareCells[squareNum-1][0] = pointNum;
        localSquareCells[squareNum-1][1] = pointNum + 1;
        localSquareCells[squareNum-1][2] = pointNum + n + 1;
        localSquareCells[squareNum-1][3] = pointNum + n + 2;

        localTriangleCells[squareNum-1][0] = pointNum;
        localTriangleCells[squareNum-1][1] = pointNum + 1;
        localTriangleCells[squareNum-1][2] = pointNum + n + 1;

        localTriangleCells[squareNum+NOSC-1][0] = pointNum + n + 2;
        localTriangleCells[squareNum+NOSC-1][1] = pointNum + n + 1;
        localTriangleCells[squareNum+NOSC-1][2] = pointNum + 1;

    }

    for (int i = 0; i < pow(SMS+1,2); i++){

        pointNum = (n+1)*floor(i/(SMS+1)) + (i % (SMS+1)) + 1 + PPSV;

        localPositions[pointNum-1][0] = ((pointNum-1) % (n+1))*meshStep;
        localPositions[pointNum-1][1] = floor((pointNum-1)/(n+1))*meshStep;


    }

    double assemA[3][3] = {};
    double assemM[3][3] = {};

    assemA[0][0] = 1./12.;
    assemA[0][1] = 1./24.;
    assemA[0][2] = 1./24.;
    assemA[1][0] = 1./24.;
    assemA[1][1] = 1./12.;
    assemA[1][2] = 1./24.;
    assemA[2][0] = 1./24.;
    assemA[2][1] = 1./24.;
    assemA[2][2] = 1./12.;

    assemM[0][0] = 1.;
    assemM[0][1] = -1./2.;
    assemM[0][2] = -1./2.;
    assemM[1][0] = -1./2.;
    assemM[1][1] = 1./2.;
    assemM[1][2] = 0.;
    assemM[2][0] = -1./2.;
    assemM[2][1] = 0.;
    assemM[2][2] = 1./2.;

    double A[(int)pow(n+1,2)][(int)pow(n+1,2)] = {};
    double M[(int)pow(n+1,2)][(int)pow(n+1,2)] = {};
    double F[(int)pow(n+1,2)] = {};

    double globalA[(int)pow(n+1,2)][(int)pow(n+1,2)];
    double globalM[(int)pow(n+1,2)][(int)pow(n+1,2)];

    double jacobianLT;
    double jacobianUT;
    double plt[3][2]={};
    double put[3][2]={};
    Gaussian GQ;
    double f_lt;
    double f_ut;

    for(int i = 0; i < pow(SMS,2); i++){

        squareNum = n*floor(i/SMS) + (i % SMS) + 1 + PSSV;



        plt[0][0]=localPositions[localTriangleCells[squareNum-1][0]-1][0];
        plt[0][1]=localPositions[localTriangleCells[squareNum-1][0]-1][1];
        plt[1][0]=localPositions[localTriangleCells[squareNum-1][1]-1][0];
        plt[1][1]=localPositions[localTriangleCells[squareNum-1][1]-1][1];
        plt[2][0]=localPositions[localTriangleCells[squareNum-1][2]-1][0];
        plt[2][1]=localPositions[localTriangleCells[squareNum-1][2]-1][1];




        put[0][0]=localPositions[localTriangleCells[squareNum+NOSC-1][0]-1][0];
        put[0][1]=localPositions[localTriangleCells[squareNum+NOSC-1][0]-1][1];
        put[1][0]=localPositions[localTriangleCells[squareNum+NOSC-1][1]-1][0];
        put[1][1]=localPositions[localTriangleCells[squareNum+NOSC-1][1]-1][1];
        put[2][0]=localPositions[localTriangleCells[squareNum+NOSC-1][2]-1][0];
        put[2][1]=localPositions[localTriangleCells[squareNum+NOSC-1][2]-1][1];

          jacobianLT = abs((plt[0][0] - plt[2][0]) * (plt[1][1] - plt[2][1]) - (plt[1][0]- plt[2][0] ) * (plt[0][1] - plt[2][1]));
          jacobianUT = abs((put[0][0] - put[2][0]) * (put[1][1] - put[2][1]) - (put[1][0]- put[2][0] ) * (put[0][1] - put[2][1]));

//        jacobianLT = abs((localPositions[localTriangleCells[squareNum-1][0]-1][0] - localPositions[localTriangleCells[squareNum-1][2]-1][0]) * (localPositions[localTriangleCells[squareNum-1][1]-1][1] - localPositions[localTriangleCells[squareNum-1][2]-1][1]) - (localPositions[localTriangleCells[squareNum-1][1]-1][0] - localPositions[localTriangleCells[squareNum-1][2]-1][0]) * (localPositions[localTriangleCells[squareNum-1][0]-1][1] - localPositions[localTriangleCells[squareNum-1][2]-1][1]));
//        jacobianUT = abs((localPositions[localTriangleCells[squareNum+NOSC-1][0]-1][0] - localPositions[localTriangleCells[squareNum+NOSC-1][2]-1][0]) * (localPositions[localTriangleCells[squareNum+NOSC-1][1]-1][1] - localPositions[localTriangleCells[squareNum+NOSC-1][2]-1][1]) - (localPositions[localTriangleCells[squareNum+NOSC-1][1]-1][0] - localPositions[localTriangleCells[squareNum+NOSC-1][2]-1][0]) * (localPositions[localTriangleCells[squareNum+NOSC-1][0]-1][1] - localPositions[localTriangleCells[squareNum+NOSC-1][2]-1][1]));

        A[localTriangleCells[squareNum-1][0]-1][localTriangleCells[squareNum-1][0]-1] += jacobianLT * assemA[0][0];
        A[localTriangleCells[squareNum-1][0]-1][localTriangleCells[squareNum-1][1]-1] += jacobianLT * assemA[0][1];
        A[localTriangleCells[squareNum-1][0]-1][localTriangleCells[squareNum-1][2]-1] += jacobianLT * assemA[0][2];
        A[localTriangleCells[squareNum-1][1]-1][localTriangleCells[squareNum-1][0]-1] += jacobianLT * assemA[1][0];
        A[localTriangleCells[squareNum-1][1]-1][localTriangleCells[squareNum-1][1]-1] += jacobianLT * assemA[1][1];
        A[localTriangleCells[squareNum-1][1]-1][localTriangleCells[squareNum-1][2]-1] += jacobianLT * assemA[1][2];
        A[localTriangleCells[squareNum-1][2]-1][localTriangleCells[squareNum-1][0]-1] += jacobianLT * assemA[2][0];
        A[localTriangleCells[squareNum-1][2]-1][localTriangleCells[squareNum-1][1]-1] += jacobianLT * assemA[2][1];
        A[localTriangleCells[squareNum-1][2]-1][localTriangleCells[squareNum-1][2]-1] += jacobianLT * assemA[2][2];

        M[localTriangleCells[squareNum-1][0]-1][localTriangleCells[squareNum-1][0]-1] += jacobianLT * assemM[0][0];
        M[localTriangleCells[squareNum-1][0]-1][localTriangleCells[squareNum-1][1]-1] += jacobianLT * assemM[0][1];
        M[localTriangleCells[squareNum-1][0]-1][localTriangleCells[squareNum-1][2]-1] += jacobianLT * assemM[0][2];
        M[localTriangleCells[squareNum-1][1]-1][localTriangleCells[squareNum-1][0]-1] += jacobianLT * assemM[1][0];
        M[localTriangleCells[squareNum-1][1]-1][localTriangleCells[squareNum-1][1]-1] += jacobianLT * assemM[1][1];
        M[localTriangleCells[squareNum-1][1]-1][localTriangleCells[squareNum-1][2]-1] += jacobianLT * assemM[1][2];
        M[localTriangleCells[squareNum-1][2]-1][localTriangleCells[squareNum-1][0]-1] += jacobianLT * assemM[2][0];
        M[localTriangleCells[squareNum-1][2]-1][localTriangleCells[squareNum-1][1]-1] += jacobianLT * assemM[2][1];
        M[localTriangleCells[squareNum-1][2]-1][localTriangleCells[squareNum-1][2]-1] += jacobianLT * assemM[2][2];

        A[localTriangleCells[squareNum+NOSC-1][0]-1][localTriangleCells[squareNum+NOSC-1][0]-1] += jacobianUT * assemA[0][0];
        A[localTriangleCells[squareNum+NOSC-1][0]-1][localTriangleCells[squareNum+NOSC-1][1]-1] += jacobianUT * assemA[0][1];
        A[localTriangleCells[squareNum+NOSC-1][0]-1][localTriangleCells[squareNum+NOSC-1][2]-1] += jacobianUT * assemA[0][2];
        A[localTriangleCells[squareNum+NOSC-1][1]-1][localTriangleCells[squareNum+NOSC-1][0]-1] += jacobianUT * assemA[1][0];
        A[localTriangleCells[squareNum+NOSC-1][1]-1][localTriangleCells[squareNum+NOSC-1][1]-1] += jacobianUT * assemA[1][1];
        A[localTriangleCells[squareNum+NOSC-1][1]-1][localTriangleCells[squareNum+NOSC-1][2]-1] += jacobianUT * assemA[1][2];
        A[localTriangleCells[squareNum+NOSC-1][2]-1][localTriangleCells[squareNum+NOSC-1][0]-1] += jacobianUT * assemA[2][0];
        A[localTriangleCells[squareNum+NOSC-1][2]-1][localTriangleCells[squareNum+NOSC-1][1]-1] += jacobianUT * assemA[2][1];
        A[localTriangleCells[squareNum+NOSC-1][2]-1][localTriangleCells[squareNum+NOSC-1][2]-1] += jacobianUT * assemA[2][2];

        M[localTriangleCells[squareNum+NOSC-1][0]-1][localTriangleCells[squareNum+NOSC-1][0]-1] += jacobianUT * assemM[0][0];
        M[localTriangleCells[squareNum+NOSC-1][0]-1][localTriangleCells[squareNum+NOSC-1][1]-1] += jacobianUT * assemM[0][1];
        M[localTriangleCells[squareNum+NOSC-1][0]-1][localTriangleCells[squareNum+NOSC-1][2]-1] += jacobianUT * assemM[0][2];
        M[localTriangleCells[squareNum+NOSC-1][1]-1][localTriangleCells[squareNum+NOSC-1][0]-1] += jacobianUT * assemM[1][0];
        M[localTriangleCells[squareNum+NOSC-1][1]-1][localTriangleCells[squareNum+NOSC-1][1]-1] += jacobianUT * assemM[1][1];
        M[localTriangleCells[squareNum+NOSC-1][1]-1][localTriangleCells[squareNum+NOSC-1][2]-1] += jacobianUT * assemM[1][2];
        M[localTriangleCells[squareNum+NOSC-1][2]-1][localTriangleCells[squareNum+NOSC-1][0]-1] += jacobianUT * assemM[2][0];
        M[localTriangleCells[squareNum+NOSC-1][2]-1][localTriangleCells[squareNum+NOSC-1][1]-1] += jacobianUT * assemM[2][1];
        M[localTriangleCells[squareNum+NOSC-1][2]-1][localTriangleCells[squareNum+NOSC-1][2]-1] += jacobianUT * assemM[2][2];

      //  Right hand side of the matrix



        f_lt=GQ.GQ2d(plt);
        f_ut=GQ.GQ2d(put);
//        cout << "Quadrature: " << f_k << endl;
        f[localTriangleCells[squareNum-1][0]] += 2*(pow(pi,2)+1)*cos(pi*plt[0][0])*sin(pi*plt[0][1]) + f_lt ;
        f[localTriangleCells[squareNum-1][1]] += 2*(pow(pi,2)+1)*cos(pi*plt[1][0])*sin(pi*plt[1][1]) + f_lt ;
        f[localTriangleCells[squareNum-1][2]] += 2*(pow(pi,2)+1)*cos(pi*plt[2][0])*sin(pi*plt[2][1]) + f_lt ;

        f[localTriangleCells[squareNum+NOSC-1][0]] += 2*(pow(pi,2)+1)*cos(pi*put[0][0])*sin(pi*put[0][1]) + f_ut ;
        f[localTriangleCells[squareNum+NOSC-1][1]] += 2*(pow(pi,2)+1)*cos(pi*put[1][0])*sin(pi*put[1][1]) + f_ut ;
        f[localTriangleCells[squareNum+NOSC-1][2]] += 2*(pow(pi,2)+1)*cos(pi*put[2][0])*sin(pi*put[2][1]) + f_ut ;


    }

//    MPI_Reduce(&A, &globalA, pow(n+1,4), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    MPI_Reduce(&M, &globalM, pow(n+1,4), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//
//****************   Test for printing positions ***************
//    if (rank==3){
//        for (int i=0;i<=2;i++) {
//                cout <<  plt[i][0] << "  " << plt[i][1]<< "   for upper   "  <<  put[i][0] << "  " << put[i][1]<<endl;}
//    }


//***************  Test for printing Matrix    ****************
//    if (rank == 0){
//
//        for(int i = 0; i < (int)pow(n+1,2); i++){
//
//            for(int j = 0; j< (int)pow(n+1,2); j++){
//
//                cout << globalA[i][j] << " ";
//
//            }
//
//            cout << endl;
//
//        }
//
//    }

    if (rank == 1){
    cout << "Rank:" << rank << endl;
        for (int i = 0; i < pow(SMS,2); i++){
        double squareNum;

        squareNum = n*floor(i/SMS) + (i % SMS) + 1 + PSSV;
        cout << i+1 << ": " << localTriangleCells[squareNum+NOSC-1][0]<< " " << localTriangleCells[squareNum+NOSC-1][1]<< " " << localTriangleCells[squareNum+NOSC-1][2] << endl;
}

    }

    MPI_Finalize();

}
