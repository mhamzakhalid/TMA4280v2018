#include <iostream>
#include <cmath>

using namespace std;
#define pi 3.141592653589793238462643383279502

class Gaussian{
    public:

         double GQ1d(double P[2][1]){
            //compute jacobian
            double J[1];
            J[0]=P[1][0] - P[0][0];
            return J[1];
        }

        double GQ2d(double P[3][2]){

            double J[2][2];
            J[0][0] = P[0][0] - P[2][0];
            J[1][0] = P[0][1] - P[2][1];
            J[0][1] = P[1][0] - P[2][0];
            J[1][1] = P[1][1] - P[2][1];
            cout << "Jacobian is: " << endl;
           for (int i=0;i<2;i++){
            for (int j=0;j<2;j++){
                            cout << J[i][j] << "  " ;
                            }
                            cout << endl;
                        }
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
            cout << "Answer is:  ";
            cout << solution << endl;
        }

};

//class determinant{           //Compute determinant
//    public:
//        double det1d(double J[1]){
//            double det;
//            det = 1./J[0];
//            return det;
//            }
//        double det2d(double J[2][2]){
//            double det;
//            det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
//            return det;
//            }
//};
//
//class inverse{            // Compute inverse
//     public:
//
//        double inv1d(double J[1],double det){
//
//            double K[1];
//            K[0] = 1.0 / det;
//            return K[2];
//            }
//
//        double inv2d(double J[2][2],double det){
//
//
//            }
//};
//
//
//class Gaussian{
//    public:
//
//        double gau2d(double P[3][2]){
//            double lambda[3][3];
//            double rho=1./3.;
//            //Making objects from Classes
//
//            Jacobian Jac;
//            determinant det;
//            inverse inv;
//            Func ftn;
//            double K;
//          //  K=abs(det.det2d(Jac.Jac2d(P)))/2;
//            double X[3][2];
//
//            X[0][0]=lambda[0][0]*P[0][0] + lambda[0][1]*P[1][0] + lambda[0][2]*P[2][0];
//            X[0][1]=lambda[0][0]*P[0][1] + lambda[0][1]*P[1][1] + lambda[0][2]*P[2][1];
//            X[1][0]=lambda[1][0]*P[0][0] + lambda[1][1]*P[1][0] + lambda[1][2]*P[2][0];
//            X[1][1]=lambda[1][0]*P[0][1] + lambda[1][1]*P[1][1] + lambda[1][2]*P[2][1];
//            X[2][0]=lambda[2][0]*P[0][0] + lambda[2][1]*P[1][0] + lambda[2][2]*P[2][0];
//            X[2][1]=lambda[2][0]*P[0][1] + lambda[2][1]*P[1][1] + lambda[2][2]*P[2][1];
//
//            double solution;
//
////            solution=K*rho*ftn.ftn2d(X[0][0],X[0][1]) + K*rho*ftn.ftn2d(X[1][0],X[1][1]) +K*rho*ftn.ftn2d(X[2][0],X[2][1]);
//            return solution;
//        }



//};


//class Gaussian{
//    public:
//      double Jacobian1d(double x[1][1]){
//            //compute jacobian
//            double J[1][1];
//            J[0][0]=x[1][0] - x[0][0];
//
//            //Compute determinant
//            double det;
//            det = 1./J[0][0];
//
//            // Compute inverse
//            double K[1][1];
//            K[0][0] = 1.0 / det;
//
//            //return J[0][0],det,K[0][0];
//        }
//
//    //Setting Gaussian quadrature for a 1d problem
//    void Gaussian1d(int a, int b, void f){
//    double z[4],w[4];
//    double result=0.;
//
//     for (int i=0;i<4;i++){
//     result+=(b-a)/2. * w[i]*f((b-a)/2.*z[i]+(a+b)/2.);
//     }
//
//    };
//
//    //for Jacobian in 2d
//    double Jacobian2d(double x[2][2]){
//
//            double J[2][2];
//            J[0][0] = x[1][0] - x[0][0];
//            J[1][0] = x[1][1] - x[0][1];
//            J[0][1] = x[2][0] - x[0][0];
//            J[1][1] = x[2][1] - x[0][1];
//
//            //Compute determinant
//            double det;
//            det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
//
//            // Compute inverse
//            double K[2][2];
//            // Compute inverse
//            K[0][0] = + J[1][1] / det;
//            K[1][0] = - J[1][0] / det;
//            K[0][1] = - J[0][1] / det;
//            K[1][1] = + J[0][0] / det;
//        }
//};

int main()
{
    double y[2][2];
    y[0][0]=0;
    y[0][1]=0;
    y[1][0]=0;
    y[1][1]=1;
    y[2][0]=1;
    y[2][1]=1;
    Gaussian GQ;

    GQ.GQ2d(y);

//    inverse inv;

    //GQ.gau2d(y);
    //cout << "J: " << J[0][0] << " " << det << endl;
    return 0;
}
