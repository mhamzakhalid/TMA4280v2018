#include <iostream>
#include <math.h>

using namespace std;

class Gaussian{
    public:
      double Jacobian1d(double x[1][1]){
            //compute jacobian
            double J[1][1];
            J[0][0]=x[1][0] - x[0][0];

            //Compute determinant
            double det;
            det = 1./J[0][0];

            // Compute inverse
            double K[1][1];
            K[0][0] = 1.0 / det;

            //return J[0][0],det,K[0][0];
        }

    //Setting Gaussian quadrature for a 1d problem
    void Gaussian1d(int a, int b, void f){
    double z[4],w[4];
    double result=0.;

     for (int i=0;i<4;i++){
     result+=(b-a)/2. * w[i]*f((b-a)/2.*z[i]+(a+b)/2.);
     }

    };

    //for Jacobian in 2d
    double Jacobian2d(double x[2][2]){

            double J[2][2];
            J[0][0] = x[1][0] - x[0][0];
            J[1][0] = x[1][1] - x[0][1];
            J[0][1] = x[2][0] - x[0][0];
            J[1][1] = x[2][1] - x[0][1];

            //Compute determinant
            double det;
            det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

            // Compute inverse
            double K[2][2];
            // Compute inverse
            K[0][0] = + J[1][1] / det;
            K[1][0] = - J[1][0] / det;
            K[0][1] = - J[0][1] / det;
            K[1][1] = + J[0][0] / det;
        }
};

int main()
{
    double y[1][1];
    y[0][0]=10;
    y[0][1]=2;

    Gaussian GQ;
    GQ.Jacobian1d(y);
    //cout << "J: " << J[0][0] << " " << det << endl;
    return 0;
}
