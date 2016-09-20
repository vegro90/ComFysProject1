#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

// object for output files
ofstream outputFile;
// Functions used
double f(double x){return 100.0*exp(-10.0*x);
}
double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

// Begin main program
int main( ){
    int exponent = 3;
    string fileName = "Result_n=10^";
    //clock_t startClockGeneralSolver, stopClockGeneralSolver, startClockTriDiagonalSolver, stopClockTriDiagonalSolver ;

    for (int i = 1; i <= exponent; i++) {                           // Looping over n = 10^n (10, 100, 1000,...)
        int  n = (int) pow(10.0,i);                                 // Construct n for matrix size (n x n)
        string fileout = fileName;                                  // Declare new file name from fileName
        string argument = to_string(i);                             // Convert the power 10^i to string
        fileout.append(argument + ".txt");                          // Add string(i) and file extention ".txt"

      double h = 1.0/(n);
      double hh = h*h;
      n = n-1;  //  shift so that only points between endpoints are studied
      mat A = zeros<mat>(n,n);
      // Set up arrays for the simple case
      vec b(n);  vec x(n);
      A(0,0) = 2.0;  A(0,1) = -1;  x(0) = h;  b(0) =  hh*f(x(0));
      x(n-1) = x(0)+(n-1)*h; b(n-1) = hh*f(x(n-1));
      for (int i = 1; i < n-1; i++){
        x(i) = x(i-1)+h;
    b(i) = hh*f(x(i));
        A(i,i-1)  = -1.0;
        A(i,i)    = 2.0;
        A(i,i+1)  = -1.0;
      }
      A(n-1,n-1) = 2.0; A(n-2,n-1) = -1.0; A(n-1,n-2) = -1.0;
  // solve Ax = b
      vec solution  = solve(A,b);
      cout << x(i) << solution(i) << exact(x(i)) << endl;
      /*
      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      for (int i = 0; i < n ;i++) {
        double RelativeError = fabs((exact(x(i))-solution(i))/exact(x(i)));
        ofile << setw(15) << setprecision(8) << x(i);
        ofile << setw(15) << setprecision(8) << solution(i);
        ofile << setw(15) << setprecision(8) << exact(x(i));
        ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
*/
    }
    return 0;
}
