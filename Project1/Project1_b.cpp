#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <time.h>

using namespace std;
ofstream outputFile;

double sourceFunction(double x) {
    return 100*exp(-10*x);
}

double sourceSolution(double x) {
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}

void gaussianEliminationGeneralSolver(int n, double* v) {       // Construct solver with forward- and backward substitution for the general case
    double h = 1.0/(n+1.0);                                     // Construct steplength, h
    double *a = new double[n+2];                                // Construct constituent a in the tridiagonal matrix A
    double *b = new double[n+2];                                // Construct constituent b in the tridiagonal matrix A
    double *c = new double[n+2];                                // Construct constituent c in the tridiagonal matrix A
    double *f = new double[n+2];

    for (int i=0; i < n+2; i++) {                               // Filling arrays: a, b, c, and bTilde
        f[i] = h*h*sourceFunction(i*h);                         // x[i] = i*h
        a[i] = -1;                                              // Filling a-array
        b[i] = 2;                                               // Filling b-array
        c[i] = -1;                                              // Filling c-array
    }

    for (int i=2; i<=n; i++) {                                  // Forward substitution
        b[i] = b[i] - ( (a[i-1] / b[i-1])*c[i-1] );
        f[i] = f[i]-( ( f[i-1]*a[i-1]) / b[i-1] );
    }

    v[n] = f[n] / b[n];                                         // Constructing the last value of v for backward substitution

    for (int i=n-1; i >= 1; i--) {                              // Backward substitution
        v[i] = ( f[i]-(c[i]*v[i+1]) ) / (b[i]);
    }

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] f;
}

void gaussianEliminationTriDiagonalSolver(int n, double* w) {   // Construct solver with forward- and backward substitution for the general case
    double h = 1.0/(n+1.0);                                     // Construct steplength, h
    double *a = new double[n+2];                                // Construct constituent a in the tridiagonal matrix A
    double *b = new double[n+2];                                // Construct constituent b in the tridiagonal matrix A
    double *c = new double[n+2];                                // Construct constituent c in the tridiagonal matrix A
    double *f = new double[n+2];

    for (int i=0; i < n+2; i++) {                               // Filling arrays: a, b, c, and bTilde
        f[i] = h*h*sourceFunction(i*h);                         // x[i] = i*h
        a[i] = -1.0;                                            // Filling a-array
        b[i] = 2.0;                                             // Filling b-array
        c[i] = -1.0;                                            // Filling c-array
    }

    for (int i=2; i<=n; i++) {                                  // Fill diagonal elements from (i+1)/i
        b[i] = (i+1.0) / i;
    }

    for (int i=2; i<=n; i++) {                                  // fill array for right hand side
        //f[i] = f[i] + 1/b[i-1] * f[i-1];
        f[i] = f[i] + (i-1.)/(i) * f[i-1];                      // 1/(b[i-1]) = 1 / ((i)/(i-1)) = (i-1)/i
    }

    w[n] = f[n] / b[n];                                         // declaring last diagonal element before backward substitution

    for (int i=n-1; i >= 1; i--) {                              // Special backward substitution
        w[i] = (f[i] + w[i+1]) / b[i];
    }

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] f;
}

void calculateRelativeError(int n, double* u, double* v, double* error) {

    for (int i=1; i<=n; i++) {
        error[i] = log10( fabs( (v[i]-u[i]) / u[i] ) );
    }
}

int main() {
    int exponent = 7;
    string fileName = "Result_n=10^";
    //clock_t startClockGeneralSolver, stopClockGeneralSolver, startClockTriDiagonalSolver, stopClockTriDiagonalSolver ;

    for (int i = 1; i <= exponent; i++) {                           // Looping over n = 10^n (10, 100, 1000,...)
        int  n = (int) pow(10.0,i);                                 // Construct n for matrix size (n x n)
        string fileout = fileName;                                  // Declare new file name from fileName
        string argument = to_string(i);                             // Convert the power 10^i to string
        fileout.append(argument + ".txt");                          // Add string(i) and file extention ".txt"

        double *x = new double[n+2];                                // Construct x-array
        double h = 1.0/(n+1.0);                                     // Construct steplength, h
        double *u = new double[n+2];                                // Defining vector u - Analytic solution
        double *v = new double[n+2];                                // Defining vector v - Numerical solution for general case
        double *w = new double[n+2];                                // Defining vector w - Numerical solution for spesial case
        double *error = new double[n+2];

        for (int i=0; i<=n+1; i++) {                                // Filling up x-array
            x[i] = i*h;
        }

        for(int i = 0; i < n+2; i++) {                              // set array values  0
            u[i] = 0;
            v[i] = 0;
            w[i] = 0;
            error[i] = 0;
        }

        for (int i=0; i < n+1; i++) {                               // Filling up vector u from sourceSolution
            u[i] = sourceSolution(x[i]);
        }

        double startClockGeneralSolver = clock();                   // Start timer for general solver

        gaussianEliminationGeneralSolver(n,v);                      // Call forward- and backward substitution for vector v with n elements

        double stopClockGeneralSolver = clock();                    // Stop timer for general solver
        double timeGeneralSolver = (( stopClockGeneralSolver - startClockGeneralSolver) / CLOCKS_PER_SEC);      // Calculating time used on calculation.

        double startClockTriDiagonalSolver = clock();

        gaussianEliminationTriDiagonalSolver(n,w);                  // Call backward substitution for vector w with n elements

        double stopClockTriDiagonalSolver = clock();
        double timeTriDiagonalSolver = ((stopClockTriDiagonalSolver - startClockTriDiagonalSolver) / CLOCKS_PER_SEC);

        calculateRelativeError(n,u,v,error);                        // Calculate relative errover between vector u and v for n elements

/*        cout << "v[n]" << "\t\t" << "w[n]" << "\t\t" << "u[n]"<< "\t\t" << "error[u,v]" << endl;

        for (int i=0; i<n+2; i++ ) {                                // Printing values of v and u
            cout << v[i] << "\t\t" << w[i] << "\t\t" << u[i] << "\t\t" << error[i] << endl;
        }

        for (int i=1; i < n; i++) {                                 // Printing relative error
            cout <<"error:  " <<  error[i] <<"  v[x]:  " <<  v[i] << endl;
        } */

        outputFile.open(fileout);
        outputFile << setiosflags(ios::showpoint | ios::uppercase);
        outputFile << "x" << "\t\t" << "u(x)" << "\t\t" << "v(x)" << "\t\t" << "w(x)" << "\t\t" << "error(u,v)" << "\t" << "time General"  << "\t" << "Time special" << endl;

        for (int i=1; i <= n; i++) {
            outputFile  << setprecision(8) << x[i];
            outputFile << "\t" << setprecision(8) << u[i];
            outputFile << "\t" << setprecision(8) << v[i];
            outputFile << "\t" << setprecision(8) << w[i];
            outputFile << "\t" << setprecision(8) << error[i];
            outputFile << "\t" << setprecision(7) << timeGeneralSolver;
            outputFile << "\t" << setprecision(7) << timeTriDiagonalSolver << endl;
        }

        delete [] x;
        delete [] u;
        delete [] v;
        delete [] w;
        delete [] error;

        outputFile.close();
    }

    return 0;
}
