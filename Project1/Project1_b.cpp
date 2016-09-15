#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;
ofstream outputFile;

double sourceFunction(double x) {
    return 100*exp(-10*x);
}

double sourceSolution(double x) {
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}

void generalSolver(int n, double* v) {              // Construct solver with forward- and backward substitution for the general case
    double h = 1.0/(n+1.0);                         // Construct steplength, h
    double *a = new double[n+2];                    // Construct constituent a in the tridiagonal matrix A
    double *b = new double[n+2];                    // Construct constituent b in the tridiagonal matrix A
    double *c = new double[n+2];                    // Construct constituent c in the tridiagonal matrix A
    double *bTilde = new double[n+2];

    for (int i=0; i < n+2; i++) {                   // Filling arrays: a, b, c, and bTilde
        bTilde[i] = h*h*sourceFunction(i*h);        // x[i] = i*h
        a[i] = -1;                                  // Filling a-array
        b[i] = 2;                                   // Filling b-array
        c[i] = -1;                                  // Filling c-array
    }

    for (int i=2; i<=n; i++) {                      // Forward substitution
        b[i] = b[i] - ( (a[i-1] / b[i-1])*c[i-1] );
        bTilde[i] = bTilde[i]-( ( bTilde[i-1]*a[i-1]) / b[i-1] );
    }

    v[n] = bTilde[n] / b[n];                        // Constructing the last value of v for backward substitution

    for (int i=n-1; i >= 1; i--) {                  // Backward substitution
        v[i] = ( bTilde[i]-(c[i]*v[i+1]) ) / (b[i]);
    }

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] bTilde;
}

int main() {
    int n = 10;                                 // Construct n for matrix size (n x n)
    double *x = new double[n+2];                // Construct x-array
    double h = 1.0/(n+1.0);                     // Construct steplength, h
    double *v = new double[n+2];                // Defining vector v
    double *u = new double[n+2];                // Defining vector u

    for (int i=0; i<=n+1; i++) {                // Filling up x-array
        x[i] = i*h;
    }

    for(int i = 0; i < n+2; i++) {              // set array values  0
        v[i] = 0;
        u[i] = 0;
    }

    generalSolver(n,v);                         // Call forward- and backward substitution vector v with n elements

    for (int i=0; i < n+2; i++) {               // Filling up vector u from sourceSolution
        u[i] = sourceSolution(x[i]);
    }

    cout << "v[n]" << "\t" << "u[n]"<< endl;

    for (int i=0; i<n+2; i++ ) {                // Printing values of v and u
        cout << v[i] << "\t" << u[i] << endl;
    }

    outputFile.open("n=10.txt");
    outputFile << setiosflags(ios::showpoint | ios::uppercase);
    outputFile << "\t" << "x" << "\t\t" << "u(x)" << "\t\t" << "v(x)" << endl;
    for (int i=1;i<=n;i++) {
        outputFile << setw(15) << setprecision(10) << x[i];
        outputFile << setw(15) << setprecision(10) << u[i];
        outputFile << setw(15) << setprecision(10) << v[i] << endl;
    }

    outputFile.close();
    delete [] x;
    delete [] v;
    delete [] u;
    return 0;
}
