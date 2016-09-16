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

int main() {
    int n = 10;                                                 // Construct n for matrix size (n x n)
    double *x = new double[n+2];                                // Construct x-array
    double h = 1.0/(n+1.0);                                     // Construct steplength, h
    double *u = new double[n+2];                                // Defining vector u
    double *v = new double[n+2];                                // Defining vector v
    double *w = new double[n+2];                                // Defining vector w


    for (int i=0; i<=n+1; i++) {                                // Filling up x-array
        x[i] = i*h;
    }

    for(int i = 0; i < n+2; i++) {                              // set array values  0
        u[i] = 0;
        v[i] = 0;
        w[i] = 0;
    }

    gaussianEliminationGeneralSolver(n,v);                      // Call forward- and backward substitution for vector v with n elements
    gaussianEliminationTriDiagonalSolver(n,w);                  // Call backward substitution for vector w with n elements

    for (int i=0; i < n+1; i++) {                               // Filling up vector u from sourceSolution
        u[i] = sourceSolution(x[i]);
    }

    cout << "v[n]" << "\t\t" << "w[n]" << "\t\t" << "u[n]"<< endl;

    for (int i=0; i<n+2; i++ ) {                                // Printing values of v and u
        cout << v[i] << "\t\t" << w[i] << "\t\t" << u[i] << endl;
    }

    outputFile.open("n=10.txt");
    outputFile << setiosflags(ios::showpoint | ios::uppercase);
    outputFile << "\t" << "x" << "\t\t" << "u(x)" << "\t\t" << "v(x)" << "\t\t" << "w(x)" << endl;
    for (int i=1;i<=n;i++) {
        outputFile << setw(15) << setprecision(10) << x[i];
        outputFile << setw(15) << setprecision(10) << u[i];
        outputFile << setw(15) << setprecision(10) << v[i];
        outputFile << setw(15) << setprecision(10) << w[i] << endl;
    }

    outputFile.close();
    delete [] x;
    delete [] u;
    delete [] v;
    delete [] w;
    return 0;
}
