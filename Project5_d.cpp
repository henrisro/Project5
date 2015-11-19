/*
**     Project5: 
**     Implementation of the three methods
**     1) Forward Euler explicit scheme
**     2) Backward Euler implicit scheme
**     3) Crank-Nicolson implicit scheme
*/

#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <time.h>

using namespace std;
using namespace arma;
ofstream ofile;

// Declaration of functions:
void Forward_Euler(vec &, int, int, double);
void Backward_Euler(vec &, int, int, double);
void Crank_Nicolson(vec &, int, int, double);
void tridiag(double *, double *, double *, vec &, vec &, int);
void Output_write(vec &, double, int);

double func(double x) {
	// Use trivial initial condition on the interior of the interval:
	return 0;
}

int main(int argc, char* argv[])
{
char *outfilename1, *outfilename2, *outfilename3;

// Read in output file, abort if there are too few command-line arguments
if( argc <= 2 ){
  cout << "Bad Usage: " << argv[0] << 
  " Read also outputfile 1, 2 and 3 on same line!" << endl;
  exit(1);
}
else{
  outfilename1 = argv[1];
  outfilename2 = argv[2];
  outfilename3 = argv[3];
}

// Apply the Forward Euler scheme to test:
int n = 59; // meaning Delta_x = 1/10 and Delta_t <= 0.005 for stability:
int t_steps = 1000;
// delta_t determined by the stability condition for the Forward Euler scheme:
double delta_t = 0.45 /(double (n+1)*(n+1));

// We apply the three methods of interest below:
// Apply Forward Euler method:
ofile.open(outfilename1);
ofile << setiosflags(ios::showpoint | ios::uppercase);
vec u_solution1(n+1);
Forward_Euler(u_solution1, n, t_steps, delta_t);
ofile.close();

// Apply Backward Euler method:
ofile.open(outfilename2);
ofile << setiosflags(ios::showpoint | ios::uppercase);
vec u_solution2(n+1);
Backward_Euler(u_solution2, n, t_steps, delta_t);
ofile.close();

// Apply Crank-Nicoloson method:
ofile.open(outfilename3);
ofile << setiosflags(ios::showpoint | ios::uppercase);
vec u_solution3(n+1);
Crank_Nicolson(u_solution3, n, t_steps, delta_t);
ofile.close();
} // End of main function.


// N is the total number of elements in the vectors y and u.
// This implementation is similar to that of Project 1.
void tridiag(double *a, double *b, double *c, vec & y, vec & u, int N) 
// A*u = y where A is tridiagonal (b on diagonal)
{
  vec diag_temp(N+1);
  //diag_temp(1) = 0.0;
  double b_temp = b[1];
  u(1) = y(1)/b_temp;
  for (int i=2; i <= N; i++) {
       // Temporary value needed also in next loop:
      diag_temp(i) = c[i-1]/b_temp;
      // Temporary diagonal element:
      b_temp = b[i] - a[i]*diag_temp(i);
      // Updating right hand side of matrix equation:
      u(i) = (y(i) - u(i-1)*a[i])/b_temp;
  }

    // Row reduction; backward substition:
  for (int i=N-1; i >= 1; i--) {
        u(i) -= diag_temp(i+1)*u(i+1);
  }
} // End of function tridiag


void Forward_Euler(vec & u, int n, int tsteps, double delta_t)
{
  // n+1 is the number of mesh points in x direction
  // Armadillo is used to handle the vectors.
  //vec u(n+1);
  vec unew(n+1);
  double delta_x = 1.0/((double) (n+1));
  double alpha = delta_t/(delta_x*delta_x);
  double t_val = 0.0;
  int result_print = tsteps/5;
  // Boundary conditions:
  u(0) = unew(0) = 1.0;
  u(n) = unew(n) = 0.0;
  // Initialize the vector acording to the initial condition in func(x):
  double x_val = 0.0;
  for (int i=1; i < n; i++) { // Initialize only interior solution
  	x_val = i*delta_x;
  	u(i) = func(x_val);
  	unew(i) = 0;
  }
  // Time integration according to the explicit scheme:
  // Loop over time:
  for (int j = 1; j <= tsteps; j++) {
  	t_val = j*delta_t;
  	// For this time, find solution for all values of x:
  	for (int i = 1; i < n; i++) {
  		unew(i) = alpha*u(i-1) + (1 - 2*alpha)*u(i) + alpha*u(i+1);
  	}
  	u = unew;
  	if (j % result_print == 0) {
  		// Write current values (all x) to txt file:
  		Output_write(unew, t_val, n);
    }
  }
} // End of function Forward_Euler function

void Backward_Euler(vec & u, int n, int tsteps, double delta_t) 
{
  vec unew(n+1);
  double delta_x = 1.0/((double) (n+1));
  double alpha = delta_t/(delta_x*delta_x);
  double t_val = 0.0;
  int result_print = tsteps/5;
  // Boundary conditions:
  u(0) = unew(0) = 1.0;
  u(n) = unew(n) = 0.0;

  // Zeroth elements not used:
  double *a = new double[n+1];
  double *b = new double[n+1];
  double *c = new double[n+1];
  for (int i=1; i<=n; i++) {
        b[i] = 1+2*alpha;
        a[i] = -alpha;
        c[i] = -alpha;
  }
  a[1] = 0.0; c[n] = 0.0;
  // Initialize the vector acording to the initial condition in func(x):
  double x_val;
  for (int i=1; i < n; i++) { // Initialize only interior solution
    x_val = i*delta_x; u(i) = func(x_val); unew(i) = 0;
  }
  // Time loop:
  for (int j = 1; j <= tsteps; j++) {
    t_val = j*delta_t;
    // For this time, find solution for all values of x.
    // Solver for tridiagonal matrix goes here.
    // Correct right hand side:
    u(1) += alpha; // A simple way to take care of the boundary conditions.
    tridiag(a, b, c, u, unew, n);
    //u.print("u after tridiag: ");
    u(0) = unew(0) = 1.0;
    u(n) = unew(n) = 0.0;
    // Update old result for next loop iteration:
    u = unew;
    // Taking care of the boundary conditions for the next iteration: 
    if (j % result_print == 0) {
      // Write current values (all x) to txt file:
      Output_write(unew, t_val, n);
    }
  }
  delete [] a;
  delete [] b;
  delete [] c;
} // end of function Backward_Euler function


void Crank_Nicolson(vec & u, int n, int tsteps, double delta_t) 
{
  vec unew(n+1);
  double delta_x = 1.0/((double) (n+1));
  double alpha = delta_t/(delta_x*delta_x);
  double t_val = 0.0;
  int result_print = tsteps/5;
  // Zeroth elements not used:
  double *a = new double[n+1];
  double *b = new double[n+1];
  double *c = new double[n+1];
  for (int i=1; i<=n; i++) {
        b[i] = 2+2*alpha;
        a[i] = -alpha;
        c[i] = -alpha;
  }
  a[1] = 0.0; c[n] = 0.0;
  // Initialize the vector acording to the initial condition in func(x):
  double x_val = 0.0;
  for (int i=1; i < n; i++) { // Initialize only interior solution
    x_val = i*delta_x;  u(i) = func(x_val);  unew(i) = 0;
  }
  u(0) = unew(0) = 1.0;
  u(n) = unew(n) = 0.0;
  vec y(n+1);
  // Time loop: 
  for (int j = 1; j <= tsteps; j++) {
    t_val = j*delta_t;
    // For this time, find solution for all values of x.
    // Solver for tridiagonal matrix goes here.
    // Calculate the right hand side in the Crank-Nicoloson scheme:
    for (int l=1; l < n; l++) {
        y(l) = alpha*u(l+1) + (2 - 2*alpha)*u(l) + alpha*u(l-1);
    }
    y(0) = (2-2*alpha)*u(0) + alpha*u(1); 
    y(n) = (2-2*alpha)*u(n) + alpha*u(n-1);
    // Impose boundary conditions on this one as well:
    y(1) += alpha;
    tridiag(a, b, c, y, unew, n);
    //if (j == 1) { unew.print("u_new after one time step: "); }
    // Impose boundary conditions to be sure:
    u(0) = unew(0) = 1.0;
    u(n) = unew(n) = 0.0;
    // Update old result for next loop iteration:
    u = unew;

    if (j % result_print == 0) {
      // Write current values (all x) to txt file:
      Output_write(unew, t_val, n);
    }
  }
  delete [] a;
  delete [] b;
  delete [] c;
} // end of function Crank_Nicolson function


void Output_write(vec & u, double t, int n) 
{
  ofile << "Solution for time: t = " << t << endl;
  for (int i = 0; i <= n; i++) {
     ofile << setw(15) << setprecision(8) << u(i) << endl;
  }
} // end output function

// double Solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}
// double f(double x) {return 100*exp(-10*x);}

// TEST THE TRIDIAG FUNCTION:
// int N = 100;
// double h = 1.0/(N+1.0);
// double *x = new double[N+2];
// vec y(N+1);

// // The constituents of the tridiagonal matrix A:
// // Zeroth element not needed, but included to make indexing easy:
// double *a = new double[N+1];
// double *b = new double[N+1];
// double *c = new double[N+1];
// // Filling up x-array:
// for (int i=0; i<=N+1; i++) {
//         // Making x[0] = 0 and x[n+1] = 1:
//         x[i] = i*h;
// }

// // Filling up b_twiddle array, i.e. right hand side of equation:
// for (int i=1; i<=N; i++) {
//         y(i) = h*h*f(x[i]);
//         b[i] = 2.0;
//         a[i] = -1.0;
//         c[i] = -1.0;
// }
// c[N] = 0;
// a[1] = 0;
// //for (int j=1; j<=N; j++) {
// //  cout << a[j] << " " << b[j] << " " << c[j] << " " << y(j) << endl;
// //}

// vec u = zeros(N+1); 
// tridiag(a, b, c, y, u, N);
// u.print(" solution: ");

//Test tridiagonal on trivial example: 
// double a[] = {0.0, 0.0, 2.0};
// double c[] = {0.0, 1.0, 0.0};
// double b[] = {0.0, 1.0, 1.0};
// vec solution(3);
// vec LHS(3);
// LHS(0) = 0.0; LHS(1) = 1.0; LHS(2) = 0.0;
// LHS.print(" LHS before: ");
// tridiag(a, b, c, LHS, solution, 2);
// solution.print(" Solution: ");
