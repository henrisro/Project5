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
void Forward_Euler(vec &, int, int, double, double, double);
void Backward_Euler(vec &, int, int, double, double&, double&, double&, double, double);
void Crank_Nicolson(vec &, int, int, double, double&, double&, double&, double, double);
void tridiag(double*, double*, double*, vec, vec &, int);
void Output_write(vec &, double, int);

double func(double x) {
	// Use trivial initial condition on the interior of the interval:
	return 0;
}

double func2(double x) {
  // changed boundary conditions for the solution v(x):
  return x-1.0;
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
int n = 50; // meaning Delta_x = 1/10 and Delta_t <= 0.005 for stability:
int t_steps = 5000;
// delta_t determined by the stability condition for the Forward Euler scheme:
double delta_t = 0.50 /((double) n*n);
// use these times for comparison:
double t1 = 0.05; double t2 = 0.5;

// We apply the three methods of interest below:
// Apply Forward Euler method:
ofile.open(outfilename1);
ofile << setiosflags(ios::showpoint | ios::uppercase);
vec u_solution1(n+1);
Forward_Euler(u_solution1, n, t_steps, delta_t, t1, t2);
ofile.close();

// Apply Backward Euler method:
ofile.open(outfilename2);
ofile << setiosflags(ios::showpoint | ios::uppercase);
vec u_solution2(n+1);
double a, b, c;
Backward_Euler(u_solution2, n, t_steps, delta_t, a, b, c, t1, t2);
ofile.close();

// Apply Crank-Nicoloson method:
ofile.open(outfilename3);
ofile << setiosflags(ios::showpoint | ios::uppercase);
vec u_solution3(n+1);
Crank_Nicolson(u_solution3, n, t_steps, delta_t, a, b, c, t1, t2);
ofile.close();
} // End of main function.


// N is the total number of elements in the vectors y and u.
// This implementation is similar to that of Project 1.
void tridiag(double a, double b, double c, vec y, vec & u, int N) 
// A*u = y where A is tridiagonal (b on diagonal)
{
  vec diag_temp(N+1);
  //diag_temp(1) = 0.0;
  double b_temp = b;
  u(1) = y(1)/b_temp;
  for (int i=2; i <= N; i++) {
       // Temporary value needed also in next loop:
      diag_temp(i) = c/b_temp;
      // Temporary diagonal element:
      b_temp = b - a*diag_temp(i);
      // Updating right hand side of matrix equation:
      u(i) = (y(i) - a*u(i-1)) / b_temp;
  }

    // Row reduction; backward substition:
  for (int i=N-1; i >= 1; i--) {
        u(i) -= diag_temp(i+1)*u(i+1);
  }
} // End of function tridiag


void Forward_Euler(vec & u, int n, int tsteps, double delta_t, double t1, double t2)
{
  // n+1 is the number of mesh points in x direction
  // Armadillo is used to handle the vectors.
  vec unew(n+1);
  double delta_x = 1.0/((double) n);
  double alpha = delta_t/(delta_x*delta_x);
  double t_val = 0.0;
  double dummy1 = t1/delta_t; double dummy2 = t2/delta_t;
  int iter1 = (int) dummy1; int iter2 = (int) dummy2;
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
  	if (j == iter1 or j == iter2) {
  		// Write current values (all x) to txt file:
  	  Output_write(unew, t_val, n);
    }
  }
} // End of function Forward_Euler function

void Backward_Euler(vec & u, int n, int tsteps, double delta_t, double &a, double &b, double &c, double t1, double t2) 
{
  vec vnew(n+1); vec unew(n+1);
  vec v(n+1);
  double delta_x = 1.0/((double) n);
  double alpha = delta_t/(delta_x*delta_x);
  double t_val = 0.0;
  double v_val;
  double dummy1 = t1/delta_t; double dummy2 = t2/delta_t;
  int iter1 = (int) dummy1; int iter2 = (int) dummy2;
  // Boundary conditions:
  v(0) = vnew(0) = 0.0;
  v(n) = vnew(n) = 0.0;

  // Zeroth elements not used:
  a = -alpha; 
  b = 1+2*alpha; 
  c = -alpha;
  // Initialize the vector acording to the initial condition in func(x):
  double x_val;
  for (int i=1; i < n; i++) { // Initialize only interior solution
    x_val = i*delta_x; v(i) = func2(x_val); vnew(i) = 0;
  }
  // Time loop:
  for (int j = 1; j <= tsteps; j++) {
    t_val = j*delta_t;
    // For this time, find solution for all values of x.
    // Solver for tridiagonal matrix goes here.
    tridiag(a, b, c, v, vnew, n-1);
    //u.print("u after tridiag: ");
    v(0) = vnew(0) = 0.0;
    v(n) = vnew(n) = 0.0;
    // Update old result for next loop iteration:
    v = vnew;
    if (j == iter1 or j == iter2) {
      for (int l=1; l<n; l++) {
        v_val = vnew(l);
        unew(l) = v_val + 1.0 - l*delta_x;
      }
      unew(0) = 1.0; unew(n) = 0.0;
      // Write current values (all x) to txt file:
      Output_write(unew, t_val, n);
    }
  }
} // end of function Backward_Euler function


void Crank_Nicolson(vec & u, int n, int tsteps, double delta_t, double &a, double &b, double &c, double t1, double t2) 
{
  vec unew(n+1); vec y(n+1);
  vec vnew(n+1); vec v(n+1);
  double delta_x = 1.0/((double) n);
  double alpha = delta_t/(delta_x*delta_x);
  double t_val = 0.0;
  double v_val;
  double dummy1 = t1/delta_t; double dummy2 = t2/delta_t;
  int iter1 = (int) dummy1; int iter2 = (int) dummy2;
  // Boundary conditions:
  v(0) = vnew(0) = 0.0;
  v(n) = vnew(n) = 0.0;

  // Zeroth elements not used:
  a = -alpha; 
  b = 2+2*alpha; 
  c = -alpha;

  double x_val;
  for (int i=1; i < n; i++) { // Initialize only interior solution
    x_val = i*delta_x; v(i) = func2(x_val); vnew(i) = 0;
  }

  // Time loop: 
  for (int j = 1; j <= tsteps; j++) {
    t_val = j*delta_t;
    // For this time, find solution for all values of x.
    // Solver for tridiagonal matrix goes here.
    for (int l=1; l < n; l++) {
         y(l) = alpha*v(l+1) + (2 - 2*alpha)*v(l) + alpha*v(l-1);
    }
    tridiag(a, b, c, y, vnew, n-1);
    //u.print("u after tridiag: ");
    v(0) = vnew(0) = 0.0;
    v(n) = vnew(n) = 0.0;
    // Update old result for next loop iteration:
    v = vnew;
    if (j == iter1 or j == iter2) {
      for (int l=1; l<n; l++) {
        v_val = vnew(l);
        unew(l) = v_val + 1.0 - l*delta_x;
      }
      unew(0) = 1.0; unew(n) = 0.0;
      // Write current values (all x) to txt file:
      Output_write(unew, t_val, n);
    }
  }
} // end of function Crank_Nicolson function


void Output_write(vec & u, double t, int n) 
{
  ofile << "Solution for time: t = " << t << endl;
  for (int i = 0; i <= n; i++) {
     ofile << setw(15) << setprecision(8) << u(i) << endl;
  }
} // end output function