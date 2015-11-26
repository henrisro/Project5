/*
**     Project 5, part g)
**     Implementation of Monte Carlo method for random walk
**     This time we draw the step length (and its direction)
**     from a Gaussian distribution with mean 0 and std 1.
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstring>
#include "lib.h"
#include <time.h>

using namespace std;
ofstream ofile;

// Declaration of functions:
void  output(int, double *, double *, int);
void  mc_sampling(int, int, double *, double *, int, int &);
double gaussian_deviate(long * idum);
// This is not used, but included here as well.
double rand_normal(double, double);

// Initial probability distribution. Area: 1/2 (for normalization in the end)
// This is not directly used in this program.
double p(double x) {
  return 2*(1-x);
}

int main(int argc, char* argv[])
{
char *outfilename1;
int n = 100;
double delta_t = 0.5 /((double) n*n);
// use these times for comparison:
double t1 = 0.05; double t2 = 0.3; // Same time steps as in previous exercise
double dummy1 = t1/delta_t; double dummy2 = t2/delta_t;
int iter1 = (int) dummy1; int iter2 = (int) dummy2; // iter1 and iter2 now contains the number of MC steps
int number_walks;
int N = 0;
double calculation_time;
cout << " iter1: " << iter1 << " and iter2: " << iter2 << endl;
iter1 = iter2;

// Read in output file, abort if there are too few command-line arguments
if( argc <= 2 ) {
  cout << "Bad Usage: " << argv[0] << 
  " Read also outputfile1 and walks on same line!" << endl;
  exit(1);
}
else{
  outfilename1 = argv[1];
  number_walks = atoi(argv[2]);
}

double *walk_position = new double [number_walks];
double *walk_variance = new double [number_walks];
for (int walks = 0; walks < number_walks; walks++){   
  walk_position[walks] = 0.0;
}

ofile.open(outfilename1);
ofile << setiosflags(ios::showpoint | ios::uppercase);

clock_t start, finish;
start = clock();
mc_sampling(iter1, number_walks, walk_position, walk_variance, n, N);
finish = clock();
calculation_time = (finish - start)/(double)CLOCKS_PER_SEC;

ofile << "Number of time steps: " << setw(10) << setprecision(5) << iter1 << endl;
ofile << "Number of walkers: " << setw(10) << setprecision(5) << number_walks << endl;
ofile << "Delta t: " << setw(10) << setprecision(8) << delta_t << endl;
cout << "Final number of escaped walkers: " << N << endl;
output(number_walks, walk_position, walk_variance, iter1);
ofile << "Final number of particles: " << setw(10) << setprecision(8) << N << endl;
ofile << " Calculation time on 'mc_sampling': " << setw(10) << setprecision(8) << calculation_time << endl;
ofile.close();
delete [] walk_position;
return 0; 

} // End of main function.
 
void mc_sampling(int tsteps, int number_walks, 
      double *walk_position, double *walk_variance,
      int n, int & N)
{
  long idum, idum2;
  double delta_t = 0.5 /((double) n*n);
  double l0 = sqrt(2.0*delta_t);
  double r, position, l;
  double eps = 1.0E-5;
  idum = -1;  // initialise random number generator ran0()
  idum2 = 1; // for the gaussian number generator
  int hit_wall;
  int flag;
  double local_average = 0.0; 
  double local_av2 = 0.0;

  N = 0;
  // outer loop over time:
  for (int trial=0; trial < number_walks; trial++){
    // inner loop over walkers:
    hit_wall = 1;
    // Checking whether the walker went outside or not:
    while (hit_wall == 1) {
      position = 0.0;
      r = ran0(&idum); // Used for the mapping.
      position =  1.0 - sqrt(1.0-r); // Mapping to initial probability distribution:
      local_av2 = local_average = 0.0;
      // Finding closest multiple of l0:
      flag = 0;
      for (int tstep = 1; tstep <= tsteps; tstep++) {
        l = l0*gaussian_deviate(&idum2); // Draw with mu = 0 and std = 1.0
        // Move right:
        if (l > 0.0 ){
          if (position + l + l0/2.0 > 1.0) {
            flag = 1; N += 1;
            break;
          }  // Hit right wall: discard!
          else {
            position += l;
          }
        }
        // Move left:
        if (l < 0.0 ){
          if (position + l - l0/2.0 < 0.0) {
             flag = 1; N += 1;
             break;
          }  // Hit left wall: discard!
          else { 
            position += l;
          }
        }
        local_average += position; local_av2 += position*position;
      }  // end of loop over walks
      if (flag == 0){ hit_wall = 0; } // Walk ended without hitting a wall!
    } // end of loop over trials
    // Updating variance array and averages:
    local_average /= ((double) tsteps); local_av2 /= ((double) tsteps);
    walk_position[trial] = position;
    walk_variance[trial] = local_av2 - local_average*local_average;
  }
}   // end mc_sampling function  


void output(int number_walks, double *walk_position, double *walk_variance, int tsteps)
{
  for( int  i = 0; i < number_walks; i++) {
    ofile << setw(15) << setprecision(8) << walk_position[i];
    ofile << setw(15) << setprecision(8) << walk_variance[i] << endl;;
  }
}  // end of function output 

// random number with gaussian distribution:
double gaussian_deviate(long * idum) {
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;

  if (idum < 0) {iset = 0;}
  if (iset == 0) {
    do {
      v1 = 2.0*ran2(idum) -1.0;
      v2 = 2.0*ran2(idum) -1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0);
    fac = sqrt(-2.0*log(rsq)/rsq);
    gset = v1*fac;
    iset = 1;
    return v2*fac;
  } else {
    iset = 0;
    return gset;
  }
}// end of gaussian_deviate function

// Alternative of the same function:
double rand_normal(double mean, double stddev) {
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached) {
        double x, y, r;
 do {
     x = 2.0*rand()/RAND_MAX - 1;
     y = 2.0*rand()/RAND_MAX - 1;

     r = x*x + y*y;
 } while (r == 0.0 || r > 1.0);
        {
        double d = sqrt(-2.0*log(r)/r);
        double n1 = x*d;
        n2 = y*d;
        double result = n1*stddev + mean;
        n2_cached = 1;
        return result;
        }
    } else {
        n2_cached = 0;
        return n2*stddev + mean;
    }
} // end rand_normal function

