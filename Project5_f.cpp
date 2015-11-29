/*
**     Project 5, part f)
**     Implementation of Monte Carlo method for random walk.
**     Uniform distribution of steps (equal probability in each
**     direction). 
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
void  mc_sampling(int, int, double, double *, double *, double, int, int &);

// Initial probability distribution. Area: 1/2 (for normalization in the end)
// This is not directly used in the program below.
double p(double x) {
  return 2*(1-x);
}

int main(int argc, char* argv[])
{
char *outfilename1;
int n = 100;
double delta_t = 0.5 /((double) n*n);
double l0 = sqrt(2*delta_t);  // step length in Diffusion problem.
// use these times for comparison:
double t1 = 0.05; double t2 = 0.3; // Same time steps as in previous exercise
double dummy1 = t1/delta_t; double dummy2 = t2/delta_t;
int iter1 = (int) dummy1; int iter2 = (int) dummy2; // iter1 and iter2 now contains the number of MC steps
double move_probability = 0.5;
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
mc_sampling(iter1, number_walks, move_probability, walk_position, walk_variance, l0, n, N);
finish = clock();
calculation_time = (finish - start)/(double)CLOCKS_PER_SEC;

ofile << "Number of time steps: " << setw(10) << setprecision(5) << iter1 << endl;
ofile << "Number of walkers: " << setw(10) << setprecision(5) << number_walks << endl;
ofile << "Step length: " << setw(10) << setprecision(8) << l0 << endl;
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
      double move_probability, double *walk_position,
      double *walk_variance, double l0, int n, int & N)
{
  long idum;
  double r, position, r2, y;
  double eps = 1.0E-5;
  idum = -1;  // initialise random number generator (ran0)
  long idum2 = -2;
  double dummy, val;
  int hit_wall;
  double *possible_pos = new double [n-1];
  int flag;
  double local_average = 0.0; 
  double local_av2 = 0.0;

  for (int i=1; i < n; i++) {
    possible_pos[i] = i*l0;
    //cout << " Possible position: " << possible_pos[i] << endl;
  }
  N = 0;
  // outer loop over walkers:
  for (int trial=0; trial < number_walks; trial++){
    // inner loop over walkers:
    hit_wall = 1;
    // Checking whether the walker went outside or not:
    while (hit_wall == 1) {
      position = 0.0;
      dummy = 2.0; // Some random (and large) number to ensure new position is found
      r2 = ran0(&idum2);
      // Mapping to initial probability distribution:
      y = 1.0 - sqrt(1.0-r2);
      local_av2 = local_average = 0.0;
      // Finding closest multiple of l0:
      for (int j=1; j<n; j++) {
        val = fabs(possible_pos[j] - y);
        if (val < dummy) { dummy = val; position = possible_pos[j]; }
      }
      flag = 0;
      for (int tstep = 1; tstep <= tsteps; tstep++) {
        r = ran0(&idum);
        // Move right:
        if (r <= move_probability){
          if (1.0 - l0 - position <= eps){flag = 1; N += 1; break; } // Hit right wall: discard!
          else{ position += l0; }
        }
        // Move left:
        if (r > move_probability){
          if (position - l0 <= eps){ flag = 1; N += 1; break; } // Hit left wall: discard!
          else { position -= l0; }
        }
        local_average += position; local_av2 += position*position;
      }  // end of loop over walks
      if (flag == 0){ hit_wall = 0; } // Walk ended without hitting a wall!
    } // end of loop over trials

    local_average /= ((double) tsteps); local_av2 /= ((double) tsteps);
    walk_position[trial] = position;
    walk_variance[trial] = local_av2 - local_average*local_average;
  }
  delete [] possible_pos;
}   // end mc_sampling function  


void output(int number_walks, double *walk_position, double *walk_variance, int tsteps)
{
  for( int  i = 0; i < number_walks; i++) {
    ofile << setw(15) << setprecision(8) << walk_position[i];
    ofile << setw(15) << setprecision(8) << walk_variance[i] << endl;
  }
}  // end of function output 