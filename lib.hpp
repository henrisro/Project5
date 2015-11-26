#ifndef LIB_LIBRARY_IS_INCLUDED
#define LIB_LIBRARY_IS_INCLUDED
 
    /*
     * The definition module                              
     *                      lib.hpp                    
     * for the library function common for all C++ programs.
     */

//    IMPORTANT: you need to change the path
//    NOTE : at the lab you need to use /site/Blitz++-0.9_64 as 
//    path where Blitz++ is stored. All PCs at the lab are 64 bits
//    machines
//   include "/site/Blitz++-0.9_64/include/blitz/array.h"

#include "/usr/include/blitz-0.9/blitz/array.h"

#include <iostream>          // Standard ANSI-C++ include files
#include <new>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream> 
#include <iomanip> 
#include <time.h>
#include <ctype.h>
#include <sys/time.h>

 using namespace std;
 using namespace blitz;

#define   NULL_PTR       (void *) 0
#define   ZERO           0
#define   D_ZERO         1.0E-10
#define   UI             unsigned int


         /* a macro used in function pythag() */

static float sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)


     /* Macro definitions for integer arguments only */

#define   SIGN(a,b) ((b) < 0 ? -fabs(a) : fabs(a))

    // ******   data declaration  ******* 



     // Function declarations 
template <typename T>
T inline scalarVectorProd(Array<T,1>& A,  Array<T,1>& B);
       /*
       ** calculates the product 
       **   scalarValue = Array<T,1> A * Array<T,1> B
       ** using Blitz array package
       */
template <typename T>
void  inline matrixVectorProd(Array<T,2>& A,  Array<T,1>& B, Array<T,1> C);
       /*
       ** calculates the product 
       **   Array<T,1> C = Array<T,2> A * Array<T,1> B
       ** using Blitz array package
       */
template <typename T>
void  inline matrixMatrixProd(Array<T,2>& A,  Array<T,2>& B, Array<T,2> C);
       /*
       ** calculates the product 
       **   Array<T,2> C = Array<T,2> A * Array<T,2> B
       ** using Blitz array package
       */


void psdes(unsigned long& lword, unsigned long& irword);
      /*
      ** performs  a "Pseudo_DES" hashing  of a 64-bit word (lword, irword).
      ** both 32-bit arguments are returned hashed on all bits.
      */

//template <typename T>        void ludcmp(T ** , int, int * , T * );
//template <typename T>        void lubksb(T **, int, int *, T *);
//template <typename T>        T     ran0(int&);


template <typename T> void ludcmp(Array<T,2>& A, int n, Array<int,1>& Index, T& d);
      /*
      ** takes as input a BLITZ matrix Array<T,2> A() of dimension n and
      ** replaces it by the LU decomposition of a rowwise permutation of
      ** itself. The results is stored in A() in the form given by 
      ** eq. (2.3.14) in "Numerical Recipe", sect. 2.3, page 45. The vector
      ** Array<T,1> Index() records the row permutation effected by the
      ** partial pivoting;
      ** The parameter d is output as +1 or -1 depending on whether the 
      ** number of row interchanges was even or odd, respectively. This
      ** routine is used in combination with the template function lubksb()
      ** to solve linear equations or invert a matrix.
      ** The function is slightly modified from the version in Numerical
      ** recipe
      */

template <typename T> void lubksb(Array<T,2>& A, int n, Array<int,1>& Index, Array<T,1>&  B);
      /*
      ** solves the set of linear equations 
      **           A X = B 
      ** of dimension n. Array<T,2> A() is input, not as the
      ** matrix A() but rather as its LU decomposition, 
      ** determined by the template function ludcmp(),
      ** Array<int,1> Index() is input as the permutation vector
      ** returned by ludcmp(). Array<T,1> B() is input as the
      **  right-hand side vector B,
      ** The solution X is returned in B(). The input data A(),
      ** n and Index() are not modified. This routine take into 
      ** account the possibility that B() will begin with many
      ** zero elements, so it is efficient for use in matrix
      ** inversion.
      ** The function is slightly modified from the version in 
      ** in Numerical recipe.
      */
void  chsone(Array<double,1>& bins, Array<double,1>& ebins, const int knstrn, 
                                        double& df, double& chsq, double& prob);
      /*
      ** takes two Blitz arrays, Array<double,1> bins() containing the observed 
      ** number of events and Array<double,1> ebins() containing expected number
      ** of events. Given the number of constrains knstrn (normally one). this 
      ** routine returns (trivally) the number of degrees of freedom df, and 
      ** (non-trivially) the chi-square chsq and the significance prob. A small 
      ** value of prob indicates a significant difference the distributions 
      ** bins() and ebins(). Note that bins() and ebins() are both arrays, 
      ** although bins will normally contain integer values.
      */   

void  chstwo(Array<double,1>& bins1, Array<double,1>& bins2, const int knstrn, 
                                        double& df, double& chsq, double& prob);
      /*
      ** takes two Blitz arrays,  Array<double,1> bins1() and Array<double,1> bins2()
      **  containing two sets of binned data. Given the number of constrains 
      ** knstrn (normally one or two). this routine returns the number of degrees of
      ** freedom df, the chi-square chsq and the significance prob. A small 
      ** value of prob indicates a significant difference between the distributions 
      ** bins1() and bins2(). Note that bins1() and bins2() are both arrays, 
      ** although they will normally contain integer values.
      */  

double gammq(const double a, const double x);
      /*
      ** returns the incomplete gamma function 
      **       Q/a,x) = 1 - P/a,x)
      */

void gser(double& gamser, const double a, const double x, double &gln);
      /*
      ** returns the incomplete gamma function P(a,x) 
      ** evaluated by its series representaion as
      ** gamser 
      ** It also returns ln(gamma(a)) as gln
      */

double gammln(const double xx);
      /*
      ** return the value ln[gamma(xx)]
      */

void gcf(double& gammcf, const double a, const double x, double &gln);
     /*
     ** returns the incomplete gamma function Q(a,x)
     ** evaluated by its continued fraction representation
     ** as gammacf.
     ** Also returns ln gamma(a)  as gln
     */

void gauleg(double x1, double x2, Array<double,1>& x, Array<double,1>& w, int n);
     /*
     ** The function 
     **              gauleg()
     ** takes the lower and upper limits of integration x1, x2, calculates
     ** and return the abcissas in x[0,...,n - 1] and the weights in 
     **w[0,...,n - 1] of length n of the Gauss--Legendre n--point 
     ** quadrature formulae.
     */

template<typename T>
inline void rot(Array<T,2>& a, const T s, const T tau, const int i,
		const int j, const int k, const int l);

template<typename T>
void jacobi(Array<T,2>& a, Array<T,1>& d, Array<T,2>& v, int &nrot);

template<typename T>
void jacobn_s(const T x, Array<T,1>& y, Array<T,1>& dfdx, Array<T,2>& dfdy);


template<typename T>
void derivs_s(const T x, Array<T,1>& y, Array<T,1>& dydx);

template<typename T>
void eigsrt(Array<T,1>& d, Array<T,2>& v);
 
template<typename T> 
void tqli(Array<T,1>& d, Array<T,1>& e, Array<T,2>& z);
   /*
   ** determine eigenvalues and eigenvectors of a real symmetric
   ** tri-diagonal matrix, or a real, symmetric matrix previously
   ** reduced by function tred2() to tri-diagonal form. On input,
   ** Array<T,1> d() contains the diagonal element and 
   ** Array<T,1>e() the sub-diagonal of the tri-diagonal matrix.
   ** On output d[] contains the eigenvalues and  Array<T,1> e()
   ** is destroyed. If eigenvectors are desired Array<T,2>z( , )
   ** on input contains the identity matrix. If eigenvectors of a 
   ** matrix reduced by tred2() are required, then Array<T,2>z( , )
   ** on input is the matrix output from tred2().
   ** On output, the k'th column returns the normalized eigenvector
   ** corresponding to Array<T,1> d(k). 
   ** The function is modified from the version in Numerical recipe.
   */

template<typename T> 
void tqli(Array<T,1>& d, Array<T,1>& e);

template<typename T>
void tred2(Array<T,2>& a, Array<T,1>& d, Array<T,1>& e);
    /*
    ** perform a Housholder reduction of a real symmetric matrix
    ** Array<T,2> a(). On output Array<T,2> a() is replaced by 
    ** the orthogonal matrix effecting the transformation. 
    ** Array<T,1> d() returns the diagonal elements of the 
    ** tri-diagonal matrix, and Array<T,1> e() the off-diagonal 
    ** elements,  with Array<T,1> e() = 0.
    ** The function is modified from the version in Numerical recipe.
    */

template<typename T>
T  pythag(T a, T b);

void rk4(Array<double,2>& y, Array<double,2>& dydx, int n, double x, double h, 
                 Array<double,2>& yout,
	 void (*derivs)(double, Array<double,2>&, Array<double,2>&));
      /*
      ** takes a set of variables Array<double,1> y(0 : n-1) for the function y(x)
      ** together with the derivatives  Array<double,1> dydx(0:n-1) and uses the
      ** fourth-order Runge-Kutta method to advance the solution over an interval
      ** h and return incremented variables as Array<double,1> yout(0:n-1), which
      ** not need to be a distinct array from y[1:n]. The users supply the routine
      ** derivs(double x,Array<double,1> y, Array<double,1> dydx), which returns 
      ** the derivatives dydx at x.
      */ 

void spline(Array<double,1>& x, Array<double,1>& y, int n, double& yp1, 
	    double yp2, Array<double,1> y2);
         /*
         ** takes as input x(0,..,n - 1) and y(0,..,n - 1) containing a tabulation
         ** y_i = f(x_i) with x_0 < x_1 < .. < x_(n - 1) together with yp_1 and yp2
         ** for first derivatives  f(x) at x_0 and x_(n-1), respectively. Then the
         ** function returns y2(0,..,n-1) which contanin the second derivatives of
         ** f(x_i)at each point x_i. If yp1 and/or yp2 is larger than the constant
         ** INFINITY the function will put corresponding second derivatives to zero.
         */ 
void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2, 
	    int n, double x, double& y);
     /*
     ** takes xa(0,..,n - 1) and y(0,..,n - 1) which tabulates a function 
     ** (with the xa(i)'s in order) and given ya(0,..,n - 1), which is the
     ** output from function spline() and with given value of x returns a 
     ** cubic--spline interpolation value y.
     */

void polint(Array<double,1>& xa, Array<double,1>& ya, int n, double x,
            Array<double,1>& y, Array<double,1>& dy);

   /*
   ** takes as input xa(0,..,n-1) and ya(0,..,n-1) together with a given value
   ** of x and returns a value y and an error estimate dy. If P(x) is a polynomial
   ** of degree N - 1 such that P(xa_i) = ya_i, i = 0,..,n-1, then the returned 
   ** value is y = P(x). 
   */


double ran0(long *idum);
/*
    ** is an "Minimal" random number generator of Park and Miller
    ** (see Numerical recipe page 279). Set or reset the input value
    ** idum to any integer value (except the unlikely value MASK)
    ** to initialize the sequence; idum must not be altered between
    ** calls for sucessive deviates in a sequence.
    ** The function returns a uniform deviate between 0.0 and 1.0.
    */
double ran1(long *idum);
/*
    ** is an "Minimal" random number generator of Park and Miller
    ** (see Numerical recipe page 280) with Bays-Durham shuffle and
    ** added safeguards. Call with idum a negative integer to initialize;
    ** thereafter, do not alter idum between sucessive deviates in a
    ** sequence. RNMX should approximate the largest floating point value
    ** that is less than 1.
    ** The function returns a uniform deviate between 0.0 and 1.0
    ** (exclusive of end-point values).
    */
double ran2(long *idum);
/*
    ** is a long periode (> 2 x 10^18) random number generator of 
    ** L'Ecuyer and Bays-Durham shuffle and added safeguards.
    ** Call with idum a negative integer to initialize; thereafter,
    ** do not alter idum between sucessive deviates in a
    ** sequence. RNMX should approximate the largest floating point value
    ** that is less than 1.
    ** The function returns a uniform deviate between 0.0 and 1.0
    ** (exclusive of end-point values).
    */ 
double ran3(long *idum);
/*
    ** returns a uniform random number deviate between 0.0 and 1.0. Set
    ** the idum to any negative value to initialize or reinitialize the
    ** sequence. Any large MBIG, and any small (but still large) MSEED
    ** can be substituted for the present values. 
    */


     /***********   end function declaration ******/

      /*
      ** The function
      **        vectorScalarProd()
      ** calculates the product 
      **   scalarValue = Array<T,1> A * Array<T,1> B
      ** using Blitz array package
      */
template <typename T>
T inline scalarVectorProd(Array<T,1>& A,  Array<T,1>& B)
{
  firstIndex   i;

  return  sum(A(i) *  B(i),i); 
} // End: function  scalarVectorProd()

      /*
      ** The function
      **        matrixVectorProd()
      ** calculates the product 
      **   Array<T,1> C = Array<T,2> A * Array<T,1> B
      ** using Blitz array package
      */
template <typename T>
void  inline matrixVectorProd(Array<T,2>& A,  Array<T,1>& B, Array<T,1> C)
{
  firstIndex   i;
  secondIndex  j;

  C = sum(A(i,j) *  B(j),j); 
} // End: function  matrixVectorProd()

      /*
      ** The function
      **        matrixMatrixProd()
      ** calculates the product 
      **   Array<T,2> C = Array<T,2> A * Array<T,2> B
      ** using Blitz array package
      */
template <typename T>
void  inline matrixMatrixProd(Array<T,2>& A,  Array<T,2>& B, Array<T,2> C)
{
  firstIndex   i;
  secondIndex  j;
  thirdIndex   k;

  C = sum(A(i,k) *  B(k,j),k); 
} // End: function  matrixMatrixProd()


     /*
     ** The function 
     **       psdes()
     ** performs  a "Pseudo_DES" hashing  of a 64-bit word (lword, irword).
     ** both 32-bit arguments are returned hashed on all bits.
     */

void psdes(unsigned long& lword, unsigned long& irword)
{
  const int    NITER = 4;
  static const unsigned long c1[NITER] =
                  {0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
  static const unsigned long c2[NITER]=
                  {0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};
  unsigned long   i,ia,ib,iswap,itmph=0,itmpl=0;

  for(i = 0; i < NITER;i++) {
    ia = (iswap=irword) ^ c1[i];
    itmpl = ia & 0xffff;
    itmph = ia >> 16;
    ib     = itmpl*itmpl + ~(itmph*itmph);
    irword = lword ^ (((ia = (ib >> 16) |
		      ((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
    lword = iswap;
  }
}// End: function psdes()

    /*
    ** The template function
    **       ludcmp()
    ** takes as input a BLITZ matrix Array<T,2> A() of dimension n and
    ** replaces it by the LU decomposition of a rowwise permutation of
    ** itself. The results is stored in A() in the form given by 
    ** eq. (2.3.14) in "Numerical Recipe", sect. 2.3, page 45. The vector
    ** Array<T,1> Index() records the row permutation effected by the
    ** partial pivoting;
    ** The parameter d is output as +1 or -1 depending on whether the 
    ** number of row interchanges was even or odd, respectively. This
    ** routine is used in combination with the template function lubksb()
    ** to solve linear equations or invert a matrix.
    ** The function is slightly modified from the version in Numerical
    ** recipe
    */

template <typename T>
void ludcmp(Array<T,2>& A, int n, Array<int,1>& Index, T& d)
{
   int         i, imax, j, k;
   T           big, dum, sum, temp;
   Array<T,1>  vv(n);

   d = 1.0;                              // no row interchange yet
   for(i = 0; i < n; i++) {     // loop over rows to get scaling information
      big = ZERO;
      for(j = 0; j < n; j++) {
         if((temp = fabs(A(i,j))) > big) big = temp;
      }
      if(big == ZERO) {
	cout << endl << endl
             << "Singular Array<T,2> A() in template function ludcmp()"
             << endl <<endl;
         exit(1);
      }               
      vv(i) = 1.0/big;                 // save scaling */
   } // end i-loop */

   for(j = 0; j < n; j++) {     // loop over columns of Crout's method
      for(i = 0; i < j; i++) {   // not i = j
         sum = A(i,j);    
	 for(k = 0; k < i; k++) sum -= A(i,k) * A(k,j);
	 A(i,j) = sum;
      }
      big = ZERO;   // initialization for search for largest pivot element
      for(i = j; i < n; i++) {
         sum = A(i,j);
	 for(k = 0; k < j; k++)  {
	   sum -=  A(i,k) * A(k,j);
	 }
	 A(i,j) = sum;
	 if((dum = vv(i)*fabs(sum)) >= big) {
  	    big = dum;
	    imax = i;
	 }
      } // end i-loop
      if(j != imax) {    // do we need to interchange rows ?
         for(k = 0; k< n; k++) {       // yes
	    dum        = A(imax,k);
	     A(imax,k) = A(j,k);
	     A(j,k)    = dum;
	 }
	 d *= -1;            // and change the parit of d
	 vv(imax) = vv(j);         // also interchange scaling factor 
      }
      Index(j) = imax;
      if(fabs(A(j,j)) < ZERO)  A(j,j) = ZERO;

        /*
        ** if the pivot element is zero the matrix is singular
        ** (at least to the precision of the algorithm). For 
        ** some application of singular matrices, it is desirable
        ** to substitute ZERO for zero,
        */

      if(j < (n - 1)) {                   // divide by pivot element 
         dum = 1.0/A(j,j);
	 for(i = j+1; i < n; i++) A(i,j) *= dum;
      }
   } // end j-loop over columns
  
   vv.free();   // release local memory

}  // End: template function ludcmp()
 
    /*
    ** The template function 
    **             lubksb()
    ** solves the set of linear equations 
    **           A X = B 
    ** of dimension n. Array<T,2> A() is input, not as the
    ** matrix A() but rather as its LU decomposition, 
    ** determined by the template function ludcmp(),
    ** Array<int,1> Index() is input as the permutation vector
    ** returned by ludcmp(). Array<T,1> B() is input as the
    **  right-hand side vector B,
    ** The solution X is returned in B(). The input data A(),
    ** n and Index() are not modified. This routine take into 
    ** account the possibility that B() will begin with many
    ** zero elements, so it is efficient for use in matrix
    ** inversion.
    ** The function is slightly modified from the version in 
    ** in Numerical recipe.
    */

template <typename T>
void lubksb(Array<T,2>& A, int n, Array<int,1>& Index, Array<T,1>&  B)
{
   int        i, ii = 0, ip, j;
   T          sum;

   for(i = 0; i < n; i++) {
      ip    = Index(i);
      sum   = B(ip);
      B(ip) = B(i);
      if(ii != 0) {
	for(j = ii - 1; j < i; j++) sum -= A(i,j) * B(j);
      }
      else if(sum != 0.0 )  ii = i + 1;
      B(i) = sum;
   }
   for(i = n - 1; i >= 0; i--) {
      sum = B(i);
      for(j = i+1; j < n; j++) sum -= A(i,j) * B(j);
      B(i) = sum/A(i,i);
   }
} // end: template function lubksb()

      /*
      ** The function 
      **        chsone()
      ** takes two Blitz arrays, Array<double,1> bins() containing the observed 
      ** number of events and Array<double,1> ebins() containing expected number
      ** of events. Given the number of constrains knstrn (normally one). this 
      ** routine returns (trivally) the number of degrees of freedom df, and 
      ** (non-trivially) the chi-square chsq and the significance prob. A small 
      ** value of prob indicates a significant difference the distributions 
      ** bins() and ebins(). Note that bins() and ebins() are both arrays, 
      ** although bins will normally contain integer values.
      */   

void  chsone(Array<double,1>& bins, Array<double,1>& ebins, const int knstrn, 
                                        double& df, double& chsq, double& prob)
{
 int j;
 double temp;

 int nbins = bins.size();
 df        = nbins - knstrn;
 chsq      = 0.0;

 for(j = 0;j < nbins; j++) {
   if(ebins(j) <= 0.0) {
     cout << endl << "Bad expected number in chsone" << endl;
     exit(1);
   }
   temp = bins(j) - ebins(j);
   chsq += temp * temp/ebins(j);
 }
 prob = gammq(0.5 * df, 0.5 * chsq);

} // End: function  chsone()

      /*
      ** The function 
      **        chstwo()
      ** takes two Blitz arrays,  Array<double,1> bins1() and Array<double,1> bins2()
      **  containing two sets of binned data. Given the number of constrains 
      ** knstrn (normally one or two). this routine returns the number of degrees of
      ** freedom df, the chi-square chsq and the significance prob. A small 
      ** value of prob indicates a significant difference between the distributions 
      ** bins1() and bins2(). Note that bins1() and bins2() are both arrays, 
      ** although they will normally contain integer values.
      */  

void  chstwo(Array<double,1>& bins1, Array<double,1>& bins2, const int knstrn, 
                                        double& df, double& chsq, double& prob)
{
  int      j;
  double   temp;

  int nbins = bins1.size();   // initialization
  df        = nbins - knstrn;
  chsq      = 0.0;
  for(j = 0; j < nbins; j++)
    if(bins1(j) == 0.0 && bins2(j) == 0.0) --df;
    else {
      temp  = bins1(j) - bins2(j);
      chsq += temp * temp/(bins1(j) + bins2(j));
    }
  prob = gammq(0.5 * df, 0.5 * chsq);
} // End: function chstwo()

    /*
    ** The function
    **       gammq()
    ** returns the incomplete gamma function 
    **       Q/a,x) = 1 - P/a,x)
    */

double gammq(const double a, const double x)
{
  double    gamser, gammcf, gln;

  if(x < 0.0 || a <= 0.0){
    cout << "Invalid arguments in routine gammq";
    exit(1);
  }
  if(x < a + 1.0) {
    gser(gamser,a,x,gln);
    return (1.0 - gamser);
  } else {
    gcf(gammcf,a,x,gln);
    return gammcf;
  }
} // End: function gammq()

    /*
    ** The function
    **       gser()
    ** returns the incomplete gamma function P(a,x) 
    ** evaluated by its series representaion as
    ** gamser 
    ** It also returns ln(gamma(a)) as gln
    */

void gser(double& gamser, const double a, const double x, double &gln)
{
  const  int    ITMAX = 100;
  const  double EPS = numeric_limits<double>::epsilon();
  int    n;
  double sum,del,ap;

  gln = gammln(a);
  if(x <= 0.0) {
    if(x < 0.0) {
      cout <<endl <<"x less than 0 in function gser" << endl;
      exit(1);
    }
    gamser = 0.0;
    return;
  } else {
    ap = a;
    del = sum = 1.0/a;
    for(n = 0; n < ITMAX; n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if(fabs(del) < fabs(sum) * EPS) {
	gamser = sum*exp(-x+a * log(x)-gln);
	return;
      }
    }
    cout << endl 
         <<"a too large, ITMAX too small in routine gser"
         << endl;
    return;
  }
} // End: function gser()


    /*
    ** The function
    **       gammln()
    ** return the value ln[gamma(xx)]
    */

double gammln(const double xx)
{
  int j;
  double  x,y,tmp,ser;
  static const 
  double cof[6] = {76.18009172947146,     -86.50532032941677,
		   24.01409824083091,      -1.231739572450155,
                    0.1208650973866179E-2, -0.5395239384953E-5};

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for(j = 0;j < 6; j++)  ser += cof[j]/++y;
  
  return (-tmp+log(2.5066282746310005*ser/x));

} //End: template function gammln() 

    /*
    ** The function
    **       gcf()
    ** returns the incomplete gamma function Q(a,x)
    ** evaluated by its continued fraction representation
    ** as gammacf.
    ** Also returns ln gamma(a)  as gln
    */

void gcf(double& gammcf, const double a, const double x, double &gln)
{
  const int     ITMAX = 100;
  const double  EPS   = numeric_limits<double>::epsilon();
  const double  FPMIN = numeric_limits<double>::min()/EPS;
  int       i;
  double    an,b,c,d,del,h;

  gln = gammln(a);
  b   = x + 1.0 - a;
  c   = 1.0/FPMIN;
  d   = 1.0/b;
  h   = d;
  for(i = 1; i <= ITMAX; i++) {
    an  = -i * (i - a);
    b  += 2.0;
    d   = an * d + b;
    if(fabs(d) < FPMIN) d = FPMIN;
    c = b + an/c;
    if(fabs(c) < FPMIN) c = FPMIN;
    d = 1.0/d;
    del = d * c;
    h *= del;
    if(fabs(del - 1.0) <= EPS) break;
  }
  if(i > ITMAX)  {
    cout << "a too large, ITMAX too small in function gcf()";
    exit(1);
  }
  gammcf = exp(-x + a * log(x) - gln)*h;

} // End: templatefunction  gcf()

       /*
       ** The function 
       **              gauleg()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

void gauleg(double x1, double x2, Array<double,1>& x, Array<double,1>& w, int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359; 

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
	   ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
	 p2 =0.0;

   	   /*
	   ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

	 for(j = 1; j <= n; j++) {
	    p3 = p2;
	    p2 = p1;
	    p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
	 }

	   /*
	   ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */
 
	 pp = n * (z * p1 - p2)/(z * z - 1.0);
	 z1 = z;
	 z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > D_ZERO);

          /* 
	  ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      x(i - 1)  = xm - xl * z;
      x(n - i)  = xm + xl * z;
      w(i - 1)  = 2.0 * xl/((1.0 - z * z) * pp * pp);
      w(n - i)  = w(i - 1);
   } // end for(i) loop

} // End_ function gauleg()

template<typename T>
inline void rot(Array<T,2>& a, const T s, const T tau, const int i,
		const int j, const int k, const int l)       
{
  T     g, h;

  g = a(i,j);
  h = a(k,l);
  a(i,j) = g - s * (h + g * tau);
  a(k,l) = h + s * (g - h * tau);
} // End: template function rot()

template<typename T>
void jacobi(Array<T,2>& a, Array<T,1>& d, Array<T,2>& v, int &nrot)
{
  int 
           i,j,ip,iq;
  T 
           tresh,theta,tau,t,sm,s,h,g,c;
  
  int n = d.size();
  Array<T,1> b(n), z(n);
  for(ip = 0; ip < n; ip++) {
    for(iq = 0; iq < n; iq++) v(ip,iq) = 0.0;
    v(ip,ip) = 1.0;
  }
  for(ip = 0; ip < n; ip++) {
    b(ip) = d(ip) = a(ip,ip);
    z(ip) = 0.0;
  }
  nrot=0;
  for(i = 1; i <= 50; i++) {
    sm = 0.0;
    for( ip = 0; ip < n-1; ip++) {
      for(iq = ip + 1; iq < n; iq++)
	sm += fabs(a(ip,iq));
    }
    if(sm == 0.0)    return;
    if(i < 4)
      tresh = 0.2 * sm/(n * n);
    else
      tresh=0.0;
    for(ip = 0; ip < n-1; ip++) {
      for(iq = ip+1; iq < n; iq++) {
	g = 100.0 * fabs(a(ip,iq));
	if(i > 4 && (fabs(d(ip)) + g) == fabs(d(ip))
	    && (fabs(d(iq)) + g) == fabs(d(iq)))
	  a(ip,iq) = 0.0;
	else if(fabs(a(ip,iq)) > tresh) {
	  h = d(iq) - d(ip);
	  if((fabs(h)+g) == fabs(h))
	    t = (a(ip,iq))/h;
	  else {
	    theta = 0.5 * h/(a(ip,iq));
	    t = 1.0/(fabs(theta) + sqrt(1.0 + theta * theta));
	    if(theta < 0.0) t = -t;
	  }
	  c = 1.0/sqrt(1 + t * t);
	  s = t * c;
	  tau = s/(1.0 + c);
	  h =t * a(ip,iq);
	  z(ip) -= h;
	  z(iq) += h;
	  d(ip) -= h;
	  d(iq) += h;
	  a(ip,iq) = 0.0;
	  for(j = 0; j < ip; j++)
	    rot(a, s, tau, j, ip, j, iq);
	  for(j = ip+1; j < iq;j++)
	    rot(a, s, tau, ip, j, j, iq);
	  for(j = iq +1; j < n; j++)
	    rot(a, s, tau, ip, j, iq, j);
	  for(j = 0; j < n; j++)
	    rot(v, s, tau, j, ip, j, iq);
	  ++nrot;
	}
      }
    }
    for(ip = 0; ip < n; ip++) {
      b(ip) += z(ip);
      d(ip) = b(ip);
      z(ip) = 0.0;
    }
  }
  cout << endl << "Too many iterations in routine jacobi";
  exit(1);
} // End: template function jacobi()

template<typename T>
void jacobn_s(const T x, Array<T,1>& y, Array<T,1>& dfdx, Array<T,2>& dfdy)
{
  int i;

  int n = y.size();
  for(i = 0; i < n; i++) dfdx(i) = 0.0;
  dfdy(0,0) = -0.013 - 1000.0 * y(2);
  dfdy(0,1) = 0.0;
  dfdy(0,2) = -1000.0 * y(0);
  dfdy(1,0) = 0.0;
  dfdy(1,1) = -2500.0 * y(2);
  dfdy(1,2) = -2500.0 * y(1);
  dfdy(2,0) = -0.013 - 1000.0 * y(2);
  dfdy(2,1) = -2500.0 * y(2);
  dfdy(2,2) = -1000.0 * y(0) - 2500.0 * y(1);

} // End: template function jacobn_s


template<typename T>
void derivs_s(const T x, Array<T,1>& y, Array<T,1>& dydx)
{
  dydx(0) = -0.013 * y(0) - 1000.0 * y(0) * y(2);
  dydx(1) = -2500.0 * y(1) * y(2);
  dydx(2) = -0.013 * y(0) - 1000.0 * y(0) * y(2) - 2500.0 * y(1) * y(2);

} // End: template function derivs_s()

template<typename T>
void eigsrt(Array<T,1>& d, Array<T,2>& v)
{
  int 
      i,j,k;
  T 
      p;

  int n = d.size();
  for(i = 0; i <n -1; i++) {
    p=d(k = i);
    for(j = i; j < n; j++)
      if(d(j) >= p) p = d(k = j);
    if(k != i) {
      d(k) = d(i);
      d(i) = p;
      for(j = 0;j < n; j++) {
	p = v(j,i);
	v(j,i) = v(j,k);
	v(j,k) = p;
      }
    }
  }
}// End: template function eigsrt()

    /*
    ** The template function 
    **            tqli()
    ** determine eigenvalues and eigenvectors of a real symmetric
    ** tri-diagonal matrix, or a real, symmetric matrix previously
    ** reduced by function tred2() to tri-diagonal form. On input,
    ** Array<T,1> d() contains the diagonal element and 
    ** Array<T,1>e() the sub-diagonal of the tri-diagonal matrix.
    ** On output d[] contains the eigenvalues and  Array<T,1> e()
    ** is destroyed. If eigenvectors are desired Array<T,2>z( , )
    ** on input contains the identity matrix. If eigenvectors of a 
    ** matrix reduced by tred2() are required, then Array<T,2>z( , )
    ** on input is the matrix output from tred2().
    ** On output, the k'th column returns the normalized eigenvector
    ** corresponding to Array<T,1> d(k). 
    ** The function is modified from the version in Numerical recipe.
    */

template<typename T> 
void tqli(Array<T,1>& d, Array<T,1>& e, Array<T,2>& z)
{
  int 
        m, l, iter, i, k;
  T
         s, r, p, g, f, dd, c, b;

  int n = d.size();
  for(i = 1; i < n; i++) e(i - 1) = e(i);
  e(n - 1) = 0.0;
  for(l = 0; l < n; l++) {
    iter = 0;
    do {
      for(m = l; m < n-1; m++) {
	dd = fabs(d(m)) + fabs(d(m + 1));
	if(fabs(e(m)) + dd == dd) break;
      }
      if(m != l) {
	if(iter++ == 30) {
	  cout << endl <<"Too many iterations in tqli";
	  exit(1);
	}
	g = (d(l + 1) - d(l))/(2.0 * e(l));
	r = pythag(g, 1.0);
	g = d(m) - d(l) + e(l)/(g + SIGN(r,g));
	s = c = 1.0;
	p = 0.0;
	for(i = m - 1;i >= l; i--) {
	  f = s * e(i);
	  b = c * e(i);
	  e(i + 1) = (r = pythag(f,g));
	  if(r == 0.0) {
	    d(i + 1) -= p;
	    e(m) = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d(i + 1) - p;
	  r = (d(i) - g) * s + 2.0 * c * b;
	  d(i + 1) = g + (p = s * r);
	  g = c * r - b;

	  // Next loop can be omitted if eigenvectors not wanted

	  for(k = 0; k < n; k++) {
	    f = z(k,i + 1);
	    z(k,i + 1) = s * z(k,i) + c * f;
	    z(k,i) = c * z(k,i) - s * f;
	  }
	}
	if(r == 0.0 && i >= l) continue;
	d(l) -= p;
	e(l) = g;
	e(m) = 0.0;
      }
    } while(m != l);
  }
} // End: template function  tqli()


template<typename T> 
void tqli(Array<T,1>& d, Array<T,1>& e)
{
  int 
        m, l, iter, i, k;
  T
         s, r, p, g, f, dd, c, b;

  int n = d.size();
  for(i = 1; i < n; i++) e(i - 1) = e(i);
  e(n - 1) = 0.0;
  for(l = 0; l < n; l++) {
    iter = 0;
    do {
      for(m = l; m < n-1; m++) {
	dd = fabs(d(m)) + fabs(d(m + 1));
	if(fabs(e(m)) + dd == dd) break;
      }
      if(m != l) {
	if(iter++ == 30) {
	  cout << endl <<"Too many iterations in tqli";
	  exit(1);
	}
	g = (d(l + 1) - d(l))/(2.0 * e(l));
	r = pythag(g, 1.0);
	g = d(m) - d(l) + e(l)/(g + SIGN(r,g));
	s = c = 1.0;
	p = 0.0;
	for(i = m - 1;i >= l; i--) {
	  f = s * e(i);
	  b = c * e(i);
	  e(i + 1) = (r = pythag(f,g));
	  if(r == 0.0) {
	    d(i + 1) -= p;
	    e(m) = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d(i + 1) - p;
	  r = (d(i) - g) * s + 2.0 * c * b;
	  d(i + 1) = g + (p = s * r);
	  g = c * r - b;
	}
	if(r == 0.0 && i >= l) continue;
	d(l) -= p;
	e(l)=g;
	e(m) =0.0;
      }
    } while(m != l);
  }
} // End: template function  tqli()

    /*
    ** The function
    **                tred2()
    ** perform a Housholder reduction of a real symmetric matrix
    ** Array<T,2> a(). On output Array<T,2> a() is replaced by 
    ** the orthogonal matrix effecting the transformation. 
    ** Array<T,1> d() returns the diagonal elements of the 
    ** tri-diagonal matrix, and Array<T,1> e() the off-diagonal 
    ** elements,  with Array<T,1> e() = 0.
    ** The function is modified from the version in Numerical recipe.
    */

template<typename T>
void tred2(Array<T,2>& a, Array<T,1>& d, Array<T,1>& e)
{
  int
         l,k,j,i;
  T
         scale, hh, h, g, f;

  int n = d.size();
  for(i = n-1; i > 0; i--) {
    l = i - 1;
    h = scale = 0.0;
    if(l > 0) {
      for(k = 0; k < l + 1; k++)
	scale += fabs(a(i,k));
      if(scale == 0.0)
	e(i) = a(i,l);
      else {
	for(k = 0; k < l + 1; k++) {
	  a(i,k) /= scale;
	  h += a(i,k) * a(i,k);
	}
	f = a(i,l);
	g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
	e(i) = scale * g;
	h -= f * g;
	a(i,l) = f - g;
	f = 0.0;
	for(j = 0; j < l+1;j++) {

	     // Next statement can be omitted if eigenvectors not wanted

	  a(j,i) = a(i,j)/h;
	  g = 0.0;
	  for(k = 0; k < j + 1; k++)
	    g += a(j,k) * a(i,k);
	  for(k = j + 1; k < l + 1; k++)
	    g += a(k,j) * a(i,k);
	  e(j) = g/h;
	  f   += e(j) * a(i,j);
	}
	hh = f/(h + h);
	for(j = 0; j < l + 1; j++) {
	  f = a(i,j);
	  e(j) = g = e(j) - hh * f;
	  for(k = 0; k < j + 1; k++)
	    a(j,k) -= (f * e(k) + g * a(i,k));
	}
      }
    } else
      e(i) = a(i,l);
    d(i) = h;
  }

        // Next statement can be omitted if eigenvectors not wanted

  d(0) = 0.0;
  e(0) = 0.0;

        // Contents of this loop can be omitted if eigenvectors not
        //	wanted except for statement d[i]=a[i][i];

  for(i = 0; i < n; i++) {
    l = i;
    if(d(i) != 0.0) {
      for(j = 0; j < l; j++) {
	g = 0.0;
	for(k = 0; k < l; k++)
	  g += a(i,k) * a(k,j);
	for(k = 0; k < l; k++)
	  a(k,j) -= g * a(k,i);
      }
    }
    d(i) = a(i,i);
    a(i,i) = 1.0;
    for(j = 0; j < l; j++) a(j,i) = a(i,j) = 0.0;
  }
} // End: template function  tred2()

template<typename T>
 T  pythag(T a, T b)
{
  T    absa, absb;

  absa = fabs(a);
  absb = fabs(b);
  if(absa > absb) return absa * sqrt(1.0 + SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa/absb)));
}
// End: function pythag(),

      /*
      ** The function
      **              rk4()
      ** takes a set of variables Array<double,1> y(0 : n-1) for the function y(x)
      ** together with the derivatives  Array<double,1> dydx(0:n-1) and uses the
      ** fourth-order Runge-Kutta method to advance the solution over an interval
      ** h and return incremented variables as Array<double,1> yout(0:n-1), which
      ** not need to be a distinct array from y[1:n]. The users supply the routine
      ** derivs(double x,Array<double,1> y, Array<double,1> dydx), which returns 
      ** the derivatives dydx at x.
      */ 

void rk4(Array<double,2>& y, Array<double,2>& dydx, int n, double x, double h, 
                 Array<double,2>& yout,
	 void (*derivs)(double, Array<double,2>&, Array<double,2>&))
{
  int      i;
  double   xh, hh, h6;
  Array<double,2> dym(n),dyt(n), yt(n);

   hh = h * 0.5;
   h6 = h/6.0;
   xh = x+hh;

   for(i = 0; i < n; i++) {                 // first step
      yt(i) = y(i) + hh * dydx(i);
   }

   (*derivs)(xh, yt, dyt);                 // second step

   for(i = 1; i < n; i++) {
      yt(i) = y(i) + hh * dyt(i);
   }

   (*derivs)(xh, yt, dym);                // third step 

   for(i = 1; i < n; i++) {
      yt(i)   = y(i) + h * dym(i);
      dym(i) += dyt(i);
   }
	
   (*derivs)(x+h,yt,dyt);                // fourth step 

        // acummulate increments with proper weights

   for(i = 0; i< n; i++) {
      yout(i) = y(i) + h6 *(dydx(i) + dyt(i) + 2.0 * dym(i));
   }
   
} // End: function rk4()

         /*
	 ** The function 
         **           spline()
         ** takes as input x(0,..,n - 1) and y(0,..,n - 1) containing a tabulation
         ** y_i = f(x_i) with x_0 < x_1 < .. < x_(n - 1) together with yp_1 and yp2
         ** for first derivatives  f(x) at x_0 and x_(n-1), respectively. Then the
         ** function returns y2(0,..,n-1) which contanin the second derivatives of
         ** f(x_i)at each point x_i. If yp1 and/or yp2 is larger than the constant
         ** INFINITY the function will put corresponding second derivatives to zero.
         */ 

void spline(Array<double,1>& x, Array<double,1>& y, int n, double& yp1, 
          double yp2, Array<double,1> y2)
{ 
   int             i,k;
   double          p,qn,sig,un;
   Array<double,1> u;

   if(yp1 > INFINITY)  y2(0) = u(0) = 0.0;
   else {
      y2(0) = -0.5;
      u(0)  = (3.0/(x(1) - x(0))) * ((y(1) - y(0))/(x(1) - x(0)) - yp1);
   }
   for(i = 1; i < (n - 1); i++) {
      sig   = (x(i) - x(i - 1))/(x(i + 1) - x(i - 1));
      p     = sig * y2(i - 1) + 2.0;
      y2(i) = (sig - 1.0)/p;
      u(i)  = (y(i + 1) - y(i))/(x(i + 1) - x(i)) - (y(i) - y(i - 1))/(x(i) - x(i - 1));
      u(i)  = (6.0 * u(i)/(x(i + 1) - x(i - 1)) - sig*u(i - 1))/p;
   }
   if(yp2 > INFINITY)  qn = un = ZERO;
   else {
      qn = 0.5;
      un = (3.0/(x(n - 1) - x(n - 2))) * (yp2 - (y(n - 1) - y(n - 2))/(x(n - 1) - x(n - 2)));
   }
   y2(n - 1) = (un - qn * u(n - 2))/(qn * y2(n - 2) + 1.0);

   for(k = n - 2; k >= 0; k--) {
      y2(k) = y2(k)*y2(k+1)+u(k);
   }

}  // End: function spline()

     /*
     ** The function 
     **           splint()
     ** takes xa[0,..,n - 1] and y[0,..,n - 1] which tabulates a function 
     ** (with the xa[i]'s in order) and given ya[0,..,n - 1], which is the
     ** output from function spline() and with given value of x returns a 
     ** cubic--spline interpolation value y.
     */

void splint(Array<double,1>& xa, Array<double,1>& ya, Array<double,1>& y2a, 
             int n, double x, double& y)
{
   int       klo,khi,k;
   double    h,b,a;

   klo = 0;
   khi = n - 1;
   while((khi - klo) > 1) {   // binary search
      k = (khi + klo) >> 1;
      if(xa(k) > x)   khi = k;
      else            klo = k;
   }
   h = xa(khi) - xa(klo);
   if(fabs(h) < D_ZERO) {
      printf("\n\n Error in function splint(): ");
      printf("\n The difference h = %4.1E -- too small\n",h);
      exit(1);
   }
   a  = (xa(khi) - x)/h;
   b  = (x - xa(klo))/h;
   y  =   a * ya(klo) + b * ya(khi) + ((a * a * a - a) * y2a(klo) 
       + (b * b * b - b) * y2a(khi)) * (h * h)/6.0;

} // End: function splint()



     /*
     ** The function
     **           ran0()
     ** is an "Minimal" random number generator of Park and Miller
     ** (see Numerical recipe page 279). Set or reset the input value
     ** idum to any integer value (except the unlikely value MASK)
     ** to initialize the sequence; idum must not be altered between
     ** calls for sucessive deviates in a sequence.
     ** The function returns a uniform deviate between 0.0 and 1.0.
     */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double ran0(long *idum)
{
   long     k;
   double   ans;

   *idum ^= MASK;
   k = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   ans=AM*(*idum);
   *idum ^= MASK;
   return ans;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef MASK

// End: function ran0() 

     /*
     ** The function
     **           ran1()
     ** is an "Minimal" random number generator of Park and Miller
     ** (see Numerical recipe page 280) with Bays-Durham shuffle and
     ** added safeguards. Call with idum a negative integer to initialize;
     ** thereafter, do not alter idum between sucessive deviates in a
     ** sequence. RNMX should approximate the largest floating point value
     ** that is less than 1.
     ** The function returns a uniform deviate between 0.0 and 1.0
     ** (exclusive of end-point values).
     */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
   int             j;
   long            k;
   static long     iy=0;
   static long     iv[NTAB];
   double          temp;

   if (*idum <= 0 || !iy) {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (*idum)/IQ;
         *idum = IA*(*idum - k*IQ) - IR*k;
         if(*idum < 0) *idum += IM;
         if(j < NTAB) iv[j] = *idum;
      }
      iy = iv[0];
   }
   k     = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   j     = iy/NDIV;
   iy    = iv[j];
   iv[j] = *idum;
   if((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran1()

     /*
     ** The function 
     **         ran2()
     ** is a long periode (> 2 x 10^18) random number generator of 
     ** L'Ecuyer and Bays-Durham shuffle and added safeguards.
     ** Call with idum a negative integer to initialize; thereafter,
     ** do not alter idum between sucessive deviates in a
     ** sequence. RNMX should approximate the largest floating point value
     ** that is less than 1.
     ** The function returns a uniform deviate between 0.0 and 1.0
     ** (exclusive of end-point values).
     */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(long *idum)
{
   int            j;
   long           k;
   static long    idum2 = 123456789;
   static long    iy=0;
   static long    iv[NTAB];
   double         temp;

   if(*idum <= 0) {
      if(-(*idum) < 1) *idum = 1;
      else             *idum = -(*idum);
      idum2 = (*idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (*idum)/IQ1;
	 *idum = IA1*(*idum - k*IQ1) - k*IR1;
	 if(*idum < 0) *idum +=  IM1;
	 if(j < NTAB)  iv[j]  = *idum;
      }
      iy=iv[0];
   }
   k     = (*idum)/IQ1;
   *idum = IA1*(*idum - k*IQ1) - k*IR1;
   if(*idum < 0) *idum += IM1;
   k     = idum2/IQ2;
   idum2 = IA2*(idum2 - k*IQ2) - k*IR2;
   if(idum2 < 0) idum2 += IM2;
   j     = iy/NDIV;
   iy    = iv[j] - idum2;
   iv[j] = *idum;
   if(iy < 1) iy += IMM1;
   if((temp = AM*iy) > RNMX) return RNMX;
   else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran2()

    /*
    ** The function
    **        ran3()
    ** returns a uniform random number deviate between 0.0 and 1.0. Set
    ** the idum to any negative value to initialize or reinitialize the
    ** sequence. Any large MBIG, and any small (but still large) MSEED
    ** can be substituted for the present values. 
    */

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3(long *idum)
{
   static int        inext, inextp;
   static long       ma[56];      // value 56 is special, do not modify
   static int        iff = 0;
   long              mj, mk;
   int               i, ii, k;

   if(*idum < 0 || iff == 0) {                 // initialization
      iff    = 1;

      mj     = MSEED - (*idum < 0 ? -*idum : *idum);
      mj    %= MBIG;
      ma[55] = mj;                            // initialize ma[55] 

      for(i = 1, mk = 1; i <= 54; i++) {      // initialize rest of table 
         ii     = (21*i) % 55;
	 ma[ii] = mk;
	 mk     = mj - mk;
	 if(mk < MZ) mk += MBIG;
	 mj = ma[ii];
      }

      for(k = 1; k <= 4; k++) {   // randimize by "warming up" the generator
         for(i = 1; i <= 55; i++) {
	    ma[i] -= ma[1 + (i + 30) % 55];
	    if(ma[i] < MZ) ma[i] += MBIG;
	 }
      }

      inext  =  0;              // prepare indices for first generator number
      inextp = 31;              // 31 is special
      *idum  = 1;
   }

   if(++inext == 56)  inext  = 1;
   if(++inextp == 56) inextp = 1;
   mj = ma[inext] - ma[inextp];
   if(mj < MZ) mj += MBIG;
   ma[inext] = mj;
   return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

// End: function ran3()


#endif
