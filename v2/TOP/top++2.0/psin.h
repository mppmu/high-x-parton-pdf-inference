/*******************************************************************************
*									       *	
*	     The Polygamma function Psi(z,k), for k=0,1,2,3,4 ,                *
*	     for Complex argument z.			                       *
*		    - Psi(z) == Psi(z,0).			               *
*	     A.K.A. Digamma function (for k=0), Trigamma (for k=1) etc.        *
*									       *	
*******************************************************************************/


#ifndef guard_psin_h
#define guard_psin_h


#include <complex>
#include <iostream>

inline int nint(double x) //returns the (same-sign) integer closest to x.
{
  return x>=0 
    ? 
    ((x - int(x) <= 0.5 ) ? int(x) : int(x) + 1) 
    : 
    ((int(x)-x <= 0.5 ) ? int(x) : int(x) - 1) ;
}

std::complex<double> psin(std::complex<double> z, 
			  int k);

inline std::complex<double> psi0(std::complex<double> z)
{
  return psin(z,0);
}

inline std::complex<double> psi1(std::complex<double> z)
{
  return psin(z,1);
}

#endif
