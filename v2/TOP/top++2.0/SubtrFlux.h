/*******************************************************************************	
*									       *	
*     A fake partonic subtraction flux.				               *
*     Implementation of http://arXiv.org/abs/hep-ph/9604351v2                  *
*     The idea is: SubtrFlux ~ Flux close to absolute treshold, for any Flux.  *
*       ETA ~ 10^-5 is a small step. It defines the distance from threshold    *
*        where SubtrFlux approximates Flux.				       *
*									       *
*******************************************************************************/

#ifndef guard_class_SubtrFlux_h
#define guard_class_SubtrFlux_h

#include <vector>
#include <stdexcept>
#include <cstddef>
#include <complex>

#include "PartonicFlux.h"
#include "lgamma.h"


class SubtrFlux 
{
 public:
  
 SubtrFlux(PartonicFlux iFlux, 
	   double irho, 
	   double iETA) :
  Flux(iFlux), 
    rho(irho), 
    ETA(iETA) 
    {}
  
  double subtrfluxgg(double x);
  double subtrfluxqqbar(double x);
  
  std::complex<double> subtrfluxggN(std::complex<double> N);
  std::complex<double> subtrfluxqqbarN(std::complex<double> N);
  
 private:
  PartonicFlux Flux;
  double rho, ETA;
  
  typedef double (PartonicFlux::*PointerToFluxMethod )(double );
  
  double CompSubtrFluxX (PointerToFluxMethod, 
			 double tau);
  std::complex<double> CompSubtrFluxN(PointerToFluxMethod, 
				      std::complex<double> N);
};

inline std::complex<double> EulerBeta(std::complex<double> x1, 
				      std::complex<double> x2)
{	
  std::complex<double> f = exp(lgamma(x1)+lgamma(x2)-lgamma(x1+x2));
  return f;
}

#endif
