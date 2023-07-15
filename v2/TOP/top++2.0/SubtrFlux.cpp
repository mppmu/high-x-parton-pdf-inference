#include "SubtrFlux.h"

using std::vector;
using std::string;
using std::domain_error;
using std::size_t;

typedef std::complex<double>  dcomp;

/*******************************************************************************	
*									       *	
*                  The four public member functions                            * 
*           (return the subtraction flux in x- and N-space)    		       *
*									       *
*******************************************************************************/

double SubtrFlux::subtrfluxgg(double x)
{
  PointerToFluxMethod pf = &PartonicFlux::fluxgg;
  return CompSubtrFluxX(pf, x);
}

double SubtrFlux::subtrfluxqqbar(double x)
{
  PointerToFluxMethod pf = &PartonicFlux::fluxqqbar;
  return CompSubtrFluxX(pf, x);
}

dcomp SubtrFlux::subtrfluxggN(dcomp N)
{
  PointerToFluxMethod pf = &PartonicFlux::fluxgg;
  return CompSubtrFluxN(pf, N);
}

dcomp SubtrFlux::subtrfluxqqbarN(dcomp N)
{
  PointerToFluxMethod pf = &PartonicFlux::fluxqqbar;
  return CompSubtrFluxN(pf, N);
}

/*******************************************************************************	
*                                                                              *
*          General result for the subtraction flux in analytical form.         *
*          Implemented as in http://arXiv.org/abs/hep-ph/9604351v2	       *	
*									       *
*******************************************************************************/

double SubtrFlux::CompSubtrFluxX(PointerToFluxMethod PFUNC, 
				 double tau)
{
  vector<double> x;//point over which SubtrFlux=Flux
  vector<double> l;//value of Flux at the points x
  
  x.push_back(0.0);// The count should be from 1 to 4; So 
  l.push_back(0.0);// add a fake 0-th element.
  for (int j=1; j<=4; ++j)
    {
      x.push_back( rho + (j-1)*ETA );
      l.push_back( (Flux.*PFUNC)(x.back()) );
    }
  double al = (log(l[1]/l[3])*log(x[1]/x[2])
	       -log(l[1]/l[2])*log(x[1]/x[3]))
    /(log(x[1]/x[2])*log((1-x[1])/(1-x[3]))
      -log(x[1]/x[3])*log((1-x[1])/(1-x[2])));
  double be = (log(l[1]/l[3])*log((1-x[1])/(1-x[2]))
	       -log(l[1]/l[2])*log((1-x[1])/(1-x[3])))
    /(log(x[1]/x[2])*log((1-x[1])/(1-x[3]))
      -log(x[1]/x[3])*log((1-x[1])/(1-x[2])));
  double B  = l[1]/( pow(1-x[1],al) * pow(x[1],-be) );
  double C  = ( B * pow(1-x[4],al) 
		* pow(x[4],-be) - l[4] )
    /( B * pow(1-x[4],al) * pow(x[4],-be) 
       * (x[4]-x[1]) * (x[4]-x[2]) * (x[4]-x[3]) );
  double f = B * pow(1-tau,al) * pow(tau,-be) 
    * (1-C*(tau-x[1])*(tau-x[2])*(tau-x[3]));
  return f;
}

dcomp SubtrFlux::CompSubtrFluxN(PointerToFluxMethod PFUNC, 
				dcomp N)
{
  vector<double> x;
  vector<double> l;
  x.push_back(0.0);
  l.push_back(0.0);
  for (int j=1; j<=4; ++j)
    {
      x.push_back( rho + (j-1)*ETA );
      l.push_back( (Flux.*PFUNC)(x.back()) );
    }
  double al = (log(l[1]/l[3])*log(x[1]/x[2])
	       -log(l[1]/l[2])*log(x[1]/x[3]))
    /(log(x[1]/x[2])*log((1-x[1])/(1-x[3]))
      -log(x[1]/x[3])*log((1-x[1])/(1-x[2])));
  double opal = 1+al;
  double be = (log(l[1]/l[3])*log((1-x[1])/(1-x[2]))
	       -log(l[1]/l[2])*log((1-x[1])/(1-x[3])))
    /(log(x[1]/x[2])*log((1-x[1])/(1-x[3]))
      -log(x[1]/x[3])*log((1-x[1])/(1-x[2])));
  double B  = l[1]/( pow(1-x[1],al) * pow(x[1],-be) );
  double C  = ( B * pow(1-x[4],al) * pow(x[4],-be) - l[4] )
    /( B * pow(1-x[4],al) * pow(x[4],-be) 
       * (x[4]-x[1]) * (x[4]-x[2]) * (x[4]-x[3]) );
  
  dcomp f = B*(EulerBeta(N-be,opal) 
	       + C*( x[1]*x[2]*x[3]*EulerBeta(N-be,opal)
		     -(x[1]*x[2]+x[1]*x[3]+x[2]*x[3])*EulerBeta(N-be+1.0,opal)
		     +(x[1]+x[2]+x[3])*EulerBeta(N-be+2.0,opal)
		     -EulerBeta(N-be+3.0,opal)));
  return f;
}
