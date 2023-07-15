/*******************************************************************************	
*									       *	
*    A class that initializes the partonic fluxes for a given set of PDF's.    *
*    The initialization computes the fluxes on a grid of size 'nmax'.	       *
*    It provides a function-interface that allows to take the value	       *
*    of the flux for any real value of 'x'; Second order finite-difference     *
*    scheme implemented.						       *
*       - tmin = min. value for the convolution variable 'x'	               *
*       - tmax = max. value for the convolution variable 'x'; typically tmax=1 *
*       - muF = Factorization scale (in GeV)				       *
*       - nmax = number of points to split the interval (tmin,tmax).	       *
*									       *
*******************************************************************************/

#ifndef guard_class_PartonicFlux_h
#define guard_class_PartonicFlux_h

#include <math.h>
#include <gsl/gsl_integration.h>
#include "LHAPDF/LHAPDF.h"
#include <vector>
#include <stdexcept>
#include <cstddef>

class PartonicFlux 
{
 public:
 PartonicFlux(bool _ColliderName, 
	      size_t _PDFmember, 
	      double _tmin, 
	      double _tmax, 
	      double _muF, 
	      std::size_t _nmax, 
	      int _prec) :
  ColliderName(_ColliderName),
    PDFmember(_PDFmember),
    tmin(_tmin),
    tmax(_tmax),
    muF(_muF),
    nmax(_nmax),
    prec(_prec)
    { 
      create(); 
    }
  
  double fluxgg(double x);
  double fluxqqbar(double x);
  double fluxgq(double x);
  double fluxqq(double x);
  double fluxqqprime(double x);
  double fluxqqbarprime(double x);
  
 private:
  bool ColliderName;//1==LHC,0==TEV
  size_t PDFmember;
  double tmin, tmax, muF, logtmin, logtmax;
  std::size_t nmax;
  int prec;
  
  std::vector<double> vfgg,vfqqbar,vfgq,vfqq,vfqqprime,vfqqbarprime;
  
  void create();
  void initflux();
  double flux(int state, double tau);
  
  static double PDFproduct(int state, 
			   double x1, 
			   double x2, 
			   double scaleF);
  
  static double f_p (double x, void * par); 
  
  static double PDFConvolution(int state, 
			       double t, 
			       double scaleF,
			       int prec);
};

struct params { int state; double t; double scaleF; };		

#endif
