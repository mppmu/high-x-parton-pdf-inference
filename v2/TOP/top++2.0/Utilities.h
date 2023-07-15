/*******************************************************************************	
*									       *	
*  - Functions that convolute the perturbative components with the fluxes      *
*  - Functions and functionalities for calculating PDF uncertainties           * 
*  - Other functions                                                           * 
*									       *	
*******************************************************************************/

#ifndef guard_Utilities_h
#define guard_Utilities_h

#include <stdio.h>
#include <gsl/gsl_statistics.h>
#include "PartonicFlux.h"
#include "SubtrFlux.h"
#include "FixedOrder.h"
#include "Resummation.h"

struct FOcomponents 
{
  PartonicFlux Flux; 
  FixedOrder FO; 
  double rho;
  std::string partonchannel;
};		

struct REScomponents 
{
  PartonicFlux Flux; 
  SubtrFlux SFlux; 
  Resummation Res; 
  double rho;
  std::string partonchannel;
};		

double ComputeFinalResultFO(PartonicFlux Flux, 
			    FixedOrder FO, 
			    double rho, 
			    int precision,
			    std::string partonchannel);

double ComputeFinalResultRES(PartonicFlux Flux, 
			     SubtrFlux SFlux, 
			     FixedOrder FO, 
			     Resummation Res, 
			     double rho, 
			     int precision,
			     std::string partonchannel);

double f_FO (double x, void * par);
double f_RES (double x, void * par);

double ConvolutionFO(PartonicFlux Flux, 
		     FixedOrder FO, 
		     double rho, 
		     int precision,
		     std::string partonchannel);

double ConvolutionRES(PartonicFlux Flux, 
		      SubtrFlux SFlux, 
		      Resummation Res, 
		      double rho, 
		      int precision,
		      std::string partonchannel);

// If Fact. and Ren. scales pass the scale variation criterion.
bool ScalesPass(std::string RestrictedScaleVariation, 
		double Fact, 
		double Ren);

// Converts OrderRES and OrderFO to integers:
int IntegerForm(std::string );

// Print a greeting message at startup
// (argument == program's version):
void Greetings(std::string );

//------------------------------ Calculation of PDF uncertainty
void deltaPDFasymmetric(const std::vector<double>& result, 
			double& smax, 
			double& smin);

void deltaPDFnnpdf(const std::vector<double>& result, 
		   double& mean, 
		   double& smax, 
		   double& smin);

void deltaPDFsymmetric(const std::vector<double>& result, 
		       double& smax, 
		       double& smin);

void deltaPDFheraVAR(const std::vector<double>& result, 
		     double& smax, 
		     double& smin);

#endif
