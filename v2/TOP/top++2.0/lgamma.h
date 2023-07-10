/*******************************************************************************
*      								               *
*	        The function Log[Gamma[z]] for Complex argument z.             *
*	       	A wrap for the corresponding function in the GSL library.      *
*									       *
*******************************************************************************/

#ifndef guard_lgamma_h
#define guard_lgamma_h

#include <gsl/gsl_sf_gamma.h>
#include <complex>

std::complex<double> lgamma(std::complex<double> z);

#endif
