#include "lgamma.h"

std::complex<double> lgamma(std::complex<double> z)
{
  double re = std::real(z);
  double im = std::imag(z);
  
  gsl_sf_result r1;
  gsl_sf_result r2;
  
  gsl_sf_lngamma_complex_e (re, im, &r1, &r2);
  std::complex<double> result(r1.val, r2.val);
  
  return result;
}
