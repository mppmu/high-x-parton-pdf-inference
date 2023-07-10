/*******************************************************************************	
*									       *	
*	  A class that evaluates the fixed order results through NNLO.	       *
*     It wraps a large number of functions in a simple interface.	       *
*	  Its initialization involves no calculations.			       *
*									       *	
*******************************************************************************/

#ifndef guard_class_FixedOrder_h
#define guard_class_FixedOrder_h

#include <math.h>
#include <vector>
#include <stdexcept>
#include <cstddef>
#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_sf_dilog.h>

class FixedOrder {
 public:
 FixedOrder(double _as,
	    double _mt, 
	    double _muF, 
	    double _muR, 
	    bool _LO, 
	    bool _NLO, 
	    bool _NNLO) :
  as(_as),
    mt(_mt),
    muF(_muF),
    muR(_muR), 
    LO(_LO),
    NLO(_NLO),
    NNLO(_NNLO)
    { 
      create(); 
    }
  
  double FOgg(double rho);
  double FOqqbar(double rho);
  double FOgq(double rho);
  double FOqq(double rho);
  double FOqqprime(double rho);
  double FOqqbarprime(double rho);
  
 private:
  static double pi;
  double as, mt, muF, muR, lg, lgfr;
  bool LO, NLO, NNLO;	
  
  void create();	
  
  double PXS_LO_gg(double rho);
  double PXS_NLO_gg(double rho);
  double PXS_NNLO_EXACT_gg(double rho);
  double PXS_NNLO_THRESHOLD_gg(double rho);
  double PXS_NNLO_FIT_gg(double rho);
  double PXS_NNLO_LO_NLOscales_gg(double rho);
  double PXS_NNLO_scale_terms_gg(double rho);
  
  double PXS_LO_qqbar(double rho);
  double PXS_NLO_qqbar(double rho);
  double PXS_NNLO_THRESHOLD_qqbar(double rho);
  double PXS_NNLO_FIT_qqbar(double rho);
  double PXS_NNLO_LO_NLOscales_qqbar(double rho);	
  double PXS_NNLO_EXACT_qqbar(double rho);
  double PXS_NNLO_identical_qqbar(double rho);
  
  double PXS_NLO_gq(double rho);
  double PXS_NNLO_gq(double rho);
  
  double PXS_NNLO_qq(double rho);
  double PXS_NNLO_qqprime(double rho);
  double PXS_NNLO_qqbarprime(double rho);  
  
  // The NNLO fits. Normalization as^4/m^2.
  double f22nl2_FIT_qqbar(double rho);
  double f22nl1_FIT_qqbar(double rho);
  double f22nl0_FIT_qqbar(double rho);
  double f21nl2_FIT_qqbar(double rho);
  double f21nl1_FIT_qqbar(double rho);
  double f21nl0_FIT_qqbar(double rho);
  double f20nl2_FIT_qqbar(double rho);
  double f20nl1_FIT_qqbar(double rho);
  double f20nl0_FIT_qqbar(double rho);
  
  double f22nl1_FIT_gg(double rho);
  double f22nl0_FIT_gg(double rho);
  double f21nl2_FIT_gg(double rho);
  double f21nl1_FIT_gg(double rho);
  double f21nl0_FIT_gg(double rho);
  double f20nl2_FIT_gg(double rho);
  double f20nl1_FIT_gg(double rho);
  double f20nl0_FIT_gg(double rho);
  
  double f20_FIT_qq(double rho);
  double f20_FIT_qqprime(double rho);
  double f20_FIT_qqbarprime(double rho);
  double f22_FIT_qqprime(double rho);
  double f21_FIT_qq(double rho);
  double f21_FIT_qqprime(double rho);
  
  double f20nl0_FIT_gq(double rho);
  double f20nl1_FIT_gq(double rho);
  double f21nl0_FIT_gq(double rho);
  double f21nl1_FIT_gq(double rho);
  double f22nl0_FIT_gq(double rho);
  double f22nl1_FIT_gq(double rho);
  
  // LO results. Normalization as^2/m^2.
  double f0gg(double rho);
  double f0qq(double rho);
  double f0ggLeading(double rho);//the O(beta) term only
  double f0qqLeading(double rho);//the O(beta) term only
  
  // NLO results. Normalization 4*Pi*as^3/m^2.
  double h1(double beta);
  double h2(double beta);
  double f1qqbar(double rho);
  double f1ggbar(double rho);
  double f1gqbar(double rho);
  double f1qq(double rho);
  double f1gg(double rho);
  double f1gq(double rho);
  
  double f1qqCM(double rho);
  double f1ggCM(double rho);
  double f1gqCM(double rho);
  double f1qqNDE(double rho);
  double f1ggNDE(double rho);
  double f1gqNDE(double rho);
};

#endif
