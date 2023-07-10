/*******************************************************************************	
*									       *	
*     A class that computes the resummed perturbative result	               *
*     Output in N-space and in x-space.				               *
*     The implementation is:					               *
*      	- automatic matching to the FO result for requested N^nLO+N^kLL        *
*      	- Only all-order resummation is performed	  	               *
*      	  (will not produce a fixed order truncation of the exponent)          *
*      	- Option to switch ON/OFF the two-loop Coulombic terms	               *
*									       *
*******************************************************************************/

#ifndef guard_class_Resummation_h
#define guard_class_Resummation_h

#include <math.h>
#include <stdexcept>
#include <complex>
#include "psin.h"
#include "lgamma.h"
#include "SubtrFlux.h"

class Resummation 
{
 public:
 Resummation(double _as,
	     SubtrFlux _SF, 
	     int _precision, 
	     double _mt, 
	     double _muF, 
	     double _muR, 
	     double _A, 
	     double _CMP,
	     int _OrderFO, 
	     int _OrderRES, 
	     bool _TwoLoopCoulombs,
	     double _H2qq, 
	     double _H2gg1, 
	     double _H2gg8) :
  as(_as),
    SF(_SF),
    precision(_precision),
    mt(_mt),
    muF(_muF),
    muR(_muR),
    A(_A),
    CMP(_CMP),
    H2qq(_H2qq),
    H2gg1(_H2gg1),
    H2gg8(_H2gg8),
    OrderFO(_OrderFO),
    OrderRES(_OrderRES),
    TwoLoopCoulombs(_TwoLoopCoulombs)
    { 
      create(); 
    }
  
  std::complex<double> resumNqqbar(std::complex<double> N);
  std::complex<double> resumNgg(std::complex<double> N);
  
  double resumqqbar(double rho);
  double resumgg(double rho);
  double resumqqbarSubFlux(double rho);
  double resumggSubFlux(double rho);
  
 private:
  static double pi;
  static double inv2pi;
  SubtrFlux SF;
  int precision;
  double mt, muF, muR, A, CMP, as, lg, lgfr;
  double H2qq, H2gg1, H2gg8;
  bool TwoLoopCoulombs;
  int OrderFO, OrderRES;
  
  void create();
  
  //--N-space expressions for Born, Soft, Coulomb.
  std::complex<double> f0qqN(std::complex<double> N);
  std::complex<double> f1qqN(std::complex<double> N);
  std::complex<double> f2qqN(std::complex<double> N);
  std::complex<double> f1qqCoulN(std::complex<double> N);
  std::complex<double> f2qqCoulN(std::complex<double> N);
  std::complex<double> qqbar(std::complex<double> N);
  std::complex<double> DqqN(std::complex<double> N);
  
  std::complex<double> f0gg1N(std::complex<double> N);
  std::complex<double> f1gg1N(std::complex<double> N);
  std::complex<double> f2gg1N(std::complex<double> N);
  std::complex<double> f1gg1CoulN(std::complex<double> N);
  std::complex<double> f2gg1CoulN(std::complex<double> N);
  std::complex<double> singlet(std::complex<double> N);
  std::complex<double> Dgg1N(std::complex<double> N);
  
  std::complex<double> f0gg8N(std::complex<double> N);
  std::complex<double> f1gg8N(std::complex<double> N);
  std::complex<double> f2gg8N(std::complex<double> N);
  std::complex<double> f1gg8CoulN(std::complex<double> N);
  std::complex<double> f2gg8CoulN(std::complex<double> N);
  std::complex<double> octet(std::complex<double> N);
  std::complex<double> Dgg8N(std::complex<double> N);
  
  std::complex<double> IN(std::complex<double> N);
  
  static double f_MI_noSF_qq (double x, void * par);
  static double f_MI_SF_qq (double x, void * par);
  static double f_MI_noSF_gg (double x, void * par);
  static double f_MI_SF_gg (double x, void * par);	
};

//-- for the Inverse Mellin transform:	
struct IMcomponents {Resummation Res; SubtrFlux SFlux; double rho; double CMP; };

#endif
