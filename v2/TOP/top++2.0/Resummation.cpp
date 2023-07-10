#include "Resummation.h"

using std::string;
using std::domain_error;
typedef std::complex<double>  dcomp;

double Resummation::pi = M_PI;
double Resummation::inv2pi = M_1_PI/2.0;

void Resummation::create()
{	
  lg   = 2*log(muR/mt);
  lgfr = 2*log(muF/muR);
  
  if (OrderRES==!1 && OrderRES==!2 && OrderRES==!3) 
    {
      throw domain_error("(from Class Resummation): unknown value for OrderRES"
			 " (LL,NLL,NNLL are implemented).");
    }
  if (OrderFO==!1 && OrderFO==!2 && OrderFO==!3) 
    {
      throw domain_error("(from Class Resummation): unknown value for OrderFO"
			 " (LO,NLO,NNLO are implemented).");
    }  
}

dcomp Resummation::resumNqqbar(dcomp N)// does all the N-space work.
{
  double g01qq = 1.921046537439473;
  double  CL1  = 0.9037021348921612, 
    CL2  = 3.517994538882637, 
    CL3  = -0.157447215011064, 
    CL4  = -1.0611493499032252, 
    CL5  = 0.001937936937423025;
  dcomp match1 = (g01qq+1.073523998265771*(lg+lgfr))*(1.0-A/(N+A)) 
    -1.2201878970378643*lgfr*(1.0-A/(N+A));
  
  double match2 = H2qq/pow(pi,2) + CL1*pow(lg,2) + CL2*lg 
    + CL3*lg*lgfr + CL4*pow(lgfr,2) + CL5*lgfr;
  dcomp qqCOUL;
  dcomp qqHARD;
  if (OrderFO == 1)//LO
    {
      qqCOUL = f0qqN(N);
      qqHARD = 1.0;
    }
  else if (OrderFO == 2)//NLO
    {
      qqCOUL = f0qqN(N) + as*f1qqCoulN(N);
      qqHARD = 1.0 + as*match1;
    }
  else if (OrderFO == 3)//NNLO
    {
      qqCOUL = f0qqN(N) + as*f1qqCoulN(N) + pow(as,2)*f2qqCoulN(N);
      qqHARD = 1.0 + as*match1 + pow(as,2)*match2;
    }
  
  dcomp qqSOFT = DqqN(N+1.0);
  dcomp exponent = qqCOUL*qqHARD*qqSOFT;
  
  dcomp subtract;	
  if (OrderFO == 1)//LO 
    {
      subtract = f0qqN(N);
    }
  else if (OrderFO == 2)//NLO 
    {
      subtract = f0qqN(N) 
	+ as*(f1qqCoulN(N) + f0qqN(N)*(f1qqN(N+1.0) + match1));
    }
  else if (OrderFO == 3)//NNLO 
    {
      subtract = f0qqN(N) 
	+ as*(f1qqCoulN(N) + f0qqN(N)*(f1qqN(N+1.0) + match1))
	+ pow(as,2)*(f2qqCoulN(N) + f1qqCoulN(N)*(f1qqN(N+1.0) + match1) 
		     +f0qqN(N)*(f2qqN(N+1.0)+f1qqN(N+1.0)*match1+match2));
    }
  
  dcomp f = exponent - subtract;
  f *= 0.38937966*pow(10.,9)*pow(as/mt,2);
  return f;
}

dcomp Resummation::resumNgg(dcomp N)
{
  const double as2 = pow(as,2);
  double  g01gg1 = 2.915123185163655,
    g01gg8 = 4.198085614965150;
  double  CL1gg1 = 0.9439276715786575, 
    CL2gg1 = 3.9032838515013357, 
    CL3gg1 = -0.12984867311560472, 
    CL4gg1 = -1.0737763446942623, 
    CL5gg1 = -1.4322131918654093,
    CL1gg8 = 0.9439276715786576, 
    CL2gg8 = 6.100350166788452, 
    CL3gg8 = -0.1298486731156053, 
    CL4gg8 = -1.0737763446942619, 
    CL5gg8 = -1.583329720375411;
  
  dcomp matchS1 = (g01gg1+1.102400715589848*(lg+lgfr))*(1.0-A/(N+A))
    - 1.2201878970378643*lgfr*(1.0-A/(N+A));
  dcomp matchO1 = (g01gg8+1.102400715589848*(lg+lgfr))*(1.0-A/(N+A))
    - 1.2201878970378643*lgfr*(1.0-A/(N+A));
  
  double matchS2 = H2gg1/pow(pi,2) + CL1gg1*pow(lg,2) + CL2gg1*lg 
    + CL3gg1*lg*lgfr + CL4gg1*pow(lgfr,2) + CL5gg1*lgfr;
  double matchO2 = H2gg8/pow(pi,2) + CL1gg8*pow(lg,2) + CL2gg8*lg 
    + CL3gg8*lg*lgfr + CL4gg8*pow(lgfr,2) + CL5gg8*lgfr;
  dcomp gg1COUL;
  dcomp gg1HARD;
  dcomp gg8COUL;
  dcomp gg8HARD;
  if (OrderFO == 1)//LO
    {
      gg1COUL = f0gg1N(N);
      gg1HARD = 1.0;
      gg8COUL = f0gg8N(N);
      gg8HARD = 1.0;
    }
  else if (OrderFO == 2)//NLO
    {
      gg1COUL = f0gg1N(N)+as*f1gg1CoulN(N);
      gg1HARD = 1.0+as*matchS1;
      gg8COUL = f0gg8N(N)+as*f1gg8CoulN(N);
      gg8HARD = 1.0+as*matchO1;
    }
  else if (OrderFO == 3)//NNLO
    {
      gg1COUL = f0gg1N(N)+as*f1gg1CoulN(N) + as2*f2gg1CoulN(N);
      gg1HARD = 1.0+as*matchS1+as2*matchS2;
      gg8COUL = f0gg8N(N)+as*f1gg8CoulN(N) + as2*f2gg8CoulN(N);
      gg8HARD = 1.0+as*matchO1+as2*matchO2;
    }
  
  dcomp gg1SOFT = Dgg1N(N+1.0);
  dcomp gg8SOFT = Dgg8N(N+1.0);
  
  dcomp exponent = gg1COUL*gg1HARD*gg1SOFT
    + gg8COUL*gg8HARD*gg8SOFT;
  
  dcomp subtract;	
  if (OrderFO == 1)//LO 
    {
      subtract = f0gg1N(N) + f0gg8N(N);
    }
  else if (OrderFO == 2)//NLO 
    {
      subtract = f0gg1N(N) + f0gg8N(N) 
	+ as*(f1gg1CoulN(N) + f1gg8CoulN(N) 
	      + f0gg8N(N)*(f1gg8N(N+1.0) + matchO1) 
	      + f0gg1N(N)*(f1gg1N(N+1.0) + matchS1));
    }
  else if (OrderFO == 3)//NNLO 
    {
      subtract = f0gg1N(N) + f0gg8N(N) 
	+ as*(f1gg1CoulN(N) + f1gg8CoulN(N) 
	      + f0gg8N(N)*(f1gg8N(N+1.0) + matchO1) 
	      + f0gg1N(N)*(f1gg1N(N+1.0) + matchS1))
	+ as2*(f2gg1CoulN(N) + f2gg8CoulN(N)+f1gg8CoulN(N)*(f1gg8N(N+1.0)+matchO1) 
	       + f0gg8N(N)*(f2gg8N(N+1.0) +f1gg8N(N+1.0)*matchO1 + matchO2)
	       +f1gg1CoulN(N)*(f1gg1N(N+1.0) + matchS1) 
	       +f0gg1N(N)*(f2gg1N(N+1.0) + f1gg1N(N+1.0)*matchS1+matchS2));
    }
  
  dcomp f = exponent - subtract;
  f *= 0.38937966*pow(10,9)*pow(as/mt,2);
  return f;
}

/*******************************************************************************	
*					  				       *	
*	      Functions related to the Inverse Mellin transform		       *
*									       *	
*******************************************************************************/

//---------------- WITHOUT subtraction flux:
double Resummation::resumqqbar(double r)
{
  IMcomponents temp_p = {*this, SF, r, CMP};
  gsl_function F;
  F.function = f_MI_noSF_qq;
  F.params = &temp_p;
  //------------------- setup for the integration routine (not meant to be modified):
  int prec = precision;
  double epsabs = 0.0;
  double epsrel = 1/pow(10,prec);
  size_t limit = 1000;
  //-------------------------------
  double result, abserr;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);	
  gsl_integration_qagi (&F, epsabs, epsrel, limit, w, &result, &abserr);
  gsl_integration_workspace_free (w);
  if (abserr/result > epsrel) 
    {
      std::cout 
	<< "Warning (Inv. Mellin Integration [qq, no_SubFlux]): relative error reached " 
	<< abserr/result 
	<< std::endl;
    }
  return result;
}

double Resummation::f_MI_noSF_qq (double x, void * par)//type required by GSL 
{
  IMcomponents * p = (IMcomponents*) par;
  dcomp N((p->CMP)-pow(x,2),x);
  dcomp Jacobian(-2*x,1);
  dcomp f = (p->Res).resumNqqbar(N);
  f *=  Jacobian * pow((p->rho),-N);
  return imag(f) * inv2pi;
}

double Resummation::resumgg(double r)
{
  IMcomponents temp_p = {*this, SF, r, CMP};
  gsl_function F;
  F.function = f_MI_noSF_gg;
  F.params = &temp_p;
  //------------------- setup for the integration routine (not meant to be modified):
  int prec = precision;
  double epsabs = 0.0;
  double epsrel = 1/pow(10,prec);
  size_t limit = 1000;
  //-------------------------------
  double result, abserr;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);	
  gsl_integration_qagi (&F, epsabs, epsrel, limit, w, &result, &abserr);
  gsl_integration_workspace_free (w);
  if (abserr/result > epsrel) 
    {
      std::cout 
	<< "Warning (Inv. Mellin Integration [gg, no_SubFlux]): relative error reached " 
	<< abserr/result 
	<< std::endl;
    }
  return result;
}

double Resummation::f_MI_noSF_gg (double x, void * par) 
{
  IMcomponents * p = (IMcomponents*) par;
  dcomp N((p->CMP)-pow(x,2),x);
  dcomp Jacobian(-2*x,1);
  dcomp f = (p->Res).resumNgg(N);
  f *=  Jacobian * pow((p->rho),-N);
  return imag(f) * inv2pi;
}

//---------------- WITH subtraction flux:
double Resummation::resumqqbarSubFlux(double r)
{
  IMcomponents temp_p = {*this, SF, r, CMP};
  gsl_function F;
  F.function = f_MI_SF_qq;
  F.params = &temp_p;
  //------------------- setup for the integration routine (not meant to be modified):
  int prec = precision;
  double epsabs = 0.0;
  double epsrel = 1/pow(10,prec);
  size_t limit = 1000;
  //-------------------------------
  double result, abserr;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);	
  gsl_integration_qagi (&F, epsabs, epsrel, limit, w, &result, &abserr);
  gsl_integration_workspace_free (w);
  if (abserr/result > epsrel) 
    {
      std::cout 
	<< "Warning (Inv. Mellin Integration [qq, with_SubFlux]): relative error reached " 
	<< abserr/result 
	<< std::endl;
    }
  return result;
}

double Resummation::f_MI_SF_qq (double x, void * par) 
{
  IMcomponents * p = (IMcomponents*) par;
  dcomp N((p->CMP)-pow(x,2) , x);
  dcomp Jacobian(-2*x,1);
  dcomp f = (p->Res).resumNqqbar(N) * (p->SFlux).subtrfluxqqbarN(N);
  f *=  Jacobian * pow((p->rho),-N);
  return imag(f) * inv2pi;
}

double Resummation::resumggSubFlux(double r)
{
  IMcomponents temp_p = {*this, SF, r, CMP};
  gsl_function F;
  F.function = f_MI_SF_gg;
  F.params = &temp_p;
  //------------------- setup for the integration routine (not meant to be modified):
  int prec = precision;
  double epsabs = 0.0;
  double epsrel = 1/pow(10,prec);
  size_t limit = 1000;
  //-------------------------------
  double result, abserr;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);	
  gsl_integration_qagi (&F, epsabs, epsrel, limit, w, &result, &abserr);
  gsl_integration_workspace_free (w);
  if (abserr/result > epsrel) 
    {
      std::cout 
	<< "Warning (Inv. Mellin Integration [gg, with_SubFlux]): relative error reached " 
	<< abserr/result 
	<< std::endl;
    }
  return result;
}

double Resummation::f_MI_SF_gg (double x, void * par) 
{
  IMcomponents * p = (IMcomponents*) par;
  dcomp N((p->CMP)-pow(x,2),x);
  dcomp Jacobian(-2*x,1);
  dcomp f = (p->Res).resumNgg(N) * (p->SFlux).subtrfluxggN(N) ;
  f *=  Jacobian * pow((p->rho),-N);
  return imag(f) * inv2pi;
}

/*******************************************************************************	
*									       *	
*	     Mathematical functions needed for N-space resummation	       *
*									       *	
*******************************************************************************/

dcomp Resummation::f0qqN(dcomp N)//Born_qq
{
  double c0 = 0.3093515553795393;
  dcomp f = c0*exp(lgamma(1.0 + N) 
		   - lgamma(7.0/2.0 + N))*(N + 2.0);
  return f;
}

dcomp Resummation::f1qqN(dcomp N)//NLO expanded soft function (qq)
{
  dcomp LN = log(N);
  dcomp f;
  if (OrderRES == 1)//LL 
    {
      f = 0.8488263631567751*pow(LN,2);
    }
  else//NLL and NNLL are equal 
    {
      f = (0.7581182049282058+0.8488263631567752*lg
	   + 0.8488263631567752*lgfr)*LN
	+0.8488263631567752*pow(LN,2);
    }
  return f;
}

dcomp Resummation::f2qqN(dcomp N)//NNLO expanded soft function (qq)
{
  dcomp LN = log(N);
  dcomp f;
  if (OrderRES == 1)//LL 
    {
      f = 0.3452425516701878*pow(LN,3) + 0.3602530973949787*pow(LN,4);
    }
  else if (OrderRES == 2)//NLL 
    {
      f = (1.2165244709655543 + 1.1613745462374336*lg + 0.3602530973949787*pow(lg,2) 
	   + 0.6435107187321518*lgfr + 0.7205061947899574*lg*lgfr 
	   + 0.3602530973949787*pow(lgfr,2))*pow(LN,2) 
	+ (0.9887532704023396 + 0.7205061947899574*lg + 0.7205061947899574*lgfr)*pow(LN,3) 
	+ 0.3602530973949787*pow(LN,4);
    }
  else if (OrderRES == 3)//NNLL 
    {
      f = (1.0136276961537627+0.9291528646437719*lg
	   + 0.25893191375264096*pow(lg,2)+0.4666295355550373*lgfr
	   - 0.25893191375264096*pow(lgfr,2))*LN
	+ (1.2165244709655547+1.1613745462374339*lg
	   +0.36025309739497874*pow(lg,2)
	   + 0.6435107187321519*lgfr+0.7205061947899575*lg*lgfr
	   + 0.36025309739497874*pow(lgfr,2))*pow(LN,2)
	+ (0.9887532704023397+0.7205061947899575*lg
	   +0.7205061947899575*lgfr)*pow(LN,3)
	+ 0.36025309739497874*pow(LN,4);
    }
  return f;
}

dcomp Resummation::f1qqCoulN(dcomp N)//NLO Coulomb function (qq)
{
  double c0 = -0.03046174197867086;
  dcomp f = c0*(3.0*N+5.0)/((N+1.0)*(N+2.0));
  return f;
}

dcomp Resummation::f2qqCoulN(dcomp N)//NNLO Coulomb function (qq)
{
  dcomp f;
  if (TwoLoopCoulombs) 
    {
      f = qqbar(N)/pow(pi,2) 
	+ 0.61009394851893212045*lg*f1qqCoulN(N);
    }
  else
    {
      f = (0.0,0.0);
    }
  return f;
}

dcomp Resummation::qqbar(dcomp N)// input to the NNLO Coulomb function
{
  dcomp f =
    (0.10228905314089776*(-5.901687825416445 + N)*(0.46730279875222375 + 0.20725313959425187*N + pow(N,2))
     *(2.600802746889236 + 3.1884118574550286*N + pow(N,2)))/
    (pow(N,2)*pow(1. + 1.*N,2)*pow(2. + 1.*N,2)) + 0.13950765657164715*exp(lgamma(N) - lgamma(0.5 + N)) 
    + 0.02682101416269378*(-1. + 1.*N)*exp(lgamma(-1.0 + N) - lgamma(1.5 + N)) - 0.0930051043810981*exp(lgamma(N) - lgamma(1.5 + N)) 
    - 1.0336743034352796*(-1. + 1.*N)*exp(lgamma(-1. + N) - lgamma(2.5 + N)) + 0.03487691414291179*exp(lgamma(N) - lgamma(2.5 + N)) 
    + 1.013558542813259*exp(lgamma(N) - lgamma(3.5 + N)) 
    - (0.3668438084603374*(1.0000000000000004 + N)*(1.9999999999999996 + N)*psi0(1. - N))/(N*pow(1. + 1.*N,2)*pow(2. + 1.*N,2)) 
    + (0.3668438084603374*(1.0000000000000004 + N)*(1.9999999999999996 + N)*psi0(2. - N))/(N*pow(1. + 1.*N,2)*pow(2. + 1.*N,2)) 
    - (0.3668438084603374*(1.0000000000000004 + N)*(1.9999999999999996 + N)*psi0(-1. + N))/(N*pow(1. + 1.*N,2)*pow(2. + 1.*N,2)) 
    - (0.550265712690506*(-0.3333333333333336 + 1.*N)*psi0(1. + N))/(N*(1. + 1.*N)) 
    - 0.7350242063324187*(-1. + 1.*N)*exp(lgamma(-1. + N) - lgamma(1.5 + N))*psi0(1.5 + N) 
    + 1.4700484126648372*(-1. + 1.*N)*exp(lgamma(-1. + N) - lgamma(2.5 + N))*psi0(2.5 + N) 
    - 0.9187802579155233*exp(lgamma(N) - lgamma(3.5 + N))*psi0(3.5 + N);
  return f;
}

//The all-order soft exponent (qq)
dcomp Resummation::DqqN(dcomp N)
{
  double b0 = 0.6100939485189321;
  dcomp L = b0*as*log(N);
  dcomp L1 = log(1.0-2.0*L);
  
  dcomp g1 = 1.391304347826087 - 1.391304347826087*L1 
    + (0.6956521739130435*L1)/L;
  
  dcomp g2 = -0.3383945767837618*L - 0.7905099622584375*L1 
    + 0.2288156488863319*pow(L1,2) - 0.6956521739130435*L1*lg 
    + 1.391304347826087*L*lgfr;
  
  dcomp g3 = (1.0003181247683397*L + 0.17113406096200867*pow(L,2) 
	      - 0.3305553698884765*L1 + 0.026572744359957722*L*L1 
	      + 0.09183453093340334*pow(L1,2) + 0.9645706884356047*L*lg 
	      - 0.2791980854239677*L1*lg + 0.4244131815783876*L*pow(lg,2) 
	      + 0.7648486543553337*L*lgfr - 1.5296973087106673*pow(L,2)*lgfr 
	      - 0.4244131815783876*L*pow(lgfr,2) 
	      + 0.8488263631567752*pow(L,2)*pow(lgfr,2))/(1.0 - 2.0*L);
  
  dcomp f;
  if      (OrderRES == 1) {f = exp(log(N)*g1); }
  else if (OrderRES == 2) {f = exp(log(N)*g1+g2); }
  else if (OrderRES == 3) {f = exp(log(N)*g1+g2+as*g3); }  
  return f;
}

//-------------------------------- gg-singlet (functions similar to the qq case; see above)

dcomp Resummation::f0gg1N(dcomp N)
{
  double c0 = 5.568327996831708;
  dcomp f = (14.0 + 20.0*N + 9.0*pow(N,2) + pow(N,3))
    /(192.0*(6.0 + 11.0*N + 6.0*pow(N,2) + pow(N,3)));
  f *= c0 * exp(lgamma(1.0 + N) - lgamma(2.5 + N));
  return f;
}

dcomp Resummation::f1gg1N(dcomp N)
{
  dcomp LN = log(N);
  dcomp f;
  if (OrderRES == 1)//LL 
    {
      f = 1.9098593171027443*pow(LN,2);
    }
  else//NLL and NNLL are equal
    {
      f = (-0.44282577065212303+1.909859317102744*lg 
	   + 1.909859317102744*lgfr)*LN+1.909859317102744*pow(LN,2);
    }
  return f;
}

dcomp Resummation::f2gg1N(dcomp N)
{
  dcomp LN = log(N);
  double lg2 = pow(lg,2);
  double lgfr2 = pow(lgfr,2);
  dcomp f;	
  if (OrderRES == 1)//LL 
    {
      f = 0.776795741257923*pow(LN,3) + 1.8237813055620804*pow(LN,4);
    }
  else if (OrderRES == 2)//NLL 
    {
      f = (0.8777984636525643 + 0.3194586879537237*lg 
	   + 1.82378130556208*pow(lg,2) - 0.8457349239331605*lgfr 
	   + 3.64756261112416*lg*lgfr + 1.82378130556208*pow(lgfr,2))*pow(LN,2) 
	+ (-0.06893918267523769 + 3.6475626111241604*lg 
	   + 3.6475626111241604*lgfr)*pow(LN,3) 
	+ 1.8237813055620804*pow(LN,4);
    }
  else if (OrderRES == 3)//NNLL 
    {
      f = (-1.4255288604709904
	   +0.7797511320757415*lg+0.5825968059434422*lg2
	   +1.0499164549988342*lgfr-0.5825968059434422*lgfr2)*LN
	+(0.8777984636525642+0.31945868795372423*lg
	  +1.82378130556208*lg2-0.8457349239331604*lgfr
	  +3.64756261112416*lg*lgfr+1.82378130556208*lgfr2)*pow(LN,2)
	+(-0.06893918267523702+3.64756261112416*lg
	  +3.64756261112416*lgfr)*pow(LN,3)
	+1.82378130556208*pow(LN,4);
    }
  return f;
}

dcomp Resummation::f1gg1CoulN(dcomp N)
{
  double c0 = 0.03426945972600472;
  dcomp  f  = c0*(2.0*IN(N)+2.0*IN(N+1.0)
		  -IN(N+2.0)-2.0/(N+1.0)-2.0/(N+2.0));
  return f;
}

dcomp Resummation::f2gg1CoulN(dcomp N)
{
  dcomp f;
  if (TwoLoopCoulombs) 
    {
      f = singlet(N)/pow(pi,2) 
	+ 0.61009394851893212045*lg*f1gg1CoulN(N);
    }
  else
    {
      f = (0.0,0.0);
    }
  return f;
}

dcomp Resummation::singlet(dcomp N)
{
  dcomp f =
    (-5.981149050983761e-7*(2.408307241524723e9 + 3.0419325750721016e9*N 
			    - 1.717903402540184e10*pow(N,2) 
			    - 4.69908784214079e10*pow(N,3) - 4.505839365414517e10*pow(N,4) - 
			    1.1583253155331848e10*pow(N,5) + 1.3884869378878479e10*pow(N,6) 
			    + 1.5493888085853691e10*pow(N,7) 
			    + 7.482276506396107e9*pow(N,8) + 2.0890941783972125e9*pow(N,9) 
			    + 3.4778172282993746e8*pow(N,10) + 3.2135436286905903e7*pow(N,11) 
			    + 1.2715287879426382e6*pow(N,12)))/
    (pow(N,2)*pow(1. + N,3)*pow(2. + N,3)*pow(3. + N,3)*(4. + N)*(5. + N)) 
    + 0.8347859153934222*exp(lgamma(N) - 1.*lgamma(0.5 + N)) - 
    0.020020907073856613*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(1.5 + N)) 
    + 0.716386280388247*exp(lgamma(N) - 1.*lgamma(1.5 + N)) - 
    1.32694186664566*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(2.5 + N)) 
    - 2.238463911460011*exp(lgamma(N) - 1.*lgamma(2.5 + N)) + 
    11.663432027295865*exp(lgamma(N) - 1.*lgamma(3.5 + N)) 
    - 9.198642375521679*exp(lgamma(N) - 1.*lgamma(4.5 + N)) 
    - 21.835333397021838*exp(lgamma(N) - 1.*lgamma(5.5 + N)) + 
    (0.7498853660372788*psi0(1. - 1.*N))/(N*(1. + N)*(2. + N)) 
    - (0.7498853660372788*psi0(2. - 1.*N))/(N*(1. + N)*(2. + N)) 
    + (0.7498853660372788*psi0(-1. + N))/(N*(1. + N)*(2. + N)) + 
    1.317126145774745*exp(lgamma(N) - 1.*lgamma(0.5 + N))*psi0(0.5 + N) 
    + 1.0972356541644999*exp(lgamma(N) - 1.*lgamma(0.5 + N))*pow(psi0(0.5 + N),2) + 
    (1.7943447152951283e-6*(-583999. + 229379.*N)*psi0(1. + N))/(N*(1. + N)) + 
    (0.9673222444161708*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(1.5 + N)) 
     + 0.020019050234414578*exp(lgamma(N) - 1.*lgamma(1.5 + N)))*psi0(1.5 + N) - 
    0.5486178270822499*exp(lgamma(N) - 1.*lgamma(1.5 + N))*pow(psi0(1.5 + N),2) + 
    (2.096448055349478*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(2.5 + N)) 
     + 0.1928821341446239*exp(lgamma(N) - 1.*lgamma(2.5 + N)))*psi0(2.5 + N) - 
    0.27430891354112497*exp(lgamma(N) - 1.*lgamma(2.5 + N))*pow(psi0(2.5 + N),2) + 
    (-0.5233169138033914*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(3.5 + N)) 
     - 10.920995260745565*exp(lgamma(N) - 1.*lgamma(3.5 + N)))*psi0(3.5 + N) + 
    0.6857722838528125*exp(lgamma(N) - 1.*lgamma(3.5 + N))*pow(psi0(3.5 + N),2) + 
    psi0(N)*((0.6100772032003436*(86. + 156.*N + 105.*pow(N,2) + 30.*pow(N,3) 
				  + 3.*pow(N,4)))/(pow(1. + N,2)*pow(2. + N,2)*pow(3. + N,2)) - 
	     1.317126145774745*exp(lgamma(N) - 1.*lgamma(0.5 + N)) 
	     - 0.41865353104271313*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(1.5 + N)) - 
	     0.020019050234414578*exp(lgamma(N) - 1.*lgamma(1.5 + N)) 
	     - 0.20932676552135657*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(2.5 + N)) - 
	     0.1928821341446239*exp(lgamma(N) - 1.*lgamma(2.5 + N)) 
	     + 0.5233169138033914*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(3.5 + N)) + 
	     0.7565142489026848*exp(lgamma(N) - 1.*lgamma(3.5 + N)) 
	     - 1.0972356541644999*exp(lgamma(N) - 1.*lgamma(0.5 + N))*psi0(0.5 + N) + 
	     0.5486178270822499*exp(lgamma(N) - 1.*lgamma(1.5 + N))*psi0(1.5 + N) 
	     + 0.27430891354112497*exp(lgamma(N) - 1.*lgamma(2.5 + N))*psi0(2.5 + N) - 
	     0.6857722838528125*exp(lgamma(N) - 1.*lgamma(3.5 + N))*psi0(3.5 + N)) 
    + (6.95890662602504*psi0(4. + N))/(N*(1. + N)*(2. + N)*(3. + N)) + 
    8.688400895129085*exp(lgamma(N) - 1.*lgamma(4.5 + N))*psi0(4.5 + N) 
    - (9.599758581586658*psi0(5. + N))/(N*(1. + N)*(2. + N)*(3. + N)*(4. + N)) + 
    14.281625237842384*exp(lgamma(N) - 1.*lgamma(5.5 + N))*psi0(5.5 + N) 
    - (10.42499924828747*psi0(6. + N))/(N*(1. + N)*(2. + N)*(3. + N)*(4. + N)*(5. + N)) + 
    0.7083788999659292*psi1(N) - 1.0972356541644999*exp(lgamma(N) - 1.*lgamma(0.5 + N))*psi1(0.5 + N) 
    - (1.8302316096010307*psi1(1. + N))/N + 
    0.5486178270822499*exp(lgamma(N) - 1.*lgamma(1.5 + N))*psi1(1.5 + N) 
    + (1.8302316096010307*psi1(2. + N))/(N*(1. + N)) + 
    0.27430891354112497*exp(lgamma(N) - 1.*lgamma(2.5 + N))*psi1(2.5 + N) 
    + (1.2201544064006873*psi1(3. + N))/(N*(1. + N)*(2. + N)) - 
    0.6857722838528125*exp(lgamma(N) - 1.*lgamma(3.5 + N))*psi1(3.5 + N) 
    - (3.6604632192020614*psi1(4. + N))/(N*(1. + N)*(2. + N)*(3. + N));
  return f;
}

dcomp Resummation::Dgg1N(dcomp N)
{
  double b0 = 0.6100939485189321;
  dcomp L = b0*as*log(N);
  dcomp L2 = pow(L,2);
  dcomp L1 = log(1.0-2.0*L);
  dcomp L12 = pow(L1,2);
  double lg2 = pow(lg,2);
  double lgfr2 = pow(lgfr,2);
  
  dcomp g1 = 3.130434782608696 
    - 3.130434782608696*L1 
    + (1.565217391304348*L1)/L;
  
  dcomp g2 = -0.761387797763464*L 
    - 0.01777784986409303*L1 
    + 0.5148352099942467*L12 
    - 1.565217391304348*L1*lg 
    + 3.130434782608696*L*lgfr;
  
  dcomp g3 = (-2.410631550362379*L 
	      + 0.3850516371645214*L2 
	      - 0.03702942851965361*L1 
	      + 0.05978867480990487*L*L1 
	      + 0.2066276946001575*L12 
	      + 0.02169231723952255*L*lg
	      - 0.6281956922039271*L1*lg 
	      + 0.954929658551372*L*lg2 
	      + 1.7209094722995002*L*lgfr 
	      - 3.4418189445990004*L2*lgfr 
	      - 0.954929658551372*L*lgfr2 
	      + 1.909859317102744*L2*lgfr2)/(1.0-2.0*L);
  
  dcomp f;
  if      (OrderRES == 1) {f = exp(log(N)*g1); }
  else if (OrderRES == 2) {f = exp(log(N)*g1+g2); }
  else if (OrderRES == 3) {f = exp(log(N)*g1+g2+as*g3); }
  return f;
}

//=============================== gg-octet (functions similar to the qq case; see above)

dcomp Resummation::f0gg8N(dcomp N)
{
  dcomp c0 = 5.568327996831708;
  dcomp f = (292.0 + 455.0*N + 262.0*pow(N,2) 
	     + 62.0*pow(N,3) + 5.0*pow(N,4))
    /(192.0*(30.0 + 67.0*N + 52.0*pow(N,2) 
	     + 17.0*pow(N,3) + 2.0*pow(N,4)));
  
  f *= c0 * exp(lgamma(1.0 + N) - lgamma(2.5 + N));
  return f;
}

dcomp Resummation::f1gg8N(dcomp N)
{
  dcomp LN = log(N);
  dcomp f;
  if (OrderRES == 1)//LL 
    {
      f = 1.9098593171027443*pow(LN,2);
    }
  else//NLL and NNLL are equal
    {
      f = (0.512103887899249
	   +1.909859317102744*lg
	   +1.909859317102744*lgfr)*LN
	+1.909859317102744*pow(LN,2);
    }
  return f;
}

dcomp Resummation::f2gg8N(dcomp N)
{
  dcomp LN = log(N);
  double lg2 = pow(lg,2);
  double lgfr2 = pow(lgfr,2);
  dcomp f;	
  if (OrderRES == 1)//LL 
    {
      f = 0.776795741257923*pow(LN,3) + 1.8237813055620804*pow(LN,4);
    }
  else if (OrderRES == 2)//NLL 
    {
      f = (1.493473134019946 + 2.1432399935158033*lg 
	   + 1.82378130556208*pow(lg,2) + 0.978046381628919*lgfr 
	   + 3.64756261112416*lg*lgfr + 1.82378130556208*pow(lgfr,2))*pow(LN,2) 
	+ (1.754842122886842 + 3.6475626111241604*lg 
	   + 3.6475626111241604*lgfr)*pow(LN,3) 
	+ 1.8237813055620804*pow(LN,4);
    }
  else if (OrderRES == 3)//NNLL 
    {
      f = (0.22166721811432488
	   +1.3623479380191839*lg
	   +0.5825968059434422*lg2
	   +1.0499164549988342*lgfr
	   -0.5825968059434422*lgfr2)*LN
	+(1.4934731340199452+2.143239993515804*lg
	  +1.82378130556208*lg2+0.9780463816289195*lgfr
	  +3.64756261112416*lg*lgfr+1.82378130556208*lgfr2)*pow(LN,2)
	+(1.754842122886842+3.64756261112416*lg
	  +3.64756261112416*lgfr)*pow(LN,3)
	+1.82378130556208*pow(LN,4);
    }
  return f;
}

dcomp Resummation::f1gg8CoulN(dcomp N)
{
  double c0 = -0.004283682465750590;
  dcomp f = c0*(14.0*IN(N)+14.0*IN(N+1.0)
		+2.0*IN(N+2.0)-26.0/(N+1.0)-29.0/(N+2.0));
  return f;
}

dcomp Resummation::f2gg8CoulN(dcomp N)
{
  dcomp f;
  if (TwoLoopCoulombs) 
    {
      f = octet(N)/pow(pi,2) 
	+ 0.61009394851893212045*lg*f1gg8CoulN(N);
    }
  else
    {
      f = (0.0,0.0);
    }
  return f;
}

dcomp Resummation::octet(dcomp N)
{
  dcomp f = 
    (-1.4952872627459403e-7*(3.727309102595295e10 + 2.3418883719642007e11*N 
			     + 6.120918386456074e11*pow(N,2) + 9.002441922189667e11*pow(N,3) 
			     + 8.498748235506926e11*pow(N,4) + 
			     5.523881395195466e11*pow(N,5) + 2.5614771724907263e11*pow(N,6) 
			     + 8.56384507179771e10*pow(N,7) + 2.0414291421397717e10*pow(N,8) 
			     + 3.3550010934929323e9*pow(N,9) + 
			     3.5651725146765655e8*pow(N,10) + 2.165455573096964e7*pow(N,11) 
			     + 556157.9432938319*pow(N,12)))/
    (pow(N,2)*pow(1. + N,3)*pow(2. + N,3)*pow(3. + N,3)*(4. + N)*(5. + N)) 
    + 0.03234397825328282*exp(lgamma(N) - 1.*lgamma(0.5 + N)) - 
    0.018423038707754948*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(1.5 + N)) 
    - 0.4956743969380271*exp(lgamma(N) - 1.*lgamma(1.5 + N)) + 
    1.1133861221926307*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(2.5 + N)) 
    + 1.0381220610749051*exp(lgamma(N) - 1.*lgamma(2.5 + N)) - 
    3.609606430559659*exp(lgamma(N) - 1.*lgamma(3.5 + N)) 
    + 2.7928101078869383*exp(lgamma(N) - 1.*lgamma(4.5 + N)) 
    - 2.085557839982852*exp(lgamma(N) - 1.*lgamma(5.5 + N)) - 
    (2.358908364792011*psi0(1. - 1.*N))/(N*(1. + N)*(2. + N)) 
    + (2.358908364792011*psi0(2. - 1.*N))/(N*(1. + N)*(2. + N)) 
    - (2.358908364792011*psi0(-1. + N))/(N*(1. + N)*(2. + N)) - 
    1.646407682218431*exp(lgamma(N) - 1.*lgamma(0.5 + N))*psi0(0.5 + N) 
    - 1.371544567705625*exp(lgamma(N) - 1.*lgamma(0.5 + N))*pow(psi0(0.5 + N),2) - 
    (4.485861788237821e-7*(-1.863512e6 + 284395.*N)*psi0(1. + N))/(N*(1. + N)) + 
    (-0.3324275924927336*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(1.5 + N)) 
     - 0.040038100468829156*exp(lgamma(N) - 1.*lgamma(1.5 + N)))*psi0(1.5 + N) + 
    1.0972356541644999*exp(lgamma(N) - 1.*lgamma(1.5 + N))*pow(psi0(1.5 + N),2) + 
    (-1.0600943175417972*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(2.5 + N)) 
     + 0.48220533536155974*exp(lgamma(N) - 1.*lgamma(2.5 + N)))*psi0(2.5 + N) - 
    0.6857722838528125*exp(lgamma(N) - 1.*lgamma(2.5 + N))*pow(psi0(2.5 + N),2) + 
    (-0.13082922845084785*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(3.5 + N)) 
     + 2.4088034520243164*exp(lgamma(N) - 1.*lgamma(3.5 + N)))*psi0(3.5 + N) + 
    0.17144307096320313*exp(lgamma(N) - 1.*lgamma(3.5 + N))*pow(psi0(3.5 + N),2) + 
    psi0(N)*((-0.07177378861180513*(319. + 600.*N + 426.*pow(N,2) + 132.*pow(N,3) 
				    + 15.*pow(N,4)))/(pow(1. + N,2)*pow(2. + N,2)*pow(3. + N,2)) + 
	     1.646407682218431*exp(lgamma(N) - 1.*lgamma(0.5 + N)) 
	     + 0.8373070620854263*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(1.5 + N)) + 
	     0.040038100468829156*exp(lgamma(N) - 1.*lgamma(1.5 + N)) 
	     - 0.5233169138033914*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(2.5 + N)) - 
	     0.48220533536155974*exp(lgamma(N) - 1.*lgamma(2.5 + N)) 
	     + 0.13082922845084785*(-1. + N)*exp(lgamma(-1. + N) - 1.*lgamma(3.5 + N)) + 
	     0.1891285622256712*exp(lgamma(N) - 1.*lgamma(3.5 + N)) 
	     + 1.371544567705625*exp(lgamma(N) - 1.*lgamma(0.5 + N))*psi0(0.5 + N) - 
	     1.0972356541644999*exp(lgamma(N) - 1.*lgamma(1.5 + N))*psi0(1.5 + N) 
	     + 0.6857722838528125*exp(lgamma(N) - 1.*lgamma(2.5 + N))*psi0(2.5 + N) - 
	     0.17144307096320313*exp(lgamma(N) - 1.*lgamma(3.5 + N))*psi0(3.5 + N)) 
    + (1.7716292083718497*psi0(4. + N))/(N*(1. + N)*(2. + N)*(3. + N)) - 
    2.192595870213594*exp(lgamma(N) - 1.*lgamma(4.5 + N))*psi0(4.5 + N) 
    + (1.402980189442108*psi0(5. + N))/(N*(1. + N)*(2. + N)*(3. + N)*(4. + N)) + 
    1.476625752922254*exp(lgamma(N) - 1.*lgamma(5.5 + N))*psi0(5.5 + N) 
    - (1.2264704997985258*psi0(6. + N))/(N*(1. + N)*(2. + N)*(3. + N)*(4. + N)*(5. + N)) + 
    0.11068420311967643*psi1(N) 
    + 1.371544567705625*exp(lgamma(N) - 1.*lgamma(0.5 + N))*psi1(0.5 + N) 
    + (1.076606829177077*psi1(1. + N))/N - 
    1.0972356541644999*exp(lgamma(N) - 1.*lgamma(1.5 + N))*psi1(1.5 + N) 
    - (1.7225709266833231*psi1(2. + N))/(N*(1. + N)) + 
    0.6857722838528125*exp(lgamma(N) - 1.*lgamma(2.5 + N))*psi1(2.5 + N) 
    + (1.4354757722361025*psi1(3. + N))/(N*(1. + N)*(2. + N)) - 
    0.17144307096320313*exp(lgamma(N) - 1.*lgamma(3.5 + N))*psi1(3.5 + N) 
    - (0.4306427316708308*psi1(4. + N))/(N*(1. + N)*(2. + N)*(3. + N));
  return f;
}

dcomp Resummation::Dgg8N(dcomp N)
{
  double b0 = 0.6100939485189321;
  dcomp L = b0*as*log(N);
  dcomp L2 = pow(L,2);
  dcomp L1 = log(1.0-2.0*L);
  dcomp L12 = pow(L1,2);
  
  dcomp g1 = 3.130434782608696 
    - 3.130434782608696*L1 
    + (1.565217391304348*L1)/L;
  
  dcomp g2 = -0.761387797763464*L 
    - 0.8003865455162669*L1 
    + 0.5148352099942467*L12 
    - 1.565217391304348*L1*lg 
    + 3.130434782608696*L*lgfr;
  
  dcomp g3 = (-0.33892162543298454*L 
	      + 0.3850516371645214*L2 
	      - 0.3511272746216172*L1 
	      + 0.05978867480990487*L*L1 
	      + 0.2066276946001575*L12 
	      + 0.9766219757908946*L*lg 
	      - 0.6281956922039271*L1*lg 
	      + 0.954929658551372*L*pow(lg,2) 
	      + 1.7209094722995002*L*lgfr 
	      - 3.4418189445990004*L2*lgfr 
	      - 0.954929658551372*L*pow(lgfr,2) 
	      + 1.909859317102744*L2*pow(lgfr,2))/(1.0-2.0*L);
  dcomp f;
  if (OrderRES == 1) { f = exp(log(N)*g1); }
  else if (OrderRES == 2) { f = exp(log(N)*g1+g2); }
  else if (OrderRES == 3) {f = exp(log(N)*g1+g2+as*g3); } 
  return f;
}

// The integral I(N) defined in the Bonciani et al paper (1998)
// (approx. Mellin transform of the gg NLO Coulomb function) 
dcomp Resummation::IN(dcomp N)
{
  double c0 = 1.772453850905516;
  double d[] = {0.0, 0.9991, -0.4828, 0.2477, -0.0712};
  
  dcomp f = c0*exp(lgamma(N+1.0)-lgamma(N+1.5)) 
    * (psi0(N+1.5)-psi0(N+1.0));
  for (int i=1; i<=4; ++i) 
    {
      f += 2*d[i]*exp(lgamma(N+1.0)
		      +lgamma(0.5*(i+1.0))-lgamma(N+0.5*(i+3.0)));
    }
  return f;
}
