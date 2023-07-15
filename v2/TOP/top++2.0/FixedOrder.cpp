#include "FixedOrder.h"

using std::string;
using std::domain_error;

double FixedOrder::pi = M_PI;

void FixedOrder::create() 
{
  lg   = 2*log(muR/mt);
  lgfr = 2*log(muF/muR);
}

/*******************************************************************************	
*									       *	
*		  The Fixed Order results for all partonic channels	       *
*									       *	
*******************************************************************************/

double FixedOrder::FOgg(double rho)
{
  double f = 0.0;	
  if (LO)   { f += PXS_LO_gg(rho) ;}
  if (NLO)  { f += PXS_NLO_gg(rho) ;}
  if (NNLO) { f += PXS_NNLO_EXACT_gg(rho) ;}
  return f * 0.38937966*pow(10,9)*pow(as/mt,2);
}

double FixedOrder::FOqqbar(double rho)
{
  double f = 0.0;	
  if (LO)   { f += PXS_LO_qqbar(rho) ;}
  if (NLO)  { f += PXS_NLO_qqbar(rho) ;}
  if (NNLO) { f += PXS_NNLO_EXACT_qqbar(rho) ;}  
  return f * 0.38937966*pow(10,9)*pow(as/mt,2);
}

double FixedOrder::FOgq(double rho)
{
  double f = 0.0;	
  if (LO)  { f += 0;}
  if (NLO) { f += PXS_NLO_gq(rho) ;}
  if (NNLO){ f += PXS_NNLO_gq(rho) ;}
  return f * 0.38937966*pow(10,9)*pow(as/mt,2);
}

double FixedOrder::FOqq(double rho)
{
  double f = 0.0;	
  if (LO)  { f += 0;}
  if (NLO) { f += 0;}
  if (NNLO){ f += PXS_NNLO_qq(rho) ;}
  return f * 0.38937966*pow(10,9)*pow(as/mt,2);
}

double FixedOrder::FOqqprime(double rho)
{
  double f = 0.0;	
  if (LO)  { f += 0;}
  if (NLO) { f += 0;}
  if (NNLO){ f += PXS_NNLO_qqprime(rho) ;}
  return f * 0.38937966*pow(10,9)*pow(as/mt,2);
}

double FixedOrder::FOqqbarprime(double rho)
{
  double f = 0.0;	
  if (LO)  { f += 0;}
  if (NLO) { f += 0;}
  if (NNLO){ f += PXS_NNLO_qqbarprime(rho) ;}
  return f * 0.38937966*pow(10,9)*pow(as/mt,2);
}

/*******************************************************************************	
*									       *	
*		  LO Partonic X-sections		                       *
*									       *	
*******************************************************************************/	

double FixedOrder::PXS_LO_qqbar(double rho)
{
  return f0qq(rho);
}

double FixedOrder::PXS_LO_gg(double rho)
{
  return f0gg(rho);
}

/*******************************************************************************	
*									       *	
*		  NLO Partonic X-sections		                       *
*									       *	
*******************************************************************************/	

double FixedOrder::PXS_NLO_gg(double rho)
{
  double f = f0gg(rho)*(-1.2201878970378643*lgfr*as) 
    + 4*pi*as*(f1gg(rho)+f1ggbar(rho)*(lg+lgfr));
  return f;
}

double FixedOrder::PXS_NLO_qqbar(double rho)
{
  double f = f0qq(rho)*(-1.2201878970378643*lgfr*as)
    +4*pi*as*(f1qq(rho)+f1qqbar(rho)*(lg+lgfr));
  return f;
}

double FixedOrder::PXS_NLO_gq(double rho)
{
  return 4*pi*as*(f1gq(rho)+f1gqbar(rho)*(lg+lgfr));
}

/*******************************************************************************	
*									       *	
*	 Exact NNLO results: gq reaction                                       *
*									       *	
*******************************************************************************/	

double FixedOrder::PXS_NNLO_gq(double rho)
{	
  double nl=5.0;
  // the functions fij are normalized as: as^n/m^2
  double f10 = 4*pi*f1gq(rho);
  double f11 = 4*pi*f1gqbar(rho);	  
  double f20 = f20nl0_FIT_gq(rho) + nl * f20nl1_FIT_gq(rho);
  double f21 = f21nl0_FIT_gq(rho) + nl * f21nl1_FIT_gq(rho);
  double f22 = f22nl0_FIT_gq(rho) + nl * f22nl1_FIT_gq(rho);
  
  double f = f20 + f21*(lg + lgfr) + f22*(pow(lg,2) + 2.*lg*lgfr + pow(lgfr,2)) 
    + f10*(-2.626056561016273*lgfr + 0.15915494309189535*lgfr*nl) 
    + f11*(-2.626056561016273*lg*lgfr - 2.626056561016273*pow(lgfr,2) 
	   + (0.15915494309189535*lg*lgfr + 0.15915494309189535*pow(lgfr,2))*nl);
  return f * pow(as,2);
}

//=======================================================================
// The NNLO gq functions (i.e. the ones without scale-logs):
//=======================================================================

double FixedOrder::f20nl0_FIT_gq(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lrho3 = pow(lrho,3);
  double lbe2 = pow(lbe,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b5 = pow(b,5);
  double b6 = pow(b,6);
  double b7 = pow(b,7);
  
  double f = 28.0998*b2 + 24.1753*b3 - 12.3211*b5 - 49.909*b7 
    + 11.7853*b3*lbe + 28.6697*b6*lbe2 - 3.82993*lrho3*rho 
    + lrho2*(-9.80339*rho - 76.7407*rho2) 
    + lrho*(-1.68957 + 30.6335*rho2);
  return f;
}

double FixedOrder::f20nl1_FIT_gq(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lrho3 = pow(lrho,3);
  double lbe2 = pow(lbe,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b7 = pow(b,7);
  
  double f = 0.363838*b2 - 1.44391*b3 + 1.1146*b7 - 0.309165*b3*lbe 
    + 0.990057*b4*lbe2 + 0.0401411*lrho3*rho + 0.362183*lrho*rho2 
    + lrho2*(0.194867*rho + 1.57274*rho2);
  return f;
}

//==========================================================================
// The scaling functions ~Log(muF^2/m^2) for gq-reaction.
//
// Implemented fits: as in ``Hathor" (Aliev et al. arXiv:1007.1327v1)
//==========================================================================

double FixedOrder::f21nl0_FIT_gq(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lbe2 = pow(lbe,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  double b6 = pow(b,6);
  double b7 = pow(b,7);
  
  double f = -0.00272708*b2 - 0.194727*b3 
    - 7.74780189742068692567*b4 - 32.9814144136482376215*b5 
    - 2168.37139107242508102*b6 + 2213.66358353965767595*b7 
    + 0.62267243714534082136754351670760392811*b3*lbe - 1.469367779446109471946*b4*lbe 
    - 82.4683103579544036772*b5*lbe - 739.731611466006394139*b6*lbe 
    - 1201.795933886655196561*b7*lbe - 0.56735789898499726176870600831868545725*b3*lbe2 
    + 215.848631557653764823*b3*lrho - 1146.873278298064968087*b4*lrho 
    + 1744.26195133738780366*b5*lrho - 1000.346862873385464604*b6*lrho 
    + 187.933944366753454585*b7*lrho + 291.663849809244439603*b3*lrho2 
    - 772.77269566517296654*b4*lrho2 + 653.225988602083032794*b5*lrho2 
    - 172.118699492821641978*b6*lrho2;
  return f;
}

double FixedOrder::f21nl1_FIT_gq(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lbe2 = pow(lbe,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  double b6 = pow(b,6);
  double b7 = pow(b,7);
  
  double f = 0.00554947127181299611557*b3 + 0.0225430311833103251235*b4 
    - 0.63483974233126651067*b5 + 9.99566887547610873574*b6 
    - 9.40032480635199279045*b7 - 0.0112980360836839281564331066591599701444*b3*lbe 
    + 0.00425431823252394833895*b4*lbe + 0.252351741613262253035*b5*lbe 
    + 2.40450791609598725401*b6*lbe + 4.53123563973969674337*b7*lbe 
    - 1.382674579044874428*b3*lrho + 6.3852062287204969368*b4*lrho 
    - 8.30987940966125596394*b5*lrho + 3.65509528520063341522*b6*lrho 
    - 0.347592328515009762391*b7*lrho - 1.263737185806346314577*b3*lrho2 
    + 3.10439060379970187549*b4*lrho2 - 2.30142593570054254576*b5*lrho2 
    + 0.460779458031490856105*b6*lrho2;
  return f;
}

double FixedOrder::f22nl0_FIT_gq(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lbe2 = pow(lbe,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  double b6 = pow(b,6);
  double b7 = pow(b,7);
  
  double f = -0.0956665797118411055*b3 + 0.0796230395754896639*b4 
    - 4.65135359743332977*b5 + 54.2282042897508657*b6 
    - 50.3664336690071789*b7 + 0.094559649830832876961451001386447576209*b3*lbe 
    + 0.01454877297063897196*b4*lbe + 1.090238571394724594*b5*lbe 
    + 12.39185329777043355*b6*lbe + 25.3325922738649971*b7*lbe 
    - 8.19094150408542609*b3*lrho + 34.1888813941803399*b4*lrho 
    - 44.3046767125567972*b5*lrho + 20.6517749275371715*b6*lrho 
    - 2.50412783544030659*b7*lrho - 6.09708063744189756*b3*lrho2 
    + 16.0087447100660202*b4*lrho2 - 12.77595728730197071*b5*lrho2 
    + 2.8643356372167124*b6*lrho2;
  return f;
}

double FixedOrder::f22nl1_FIT_gq(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lbe2 = pow(lbe,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  double b6 = pow(b,6);
  double b7 = pow(b,7);
  
  double f = 0.00282450625811508503*b3 + 0.0001125529801773163916*b4 
    - 0.0325005325377584271*b5 + 0.171745428117885934*b6 
    - 0.1362549485229853436*b7 + 0.0000160504261251611439*b4*lbe 
    + 0.00258676966103027447*b5*lbe + 0.0358936097435228275*b6*lbe 
    + 0.0721640494225125838*b7*lbe - 0.0404568351507943864*b3*lrho 
    + 0.1116658960293489047*b4*lrho - 0.1145615544900474319*b5*lrho 
    + 0.0537148479428972514*b6*lrho - 0.01036242513667345341*b7*lrho 
    - 0.0253066537137784997*b3*lrho2 + 0.0551710547372787562*b4*lrho2 
    - 0.0403074509370066612*b5*lrho2 + 0.01044304642981058897*b6*lrho2;
  return f;
}

/*******************************************************************************	
*									       *	
*	 Exact NNLO results: qq, qqprime and qqbarprime reactions              *
*									       *	
*******************************************************************************/	

double FixedOrder::PXS_NNLO_qq(double rho)
{
  double LOG = lg+lgfr;
  double LOG2 = pow(LOG,2);
  
  double f = f20_FIT_qq(rho) 
    + LOG2 * f22_FIT_qqprime(rho) 
    + LOG  * f21_FIT_qq(rho);
  return f * pow(as,2);
}

double FixedOrder::PXS_NNLO_qqprime(double rho)
{
  double LOG = lg+lgfr;
  double LOG2 = pow(LOG,2);
  
  double f = f20_FIT_qqprime(rho)
    + LOG2 * f22_FIT_qqprime(rho) 
    + LOG  * f21_FIT_qqprime(rho);
  return f * pow(as,2);
}

double FixedOrder::PXS_NNLO_qqbarprime(double rho)
{
  double LOG = lg+lgfr;
  double LOG2 = pow(LOG,2);
  
  double f = f20_FIT_qqbarprime(rho)
    + LOG2 * f22_FIT_qqprime(rho) 
    + LOG  * f21_FIT_qqprime(rho);
  return f * pow(as,2);
}

//=======================================================================
// The NNLO functions (no scale-logs) for qq-, qq'- and qqbar'-reactions
//=======================================================================

double FixedOrder::f20_FIT_qq(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  
  double fexp = -0.740558 - 22.8129*pow(b,2) - 0.191648*pow(b,3) 
    - 6.58031*pow(b,4) - 0.537669*pow(b,5) + 31.7872*pow(b,6) 
    - 21.0783*lrho*rho - 3.25313*pow(lrho,2)*rho - 10.8176*lrho*pow(rho,2) 
    + 15.8988*pow(lrho,2)*pow(rho,2) + 8.64557*lrho*pow(rho,3);
  return -0.4768323995789214*lrho - pow(b,2)*exp(fexp);	
}

double FixedOrder::f20_FIT_qqprime(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  
  double fexp = -0.740558 - 23.4518*pow(b,2) - 0.193073*pow(b,3) 
    - 5.97215*pow(b,4) - 0.541402*pow(b,5) + 31.8227*pow(b,6) 
    - 21.3725*lrho*rho - 3.29162*pow(lrho,2)*rho - 11.1642*lrho*pow(rho,2) 
    + 15.9932*pow(lrho,2)*pow(rho,2) + 8.64746*lrho*pow(rho,3);
  return -0.4768323995789214*lrho - pow(b,2)*exp(fexp);
}

double FixedOrder::f20_FIT_qqbarprime(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  
  double fexp = -0.740572 - 31.2117*pow(b,2) - 0.31495*pow(b,3) 
    + 15.8601*pow(b,4) - 1.64639*pow(b,5) + 18.9767*pow(b,6) 
    - 19.6977*lrho*rho - 3.16565*pow(lrho,2)*rho - 16.1386*lrho*pow(rho,2) 
    + 12.3828*pow(lrho,2)*pow(rho,2) + 4.17707*lrho*pow(rho,3);
  return -0.4768323995789214*lrho - pow(b,2)*exp(fexp);
}

//==========================================================================
// The scaling functions ~Log(muF^2/m^2) for qq-, qq'- and qqbar'-reactions
//==========================================================================

double FixedOrder::f22_FIT_qqprime(double rho)// same for all 3 reactions
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  
  double fexp = -2.20937 - 99.5512*pow(b,2) + 0.394393*pow(b,3) 
    + 92.3169*pow(b,4) + 3.9297*pow(b,5) + 4.81762*pow(b,8) 
    - 43.9374*lrho*rho - 6.31738*pow(lrho,2)*rho - 48.7898*lrho*pow(rho,2) 
    + 25.8588*pow(lrho,2)*pow(rho,2) - 7.30711*lrho*pow(rho,3);	
  double f = -0.157185*pow(b,2) - 0.047419*lrho + pow(b,2)*pow(rho,0.72854)*exp(fexp);
  return f;
}

double FixedOrder::f21_FIT_qqprime(double rho)//same for qq' and qqbar'
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  
  double fexp = -0.07924228 + 211.8104*pow(b,2) - 0.4464898*pow(b,3) - 196.6873*pow(b,4) 
    - 5.141548*pow(b,5) - 9.448772*pow(b,8) + 90.67104*lrho*rho + 11.02778*pow(lrho,2)*rho 
    + 104.2288*lrho*pow(rho,2) - 55.31094*pow(lrho,2)*pow(rho,2) + 16.62298*lrho*pow(rho,3);
  double f = 1.17852*pow(b,2) + 0.254683*lrho - 1.*pow(b,2)*pow(rho,0.4)*exp(fexp);
  return f;
}

double FixedOrder::f21_FIT_qq(double rho)//for qq only
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  
  double fexp = -0.07923334 + 209.7125*pow(b,2) - 0.3786914*pow(b,3) - 195.7483*pow(b,4) 
    - 4.500831*pow(b,5) - 9.000056*pow(b,8) + 89.56873*lrho*rho + 10.92725*pow(lrho,2)*rho 
    + 103.6116*lrho*pow(rho,2) - 54.04722*pow(lrho,2)*pow(rho,2) + 16.25028*lrho*pow(rho,3);
  double f = 1.17852*pow(b,2) + 0.254683*lrho - 1.*pow(b,2)*pow(rho,0.4)*exp(fexp);
  return f;
}

/*******************************************************************************	
*									       *	
*     NNLO Partonic X-sections	for gg-reaction. It is a sum of 4 terms:       *
*       - exact_threshold term (without scale dependence)                      *
*       - fit for the power suppressed terms (including the "constant" term)   *
*       - Full NNLO scale dependence (not split into threshold and powers)     *
*       - scale dependent contributions from muF=/=muR originating at LO/NLO   *
*									       *	
*******************************************************************************/	

double FixedOrder::PXS_NNLO_EXACT_gg(double rho)
{
  double f = 
    + PXS_NNLO_THRESHOLD_gg(rho) 
    + PXS_NNLO_FIT_gg(rho)
    + PXS_NNLO_LO_NLOscales_gg(rho)
    + PXS_NNLO_scale_terms_gg(rho);
  return f;
}

//=========================================================================
//   threshold expansion of the NNLO function f20 
//   (with exact BORN factored out) from Beneke et al (2009); C2gg=0. 
//   The constant C2gg is implicitly contained in the fits below
//   (but note it also receives contrubution from: Born x 1/beta^2)
//=========================================================================

double FixedOrder::PXS_NNLO_THRESHOLD_gg(double rho)
{
  double nl = 5.0;// number of light flavors
  double b = sqrt(1-rho);
  double b2 = pow(b,2);
  double lbe = log(b);
  double lbe2 = pow(lbe,2);
  double lbe3 = pow(lbe,3);
  double lbe4 = pow(lbe,4);
  
  double thresh20nl0 = -0.024007157914810687/b + 0.43407982319605976/b2 
    + 14.861847493860669*lbe + (1.8153673293788755*lbe)/b 
    - 1.9983842014515505*lbe2 + (3.142857142857143*lbe2)/b 
    - 14.701583665657568*lbe3 + 29.18050088899328*lbe4;
  
  double thresh20nl1 = -0.0061192368274098/b + 0.13912366816397184*lbe 
    + (0.04365079365079365*lbe)/b - 0.7558262174090832*lbe2 
    + 0.5403796460924681*lbe3;
  
  double f = thresh20nl0 + nl * thresh20nl1;
  return pow(as,2) * f0gg(rho) * f;
}

//=========================================================================
//   Scale dependent contributions for muF=/=muR originating at LO and NLO
//=========================================================================

double FixedOrder::PXS_NNLO_LO_NLOscales_gg(double rho)
{
  double nl=5.0;
  double nl2=pow(nl,2);
  double lgfr2=pow(lgfr,2);
  
  // the functions fij are normalized as: as^n/m^2
  double f10 = 4*pi*f1gg(rho);
  double f11 = 4*pi*f1ggbar(rho);
  double f00 = f0gg(rho);
  
  double f = 
    f10*(-2.626056561016273*lgfr + 0.15915494309189535*lgfr*nl) 
    + f11*(-2.626056561016273*lg*lgfr - 2.626056561016273*lgfr2 
	   + 0.15915494309189535*lg*lgfr*nl + 0.15915494309189535*lgfr2*nl) 
    + f00*(-1.2918450914398067*lgfr + 2.298724353885538*lgfr2 
	   + 0.16042520743370148*lgfr*nl - 0.2786332550164289*lgfr2*nl 
	   + 0.008443431970194815*lgfr2*nl2);
  return f * pow(as,2);
}

//==========================================================================
// The scaling functions ~Log(muF^2/m^2) for gg-reaction.
//
// Some scaling functions from ``Hathor" (Aliev et al. arXiv:1007.1327v1)
//==========================================================================

double FixedOrder::PXS_NNLO_scale_terms_gg(double rho)
{
  double nl=5.0;
  double nl2=pow(nl,2);
  double LOG = lg+lgfr;
  double LOG2 = pow(LOG,2);
  
  double f22nl0 = f22nl0_FIT_gg(rho);
  double f22nl1 = f22nl1_FIT_gg(rho);
  double f21nl0 = f21nl0_FIT_gg(rho) + 20 * f21nl2_FIT_gg(rho);
  double f21nl1 = f21nl1_FIT_gg(rho) -  9 * f21nl2_FIT_gg(rho);
  double f21nl2 = f21nl2_FIT_gg(rho);
  
  double f = 
    + LOG *  (f21nl0 + nl * f21nl1 + nl2 * f21nl2)
    + LOG2 * (f22nl0 + nl * f22nl1);
  return f * pow(as,2);
}

double FixedOrder::f22nl1_FIT_gg(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lbe2 = pow(lbe,2);
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  
  double f = -0.00356104167821861222*b - 0.600682007963103239*b2 
    + 13.82970988281562116*b3 - 16.9635082827441165*b4 
    + 3.76261923077192499*b5 + 0.0116050479337840349*b*lbe 
    - 0.1575348075472875601*b2*lbe + 4.63151511439960648*b3*lbe 
    + 7.74747644911347943*b4*lbe - 2.16931397003164326*b5*lbe 
    - 0.01103205129641791145*b2*lbe2 + 1.039657934376285962*b3*lbe2 
    + 3.29932293902825979*b*lrho - 10.34880689568112313*b2*lrho 
    + 10.35776233604043958*b3*lrho - 3.30838403677087152*b4*lrho 
    + 0.782326666976654504*b*lrho2 - 1.65278029222162306*b2*lrho2 
    + 0.870448931607740223*b3*lrho2;
  return f;
}

double FixedOrder::f22nl0_FIT_gg(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lbe2 = pow(lbe,2);
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  
  double f = -0.0578389*b + 1.93094094141429693*b2 
    - 122.9455571577786397*b3 + 211.53435937689502*b4 
    - 93.3402588876750071*b5 - 0.704273292570916736*b*lbe 
    + 0.392152298237381979*b2*lbe - 37.6264268090281144*b3*lbe 
    - 61.3793218095110394*b4*lbe + 44.7561338172871625*b5*lbe 
    + 0.835563451232450513*b*lbe2 + 0.0171747976559594831*b2*lbe2 
    - 5.34358073726912582*b3*lbe2 - 45.9071083443894262*b*lrho 
    + 156.2229978854608434*b2*lrho - 167.54059459682459*b3*lrho 
    + 56.7473196239219973*b4*lrho - 18.1184678733560682*b*lrho2 
    + 39.3190338078504655*b2*lrho2 - 21.2004447764154102*b3*lrho2;
  return f;
}

double FixedOrder::f21nl1_FIT_gg(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lbe2 = pow(lbe,2);
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  
  double f = -0.00249982 - 0.0272733*b 
    + 19.4340755476554982*b2 - 435.539809500685645*b3 
    + 504.912020988537479*b4 - 88.799664693707602*b5 
    + 0.0432851214688035388*b*lbe + 5.11647907859189718*b2*lbe 
    - 146.1361340546396188*b3*lbe - 248.168773880287836*b4*lbe 
    + 55.4431105669571417*b5*lbe - 0.0464201917351361396*b*lbe2 
    + 0.35994235337967634*b2*lbe2 - 33.2094866285480819*b3*lbe2 
    - 100.5727014630313043*b*lrho + 302.705056562542737*b2*lrho 
    - 293.963282111113411*b3*lrho + 91.8347828939951401*b4*lrho 
    - 20.4652573057790352*b*lrho2 + 42.250982499345409*b2*lrho2 
    - 21.7855563762577368*b3*lrho2;
  return f;
}

double FixedOrder::f21nl0_FIT_gg(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double lbe2 = pow(lbe,2);
  double lbe3 = pow(lbe,3);
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  
  double f = 0.041247 + 0.0123377*b 
    - 661.55005204355904*b2 + 13004.29562244225291*b3 
    - 13876.36686145908297*b4 + 1547.963484189548846*b5 
    - 0.179987*lbe + 0.0505162*b*lbe - 177.28738825731523*b2*lbe 
    + 4661.3434534046446*b3*lbe + 7636.87990302021721*b4*lbe 
    - 1115.28388415014191*b5*lbe + 2.12283040860239304*b*lbe2 
    - 12.67292918310723623*b2*lbe2 + 1107.754574260264376*b3*lbe2 
    - 3.34225380492980205*b*lbe3 + 2369.63380918951308*b*lrho 
    - 7425.35627047328624*b2*lrho + 7426.27018972100489*b3*lrho 
    - 2368.07297895804527*b4*lrho + 606.613428260526799*b*lrho2 
    - 1272.126282331000418*b2*lrho2 + 665.508296075487579*b3*lrho2;
  return f;
}

double FixedOrder::f21nl2_FIT_gg(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lbe = log(b);
  double lrho2 = pow(lrho,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  
  double f = 0.00005344*b4*lbe + 0.00010271*b5*lbe 
    + 0.000206844*b2*rho - 0.00013099*b3*rho 
    + 2.21e-6*b5*rho + 0.000207078*lrho*rho2;
  return f;
}

//===========================================================================
//    The power suppressed terms (no scales). 
//    Includes the gg "constant" term.
//===========================================================================

double FixedOrder::PXS_NNLO_FIT_gg(double rho)
{
  double nl=5.0;
  double nl2 = pow(nl,2);
  // the functions fij are normalized in as^n/m^2
  double f20 = 
    + nl2 * f20nl2_FIT_gg(rho) 
    + nl  * f20nl1_FIT_gg(rho) 
    +       f20nl0_FIT_gg(rho);
  return f20 * pow(as,2);
}

double FixedOrder::f20nl2_FIT_gg(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double lrho2 = pow(lrho,2);
  double lbe = log(b);
  double rho2 = pow(rho,2);	
  double rho3 = pow(rho,3);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  
  double f = (6.44022*b*rho - 4.8664*b2*rho 
	      - 0.0324653*lrho2*rho - 13.8424*b*rho2 
	      + 4.7366*b2*rho2 - 2.91398*lrho*rho2 
	      + 8.43828*b*rho3 - 2.78748*b2*rho3 
	      + 2.38971*b3*rho3)/10000.0 ;
  return f;
}

double FixedOrder::f20nl1_FIT_gg(double rho)
{
  double b = sqrt(1-rho);
  double lbe = log(b);
  double lbe2 = pow(lbe,2);
  double lrho = log(rho);
  double lrho2 = pow(lrho,2);
  double lrho3 = pow(lrho,3);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b5 = pow(b,5);
  double b6 = pow(b,6);
  double b7 = pow(b,7);
  
  double f = -0.0195046*b - 1.4717*b2 
    - 0.223616*b3 + 0.499196*b5 
    + 1.32756*b7 + 0.00466872*b3*lbe 
    + 0.0321469*b6*lbe2 + 0.579781*lrho2*rho 
    + 0.166646*lrho3*rho - 1.36644*lrho*rho2 
    + 2.24909*lrho2*rho2 ;
  
  return f;
}

double FixedOrder::f20nl0_FIT_gg(double rho)
{
  double b = sqrt(1-rho);
  double lbe = log(b);
  double lrho = log(rho);
  double lrho2 = pow(lrho,2);
  double rho2 = pow(rho,2);	
  double rho3 = pow(rho,3);	
  double rho4 = pow(rho,4);	
  double rho5 = pow(rho,5);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  
  double f = 581.27542*b + 1251.4057*b2 
    - 60.478096*b3 + 1101.2272*b4 
    - 2905.3858*b5 + 629.9128*b4*lbe 
    - 5.189107503169016*lrho 
    + 1200.741*lrho*rho + 162.50333*lrho2*rho 
    - 1810.2849*b*rho2 + 36.074524*lrho*rho2 
    - 1192.8918*lrho2*rho2 + 1568.7591*b*rho3 
    - 461.21326*b*rho4 + 121.6379*b*rho5 ;
  
  return f;
}

/*******************************************************************************	
*									       *	
*     NNLO Partonic X-sections	for qqbar-reaction (sum of 4 terms)	       *
*       - exact_threshold term (with its own scale dependence)                 *
*       - fit for the power suppressed terms (including the "constant" term)   *
*       - the contribution with ID final state not icluded above               * 
*         (it is a sum of "pure interference" + qqbar' -> ttbar + qqbar')      *
*       - scale dependent contributions from muF=/=muR originating at LO/NLO   *
*     NOTE: the term ~nl^2 has no threshold, and is known analytically         *
*									       *	
*******************************************************************************/	

double FixedOrder::PXS_NNLO_EXACT_qqbar(double rho)
{
  double f = 
    + PXS_NNLO_THRESHOLD_qqbar(rho) 
    + PXS_NNLO_FIT_qqbar(rho)
    + PXS_NNLO_LO_NLOscales_qqbar(rho)
    + PXS_NNLO_identical_qqbar(rho);
  return f;
}

//================================================================================
// The part of sigma(q qbar -> t tbar q qbar) not included in PXS_NNLO_FIT_qqbar.
// It is a sum of a pure interference term and sigma(q qbar' -> t tbar q qbar').
//================================================================================

double FixedOrder::PXS_NNLO_identical_qqbar(double rho)
{
  double b = sqrt(1-rho);
  double lbe = log(b);
  double lrho = log(rho);
  // the fit of the pure interference contribution:
  double f = -1.26644*pow(b,20)*lbe + 1.53647*pow(b,3)*rho 
    + 10.7411*pow(b,4)*rho + 3.86236*lrho*rho + 0.327143*pow(lrho,2)*rho 
    - 24.3298*pow(b,4)*pow(rho,2) - 21.332*lrho*pow(rho,2) 
    - 10.7669*pow(lrho,2)*pow(rho,2) - 4.50719*pow(b,3)*pow(rho,3) 
    + 15.4975*pow(b,4)*pow(rho,3) + 17.4705*lrho*pow(rho,3) 
    + 2.90068*pow(b,3)*pow(rho,4) - 4.98808*pow(b,4)*pow(rho,4);
  
  f += f20_FIT_qqbarprime(rho);	 
  return f * pow(as,2);
}

//=========================================================================
//   threshold expansion of the NNLO functions f22, f21 and f20 
//   (with exact BORN factored out) from Beneke et al (2009); C2qq=0. 
//   The actual constant is set to zero here and is contained in the fits.
//=========================================================================

double FixedOrder::PXS_NNLO_THRESHOLD_qqbar(double rho)
{
  double nl = 5.0;// number of light flavors
  double b = sqrt(1-rho);
  double lbe = log(b);
  double lbe2 = pow(lbe,2);
  double lbe3 = pow(lbe,3);
  double lbe4 = pow(lbe,4);
  
  double thresh20nl0 = 0.022846306484003143/pow(b,2) 
    - 0.033390549389716966/b + 1.5810883840639216*lbe 
    + (0.3422025061689375*lbe)/b + 6.62693178820539*lbe2 
    - (0.8888888888888888*lbe2)/b - 9.531529857990602*lbe3 
    + 5.76404955831966*lbe4;
  double thresh20nl1 = 0.01168217939778234/b + 0.35320738606771657*lbe 
    - (0.027777777777777776*lbe)/b - 0.5752392118254954*lbe2 + 0.2401687315966525*lbe3;
  double thresh21nl0 = 2.963036815826603 - 0.5208333333333334/b - 5.322595839415953*lbe 
    + (0.4444444444444444*lbe)/b + 11.307833327841358*lbe2 - 5.76404955831966*lbe3;
  double thresh21nl1 = -0.40917468164069626 + 0.041666666666666664/b 
    + 0.46164291174543004*lbe - 0.5403796460924681*lbe2;
  double thresh22nl0 = 1.7154768057697534 - 3.5187082038820767*lbe + 1.441012389579915*lbe2;
  double thresh22nl1 = -0.26328935946926935 + 0.22515818587186173*lbe;
  
  double f20 = thresh20nl0 + nl*thresh20nl1;
  double f21 = thresh21nl0 + nl*thresh21nl1;
  double f22 = thresh22nl0 + nl*thresh22nl1;
  
  double f = f20 + f21*(lg + lgfr) + f22*(pow(lg,2) + 2.*lg*lgfr + pow(lgfr,2));
  return pow(as,2) * f0qq(rho) * f;
}

//=========================================================================
//   Scale dependent contributions for muF=/=muR originating at LO and NLO
//=========================================================================

double FixedOrder::PXS_NNLO_LO_NLOscales_qqbar(double rho)
{
  double nl=5.0;
  // the functions fij are normalized as: as^n/m^2
  double f10 = 4*pi*f1qq(rho);
  double f11 = 4*pi*f1qqbar(rho);
  double f00 = f0qq(rho);
  
  double f =
    + f10*(-2.626056561016273*lgfr + 0.15915494309189535*lgfr*nl) 
    + f11*(-2.626056561016273*lg*lgfr - 2.626056561016273*pow(lgfr,2) 
	   + 0.15915494309189535*lg*lgfr*nl + 0.15915494309189535*pow(lgfr,2)*nl) 
    + f00*(-1.2918450914398067*lgfr + 2.298724353885538*pow(lgfr,2) + 0.16042520743370148*lgfr*nl 
	   - 0.2786332550164289*pow(lgfr,2)*nl + 0.008443431970194815*pow(lgfr,2)*pow(nl,2));
  return f * pow(as,2);
}

//===========================================================================
//   Fit for the power suppressed terms (contains the qqbar "constant" term)
//   NOTE: the term ~nl^2 has no "threshold", and is known analytically
//===========================================================================

double FixedOrder::PXS_NNLO_FIT_qqbar(double rho)
{
  double nl=5.0;
  // the functions fij are normalized in as^n/m^2
  double f20 = pow(nl,2)*f20nl2_FIT_qqbar(rho) + nl*f20nl1_FIT_qqbar(rho) + f20nl0_FIT_qqbar(rho);
  double f21 = pow(nl,2)*f21nl2_FIT_qqbar(rho) + nl*f21nl1_FIT_qqbar(rho) + f21nl0_FIT_qqbar(rho);
  double f22 = pow(nl,2)*f22nl2_FIT_qqbar(rho) + nl*f22nl1_FIT_qqbar(rho) + f22nl0_FIT_qqbar(rho);
  
  double f = f20 + f21*(lg + lgfr) + f22*(pow(lg,2) + 2.*lg*lgfr + pow(lgfr,2));
  return f * pow(as,2);
}

//==========================================================================
//   The fits for qqbar -> ttbar:
//==========================================================================

double FixedOrder::f22nl2_FIT_qqbar(double rho)
{
  return f0qq(rho)/(12.*pow(pi,2));
}

double FixedOrder::f22nl1_FIT_qqbar(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double f = 0.00239658*b*rho - 0.0220604*pow(b,2)*rho 
    + 0.0188993*pow(b,3)*rho - 0.0131116*lrho*rho 
    - 0.000963099*b*pow(rho,2) + 0.0141569*pow(b,2)*pow(rho,2) 
    - 0.00104753*pow(b,3)*pow(rho,2) - 0.000699993*lrho*pow(rho,2) 
    - 0.0014341*b*pow(rho,3) + 0.00610629*lrho*pow(rho,3) 
    - 0.000214179*lrho*pow(rho,4);
  return f;
}

double FixedOrder::f22nl0_FIT_qqbar(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double f = 0.191165*pow(b,2) - 13.2537*b*rho 
    + 195.393*pow(b,2)*rho - 139.565*pow(b,3)*rho + 10.4168*lrho*rho 
    + 16.2786*b*pow(rho,2) + 165.851*pow(b,2)*pow(rho,2) 
    + 235.88*pow(b,3)*pow(rho,2) + 213.556*lrho*pow(rho,2) 
    - 3.03232*b*pow(rho,3) + 835.447*pow(b,2)*pow(rho,3) 
    - 90.1477*pow(b,3)*pow(rho,3) + 677.134*lrho*pow(rho,3) 
    + 295.537*lrho*pow(rho,4);
  return f;
}

double FixedOrder::f21nl2_FIT_qqbar(double rho)
{
  double lrho = log(rho);
  return f0qq(rho) * (5 + 3*lrho - 6*log(2))/(18.*pow(pi,2));
}

double FixedOrder::f21nl1_FIT_qqbar(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double f = - 0.000188644*pow(b,2) + 0.371362*b*rho 
    - 3.14887*pow(b,2)*rho + 3.50366*pow(b,3)*rho 
    + 0.102822*lrho*rho - 0.328962*b*pow(rho,2) 
    + 8.58138*pow(b,2)*pow(rho,2) - 3.93341*pow(b,3)*pow(rho,2) 
    + 1.55362*lrho*pow(rho,2) - 0.0428876*b*pow(rho,3) 
    - 1.31531*pow(b,2)*pow(rho,3) + 0.896763*pow(b,3)*pow(rho,3) 
    + 2.39531*lrho*pow(rho,3) + 0.0208085*lrho*pow(rho,4);
  return f;
}

double FixedOrder::f21nl0_FIT_qqbar(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double f = -0.720681*pow(b,2) + 46.5167*b*rho - 761.601*pow(b,2)*rho 
    + 503.018*pow(b,3)*rho - 50.2996*lrho*rho - 61.2529*b*pow(rho,2) 
    - 1148.19*pow(b,2)*pow(rho,2) - 914.708*pow(b,3)*pow(rho,2) 
    - 1046.69*lrho*pow(rho,2) + 14.8024*b*pow(rho,3) 
    - 3847.6*pow(b,2)*pow(rho,3) + 377.917*pow(b,3)*pow(rho,3) 
    - 3269.96*lrho*pow(rho,3) - 1386.84*lrho*pow(rho,4);
  return f;
}

double FixedOrder::f20nl2_FIT_qqbar(double rho)
{
  double lrho = log(rho);
  double f = (25 - 3*pow(pi,2) + 30*(lrho - 2*log(2)) 
	      + 9*pow(lrho - 2*log(2),2))/(108.*pow(pi,2));
  return f * f0qq(rho);	
}

double FixedOrder::f20nl1_FIT_qqbar(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double f = 0.90756*b*rho - 6.75556*pow(b,2)*rho + 9.18183*pow(b,3)*rho 
    + 1.3894*lrho*rho + 0.137881*pow(lrho,2)*rho - 0.99749*b*pow(rho,2) 
    + 27.7454*pow(b,2)*pow(rho,2) - 12.9055*pow(b,3)*pow(rho,2) 
    + 6.13693*lrho*pow(rho,2) - 0.0077383*b*pow(rho,3) 
    - 4.49375*pow(b,2)*pow(rho,3) + 3.86854*pow(b,3)*pow(rho,3) 
    + 8.78276*lrho*pow(rho,3) - 0.380386*pow(b,4)*pow(rho,4) 
    - 0.0504095*lrho*pow(rho,4);
  return f;
}

double FixedOrder::f20nl0_FIT_qqbar(double rho)
{
  double b = sqrt(1-rho);
  double lrho = log(rho);
  double f = -2.32235*b*rho + 44.3522*pow(b,2)*rho - 24.6747*pow(b,3)*rho 
    + 3.11129*lrho*rho + 2.92101*b*pow(rho,2) + 224.311*pow(b,2)*pow(rho,2) 
    + 21.5307*pow(b,3)*pow(rho,2) + 100.125*lrho*pow(rho,2) 
    + 2.05531*b*pow(rho,3) + 945.506*pow(b,2)*pow(rho,3) 
    + 36.1101*pow(b,3)*pow(rho,3) - 176.632*pow(b,4)*pow(rho,3) 
    + 563.1*lrho*pow(rho,3) + 7.68918*pow(b,4)*pow(rho,4) 
    + 568.023*lrho*pow(rho,4);
  return f;
}

/*******************************************************************************	
*									       *	
*	                          LO results	                               *
*									       *	
*******************************************************************************/	

double FixedOrder::f0qq(double rho)
{
  double beta = sqrt(1-rho);
  return pi*beta*rho/27.*(2+rho);
}

double FixedOrder::f0qqLeading(double rho)//only the O(beta) term (for testing)
{
  double beta = sqrt(1-rho);
  return 0.111111*beta*pi;
}

double FixedOrder::f0gg(double rho)
{
  double beta = sqrt(1-rho);
  double f = pi*beta*(rho/192.)
    *((pow(rho,2) + 16*rho+16)
      *log((1+beta)/(1-beta))/beta-28-31*rho);	
  return f;
}

double FixedOrder::f0ggLeading(double rho)//only the O(beta) term (for testing)
{
  double beta = sqrt(1-rho);
  return 0.0364583*beta*pi;
}

/*******************************************************************************	
*									       *	
*                                 NLO results	                               *
*									       *	
*******************************************************************************/	

//==========================================================================
//   Scale-dependent terms (as in: Nason, Dawson, Ellis '88)
//==========================================================================

double FixedOrder::h1(double beta)
{
  double f = -2*gsl_sf_dilog((1 - beta)/2.) 
    + 2*gsl_sf_dilog((1 + beta)/2.) 
    - pow(log((1 - beta)/2.),2) 
    + pow(log((1 + beta)/2.),2);
  return f;
}

double FixedOrder::h2(double beta)
{
  double f = -gsl_sf_dilog((-2*beta)/(1 - beta)) 
    + gsl_sf_dilog((2*beta)/(1 + beta));
  return f;
}

double FixedOrder::f1qqbar(double rho)
{
  double beta = sqrt(1-rho);
  double nlf = 5;
  double f = (2*rho*log((1 + beta)/(1 - beta)))/(81.*pi) 
    + f0qq(rho)*(-(-127 + 6*nlf)/(72.*pow(pi,2)) 
		 + (2*log(rho/(4.*pow(beta,2))))/(3.*pow(pi,2)));
  return f;
}

double FixedOrder::f1ggbar(double rho)
{
  double beta = sqrt(1-rho);	
  double f = (-181*beta)/(1440.*pi) 
    + (26*beta*rho)/(45.*pi) - (2483*beta*pow(rho,2))/(1920.*pi) 
    + (-rho/(8.*pi) + pow(rho,2)/(16.*pi) - pow(rho,3)/(256.*pi))*h1(beta) 
    + (rho/(8.*pi) + pow(rho,2)/(8.*pi) + pow(rho,3)/(128.*pi))*h2(beta) 
    + ((-3*rho)/(8.*pi) + (33*pow(rho,2))/(128.*pi) 
       + (59*pow(rho,3))/(768.*pi))*log((1 + beta)/(1 - beta)) 
    + (3*f0gg(rho)*log(rho/(4.*pow(beta,2))))/(2.*pow(pi,2));
  return f;
}

double FixedOrder::f1gqbar(double rho)
{
  double beta = sqrt(1-rho);	
  double f = (-181*beta)/(6480.*pi) + (289*beta*rho)/(2160.*pi) 
    - (1319*beta*pow(rho,2))/(25920.*pi) 
    + (-rho/(72.*pi) + pow(rho,2)/(144.*pi))*h1(beta) 
    + ((-17*rho)/(432.*pi) + pow(rho,2)/(128.*pi) 
       + (7*pow(rho,3))/(1728.*pi))*log((1 + beta)/(1 - beta));
  return f;
}

//==========================================================================
//   Scale-independent terms at NLO (starting from ver2.0 use CM'08)
//==========================================================================

double FixedOrder::f1qq(double rho)
{
  return f1qqCM(rho);
}

double FixedOrder::f1gg(double rho)
{
  return f1ggCM(rho);
}

double FixedOrder::f1gq(double rho)
{
  return f1gqCM(rho);
}

//==========================================================================
//   High-quality parameterization exact NLO results (Czakon, Mitov '08)
//   Fits from HATHOR: Aliev et al arXiv:1007.1327
//==========================================================================

double FixedOrder::f1qqCM(double rho)
{
  double nl = 5.0;
  double b = sqrt(1-rho);
  double lbe = log(b);
  double lbe2 = pow(lbe,2);
  double lrho = log(rho);
  double lrho2 = pow(lrho,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  
  double f10nl0 = 
    + 0.894801363487060697*b2 - 15.9806533566333508*b3 
    + 15.5948076822314514*b4 - 0.508993720940668306*b5 
    + 0.25808023534387927*b2*lbe - 3.99149865104188887*b3*lbe 
    - 8.97737571174816865*b4*lbe + 0.147026822789739241*b5*lbe 
    + 0.0187135812730567129*b2*lbe2 - 1.81602862945435741*b3*lbe2 
    - 1.74752529849853127*b*lrho + 6.20751624259016067*b2*lrho 
    - 6.75340649699489205*b3*lrho + 2.29340856664909682*b4*lrho 
    + 0.135309206882346187*b*lrho2 - 0.0712992055892768895*b2*lrho2 
    - 0.0640103309259597209*b3*lrho2 - 0.0913852259360125798*rho 
    + 0.44171954644071765*b*rho - 0.57251372837945371*b*lbe*rho 
    + 1.185185185185185185*b*lbe2*rho + 0.0987654320987654321*b*lrho*rho 
    + 0.0138455459529272122*b*rho2 + 0.049382716049382716*b*lrho*rho2;
  
  double f10nl1 = (b*(3 - 4*b2 + b4)*(-5 - 3*lrho + log(64)))/243.;
  return (f10nl0 + nl*f10nl1)/(4.0*pi);//CM normalization: as |-> NDE's: 4*pi*as
}

double FixedOrder::f1ggCM(double rho)
{
  double nl = 5.0;
  double b = sqrt(1-rho);
  double lbe = log(b);
  double lbe2 = pow(lbe,2);
  double lrho = log(rho);
  double lrho2 = pow(lrho,2);
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  double logxconf = log((1.0-b)/(1.0+b));
  
  double f10nl0 = 
    + 0.0471205071232564865 + 0.235655553965240606*b 
    - 112.1628025020790127*b2 + 1883.77093917245376*b3 
    - 1766.27896687670163*b4 - 4.28709466299521533*b5 
    - 0.086738651030143563*b*lbe - 30.2912153871041909*b2*lbe 
    + 687.805453265266195*b3*lbe + 1142.477618105764338*b4*lbe 
    - 61.3742807324254932*b5*lbe + 0.875*b*lbe2 
    - 2.19494021495237345*b2*lbe2 + 169.273238020325875*b3*lbe2 
    + 284.814617211531812*b*lrho - 849.780999197478008*b2*lrho 
    + 817.825794348102379*b3*lrho - 252.863693983687393*b4*lrho 
    + 57.8966224150269576*b*lrho2 - 121.9429831242445605*b2*lrho2 
    + 64.0461741451940366*b3*lrho2 + 0.03125*b*rho2 
    + 0.015625*logxconf*rho2;
  
  double f10nl1 = -((2*b + logxconf)*rho2)/256.;
  return (f10nl0 + nl*f10nl1)/(4.0*pi);//CM normalization: as |-> NDE's: 4*pi*as
}	

double FixedOrder::f1gqCM(double rho)
{
  double b = sqrt(1-rho);
  double lbe = log(b);
  double lbe2 = pow(lbe,2);
  double lrho = log(rho);
  double lrho2 = pow(lrho,2);
  double lrho3 = lrho*lrho2;
  double lrho4 = pow(lrho2,2);
  double lrho5 = lrho*lrho4;
  double rho2 = pow(rho,2);	
  double b2 = pow(b,2);
  double b3 = pow(b,3);
  double b4 = pow(b,4);
  double b5 = pow(b,5);
  double b6 = pow(b,6);
  
  double f10nl0 = 
    -0.02457581886482620869*b3 - 3.28032159897403943*b4 
    + 3.7941231153033735*b5 - 0.189185108551928629*b6 
    + 0.1388888888888888889*b3*lbe - 0.0178630821547912901*b4*lbe 
    - 0.585680726926105217*b5*lbe - 1.8961443995961806*b6*lbe 
    - 3.19157673250705547*b2*lrho*rho + 5.01130734325825868*b3*lrho*rho 
    - 1.78170350247765988*b4*lrho*rho - 0.1255542812553638502*b2*lrho2*rho 
    - 0.307143757144090147*b3*lrho2*rho + 0.23465015682704749*b4*lrho2*rho 
    + 0.0299903911910146112*b2*lrho3*rho - 0.000427110882824291123*b2*lrho4*rho 
    - 0.00001115993665476662179*b2*lrho5*rho;
  return f10nl0/(4.0*pi);//CM normalization: as |-> NDE's: 4*pi*as
}

//==========================================================================
//   NDE parameterization of NLO results (Nason, Dawson, Ellis '88)
//==========================================================================

double FixedOrder::f1qqNDE(double rho)
{
  double beta = sqrt(1-rho);
  double nlf = 5;
  double a0 = 0.180899,
    a1 = 0.101949,
    a2 = -0.234371,
    a3 = -0.0109950,
    a4 = -0.0185575,
    a5 = 0.00907527,
    a6 = 0.0160367,
    a7 = 0.00786727;
  
  double f = (a0*beta + a2*pow(beta,3) + a4*pow(beta,5) - pi/432.)*rho 
    + (a1*pow(beta,3) + a3*pow(beta,5) + a5*pow(beta,7) 
       - (41*beta)/(108.*pi))*rho*log(8*pow(beta,2)) 
    + (2*beta*rho*pow(log(8*pow(beta,2)),2))/(27.*pi) 
    + f0qq(rho)*((-5*(-4 + nlf))/(36.*pow(pi,2)) 
		 + ((-4 + nlf)*log(4/rho))/(12.*pow(pi,2))) 
    + a6*beta*rho*log(rho) + a7*beta*rho*pow(log(rho),2);
  return f;
}

double FixedOrder::f1ggNDE(double rho)
{
  double beta = sqrt(1-rho);
  double nlf = 5;	
  double a0 = 0.108068,
    a1 = -0.114997,
    a2 = 0.0428630,
    a3 = 0.131429,
    a4 = 0.0438768,
    a5 = -0.0760996,
    a6 = -0.165878,
    a7 = -0.158246;
  
  double f = a0*beta + a2*pow(beta,3) + (11*pi)/9216. 
    - (beta*(-4 + nlf)*pow(rho,2))/(512.*pi)
    + (a1*pow(beta,3) + a3*pow(beta,5) - (61*beta)/(256.*pi))*log(8*pow(beta,2)) 
    + (7*beta*pow(log(8*pow(beta,2)),2))/(128.*pi) 
    + ((-4 + nlf)*pow(rho,2)*log((1 + beta)/(1 - beta)))/(1024.*pi) 
    + (a6*beta*rho + a4*beta*pow(rho,2))*log(rho) 
    + (a7*beta*rho + a5*beta*pow(rho,2))*pow(log(rho),2);
  return f;
}

double FixedOrder::f1gqNDE(double rho)
{
  double beta = sqrt(1-rho);
  double a0 = 0.0110549,
    a1 = -0.426773,
    a2 = -0.00103876,
    a3 = 0.450626,
    a4 = -0.227229,
    a5 = 0.0472502,
    a6 = -0.197611,
    a7 = -0.0529130;
  
  double f = a1*pow(beta,3) + a3*pow(beta,5) 
    + (a0*pow(beta,3) + a2*pow(beta,5))*log(beta) 
    + (a6*beta*rho + a4*beta*pow(rho,2))*log(rho) 
    + (a7*beta*rho + a5*beta*pow(rho,2))*pow(log(rho),2);
  return f;
}
