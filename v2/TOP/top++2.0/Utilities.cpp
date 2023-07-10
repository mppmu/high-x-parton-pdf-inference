#include <algorithm>
#include "Utilities.h"

using std::string;
using std::vector;
using std::domain_error;
using std::max;
using std::cout;
using std::endl;
using std::max_element;

/*******************************************************************************	
*									       *	
*  Compute the final result when Fixed Order Calculation is requested          *
*									       *	
*******************************************************************************/

double ComputeFinalResultFO(PartonicFlux Flux,
			    FixedOrder FO, 
			    double rho, 
			    int precision,
			    string partonchannel)
{
  return ConvolutionFO(Flux,FO, rho, precision, partonchannel);
}

double f_FO (double x, 
	     void * par) 
{
  FOcomponents * p = (FOcomponents*) par;
  double f;
  
  if ((p->partonchannel)=="qqbar")
    {
      f = (p->Flux).fluxqqbar(x) * (p->FO).FOqqbar((p->rho)/x) ;
    }
  else if ((p->partonchannel)=="gg")
    {
      f = (p->Flux).fluxgg(x) * (p->FO).FOgg((p->rho)/x) ;
    }
  else if ((p->partonchannel)=="qg")
    {
      f = (p->Flux).fluxgq(x) * (p->FO).FOgq((p->rho)/x) ;
    }
  else if ((p->partonchannel)=="qq")
    {
      f = (p->Flux).fluxqq(x) * (p->FO).FOqq((p->rho)/x) ;
    }
  else if ((p->partonchannel)=="qqprime")
    {
      f = (p->Flux).fluxqqprime(x) * (p->FO).FOqqprime((p->rho)/x) ;
    }
  else if ((p->partonchannel)=="qqbarprime")
    {
      f = (p->Flux).fluxqqbarprime(x) * (p->FO).FOqqbarprime((p->rho)/x) ;
    }
  else//all channels combined
    {
      f = (p->Flux).fluxqqbar(x)    * (p->FO).FOqqbar((p->rho)/x)      +
	(p->Flux).fluxgg(x)         * (p->FO).FOgg((p->rho)/x)         +
	(p->Flux).fluxgq(x)         * (p->FO).FOgq((p->rho)/x)         +
	(p->Flux).fluxqq(x)         * (p->FO).FOqq((p->rho)/x)         +
	(p->Flux).fluxqqprime(x)    * (p->FO).FOqqprime((p->rho)/x)    +
	(p->Flux).fluxqqbarprime(x) * (p->FO).FOqqbarprime((p->rho)/x) ;
    }
  return f/x;
}

double ConvolutionFO(PartonicFlux Flux, 
		     FixedOrder FO, 
		     double rho, 
		     int precision,
		     string partonchannel)
{
  FOcomponents temp_p = {Flux, FO , rho, partonchannel};
  gsl_function F;
  F.function = f_FO;
  F.params = &temp_p;
  //------------------- setup for the integration routine (not meant to be modified):
  int prec = precision;
  double a = rho;
  double b = 1.0;
  double epsabs = 0.0;
  double epsrel = 1/pow(10,prec);
  size_t limit = 1000;
  int key = 6;
  //-------------------------------
  double result, abserr;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);	
  gsl_integration_qag (&F, a, b, epsabs, epsrel,limit, key, w, &result, &abserr); 
  gsl_integration_workspace_free (w);
  
  if (abserr/result > epsrel) 
    {
      cout 
	<< "Warning (Fixed Order cross-section convolution): relative error reached " 
	<< abserr/result 
	<< endl;
    }
  return result;
}

/*******************************************************************************	
*									       *	
*        Compute the final result for mixed F.O. with Resummation              *
*									       *	
*******************************************************************************/

double ComputeFinalResultRES(PartonicFlux Flux,
			     SubtrFlux SFlux, 
			     FixedOrder FO, 
			     Resummation Res, 
			     double rho, 
			     int precision,
			     string partonchannel)
{	
  double f1 = ConvolutionRES(Flux, SFlux, Res, rho, precision, partonchannel);
  double f2;
  if (partonchannel=="qqbar")
    {
      f2 = Res.resumqqbarSubFlux(rho) ;
    }
  else if (partonchannel=="gg")
    {
      f2 = Res.resumggSubFlux(rho) ;
    }
  else if (partonchannel=="qg" ||
	   partonchannel=="qqbarprime" ||
	   partonchannel=="qq" ||
	   partonchannel=="qqprime")
    {
      f2 = 0 ;
    }
  else//all channels combined
    {
      f2 = Res.resumqqbarSubFlux(rho) + Res.resumggSubFlux(rho);	
    }
  return f1+f2;
}

double ConvolutionRES(PartonicFlux Flux, 
		      SubtrFlux SFlux, 
		      Resummation Res, 
		      double rho, 
		      int precision,
		      string partonchannel)
{
  REScomponents temp_p = {Flux, SFlux, Res, rho, partonchannel};
  gsl_function F;
  F.function = f_RES;
  F.params = &temp_p;
  //------------------- setup for the integration routine (not meant to be modified):
  int prec = precision;
  double a = rho;
  double b = 1.0;
  double epsabs = 0.0;
  double epsrel = 1/pow(10,prec);
  size_t limit = 1000;
  int key = 6;
  //-------------------------------
  double result, abserr;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);	
  gsl_integration_qag (&F, a, b, epsabs, epsrel,limit, key, w, &result, &abserr); 
  gsl_integration_workspace_free (w);
  if (abserr/result > epsrel) 
    {
      cout 
	<< "Warning (Resummed cross-section x-space convolution): relative error reached " 
	<< abserr/result 
	<< endl;
    }
  return result;
}

double f_RES (double x, void * par) 
{
  REScomponents * p = (REScomponents*) par;	
  double f;
  
  if ((p->partonchannel)=="qqbar")
    {
      f = ((p->Flux).fluxqqbar(x) - (p->SFlux).subtrfluxqqbar(x) ) 
	* (p->Res).resumqqbar((p->rho)/x) ; 
    }
  else if ((p->partonchannel)=="gg")
    {
      f = ((p->Flux).fluxgg(x) - (p->SFlux).subtrfluxgg(x) )  
	* (p->Res).resumgg((p->rho)/x) ; 
    }
  else if ((p->partonchannel)=="qg" ||
	   (p->partonchannel)=="qqbarprime" ||
	   (p->partonchannel)=="qq" ||
	   (p->partonchannel)=="qqprime")
    {
      f = 0 ;
    }
  else//all channels combined
    {
      f = ((p->Flux).fluxqqbar(x) - (p->SFlux).subtrfluxqqbar(x) ) 
	* (p->Res).resumqqbar((p->rho)/x) 
	+ ((p->Flux).fluxgg(x) - (p->SFlux).subtrfluxgg(x) )  
	* (p->Res).resumgg((p->rho)/x) ;
    }
  return f/x;
}

/*******************************************************************************	
*									       *	
*        If Fact. and Ren. scales pass the scale variations criterion.         *
*									       *	
*******************************************************************************/

bool ScalesPass(string RestrictedScaleVariation, 
		double F, 
		double R)
{// At this point RestrictedScaleVariation is either "NO" or a "double >= 1.0"
  if (RestrictedScaleVariation=="NO") return 1;
  else {
    double x = atof(RestrictedScaleVariation.c_str());
    if (F/R <= x && F/R >= 1/x) return 1;
    else return 0;
  }
}

/*******************************************************************************	
*									       *	
*        Converts the options OrderRES and OrderFO to integers:                *
*									       *	
*******************************************************************************/

int IntegerForm(string s)
{
  if      (s=="LO"   || s=="LL")   {return 1;}
  else if (s=="NLO"  || s=="NLL")  {return 2;}
  else if (s=="NNLO" || s=="NNLL") {return 3;}
  else 
    {
      cout 
	<< "Warning: unknown option for OrderRES or OrderFO (LO/LL,NLO/NLL,NNLO/NNLL expected)." 
	<< endl;
      return 0;
    }
}

/*******************************************************************************	
*									       *	
*               Print a greeting at the beginning of the program               *
*									       *	
*******************************************************************************/

void Greetings(string ver)// ver == program's version
{
  const string::size_type pad = 7;//number of L/R blanks in longest line
  
  //characters for the frame:
  char up = '=';
  char le = '=';
  char ri = '=';
  char dw = '=';
  
  string text0 = "Top++ version";
  string text1 = "Program for the calculation of the inclusive top pair";
  string text2 = "production cross-section at hadron colliders";
  string text3 = "Authors: Michal Czakon and Alexander Mitov";
  
  string line1 = le+string(pad, ' ')+text1+string(pad, ' ')+ri;
  string line0 = text0+' '+ver;
  const string::size_type L1 = line1.size();// the longest text
  const string::size_type left0 = (L1-2-line0.size())/2;
  const string::size_type right0 = L1-2-line0.size()-left0;
  const string::size_type left2 = (L1-2-text2.size())/2;
  const string::size_type right2 = L1-2-text2.size()-left2;
  const string::size_type left3 = (L1-2-text3.size())/2;
  const string::size_type right3 = L1-2-text3.size()-left3;
  
  cout 
    << string(L1, up)+"\n"
    << le+string(L1-2, ' ')+ri+"\n"
    << le+string(left0, ' ')+line0+string(right0, ' ')+ri+"\n"
    << le+string(L1-2, ' ')+ri+"\n"
    << line1+"\n" 
    << le+string(left2, ' ')+text2+string(right2, ' ')+ri+"\n"
    << le+string(L1-2, ' ')+ri+"\n"
    << le+string(left3, ' ')+text3+string(right3, ' ')+ri+"\n"
    << le+string(L1-2, ' ')+ri+"\n"
    << string(L1, dw)+"\n"
    << endl;
}

/*******************************************************************************	
*									       *	
*		 PDF uncertainties and related functionalities                 *
*									       *	
*******************************************************************************/	

//=======================================================================
//   1)    The asymmetric prescription used by MSTW and others.
//         P. Nadolsky and Z. Sullivan [arXiv:hep-ph/0110378], p. P510
//         (see Cacciari et al. http://arxiv.org/abs/0804.2800)  
//=======================================================================

void deltaPDFasymmetric(const vector<double>& result, 
			double& smax, 
			double& smin)
{
  double central = result[0];
  if (LHAPDF::numberPDF()%2 != 0) 
    {
      throw domain_error("The number of pdf eigenvectors is not even. Don't know how to proceed.");
    }
  int Imax = LHAPDF::numberPDF()/2; //int division OK; numberPDF is even
  
  double odd, even;
  double tp,tm;
  vector<double> plus,minus;
  for (int i = 1; i<=Imax; ++i ) 
    {
      odd  = result[2*i-1] - central;
      even = result[2*i] - central;
      
      tp = max( max(even,odd) ,0.0);
      tm = max( max(-even,-odd) ,0.0);
      
      plus.push_back(  tp );
      minus.push_back( tm );
    }
  
  double deltaP=0.0;
  for(vector<double>::iterator it = plus.begin(); it != plus.end(); ++it) 
    {
      deltaP += pow( *it ,2);
    }
  double deltaM=0.0;
  for(vector<double>::iterator it = minus.begin(); it != minus.end(); ++it) 
    {
      deltaM += pow( *it ,2);
    }  
  smax = central + sqrt(deltaP);
  smin = central - sqrt(deltaM);
}

//=======================================================================
//  2)     The prescription used by the NNPDF family of sets.       
//=======================================================================

void deltaPDFnnpdf(const vector<double>& result, 
		   double& mean, 
		   double& smax, 
		   double& smin)
{
  // The first element of 'result' has to be excluded;
  // it is an approximation to 'mean' to be derived below.	
  
  int size = result.size()-1;
  double ar[size];
  for(int i = 0; i != size; ++i)
    {
      ar[i] = result[i+1];
    }  
  mean = gsl_stats_mean (ar, 1, size);
  double sd = gsl_stats_sd_m (ar, 1, size, mean);
  smax = mean + sd;
  smin = mean - sd;
}

//=======================================================================
//  3)   The symmetric prescription used by the Alekhin family of sets  
//       (including ABKM09 and ABM11)                     
//=======================================================================

void deltaPDFsymmetric(const vector<double>& result, 
		       double& smax, 
		       double& smin)
{	
  double central = result[0];
  int size = result.size();
  double deltasq = 0;
  for(int i = 0; i != size; ++i)
    {
      deltasq += pow(central-result[i],2);
    }
  double delta = sqrt(deltasq);
  smax = central + delta;
  smin = central - delta;
}

//=======================================================================
//  4)   The prescription for the HERA's *_VAR family of sets  
//=======================================================================

void deltaPDFheraVAR(const vector<double>& result, 
		     double& smax, 
		     double& smin)
{
  double central = result[0];
  double pos1=0, neg1=0;	
  for(int i = 1; i != 9; ++i)
    {
      double x = result[i]-central;
      if (x < 0) 
	{ 
	  neg1 += x*x;
	}
      else 
	{ 
	  pos1 += x*x;
	}
    }
  vector<double> vp,vn;
  for(int i = 9; i != result.size(); ++i)
    {
      double x = result[i]-central;
      if (x < 0) 
	{ 
	  vn.push_back(-x);
	}
      else 
	{ 
	  vp.push_back(x);
	}
    }  
  double pos2=0, neg2=0;
  if (!vp.empty()) 
    {
      pos2 = *max_element(vp.begin(), vp.end());
    }
  if (!vn.empty()) 
    {
      neg2 = *max_element(vn.begin(), vn.end());
    }
  smax = central + sqrt(pos1 + pow(pos2,2));
  smin = central - sqrt(neg1 + pow(neg2,2));
}

//=======================================================================
//  n)  More prescriptions can be defined by the user in the following                     
//=======================================================================

