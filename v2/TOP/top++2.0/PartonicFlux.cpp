#include "PartonicFlux.h"

using std::vector;
using std::domain_error;
using std::size_t;

void PartonicFlux::create()
{
  logtmin = log(tmin);
  logtmax = log(tmax);
  if (nmax <= 0) 
    {
      throw domain_error("Number of grid points nmax<=0.");
    }
  
  //initialize one specific member of the pdf set:
  LHAPDF::initPDF(PDFmember);
  
  //initialize the flux (i.e. compute it once and for all)
  initflux();	
}

/*******************************************************************************	
*     Parton flux as a standard convolution of two                             *
*     PDFs at the point t = x1*x2, multiplied by t.                            *
*     Implemented as standard convolution of x*pdf                             *
*     at the point t=x1*x2 but without a factor of t.                          *
*        state can be either                                                   *
*           0 - gluon flux                                                     *
*           1 - q-qbar flux                                                    *
*           2 - gluon-quark flux                                               *
*           3 - quark-quark flux                                               *
*           4 - quark-quark' flux (with q=/=q')                                *
*           5 - quark-quark_bar' flux (with q=/=qbar' or q=/=q')               *
*******************************************************************************/

double PartonicFlux::fluxgg(double x)
{
  return flux(0, x);
}

double PartonicFlux::fluxqqbar(double x)
{
  return flux(1, x);
}

double PartonicFlux::fluxgq(double x)
{
  return flux(2, x);
}

double PartonicFlux::fluxqq(double x)
{
  return flux(3, x);
}

double PartonicFlux::fluxqqprime(double x)
{
  return flux(4, x);
}

double PartonicFlux::fluxqqbarprime(double x)
{
  return flux(5, x);
}

/*******************************************************************************	
*      This function 'changes' the flux from discrete -> continuous.           *
*******************************************************************************/

double PartonicFlux::flux(int state, 
			  double t)
{
  if (t>tmax||t<tmin) 
    {
      throw domain_error("Value of x passed to \"flux\" is out of range.");
    }
  double lt0 = log(t)-logtmin;
  double dt0 = (logtmax-logtmin)/nmax;
  int i = int(lt0/dt0)+1;
  if (i <= 0) i = 1;
  if (i >= nmax) i = nmax-1;
  
  double f0,f1,f2;
  if (state == 0) 
    {
      f0 = vfgg[i-1];
      f1 = vfgg[i];
      f2 = vfgg[i+1];
    }
  else if (state == 1) 
    {
      f0 = vfqqbar[i-1];
      f1 = vfqqbar[i];
      f2 = vfqqbar[i+1];
    }
  else if (state == 2) 
    {
      f0 = vfgq[i-1];
      f1 = vfgq[i];
      f2 = vfgq[i+1];
    }
  else if (state == 3) 
    {
      f0 = vfqq[i-1];
      f1 = vfqq[i];
      f2 = vfqq[i+1];
    }
  else if (state == 4) 
    {
      f0 = vfqqprime[i-1];
      f1 = vfqqprime[i];
      f2 = vfqqprime[i+1];
    }
  else if (state == 5) 
    {
      f0 = vfqqbarprime[i-1];
      f1 = vfqqbarprime[i];
      f2 = vfqqbarprime[i+1];
    }	
  else 
    {
      throw domain_error("Incorrect state passed to function flux.");
    }
  double dt = lt0-dt0*i;
  double iflux= (f0*dt*(dt-dt0)/2+f1*(pow(dt0,2)-pow(dt,2))
		 + f2*dt*(dt+dt0)/2)/pow(dt0,2);
  return iflux;
}

/*******************************************************************************	
*     Compute the fluxes on a grid and store them in vectors.                  *
*******************************************************************************/

void PartonicFlux::initflux()
{
  vfqqbar.clear();
  vfgg.clear();
  vfgq.clear();
  vfqq.clear();
  vfqqprime.clear();
  vfqqbarprime.clear();
  
  for (size_t i=0; i<nmax; i++) 
    {
      double logt = (i*(logtmax-logtmin))/nmax+logtmin;
      double t = exp(logt);
      
      vfgg.push_back( PDFConvolution(0,t,muF,prec) );
      if (ColliderName)//LHC 
	{
	  vfqqbar.push_back( PDFConvolution(1,t,muF,prec) );
	  vfqq.push_back( PDFConvolution(4,t,muF,prec) );
	  vfqqprime.push_back( PDFConvolution(6,t,muF,prec) );
	  vfqqbarprime.push_back( PDFConvolution(8,t,muF,prec) );
	}
      else//TEV 
	{
	  vfqqbar.push_back( PDFConvolution(2,t,muF,prec) );
	  vfqq.push_back( PDFConvolution(5,t,muF,prec) );
	  vfqqprime.push_back( PDFConvolution(7,t,muF,prec) );
	  vfqqbarprime.push_back( PDFConvolution(9,t,muF,prec) );
	}
      vfgq.push_back( PDFConvolution(3,t,muF,prec) );
    }
  
  // the flux vanishes at the end-point nmax, but the PDFs are too unstable to integrate:
  vfgg.push_back(0.0);
  vfqqbar.push_back(0.0);
  vfgq.push_back(0.0);
  vfqq.push_back(0.0);
  vfqqprime.push_back(0.0);
  vfqqbarprime.push_back(0.0);
}

/*******************************************************************************	
*     PDFConvolution.                                                          *
*        state can be either:                                                  *
*            0 - gluon flux                                                    *
*            1 - q-qbar flux in proton-proton collisions                       *
*            2 - q-qbar flux in proton-anti-proton collisions                  *
*            3 - gluon-quark flux                                              *
*            4 - q-q flux in proton-proton collisions                          *
*            5 - q-q flux in proton-anti-proton collisions                     *
*            6 - q-q' flux in proton-proton collisions                         *
*            7 - q-q' flux in proton-anti-proton collisions                    *
*            8 - q-qbar' flux in proton-proton collisions                      *
*            9 - q-qbar' flux in proton-anti-proton collisions                 *
*******************************************************************************/

double PartonicFlux::PDFConvolution(int state, 
				    double t, 
				    double scaleF, 
				    int prec)
{
  params temp_p = {state, t, scaleF};
  gsl_function F;
  F.function = f_p;
  F.params = &temp_p;
  //------------------- setup for the integration routine (not meant to be modified):
  double a = t;
  double b = 1.0;
  double epsabs = 0.0;
  double epsrel = 1/pow(10,prec);
  size_t limit = 1000;
  int key = 6;
  //-------------------
  double result, abserr;
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);	
  gsl_integration_qag (&F, a, b, epsabs, epsrel,limit, key, w, &result, &abserr); 
  gsl_integration_workspace_free (w);
  
  if (abserr/result > epsrel) 
    {
      std::cout << "Warning (Flux initiation): relative error reached " 
		<< abserr/result 
		<< std::endl;
    }
  return result;
}

/*******************************************************************************	
*        A special GSL wrap for the integrand of PDFconvolutionc               *
*******************************************************************************/

double PartonicFlux::f_p (double x, void * par)
{
  params * p = (params*) par;
  double f = PDFproduct((p->state), x, (p->t)/x, (p->scaleF));
  return f;
}

/*******************************************************************************	
*                The integrand of PDFconvolution                               *
*******************************************************************************/

double PartonicFlux::PDFproduct(int state, 
				double x1, 
				double x2, //x2=t/x1 
				double scaleF)
{
  double PDFKrnl = 0.0;
  if (state == 0)
    {
      PDFKrnl += LHAPDF::xfx(x1,scaleF,0)*LHAPDF::xfx(x2,scaleF,0);
    }
  else
    {
      for (int i = -5; i<6; i++) 
	{
	  if (i!=0) 
	    {
	      if (state == 1 || state == 5) 
		{
		  PDFKrnl += LHAPDF::xfx(x1,scaleF,+i)*LHAPDF::xfx(x2,scaleF,-i);
		}
	      else if (state == 2 || state == 4) 
		{
		  PDFKrnl += LHAPDF::xfx(x1,scaleF,+i)*LHAPDF::xfx(x2,scaleF,+i);
		}
	      else if (state == 3) 
		{
		  PDFKrnl += LHAPDF::xfx(x1,scaleF,0)*LHAPDF::xfx(x2,scaleF,+i) 
		    + LHAPDF::xfx(x1,scaleF,+i)*LHAPDF::xfx(x2,scaleF,0);
		}
	      //-------------------------------
	      else if (state == 6 || state == 9) 
		{
		  for (int j = -5; j<6; j++)
		    {
		      if (j != 0 && j != i && j*i > 0)
			{
			  PDFKrnl += LHAPDF::xfx(x1,scaleF,+i)*LHAPDF::xfx(x2,scaleF,+j);
			}
		    }
		}
	      else if (state == 7 || state == 8) 
		{
		  for (int j = -5; j<6; j++)
		    {
		      if (j != 0 && j != i && j*i > 0)
			{
			  PDFKrnl += LHAPDF::xfx(x1,scaleF,+i)*LHAPDF::xfx(x2,scaleF,-j);
			}
		    }
		}
	    }
	}		
    }
  PDFKrnl /= x1;
  return PDFKrnl;
}
