/*******************************************************************************	
*									       *	
*     Program: Top++ (ver. 2.0)                                                *
*                                                                              *
*     Authors: Michal Czakon (Aahen) and Alexander Mitov (CERN)                *
*     Preprint number: CERN-PH-TH/2011-303, TTK-11-58                          *
*                                                                              *
*     Licence: GNU Public License. No warranty given or implied.               *
*     Website: http://www.alexandermitov.com/software                          *
*     Compiler: Developed and tested with GCC’s C++ compiler.                   *
*     Operating system: Linux; could be adapted for Mac OS X.                  *
*     Program Language: C++                                                    *
*     External libraries: 1) GNU Scientific Library (GSL)                      *
*                         2) the Les Houches Accord PDF Interface (LHAPDF)     *
*                                                                              *
*     What is this program for?     			                       *
*           		                                                       *	
*     This program calculates the top-pair total inclusive cross-section in    *
*     hadron collisions. The program is able to calculate the total            *
*     top pair cross-section in both pure fixed order perturbation theory      *
*     through exact NNLO and by including soft-gluon resummation through       *
*     next-to-next-to-leading logarithmic order (NNLL).                        *
*                                                                              *
* [1] The implementation of soft-gluon resummation  is described in:           *
*        M. Cacciari, M. Czakon, M. L. Mangano, A. Mitov and P. Nason,         *
*        arXiv:1111.5869 [hep-ph].                                             *
* [2] The (dominant) NNLO result for qqbar->tt+X is computed in:               *
*        P. Baernreuther, M. Czakon and A. Mitov, arXiv:1204.5201 [hep-ph].    *
* [3] The NNLO corrections to qq->tt+X, qq'->tt+X, qqbar'->tt+X,               *
*     together with the remaining contribution from qqbar->tt+X                *
*     are computed in the paper:                                               *
*        M.Czakon and A.Mitov, arXiv:1207.0236 [hep-ph].                       *
* [4] The NNLO corrections to gq->tt+X                                         *
*     are computed in the paper:                                               *
*        M.Czakon and A.Mitov, arXiv:1210.6832 [hep-ph].                       * 
* [5] The NNLO corrections to gg-> tt+X                                        *
*     are computed  in the paper:                                              *
*        M.Czakon, P.Fiedler and A.Mitov, arXiv:1303.6254 [hep-ph].            *
*                                                                              *
*******************************************************************************/

#include <map>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <time.h>
#include <complex>
#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_errno.h>

#include "PartonicFlux.h"
#include "SubtrFlux.h"
#include "FixedOrder.h"
#include "Resummation.h"
#include "Utilities.h"
#include "psin.h"
#include "lgamma.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::map;
using std::domain_error;
using std::max;
using std::ofstream;
using std::ifstream;
using std::stringstream;

typedef std::complex<double> dcomp;

int main() {	
  
  gsl_set_error_handler_off();//Supress the GSL-default abort when integration problem occurs. Report it instead.
  time_t start,end, calcstart,calcend; 
  double dif,calcdif;		  
  vector<vector<vector<double> > > resultGlobal;// Stores the global result.	
  vector<vector<double> > resultGlobalSummary;// Stores the summary of the global result.	
  
  time (&calcstart);	
  
  Greetings("2.0");//Startup greeting on the screen	
  
/*******************************************************************************	
*									       *	
*		     Default and user_provided parameters		       *
*									       *	
*******************************************************************************/		
  
  map<string, vector<string> > options;
  const int NumberParameters = 29;//Number of param's !!!!!!!!! Must be modified if new options are added !!!!!!!
  
  // default values for all options
  
  string defaults[NumberParameters*2] = {
    //================================= General setup (Collider, pdf, F.O. vs RES.)
    "Collider", "LHC",
    "WithResummation", "YES",
    "PDFuncertainty", "NO",
    "RestrictedScaleVariation", "2.0",  
    "PDFset", "MSTW2008nnlo68cl",
    "PDFmember", "0",
    //================================= mt(GeV); muF and muR (in units of mt)
    "Mtop", "173.3",
    "MtopLimit", "-1.0",//unattainable value. Do NOT modify!
    "MtopStep", "1",
    "muF", "0.5 1.0 2.0",
    "muR", "0.5 1.0 2.0",
    //================================= Resummation
    "OrderFO", "NNLO",
    "OrderRES", "NNLL",		
    "A", "0.0",
    "TwoLoopCoulombs", "YES",
    "H2gg1", "53.17",//normalization as/(pi)
    "H2gg8", "96.34",
    "H2qq", "84.81",
    //================================= Fixed Order
    "LO", "YES",
    "NLO", "YES",
    "NNLO", "YES",
    //================================= Setup parameters
    "ECMLHC", "8000",
    "ECMTEV", "1960",
    "Precision", "2",
    "NPdfGrid", "100",
    "ETA", "1e-5",
    "CMP", "2.7",
    "PdfFileType", "LHgrid",
    "PartonChannel", "ALL"  
  };
  
  for (int i = 0; i < NumberParameters; ++i)
    {
      string key=defaults[2*i];
      string buffer=defaults[2*i+1];
      stringstream line(buffer);
      string val;
      vector<string> values;		
      while (line)
	{
	  val.clear();
	  line >> val;
	  if (val.empty()) continue;
	  values.push_back(val);
	}      		
      options[key] = values;		
    }
  
  // Override the default options by those from the user (external file). 
  // Every option is a string at this point. 
  // Comments start with a "/"
  
  ifstream config("top++.cfg");
  
  while (config)
    {
      string buffer;
      getline(config, buffer);
      if (buffer.empty()) continue;
      
      stringstream line(buffer);
      string key, val;
      vector<string> values;
      line >> key;
      while (line)
	{
	  val.clear();
	  line >> val;
	  if (val.empty()) continue;
	  if (val[0] == '/') break;
	  values.push_back(val);
	}
      options[key] = values;
    }
  
  // assignment of the actual options
  
  double A = atof(options["A"][0].c_str());
  double CMP = atof(options["CMP"][0].c_str());
  double ECMLHC = atof(options["ECMLHC"][0].c_str());
  double ECMTEV = atof(options["ECMTEV"][0].c_str());
  double ETA = atof(options["ETA"][0].c_str());
  double H2gg1 = atof(options["H2gg1"][0].c_str());
  double H2gg8 = atof(options["H2gg8"][0].c_str());
  double H2qq = atof(options["H2qq"][0].c_str());
  double Mtop = atof(options["Mtop"][0].c_str());
  double MtopLimit   = atof(options["MtopLimit"][0].c_str());
  double MtopStep = atof(options["MtopStep"][0].c_str());
  int NPdfGrid = atoi(options["NPdfGrid"][0].c_str());
  int PDFmember = atoi(options["PDFmember"][0].c_str());
  int Precision = atoi(options["Precision"][0].c_str());
  
  string LO = options["LO"][0];
  string NLO = options["NLO"][0];
  string NNLO = options["NNLO"][0];
  string OrderFO = options["OrderFO"][0];
  string OrderRES = options["OrderRES"][0];
  string PDFuncertainty = options["PDFuncertainty"][0];
  string WithResummation = options["WithResummation"][0];
  string Collider = options["Collider"][0];
  string PDFset = options["PDFset"][0];
  string TwoLoopCoulombs = options["TwoLoopCoulombs"][0];
  string RestrictedScaleVariation = options["RestrictedScaleVariation"][0];
  string PdfFileType = options["PdfFileType"][0];
  string PartonChannel = options["PartonChannel"][0];
  
  // Check if this option's value is allowed (either "NO" or a "double >= 1.0"):
  // Recall: atof(ptr-to-non-double) == 0.0
  if (RestrictedScaleVariation != "NO" && 
      atof(RestrictedScaleVariation.c_str()) < 1.0) 
    {
      throw domain_error("Incorrect option RestrictedScaleVariation: (must be either NO or a double >= 1.0).");
    }
  
  // Set default value MtopLimit=Mtop.
  if (MtopLimit == -1.0){ MtopLimit = Mtop;}
  
  // Guard the range of this option to avoid infinite loops:
  if (MtopStep <= 0) 
    {
      throw domain_error("Incorrect option MtopStep: (must be a positive number).");
    }
  
  if (PdfFileType!="LHgrid" && PdfFileType!="LHpdf")
    {
      throw domain_error("Unknown type of pdf file (possibilities: LHgrid or LHpdf).");
    }
  
  vector<double> vmuF, vmuR;
  if (PDFuncertainty == "YES") 
    {//For pdf variation, use oly central scales.
      vmuF.push_back(1.0);
      vmuR.push_back(1.0);
    }	
  else 
    {
      for (vector<string>::iterator mu = options["muR"].begin();
	   mu != options["muR"].end(); ++mu) vmuR.push_back(atof(mu->c_str()));
      for (vector<string>::iterator mu = options["muF"].begin();
	   mu != options["muF"].end(); ++mu) vmuF.push_back(atof(mu->c_str()));	
      
      if (vmuF.size() == 0) 
	{
	  throw domain_error("Please specify values for muF.");
	}
      if (vmuR.size() == 0) 
	{
	  throw domain_error("Please specify values for muR.");
	}
    }
  
  // When with 'RESUMMATION', set the F.O. parameters according to the RES ones:
  if (WithResummation=="YES") 
    {      
      if (OrderFO=="LO") 
	{	
	  LO =   "YES";
	  NLO =  "NO";
	  NNLO = "NO";
	}
      else if (OrderFO=="NLO") 
	{	
	  LO =   "YES";
	  NLO =  "YES";
	  NNLO = "NO";
	}
      else if (OrderFO=="NNLO") 
	{	
	  LO =   "YES";
	  NLO =  "YES";
	  NNLO = "YES";
	}
      else 
	{
	  throw domain_error("Unknown option for Resummation matching (Expected LO, NLO or NNLO)");
	}
    }	
  
  double ECM;//Collider c.m. energy
  if (Collider == "TEV") { ECM=ECMTEV; }
  else { ECM=ECMLHC; }
  
/*******************************************************************************	
*									       *	
*	      Check for consistency and completeness all string                *
*             options that are passed as bool to the classes.                  *
*									       *	
*******************************************************************************/		
  
  if (Collider!="LHC" && Collider!="TEV")
    {
      throw domain_error("Unknown type of collider (possibilities: TEV or LHC).");
    }
  
  if (LO!="YES" && LO!="NO")
    {
      throw domain_error("Unknown option for the LO result (possibilities: YES or NO).");
    }
  
  if (NLO!="YES" && NLO!="NO")
    {
      throw domain_error("Unknown option for the NLO result (possibilities: YES or NO).");
    }
  
  if (NNLO!="YES" && NNLO!="NO")
    {
      throw domain_error("Unknown option for the NNLO result (possibilities: YES or NO).");
    }
  
  if (TwoLoopCoulombs!="YES" && TwoLoopCoulombs!="NO")
    {	
      throw domain_error("Unknown option about Two-loop Coulombs (possibilities are: YES or NO).");
    }
  
  if (OrderRES!="LL" && OrderRES!="NLL" && OrderRES!="NNLL")
    {	
      throw domain_error("Unknown option with resummation order (possibilities: LL, NLL, NNLL).");
    }
  
  if (OrderFO!="LO" && OrderFO!="NLO" && OrderFO!="NNLO")
    {	
      throw domain_error("Unknown option with resummation matching to F.O. (possibilities: LO, NLO, NNLO).");
    }
  	
/*******************************************************************************	
*									       *	
*                 Setup for the computation of pdf uncertainties               *
*									       *	
*******************************************************************************/		
  
  // ALL prescriptions for comp. PDF uncert's MUST be included in the next 2 lines:
  const int NumPdfPrescr = 4;
  string arrpdf[NumPdfPrescr] = {"Asymmetric","NNPDF","Symmetric","HERA_VAR"};
  
  vector<string> ListKnownPdfPrescriptions;
  for (int i = 0; i < NumPdfPrescr; ++i)
    ListKnownPdfPrescriptions.push_back( arrpdf[i] );
  
  // Store all user-defined pairs (pdf set, pdf uncert. prescription)
  map<string, vector<string> > PdfPairs;			
  {
    ifstream readpdf("pdf.cfg");
    while (readpdf)
      {
	string buffer;
	getline(readpdf, buffer);
	if (buffer.empty()) continue;
	
	stringstream line(buffer);
	string key, val;
	vector<string> values;
	line >> key;			
	while (line)
	  {
	    val.clear();
	    line >> val;
	    if (val.empty()) continue;
	    if (val[0] == '/') break;
	    values.push_back(val);
	  }			
	PdfPairs[key] = values;
      }
  }
  
  // assignment of the actual option for Computing Pdf Uncertainty
  string PdfUncertMethod = "Unknown";
  
  if (PDFuncertainty == "YES") 
    {
      if (PdfPairs.count(PDFset) > 0)//pdf set is in the list: 
	{
	  bool flag = true;
	  
	  if (PdfPairs[PDFset].size() == 0) 
	    {
	      flag = false;
	      PdfPairs[PDFset].push_back("Unknown");
	    }
	  if (PdfPairs[PDFset].size() != 1) 
	    {
	      cout 
		<< "\n!!!!!!!!!!!!!!!!!!!!!!!!!   WARNING   !!!!!!!!!!!!!!!!!!!!!!!!!\n"
		" Too many pdf prescriptions for pdf set "
		<< PDFset 
		<< " (must be one).\n"
		" The program will proceed using the prescription \""
		<< PdfPairs[PDFset][0]
		<< "\".\n" 
		<< "!!!!!!!!!!!!!!!!!!!!!!!   END WARNING   !!!!!!!!!!!!!!!!!!!!!!!\n"
		<< endl;
	    }			
	  
	  // check if the method requested by the user is known:
	  vector<string>::iterator it = find(ListKnownPdfPrescriptions.begin() , 
					     ListKnownPdfPrescriptions.end() , 
					     PdfPairs[PDFset][0]);			
	  
	  if ( it == ListKnownPdfPrescriptions.end() && flag != false) //unknown method
	    {
	      cout 
		<< "\n!!!!!!!!!!!!!!!!!!!!!!!!!   WARNING   !!!!!!!!!!!!!!!!!!!!!!!!!\n"
		" The prescription "
		<< PdfPairs[PDFset][0]
		<< " for computing pdf uncertainty is unknown.\n" 
		" The program will continue without computing pdf uncertainty.\n"
		" The results for the individual pdf members will be displayed.\n"
		<< "!!!!!!!!!!!!!!!!!!!!!!!   END WARNING   !!!!!!!!!!!!!!!!!!!!!!!\n"
		<< endl;
	    }
	  else if (flag == false) //no method was specified 
	    {
	      cout 
		<< "\n!!!!!!!!!!!!!!!!!!!!!!!!!   WARNING   !!!!!!!!!!!!!!!!!!!!!!!!!\n"
		" No prescription for computing pdf uncertainty with set "
		<< PDFset
		<< " is specified.\n" 
		" The program will continue without computing pdf uncertainty.\n"
		" The results for the individual pdf members will be displayed.\n"
		<< "!!!!!!!!!!!!!!!!!!!!!!!   END WARNING   !!!!!!!!!!!!!!!!!!!!!!!\n"
		<< endl;				
	    }
	  else //the method PdfPairs[PDFset][0] is known
	    {
	      PdfUncertMethod = PdfPairs[PDFset][0];
	      cout 
		<< "NOTE: Pdf uncertainty will be computed with prescription: "
		<< PdfUncertMethod
		<< "\n"
		<< endl;
	    }
	}
      else //pdf set not in the list:
	{
	  cout 
	    << "\n!!!!!!!!!!!!!!!!!!!!!!!!!   WARNING   !!!!!!!!!!!!!!!!!!!!!!!!!\n"
	    " No prescription for computing pdf uncertainty with set "
	    << PDFset
	    << " is specified.\n" 
	    " The program will continue without computing pdf uncertainty.\n"
	    " The results for the individual pdf members will be displayed.\n"
	    << "!!!!!!!!!!!!!!!!!!!!!!!   END WARNING   !!!!!!!!!!!!!!!!!!!!!!!\n"
	    << endl;
	}
    }
  
/*******************************************************************************	
*									       *	
*       		     Initialize the PDF set			       *
*									       *	
*******************************************************************************/	
  
  time (&start);	
  if (PdfFileType == "LHgrid") 
    {
      initPDFSetByName(PDFset,LHAPDF::LHGRID);
    } 
  else 
    {
      initPDFSetByName(PDFset,LHAPDF::LHPDF);
    }
  time (&end);
  dif = difftime (end,start);
  cout 
    << "................. PDF initialized (it took " 
    << dif 
    << " [sec.]).\n" 
    << endl;
  
/*******************************************************************************	
*									       *	
*		       BEGIN: loop over the top mass mt			       *
*									       *	
*******************************************************************************/		
  
  //check user's input for consistency	
  if (MtopLimit < Mtop) 
    {
      MtopLimit = Mtop;
      cout 
	<< "\n!!!!!!!!!!!!!!!!!!!!!!!!!   WARNING   !!!!!!!!!!!!!!!!!!!!!!!!!\n"
	" MtopLimit is smaller than Mtop. \n"
	" The program will continue with MtopLimit=Mtop.\n"
	<< "!!!!!!!!!!!!!!!!!!!!!!!   END WARNING   !!!!!!!!!!!!!!!!!!!!!!!\n"
	<< endl;
    }	
  
  bool WillOutputSummary = true;//if true at the end, will output in file summary of each run
  
  for (double mt=Mtop; mt <=MtopLimit; mt+=MtopStep) {
    double Rhad = 4*pow(mt/ECM,2);
    
/*******************************************************************************	
*									       *	
*        -  Initialize the partonic fluxes			               *
*        -  Initialize alpha_s (for each pdf member or for each muR)           * 
*									       *	
*******************************************************************************/	
    
    time (&start);
    vector<PartonicFlux> vFlux;
    vector<double> AlphaS;
    
    if (PDFuncertainty == "YES") //all pdf members; a single muF value (=mt)
      {
	cout 
	  << "Initializing fluxes and alpha_s for each pdf member ... \n" 
	  << endl;	
	for(int i = 0; i <= LHAPDF::numberPDF(); ++i) 
	  {
	    PartonicFlux flux(Collider=="LHC" ? 1:0 , 
			      i, 
			      Rhad, 
			      1.0, //Upper integration limit
			      vmuF[0] * mt, 
			      NPdfGrid, 
			      Precision);
	    vFlux.push_back( flux );
	    //alphasPDF always calling current pdf member (i.e. "i"):
	    AlphaS.push_back( LHAPDF::alphasPDF( vmuR[0] * mt ) );  		  
	    cout
	      << i 
	      << " out of "
	      << 1+LHAPDF::numberPDF() 
	      << ": alpha_s(muR="
	      << vmuR[0] * mt 
	      << ")=" 
	      << LHAPDF::alphasPDF(vmuR[0] * mt)
	      << ", alpha_s(mZ=91.1876)=" 
	      << LHAPDF::alphasPDF(91.1876) 
	      << endl;
	  }
      }
    else // all user's muF values; A single pdf member
      {
	cout 
	  << "Initializing fluxes and alpha_s for pdf member " 
	  << PDFmember 
	  <<" ...\n" 
	  << endl;	
	for(vector<double>::iterator it = vmuF.begin(); it != vmuF.end(); ++it) 
	  {
	    PartonicFlux flux(Collider=="LHC" ? 1:0, 
			      PDFmember, 
			      Rhad, 
			      1.0, 
			      (*it)*mt, 
			      NPdfGrid, 
			      Precision);
	    vFlux.push_back( flux );
	  }
	cout 
	  << "alpha_s(mZ=91.1876)=" 
	  << LHAPDF::alphasPDF(91.1876) 
	  << endl;		  		  
	for(vector<double>::iterator it = vmuR.begin(); it != vmuR.end(); ++it) 
	  {
	    AlphaS.push_back( LHAPDF::alphasPDF( (*it)*mt ) );
	    cout 
	      << "alpha_s(muR="
	      << (*it)*mt 
	      << ")=" 
	      << LHAPDF::alphasPDF((*it)*mt) 
	      << endl;
	  }
      }
    time (&end);
    dif = difftime (end,start);
    cout 
      << "\n.............. Fluxes initialized (it took " 
      << dif 
      << " [sec.]).\n" 
      << endl;
    
/*******************************************************************************	
*									       *	
*		       	Initialize the FixedOrder result		       *
*									       *	
*******************************************************************************/
    
    time (&start);
    cout << "\nInitializing Fixed Order ..." 
	 << endl;
    vector<FixedOrder> vFO;
    
    if (PDFuncertainty == "YES") 
      {
	for(int i = 0; i <= LHAPDF::numberPDF(); ++i) 
	  {
	    FixedOrder FO(AlphaS[i],
			  mt, 
			  vmuF[0] * mt, 
			  vmuR[0] * mt, 
			  LO=="YES" ? 1:0 , 
			  NLO=="YES" ? 1:0 , 
			  NNLO=="YES" ? 1:0);
	    vFO.push_back( FO );
	  }	
      }
    else 
      {
	for(int i=0; i != vmuF.size(); ++i) 
	  {
	    for(int j=0; j != vmuR.size(); ++j) 
	      {
		FixedOrder FO(AlphaS[j],
			      mt, 
			      vmuF[i] * mt, 
			      vmuR[j] * mt, 
			      LO=="YES" ? 1:0 , 
			      NLO=="YES" ? 1:0 , 
			      NNLO=="YES" ? 1:0);
		vFO.push_back( FO );
	      }
	  }		  
      }
    time (&end);
    dif = difftime (end,start);
    cout 
      << "......... Fixed Order initialized (it took " 
      << dif 
      << " [sec.]).\n" 
      << endl;  
    
/*******************************************************************************	
*      					                                       *	
*                        Initialize the Subtraction flux 		       *
*								               *	
*******************************************************************************/
    
    vector<SubtrFlux> vSF;
    if (WithResummation == "YES") 
      {					
	time (&start);
	cout << "Initializing the subtraction fluxes ... " 
	     << endl;	
	
	for(vector<PartonicFlux>::iterator it = vFlux.begin(); it != vFlux.end(); ++it) 
	  {
	    SubtrFlux SF(*it , Rhad, ETA);
	    vSF.push_back( SF );
	  }
	time (&end);
	dif = difftime (end,start);
	cout 
	  << "...... Subtraction Fluxes initialized (it took "
	  << dif 
	  << " [sec.]).\n" 
	  << endl;
      }
    
/*******************************************************************************	
*					 			               *	
*			      Initialize Resummation			       *
*								       	       *	
*******************************************************************************/
    
    vector<Resummation> vRES;
    if (WithResummation == "YES") 
      {		
	time (&start);
	cout 
	  << "Initializing resummation ... " 
	  << endl;
	if (PDFuncertainty == "YES") 
	  {
	    for(int i = 0; i <= LHAPDF::numberPDF(); ++i) 
	      {	
		Resummation RES(AlphaS[i],
				vSF[i], 
				Precision, 
				mt, 
				vmuF[0] * mt, 
				vmuR[0] * mt, 
				A, 
				CMP, 
				IntegerForm(OrderFO), 
				IntegerForm(OrderRES), 
				TwoLoopCoulombs=="YES" ? 1:0 , 
				H2qq,H2gg1,H2gg8);
		vRES.push_back( RES );
	      }
	  }
	else 
	  {
	    for(int i=0; i != vmuF.size(); ++i) 
	      {
		for(int j=0; j != vmuR.size(); ++j) 
		  {	
		    Resummation RES(AlphaS[j],
				    vSF[i], 
				    Precision, 
				    mt, 
				    vmuF[i] * mt, 
				    vmuR[j] * mt, 
				    A, 
				    CMP, 
				    IntegerForm(OrderFO),
				    IntegerForm(OrderRES), 
				    TwoLoopCoulombs=="YES" ? 1:0 , 
				    H2qq,H2gg1,H2gg8);
		    vRES.push_back( RES );
		  }
	      }
	  }
	time (&end);
	dif = difftime (end,start);
	cout 
	  << "...... Resummation initialization done (it took " 
	  << dif 
	  << " [sec.]).\n" 
	  << endl;
      }
    
/*******************************************************************************	
*		 			                                       *	
*			  Compute the final result			       *
*								               *	
*******************************************************************************/
    
    time (&start);
    cout 
      << "\nComputing the cross-section:\n" 
      << endl;
    
    double central = -1.0;//an unattainable value for sigma_tot
    vector<double> result;// all values of x-section for fixed mt
    
    vector<vector<double> > resultPrint;// with pdf uncert  : [mt, pdf member, muF, muR, alpha_s, sigma_tot]
    //  with scale uncert: [mt, muF, muR, alpha_s, sigma_tot]
    
    if (PDFuncertainty == "YES") //compute sigma for all pdf members, and for a single (muF,muR) value 
      {
	if (vFO.size() != 1+LHAPDF::numberPDF() || vmuF.size() != 1 || vmuR.size() != 1) 
	  {
	    throw domain_error("Inconsistent number of (muF,muR)"
			       " pairs with PDF uncertainty determination."
			       " (should be only one), " 
			       " Or inconsistent number of pdf members");
	  }
	bool pass = ScalesPass(RestrictedScaleVariation,vmuF[0],vmuR[0]);
	for(int i = 0; i <= LHAPDF::numberPDF(); ++i) 
	  {
	    vector<double> element;
	    if (pass)
	      {
		element.push_back(mt);
		element.push_back(i);
		element.push_back(vmuF[0]);
		element.push_back(vmuR[0]);
		element.push_back(AlphaS[i]);
		if (WithResummation == "YES") 
		  {			
		    result.push_back(  ComputeFinalResultFO (vFlux[i],vFO[i],Rhad,Precision,PartonChannel)
				       + ComputeFinalResultRES(vFlux[i], vSF[i], vFO[i], vRES[i], Rhad, Precision, PartonChannel) );
		    element.push_back(result.back());
		  }
		else 
		  {
		    result.push_back(  ComputeFinalResultFO (vFlux[i],vFO[i],Rhad,Precision,PartonChannel) );
		    element.push_back(result.back());
		  }
		if (i == 0) 
		  {
		    central = result.back();
		  }
		cout 
		  << "m_top="
		  << mt
		  << ", Pdf member " 
		  << i 
		  << ", muF="
		  << vmuF[0]
		  << ", muR="
		  << vmuR[0]
		  << ", sigma_tot=" 
		  << result.back()
		  << " [pb]." 	
		  << endl;
		
		resultPrint.push_back(element);
	      }
	  }
	time (&end);
	dif = difftime (end,start);
	
	if (result.empty()) 
	  {
	    throw domain_error("The pair (muF,muR) does not meet the restricted scale variation criteria."
			       " Please review muF and muR.");
	  }
	else 
	  {
	    cout 
	      << "\n......... Cross-section computed. (it took " 
	      << dif 
	      << " [sec.]).\n" 
	      << endl;
	  }
      }	
    else //compute sigma for the central pdf member, and for all user supplied (muF,muR) values
      {
	for(int i = 0; i != vmuF.size(); ++i) 
	  {
	    for(int j = 0; j != vmuR.size(); ++j) 
	      {
		vector<double> element;
		
		bool pass = ScalesPass(RestrictedScaleVariation,vmuF[i],vmuR[j]);
		if (pass)
		  {
		    element.push_back(mt);
		    element.push_back(vmuF[i]);
		    element.push_back(vmuR[j]);
		    element.push_back(AlphaS[j]);
		    if (WithResummation == "YES") 
		      {				
			result.push_back(  ComputeFinalResultFO (vFlux[i],vFO[ i * vmuR.size() + j ],
								 Rhad,Precision,PartonChannel)
					   + ComputeFinalResultRES(vFlux[i], vSF[i], vFO[ i * vmuR.size() + j ], 
								   vRES[ i * vmuR.size() + j ], Rhad, Precision, PartonChannel) );
			element.push_back(result.back());
		      }
		    else 
		      {
			result.push_back(  ComputeFinalResultFO (vFlux[i],vFO[ i * vmuR.size() + j ],Rhad,Precision,PartonChannel) );  
			element.push_back(result.back());
		      }
		    if (vmuF[i] == 1.0 && vmuR[j] == 1.0) 
		      {
			central = result.back();
		      }
		    cout 
		      << "m_top="
		      << mt
		      << ", muF="
		      << vmuF[i]
		      << ", muR="
		      << vmuR[j]
		      << ", sigma_tot=" 
		      << result.back() 
		      << " [pb]." 
		      << endl;
		    
		    resultPrint.push_back(element);
		  }
	      }
	  }
	time (&end);
	dif = difftime (end,start);
	
	if (result.empty()) 
	  {
	    throw domain_error("Not a single pair (muF,muR) meets the restricted scale-variation criteria."
			       " Please review your muF and muR or the option RestrictedScaleVariation.");
	  }
	else 
	  {
	    cout 
	      << "\n......... Cross-section computed! (it took " 
	      << dif 
	      << " [sec.]).\n" 
	      << endl;
	  }
      }
    
/*******************************************************************************	
*									       *	
*                Compute scale/pdf variation for a fixed mtop                  *
*									       *	
*******************************************************************************/	
    
    double s0 = central;
    double mean=-1.0, smax=-1.0, smin=-1.0;//unattainable values
    
    if (PDFuncertainty == "YES")
      {
	if (PdfUncertMethod == "NNPDF")
	  {
	    deltaPDFnnpdf(result, mean, smax, smin);
	  }
	else if (PdfUncertMethod == "Symmetric")
	  {
	    deltaPDFsymmetric(result, smax, smin);
	  }		  
	else if (PdfUncertMethod == "HERA_VAR")
	  {
	    deltaPDFheraVAR(result, smax, smin);
	  }		  		  
	else if (PdfUncertMethod == "Asymmetric")
	  {
	    deltaPDFasymmetric(result, smax, smin);
	  }
      }
    else//scale uncertainty
      {
	smax = *max_element(result.begin(),result.end());
	smin = *min_element(result.begin(),result.end());
      }
    
    resultGlobal.push_back(resultPrint);
    
/*******************************************************************************	
*									       *	
*                         Output on the screen for fixed mtop.                 *
*									       *	
*******************************************************************************/	
	  	  
    bool WillOutput = true;
    
    if (s0 == -1.0) 
      {
	WillOutput = false;
	
	cout 
	  << "\n!!!!!!!!!!!!!!!!!!!!!!!!!   WARNING   !!!!!!!!!!!!!!!!!!!!!!!!!\n"
	  " Central value cannot be determined.\n"
	  " Please include the point muF=muR=mt.\n" 
	  << "!!!!!!!!!!!!!!!!!!!!!!!   END WARNING   !!!!!!!!!!!!!!!!!!!!!!!\n"
	  << endl;
      }
    
    if (PDFuncertainty == "YES" && PdfUncertMethod == "Unknown")
      {
	WillOutput = false;
      }	  
    
    if (WillOutput)
      {
	if (PDFuncertainty == "YES") 
	  {
	    cout 
	      << "****************    Final result (PDF uncertainty):   ****************\n" 
	      << endl;		
	  }
	else 
	  {
	    cout 
	      << "****************    Final result (scale variation):   ****************\n" 
	      << endl;
	  }
	
	if (PartonChannel=="gg" || 
	    PartonChannel=="qg" || 
	    PartonChannel=="qqbar" || 
	    PartonChannel=="qqbarprime" || 
	    PartonChannel=="qq" || 
	    PartonChannel=="qqprime") 
	  {
	    cout 
	      << " (Calculation of a single partonic channel: "
	      << PartonChannel 
	      << " -> ttbar+X)\n"
	      << endl;
	  }
	cout 
	  << " sigma_tot = " 
	  << s0 
	  << " + " 
	  << smax - s0 
	  << " (" 
	  << (smax - s0)/s0*100.0 
	  << " \%) - " 
	  << s0 - smin 
	  << " (" 
	  << (s0 - smin)/s0*100.0 
	  << " \%)"
	  << " [pb].\n"; 
	if (PDFuncertainty == "YES" && PdfUncertMethod == "NNPDF")
	  {
	    cout 
	      << " Statistical mean = "
	      << mean
	      << " [pb] (Note: this is NNPDF-specific info)\n";
	  }	  
	cout
	  << "**********************************************************************\n" 
	  << endl;
	
	vector<double> summary;// summary of the result; for fixed m_top
	summary.push_back(mt); 
	summary.push_back(s0); 
	summary.push_back(smin); 
	summary.push_back(smax); 		  
	
	resultGlobalSummary.push_back(summary);
      }
    
    if (!WillOutput) //if for any mt WillOutput=0, then no Summary for any mt
      {
	WillOutputSummary = false;
	
	if (PartonChannel=="gg" || 
	    PartonChannel=="qg" || 
	    PartonChannel=="qqbar" || 
	    PartonChannel=="qqbarprime" || 
	    PartonChannel=="qq" || 
	    PartonChannel=="qqprime") 
	  {
	    cout 
	      << "\n>>>>> NOTE: calculation of a single partonic channel: "
	      << PartonChannel 
	      << " -> ttbar+X <<<<<\n"
	      << endl;
	  }		  
      }
    
  }//END loop over top mass mt.
  
  time (&calcend);
  calcdif = difftime (calcend,calcstart);
  
  cout 
    << "\n\tTotal time used: " 
    << calcdif 
    << " [sec.].\n" 
    << endl;
  	
/*******************************************************************************	
*									       *	
*       		  Print the results in a file	                       *
*									       *	
*******************************************************************************/
  
  ofstream file("top++.res");
  
  if (PartonChannel=="gg" || 
      PartonChannel=="qg" || 
      PartonChannel=="qqbar" || 
      PartonChannel=="qqbarprime" || 
      PartonChannel=="qq" || 
      PartonChannel=="qqprime") 
    {
      file 
	<< "#>>>>> Calculation of a single partonic channel: "
	<< PartonChannel 
	<< " -> ttbar+X <<<<<"
	<< endl;
    }

  file 
    << "#--------------------------- Input parameters ---------------------------#"
    << endl;
  
  if (PDFuncertainty == "YES")
    {
      file 
	<< "# PDF variation for sigma_tot" 
	<< "\n# m_top in the range (" 
	<< Mtop 
	<< "," 
	<< MtopLimit 
	<< ")" 
	<< endl;
      
      file 
	<< "# Collider: " 
	<< Collider 
	<< "; c.m. energy " 
	<< ECM << " GeV" 
	<< endl;	
      
      if (WillOutputSummary) 
	{
	  file 
	    << "# PDF set: " 
	    << PDFset
	    << "; alpha_s(PDFmember "
	    << PDFmember
	    << ", mZ=91.1876)=" 
	    << LHAPDF::alphasPDF(91.1876)
	    << "\n# pdf uncertainty prescription: "
	    << PdfUncertMethod
	    << endl;
	}
      else
	{
	  file 
	    << "# PDF set: " 
	    << PDFset
	    << "; alpha_s(PDFmember "
	    << PDFmember
	    << ", mZ=91.1876)=" 
	    << LHAPDF::alphasPDF(91.1876)		
	    << endl;			
	}
      if (WithResummation == "YES")
	{
	  file 
	    << "# Resummed calculation: " 
	    << OrderFO 
	    << " + " 
	    << OrderRES  
	    << "\n# A=" 
	    << A 
	    << ", TwoLoopCoulombs="
	    << TwoLoopCoulombs
	    << ", H2gg1="
	    << H2gg1
	    << ", H2gg8="
	    << H2gg8
	    << ", H2qq="
	    << H2qq
	    << endl;	
	}
      else 
	{
	  file 
	    << "# Fixed Order calculation: "
	    "LO(" 
	    << LO 
	    << "), NLO(" 
	    << NLO 
	    << "), NNLO(" 
	    << NNLO 
	    <<")" 
	    << endl;	
	}
      file 
	<< "# Precision=" 
	<< Precision
	<< "; NPdfGrid=" 
	<< NPdfGrid
	<< "; ETA=" 
	<< ETA
	<< "; CMP=" 
	<< CMP
	<< endl;					
      
      if (WillOutputSummary) 
	{
	  file 
	    << "\n#---------------------- Summary of the calculation ----------------------#\n"
	    "# mt   central   pdf_min   pdf_max"
	    << endl;
	  
	  for(vector<vector<double> >::iterator imt = resultGlobalSummary.begin(); imt != resultGlobalSummary.end(); ++imt) 
	    {
	      for(vector<double>::iterator r = (*imt).begin(); r != (*imt).end(); ++r) 
		{	 
		  file <<  *r << string(2,' ');
		}
	      file << endl;
	    }
	}
      
      file 
	<< "\n#---------------- Detailed results from the calculation -----------------#\n"
	"# mt   pdf member   muF   muR   alpha_s   sigma_tot"		
	<< endl;
      for(vector<vector<vector<double> > >::iterator a = resultGlobal.begin(); a != resultGlobal.end(); ++a)
	{
	  for(vector<vector<double> >::iterator b = (*a).begin(); b != (*a).end(); ++b) 
	    {	 
	      for(vector<double>::iterator c = (*b).begin(); c != (*b).end(); ++c) 
		{	 
		  file <<  *c << string(8,' ');
		}
	      file << endl;
	    }
	}
    }	
  else 
    {
      if (RestrictedScaleVariation == "NO")
	{
	  file 
	    << "# Unrestricted scale variation for sigma_tot.\n"
	    "# m_top in the range (" 
	    << Mtop 
	    << "," 
	    << MtopLimit 
	    << ")." 
	    << endl;
	}
      else
	{
	  file 
	    << "# Restricted scale variation for sigma_tot with "
	    "RestrictedScaleVariation = " 
	    << RestrictedScaleVariation 
	    << ".\n"
	    "# m_top in the range (" 
	    << Mtop 
	    << "," 
	    << MtopLimit 
	    << ")." 
	    << endl;
	}
      
      file 
	<< "# Collider: " 
	<< Collider 
	<< "; c.m. energy " 
	<< ECM << " GeV" 
	<< endl;	
      
      file 
	<< "# PDF set: " 
	<< PDFset
	<< "; PDF member: " 
	<< PDFmember
	<< "; alpha_s(mZ=91.1876)=" 
	<< LHAPDF::alphasPDF(91.1876)
	<< endl;	
      
      if (WithResummation == "YES")
	{
	  file 
	    << "# Resummed calculation: " 
	    << OrderFO 
	    << " + " 
	    << OrderRES  
	    << "\n# A=" 
	    << A 
	    << ", TwoLoopCoulombs="
	    << TwoLoopCoulombs
	    << ", H2gg1="
	    << H2gg1
	    << ", H2gg8="
	    << H2gg8
	    << ", H2qq="
	    << H2qq
	    << endl;
	}
      else 
	{
	  file 
	    << "# Fixed Order calculation: "
	    "LO(" 
	    << LO 
	    << "), NLO(" 
	    << NLO 
	    << "), NNLO(" 
	    << NNLO
	    << ")"
	    << endl;	
	}
      file 
	<< "# Precision=" 
	<< Precision
	<< "; NPdfGrid=" 
	<< NPdfGrid
	<< "; ETA=" 
	<< ETA
	<< "; CMP=" 
	<< CMP
	<< endl;			
      
      if (WillOutputSummary) 
	{
	  file 
	    << "\n#---------------------- Summary of the calculation ----------------------#\n"
	    "# mt   central   scale_min   scale_max"
	    << endl;
	  
	  for(vector<vector<double> >::iterator imt = resultGlobalSummary.begin(); imt != resultGlobalSummary.end(); ++imt) 
	    {
	      for(vector<double>::iterator r = (*imt).begin(); r != (*imt).end(); ++r) 
		{	 
		  file <<  *r << string(2,' ');
		}
	      file << endl;
	    }
	}
      
      file
	<< "\n#---------------- Detailed results from the calculation -----------------#\n"
	"# mt      muF      muR      alpha_s      sigma_tot"		
	<< endl;
      for(vector<vector<vector<double> > >::iterator a = resultGlobal.begin(); a != resultGlobal.end(); ++a)
	{
	  for(vector<vector<double> >::iterator b = (*a).begin(); b != (*a).end(); ++b) 
	    {	 
	      for(vector<double>::iterator c = (*b).begin(); c != (*b).end(); ++c) 
		{	 
		  file <<  *c << string(8,' ');
		}
	      file << endl;
	    }
	}
    }
  
  return 0;
}
