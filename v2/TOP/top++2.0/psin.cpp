#include "psin.h"

typedef std::complex<double> dcomp;

dcomp psin(dcomp z, int k)
{
  double a,b,y,t;
  dcomp x,u,v,h,r,p;
  
  dcomp Cunit(1.0,0.0);
  double delta = 5e-13;                                                 
  double R1 = 1.0, HF = R1/2.0;                                             
  double PI = 3.14159265358979324;
  double C1 = pow(PI,2), C2 = 2*pow(PI,3), C3 = 2*pow(PI,4), C4 = 8*pow(PI,5);          
  
  double sgn[]= {-1.0,+1.0,-1.0,+1.0,-1.0};
  double fct[]= {0.0,1.0,1.0,2.0,6.0,24.0};
  
  double C[7][5];
  
  C[0][0] = 0.0;
  C[1][0] =  8.33333333333333333e-2;                                  
  C[2][0] = -8.33333333333333333e-3;                                  
  C[3][0] =  3.96825396825396825e-3;                                  
  C[4][0] = -4.16666666666666667e-3;                                  
  C[5][0] =  7.57575757575757576e-3;                                  
  C[6][0] = -2.10927960927960928e-2;                                  
  
  C[0][1] = 0.0;
  C[1][1] =  1.66666666666666667e-1;                                  
  C[2][1] = -3.33333333333333333e-2;                                  
  C[3][1] =  2.38095238095238095e-2;                                  
  C[4][1] = -3.33333333333333333e-2;                                  
  C[5][1] =  7.57575757575757576e-2;                                  
  C[6][1] = -2.53113553113553114e-1;                                  
  
  C[0][2] = 0.0;
  C[1][2] =  5.00000000000000000e-1;                                  
  C[2][2] = -1.66666666666666667e-1;                                  
  C[3][2] =  1.66666666666666667e-1;                                  
  C[4][2] = -3.00000000000000000e-1;                                  
  C[5][2] =  8.33333333333333333e-1;                                  
  C[6][2] = -3.29047619047619048e+0;                                  
  
  C[0][3] = 0.0;
  C[1][3] =  2.00000000000000000e+0;                                  
  C[2][3] = -1.00000000000000000e+0;                                  
  C[3][3] =  1.33333333333333333e+0;                                  
  C[4][3] = -3.00000000000000000e+0;                                  
  C[5][3] =  1.00000000000000000e+1;                                  
  C[6][3] = -4.60666666666666667e+1;                                  
	
  C[0][4] = 0.0;
  C[1][4]= 10.0;
  C[2][4]= -7.0;
  C[3][4]= 12.0;
  C[4][4]= -33.0;
  C[5][4]= 130.0;
  C[6][4]= -691.0;
  
  u=z;                                                                       
  x=u;                                                                       
  a=abs(x); 
  
  if (k < 0 || k > 4 || k - int(k) !=0) 
    {                                           
      h=1e+278;                                                                     
      std::cout << "\tError: Psi(k,z) is defined only for integer k=0,1,2,3,4. Please correct k." 
		<< std::endl;  
    }
  else if (abs(imag(u)*Cunit) < delta && abs(x + nint(a)*1.0) < delta) 
    {
      h=1e+278;                                                                      
      std::cout << "\tError: Psi(k,z) called with argument z that is a non-possitive integer. " 
		<< std::endl;
    }
  else 
    {
      double k1=k+1;                                                                   
      if (real(x) < 0) u=-u;
      v=u;                                                                      
      h=0;
      if (a < 15) 
	{                                                      
	  h=pow(v,-k1);
	  for (int i = 1; i<=14-int(a); ++i)
	    {                                                    
	      v=v+1.0;                                                                   
	      h+= + pow(v,-k1);                                                             
	    }
	  v=v+1.0;
	}                                                                   
      r = pow(v,-2);                                                                 
      p=r*C[6][k];                                                               
      for( int i = 5; i>=1; --i)
	{                                                          
	  p=r*(C[i][k]+p);
	}
      h=sgn[k]*(fct[k+1]*h+(v*(fct[k]+p)+HF*fct[k+1])/pow(v,k1));                     
      if (k == 0) h+=log(v);                                                  
      
      if( real(x) < 0)
	{
	  v=PI*u;                                                                  
	  x=v;                                                                     
	  y=imag(v);                                                               
	  a=sin(real(x));                                                                
	  b=cos(real(x));                                                                
	  t=tanh(y);                   
	  dcomp num(b,-a*t);
	  dcomp den(a,b*t);
	  p=num/den;                                          
	  
	  if(k == 0)
	    {                                                       
	      h+= pow(u,-1)+PI*p;                                                           
	    }
	  else if(k == 1)
	    {                                                   
	      h=-h+pow(u,-2)+C1*(pow(p,2)+1.0); 
	    }
	  else if(k == 2)
	    {                                                   
	      h+=2.0*pow(u,-3)+C2*p*(pow(p,2)+1.0);
	    }
	  else if(k == 3)
	    {                                                   
	      r=pow(p,2);                                                                 
	      h=-h+6.0*pow(u,-4)+C3*((3.0*r+4.0)*r+1.0);                                           
	    }
	  else if(k == 4)
	    {                                                   
	      r=pow(p,2);                                                                 
	      h+=24.0*pow(u,-5)+C4*p*((3.0*r+5.0)*r+2.0);                                         
	    }
	}
    }
  return h;
}
