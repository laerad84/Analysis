#include "E14WaveFitterMinimum.h"
TSpline3* E14WaveFitterMinimum::m_spl = new TSpline3();
E14WaveFitterMinimum::E14WaveFitterMinimum(){
  //m_FuncFlag = SelectFunction;
  MakeFunction();
}

E14WaveFitterMinimum::~E14WaveFitterMinimum(){
  ;
}

int E14WaveFitterMinimum::MakeFunction(){
  m_FitFunc = new TF1("fitFuncTemplete",fTemplateFunction,
			0., 500., 3);
  return 1;
}    

int E14WaveFitterMinimum::Fit( TGraph* gr , char* goption, char* foption,
			       Axis_t xmin, Axis_t xmax){
  int rst = gr->Fit(this->m_FitFunc, goption, foption, xmin, xmax);
  return rst;
}

double E14WaveFitterMinimum::fTemplateFunction( double* x ,double* par){

  double t=x[0];
  double height = par[0];
  double mean   = par[1];
  double ped    = par[2];

  double t_fcn  = 0;
  if( t - mean  <= -150 || t- mean >= 300 ){ return ped;}
  //if( t -mean < -80 || t-mean >= 90 ){ return ped ;}
  if( t == 1./0 ){ return ped; }
  t_fcn         = height*E14WaveFitterMinimum::m_spl->Eval(t - mean) + ped;
  return t_fcn;
}

double E14WaveFitterMinimum::fAsymmetricGaussian(double* x, double* par){  
  double t=x[0];
  double height=par[0];
  double mean=par[1];
  double ped=par[2];
  double b=par[3];
  double a=par[4];
  double t_fcn=0.;  
  // asynmetric gaussian
  double sigma=a*(t-mean)+b;
  if ( sigma<0 ){ return ped; }
  t_fcn=height*TMath::Gaus(t,mean,sigma)+ped;  
  return t_fcn;
}

bool E14WaveFitterMinimum::Clear(){
  m_FitFunc->SetParameters(0,0,0);
  return true;
}
		    