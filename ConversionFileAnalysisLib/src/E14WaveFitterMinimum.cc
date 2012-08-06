#include "E14WaveFitterMinimum.h"
TSpline3* E14WaveFitterMinimum::m_spl = new TSpline3();
E14WaveFitterMinimum::E14WaveFitterMinimum(int SelectFunction){
  m_FuncFlag = SelectFunction;
  MakeFunction();
}

E14WaveFitterMinimum::~E14WaveFitterMinimum(){
  ;
}

int E14WaveFitterMinimum::MakeFunction(){
  if(m_FuncFlag == 0){
    m_FitFunc = new TF1("fitFuncTemplete",fTempleteFunction,
			0., 500., 3);
  }else{
    m_FitFunc = new TF1("fitFuncAsymmetric",fAsymmetricGaussian,
			0., 500., 5);
  }
  return m_FuncFlag;
}    
int E14WaveFitterMinimum::Fit( TGraph* gr , char* goption, char* foption,
			       Axis_t xmin, Axis_t xmax){
  int rst = gr->Fit(this->m_FitFunc, goption, foption, xmin, xmax);
  return rst;
}

double E14WaveFitterMinimum::fTempleteFunction( double* x ,double* par){
  double t=x[0];
  double height = par[0];
  double mean   = par[1];
  double ped    = par[2];

  double t_fcn  = 0;
  //std::cout<< t << " : " << spl->Eval(t+mean) << std::endl;
  if( t - mean  < -100 || t- mean >200 ){ return ped;}
  t_fcn         = height* m_spl->Eval(t - mean) + ped;
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
  if( m_FuncFlag == 0){
    m_FitFunc->SetParameters(0,0,0);
  }else{
    m_FitFunc->SetParameters(0,0,0,0,0);
  }
  return true;
}
		    
