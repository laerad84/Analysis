#include "E14WaveFitterMinimum.h"
TSpline3* E14WaveFitterMinimum::m_spl = new TSpline3();
E14WaveFitterMinimum::E14WaveFitterMinimum(){
  m_spl = new TSpline3();
  m_FitFunc = new TF1("fitFunc",fTemplateFunction,
		      0, 500, 5);
}

E14WaveFitterMinimum::~E14WaveFitterMinimum(){
  ;
}

int E14WaveFitterMinimum::Fit( TGraph* gr , char* goption, char* foption,
			       double xmin, double xmax){
  int rst = gr->Fit(this->m_FitFunc, goption, foption, xmin, xmax);
  return rst;
}

double E14WaveFitterMinimum::fTemplateFunction( double* x ,double* par){
  double t=x[0];
  double height = par[0];
  double mean   = par[1];
  double ped    = par[2];
  /*
    double b      = par[2];
    double a      = par[3];
  */
  double t_fcn  = 0;
  //std::cout<< t << " : " << spl->Eval(t+mean) << std::endl;
  if( t - mean  < -100 || t- mean >150 ){ return ped;}
  t_fcn         = height* m_spl->Eval(t - mean) + ped;
  return t_fcn;
}

bool E14WaveFitterMinimum::Clear(){
  m_FitFunc->SetParameters(0,0,0,0,0);
}
		    
