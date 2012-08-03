#ifndef E14WAVEFITTER__H__
#include "E14WaveFitterMimimum.h"


E14WaveFitterMinimum::E14WaveFitterMinimum(){
  //  m_FitFunc = new TF1("fitFunc",TemplateFunction,
  //0, 500, 3);
}

E14WaveFitterMinimum::~E14WaveFitterMinimum(){
  ;
}

int E14WaveFitterMinimum::Fit( TGraph* gr , char* goption, char* foption,
			       double xmin, double xmax){
  int rst = gr->Fit(this->m_FitFunc, goption, foption, xmin, xmax);
  return rst;
}

double E14WaveFitterMinimum::TemplateFunction( double* x ,double* par){
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
  t_fcn         = height*m_spl->Eval(t - mean) + ped;
  return t_fcn;

}

		    
