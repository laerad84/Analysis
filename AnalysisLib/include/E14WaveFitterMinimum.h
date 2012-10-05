#ifndef E14WAVEFITTERMINIMUM__H__
#define E14WAVEFITTERMINIMUM__H__

#include <iostream>
#include <vector>
#include <string>

#include "TROOT.h"
#include "TF1.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TMath.h"

class E14WaveFitterMinimum {
 public:

  static TSpline3* m_spl;
  static double    fTemplateFunction( double* x, double* par );
  static double    fTemplateFunctionDoublePeak( double* x, double* par);
  static double    fAsymmetricGaussian( double* x, double* par );

  TF1* m_FitFunc;
  TF1* m_FitFuncDoublePeak;

  //int  m_FuncFlag;
  E14WaveFitterMinimum();
  virtual ~E14WaveFitterMinimum();

  virtual int   MakeFunction();
  virtual TF1*  GetFunction( void ){ return m_FitFunc; }
  virtual void  SetWaveform( TSpline3* spl ){ E14WaveFitterMinimum::m_spl = spl; }
  virtual int   Fit(TGraph* gr, char* goption="", char* foption="", Axis_t xmin = 0, Axis_t xmax = 0);
  virtual bool  Clear(); 

};

#endif
