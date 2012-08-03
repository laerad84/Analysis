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
  static double fTemplateFunction( double* x, double* par );
  static double fAsymmetricGaussian( double* x, double* par );
  TF1* m_FitFunc;
 
  E14WaveFitterMinimum();
  virtual ~E14WaveFitterMinimum();
  
  TF1* GetFunction( void ){ return m_FitFunc; }
  void SetWaveform( TSpline3* spl ){ m_spl = spl; }
  int  Fit(TGraph* gr, char* goption, char* foption, double xmin, double xmax);
  bool Clear(); 
};

#endif
