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
  static double fTempleteFunction( double* x, double* par );
  static double fAsymmetricGaussian( double* x, double* par );
  TF1* m_FitFunc;
  int  m_FuncFlag;
  E14WaveFitterMinimum(int SelectFunction = 0);
  virtual ~E14WaveFitterMinimum();
  int MakeFunction();
  TF1* GetFunction( void ){ return m_FitFunc; }
  void SetWaveform( TSpline3* spl ){ E14WaveFitterMinimum::m_spl = spl; }
  int  Fit(TGraph* gr, char* goption="", char* foption="", Axis_t xmin = 0, Axis_t xmax = 0);
  bool Clear(); 
};

#endif
