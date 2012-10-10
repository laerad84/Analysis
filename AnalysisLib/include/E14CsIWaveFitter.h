#ifndef E14CSIWAVEFITTER__H__
#define E14CSIWAVEFITTER__H__
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include "TROOT.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"
#include "TSpline.h"


class E14WaveFitter {
 private:
  TF1*     m_FitFunc;
  TSpline3 m_spl;
  double   m_height;
  double   m_peakTime;
  double   m_HHTime;
  double   m_splTime;

  Double_t TemplateFunction( double *x, double *par );
  Double_t AsymmericGaussian( double *x, double *par );

 public:
  E14WaveFitter();
  virtual ~E14WaveFitter( void );
  TF1*     GetFunction(void){ return m_FitFunc;}
  void     SetWaveform( TSpline3* spl ){m_spl = spl;}
  void     SetParameter( int parNum, double value );
  void     SetParLimit ( int parNum, double lowlimit, double highlimit);
  Double_t GetParameter( int parNum );
  void     GetFitResult( );
  bool     Fit( TGraph* gr );  
  bool     Approx( TGraph* gr );
  bool     Clear( void );
  bool     FitResult();
};


#endif
