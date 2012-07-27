#ifndef Waveform_h
#define Waveform_h
#define WAVEFORM_VERSION "WAVEFORM v0.0.2_201005XX"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <TROOT.h>
#include <TObject.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMath.h>

const double f_sigma=2.;

double AsymmetricGaussian(double* x, double* par);
double LinearGaussian(double* x, double* par);

class WaveformFitter{
public:
  WaveformFitter(int Nsamples=64, bool fixed=false, int PedSmpl=8);
  virtual ~WaveformFitter(void);
  
  TF1* GetFunction(void){ return m_fitfunc; }
  void SetParameters(double width, double asymm){ m_width=width; m_asymm=asymm; }
  double GetParameter(int parNum);
  bool Approx(TGraph* gr);
  bool Fit(TGraph* gr);
  bool Clear();
protected:
  int    m_Nsamples;
  bool   m_fixed;
  TF1*   m_fitfunc;
  int    m_pedsmpl;
  
  double m_width;
  double m_asymm;
  
  ClassDef(WaveformFitter,1)
};
#endif
