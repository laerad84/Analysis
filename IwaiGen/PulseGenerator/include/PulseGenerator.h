#ifndef PULSEGENERATOR__H__
#define PULSEGENERATOR__H__

#include <iostream>
#include <string>
#include "TGraph.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TObject.h"

static const int    nSample   = 64;
static const int    fSeg      = 10; 
static const double pdfPar[5]={3.28113, 0.541526, -0.00582465, 4.53179e-05, -1.16776e-07};  
static const double corPar[5]={-0.197912, 0.0695365, 0.161864, 0.193783, 0.232639};
//static double pFuncArray[8*nSample*fSeg];
static double pFuncArray[2560];

class PulseGenerator {
 public:

  static const int    peMean    = 160;
  static const double absLY     = 12.7;
  static const double gndSgm    = 2.4; 
  static const double normF     = 1.652;
  static const double lambdaPMT = 14.; 

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Asymmetric Gaussian //   
  // p3-p7, valid in range [p1-8,p1+160] (p1:peak time)
  ////////////////////////////////////////////////////////////////////////////////////////////////////


  TF1* pePDF;
  TF1* pfPDF;
  TF1* Waveform;


  PulseGenerator();
  virtual ~PulseGenerator();
  void MakeFunction();
  void SetLightYield( double E_to_PE );  
  TF1* GetWaveform( double Energy, double SignalTime, double E_to_PE);
  TF1* GetWaveform( std::vector<double> EnergyArr, std::vector<double> SignalTimeArr, double E_to_PE);
  static double agaus( double* x, double* par);
  static double pulseFunc( double* x, double* par);
  static double Asymmetric( double* x, double* par);
  ClassDef( PulseGenerator, 1 )
};

#endif
