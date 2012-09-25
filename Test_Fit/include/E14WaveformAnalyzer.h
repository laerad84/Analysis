#ifndef E14WAVEFORMANALYZER__H__
#define E14WAVEFORMANALYZER__H__

#include "TGraph.h"
#include "TMath.h"
#include "TROOT.h"
#include "TF1.h"
#include "TLine.h"
Double_t pol2func(Double_t *x, Double_t* par);

class E14WaveformAnalyzer {
 private:
  static const Int_t    m_CnHead    = 4;
  static const Int_t    m_CnTail    = 6; 
  static const Double_t m_TimeWidth = 8.;
  Int_t  m_nPoint;
  Bool_t m_fPedestal;
  Bool_t m_fHeight;
  Bool_t m_fMaximum;
  Bool_t m_fMinimum;
  Bool_t m_fWidth;
  Bool_t m_fmHead;
  Bool_t m_fmTail; 
  Bool_t m_fBoundary;
  Bool_t m_fADC; 

  Double_t m_MeanHead;
  Double_t m_MeanTail;
  Double_t m_SigmaHead; 
  Double_t m_SigmaTail;
  Double_t m_Pedestal;
  Double_t m_PedestalSigma;
  Double_t m_Height;
  Double_t m_PeakMaximum;
  Double_t m_PeakMinimum;
  Double_t m_TimeMaximum;
  Double_t m_TimeMinimum;
  Double_t m_BoundaryHead;
  Double_t m_BoundaryTail;
  Double_t m_Width; 
  Double_t m_ADC;

 public:
  TF1*     m_peakFunc;
  TGraph*  m_peakGraph;
  
  E14WaveformAnalyzer( int nPoints  = 48 );
  virtual ~E14WaveformAnalyzer();
  
  virtual Bool_t _GetMeanHead( Double_t* Waveform );
  virtual Bool_t _GetMeanTail( Double_t* Waveform );
  virtual Bool_t _GetMaximum ( Double_t* Waveform );
  virtual Bool_t _GetMinimum ( Double_t* Waveform );

  virtual Bool_t _GetPedestal( Double_t* Waveform );
  virtual Bool_t _GetWidth   ( Double_t* Waveform );
  virtual Bool_t _GetADC     ( Double_t* Waveform , Int_t StartTime = 0, Int_t EndTime = 48);
  
  virtual Bool_t _GetMeanHead( Double_t* Waveform , Double_t& MeanHead , Double_t& SigmaHead );
  virtual Bool_t _GetMeanTail( Double_t* Waveform , Double_t& MeanTail , Double_t& SigmaTail );  
  virtual Bool_t _GetMaximum ( Double_t* Waveform , Double_t& PeakMaximum   , Double_t& TimeMaximum   );
  virtual Bool_t _GetMinimum ( Double_t* Waveform , Double_t& PeakMinimum   , Double_t& TimeMinimum   );

  virtual Bool_t _GetPedestal( Double_t* Waveform , Double_t& Pedestal      , Double_t& PedestalSigma );
  virtual Bool_t _GetWidth   ( Double_t* Waveform , Double_t& Width);

  virtual Bool_t _GetADC     ( Double_t* Waveform , Double_t& ADC      , Int_t StartTime = 0, Int_t EndTime = 48);

  virtual void   _Clear();
  virtual void   _Draw(Double_t* Waveform);
  
};
#endif 
