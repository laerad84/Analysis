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
  Bool_t m_fSlope;
  Bool_t m_fOverflow;
  Bool_t m_fUnderflow;
  Bool_t m_fWaveformAnalyzed;
  //// Bit Information
  //// Underflow     = 0; less than 1
  //// Overflow      = 1; bigger than 16000
  //// Peak in Head  = 2; MeanHead - Pedestal > 10;
  //// Peak in Tail  = 3; MeanTail - Pedestal > 10;
  //// Width         = 4; Width > 100; 
  //// Peak Position = 5; Peak Position > 250;
  Int_t  m_WaveformState;

  Double_t m_MeanHead;
  Double_t m_MeanTail;
  Double_t m_SigmaHead; 
  Double_t m_SigmaTail;
  Double_t m_Pedestal;
  Double_t m_PedestalSigma;
  Double_t m_Height;
  Double_t m_PeakMaximum;
  Double_t m_MinimumSigma;

  Double_t m_AllPointMaximum;
  Double_t m_AllPointMinimum;

  Int_t    m_AllPointMin;
  Int_t    m_AllPointMax;

  Double_t m_PeakMinimum;
  Double_t m_TimeMaximum;
  Double_t m_TimeMinimum;
  Double_t m_BoundaryHead;
  Double_t m_BoundaryTail;
  Double_t m_Width; 
  Double_t m_ADC;

  Int_t    m_SlopeStart; 
  Double_t m_SlopeDelta;

  Double_t m_PedestalFront;
  Double_t m_PedestalRear;
  Double_t m_ChisqNDFPedestalFront;
  Double_t m_ChisqNDFPedestalRear;
  
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
  virtual Bool_t _CheckOVF   ( Double_t* Waveform );
  virtual Bool_t _GetMaximumWithRange( Double_t* Waveform , Double_t& Maximum, Double_t& MaxTime, Int_t Init_Point = 0, Int_t End_Point = 48 );
  virtual Bool_t _GetMeanHead( Double_t* Waveform , Double_t& MeanHead , Double_t& SigmaHead );
  virtual Bool_t _GetMeanTail( Double_t* Waveform , Double_t& MeanTail , Double_t& SigmaTail );  
  virtual Bool_t _GetMaximum ( Double_t* Waveform , Double_t& PeakMaximum   , Double_t& TimeMaximum   );
  virtual Bool_t _GetMinimum ( Double_t* Waveform , Double_t& PeakMinimum   , Double_t& TimeMinimum   );
  virtual Bool_t _GetPedestal( Double_t* Waveform , Double_t& Pedestal      , Double_t& PedestalSigma );
  virtual Bool_t _GetWidth   ( Double_t* Waveform , Double_t& BoundaryHead, Double_t& BoundaryTail);
  virtual Bool_t _GetADC     ( Double_t* Waveform , Double_t& ADC      , Int_t StartTime = 0, Int_t EndTime = 48);
  virtual Bool_t _GetSumSlope( Double_t* Waveform );
  virtual Bool_t _GetSumSlope( Double_t* Waveform , Int_t& StartPoint , Double_t& DeltaMaximum);
  virtual void   _Clear();
  virtual void     Draw(Double_t* Waveform);  
  virtual Double_t GetSlopeDelta( Double_t* Waveform, Int_t StartPoint ); 
  virtual Int_t    AnalyzeWaveform( Double_t* Waveform );
  virtual Bool_t   IsAnalyzed() const { return m_fWaveformAnalyzed; }
  virtual Bool_t   GetParameters( Double_t& Pedestal, Double_t& Height, Double_t& PeakTime ) const;
  /*
  virtual Double_t   GetPeakPointMaximum const () { return m_PeakPointMaximum;}
  virtual Double_t   GetMeanHead const () { return m_MeanHead; }  
  virtual Double_t   GetMeanTail const () { return m_MeanTail; }
  virtual Double_t   GetPedestal const () { return m_Pedestal; }
  virtual Double_t   GetWidth    const () { return m_BoundaryTail-m_BoundaryHead; }
  virtual Double_t   GetADC      const () { return m_ADC; }
  */
  
};
#endif 
