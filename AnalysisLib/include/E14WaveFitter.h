#ifndef E14WAVEFITTER__H__
#define E14WAVEFITTER__H__
#include "E14WaveFitterMinimum.h"
#include "E14WaveformAnalyzer.h"
#include "TF1.h"
class E14WaveFitter : public E14WaveFitterMinimum {
 private:

  double    m_height;
  double    m_peakTime;
  double    m_HHTime;
  double    m_splTime;
  double    m_gnd;
  double    m_slope;

  int      m_PeakFlag;  

  double   m_PedLimitHigh;
  double   m_PedLimitLow;
  double   m_HeightLimitLow;
  double   m_HeightLimitHigh;


 public:
  TF1*     m_linfunc;

  E14WaveFitter(int nPoints = 48);
  virtual ~E14WaveFitter( void );

  virtual  void     SetParameter( Int_t parNum, Double_t value );
  virtual  void     SetParLimits( int parNum, double lowlimit, double highlimit);
  virtual  Bool_t   SetParameter( E14WaveformAnalyzer* WaveAna);
  virtual  Double_t GetParameter( int parNum );
  virtual  Double_t GetConstantFraction();
  virtual  Double_t GetFitResult( );
  virtual  Double_t GetChisquare();
  virtual  Int_t    GetNDF();

  virtual  void     InitPar();
  virtual  bool     Fit               ( TGraph* gr ,double minX = 0, double maxX = 48*8);  
  virtual  bool     FitTime           ( TGraph* gr );
  virtual  bool     ClearFunction     ( TGraph* gr );
  virtual  bool     Draw( TGraph* gr );
};




#endif
