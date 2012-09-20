#ifndef E14WAVEFITTER__H__
#define E14WAVEFITTER__H__
#include "E14WaveFitterMinimum.h"
#include "TF1.h"

class E14WaveFitter : public E14WaveFitterMinimum {
 private:

  double    m_height;
  double    m_peakTime;
  double    m_HHTime;
  double    m_splTime;
  double    m_gnd;
  double    m_slope;

  int      m_pedsmpl;
  int      m_peakpoint;  
  double   m_gndrms;
  double   m_tempheight;
  double   m_tempgnd;

  int      m_PeakFlag;  
  TF1*     m_linfunc;


  double   m_PedLimitHigh;
  double   m_PedLimitLow;
  double   m_HighestPoint;
  double   m_LowestPoint;
  double   m_HighestTime;
  double   m_LowestTime;

  double   m_FitLimitLow;
  double   m_FitLimitHigh;

 public:

  E14WaveFitter(int SelectFunction = 0,int nPedestal = 7);
  virtual ~E14WaveFitter( void );

  virtual  void     SetParameter( Int_t parNum, Double_t value );
  virtual  void     SetParLimits( int parNum, double lowlimit, double highlimit);
  virtual  Double_t GetParameter( int parNum );
  virtual  Double_t GetConstantFraction();
  virtual  Double_t GetFitResult( );
  virtual  Double_t GetChisquare();
  virtual  Int_t    GetNDF();

  virtual  void     InitPar();
  virtual  int      CheckWaveform     ( TGraph* gr );
  virtual  bool     Fit               ( TGraph* gr ,double minX = 0, double maxX = 48*8);  
  virtual  bool     FitTime           ( TGraph* gr );
  virtual  bool     Approx            ( TGraph* gr );
  virtual  Double_t GetADC            ( TGraph* gr );
  virtual  Double_t GetADCPeak        ( TGraph* gr );
  virtual  Double_t GetFuncADC        ( TGraph* gr );
  virtual  Double_t GetFuncADCPeak    ( TGraph* gr );
  virtual  Double_t GetFitResultAspect( TGraph* gr );
  virtual  bool     ClearFunction     ( TGraph* gr );

 private:
  virtual  bool     GetLimits   ( TGraph* gr );
  virtual  bool     GetPeakPoint( TGraph* gr ,Double_t MinimumLimit = 0., Double_t MaximumLimit  = 384. );
  virtual  bool     GetPedestal ( TGraph* gr );
  virtual  bool     GetHeight   ( TGraph* gr );  
};


#endif
