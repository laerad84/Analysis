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

 public:

  E14WaveFitter(int SelectFunction = 0,int nPedestal = 7);
  virtual ~E14WaveFitter( void );

  void     SetParameter( Int_t parNum, Double_t value );
  void     SetParLimits( int parNum, double lowlimit, double highlimit);
  Double_t GetParameter( int parNum );
  Double_t GetConstantFraction();
  Double_t GetFitResult( );
  Double_t GetChisquare();
  Int_t    GetNDF();

  void     InitPar();
  int      CheckWaveform( TGraph* gr );
  bool     Fit    ( TGraph* gr ,double minX = 0, double maxX = 48*8);  
  bool     FitTime( TGraph* gr );
  bool     Approx ( TGraph* gr );
  //bool     Clear  ( void );
  Double_t GetADC ( TGraph* gr );

 private:
  bool     GetLimits  ( TGraph* gr );
  bool     GetPeakPoint( TGraph* gr ,Double_t MinimumLimit = 0., Double_t MaximumLimit  = 384. );
  bool     GetPedestal ( TGraph* gr );
  bool     GetHeight   ( TGraph* gr );  
};


#endif
