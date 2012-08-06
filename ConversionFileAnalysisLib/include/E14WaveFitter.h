#ifndef E14WAVEFITTER__H__
#define E14WAVEFITTER__H__
#include "E14WaveFitterMinimum.h"

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

  int      m_nPeakFlag;
  int      m_PeakID[4000];
  int      m_PeakFlag[4000];
  
 public:

  E14WaveFitter(int SelectFunction = 0,int nPedestal = 7);
  virtual ~E14WaveFitter( void );

  void     SetParameter( int parNum, double value );
  void     SetParLimit ( int parNum, double lowlimit, double highlimit);
  Double_t GetParameter( int parNum );
  void     GetFitResult( );

  void     InitPar();
  bool     CheckWaveform( TGraph* gr , int chID  = -1);
  bool     Fit   ( TGraph* gr ,double minX = 0, double maxX = 48*8);  
  bool     Approx( TGraph* gr );
  bool     Clear ( void );

 private:
  void     GetPeakPoint( TGraph* gr );
  void     GetPedestal ( TGraph* gr );
  void     GetHeight   ( TGraph* gr );  
  
};


#endif