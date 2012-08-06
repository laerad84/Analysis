#ifndef E14WAVEFITTER__H__
#include "E14WaveFitter.h"
#endif 

E14WaveFitter::E14WaveFitter(int SelectFunction, int nPedestal){
  m_FuncFlag = SelectFunction;
  m_pedsmpl = nPedestal;
  MakeFunction();
}
E14WaveFitter::~E14WaveFitter(){

}
bool E14WaveFitter::Fit( TGraph* gr , Double_t minX, Double_t maxX){
  if( !Approx( gr ) ){ return false; }
  int rst = gr->Fit( m_FitFunc, "QR");
  if( rst == -1 ){
    return false;
  }else{
    return true;
  }
}
bool E14WaveFitter::FitTime( TGraph* gr, Double_t minX, Double_t MaxX){

}

bool E14WaveFitter::Approx( TGraph* gr ){
  this->GetPeakPoint( gr );
  this->GetPedestal( gr );
  this->GetHeight( gr );
  //if( m_gndrms > 4 ){ return false; }   
  if( m_height < 8 ){ return false; }
  
  return true;
}

void E14WaveFitter::GetPeakPoint( TGraph* gr ){
  m_peakpoint  = 0;
  m_tempheight = 0;
  for( int ipoint = 0; ipoint < gr->GetN(); ipoint++){
    if( gr->GetY()[ipoint] > m_tempheight ){
      m_tempheight = gr->GetY()[ipoint];
      m_peakpoint  = ipoint;
    }
  }
}
     
void E14WaveFitter::GetPedestal( TGraph* gr ){
  double MinimumSigma = 1000000;
  int    MinimumPoint = 0;
  double localMaximum;
  double localSigma;
  double localMean;
  // Local Maximum RMS method /// 
  /******************************************************************************
  for( int ipoint = 0; ipoint < m_peakpoint - m_pedsmpl; ipoint++){
    localMaximum = 0;
    for( int subpoint = 0; subpoint < m_pedsmpl; subpoint++){
      if( localMaximum < gr->GetY()[ipoint + subpoint] ){
	localMaximum = gr->GetY()[ ipoint + subpoint ];
      }
    }
    localRMS = 0;
    for( int subpoint  =0; subpoint < m_pedsmpl; subpoint++){
      localRMS += TMath::Power(localMaximum - gr->GetY()[ipoint + subpoint],2);
    }
    localRMS = TMath::Sqrt( localRMS )/m_pedsmpl;
    if( localRMS < MinimumRMS ){
      MinimumRMS   = localRMS;
      MinimumPoint = ipoint;
    }
  }
  ******************************************************************************/
  /// Minimum RMS method ///
  /******************************************************************************/
  for( int ipoint = 0; ipoint < m_peakpoint - m_pedsmpl; ipoint++){
    localMean  = 0.;
    localSigma = 0.;
    for( int isubpoint  = 0; isubpoint < m_pedsmpl; isubpoint++){
      localMean += gr->GetY()[ipoint+ isubpoint];
      localSigma+= TMath::Power(gr->GetY()[ ipoint + isubpoint ], 2 );
    }
    localMean /= m_pedsmpl;
    localSigma = TMath::Sqrt(localSigma/m_pedsmpl - localMean*localMean );
    if( localSigma < MinimumSigma ){
      MinimumSigma = localSigma; 
      MinimumPoint = ipoint;
    }
  }
  /******************************************************************************/
  double gnd = 0;
  for( int ipoint = MinimumPoint; ipoint < MinimumPoint + m_pedsmpl; ipoint++){
    gnd += gr->GetY()[ipoint];
  }
  gnd = gnd / m_pedsmpl;
  m_gndrms = MinimumSigma;
  m_gnd    = gnd;
}

void E14WaveFitter::GetHeight( TGraph* gr ){ 
  m_height = gr->GetY()[ m_peakpoint ];
  m_height = m_height - m_gnd;
}
void E14WaveFitter::InitPar(){
  m_height   = 0.;

  m_peakTime   = 0.;
  m_HHTime     = 0.;
  m_splTime    = 0.; 

  m_gnd        = 0.;
  m_gndrms     = 0.;

  m_tempheight = 0.;
  m_tempgnd    = 0.; 
  m_peakpoint  = 0;
  
  for( int i = 0; i< 4000; i++){
    m_PeakID[i]   = -1;
    m_PeakFlag[i] = 0;
  } 
  m_nPeakFlag = 0;
  
}  
  
bool E14WaveFitter::CheckWaveform( TGraph* gr , int chID){
  if( m_nPeakFlag == 3999 ){
    std::cerr << "Flag Array is Full. Please use InitPar() start of event" << std::endl;
    return false;
  }
  GetPeakPoint(gr);
  GetPedestal(gr);
  GetHeight(gr);
  if( m_peakpoint < 15 ){ 
    m_nPeakFlag[m_nPeakFlag] |= 1;
  }
  if( m_gndrms >= 3. ){
    m_nPeakFlag[m_nPeakFlag] |= 2;
  }
  if(m_height <= 10. ){
    m_nPeakFlag[m_nPeakFlag] |= 4;
  }
  if( m_nPeakFlag[m_nPeakFlag] != 0){
    m_PeakID[m_nPeakFlag] = chID;   
    m_nPeakFlag++;
    return false; 
  }else{
    return true;
  }
}
