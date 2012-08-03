#ifndef E14WAVEFITTER__H__
#include "E14WaveFitter.h"
#endif 

E14WaveFitter::E14WaveFitter(int nPedestal){
  m_pedsmpl = nPedestal;

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

bool E14WaveFitter::Approx( TGraph* gr ){
  this->GetPeakPoint( gr );
  this->GetPedestal( gr );
  this->GetHeight( gr );
  if( m_gndrms > 4 ){ return false; }   
  if( m_height < 8 ){ return false; }
  
  return true;
}

void E14WaveFitter::GetPeakPoint( TGraph* gr ){
  m_peakpoint =0;
  m_tempheight = 0;
  for( int ipoint  =0; ipoint < gr->GetN(); ipoint++){
    if( gr->GetY()[ipoint] > m_tempheight ){
      m_tempheight = gr->GetY()[ipoint];
      m_peakpoint  = ipoint;
    }
  }
}
     
void E14WaveFitter::GetPedestal( TGraph* gr ){
  double MinimumRMS   = 1000000;
  int    MinimumPoint = 0;
  double localMaximum;
  double localRMS;
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
  double gnd = 0;
  for( int ipoint = MinimumPoint; ipoint < MinimumPoint + m_pedsmpl; ipoint++){
    gnd += gr->GetY()[ipoint];
  }
  gnd = gnd / m_pedsmpl;
  m_gndrms = MinimumRMS;
  m_gnd = gnd;
}

void E14WaveFitter::GetHeight( TGraph* gr ){ 
  m_height = gr->GetY()[ m_peakpoint ];
  m_height = m_height - m_gnd;
}
