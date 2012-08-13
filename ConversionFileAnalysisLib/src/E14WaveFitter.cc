#include "E14WaveFitter.h"

E14WaveFitter::E14WaveFitter(int SelectFunction, int nPedestal){
  std::cout<< __FUNCTION__ << std::endl; 
  m_FuncFlag = SelectFunction;
  m_pedsmpl = nPedestal;
  MakeFunction();
}
E14WaveFitter::~E14WaveFitter(){
  ;
}

void E14WaveFitter::SetParameter( Int_t parNum , Double_t value ){
  m_FitFunc->SetParameter( parNum, value );
}

void E14WaveFitter::SetParLimits( int parNum, double lowlimit, double highlimit ){
  m_FitFunc->SetParLimits( parNum, lowlimit, highlimit );
}

Double_t E14WaveFitter::GetParameter( int parNum ){
  Double_t rst = m_FitFunc->GetParameter( parNum );
  return rst;
}
Double_t E14WaveFitter::GetConstantFraction( ){
  return m_HHTime;
}

Double_t E14WaveFitter::GetFitResult(){
  Double_t chisq = m_FitFunc->GetChisquare();
  Int_t    ndf   = m_FitFunc->GetNDF();
  if( ndf == 0) { return -1; }
  return chisq/ndf; 
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
  
  m_PeakFlag   = 0;
}  
  
int  E14WaveFitter::CheckWaveform( TGraph* gr ){
  
  GetPeakPoint(gr);
  GetPedestal(gr);
  GetHeight(gr);

  //std::cout << typeid(*this).name() << ":" << __FUNCTION__ << std::endl;
  //std::cout <<  m_height << " : " << m_peakTime << " : " << m_gnd << std::endl; 
  if( m_peakpoint < 15 ){ 
    m_PeakFlag |= 1;
  }
  if( m_gndrms >= 3. ){
    m_PeakFlag |= 2;
  }
  if(m_height <= 10. ){
    m_PeakFlag |= 4;
  }
  return m_PeakFlag;
}

bool E14WaveFitter::Fit( TGraph* gr , double minX, double maxX ){
  //if( !Approx( gr ) ){ return false; }
#ifdef ALALYSIS_DEBUG
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
#endif
  if( CheckWaveform( gr ) != 0  ){ return false; }  

  m_FitFunc->SetParameter(2,m_gnd);
  m_FitFunc->SetParLimits(2,m_gnd,m_gnd);
  m_FitFunc->SetParameter(0,m_height);
  m_FitFunc->SetParLimits(0,m_height*0.9, m_height*2);
  m_FitFunc->SetParameter(1,m_peakTime);
  m_FitFunc->SetParLimits(1,m_peakTime-32,m_peakTime+32);
  int rst  = gr->Fit( m_FitFunc, "","",minX,maxX);
  if( rst == 0 ){
    return false;
  }else{
    //std::cout<< "FitTime" << std::endl;
    //FitTime( gr );
   return true;
  }
}

bool E14WaveFitter::FitTime( TGraph* gr){
  m_splTime = m_FitFunc->GetX( m_FitFunc->GetParameter(0)*0.5 + m_FitFunc->GetParameter(2),
			       m_FitFunc->GetParameter(1)-60,   m_FitFunc->GetParameter(1));
  //std::cout<< m_splTime << std::endl; 
  gr->Fit( "pol1","Q","",m_splTime - 16, m_splTime + 16 );
  m_linfunc = gr->GetFunction("pol1");
  double slope  = m_linfunc->GetParameter(1);
  double offset = m_linfunc->GetParameter(0);
  double xpos   = ((m_height*0.5 +m_gnd) - offset)/slope;  
  if( m_linfunc->GetParameter(1) < 0 ){ m_HHTime = -1.; }
  if( abs( xpos - m_splTime ) > 8 ){ m_HHTime = -1.; }
  //m_linfunc->Clear();
  m_HHTime = xpos;   
  if( m_HHTime > 0 ){
    return true;
  }else{
    return false;
  }
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
  m_peakTime = m_peakpoint*8;
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
