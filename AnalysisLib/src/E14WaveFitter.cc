#include "E14WaveFitter.h"

E14WaveFitter::E14WaveFitter(int SelectFunction, int nPedestal, int nPoints){
  //std::cout<< __FUNCTION__ << std::endl; 
  m_FuncFlag = SelectFunction;
  m_pedsmpl = nPedestal;
  MakeFunction();

  m_PeakFunc = new TF1("m_PeakFunc",Pol2Function, 0, 500., 3);
  m_PeakGr   = new TGraph();

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
Double_t E14WaveFitter::GetChisquare(){
  Double_t chisq = m_FitFunc->GetChisquare();
  return chisq;
}
Int_t  E14WaveFitter::GetNDF(){
  Int_t ndf = m_FitFunc->GetNDF();
  return ndf;
}
void E14WaveFitter::InitPar(){
  /*
  if( m_linfunc != NULL ){
    m_linfunc->Delete();
  }
  if( m_FitFunc != NULL ){
    m_FitFunc->Delete();
  }
  */
  m_height   = 0.;

  m_peakTime   = 0.;
  m_HHTime     = 0.;
  m_splTime    = 0.; 
  m_slope      = 0.;

  m_gnd        = 0.;
  m_gndrms     = 0.;

  m_tempheight = 0.;
  m_tempgnd    = 0.; 
  m_peakpoint  = 0;
  
  m_PeakFlag   = 0;

  m_PedLimitHigh = 0;
  m_PedLimitLow  = 0;
  m_HighestPoint = 0;
  m_LowestPoint  = 0;
  m_HighestTime  = 0;
  m_LowestTime   = 0;

  m_FitLimitLow  = 0;
  m_FitLimitHigh = 0;
}  
  
int  E14WaveFitter::CheckWaveform( TGraph* gr ){
  
  GetPeakPoint(gr);
  GetPedestal(gr);
  GetHeight(gr);

  //std::cout << typeid(*this).name() << ":" << __FUNCTION__ << std::endl;
  //std::cout <<  m_height << " : " << m_peakTime << " : " << m_gnd << std::endl; 
  if( m_peakpoint < 8 ){ 
    m_PeakFlag |= 1;
  }
  if( m_gndrms >= 5. ){
    m_PeakFlag |= 2;
  }
  if(m_height <= 10. ){
    m_PeakFlag |= 4;
  }

  return m_PeakFlag;
}

bool E14WaveFitter::Fit( TGraph* gr , double minX, double maxX ){
  //if( !Approx( gr ) ){ return false; }
#ifdef DEBUG
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
#endif
  if( CheckWaveform( gr ) != 0  ){ return false; }  

  m_FitFunc->SetParameter(2,m_gnd);
  m_FitFunc->SetParLimits(2,m_PedLimitLow,m_PedLimitHigh);
  m_FitFunc->SetParameter(0,m_height);
  m_FitFunc->SetParLimits(0,m_height*0.8, m_height*2);
  m_FitFunc->SetParameter(1,m_peakTime);
  m_FitFunc->SetParLimits(1,m_peakTime-16,m_peakTime+16);

  //std::cout <<m_peakTime << std::endl;
  int rst  = gr->Fit( m_FitFunc, "Q","",m_peakTime-150,m_peakTime+35);
  if( rst == 0 ){
    FitTime( gr );
    return true;
  }else{
    //std::cout<< "FitTime" << std::endl;
   return false;
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
  m_slope = m_linfunc->GetParameter(1);
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

bool E14WaveFitter::GetLimits( TGraph* gr ){
  m_HighestPoint = 0;
  m_LowestPoint  = 20000;
  for( int ipoint = 0; ipoint < gr->GetN(); ipoint++ ){
    double point = gr->GetY()[ipoint];
    double t     = gr->GetX()[ipoint]; 
    if( m_HighestPoint < point ) { m_HighestPoint = point;m_HighestTime =t; }
    if( m_LowestPoint  > point ) { m_LowestPoint  = point;m_LowestTime  =t; }    
  }
  // Check graph is exist or not (./.) 
  if( m_HighestPoint <= m_LowestPoint ){ return false; } 
  if( m_HighestTime  == m_LowestTime ){ return false; }
  return true;
}

bool E14WaveFitter::GetPeakPoint( TGraph* gr ,Double_t MinimumLimit, Double_t MaximumLimit){
  m_peakpoint  = 0;
  m_tempheight = 0;
  if( MinimumLimit > MaximumLimit ){ return false;}
  /// Maximum Point method (./.) /// 
  /**************************************************************************************
  for( int ipoint = 0; ipoint < gr->GetN(); ipoint++){
    if( gr->GetX()[ipoint] <= MinimumLimit || gr->GetX()[ipoint] > MaximumLimit ){
      continue;
    }
    if( gr->GetY()[ipoint] > m_tempheight ){
      m_tempheight = gr->GetY()[ipoint];
      m_peakpoint  = ipoint;
    }
  }
  **************************************************************************************/
  /// Maximum in 3point method (./.) ///
  //**************************************************************************************//
  /*
  double sum3heightest       = 0;
  double sum3heightestpoint  = 0;
  double sum3                = 0;
  double sum3Point           = 0;
  for( int ipoint = 1; ipoint < gr->GetN()-1; ipoint++){
    if( gr->GetX()[ipoint] < MinimumLimit || gr->GetX()[ipoint]> MaximumLimit ){
      continue;
    }

    sum3      = 0.;
    sum3Point = 0.;

    sum3 += gr->GetY()[ipoint - 1]; 
    sum3 += gr->GetY()[ipoint    ]; 
    sum3 += gr->GetY()[ipoint + 1]; 

    sum3Point += gr->GetY()[ipoint -1]*gr->GetX()[ipoint-1];
    sum3Point += gr->GetY()[ipoint   ]*gr->GetX()[ipoint  ];
    sum3Point += gr->GetY()[ipoint +1]*gr->GetX()[ipoint+1];
    sum3Point /= sum3;

    if( sum3 > sum3heightest ){ 
      sum3heightest      = sum3;
      sum3heightestpoint = sum3Point;
      m_peakpoint        = ipoint;
      m_tempheight       = gr->GetY()[ipoint];
    }    
  }
  */

  const int nSumPoint   = 6; 
  int    PeakFrontPoint = 0;
  double HighestSum     = 0.;
  double HighestPoint   = 0.;
  double Sum            = 0.;
  double SumPoint       = 0.;
  
  for( int ipoint = 0; ipoint < gr->GetN() - (nSumPoint-1) ; ipoint++){
    if( gr->GetX()[ipoint]  < MinimumLimit || gr->GetX()[ipoint] > MaximumLimit ){
      continue; 
    }
    Sum      = 0.;
    SumPoint = 0.;    
    for( int isubPoint = 0; isubPoint < nSumPoint -1 ; isubPoint++){
      Sum     += gr->GetY()[ipoint + isubPoint];
      SumPoint+= gr->GetY()[ipoint + isubPoint]*gr->GetX()[ipoint + isubPoint];
    }
    SumPoint  = SumPoint / Sum; 
    int PeakPoint  = 0; 
    if( nSumPoint%2 == 0){
      PeakPoint = ipoint + nSumPoint/2 -1;
    }else{
      PeakPoint = ipoint + nSumPoint/2;
    }
    if( gr->GetY()[ ipoint ] > gr->GetY()[ PeakPoint ]){ // front peak exception // 
      Sum      = 0.;
      SumPoint = 0.;
      continue; 
    }

    if( Sum > HighestSum ){
      PeakFrontPoint = ipoint;
      HeighestSum    = Sum;
      HeighestPoint  = SumPoint;
      //Test version of peak search // 
      /*************************************************************************/
      if( (nSumPoint % 2) == 0 ){ 
	m_peakpoint  = (ipoint + nSumPoint/2 + 0.5);
	if( gr->GetY()[ipoint+nSumPoint/2] > gr->GetY()[ipoint+nSumPoint/2 +1]){
	  m_tempheight = gr->GetY()[ipoint+nSumPoint/2];
	}else{
	  m_tempheight = gr->GetY()[ipoint+nSumPoint/2+1];
	}
      }else{
	m_peakpoint = ( ipoint + nSumPoint/2 ); 
	m_tempheight = gr->GetY()[ipoint+nSumPoint/2];
      } 
      /*************************************************************************/      
    }    
  }
  m_peakTime = m_peakpoint*8;

  if( Sum < 1 ){    
    m_peakTime = m_peakpoint*8;
    return false;
  }

  m_PeakGr->Set(0);
  for( int ipoint  = PeakFrontPoint; ipoint < PeakFrontPoint+ nSumPoint; ipoint++){
    m_PeakGr->SetPoint( ipoint-PeakFrontPoint,
			gr->GetX()[ipoint],
			TMath::Log( gr->GetY()[ipoint] ));
  }
  m_PeakGr->Fit( m_PeakFunc, "Q", "", HeightestPoint*8-nSumPoint*4, HeightestPoint*8+nSumPoint*4);
  m_peakTime = m_PeakFunc->GetParameter(2);

  //**************************************************************************************//
  return true;
}
     
bool E14WaveFitter::GetPedestal( TGraph* gr ){
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
  /******************************************************************************
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
  double gnd = 0;
  for( int ipoint = MinimumPoint; ipoint < MinimumPoint + m_pedsmpl; ipoint++){
    gnd += gr->GetY()[ipoint];
  }
  gnd = gnd / m_pedsmpl;
  m_gndrms = MinimumSigma;
  m_gnd    = gnd;
  ******************************************************************************/
  /// Set RMS Region and Minimum value method ///
  ///******************************************************************************
  // All Csi Output ///
  double minimumPeak = 20000;

  m_PedLimitHigh = 0;
  m_PedLimitLow  = 20000;
  
  for( int ipoint  = 0; ipoint < m_peakpoint - 8; ipoint++ ){    
    if( minimumPeak > gr->GetY()[ipoint] ){
      minimumPeak  = gr->GetY()[ipoint];
      MinimumPoint = ipoint;      
    }
    if( m_PedLimitHigh < gr->GetY()[ipoint] ){
      m_PedLimitHigh = gr->GetY()[ipoint];
    }
    if( m_PedLimitLow > gr->GetY()[ipoint] ){
      m_PedLimitLow = gr->GetY()[ipoint];
    }
  }
  if( m_PedLimitHigh - m_PedLimitLow > 10 ){
    m_PedLimitHigh  = m_PedLimitLow + 10;
    m_gnd = m_PedLimitLow;
    return false;
  }
  m_gnd = m_PedLimitLow;
  //******************************************************************************/
  return true;
}

bool E14WaveFitter::GetHeight( TGraph* gr ){ 
  m_height = gr->GetY()[ m_peakpoint ];
  m_height = m_height - m_gnd;
  return true;
}

Double_t E14WaveFitter::GetADC( TGraph* gr ){
  Double_t ADC  = 0.;
  if( m_peakpoint >= 8 && m_PeakFlag == 0){
    for( int ipoint  = 0; ipoint < gr->GetN() ; ipoint++){
      if( TMath::Abs(gr->GetX()[ipoint] - m_peakTime) < 50 ){      
	ADC += gr->GetY()[ipoint] - m_gnd;
      }
    }
    return ADC;
  }else{
    return 0.;
  }
}

Double_t E14WaveFitter::GetADCPeak( TGraph* gr ){
  Double_t ADC = 0;
  m_FitLimitLow = m_FitFunc->GetParameter(1) - 80;
  m_FitLimitHigh= m_FitFunc->GetParameter(1) + 24;
  
  if( m_peakpoint >= 8 && m_PeakFlag == 0 ){
    for( int ipoint = 0 ; ipoint < gr->GetN(); ipoint++){
      if( gr->GetX()[ipoint] > m_FitLimitLow && gr->GetX()[ipoint] < m_FitLimitHigh ){
	ADC += gr->GetY()[ipoint] - m_gnd;
      }
    }
    return ADC;
  }else{
    return 0;
  }
}

Double_t E14WaveFitter::GetFuncADC( TGraph* gr ){
  Double_t ADC = 0;
  if( m_peakpoint >= 8 && m_PeakFlag  == 0 ){
    for( int ipoint = 0; ipoint < gr->GetN() ; ipoint++ ){      
      ADC += m_FitFunc->Eval( gr->GetX()[ipoint] ) - m_gnd;
    }
  }else{
    return 0.;
  }
  return ADC; 
}

Double_t E14WaveFitter::GetFuncADCPeak( TGraph* gr ){
  Double_t ADC = 0;
  m_FitLimitLow = m_FitFunc->GetParameter(1) - 80;
  m_FitLimitHigh= m_FitFunc->GetParameter(1) + 24;
  
  if( m_peakpoint >= 8 && m_PeakFlag == 0 ){
    for( int ipoint = 0 ; ipoint < gr->GetN(); ipoint++){
      if( gr->GetX()[ipoint] > m_FitLimitLow && gr->GetX()[ipoint] < m_FitLimitHigh ){
	ADC += m_FitFunc->Eval( gr->GetX()[ipoint] ) - m_gnd;
      }
    }
    return ADC;
  }else{
    return 0;
  }
}


Double_t E14WaveFitter::GetFitResultAspect( TGraph* gr){
  
  Double_t  fOut = this->GetFuncADC( gr );
  Double_t  gOut = this->GetFuncADC( gr ); 
  
  if( fOut <= 0 || gOut <= 0 ){
    return 0;
  }else{
    return fOut/gOut; 
  }
}

bool E14WaveFitter::ClearFunction( TGraph* gr ){
  gr->GetListOfFunctions()->Delete();
  return true;
}

Double_t Pol2Function( Double_t *x, Double_t *par ){
  Double_t xpos  = x[0];
  Double_t Sigma = par[1];
  Double_t Mean  = par[2];
  Double_t rtn   = par[0] -1*(x[0]-par[2])*(x[0]-par[2])/(2*par[1]*par[1]);
  return rtn; 
}
