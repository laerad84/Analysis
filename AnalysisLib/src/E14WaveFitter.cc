#include "E14WaveFitter.h"

E14WaveFitter::E14WaveFitter(int nPoints){
  //std::cout<< __FUNCTION__ << std::endl; 
  MakeFunction();

}
E14WaveFitter::~E14WaveFitter(){
  ;
}
void     E14WaveFitter::SetParameter( Int_t parNum , Double_t value ){
  m_FitFunc->SetParameter( parNum, value );
}
void     E14WaveFitter::SetParLimits( int parNum, double lowlimit, double highlimit ){
  m_FitFunc->SetParLimits( parNum, lowlimit, highlimit );
}
Bool_t   E14WaveFitter::SetParameter( E14WaveformAnalyzer* WaveAna ){
  WaveAna->GetParameters( m_gnd, m_height, m_peakTime);
  return true; 
}
Double_t E14WaveFitter::GetParameter( int parNum ){
  Double_t rst = m_FitFunc->GetParameter( parNum );
  return rst;
}
Double_t E14WaveFitter::GetConstantFraction(){
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
Int_t    E14WaveFitter::GetNDF(){
  Int_t ndf = m_FitFunc->GetNDF();
  return ndf;
}
void     E14WaveFitter::InitPar(){


  m_height   = 0.;

  m_peakTime   = 0.;
  m_HHTime     = 0.;
  m_splTime    = 0.; 
  m_slope      = 0.;

  m_gnd        = 0.;

}  
bool     E14WaveFitter::Fit( TGraph* gr , double minX, double maxX ){
  //if( !Approx( gr ) ){ return false; }
#ifdef DEBUG
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
#endif

  if( m_height*0.1 < 4 ){
    if( m_height -4 > 0){ m_HeightLimitLow = m_height -4; }
    else{ m_HeightLimitLow = 0; }
  }else{
    m_HeightLimitLow = m_height*0.9;
  }
  if( m_height*0.5 < 4){
    m_HeightLimitHigh = m_height +4;
  }else{
    m_HeightLimitHigh = m_height*1.5;
  }

  m_PedLimitLow  = m_gnd -4;// 4 is 2*ped RMS // 
  m_PedLimitHigh = m_gnd +4;

  m_FitFunc->SetParameter(2,m_gnd);
  m_FitFunc->SetParLimits(2,m_PedLimitLow,m_PedLimitHigh);
  m_FitFunc->SetParameter(0,m_height);
  m_FitFunc->SetParLimits(0,m_HeightLimitLow, m_HeightLimitHigh);
  m_FitFunc->SetParameter(1,m_peakTime);
  m_FitFunc->SetParLimits(1,m_peakTime-16,m_peakTime+16);

  //std::cout <<m_peakTime << std::endl;
  if( minX != 0 || maxX != 48*8 ){
    int rst  = gr->Fit( m_FitFunc, "Q","", minX, maxX);
    if( rst == 0 ){
      FitTime( gr );
      return true;
    }else{
      //std::cout<< "FitTime" << std::endl;
      return false;
    }
  }else{
    int rst = gr->Fit( m_FitFunc, "Q", "", m_peakTime-120, m_peakTime+50 );
    if( rst == 0 ){
      FitTime( gr );
      return true;
    }else{
      //std::cout<< "FitTime" << std::endl;
      return false;
    }
  }  
}
bool     E14WaveFitter::FitTime( TGraph* gr){
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
bool     E14WaveFitter::ClearFunction( TGraph* gr ){
  gr->GetListOfFunctions()->Delete();
  return true;
}
bool     E14WaveFitter::Draw(TGraph* gr){
  gr->SetMarkerStyle(6);
  gr->Draw("AP");
  m_FitFunc->SetLineColor(2);
  m_FitFunc->SetLineWidth(2);
  m_linfunc->SetLineColor(4);
  m_linfunc->SetLineWidth(2);
  m_FitFunc->Draw("same");
  m_linfunc->Draw("same");
  std::cout << "////////////////////////////////////////////////////\n";
  std::cout << "PeakHeight :" << m_FitFunc->GetParameter(0) << "\n";
  std::cout << "PeakTime   :" << m_FitFunc->GetParameter(1) << "\n";
  std::cout << "Pedestal   :" << m_FitFunc->GetParameter(2) << "\n"; 
  std::cout << "FitResult  :" << GetFitResult() << "\n";

  

  
  return true;
}
