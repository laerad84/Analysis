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
Double_t E14WaveFitter::GetChisqNDF( double& chisqPed, double& chisqHead, double& chisqTail){
  
  Double_t chisq = m_FitFunc->GetChisquare() / m_FitFunc->GetNDF();
  chisqPed = m_ChisqNDF_Pedestal;
  chisqHead = m_ChisqNDF_Head;
  chisqTail = m_ChisqNDF_Tail;

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
  
  m_fFit       = false;
  m_ChisqNDF_Pedestal = (double)0xFFFF;
  m_ChisqNDF_Head     = (double)0xFFFF;
  m_ChisqNDF_Tail     = (double)0xFFFF;

}  

bool     E14WaveFitter::Fit( TGraph* gr , double minX, double maxX ){
  //if( !Approx( gr ) ){ return false; }
#ifdef DEBUG
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
#endif
  
  m_fFit = true;

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
      CalChisqNDF( gr );
      return true;
    }else{
      //std::cout<< "FitTime" << std::endl;      
      return false;
    }
  }else{
    int rst = gr->Fit( m_FitFunc, "Q", "", m_peakTime-120, m_peakTime+50 );
    if( rst == 0 ){
      FitTime( gr );
      CalChisqNDF( gr );
      return true;
    }else{
      //std::cout<< "FitTime" << std::endl;
      return false;
    }
  }  
}
bool     E14WaveFitter::CalChisqNDF( TGraph* gr ){
  if( !m_fFit ){ 
    m_ChisqNDF_Pedestal = (double)0xFFFF;
    m_ChisqNDF_Head     = (double)0xFFFF;
    m_ChisqNDF_Tail     = (double)0xFFFF;
    return false;
  }
  m_ChisqNDF_Pedestal = 0;
  m_ChisqNDF_Head     = 0;
  m_ChisqNDF_Tail     = 0;
  int nPointsPedestal = 0;
  int nPointsHead     = 0;
  int nPointsTail     = 0;
  double ChisqPedestal = 0.;
  double ChisqHead     = 0.;
  double ChisqTail     = 0.;

  double x;
  double y;
  double peakX = m_FitFunc->GetParameter(1);
  for( int i = 0; i< gr->GetN(); i++){
    x = gr->GetX()[i];
    y = gr->GetY()[i]; 

    if( x >= peakX - 120 &&
        x <= peakX - 96  ){
      ChisqPedestal += (m_FitFunc->Eval(x)-y)*(m_FitFunc->Eval(x) -y);
      nPointsPedestal++; 
    }else if( x > peakX - 96 &&
	      x <= peakX     ){
      ChisqHead     += (m_FitFunc->Eval(x)-y)*(m_FitFunc->Eval(x) -y);
      nPointsHead++; 
    }else if( x > peakX       &&
	      x <= peakX + 96 ){
      ChisqTail     += (m_FitFunc->Eval(x)-y)*(m_FitFunc->Eval(x) - y );
      nPointsTail++;
    }else{ 
      continue; 
    }
  }

  if( peakX < 96 || nPointsPedestal == 0 ){ 
    ChisqPedestal =(double)0xFFFF;
  }else{ 
    ChisqPedestal = ChisqPedestal/nPointsPedestal;
  }

  if( peakX < 0 && nPointsHead == 0 ){
    ChisqHead = ( double )0xFFFF;
  }else{ 
    ChisqHead = ChisqHead / nPointsHead;
  }

  if( peakX > gr->GetN()*8 || nPointsTail ==0 ){
    ChisqTail = ( double )0xFFFF;
  }else{
    ChisqTail = ChisqTail / nPointsTail;
  }

  m_ChisqNDF_Pedestal = ChisqPedestal;
  m_ChisqNDF_Head     = ChisqHead;
  m_ChisqNDF_Tail     = ChisqTail;
  if( m_ChisqNDF_Pedestal >= 0xFFFF ||
      m_ChisqNDF_Head     >= 0xFFFF ||
      m_ChisqNDF_Tail     >= 0xFFFF ){
    return false;
  }else{
    return true;
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
  std::cout << "FitChisq   :" 
	    << m_ChisqNDF_Pedestal << "\t" 
	    << m_ChisqNDF_Head     << "\t" 
	    << m_ChisqNDF_Tail     << "\n";
  return true;
}

TGraph* E14WaveFitter::GetDeltaFit( TGraph* gr ){
  if( !m_fFit ){
    std::cerr << "Please Execute Fit before get DeltaFit graph. " << std::endl;
    return (TGraph*)NULL;
  }
  
  TGraph* grTemp = new TGraph();
  double x;
  double y;
  for( int i = 0; i < gr->GetN(); i++){
    x = gr->GetX()[i];
    y = gr->GetY()[i];
    grTemp->SetPoint( i, x, y - m_FitFunc->Eval(x));
  }
  return grTemp;
}
