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
Double_t E14WaveFitter::GetChisqNDF( double& chisqPed, double& chisqFront, double& chisqRear, double& chisqTail){
  
  Double_t chisq = m_FitFunc->GetChisquare() / m_FitFunc->GetNDF();
  chisqPed       = m_ChisqNDF_Pedestal;
  chisqFront     = m_ChisqNDF_Front;
  chisqRear      = m_ChisqNDF_Rear;
  chisqTail      = m_ChisqNDF_Tail;
  
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
  m_ChisqNDF_Front    = (double)0xFFFF;
  m_ChisqNDF_Rear     = (double)0xFFFF;
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
    m_ChisqNDF_Front    = (double)0xFFFF;
    m_ChisqNDF_Rear     = (double)0xFFFF;
    m_ChisqNDF_Tail     = (double)0xFFFF;
    return false;
  }
  m_ChisqNDF_Pedestal = 0;
  m_ChisqNDF_Front    = 0;
  m_ChisqNDF_Rear     = 0;
  m_ChisqNDF_Tail     = 0;
  int nPointsPedestal = 0;
  int nPointsFront    = 0;
  int nPointsRear     = 0;
  int nPointsTail     = 0;

  double ChisqPedestal = 0.;
  double ChisqFront    = 0.;
  double ChisqRear     = 0.;
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
      ChisqFront     += (m_FitFunc->Eval(x)-y)*(m_FitFunc->Eval(x) -y);
      nPointsFront++; 
    }else if( x > peakX       &&
	      x <= peakX + 96 ){
      ChisqRear     += (m_FitFunc->Eval(x)-y)*(m_FitFunc->Eval(x) - y );
      nPointsRear++;
    }else if( x > peakX + 96 && 
	      x <= peakX+150 ){
      ChisqTail     += (m_FitFunc->Eval(x)-y)*(m_FitFunc->Eval(x) - y);
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

  if( peakX < 0 || nPointsFront == 0 ){
    ChisqFront = ( double )0xFFFF;
  }else{ 
    ChisqFront = ChisqFront / nPointsFront;
  }

  if( peakX > gr->GetN()*8 || nPointsRear ==0 ){
    ChisqRear = ( double )0xFFFF;
  }else{
    ChisqRear = ChisqRear / nPointsRear;
  }

  if( peakX + 96 > gr->GetN()*8 || nPointsTail ==0 ){
    ChisqTail = ( double )0xFFFF;
  }else{
    ChisqTail = ChisqTail / nPointsTail;
  }


  m_ChisqNDF_Pedestal = ChisqPedestal;
  m_ChisqNDF_Front     = ChisqFront;
  m_ChisqNDF_Tail     = ChisqTail;
  m_ChisqNDF_Rear     = ChisqRear;
  if( m_ChisqNDF_Pedestal >= 0xFFFF ||
      m_ChisqNDF_Front    >= 0xFFFF ||
      m_ChisqNDF_Rear     >= 0xFFFF ||
      m_ChisqNDF_Tail     >= 0xFFFF ){
    return false;
  }else{
    return true;
  }
}

bool     E14WaveFitter::FitTime( TGraph* gr){
  m_splTime = m_FitFunc->GetX( m_FitFunc->GetParameter(0)*0.5 + m_FitFunc->GetParameter(2),
			       m_FitFunc->GetParameter(1)-80,   m_FitFunc->GetParameter(1));
  //std::cout<< m_splTime << std::endl; 
  //std::cout<< "Fit" << std::endl;
  if( m_splTime - 16 < 0 ){ 
    m_slope  = 0;
    m_HHTime = 0;
    return false;
  }
  
  if( m_splTime + 16 > gr->GetX()[gr->GetN()-1]){
    m_slope = 0;
    m_HHTime = 9999;
    return false;
  }
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
	    << m_ChisqNDF_Front    << "\t" 
	    << m_ChisqNDF_Rear     << "\t" 
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

