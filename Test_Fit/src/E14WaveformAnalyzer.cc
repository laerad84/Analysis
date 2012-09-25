#include "E14WaveformAnalyzer.h"

Double_t pol2func( Double_t* x, Double_t* par){
  Double_t p0 = par[0];
  Double_t p1 = par[1];
  Double_t p2 = par[2];
  Double_t x0 = x[0];  
  Double_t Value = p0*( x0 - p1 )*( x0 - p1 ) + p2;
  return Value; 
}

E14WaveformAnalyzer::E14WaveformAnalyzer( int nPoints ){
  m_nPoint = nPoints; 
  m_peakGraph = new TGraph();
  m_peakFunc  = new TF1("peakFunc",pol2func,0,m_nPoint*8,2);
  _Clear();
}
E14WaveformAnalyzer::~E14WaveformAnalyzer(){
  _Clear();
  m_peakFunc->Delete();
  m_peakGraph->Delete();
}
void   E14WaveformAnalyzer::_Clear(){
  if( m_fMaximum ){
    m_peakGraph->GetListOfFunctions()->Delete();
  }

  m_fPedestal = false;
  m_fHeight   = false;  
  m_fMaximum  = false;
  m_fMinimum  = false;
  m_fWidth    = false;
  m_fmHead    = false;
  m_fmTail    = false;
  m_fBoundary = false;
  m_fADC      = false;

  m_MeanHead = 0;
  m_MeanTail = 0;
  m_SigmaHead = 0; 
  m_SigmaTail = 0;
  m_Pedestal = 0; 
  m_PedestalSigma = 0; 
  m_Height = 0;
  m_PeakMaximum = 0;
  m_PeakMinimum = 0; 
  m_TimeMaximum = 0;
  m_TimeMinimum = 0;
  m_BoundaryHead = 0;
  m_BoundaryTail = 0;
  m_Width = 0;
  m_ADC = 0; 

}
Bool_t E14WaveformAnalyzer::_GetMeanHead( Double_t* Waveform ){
  // GetMean 4points of Waveform head // 
  m_MeanHead  = 0.; 
  m_SigmaHead = 0.; 
  for( int ipoint  =0;  ipoint < m_CnHead; ipoint++){
    m_MeanHead  += Waveform[ipoint];
    m_SigmaHead += Waveform[ipoint]*Waveform[ipoint];
  }
  m_MeanHead  = m_MeanHead / m_CnHead;
  m_SigmaHead = TMath::Sqrt(m_SigmaHead/m_CnHead - m_MeanHead*m_MeanHead);
  m_fmHead    = true;
  return      m_fmHead;
}
Bool_t E14WaveformAnalyzer::_GetMeanHead( Double_t* Waveform , Double_t& MeanHead, Double_t& SigmaHead ){
  // GetMean 4points of Waveform head // 
  if( !m_fmHead ){
    _GetMeanHead( Waveform );
  }

  MeanHead    = m_MeanHead;
  SigmaHead   = m_SigmaHead; 
  m_fmHead    = true;
  return      m_fmHead;
}
Bool_t E14WaveformAnalyzer::_GetMeanTail( Double_t* Waveform ){
  m_MeanTail  = 0.;
  m_SigmaTail = 0.;
  for( int ipoint = m_nPoint-m_CnTail; ipoint < m_nPoint; ipoint++){
    m_MeanTail  += Waveform[ ipoint ]; 
    m_SigmaTail += Waveform[ ipoint ]*Waveform[ipoint];
  }
  m_MeanTail  = m_MeanTail/ m_CnTail;
  m_SigmaTail = TMath::Sqrt( m_SigmaTail/m_CnTail - m_MeanTail*m_MeanTail );
  m_fmTail    = true;
  return      m_fmTail; 
}
Bool_t E14WaveformAnalyzer::_GetMeanTail( Double_t* Waveform, Double_t& MeanTail , Double_t& SigmaTail ){
  if( !m_fmTail ){
    _GetMeanTail(Waveform);
  }
  MeanTail    = m_MeanTail;
  SigmaTail   = m_SigmaTail; 
  m_fmTail    = true;
  return      m_fmTail; 
}
Bool_t E14WaveformAnalyzer::_GetMaximum ( Double_t* Waveform ){
  // Taking 4point Mean maximum //
  m_PeakMaximum = 0;
  m_TimeMaximum = 0; 
  Double_t localMaximum;
  Double_t localMaxTime;    
  Int_t    MaxPoint; 
  for( int ipoint = m_CnHead; ipoint < m_nPoint - m_CnTail -4 ; ipoint++ ){
    localMaximum = 0; 
    localMaxTime = 0;
    for( int iSubPoint  =0; iSubPoint <4; iSubPoint++){
      localMaximum += Waveform[ipoint+iSubPoint]; 
      localMaxTime += Waveform[ipoint+iSubPoint]*ipoint;       
    }
    localMaximum = localMaximum / 4;
    localMaxTime = localMaxTime / 4 / localMaximum; 
    if( localMaximum > m_PeakMaximum ){
      MaxPoint      = ipoint;
      m_PeakMaximum = (Waveform[MaxPoint+1] +Waveform[MaxPoint+2])/2;
      m_TimeMaximum = localMaxTime*m_TimeWidth; 
    }
  }

  // Fit with pol2 
  m_peakGraph->Set(0);
  for( int ipoint = MaxPoint-2; ipoint < MaxPoint + 4 + 2 ; ipoint++){ 
    m_peakGraph->SetPoint( m_peakGraph->GetN(), ipoint*8, Waveform[ipoint]);
  }

  m_peakFunc->SetParLimits(1,m_TimeMaximum - 8   , m_TimeMaximum +8  );
  m_peakFunc->SetParLimits(2,m_PeakMaximum * 0.9 , m_PeakMaximum*1.2 );
  m_peakFunc->SetParameter(0,-100);
  m_peakFunc->SetParameter(1,m_TimeMaximum);
  m_peakFunc->SetParameter(2,m_PeakMaximum*1.05);

  m_peakGraph->Fit(m_peakFunc,"Q","",( MaxPoint-2)*8, (MaxPoint+4+2)*8);
  
  std::cout << m_peakFunc->GetParameter(2) << std::endl;
  std::cout << m_peakFunc->GetParameter(1) << std::endl;
  m_PeakMaximum = m_peakFunc->GetParameter(2);
  m_TimeMaximum = m_peakFunc->GetParameter(1);
  return m_fMaximum; 
}

Bool_t E14WaveformAnalyzer::_GetMaximum ( Double_t* Waveform , Double_t& PeakMaximum, Double_t& TimeMaximum){
  // Taking 4point Mean maximum //
  if( !m_fMaximum ){
    _GetMaximum( Waveform );
  }

  PeakMaximum = m_PeakMaximum;
  TimeMaximum = m_TimeMaximum; 
  return m_fMaximum; 
}

Bool_t E14WaveformAnalyzer::_GetMinimum     ( Double_t* Waveform ){
  m_PeakMinimum = 0xFFFF;
  m_TimeMinimum = 0;
  Double_t localMinimum;
  Double_t localMinTime;
  if( !m_fMaximum ){
    _GetMaximum( Waveform );
  }
  for( int ipoint = m_CnHead; ipoint < m_nPoint - m_CnTail - 4; ipoint++){
    localMinimum = 0;
    localMinTime = 0;
    for( int iSubPoint = 0; iSubPoint < 4; iSubPoint++){
      localMinimum += Waveform[ ipoint + iSubPoint];
      localMinTime += Waveform[ ipoint + iSubPoint]*(ipoint + iSubPoint);
    }
    localMinimum = localMinimum / 4.;   
    localMinTime = localMinTime / 4. / localMinimum;

    if( localMinTime >=m_TimeMaximum ){
      break;
    }
    if( localMinimum < m_PeakMinimum ){
      m_PeakMinimum = localMinimum;
      m_TimeMinimum = localMinTime;
    }
  } 
  m_fMinimum = true;
  return m_fMinimum;
}
Bool_t E14WaveformAnalyzer::_GetMinimum     ( Double_t* Waveform , Double_t& PeakMinimum, Double_t& TimeMinimum ){
  if( !m_fMinimum ){
    _GetMinimum( Waveform );
  }
  PeakMinimum = m_PeakMinimum;
  TimeMinimum = m_PeakMinimum;
  m_fMinimum = true;
  return m_fMinimum;
}
Bool_t E14WaveformAnalyzer::_GetPedestal( Double_t* Waveform ){
  if( !m_fMaximum ){
    _GetMaximum( Waveform ); 
  }
  if( !m_fMinimum ){
  _GetMinimum( Waveform );   
  }
  Int_t MinimumPoint = (int)(m_TimeMinimum/8);
  Int_t MaximumPoint = (int)(m_TimeMaximum/8);
  Int_t nPedPoint    = 0; 
  m_Pedestal = 0.;
  m_PedestalSigma = 0.; 
  for( int ipoint = MinimumPoint - 3; ipoint < MinimumPoint + 3; ipoint++){       
    if( ipoint < MaximumPoint - 10 ) { break ; } // All Signal has width ~ 80ns in front regeon
    m_Pedestal      += Waveform[ipoint];
    m_PedestalSigma += Waveform[ipoint]* Waveform[ipoint];
    nPedPoint++;
  }
  if( nPedPoint ==0 ){
    m_Pedestal = 0;
    m_PedestalSigma = 0; 
    m_fPedestal = false;
  }else{
    m_Pedestal /= nPedPoint;
    m_PedestalSigma = TMath::Sqrt( m_Pedestal/nPedPoint - m_Pedestal*m_Pedestal );
    m_fPedestal  =true; 
  }
  return m_fPedestal; 
} 
Bool_t E14WaveformAnalyzer::_GetPedestal( Double_t* Waveform , Double_t& Pedestal , Double_t& PedestalSigma ){
  if( !m_fMaximum ){
    _GetMaximum( Waveform ); 
  }
  if( !m_fMinimum ){
  _GetMinimum( Waveform );   
  }
  Int_t MinimumPoint = (int)(m_TimeMinimum/8);
  Int_t MaximumPoint = (int)(m_TimeMaximum/8);
  Int_t nPedPoint    = 0; 
  m_Pedestal = 0.;
  m_PedestalSigma = 0.; 
  for( int ipoint = MinimumPoint - 3; ipoint < MinimumPoint + 3; ipoint++){       
    if( ipoint < MaximumPoint - 10 ) { break ; } // All Signal has width ~ 80ns in front regeon
    m_Pedestal      += Waveform[ipoint];
    m_PedestalSigma += Waveform[ipoint]* Waveform[ipoint];
    nPedPoint++;
  }
  if( nPedPoint ==0 ){
    m_Pedestal = 0;
    m_PedestalSigma = 0; 
    m_fPedestal = false;
  }else{
    m_Pedestal /= nPedPoint;
    m_PedestalSigma = TMath::Sqrt( m_Pedestal/nPedPoint - m_Pedestal*m_Pedestal );
    m_fPedestal  =true; 
  }
  
  Pedestal  = m_Pedestal;
  PedestalSigma = m_PedestalSigma; 

  return m_fPedestal; 
}
Bool_t E14WaveformAnalyzer::_GetWidth   ( Double_t* Waveform ){  
  if( !m_fPedestal ){
    _GetPedestal( Waveform );
  }
  if( !m_fMaximum ){
    _GetMaximum( Waveform );
  }
  
  m_Height = m_PeakMaximum - m_Pedestal;
  m_BoundaryHead = 0;
  m_BoundaryTail = m_nPoint*8;
  Int_t MaxPoint = (int)(m_TimeMaximum/8);
  Int_t MinPoint = (int)(m_TimeMinimum/8);
  for( int ipoint = MaxPoint; ipoint > MinPoint; ipoint-- ){
    if( Waveform[ ipoint - 1 ] - m_Pedestal  < m_Height / 2 &&
	Waveform[ ipoint ]     - m_Pedestal  > m_Height / 2 ){
      m_BoundaryHead = 8*(ipoint - 1) + 8*( m_Height/2 - (Waveform[ipoint -1] - m_Pedestal ))/(Waveform[ipoint] -Waveform[ipoint -1]);
      break; 
    }else if(( Waveform[ipoint] - m_Pedestal ) == m_Height/2 ){
      m_BoundaryHead = ipoint*8;
      break;
    }
  }
  for( int ipoint = MaxPoint; ipoint < m_nPoint - m_CnTail; ipoint++){
    if( Waveform[ ipoint +1 ] - m_Pedestal < m_Height / 2 &&
	Waveform[ ipoint ]    - m_Pedestal > m_Height / 2 ){
      m_BoundaryTail = 8*(ipoint + 1) - 8*( m_Height/2 - (Waveform[ipoint + 1] -m_Pedestal ))/(Waveform[ipoint] -Waveform[ipoint +1]);
      break;
    }else if( Waveform[ipoint] - m_Pedestal  == m_Height/2 ){
      m_BoundaryTail = ipoint*8;
      break;
    }
  }
  
  if( m_BoundaryHead == 0 ||
      m_BoundaryTail == m_nPoint*8 ){
    m_fBoundary = false;
  }else{
    m_Width = m_BoundaryTail - m_BoundaryHead;
    m_fBoundary = true;
  }
  return m_fBoundary; 				     
}
Bool_t E14WaveformAnalyzer::_GetWidth   ( Double_t* Waveform ,Double_t& Width){  
  if( m_fBoundary ){
    _GetWidth( Waveform );
  }
  if( m_fBoundary ){
    Width == m_Width;
  }else{
    Width == m_nPoint*8;
  }
  return m_fBoundary; 				     

}

Bool_t E14WaveformAnalyzer::_GetADC( Double_t* Waveform , Int_t StartTime, Int_t EndTime ){
  if( !m_fPedestal ){
    _GetPedestal( Waveform );
  }
  m_ADC = 0;
  for( int ipoint =  StartTime ; ipoint < EndTime; ipoint ++){
    m_ADC += Waveform[ipoint] - m_Pedestal; 
  }
  m_fADC = true; 
  return m_fADC;
}

Bool_t E14WaveformAnalyzer::_GetADC( Double_t* Waveform , Double_t& ADC, Int_t StartTime , Int_t EndTime ){
  if( !m_fADC ){
    _GetADC( Waveform, StartTime, EndTime );
  }
  ADC = m_ADC;
  return m_fADC;
}

void E14WaveformAnalyzer::_Draw( Double_t* Waveform ){
  TGraph* grWave = new TGraph();
  for( int ipoint  =0 ; ipoint < m_nPoint; ipoint++){
    grWave->SetPoint( ipoint, ipoint*8, Waveform[ipoint]);
  }
  if( !m_fmHead ) _GetMeanHead( Waveform );
  if( !m_fmTail ) _GetMeanTail( Waveform );
  if( !m_fMaximum ) _GetMaximum( Waveform );
  if( !m_fMinimum ) _GetMinimum( Waveform );
  if( !m_fPedestal ) _GetPedestal( Waveform );
  if( !m_fWidth ) _GetWidth( Waveform );
  if( !m_fADC ) _GetADC( Waveform );
  grWave->SetMarkerStyle(6);
  grWave->Draw("AP");
  TLine* LPedestal = new TLine();
  TLine* LHeight   = new TLine();
  TLine* LBoundaryHead = new TLine();
  TLine* LBoundaryTail = new TLine();
  TLine* LHalfHeight   = new TLine();
  LPedestal    ->SetLineColor(1);
  LHeight      ->SetLineColor(2);
  LBoundaryHead->SetLineColor(3);
  LBoundaryTail->SetLineColor(4);
  LHalfHeight  ->SetLineColor(5);
  LPedestal    ->DrawLine( 0, m_Pedestal   , 8*m_nPoint, m_Pedestal);
  LHeight      ->DrawLine( 0, m_PeakMaximum, 8*m_nPoint, m_PeakMaximum);
  LHalfHeight  ->DrawLine( 0, (m_PeakMaximum + m_PeakMinimum)/2, 8*m_nPoint, ( m_PeakMaximum + m_PeakMinimum )/2);
  LBoundaryHead->DrawLine( m_BoundaryHead, m_PeakMinimum, m_BoundaryHead, m_PeakMaximum);
  LBoundaryTail->DrawLine( m_BoundaryTail, m_PeakMinimum, m_BoundaryTail, m_PeakMaximum);

  std::cout << "//////////////////////////////////////////////////////////\n"
	    << "***Parameter***\n"
	    << "MeanHead     :" << m_MeanHead      << "\n"
	    << "MeanTail     :" << m_MeanTail      << "\n"
	    << "Maximum      :" << m_PeakMaximum   << "\n"
	    << "Minimum      :" << m_PeakMinimum   << "\n"
	    << "Pedestal     :" << m_Pedestal      << "\n"
	    << "PedestalSigma:" << m_PedestalSigma << "\n" 
	    << "Height       :" << m_Height        << "\n"
	    << "BoundaryHead :" << m_BoundaryHead  << "\n"
	    << "BoundaryTail :" << m_BoundaryTail  << "\n"
	    << std::endl; 
}
