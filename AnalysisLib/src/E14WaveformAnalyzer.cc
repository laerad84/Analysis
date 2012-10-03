#include "E14WaveformAnalyzer.h"

Double_t pol2func( Double_t* x, Double_t* par){
  Double_t p0 = par[0];
  Double_t p1 = par[1];
  Double_t p2 = par[2];
  Double_t x0 = x[0];  
  Double_t Value = -1*p0*( x0 - p1 )*( x0 - p1 ) + p2;
  return Value; 
}

E14WaveformAnalyzer::E14WaveformAnalyzer( int nPoints ){
  m_nPoint = nPoints; 
  m_peakGraph = new TGraph();
  m_peakFunc  = new TF1("peakFunc",pol2func,0,m_nPoint*8,3);
  _Clear();
}
E14WaveformAnalyzer::~E14WaveformAnalyzer(){
  _Clear();
  m_peakFunc->Delete();
  m_peakGraph->Delete();
}
void   E14WaveformAnalyzer::_Clear(){
  if( m_peakGraph != NULL ){
    if( m_peakGraph->GetListOfFunctions()->GetEntries() != 0){
      m_peakGraph->GetListOfFunctions()->Delete();
    }
  }

  m_fPedestal  = false;
  m_fHeight    = false;  
  m_fMaximum   = false;
  m_fMinimum   = false;
  m_fWidth     = false;
  m_fmHead     = false;
  m_fmTail     = false;
  m_fBoundary  = false;
  m_fADC       = false;
  m_fSlope     = false;
  m_fOverflow  = false;
  m_fUnderflow = false;
  m_fWaveformAnalyzed = false;
  m_WaveformState = 0;

  m_MeanHead      = 0;
  m_MeanTail      = 0;
  m_SigmaHead     = 0; 
  m_SigmaTail     = 0;
  m_Pedestal      = 0; 
  m_PedestalSigma = 0; 
  m_MinimumSigma  = 0;
  m_Height        = 0;
  m_PeakMaximum   = 0;
  m_PeakMinimum   = 0; 
  m_TimeMaximum   = 0;
  m_TimeMinimum   = 0;
  m_BoundaryHead  = 0;
  m_BoundaryTail  = 0;
  m_Width         = 0;
  m_ADC           = 0; 
  m_SlopeStart    = 0;
  m_SlopeDelta    = -1*0xFFFF;
  
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
  m_PeakMaximum = -1*0xFFFF;
  m_TimeMaximum = 0; 
  Double_t localMaximum = 0;
  Double_t localMaxTime = 0;
  Double_t localPointMaximum = 0;

  m_AllPointMaximum  = 0;
  m_AllPointMinimum  = 0; 
  m_AllPointMax      = 0;
  m_AllPointMin      = 0xFFFF;
  m_PeakPointMaximum = 0;

  // Check Underflow && Overflow 
  for( int ipoint = 0; ipoint < m_nPoint;  ipoint++){
    if( m_AllPointMaximum < Waveform[ ipoint ] ){
      m_AllPointMaximum  = Waveform[ ipoint ];
      m_AllPointMax = ipoint;
    }
    if( m_AllPointMinimum > Waveform[ ipoint ] ){
      m_AllPointMinimum = Waveform[ ipoint ];
      m_AllPointMin = ipoint;
    }
  } 
  if( m_AllPointMaximum > 16000 ){ m_fOverflow = true; }
  if( m_AllPointMinimum < 1     ){ m_fUnderflow = true; }

  for( int ipoint = m_CnHead+4; ipoint < m_nPoint - m_CnTail - 4 ; ipoint++ ){
    localPointMaximum  = 0;
    localMaximum       = 0; 

    for( int iSubPoint = 0; iSubPoint <4; iSubPoint++){
      localMaximum += Waveform[ipoint+iSubPoint]; 
      localMaxTime += Waveform[ipoint+iSubPoint]*(ipoint+iSubPoint);       
      if( Waveform[ ipoint +iSubPoint ] > localPointMaximum ){
	localPointMaximum = Waveform[ ipoint + iSubPoint ];
      }
    }

    localMaxTime = ((Double_t)ipoint+1.5)*m_TimeWidth;
    localMaximum = localMaximum / 4;
    if( localMaximum > m_PeakMaximum ){
      m_PeakMaximum = localMaximum;
      m_TimeMaximum = localMaxTime;
      m_PeakPointMaximum = localPointMaximum;
    }
    
    /*
    if( localMaximum > m_PeakMaximum ){
      MaxPoint      = ipoint;
      m_PeakMaximum = (Waveform[MaxPoint+1] + Waveform[MaxPoint+2])/2;
      m_TimeMaximum = localMaxTime*m_TimeWidth; 
      m_PeakPointMaximum = localPointMaximum;
    }
    */
  }

  //m_PeakMaximum = m_PeakPointMaximum; 

  /* 
  // Fit with pol2 -> Abandon
  m_peakGraph->Set(0);
  // Check Local Maximum
  // if Maximum point exists on edge Fitting became craze... //  
  for( int ipoint = MaxPoint-2; ipoint < MaxPoint + 4 + 3 ; ipoint++){ 
    m_peakGraph->SetPoint( m_peakGraph->GetN(), ipoint*8, Waveform[ipoint]);
  }
  m_peakFunc->SetParameter(0, (2*m_PeakMaximum - (Waveform[MaxPoint-4]+Waveform[MaxPoint+4])/1024));
  m_peakFunc->SetParameter(1, m_TimeMaximum);
  m_peakFunc->SetParameter(2, m_PeakMaximum);
  m_peakFunc->SetParLimits(0, 0.001, 4);
  m_peakFunc->SetParLimits(1, m_TimeMaximum - 16, m_TimeMaximum + 16);
  m_peakFunc->SetParLimits(2, m_PeakMaximum*0.8 , m_PeakMaximum*1.2 );  
  m_peakGraph->Fit(m_peakFunc,"Q","",( MaxPoint-2)*8, (MaxPoint+4+2)*8);  
  m_PeakMaximum = m_peakFunc->GetParameter(2);
  m_TimeMaximum = m_peakFunc->GetParameter(1);
  */

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
  Double_t localMinSigma;
  Double_t localMinimum;
  Double_t localMinTime;
  if( !m_fMaximum ){
    _GetMaximum( Waveform );
  }
  for( int ipoint = 0; ipoint*8 < m_TimeMaximum ; ipoint++){
    localMinimum = 0;
    localMinTime = 0;
    for( int iSubPoint = 0; iSubPoint < 4; iSubPoint++){
      localMinimum  += Waveform[ ipoint + iSubPoint ];
      localMinSigma += Waveform[ ipoint + iSubPoint ]*Waveform[ ipoint + iSubPoint ];
      //localMinTime += Waveform[ ipoint + iSubPoint]*(ipoint + iSubPoint);
    }
    localMinimum  = localMinimum / 4.;   
    localMinSigma = TMath::Sqrt(localMinSigma/4 - localMinimum*localMinimum);
    //localMinTime = localMinTime / 4. / localMinimum;
    localMinTime = ((Double_t)ipoint +1.5)*8;
    if( localMinTime >=m_TimeMaximum ){
      break;
    }
    if( localMinimum < m_PeakMinimum ){
      m_PeakMinimum  = localMinimum;
      m_TimeMinimum  = localMinTime;
      m_MinimumSigma = localMinSigma;
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

  Int_t nPedPoint    = 0; 
  m_Pedestal         = 0.;
  m_PedestalSigma    = 0.; 
  Double_t startTime = 0;
  Double_t endTime   = 0;

  // Take 6 Point of Pedestal  in front of Peak // 
  endTime   = m_TimeMaximum - 88;
  startTime = m_TimeMaximum - 136;
  if( startTime < 0 ){
    startTime = 0;
  }
  Int_t StartPoint = (Int_t)(m_TimeMinimum/8);
  for( int ipoint = StartPoint; ipoint*8 < endTime; ipoint++){       
    m_Pedestal      += Waveform[ipoint];
    m_PedestalSigma += Waveform[ipoint]* Waveform[ipoint];
    nPedPoint++;
  }

  if( nPedPoint                     < 4                 ||
      m_Pedestal                    < m_PeakMinimum     ||
      m_Pedestal    - m_PeakMinimum > m_MinimumSigma*3  ||
      m_TimeMaximum - m_TimeMinimum < 12*m_TimeWidth    ){
    m_Pedestal      = m_PeakMinimum;
    m_PedestalSigma = m_MinimumSigma;
    m_fPedestal     = true;
  }else{
    m_Pedestal /= nPedPoint;
    m_PedestalSigma = TMath::Sqrt( m_PedestalSigma/nPedPoint - m_Pedestal*m_Pedestal );
    m_fPedestal  =true; 
  }  
  return m_fPedestal; 
}
Bool_t E14WaveformAnalyzer::_GetPedestal( Double_t* Waveform , Double_t& Pedestal , Double_t& PedestalSigma ){
  if( !m_fPedestal ){
    _GetPedestal( Waveform );
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

  if( m_Height < 0 ){
    m_BoundaryTail = m_nPoint*8;
    m_BoundaryHead = 0;
    m_fBoundary = false;
    return m_fBoundary;
  }

  Int_t MaxPoint = (int)(m_TimeMaximum/8);
  Int_t MinPoint = (int)(m_TimeMinimum/8);  
  Int_t MaxHeadPoint = (int)(m_TimeMaximum/8 - 1.5);
  Int_t MinHeadPoint = (int)(m_TimeMinimum/8 - 1.5);

  for( int ipoint = MaxHeadPoint+3; ipoint > MinPoint; ipoint-- ){
    if( Waveform[ ipoint - 1 ] - m_Pedestal  < m_Height / 2 &&
	Waveform[ ipoint ]     - m_Pedestal  > m_Height / 2 ){
      m_BoundaryHead = 8*(ipoint - 1) + 8*( m_Height/2 - (Waveform[ipoint -1] - m_Pedestal ))/(Waveform[ipoint] -Waveform[ipoint -1]);
      break; 
    }else if(( Waveform[ipoint] - m_Pedestal ) == m_Height/2 ){
      m_BoundaryHead = ipoint*8;
      break;
    }
  }

  for( int ipoint = MaxHeadPoint; ipoint < m_nPoint - m_CnTail; ipoint++){
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

Bool_t E14WaveformAnalyzer::_GetWidth   ( Double_t* Waveform ,Double_t& BoundaryHead,Double_t& BoundaryTail){  
  if( !m_fBoundary ){
    _GetWidth( Waveform );
  }
  BoundaryHead = m_BoundaryHead;
  BoundaryTail = m_BoundaryTail;    
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



Bool_t   E14WaveformAnalyzer::_GetSumSlope( Double_t* Waveform ){
  if( m_nPoint > 1024){
    std::cerr << "FX TOO LONG" << std::endl;    
    m_fSlope = false; 
    return m_fSlope;
  }
  if( !m_fPedestal ){
    _GetPedestal( Waveform );
  }

  Double_t Sum[1024] = {0}; // Size is Must bigger than m_nPoints // 
  for( int i = 0; i < m_nPoint; i++){
    if( i == 0 ){
      Sum[i] = Waveform[i] - m_Pedestal;
    }else{
      Sum[i] = Sum[i-1] + Waveform[i] - m_Pedestal;
    }
  }

  Double_t DeltaMaximum  = -1*0xFFFF;
  Int_t    StartPoint    = 0;   

  // 4 front , 8 middle, 4 tail compare // 
  for( int i = 0; i < m_nPoint - 15; i++){
    Double_t HeadSumMean   = 0;
    Double_t TailSumMean   = 0;
    Double_t MiddleSumMean = 0; 
    Double_t RegionDelta   = 0;
    for( int isubpoint  = 0; isubpoint < 4; isubpoint++){
      HeadSumMean   += ( Sum[ i + isubpoint ] / 4.) ;      
    }
    for( int isubpoint =  4; isubpoint < 12; isubpoint++){
      MiddleSumMean += ( Sum[ i + isubpoint ] / 12 );
    }
    for( int isubpoint =  12; isubpoint < 16; isubpoint++){
      TailSumMean   += ( Sum[ i + isubpoint ] / 4.);
    }
    RegionDelta  = TailSumMean - HeadSumMean;
    if( RegionDelta  >  DeltaMaximum ){
      DeltaMaximum = RegionDelta;
      StartPoint   = i;
    }
  }

  m_SlopeStart = StartPoint;
  m_SlopeDelta = DeltaMaximum;
  m_fSlope = true; 
  return m_fSlope; 
}
Bool_t   E14WaveformAnalyzer::_GetSumSlope( Double_t* Waveform , Int_t& StartPoint, Double_t& DeltaMaximum ){
  if( !m_fSlope ){
    _GetSumSlope( Waveform );    
  }
  StartPoint = m_SlopeStart;
  DeltaMaximum = m_SlopeDelta;
  return m_fSlope;
}

Double_t E14WaveformAnalyzer::GetSlopeDelta( Double_t* Waveform , Int_t StartPoint){
  if(  (StartPoint + 16) - m_nPoint <= 0 || m_nPoint > 1024 ){
    return -1*0xFFFF;
  } 
  if( !m_fPedestal ){
    _GetPedestal( Waveform );
  }
  
  Double_t Sum[1024];
  for( int i = 0; i< m_nPoint; i++){
    if( i ==0 ){
      Sum[i] = Waveform[i] - m_Pedestal ;
    }else{
      Sum[i] = Sum[i-1] + Waveform[i] - m_Pedestal; 
    }
  }
  
  Double_t HeadSumMean = 0;
  Double_t TailSumMean = 0;
  for( int iSubPoint  = 0; iSubPoint < 4; iSubPoint++ ){
    HeadSumMean += ( Sum[ StartPoint + iSubPoint ]      /4.);    
    TailSumMean += ( Sum[ StartPoint + iSubPoint + 12 ] /4.);
  }
  
  return TailSumMean - HeadSumMean ;
}
Int_t    E14WaveformAnalyzer::AnalyzeWaveform( Double_t *Waveform ){
  _Clear();
  m_fWaveformAnalyzed = false;
  m_WaveformState = 0;
  _GetMeanHead( Waveform );
  _GetMeanTail( Waveform );
  _GetMaximum ( Waveform );
  _GetMinimum ( Waveform );
  _GetPedestal( Waveform );
  _GetWidth   ( Waveform );
  _GetSumSlope( Waveform );  
  m_fWaveformAnalyzed = true;
  if( m_PeakMaximum > 16000)                { m_WaveformState |= 1;  }
  if( m_PeakMinimum < 1    )                { m_WaveformState |= 2;  }
  if( m_MeanHead - m_Pedestal > 10 )          { m_WaveformState |= 4;  }
  if( m_MeanTail - m_Pedestal > 10 )          { m_WaveformState |= 8;  }
  if( m_TimeMaximum > (m_nPoint-15)*m_TimeWidth){ m_WaveformState |= 16; }
  return m_WaveformState; 
}
Bool_t   E14WaveformAnalyzer::GetParameters( Double_t& Pedestal, Double_t& Height, Double_t& PeakTime) const {
  if( !m_fWaveformAnalyzed ) { 
    std::cerr << "Please Analyze Waveform beform GetParameter" << std::endl;
    Pedestal = 0;
    Height   = 0;
    PeakTime = 0;
    return false;
  }
  Pedestal = m_Pedestal; 
  Height   = m_PeakMaximum - m_Pedestal;
  PeakTime = m_TimeMaximum;
  return true;
}
void     E14WaveformAnalyzer::Draw( Double_t* Waveform ){
  /*
  if( !m_fmHead    ) _GetMeanHead( Waveform );
  if( !m_fmTail    ) _GetMeanTail( Waveform );
  if( !m_fMaximum  ) _GetMaximum( Waveform );
  if( !m_fMinimum  ) _GetMinimum( Waveform );
  if( !m_fPedestal ) _GetPedestal( Waveform );
  if( !m_fWidth    ) _GetWidth( Waveform );
  if( !m_fADC      ) _GetADC( Waveform );
  if( !m_fSlope    ) _GetSumSlope( Waveform );
  */

  if( !m_fWaveformAnalyzed ){ AnalyzeWaveform( Waveform );}
  TGraph* grWave    = new TGraph();
  for( int ipoint  =0 ; ipoint < m_nPoint; ipoint++){
    grWave->SetPoint( ipoint, ipoint*8, Waveform[ipoint]);
  }

  grWave->SetMarkerStyle(6);
  grWave->Draw("AP");
  TLine* LPedestal = new TLine();
  TLine* LHeight   = new TLine();
  TLine* LBoundaryHead = new TLine();
  TLine* LBoundaryTail = new TLine();
  TLine* LHalfHeight   = new TLine();
  TLine* LMinimum      = new TLine();
  TLine* LMaximum      = new TLine();

  TLine* LPeakTime     = new TLine();
  TLine* LSumTimeHead  = new TLine();
  TLine* LSumTimeTail  = new TLine();


  LPedestal    ->SetLineColor(1);
  LHeight      ->SetLineColor(2);
  LBoundaryHead->SetLineColor(3);
  LBoundaryTail->SetLineColor(4);
  LHalfHeight  ->SetLineColor(5);
  LMinimum     ->SetLineColor(6);
  LMaximum     ->SetLineColor(7);
  LSumTimeHead->SetLineWidth(2);
  LSumTimeTail->SetLineWidth(2);
  LSumTimeHead->SetLineColor(2);
  LSumTimeTail->SetLineColor(2);
  LSumTimeHead->DrawLine( (m_SlopeStart + 4 )*8, m_PeakMinimum , ( m_SlopeStart + 4 )*8, m_PeakMaximum );
  LSumTimeTail->DrawLine( (m_SlopeStart + 12)*8, m_PeakMinimum , ( m_SlopeStart + 12)*8, m_PeakMaximum ); 
  
  LPeakTime    ->DrawLine( m_TimeMaximum, m_PeakMinimum, m_TimeMaximum, m_PeakMaximum);
  LPedestal    ->DrawLine( 0, m_Pedestal   , 8*m_nPoint, m_Pedestal);
  LHalfHeight  ->DrawLine( 0, (m_PeakMaximum + m_Pedestal)/2, 8*m_nPoint, ( m_PeakMaximum + m_Pedestal)/2);
  LBoundaryHead->DrawLine( m_BoundaryHead, m_PeakMinimum, m_BoundaryHead, m_PeakMaximum);
  LBoundaryTail->DrawLine( m_BoundaryTail, m_PeakMinimum, m_BoundaryTail, m_PeakMaximum);
  LMinimum     ->DrawLine( 0, m_PeakMinimum, 8*m_nPoint, m_PeakMinimum);
  LMaximum     ->DrawLine( 0, m_PeakMaximum, 8*m_nPoint, m_PeakMaximum);

  std::cout << "//////////////////////////////////////////////////////////\n"
	    << "***Parameter***\n"
	    << "MeanHead     :" << m_MeanHead      << "\n"
	    << "MeanTail     :" << m_MeanTail      << "\n"
	    << "Maximum      :" << m_PeakMaximum   << "\n"
	    << "Minimum      :" << m_PeakMinimum   << "\n"
	    << "PeakTime     :" << m_TimeMaximum   << "\n"
	    << "PedTime      :" << m_TimeMinimum   << "\n"
	    << "Pedestal     :" << m_Pedestal      << "\n"
	    << "PedestalSigma:" << m_PedestalSigma << "\n" 
	    << "Height       :" << m_Height        << "\n"
	    << "BoundaryHead :" << m_BoundaryHead  << "\n"
	    << "BoundaryTail :" << m_BoundaryTail  << "\n"
	    << "Slope Start  :" << m_SlopeStart    << "\n"
	    << "Slope Delta  :" << m_SlopeDelta    << "\n"    
	    << std::endl; 
}
