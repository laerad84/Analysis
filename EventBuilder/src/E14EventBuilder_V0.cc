#include "E14EventBuilder_V0.h"

E14EventBuilder_V0::E14EventBuilder_V0(TTree* trout, Int_t RunNumber){
  m_trOut     = trout;
  m_RunNumber = RunNumber;
  m_TimeClusterHist = new TH1D("TimeClusterHist","TimeClusterHist",32, 0, 64*8);
  Init();
}
E14EventBuilder_V0::~E14EventBuilder_V0(){
  ;
}
bool E14EventBuilder_V0::Init(){
  InitEnvironment();
  InitIOFile();
  InitTemplate();
  InitTimeOffset();
  InitTrigger();
  for( int i = 0; i< 20 ;i++){
    COSMIC_THRESHOLD[i] = 1000;
  }
  return true;
}
bool E14EventBuilder_V0::InitEnvironment(){
  ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  wavFitter   = new WaveformFitter( 48, kFALSE );
  Fitter      = new E14WaveFitter();
  wavAnalyzer = new E14WaveformAnalyzer();
  Converter   = new EnergyConverter();
  idHandler   = new IDHandler();
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",
					  ANALIBDIR.c_str()));
  return true;     
}
bool E14EventBuilder_V0::InitIOFile(){
  // Init IO File // 
  int nFilesOpened = 0;
  for( int icrate = 0; icrate < nCrateFeb; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, m_RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
    if( !tf[icrate]->IsOpen() ){
      std::cout << "//////////////////////////////////////////////////\n";
      std::cout << "//////////////////////////////////////////////////\n";
      std::cout << "Warning\n";
      std::cout << "File is not exist : "<< tf[icrate]->GetName() << "\n";
      std::cout << "//////////////////////////////////////////////////\n";
      std::cout << "//////////////////////////////////////////////////\n";
    }else{
      nFilesOpened++;
    }      
  }
  std::cout << Form("%s/Sum%d.root",SUMFILEDIR.c_str(),m_RunNumber) << std::endl;  
  wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),m_RunNumber),
					    m_trOut);
  std::cout<< "Setting Map" << std::endl;
  {
    wConv->AddModule("Csi");
    wConv->AddModule("CC03");
    wConv->AddModule("OEV");
    wConv->AddModule("CV");
    wConv->AddModule("Cosmic");
    wConv->AddModule("Laser");
    wConv->AddModule("Etc");
    CsiModuleID    = 0;    
    CC03ModuleID   = 1;
    OEVModuleID    = 2;
    CVModuleID     = 3;
    CosmicModuleID = 4;
    LaserModuleID  = 5;
    wConv->Set();
    wConv->SetMap();
    wConv->Branch();
    int nentries = conv[0]->GetEntries();  
    for( int icrate = 1; icrate < nCrateFeb; icrate++){
      if( nentries != conv[icrate]->GetEntries() ){
	std::cout << "Entries is Different" << std::endl;
      }
    }
  }
  std::cout<< "Number of Module " << wConv->GetNmodule() << std::endl;
  m_Entries = conv[0]->GetEntries();
  grCsI = new TGraph();
  if( nFilesOpened == nCrateFeb && m_trOut != NULL){
    return true;
  }else{
    return false;
  }
}
bool E14EventBuilder_V0::InitTemplate(){
  // Set Template 
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  TFile* tfTemplate = new TFile(Form("%s/Data/Template/TEMPLATE_COMPLETE_GRAPH_3pi0RunList_200_400.root",ANALIBDIR.c_str()));
  for( int ispline = 0; ispline< 2716; ispline++){
    tempGr[ispline]  = NULL;
    tempSpl[ispline] = NULL;
    //tempGr[ispline]  = (TGraphErrors*)tfTemplate->Get(Form("Waveform_Height_%d_0",ispline));
    tempGr[ispline]  = (TGraph*)tfTemplate->Get(Form("Template_%d",ispline));
    if( tempGr[ispline]->GetN()!= 0){
      tempSpl[ispline] = new TSpline3(Form("splTemplate_%d",ispline),tempGr[ispline]);
    }else{
      std::cout<< "Non Exist channel:" << ispline << std::endl;
    }
  }
  tfTemplate->Close();
  return true;
}
bool E14EventBuilder_V0::InitTimeOffset(){
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  std::string    TimeOffsetFile = Form("%s/Data/TimeOffset/testNewWORKCompileOffset.txt",ANALIBDIR.c_str());
  std::ifstream  ifsTimeOffset(TimeOffsetFile.c_str());
  for( int i = 0; i < 2716; i++ ){
    TimeOffset[i] = 0.;
  }
  
  Int_t    tempID;
  Double_t tempOffset;
  Double_t tempChisq;
  while(  ifsTimeOffset >> tempID >> tempOffset >> tempChisq ){
    TimeOffset[tempID] = tempOffset;
  }
  return true;
}
bool E14EventBuilder_V0::InitTrigger(){
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  grTrigger = new TGraph();
  // Init
  for( int i  = 0 ; i< 3 ; i++ ){
    LaserCFC[i] = 9999;
    for( int j = 0 ; j < 10 ; j++ ){
      CVCFC[j][i] = 9999;
    }
    for( int j = 0; j < 20 ; j++ ){
      CosmicCFC[j][i] = 9999;
    }
  }
  if((wConv->ModMap[CosmicModuleID]).nMod != 20)
    { 
      std::cout<< "Cosmic nModule is not equal" << std::endl;
    }
  if((wConv->ModMap[CVModuleID]).nMod     != 10)
    { 
      std::cout<< "CV nModule is not equal"     << std::endl;
    }
  // Set 
  for( int subID = 0; subID < 3; subID++ ){    
    LaserCFC[subID] = (wConv->ModMap[LaserModuleID]).Map[0][subID];
  }
  for( int icv = 0; icv < (wConv->ModMap[CVModuleID]).nMod; icv++){
    for( int subID = 0; subID < 3; subID++){
      CVCFC[icv][subID] = (wConv->ModMap[LaserModuleID]).Map[0][subID];
    }
  }
  for( int icosmic = 0; icosmic < (wConv->ModMap[CosmicModuleID]).nMod; icosmic++){
    for( int subID = 0; subID < 3; subID++){
      CosmicCFC[icosmic][subID] = (wConv->ModMap[LaserModuleID]).Map[0][subID];
    }
  }
  return true;
}
int  E14EventBuilder_V0::TriggerDicision(){
#ifdef DEBUG
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
#endif
  TotalTriggerFlag  = 0;
  for( int iMod = 0; iMod < wConv->GetNmodule(); iMod++ ){            
    //std::cout << iMod << std::endl;
    if( iMod == wConv->GetModuleID("Csi") ){ continue; }
    int nSubModule = wConv->GetNsubmodule( iMod );
    if( nSubModule <= 0 ){ continue ;}      
    for( int iSubMod = 0; iSubMod < nSubModule; iSubMod++){	
      grTrigger->Set(0);
      if( wConv->SetGraph( iMod, iSubMod ,conv , grTrigger ) == 0 ){ continue; }	             
      bool fit = wavFitter->Fit( grTrigger ); 
      int chIndex  = (wConv->mod[iMod])->m_nDigi;	 
      if( fit ){
	TF1*   fitFunc    = wavFitter->GetFunction();	    
	double halfHeight = fitFunc->GetParameter(0)/2 + fitFunc->GetParameter(4);
	double halfTiming = fitFunc->GetX( halfHeight,
					   fitFunc->GetParameter(1)-48, fitFunc->GetParameter(1));
	wConv->mod[iMod]->m_FitHeight[chIndex]= 1;
	wConv->mod[iMod]->m_ID[chIndex]       = iSubMod;
	wConv->mod[iMod]->m_Pedestal[chIndex] = fitFunc->GetParameter(4);
	wConv->mod[iMod]->m_Signal[chIndex]   = fitFunc->GetParameter(0);
	wConv->mod[iMod]->m_Time[chIndex]     = fitFunc->GetParameter(1);
	wConv->mod[iMod]->m_HHTime[chIndex]   = halfTiming;
	wConv->mod[iMod]->m_ParA[chIndex]     = fitFunc->GetParameter(3);
	wConv->mod[iMod]->m_ParB[chIndex]     = fitFunc->GetParameter(2);
	wConv->mod[iMod]->m_nDigi++;	      	    
	
	if( halfTiming - 12 > 0 && halfTiming +12 < grTrigger->GetX()[grTrigger->GetN()-1] ){
	  TF1* linearFunction = new TF1("func","pol1",halfTiming - 12, halfTiming + 12);
	  grTrigger->Fit( linearFunction, "Q", "", halfTiming -12, halfTiming +12 );
	  double halfFitTiming = linearFunction->GetX( halfHeight, halfTiming -12, halfTiming +12);
	  wConv->mod[iMod]->m_FitTime[chIndex]= halfFitTiming;
	  delete linearFunction;
	}else{
	  wConv->mod[iMod]->m_FitTime[chIndex] = 0;
	}
	wavFitter->Clear();
      }
      grTrigger->GetListOfFunctions()->Delete();
    }	
    
    if( iMod == LaserModuleID ){
      if( wConv->mod[iMod]->m_nDigi > 0){
	if( wConv->mod[iMod]->m_Signal[0] > 100 ){
	  wConv->m_LaserTrig = 1;
	  wConv->m_TrigFlag |= 1;
	  TotalTriggerFlag  |= 1;
	}
      }
    }else if( iMod == CosmicModuleID ){
      //int nSubMod = (wConv->ModMap[iMod]).nMod;
      int nSubMod = wConv->mod[iMod]->m_nDigi;
      for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++ ){	    
	CosmicSignal[CosmicArr[wConv->mod[iMod]->m_ID[iSubMod]]] = wConv->mod[iMod]->m_Signal[iSubMod];
	CosmicTime[CosmicArr[wConv->mod[iMod]->m_ID[iSubMod]]]   = wConv->mod[iMod]->m_Time[iSubMod];
      }
      //std::cout<< __LINE__ << std::endl;
      for( int iCosmic = 0; iCosmic < 5; iCosmic++){
	if( CosmicSignal[ iCosmic    ] > COSMIC_THRESHOLD[ iCosmic     ] ||
	    CosmicSignal[ iCosmic+10 ] > COSMIC_THRESHOLD[ iCosmic +10 ]){
	  wConv->m_CosmicTrigFlagUp |= 1 << iCosmic;
	}
	if( CosmicSignal[ iCosmic+5  ] > COSMIC_THRESHOLD[ iCosmic +5  ] ||
	    CosmicSignal[ iCosmic+15 ] > COSMIC_THRESHOLD[ iCosmic +15 ]){
	  wConv->m_CosmicTrigFlagDn |= 1 << iCosmic;
	  }	  
      }
      //std::cout<< __LINE__ << std::endl;
      if( wConv->m_CosmicTrigFlagUp && 
	  wConv->m_CosmicTrigFlagDn ){
	wConv->m_CosmicTrig = 1; 
	wConv->m_TrigFlag  |= 2;
	TotalTriggerFlag   |= 2;
      }
    }else if( iMod == CVModuleID ){
      //int nSubMod = (wConv->ModMap[iMod]).nMod;
      int nSubMod = wConv->mod[iMod]->m_nDigi;
      for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++ ){	    
	CVSignal[wConv->mod[iMod]->m_ID[iSubMod]] = wConv->mod[iMod]->m_Signal[iSubMod];
	CVTime[wConv->mod[iMod]->m_ID[iSubMod]]   = wConv->mod[iMod]->m_Time[iSubMod];
      }
      for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++){
	  if( wConv->mod[iMod]->m_Signal[iSubMod] > 500 ){
	    wConv->m_CVTrig    = 1;
	    wConv->m_TrigFlag |= 4; 
	    TotalTriggerFlag  |= 4;
	  }
      }
    }else{
      continue;
    } 
  }  
  return TotalTriggerFlag;
}
int  E14EventBuilder_V0::AnalyzeCsIData(){
#ifdef DEBUG
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
#endif
  m_TimeClusterHist->Reset();
  int nSubCsiModule = wConv->GetNsubmodule( CsiModuleID );
  for( int iSubMod = 0; iSubMod < nSubCsiModule; iSubMod++){	
    int iCrate = 9999;
    int iSlot  = 9999;
    int iCh    = 9999;
    iCrate = (wConv->ModMap[CsiModuleID]).Map[iSubMod][0];
    iSlot  = (wConv->ModMap[CsiModuleID]).Map[iSubMod][1];
    iCh    = (wConv->ModMap[CsiModuleID]).Map[iSubMod][2];
    
    // Ignore unmapped channel // 
    if( iCrate == 9999 || iSlot == 9999 || iCh == 9999 ) continue;       
    grCsI->Set(0);
    if( wConv->SetGraph( CsiModuleID, iSubMod ,conv , grCsI ) == 0 ){ continue; }		
    //////////////////////////////////////////////////////
    //// Different Analysis for each Different Module //// 
    //// For CsI Using Templete Fitting               //// 
    //// For other module Using Function fitting      //// 
    //////////////////////////////////////////////////////
    if( tempSpl[iSubMod] == NULL ){ 
      std::cout<< "Spline Pointer is NULL. Channel: "<< iSubMod  << std::endl;
    }else{
      // Init Waveform Analyzer // 
      wavAnalyzer->_Clear();
      Fitter->InitPar();

      int WavAnaRst = wavAnalyzer->AnalyzeWaveform( grCsI->GetY() );

      const double heightThreshold   = 30.;// Approx 3MeV // 
      const double widthThreshold    = 32.;// At least 3 point in peak 
      const double heightMinimum     = 5.;
      const double slopeDeltaMinimum = 50;
      if( wavAnalyzer->m_Height < heightMinimum ){
	// Abort too low signal // 
	continue; 
      }else if( wavAnalyzer->m_Height > heightThreshold ){
	// Selection for high signals 
	if( wavAnalyzer->m_Width < widthThreshold ){// case of spike data  || Tail/Head Peak Data// 
	  continue;
	}
      }else {
	// Selection for low signals 
	// Abort data too small in ADC 
	if( wavAnalyzer->m_SlopeDelta < slopeDeltaMinimum ){ continue; }
	// Abort data peak is not in Maximum SlopeRegion
	if(( wavAnalyzer->m_TimeMaximum < (wavAnalyzer->m_SlopeStart+3 )*8) ||
	   ( wavAnalyzer->m_TimeMaximum > (wavAnalyzer->m_SlopeStart+11)*8) ){
	  continue;
	}
      }

      Fitter->SetWaveform(tempSpl[iSubMod]);
      Fitter->SetParameter( wavAnalyzer );
      bool fit = Fitter->Fit(grCsI);      
      int chIndex = (wConv->mod[CsiModuleID])->m_nDigi;
      
      if( fit ){ 
	Double_t ResultHeight = Fitter->GetParameter(0);
	Double_t ResultEnergy = Converter->ConvertToEnergy( iSubMod , ResultHeight );
	//Abort small signal && Small Energy event 
	if( ResultHeight  < 5. || ResultEnergy < 1. ){ continue; }

	Double_t x,y;
	idHandler->GetMetricPosition( iSubMod , x, y ); 		
	int moduleIndex = wConv->mod[CsiModuleID]->m_nDigi;	

	wConv->mod[CsiModuleID]->m_FitHeight[moduleIndex]  = 1;
	wConv->mod[CsiModuleID]->m_ID[moduleIndex]         = iSubMod;
	wConv->mod[CsiModuleID]->m_Signal[moduleIndex]     = Fitter->GetParameter(0);
	wConv->mod[CsiModuleID]->m_Time[moduleIndex]       = Fitter->GetParameter(1);
	wConv->mod[CsiModuleID]->m_Pedestal[moduleIndex]   = Fitter->GetParameter(2);
	wConv->mod[CsiModuleID]->m_HHTime[moduleIndex]     = Fitter->GetConstantFraction();
	wConv->mod[CsiModuleID]->m_Chisq[moduleIndex]      = Fitter->GetChisquare();
	wConv->mod[CsiModuleID]->m_NDF[moduleIndex]        = Fitter->GetNDF();
	  //wConv->mod[CsiModuleID]->m_ADC[moduleIndex]    = Fitter->GetADC(gr);
	wConv->mod[CsiModuleID]->m_Energy[moduleIndex]   = Converter->ConvertToEnergy(iSubMod,Fitter->GetParameter(0));

	wConv->mod[CsiModuleID]->m_wav_SlopeDelta[moduleIndex]   = wavAnalyzer->m_SlopeDelta;
	wConv->mod[CsiModuleID]->m_wav_PeakTime[moduleIndex]     = wavAnalyzer->m_TimeMaximum;
	wConv->mod[CsiModuleID]->m_wav_Width[moduleIndex]        = wavAnalyzer->m_Width;
	wConv->mod[CsiModuleID]->m_wav_Height[moduleIndex]       = wavAnalyzer->m_Height;
	wConv->mod[CsiModuleID]->m_wav_FrontHalfTime[moduleIndex]= wavAnalyzer->m_BoundaryHead;
	wConv->mod[CsiModuleID]->m_wav_RearHalfTime[moduleIndex] = wavAnalyzer->m_BoundaryTail;
	wConv->mod[CsiModuleID]->m_wav_Pedestal[moduleIndex]     = wavAnalyzer->m_Pedestal;

	wConv->mod[CsiModuleID]->m_Fit_Pedestal[moduleIndex]     = Fitter->GetParameter(2);
	wConv->mod[CsiModuleID]->m_Fit_Height[moduleIndex]       = Fitter->GetParameter(0);
	wConv->mod[CsiModuleID]->m_Fit_Time[moduleIndex]         = Fitter->GetParameter(1);
	wConv->mod[CsiModuleID]->m_Fit_HHTime[moduleIndex]       = Fitter->GetConstantFraction();
	wConv->mod[CsiModuleID]->m_Fit_ChisqNDF[moduleIndex]     
	  = Fitter->GetChisqNDF(wConv->mod[CsiModuleID]->m_Fit_ChisqPed[moduleIndex],
				wConv->mod[CsiModuleID]->m_Fit_ChisqFront[moduleIndex],				
				wConv->mod[CsiModuleID]->m_Fit_ChisqRear[moduleIndex],
				wConv->mod[CsiModuleID]->m_Fit_ChisqTail[moduleIndex]);
	wConv->mod[CsiModuleID]->m_Conv_Energy[moduleIndex]      = ResultEnergy;
	wConv->mod[CsiModuleID]->m_Conv_Time[moduleIndex]        = Fitter->GetParameter(1);// Must change // 

	// Must Edit after // 
	// Large PMT signal delay ~ 16 ns than Small PMT
	if( wConv->mod[CsiModuleID]->m_Conv_Energy[moduleIndex] > 10 ){
	  if( wConv->mod[CsiModuleID]->m_ID[moduleIndex] <  2240 ){
	    m_TimeClusterHist->Fill( wConv->mod[CsiModuleID]->m_Conv_Time[ moduleIndex ]      ,
				     wConv->mod[CsiModuleID]->m_Conv_Energy[ moduleIndex ]);
	  }else{
	    m_TimeClusterHist->Fill( wConv->mod[CsiModuleID]->m_Conv_Time[ moduleIndex ] - 16 , 
				     wConv->mod[CsiModuleID]->m_Conv_Energy[ moduleIndex ]);
	  }
	}
	wConv->mod[CsiModuleID]->m_nDigi++;	    
	
	
      }else{
	;
      }
      Fitter->Clear();	
      //std::cout << wConv->mod[CsiModuleID]->m_nDigi << std::endl;
    }
  }
  
  // Clustering with Time // 
  wConv->mod[CsiModuleID]->m_nTimeCluster = 0;
  for( int i = 0; i< 32; i++){    
    wConv->mod[CsiModuleID]->m_TotalEnergyInTimeCluster[i] = 0;
    wConv->mod[CsiModuleID]->m_TimeClusterHead[i]       = -1;
    wConv->mod[CsiModuleID]->m_TimeClusterTail[i]       = -1;
  }
  for( int ibin = 0; ibin <= m_TimeClusterHist->GetNbinsX(); ibin++ ){
    // Set Total Energy of Time Cluster 
    if( m_TimeClusterHist->GetBinContent(ibin) != 0){
      wConv->mod[CsiModuleID]->m_TotalEnergyInTimeCluster[wConv->mod[CsiModuleID]->m_nTimeCluster] += m_TimeClusterHist->GetBinContent(ibin);      
    }
    // Set Cluster Start
    if( m_TimeClusterHist->GetBinContent( ibin -1 ) == 0 &&
	m_TimeClusterHist->GetBinContent( ibin    ) != 0 ){
      wConv->mod[CsiModuleID]->m_TimeClusterHead[ wConv->mod[CsiModuleID]->m_nTimeCluster ] = m_TimeClusterHist->GetBinLowEdge(ibin);
    }
    // Set Cluster End 
    if( m_TimeClusterHist->GetBinContent( ibin ) != 0 &&
	m_TimeClusterHist->GetBinContent( ibin + 1 ) == 0 ){
      wConv->mod[CsiModuleID]->m_TimeClusterTail[ wConv->mod[CsiModuleID]->m_nTimeCluster ] = m_TimeClusterHist->GetBinLowEdge(ibin) + m_TimeClusterHist->GetBinWidth( ibin );
      wConv->mod[CsiModuleID]->m_nTimeCluster++;
    }
  }
  for( int moduleIndex  = 0; moduleIndex < wConv->mod[CsiModuleID]->m_nDigi; moduleIndex++){
    wConv->mod[CsiModuleID]->m_TimeClusterID[moduleIndex] = -1;
    if( wConv->mod[CsiModuleID]->m_ID[moduleIndex] < 2240){
      for( int timeClusterIndex = 0; timeClusterIndex < wConv->mod[CsiModuleID]->m_nTimeCluster; timeClusterIndex++){
	if( wConv->mod[CsiModuleID]->m_Fit_Time[moduleIndex] > wConv->mod[CsiModuleID]->m_TimeClusterHead[timeClusterIndex]-8&&
	    wConv->mod[CsiModuleID]->m_Fit_Time[moduleIndex] < wConv->mod[CsiModuleID]->m_TimeClusterTail[timeClusterIndex]+8){
	  wConv->mod[CsiModuleID]->m_TimeClusterID[moduleIndex] = timeClusterIndex;
	  break;
	}
      }
    }else{
      for( int timeClusterIndex = 0; timeClusterIndex < wConv->mod[CsiModuleID]->m_nTimeCluster; timeClusterIndex++){
	if( wConv->mod[CsiModuleID]->m_Fit_Time[moduleIndex] > wConv->mod[CsiModuleID]->m_TimeClusterHead[timeClusterIndex]-8 + 16&&
	    wConv->mod[CsiModuleID]->m_Fit_Time[moduleIndex] < wConv->mod[CsiModuleID]->m_TimeClusterTail[timeClusterIndex]+8 + 16){
	  wConv->mod[CsiModuleID]->m_TimeClusterID[moduleIndex] = timeClusterIndex;
	  break;
	}
      }
    }
  }
  return (wConv->mod[CsiModuleID])->m_nDigi;
}
int  E14EventBuilder_V0::EventProcess(int ievent){
  m_EventNumber = ievent;
  wConv->InitData();
  for( int icosmic = 0; icosmic < 20; icosmic++){
    CosmicSignal[icosmic] = 0;
    CosmicTime[icosmic]   = 0;    
  }
  for( int icv = 0; icv < 10; icv++){
    CVSignal[icv] = 0;
    CVTime[icv]   = 0;
  }

  wConv->m_RunNo = m_RunNumber;
  wConv->m_EventNo = ievent;

  if( ievent % 100 == 0 && ievent ){ std::cerr << ievent << "/" << conv[0]->GetEntries() << "\n";  }

  for( int icrate = 0; icrate < nCrateFeb; icrate++){
    conv[icrate]->GetEntry(ievent);
  }
  Int_t Trigger  = TriggerDicision();
  Int_t nChannel = AnalyzeCsIData();
  std::cout << "TRIGGER:" << Trigger << "\tnChannel" << nChannel << std::endl;
  m_trOut->Fill();
  return nChannel;
}
int  E14EventBuilder_V0::LoopAll(){
  for( int ievent  = 0; ievent < m_Entries; ievent++){
  //for( int ievent  = 0; ievent < 1000; ievent++){
    if( (ievent % 100)  == 0 && ievent ){
      std::cout<< "Process:" << ievent << " / " << m_Entries << "\n";
    } 
    EventProcess( ievent );    
  }
  return m_Entries; 
}
void E14EventBuilder_V0::Clear(){
  TotalTriggerFlag = 0;
}
void E14EventBuilder_V0::DrawEvent(TCanvas* can){
  CsI_Module* csiEnergy = new CsI_Module("Energy");
  CsI_Module* csiTime   = new CsI_Module("Time");
  CsI_Module* csiEnergyCluster[2];
  for( int i  =0; i< 2; i++){
    csiEnergyCluster[i] = new CsI_Module(Form("Cluster%d",i));
  }

  //TH1D* hisTimeDistrib  = new TH1D("hisTimeDistrib","hisTimeDistrib",64,0,8*64);
  
  TH2D* hisEnergyTimeDistrib = new TH2D("hisEnergyTimeDistrib","EnergyTimeDistrib",64,0,8*64,160,0,1000);
  Int_t BigEnergyClusterID[2] = {-1};
  Double_t BigEnergyCluster[2] = {0};
  for( int iCluster  =0; iCluster < wConv->mod[CsiModuleID]->m_nTimeCluster; iCluster++){
    if( wConv->mod[CsiModuleID]->m_TotalEnergyInTimeCluster[ iCluster ] > BigEnergyCluster[0] ){
      BigEnergyCluster[0]   = wConv->mod[CsiModuleID]->m_TotalEnergyInTimeCluster[ iCluster ];
      BigEnergyClusterID[0] = iCluster;
    }else{
      if( wConv->mod[CsiModuleID]->m_TotalEnergyInTimeCluster[ iCluster ] > BigEnergyCluster[1] ){
	BigEnergyCluster[1]   = wConv->mod[CsiModuleID]->m_TotalEnergyInTimeCluster[ iCluster ];
	BigEnergyClusterID[1] = iCluster;
      }
    }
  }
  std::cout << BigEnergyClusterID[0] << "\t" << BigEnergyClusterID[1] << std::endl;
  for( int ich = 0; ich < wConv->mod[CsiModuleID]->m_nDigi; ich++){
    if( wConv->mod[CsiModuleID]->m_Energy[ich] > 1 ){
      /*
    std::cout << wConv->mod[CsiModuleID]->m_ID[ich] << "\t"
		<< wConv->mod[CsiModuleID]->m_TimeClusterID[ ich ] << "\t"
		<< wConv->mod[CsiModuleID]->m_Time[ich] << "\t" 
		<< wConv->mod[CsiModuleID]->m_Signal[ich] << "\t"
		<< wConv->mod[CsiModuleID]->m_Energy[ich] << "\n";
      */
      csiEnergy->Fill( wConv->mod[CsiModuleID]->m_ID[ich], wConv->mod[CsiModuleID]->m_Energy[ich]);
      csiTime  ->Fill( wConv->mod[CsiModuleID]->m_ID[ich], wConv->mod[CsiModuleID]->m_Time[ich]);

      if( wConv->mod[CsiModuleID]->m_TimeClusterID[ich] == BigEnergyClusterID[0] ){
	csiEnergyCluster[0]->Fill( wConv->mod[CsiModuleID]->m_ID[ich], wConv->mod[CsiModuleID]->m_Energy[ich]);
      }else if( wConv->mod[CsiModuleID]->m_TimeClusterID[ich] == BigEnergyClusterID[1] ){
	csiEnergyCluster[1]->Fill( wConv->mod[CsiModuleID]->m_ID[ich], wConv->mod[CsiModuleID]->m_Energy[ich]);
      }
      //if( wConv->mod[CsiModuleID]->m_Signal[ich] > 100 ){
      //hisTimeDistrib->Fill(wConv->mod[CsiModuleID]->m_Time[ich] - TimeOffset[wConv->mod[CsiModuleID]->m_ID[ich]]);
      //}
      hisEnergyTimeDistrib->Fill(wConv->mod[CsiModuleID]->m_Time[ich],wConv->mod[CsiModuleID]->m_Energy[ich]);
    }
  }
  can->Divide(2,3);
  can->cd(1);
  gPad->SetLogz();
  csiEnergy->Draw("colz");
  can->cd(2);
  csiTime->Draw("colz");
  can->cd(3);
  m_TimeClusterHist->Draw();
  //hisTimeDistrib->Draw();
  can->cd(4);
  hisEnergyTimeDistrib->Draw("col");
  can->cd(5);
  gPad->SetLogz();
  csiEnergyCluster[0]->Draw("colz");
  std::cout<< wConv->mod[CsiModuleID]->m_TimeClusterHead[BigEnergyClusterID[0]] << "\t"
	   << wConv->mod[CsiModuleID]->m_TimeClusterTail[BigEnergyClusterID[0]] << "\n";
  can->cd(6);
  gPad->SetLogz();
  csiEnergyCluster[1]->Draw("colz");
  std::cout<< wConv->mod[CsiModuleID]->m_TimeClusterHead[BigEnergyClusterID[1]] << "\t"
	   << wConv->mod[CsiModuleID]->m_TimeClusterTail[BigEnergyClusterID[1]] << "\n";
  
}
