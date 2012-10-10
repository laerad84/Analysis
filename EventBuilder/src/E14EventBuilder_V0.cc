#include "E14EventBuilder_V0.h"

E14EventBuilder_V0::E14EventBuilder_V0(TTree* trout, Int_t RunNumber){
  m_trOut     = trout;
  m_RunNumber = RunNumber;
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
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),m_RunNumber),
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
  
  TFile* tfTemplate = new TFile("TEMPLATE_SPLINE_250_500.root");
  for( int ispline = 0; ispline< 2716; ispline++){
    std::cout<< ispline << std::endl;
    tempGr[ispline]  = NULL;
    tempSpl[ispline] = NULL;
    //tempGr[ispline]  = (TGraphErrors*)tfTemplate->Get(Form("Waveform_Height_%d_0",ispline));
    tempGr[ispline]  = (TGraph*)tfTemplate->Get(Form("Template_graph_%d",ispline));
    std::cout<< ispline << std::endl;
    if( tempGr[ispline]->GetN()!= 0){
      tempSpl[ispline] = new TSpline3(Form("waveformSpl_%d",ispline),tempGr[ispline]);
    }else{
      std::cout<< "Non Exist channel:" << ispline << std::endl;
    }
  }
  tfTemplate->Close();
  return true;
}
bool E14EventBuilder_V0::InitTimeOffset(){
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
  
  if((wConv->ModMap[CVModuleID]).nMod     != 10)
    { 
      std::cout<< "CV nModule is not equal"     << std::endl;
    }
  if((wConv->ModMap[CosmicModuleID]).nMod != 20)
    { 
      std::cout<< "Cosmic nModule is not equal" << std::endl;
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
	
	TF1* linearFunction = new TF1("func","pol1",halfTiming - 12, halfTiming + 12);
	grTrigger->Fit( linearFunction, "Q", "", halfTiming -12, halfTiming +12 );
	double halfFitTiming = linearFunction->GetX( halfHeight, halfTiming -12, halfTiming +12);
	wConv->mod[iMod]->m_FitTime[chIndex]= halfFitTiming;
	delete linearFunction;
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

  int nSubCsiModule = wConv->GetNsubmodule( CsiModuleID );
  for( int iSubMod = 0; iSubMod < nSubCsiModule; iSubMod++){	
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
      Fitter->SetWaveform(tempSpl[iSubMod]);
      Fitter->InitPar();
      //std::cout<< "Fit" << std::endl;
      bool fit = Fitter->Fit(grCsI);
      int chIndex = (wConv->mod[CsiModuleID])->m_nDigi;
      if( fit ){ 
	if( Fitter->GetParameter(0) < 5 ){ continue ; }
	for( int ipoint = 0; ipoint < grCsI->GetN(); ipoint++){
	  Double_t x,y;
	  idHandler->GetMetricPosition( iSubMod , x, y ); 		
	  std::cout<< "Fit Result : " << Fitter->GetFitResult() << std::endl;
	  
	  wConv->mod[CsiModuleID]->m_FitHeight[wConv->mod[CsiModuleID]->m_nDigi]= 1;
	  wConv->mod[CsiModuleID]->m_ID[wConv->mod[CsiModuleID]->m_nDigi]       = iSubMod;
	  wConv->mod[CsiModuleID]->m_Signal[wConv->mod[CsiModuleID]->m_nDigi]   = Fitter->GetParameter(0);
	  wConv->mod[CsiModuleID]->m_Time[wConv->mod[CsiModuleID]->m_nDigi]     = Fitter->GetParameter(1);
	  wConv->mod[CsiModuleID]->m_Pedestal[wConv->mod[CsiModuleID]->m_nDigi] = Fitter->GetParameter(2);
	  wConv->mod[CsiModuleID]->m_HHTime[wConv->mod[CsiModuleID]->m_nDigi]   = Fitter->GetConstantFraction();
	  wConv->mod[CsiModuleID]->m_Chisq[wConv->mod[CsiModuleID]->m_nDigi]    = Fitter->GetChisquare();
	  wConv->mod[CsiModuleID]->m_NDF[wConv->mod[CsiModuleID]->m_nDigi]      = Fitter->GetNDF();
	  //wConv->mod[CsiModuleID]->m_ADC[wConv->mod[CsiModuleID]->m_nDigi]      = Fitter->GetADC(gr);
	  wConv->mod[CsiModuleID]->m_nDigi++;	    
	}
      }else{
	;
      }
      Fitter->Clear();	
      //std::cout << wConv->mod[CsiModuleID]->m_nDigi << std::endl;
    }
  }
  return wConv->mod[CsiModuleID]->m_nDigi;
}
int  E14EventBuilder_V0::EventProcess(int ievent){
  m_EventNumber = ievent;
  for( int icrate = 0; icrate < nCrateFeb; icrate++){
    conv[icrate]->GetEntry(ievent);
  }
  TriggerDicision();
  Int_t nChannel = AnalyzeCsIData();
  return nChannel;
}
int  E14EventBuilder_V0::LoopAll(){
  for( int ievent  = 0; ievent < m_Entries; ievent++){
    if( (ievent % 100)  == 0 && ievent ){
      std::cout<< "Process:" << ievent << " / " << m_Entries << "\n";
    } 
    EventProcess( ievent );    
  }
}
void E14EventBuilder_V0::Clear(){
  TotalTriggerFlag = 0;
}
