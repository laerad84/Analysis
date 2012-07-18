#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TFile.h"

#include "GeneralTypes.h"
#include "GeneralMacros.h"

#include "E14ReadSumFile.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "HoughCsI.h"
#include "Chisq_cosmic.h"

#include "E14MapReader.h"
#include "Structs.h"
#include "E14ConvReader.h"
#include "E14ConvWriter.h"
#include "E14CosmicAnalyzer.h"


const Double_t COSMIC_THRESHOLD[20] = {1000,1000,1000,1000,1000,
				       1000,1000,1000,1000,1000,
				       1000,1000,1000,1000,1000,
				       1000,1000,1000,1000,1000};

int
main(int argc, char** argv){
  
  gStyle->SetOptFit ( 11111111  );
  gStyle->SetOptStat( "neMRiuo" );  
  gStyle->SetPalette(1);

  // ARGV[0]  <InputROOTFile> <OutputROOTFile> [<Vision>]

  const int CosmicArr[20]= {4 ,5 ,2 ,3 ,6 ,7 ,0 ,1 ,12,13,
			    10,11,14,15,8 ,9 ,16,17,18,19};
  int RunNumber;
  std::string OutputFile;

  Double_t ADCThreshold    = 300;//in case of rawData 300;
  Double_t HeightThreshold = 30;// in case of Signal Height;~200 in ADC
  Double_t energyThreshold = 3;

  if( argc == 2){
    RunNumber = atoi(argv[1]);
    //OutputFile = argv[2];
  }else{
    std::cerr << "ARGUEMENT ERROR" << std::endl;
    return -1;
  }
  
  // ReadEnvironment
  std::cout << "ENVIRONMENT" << std::endl;
  std::string ANALIBDIR     = std::getenv("ANALYSISLIB"    );
  std::string CONVFILEDIR   = std::getenv("ROOTFILE_CONV"  );
  std::string WAVEFILEDIR   = std::getenv("ROOTFILE_WAV"   );
  std::string SUMFILEDIR    = std::getenv("ROOTFILE_SUMUP" );
  std::string COSMICFILEDIR = std::getenv("ROOTFILE_COSMIC");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;

  // Set InputFile
  std::cout << "FILE SETTING" << std::endl;
  TFile* tfRead = new TFile(Form("%s/run%d_wav.root",WAVEFILEDIR.c_str(), RunNumber));
  TTree* trRead = (TTree*)tfRead->Get("WFTree");
  E14ConvWriter* wConv = new E14ConvWriter(Form("%s/Sum%d.root",SUMFILEDIR.c_str(), RunNumber), trRead);
  {
    
    //  ID 0 : Csi; 1 : Cosmic; 2 : Laser
    wConv->AddModule("Csi");
    wConv->AddModule("Cosmic"); 
    wConv->AddModule("Laser");
    wConv->Set();
    wConv->SetMap();
    wConv->SetBranchAddress();
    
  }
  
  // Set OutputFile
  TFile* tfOut  = new TFile(Form("%s/run%d_cosmic.root",COSMICFILEDIR.c_str(),RunNumber),"Recreate");
  TTree* trOut  = new TTree("CosmicOut","");  

  TH2D* hisTriggerHitPosition
    = new TH2D("hisTriggerHitposition",
	       "TriggerHitPosition;XPosition[mm];ZPosition;nHit",
	       200,-1000,1000,
	       10 ,0    ,10);  
  
  Int_t    CosmicTrigger[20]={0};
  Int_t    UpperID[5];
  Int_t    DownID[5];
  Int_t    nHitUp         =  0;
  Int_t    nHitDn         =  0;
  Bool_t   bCosmicTrigger =  false;
  Int_t    TriggerIndex   = -1;

  Double_t CsIdepE[N_TOTAL_CSI];
  Double_t CsITiming[N_TOTAL_CSI];
  Double_t CsIHHTiming[N_TOTAL_CSI];
  Double_t CsISplTiming[N_TOTAL_CSI];
  Double_t CsIFitTiming[N_TOTAL_CSI];

  Int_t    CsIID[N_TOTAL_CSI];
  Double_t CsIADC[N_TOTAL_CSI];
  Double_t PathLength[N_TOTAL_CSI];
  Double_t CalFactor;
  Int_t    nDigi;
  Int_t    Trigger;
  Double_t roh;
  Double_t theta;
  Int_t    CosmicFit;
  Int_t    CosmicBoolUp;
  Int_t    CosmicBoolDn;

  Int_t    HitCoinUp;
  Int_t    HitCoinDn;
  Int_t    HitUp;
  Int_t    HitDn;
  {

    trOut->Branch("nDigi"    ,&nDigi    ,"nDigi/I");
    trOut->Branch("CsIID"    ,CsIID     ,"CsIID[nDigi]/I");//nDigi;
    trOut->Branch("CsIADC"   ,CsIADC    ,"CsIADC[nDigi]/I");//nDigi;
    trOut->Branch("CsIdepE"  ,CsIdepE   ,"CsidepE[nDigi]/D");//nDigi;
    trOut->Branch("CsITiming"    ,CsITiming      ,"CsITiming[nDigi]/D");//nDigi;
    trOut->Branch("CsIHHTiming"  ,CsIHHTiming    ,"CsIHHTiming[nDigi]/D");//nDigi;
    trOut->Branch("CsIFitTiming" ,CsIFitTiming   ,"CsIFitTiming[nDigi]/D");//nDigi;
    trOut->Branch("CsISplTiming" ,CsISplTiming   ,"CsISplTiming[nDigi]/D");//nDigi;

    trOut->Branch("nHitUp"   ,&nHitUp   ,"nHitUp/I");
    trOut->Branch("nHitDn"   ,&nHitDn   ,"nHitDn/I");
    
    trOut->Branch("HitCoinUp",&HitCoinUp,"HitCoinUp/I");
    trOut->Branch("HitCoinDn",&HitCoinDn,"HitCoinDn/I");
    
    trOut->Branch("HitUp"   ,&HitUp   ,"HitUp/I");
    trOut->Branch("HitDn"   ,&HitDn   ,"HitDn/I");
    
    trOut->Branch("UpperID"  ,&UpperID  ,"UpperID[nHitUp]/I");//nHitUp;
    trOut->Branch("DownID"   ,&DownID   ,"DownID[nHitDn]/I");//nHitDn;

    trOut->Branch("Trigger"  ,&Trigger  ,"Trigger/I");
    trOut->Branch("CalFactor",&CalFactor,"CalFactor/D");
    trOut->Branch("roh"      ,&roh      ,"roh/D");
    trOut->Branch("theta"    ,&theta    ,"theta/D");
    trOut->Branch("CosmicFit",&CosmicFit,"CosmicFit/I");
    trOut->Branch("CosmicBoolUp",&CosmicBoolUp,"CosmicBoolUp/I");
    trOut->Branch("CosmicBoolDn",&CosmicBoolDn,"CosmicBoolDn/I");
    //trOut->Branch("PathLength"  ,&PathLength  ,"PathLength[nDigi]/D");
  }


  int    CsiNumber;
  int    CsiID[4096];
  double CsiSignal[4096];
  double CsiTime[4096];
  double CsiHHTime[4096];
  double CsiFitTime[4096];
  double CsiSplTime[4096];
  int    CosmicNumber;
  int    CosmicID[4096];
  double CosmicSignal[4096];
  double CosmicTime[4096];
  int    LaserNumber;
  int    LaserID[4096];
  double LaserSignal[4096];
  double LaserTime[4096];

  double CosmicOut[20];

  IDHandler* handler       = new IDHandler();
  /*
    HoughCsI*  hough         = new HoughCsI();
    Chisq_cosmic* chi2Cosmic = new Chisq_cosmic();
  */
  E14CosmicAnalyzer* CosmicAna = new E14CosmicAnalyzer();
  
  long nEntries = wConv->GetEntries();
  std::cout << "TotalEntries::" << nEntries  << std::endl;  
  
  TGraph* gr  = new TGraph();
  for( int iEntry = 0; iEntry < nEntries; iEntry++){
    bool Trigger=false;
    bool CoinTrigger=false;
    
    wConv->GetEntry(iEntry);
    if( iEntry && iEntry % 100 == 0){
      std::cout<< "\r" << iEntry << "/" << nEntries << std::endl;
      std::cout<< std::flush; 
    }

    //Init
    { 
      CosmicBoolUp = 0;
      CosmicBoolDn = 0;
      nDigi        = 0; 
      nHitUp       = 0; 
      nHitDn       = 0; 
      Trigger      = 0; 
      CalFactor    = 0;
      roh          = 0;
      theta        = 0;
      for( int i = 0; i< 5; i++){
	UpperID[i] = 0;
	DownID[i]  = 0; 
      }

      HitUp = 0;
      HitDn = 0; 
      HitCoinUp = 0;
      HitCoinDn = 0; 

      for( int i = 0;  i < N_TOTAL_CSI; i++ ){
	CsIdepE[i]      = 0.;
	CsIID[i]        = -1; 
	CsITiming[i]    = -9999;
	CsIHHTiming[i]  = -9999;
	CsIFitTiming[i] = -9999;
	CsISplTiming[i] = -9999;
	
      }
    }


    for( int i = 0; i< 20; i++){
      CosmicSignal[i] = 0;
      CosmicTime[i]   = 0;
    }
    
    // Cosmic Trigger Judge // 
    for( int iMod = 0; iMod < 3; iMod++ ){
      int nSubMod = (wConv->mod[iMod])->m_nDigi;      
      switch( iMod ){
      case 0:
	// Csi
	CsiNumber = nSubMod;
	for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++ ){	    
	  CsiID[iSubMod]     = wConv->mod[iMod]->m_ID[iSubMod];
	  CsiSignal[iSubMod] = wConv->mod[iMod]->m_Signal[iSubMod];
	  CsiTime[iSubMod]      = wConv->mod[iMod]->m_Timing[iSubMod];
	  CsiFitTime[iSubMod]   = wConv->mod[iMod]->m_FitTiming[iSubMod];
	  CsiSplTime[iSubMod]   = wConv->mod[iMod]->m_SplTiming[iSubMod];
	  CsiHHTime[iSubMod]    = wConv->mod[iMod]->m_HHTiming[iSubMod];

	}
	break;
      case 1:
	CosmicNumber = nSubMod;
	for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++ ){	    
	  //CosmicID[iSubMod]     = wConv->mod[iMod]->m_ID[iSubMod];
	  CosmicSignal[CosmicArr[wConv->mod[iMod]->m_ID[iSubMod]]] = wConv->mod[iMod]->m_Signal[iSubMod];
	  CosmicTime[CosmicArr[wConv->mod[iMod]->m_ID[iSubMod]]]   = wConv->mod[iMod]->m_Timing[iSubMod];
	}
	break;
      case 2:
	LaserNumber = nSubMod;
	for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++ ){	    
	  LaserID[iSubMod]     = wConv->mod[iMod]->m_ID[iSubMod];
	  LaserSignal[iSubMod] = wConv->mod[iMod]->m_Signal[iSubMod];
	  LaserTime[iSubMod]   = wConv->mod[iMod]->m_Timing[iSubMod];
	}
	break;
      default:
	break;
      }

    }

    
    // Analysis Code 
    gr->Set(0);

    // Set Activated Crystal to graph    
    
    nDigi = 0;
    for ( int idigi = 0; idigi < CsiNumber; idigi++){
      if( CsiSignal[idigi] < HeightThreshold ){ continue; }
      double x,y;
      handler->GetMetricPosition(CsiID[idigi], x,y);
      gr->SetPoint(gr->GetN(), x, y);

      CsIID[nDigi]   = CsiID[idigi];
      CsIdepE[nDigi]      = CsiSignal[idigi];
      CsITiming[nDigi]    = CsiTime[idigi];
      CsIHHTiming[nDigi]  = CsiHHTime[idigi];
      CsIFitTiming[nDigi] = CsiFitTime[idigi];
      CsISplTiming[nDigi] = CsiSplTime[idigi];

      nDigi++;
    }
    
    // Trigger Judgement // 
    nHitUp = 0; 
    nHitDn = 0;
    if( LaserNumber == 0 ){
      for( int iCosmic = 0; iCosmic< 5; iCosmic++){
	if( CosmicSignal[ iCosmic    ] > COSMIC_THRESHOLD[ iCosmic     ] &&
	    CosmicSignal[ iCosmic+10 ] > COSMIC_THRESHOLD[ iCosmic +10 ]){
	  HitCoinUp |= 1 << iCosmic;
	}
	if( CosmicSignal[ iCosmic    ] > COSMIC_THRESHOLD[ iCosmic     ] ||
	    CosmicSignal[ iCosmic+10 ] > COSMIC_THRESHOLD[ iCosmic +10 ]){
	  HitUp     |= 1 << iCosmic;
	}
	if( CosmicSignal[ iCosmic+5  ] > COSMIC_THRESHOLD[ iCosmic +5  ] &&
	    CosmicSignal[ iCosmic+15 ] > COSMIC_THRESHOLD[ iCosmic +15 ]){
	  HitCoinDn |= 1 << iCosmic;
	}
	if( CosmicSignal[ iCosmic+5  ] > COSMIC_THRESHOLD[ iCosmic +5  ] ||
	    CosmicSignal[ iCosmic+15 ] > COSMIC_THRESHOLD[ iCosmic +15 ]){
	  HitDn     |= 1 << iCosmic;
	}       	  
      }
    }
    if( HitUp != 0 && HitDn != 0 ){ Trigger = true; }
    if( HitCoinUp != 0 && HitCoinDn != 0 ){ CoinTrigger = true ; } 
    if( Trigger ){ 
      CosmicAna->Reset();
      if( CosmicAna->GetResult( gr, roh, theta ) ){
	CalFactor = CosmicAna->mc_chi2Cosmic->GetCalibrationFactor();      
	CosmicFit = 1; 
      }else{
	CosmicFit = 0;
	CalFactor = 0;
      }
    }else{	
      CalFactor = 0;
      nDigi     = 0;
      roh       = 0;
      theta     = 0;
    }  
    trOut->Fill();
  }
  

  std::cout << "Analysis cosmic ray event is over" << std::endl;
  hisTriggerHitPosition->Write();
  trOut->Write();
  tfOut->Close();  
}



