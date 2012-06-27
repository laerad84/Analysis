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
#include "Environment.h"

#include "E14ConvReader.h"
#include "E14ConvWriter.h"

const Double_t COSMIC_THRESHOLD[20] = {1000,1000,1000,1000,1000,
				       1000,1000,1000,1000,1000,
				       1000,1000,1000,1000,1000,
				       1000,1000,1000,1000,1000};

int
main(int argc, char** argv){
  gStyle->SetPalette(1);
  // ARGV[0]  <InputROOTFile> <OutputROOTFile> [<Vision>]
  const int CosmicArr[20] = {0 ,5 ,2 ,1 ,6 ,7 ,4 ,13,12,11,
  			     14,3 ,10,15,8 ,9 ,16,17,18,19};
  
  int RunNumber;
  std::string OutputFile;

  Double_t ADCThreshold    = 300;//in case of rawData 300;
  Double_t HeightThreshold = 20;// in case of Signal Height;~200 in ADC
  Double_t energyThreshold = 3;

  if( argc == 3){
    RunNumber = atoi(argv[1]);
    OutputFile = argv[2];
  }else{
    std::cerr << "ARGUEMENT ERROR" << std::endl;
    return -1;
  }
  
  // ReadEnvironment 
  int envRtn = GetEnvironment();
  PrintEnvironment();
  
  // Set InputFile
  TFile* tfRead = new TFile(Form("%s/run%d_wav.root",waveAnaFileDir,RunNumber));
  TTree* trRead = (TTree*)tfRead->Get("WFTree");
  E14ConvWriter* wConv = new E14ConvWriter(Form("%s/Sum%d.root",sumFileDir, RunNumber),
					   trRead);
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
  TFile* tfOut  = new TFile(OutputFile.c_str(),"Recreate");
  TTree* trOut  = new TTree("CosmicOut","");  

  TH2D* hisTriggerHitPosition
    = new TH2D("hisTriggerHitposition",
	       "TriggerHitPosition;XPosition[mm];ZPosition;nHit",
	       200,-1000,1000,
	       10,0,10);  
  
  Int_t    CosmicTrigger[20]={0};
  Int_t    UpperID[5];
  Int_t    DownID[5];
  Int_t    nHitUp  =  0;
  Int_t    nHitDn  =  0;
  Bool_t   bCosmicTrigger =  false;
  Int_t    TriggerIndex   = -1;
  Double_t CsIdepE[N_TOTAL_CSI];
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
  HoughCsI*  hough         = new HoughCsI();
  Chisq_cosmic* chi2Cosmic = new Chisq_cosmic();

  

  long nEntries = wConv->GetEntries();
  std::cout << "TotalEntries::" << nEntries  << std::endl;  
  
  TGraph* gr  = new TGraph();
  for( int iEntry = 0; iEntry < nEntries; iEntry++){
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
	  CsiTime[iSubMod]   = wConv->mod[iMod]->m_Timing[iSubMod];
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
    hough->Reset();

    // Set Activated Crystal to graph    
    for ( int idigi = 0; idigi < CsiNumber; idigi++){
      if( CsiSignal[idigi] < 20 ){ continue; }
      double x,y;
      handler->GetMetricPosition(CsiID[idigi], x,y);
      gr->SetPoint(gr->GetN(), x, y);      
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
	  HitUp |= 1 << iCosmic;
	}
	if( CosmicSignal[ iCosmic+5  ] > COSMIC_THRESHOLD[ iCosmic +5  ] &&
	    CosmicSignal[ iCosmic+15 ] > COSMIC_THRESHOLD[ iCosmic +15 ]){
	  HitCoinDn |= 1 << iCosmic;
	}
	if( CosmicSignal[ iCosmic+5  ] > COSMIC_THRESHOLD[ iCosmic +5  ] ||
	    CosmicSignal[ iCosmic+15 ] > COSMIC_THRESHOLD[ iCosmic +15 ]){
	  HitDn |= 1 << iCosmic;
	}       	  
      }
    }
  

  /*
    ///AnalysisCode    
    TGraph* gr = new TGraph();
    hough->Reset();
    
    for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
      if( Reader->CsiIntegratedADC[idigi] > ADCThreshold){//MeV
	double x,y;
	//std::cout << "ID::" << Reader->CsiModID[idigi] << std::endl;
	handler->GetMetricPosition(Reader->CsiModID[idigi], x, y);
	gr->SetPoint( gr->GetN(), x,y);	
      }
    }    
    nHitUp = 0; 
    nHitDn = 0; 
    
    //Trigger Judgement
    if( gr->GetN() < 350 ){      
      for( int icosmic = 0; icosmic< 5 ; icosmic++){	
	if( Reader->CosmicEne[icosmic]    > COSMIC_THRESHOLD[icosmic]  ||
	    Reader->CosmicEne[icosmic+10] > COSMIC_THRESHOLD[icosmic+10]){
	  UpperID[nHitUp] = icosmic;
	  CosmicBoolUp |= 1 << icosmic;
	  nHitUp++;	  
	}
	
	if( Reader->CosmicEne[icosmic+5]  > COSMIC_THRESHOLD[icosmic+5] ||
	    Reader->CosmicEne[icosmic+15] > COSMIC_THRESHOLD[icosmic+15]){	
	  DownID[nHitDn]  = icosmic%5;
	  nHitDn++;	  
	  CosmicBoolDn |= 1<< icosmic;
	}	
      }      
    }
    
    if( nHitDn  >= 1 && nHitUp >= 1){
      bCosmicTrigger = true;
      Trigger = 1;
    }else{
      bCosmicTrigger = false;
      Trigger = 0;
    }    
    nDigi = 0; 
    if( bCosmicTrigger){
      if(hough->CosmicJudgment(gr)){	
	gr->Set(0);
	CosmicFit = 1;
	roh       = hough->GetRoh();
	theta     = hough->GetTheta();
	
	for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){	
	  if( Reader->CsiIntegratedADC[idigi] > ADCThreshold ){
	    double x,y;
	    //std::cout << "ID::" << Reader->CsiModID[idigi] << std::endl;
	    handler->GetMetricPosition(Reader->CsiModID[idigi], x, y);
	    if( hough->CalDistance(x,y) <= 50 ){
	      gr->SetPoint( gr->GetN(), x,y);	    
	    }	
	  }
	}	
	
	chi2Cosmic->Reset();
	chi2Cosmic->SetFunction(gr);
	chi2Cosmic->SetRange(hough->GetRoh(),hough->GetTheta());
	chi2Cosmic->CalChisq();
	//std::cout << hough->GetRoh() << "\t" 
	//<< hough->GetTheta() << std::endl;
	roh = chi2Cosmic->GetRoh();
	theta = chi2Cosmic->GetTheta();
	
	if( nHitUp  > 0 && nHitDn > 0){
	  Double_t xPosOnScinti[2];//0:Up // 1:Down
	  Double_t RadTheta = TMath::Pi()/180*theta;
	  xPosOnScinti[0] = (roh - 1000*TMath::Sin(RadTheta))/TMath::Cos(RadTheta);
	  xPosOnScinti[1] = (roh + 1000*TMath::Sin(RadTheta))/TMath::Cos(RadTheta);
	  for( int ipos = 0; ipos < nHitUp ; ipos++){	    
	    hisTriggerHitPosition->Fill(xPosOnScinti[0],UpperID[ipos]);
	  }
	  for( int ipos = 0; ipos < nHitDn; ipos++){	    
	    hisTriggerHitPosition->Fill(xPosOnScinti[1],DownID[ipos]+5);
	  }
	}


	for( int idigi = 0; idigi<Reader->CsiNumber; idigi++){	  
	  double x,y;
	  if( Reader->CsiIntegratedADC[idigi] > ADCThreshold ){
	    handler->GetMetricPosition(Reader->CsiModID[idigi],x, y);
	    Double_t length;
	    if( Reader->CsiModID[idigi] < 2240){
	      length = 25.;
	    }else{
	      length = 50.;
	    }
	    if( chi2Cosmic->GetDistance(x,y) < length/2.*(TMath::Cos(theta*TMath::Pi()/180.)-TMath::Sin(theta*TMath::Pi()/180.)) && TMath::Abs(theta) < 45){
	      CsIID[nDigi]   = Reader->CsiModID[idigi];
	      CsIdepE[nDigi] = Reader->CsiEne[idigi];
	      CsIADC[nDigi]  = Reader->CsiIntegratedADC[idigi];
	      //hisCosmicHit[Reader->CsiModID[idigi]]->Fill(Reader->CsiIntegratedADC[idigi]);
	      nDigi++;
	    }
	  }
	}
	CalFactor = chi2Cosmic->GetCalibrationFactor();      

      }else{	
	CalFactor = 0;
	nDigi     = 0;
	roh       = 0;
	theta     = 0;
      }
      
    }else{
      CalFactor = 0;
      nDigi     = 0;
      roh       = 0; 
      theta     = 0;
    }
    trOut->Fill();
  }

  */

  std::cout << "Analysis cosmic ray event is over" << std::endl;
  hisTriggerHitPosition->Write();
  trOut->Write();
  tfOut->Close();  
}
