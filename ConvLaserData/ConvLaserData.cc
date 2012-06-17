#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"


int
main( int argc, char **argv){
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("neRMIuo");

  if( argc != 2){
    std::cout<< argv[0] << " [ RUNNUMBER ] "  << std::endl;
    return -1;
  }
  int RunNumber = atoi(argv[1]);//3893;

  std::string LaserFilenameStr = Form( "LaserData/laserdata_%d.root", RunNumber);
  std::string OutFilenameStr   = Form( "ConversionedLaserData/LaserOut%d.root",RunNumber);
  
  ///////////////////////////////////////////////////////////////////////////////
  // INPUT
  ///////////////////////////////////////////////////////////////////////////////

  TFile* laserRootFile = new TFile(LaserFilenameStr.c_str());
  TTree* laserTree     = (TTree*)laserRootFile->Get("laserEventData");
  
  Double_t T;
  Int_t ADCData[16];
  Int_t scaleData[16];
  Int_t TDCData[16];
  Int_t PINData[4];

  laserTree->SetBranchAddress("Time",&T);
  laserTree->SetBranchAddress("ADCData",ADCData);
  laserTree->SetBranchAddress("TDCData",TDCData);
  laserTree->SetBranchAddress("scaleData",scaleData);
  laserTree->SetBranchAddress("PINData",PINData);

  ///////////////////////////////////////////////////////////////////////////////
  // OUTPUT
  ///////////////////////////////////////////////////////////////////////////////
  
  TFile* laserOutFile = new TFile(OutFilenameStr.c_str(),"recreate");
  TTree* laserRecTree = new TTree("LaserEvent","Processed LaserData");
  Double_t Time;
  Double_t PIN[4];
  Double_t ADC[16];
  Double_t TDC[16];
  Int_t    SCL[16];
  Int_t    Spill;
  Int_t    SpillFlag;
  Int_t    SyncFlag;
  laserRecTree->Branch("Time",&Time,"Time/D");
  laserRecTree->Branch("PIN",PIN,"PIN[4]/D");
  laserRecTree->Branch("ADC",ADC,"ADC[16]/D");
  laserRecTree->Branch("TDC",TDC,"TDC[16]/D");
  laserRecTree->Branch("SCL",SCL,"SCL[16]/I");
  laserRecTree->Branch("Spill",&Spill,"Spill/I");
  laserRecTree->Branch("SpillFlag",&SpillFlag,"SpillFlag/I");
  laserRecTree->Branch("SyncFlag",&SyncFlag,"SyncFlag/I");

  ///////////////////////////////////////////////////////////////////////////////
  // Convert
  ///////////////////////////////////////////////////////////////////////////////

  Double_t timeOut[2];
  Double_t pinOut[2][4];
  Double_t ADCOut[2][16];
  Int_t sclOut[2][16];
  Double_t pinSum[2];
  for( int i = 0; i<2; ++i){
    for( int j = 0; j<4; ++j){
      pinOut[i][j] = 0; 
    }
    for( int j = 0; j< 16; ++j){
      ADCOut[i][j] = 0;
      sclOut[i][j] = 0;
    }
  }    


  long nentries = laserTree->GetEntries();
  Spill = 0; 
  SpillFlag = 0;
  Int_t LastEvent =0; 
  Int_t ScaleData = 0;
  for( int ientry = 0; ientry < nentries ; ++ientry){
    //------------------------------------
    //init Variables //
    //------------------------------------
    pinSum[0] = 0;
    pinSum[1] = 0;
    if( ientry <= 2 || ientry + 1 == nentries){
      continue; 
    }
    laserTree->GetEntry(ientry);
    timeOut[0] = T;
    for( int ipin = 0; ipin< 4; ipin++){
      pinOut[0][ipin] = PINData[ipin];
      pinSum[0]+=PINData[ipin];
    }
    for( int ich  = 0; ich < 16; ++ich){
      ADCOut[0][ich] = ADCData[ich];
      sclOut[0][ich] = scaleData[ich];
    }
    laserTree->GetEntry(ientry+1);
    timeOut[1] = T;
    for( int ipin = 0; ipin< 4; ipin++){
      pinOut[1][ipin] = PINData[ipin];
      pinSum[1] += PINData[ipin];
    }
    for( int ich  = 0; ich < 16; ++ich){
      ADCOut[1][ich] = ADCData[ich];
      sclOut[1][ich] = scaleData[ich];
    }

    //----------------------------------
    // Define Case
    //----------------------------------
    SyncFlag = -2;
    if( pinSum[0] > pinSum[1] ){
      if( ADCOut[0][1] > ADCOut[1][1] ){ // Good Sync Case
	Time = timeOut[0];
	SyncFlag = 0;
	for( int ipin = 0; ipin < 4; ++ipin){
	  PIN[ipin] = pinOut[0][ipin] - pinOut[1][ipin];
	}
	for( int ich  = 0; ich < 16; ++ich){
	  ADC[ich] = ADCOut[0][ich] - ADCOut[1][ich];
	  SCL[ich] = sclOut[0][ich];
	}
	
      }else{// Odd Sync Case	
	Time = timeOut[0];
	SyncFlag = -1;
	for( int ipin = 0; ipin < 4; ++ipin){
	  PIN[ipin] = pinOut[0][ipin] - pinOut[1][ipin];
	}
	for( int ich  = 0; ich < 16; ++ich){
	  ADC[ich] = ADCOut[1][ich] - ADCOut[0][ich];
	  SCL[ich] = sclOut[1][ich];
	}
      }
      if( SyncFlag == 0){
	++ientry;// Because Find SyncData.
      }
    }else{// Not Using
      continue;
    }
    if( LastEvent == SCL[3] ){
      SpillFlag = 0;
    }else{
      if(SpillFlag ==0){++Spill;}
      SpillFlag = 1;
    }

    // mask some unuse channel
    for( int ich = 0; ich < 16; ++ich){
      if( ich != 1){
	ADC[ich] = 0; 
      }
      if( ich > 4 ){
	SCL[ich] = 0; 
      }
      TDC[ich] = 0; 
    }

    //-----------------------------------
    // Using Only PIN ADC SCL TDC //


    


    //Doing After Event Processing
    laserRecTree->Fill();
    LastEvent = SCL[3];
  }
  laserRecTree->Write();
  laserOutFile->Close();
}

  

  



  


