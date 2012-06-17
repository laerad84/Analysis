#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

#include "E14ReadSumFile.h"
#include "IDHandler.h"

int
main( int argc, char **argv){
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("neRMIuo");
  
  if( argc != 2){ std::cerr << argv[0] << " [RunNUmber]" << std::endl; }

  int RunNumber = atoi(argv[1]);;


  std::string LaserFilenameStr = Form("ConversionedLaserData/LaserOut%d.root", 
				      RunNumber);
  std::string SumFilenameStr   = Form("sumup_data/Sum%d.root",RunNumber);
  std::string OutFilenameStr   = Form("mergeData/LaserMerged%d.root",RunNumber);
  // Condition 
  double  LaserThreshold = 1000.;

  IDHandler* handler  = new IDHandler();

  
  TFile* laserRootFile = new TFile(LaserFilenameStr.c_str());
  TTree* laserRootTree = (TTree*)laserRootFile->Get("LaserEvent");
  double Time;
  double PIN[4];
  double ADC[16];
  int    SCL[16];
  int    SpillFlag;

  laserRootTree->SetBranchAddress("Time",&Time);
  laserRootTree->SetBranchAddress("PIN",PIN);
  laserRootTree->SetBranchAddress("ADC",ADC);
  laserRootTree->SetBranchAddress("SCL",SCL);
  laserRootTree->SetBranchAddress("SpillFlag",&SpillFlag);

  //TFile* SumRootFile   = new TFile(SumFilenameStr.c_str());
  E14ReadSumFile* Reader = new E14ReadSumFile();
  Reader->Add(SumFilenameStr.c_str());


  TFile* tfout = new TFile(OutFilenameStr.c_str(), "recreate");
  TTree* trOut = new TTree("LaserOut","LaserOut vs Sum"); 
  int    nDigi;
  int    SumEntry;
  int    LaserEntry;
  int    ID[N_TOTAL_CSI];
  double SumOut[N_TOTAL_CSI];
  double PINOut[N_TOTAL_CSI];

  trOut->Branch("nDigi",&nDigi,"nDigi/I");
  trOut->Branch("SumEntry",&SumEntry,"SumEntry/I");
  trOut->Branch("LaserEntry",&LaserEntry,"LaserEntry/I");
  trOut->Branch("ID",ID,"ID[nDigi]/I");//nDigi
  trOut->Branch("SumOut",SumOut,"SumOut[nDigi]/D");//nDigi
  trOut->Branch("PINOut",PINOut,"PINOut[nDigi]/D");//nDigi

  
  long LaserEntries = laserRootTree->GetEntries();
  long SumEntries   = Reader->GetEntries(); 
  long preSumEntry = -1;
  for( int lentry = 0; lentry < LaserEntries; ++lentry){
    laserRootTree->GetEntry(lentry);

    LaserEntry = lentry;
    if( SpillFlag == 0 ){
      preSumEntry = SCL[3] -1;
      continue;
    }

    int tempSumEntry = SCL[3] -1;
    do{
      Reader->GetEntry(tempSumEntry);      
      SumEntry = tempSumEntry;
      --tempSumEntry;
      if( tempSumEntry < SCL[3] - 1 -5 ){ break;}
    }while( Reader->LaserEne[0] < LaserThreshold && tempSumEntry>preSumEntry);
    preSumEntry = SCL[3]-1;
    if( Reader->LaserEne[0] < LaserThreshold ){ 
      continue;
    }

    nDigi = Reader->CsiNumber;
    for( int idigi = 0; idigi < nDigi; ++idigi){

      SumOut[idigi] = Reader->CsiIntegratedADC[idigi];
      ID[idigi]     = Reader->CsiModID[idigi];
      double x,y;
      handler->GetMetricPosition(ID[idigi], x,y);
      if( x < 0 ){
	if( y> 0 ){
	  PINOut[idigi] = PIN[0];
	}else{
	  PINOut[idigi] = PIN[1];
	}
      }else{
	if( y> 0){
	  PINOut[idigi] = PIN[2];
	}else{
	  PINOut[idigi] = PIN[3];
	}	
      }
    }
    // Laser Data and Sum Data is matched //
    
    trOut->Fill();
  }
  trOut->Write();
  tfout->Close();
  return 0;
}

  

  



  


