#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"


int
main( int argc, char** argv){
  int RunNumber = atoi( argv[1] );

  std::string ROOTFILEWAV = std::getenv("ROOTFILE_WAV");
  TFile* tf = new TFile(Form("%s/TEMPLATE_FIT_RESULT_DUMP_%d.root",ROOTFILEWAV.c_str(), RunNumber));
  TTree* tr =(TTree*)tf->Get("Waveform");
  
  Double_t Waveform[48];
  Short_t  TimeInfo[48];
  Int_t    ModuleNumber;
  Int_t    EventNumber;
  Double_t ChisqNDF;
  Double_t PeakTime;
  Double_t HHTime;
  Double_t Height;
  Double_t Pedestal;

  tr->SetBranchAddress("Waveform",Waveform);
  tr->SetBranchAddress("TimeInfo",TimeInfo);
  tr->SetBranchAddress("ModuleNumber",&ModuleNumber);
  tr->SetBranchAddress("EventNumber",&EventNumber);
  tr->SetBranchAddress("ChisqNDF",&ChisqNDF);
  tr->SetBranchAddress("PeakTime",&PeakTime);
  tr->SetBranchAddress("HHTime",&HHTime);
  tr->SetBranchAddress("Height",&Height);
  tr->SetBranchAddress("Pedestal",&Pedestal);

  TGraph* grTemp = new TGraph();
  
  Double_t Chisq[2716];
  Int_t    EventFirstRunNumber;
  Int_t    EventLastRunNumber;
  
  TFile* tfout = new TFile(Form("%s/TEMPALTE_FIT_REUSLT_REBUILD_%d.root",ROOTFILEWAV.c_str(), RunNumber),"RECREATE");
  TTree* trout = new TTree("trout","");
  Int_t    nWave;
  Int_t    ID[2716];
  Double_t WaveData[2716][48];
  Short_t  TimeData[2716][48];
  Double_t Output[2716];
  Double_t Peak[2716];
  Double_t HHT[2716];
  Double_t Ped[2716];

  trout->Branch("nWave"   ,&nWave  ,"nWave/I");
  trout->Branch("ID"      ,ID      ,"ID/I");//nWave
  trout->Branch("WaveData",WaveData,"WaveData[nWave][48]/D");//nWave
  trout->Branch("TimeData",TimeData,"TimeData[nWave][48]/S");//nWave
  trout->Branch("Output"  ,Output  ,"Output[nWave]/D");//nWave
  trout->Branch("Ped"     ,Ped     ,"Ped[nWave]/D");//nWave
  trout->Branch("Peak"    ,Peak    ,"Peak[nWave]/D");//nWave
  trout->Branch("HHT"     ,HHT     ,"HHT[nWave]/D");//nWave
  
 
  Int_t OldEventNumber = -1;
  Int_t nEntries = tr->GetEntries();
  std::cout << nEntries << std::endl;
  for( int ievent = 0; ievent < nEntries; ievent++){
    tr->GetEntry( ievent );
    grTemp->Set(0);
    
    if( EventNumber    != OldEventNumber ||
	OldEventNumber <  0              ){
      if( OldEventNumber != OldEventNumber){
	trout->Fill();
      }
      nWave = 0;      
      for( int ich = 0; ich < 2716; ich++){
	for( int ipoint = 0; ipoint < 48; ipoint++){
	  WaveData[ich][ipoint] = 0;
	  TimeData[ich][ipoint] = 0;
	}
	Output[ich] = 0;
	Peak[ich]   = 0;
	HHT[ich]    = 0;
	Ped[ich]    = 0;
	ID[ich]     = 0;
      }
    }
    
    ID[nWave] = ModuleNumber;
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      grTemp->SetPoint( ipoint , TimeInfo[ipoint], Waveform[ipoint]);
      WaveData[nWave][ipoint] = Waveform[ipoint];
      TimeData[nWave][ipoint] = TimeInfo[ipoint];      
    }
    
    nWave++;
    OldEventNumber = EventNumber;
  }

  trout->Fill();

  trout->Write();
  tfout->Close();
  return 0; 
}
  
  
