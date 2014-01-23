#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "csimap/CsiMap.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "CosmicTriggerTree.h"

int main( int argc, char** argv ){

  if( argc != 2 ) { 
    std::cout<< "Wrong argument" << std::endl;
    return -1; 
  }
  Int_t RunNumber = std::atoi( argv[1]);

  std::string IFDIR = std::getenv("ROOTFILEWAV");
  std::string OFDIR = std::getenv("ROOTFILEBEAM");

  char* InputFilename = Form("%s/run_wav_%d.root",IFDIR.c_str(),RunNumber);
  char* OutputFilename= Form("%s/run_wav_Beam_%d.root",OFDIR.c_str(),RunNumber);

  CsiMap* map = CsiMap::getCsiMap();


  CosmicTriggerTree* cosmicTrig = new CosmicTriggerTree();

  TFile* tf = new TFile(InputFilename);
  TTree* tr = (TTree*)tf->Get("Tree");
  cosmicTrig->SetBranchAddress(tr);
  Int_t EventNo;

  Int_t CsiNumber;
  Double_t CsiSignal[2716];
  Double_t CsiEne[2716];
  Double_t CsiTime[2716];
  Double_t CsiHHTime[2716];
  Double_t CsiChisq[2716];
  Short_t  CsiID[2716];
  Short_t  CsiNDF[2716];
  Short_t  CsiPosID[2716];
  Short_t  CsiGB[2716];
  Short_t  CsiCrate[2716];

  Int_t LaserNumber;
  Double_t LaserSignal[10];
  Double_t LaserHHTime[10];
  Double_t LaserTime[10];
  Double_t LaserChisq[10];
  Short_t  LaserID[10];
  Short_t  LaserNDF[10];
  
  Int_t OEVNumber;
  Double_t OEVSignal[50];
  Double_t OEVTime[50];
  Double_t OEVChisq[50];
  Short_t  OEVID[50];
  Short_t  OEVNDF[50];

  Int_t CC03Number;
  Double_t CC03Signal[50];
  Double_t CC03Time[50];
  Double_t CC03Chisq[50];
  Short_t  CC03ID[50];
  Short_t  CC03NDF[50];

  Int_t CosmicNumber;
  Double_t CosmicSignal[20];
  Double_t CosmicTime[20];
  Double_t CosmicChisq[20];
  Short_t  CosmicID[20];
  Short_t  CosmicNDF[20];
  
  Int_t EtcNumber;
  Double_t EtcSignal[10];
  Double_t EtcTime[10];
  Double_t EtcHHTime[10];
  Double_t EtcChisq[10];
  Short_t  EtcID[10];
  Short_t  EtcNDF[10];

  tr->SetBranchAddress("EventNo",&EventNo);

  tr->SetBranchAddress("CsiNumber",&CsiNumber);  
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  tr->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  tr->SetBranchAddress("CsiHHTime",CsiHHTime);//CsiNumber
  tr->SetBranchAddress("CsiChisq",CsiChisq);//CsiNumber
  tr->SetBranchAddress("CsiID",CsiID);//CsiNumber
  tr->SetBranchAddress("CsiNDF",CsiNDF);//CsiNumber

  tr->SetBranchAddress("LaserNumber",&LaserNumber);
  tr->SetBranchAddress("LaserSignal",LaserSignal);//LaserNumber
  tr->SetBranchAddress("LaserHHTime",LaserHHTime);//LaserNumber
  tr->SetBranchAddress("LaserTime",LaserTime);//LaserNumber
  tr->SetBranchAddress("LaserChisq",LaserChisq);//LaserNumber
  tr->SetBranchAddress("LaserID",LaserID);//LaserNumber
  tr->SetBranchAddress("LaserNDF",LaserNDF);//LaserNumber

  tr->SetBranchAddress("EtcNumber",&EtcNumber);
  tr->SetBranchAddress("EtcSignal",EtcSignal);//EtcNumber
  tr->SetBranchAddress("EtcHHTime",EtcHHTime);//EtcNumber
  tr->SetBranchAddress("EtcTime",EtcTime);//EtcNumber
  tr->SetBranchAddress("EtcChisq",EtcChisq);//EtcNumber
  tr->SetBranchAddress("EtcID",EtcID);//EtcNumber
  tr->SetBranchAddress("EtcNDF",EtcNDF);//EtcNumber

  TFile* tfOut = new TFile(OutputFilename,"recreate");
  TTree* trOut = new TTree("EventTree","Beam Data");
  cosmicTrig->Branch(trOut);
  trOut->Branch("RunNo",&RunNumber,"RunNumber/I");
  trOut->Branch("EventNo",&EventNo,"EventNo/I");

  trOut->Branch("CsiNumber",&CsiNumber,"CsiNumber/I");
  trOut->Branch("CsiSignal",CsiSignal,"CsiSignal[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiEne",CsiEne,"CsiEne[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiTime",CsiTime,"CsiTime[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiHHTime",CsiHHTime,"CsiHHTime[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiChisq",CsiChisq,"CsiChisq[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiNDF",CsiNDF,"CsiNDF[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiID",CsiID,"CsiID[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiPosID",CsiPosID,"CsiPosID[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiGB",CsiGB,"CsiGB[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiCrate",CsiCrate,"CsiCrate[CsiNumber]/S");//CsiNumber

  trOut->Branch("LaserNumber",&LaserNumber,"LaserNumber/I");
  trOut->Branch("LaserSignal",LaserSignal,"LaserSignal[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserTime",LaserTime,"LaserTime[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserHHTime",LaserHHTime,"LaserHHTime[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserChisq",LaserChisq,"LaserChisq[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserID",LaserID,"LaserID[LaserNumber]/S");//LaserNumber
  trOut->Branch("LaserNDF",LaserNDF,"LaserNDF[LaserNumber]/S");//LaserNumber  

  trOut->Branch("EtcNumber",&EtcNumber,"EtcNumber/I");
  trOut->Branch("EtcSignal",EtcSignal,"EtcSignal[EtcNumber]/D");//EtcNumber
  trOut->Branch("EtcTime",EtcTime,"EtcTime[EtcNumber]/D");//EtcNumber
  trOut->Branch("EtcHHTime",EtcHHTime,"EtcHHTime[EtcNumber]/D");//EtcNumber
  trOut->Branch("EtcChisq",EtcChisq,"EtcChisq[EtcNumber]/D");//EtcNumber
  trOut->Branch("EtcID",EtcID,"EtcID[EtcNumber]/S");//EtcNumber
  trOut->Branch("EtcNDF",EtcNDF,"EtcNDF[EtcNumber]/S");//EtcNumber  


  std::cout<< tr->GetEntries() << std::endl;
  for( int i = 0; i< tr->GetEntries(); i++){
    tr->GetEntry(i);
    
    //if( LaserNumber > 0 && LaserSignal[0] > 200 ){ continue; }// Except Laser Event
    if( LaserNumber !=0 && LaserSignal[0] > 200 && LaserID[0] == 0){ continue; }
    cosmicTrig->TriggerDecision();
    trOut->Fill();
  }
  tfOut->Close();
}
