#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
int main( int argc, char** argv){
  std::string ROOTFILELASER=std::getenv("ROOTFILELASER");
  TFile* tfOffset =  new TFile("LaserTimingOffset.root");
  TTree* trLaserTimeOffset=(TTree*)tfOffset->Get("TreeTimeOffset");
  Int_t ID[2716];
  Double_t Offset[2716];
  Double_t Error[2716];
  trLaserTimeOffset->SetBranchAddress("ID",ID);
  trLaserTimeOffset->SetBranchAddress("Offset",Offset);
  trLaserTimeOffset->SetBranchAddress("Error",Error);
  trLaserTimeOffset->GetEntry(0);

  TFile* tf = new TFile(Form("%s/run_wav_Laser_4747.root",ROOTFILELASER.c_str()));
  TTree* tr = (TTree*)tf->Get("LaserTree");
  Int_t     CsiNumber;
  Short_t   CsiID[2716];
  Double_t  CsiSignal[2716];
  Double_t  CsiTime[2716];
  Int_t     LaserNumber;
  Short_t   LaserID[5];
  Double_t  LaserSignal[5];
  Double_t  LaserTime[5];

  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiID"    ,CsiID);//CsiNumber
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiTime"  ,CsiTime);//CsiNumber
  tr->SetBranchAddress("LaserNumber",&LaserNumber);
  tr->SetBranchAddress("LaserID",LaserID);//LaserNumber
  tr->SetBranchAddress("LaserSignal",LaserSignal);//LaserNumber
  tr->SetBranchAddress("LaserTime"  ,LaserTime);//LaserNumber

  TFile* tfLaserTimeDelay = new TFile("LaserTimeDelay.root","recreate");
  TH2D* hisLaserTimingDelay[2716];
  TTree* trOut = new TTree("LaserTree","Correct Offset");
  trOut->Branch("CsiNumber"  ,&CsiNumber,"CsiNumber/I");
  trOut->Branch("CsiID"      ,CsiID,"CsiID[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiSignal"  ,CsiSignal,"CsiSignal[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiTime"    ,CsiTime,"CsiTime[CsiNumber]/D");//CsiNumber
  trOut->Branch("LaserNumber",&LaserNumber,"LaserNumber/I");
  trOut->Branch("LaserID"    ,LaserID,"LaserID[LaserNumber]/S");//LaserNumber
  trOut->Branch("LaserSignal",LaserSignal,"LaserSignal[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserTime"  ,LaserTime,"LaserTime[LaserNumber]/D");//LaserNumber

  for( int i = 0; i< 2716; i++){
    hisLaserTimingDelay[i] = new TH2D(Form("hisLaserTimingDelay_%d",i),Form("hisLaserTimingDelay_%d",i),400,0,16000,500,-25,25);    
  }

  for( int ievt = 0; ievt < tr->GetEntries(); ievt++){
    tr->GetEntry( ievt );
    for( int i  =0; i< CsiNumber; i++){
      hisLaserTimingDelay[CsiID[i]]->Fill(CsiSignal[i],CsiTime[i] - LaserTime[0]-Offset[CsiID[i]]);
      CsiTime[i] =CsiTime[i] - LaserTime[0]-Offset[CsiID[i]]; 
    }
    trOut->Fill();
  }


  for( int i = 0; i< 2716; i++){
    hisLaserTimingDelay[i]->Write();
  }
  trOut->Write();
  tfLaserTimeDelay->Close();
}
