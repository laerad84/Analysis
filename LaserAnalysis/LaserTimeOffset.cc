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
  TFile* tf = new TFile(Form("%s/run_wav_Laser_4747.root",ROOTFILELASER.c_str()));
  TTree* tr = (TTree*)tf->Get("LaserTree");
  Int_t     CsiNumber;
  Short_t   CsiID[2716];
  Double_t  CsiSignal[2716];
  Double_t  CsiTime[2716];
  Int_t     LaserNumber;
  Double_t  LaserSignal[5];
  Double_t  LaserTime[5];

  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiID"    ,CsiID);//CsiNumber
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiTime"  ,CsiTime);//CsiNumber
  tr->SetBranchAddress("LaserNumber",&LaserNumber);
  tr->SetBranchAddress("LaserSignal",LaserSignal);//LaserNumber
  tr->SetBranchAddress("LaserTime"  ,LaserTime);//LaserNumber

  TFile* tfOut = new TFile("LaserTimingOffset.root","recreate");
  TH2D* hisLaserTimingOffset2D = new TH2D("hisLaserTimingOffset2D","hisLaserTimingOffset2D",2716,0,2716,500,0,50);
  
  TH1D* hisLaserTimingOffset[2716];
  for( int i = 0; i< 2716; i++){
    hisLaserTimingOffset[i] = new TH1D(Form("hisLaserTimingOffset_%d",i),Form("hisLaserTimingOffset_%d",i),500,0,50);    
  }
  TTree* trLaserTimeOffset= new TTree("TreeTimeOffset","TimeOffset:Laser");
  Int_t ID[2716];
  Double_t Offset[2716];
  Double_t Error[2716];
  trLaserTimeOffset->Branch("ID",ID,"ID[2716]/I");
  trLaserTimeOffset->Branch("Offset",Offset,"Offset[2716]/D");
  trLaserTimeOffset->Branch("Error",Error,"Error[2716]/D");

  for( int ievt = 0; ievt < tr->GetEntries(); ievt++){
    tr->GetEntry( ievt );
    for( int i  =0; i< CsiNumber; i++){
      if( CsiSignal[i] > 1000 && CsiSignal[i] < 2000 ){
	hisLaserTimingOffset[CsiID[i]]->Fill(CsiTime[i] - LaserTime[0]);
      }
      hisLaserTimingOffset2D->Fill(CsiID[i], CsiTime[i]-LaserTime[0]);
    }
  }

  std::ofstream ofs("LaserTimeOffset.txt");
  for( int i = 0; i< 2716; i++){
    ID[i] = i;
    if( hisLaserTimingOffset[i]->GetEntries() < 100 ){
      Offset[i] = 0;
      Error[i]  = 100;
    }else{
      Offset[i] = hisLaserTimingOffset[i]->GetMean();
      Error[i]  = hisLaserTimingOffset[i]->GetMeanError();
    }    
    ofs << ID[i] << "\t" << Offset[i] << "\t" << Error[i] << std::endl;

  }
  ofs.close();

  trLaserTimeOffset->Fill();

  hisLaserTimingOffset2D->Write();
  for( int i = 0; i< 2716; i++){
    hisLaserTimingOffset[i]->Write();
  }


  trLaserTimeOffset->Write();
  tfOut->Close();
}
