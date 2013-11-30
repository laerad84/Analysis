#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"

int main( int argc, char** argv ){
  
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  int RunNumber = std::atoi( argv[1] );

  TChain* ch = new TChain("Tree");
  ch->Add(Form("%s/run_wav_%d.root",ROOTFILE_WAV.c_str(),RunNumber));
  //ch->Add(Form("%s/run_wav_%d.root",ROOTFILE_WAV.c_str(),4748));
  
  Int_t CsiNumber;
  Short_t CsiID[2716];
  Double_t CsiTime[2716];
  Double_t CsiSignal[2716];
  
  Int_t LaserNumber;
  Short_t LaserID[5];
  Double_t LaserSignal[5];
  Double_t LaserTime[5];
  
  ch->SetBranchAddress("CsiNumber",&CsiNumber);
  ch->SetBranchAddress("CsiID",CsiID);//CsiNumber
  ch->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  ch->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  ch->SetBranchAddress("LaserNumber",&LaserNumber);
  ch->SetBranchAddress("LaserID",LaserID);//LaserNumber
  ch->SetBranchAddress("LaserTime",LaserTime);//LaserNumber
  ch->SetBranchAddress("LaserSignal",LaserSignal);//LaserNumber

  TFile* tfOut = new TFile(Form("LaserHeightTiming_%d.root",RunNumber),"RECREATE");

  TH1D* hisLaserHeight[2716];
  TH1D* hisLaserT0[2716];
  TH2D* hisLaserHT[2716];
  for( int i = 0; i < 2716; i++){
    hisLaserHeight[i] = new TH1D(Form("hisLaserHeight_%d",i),Form("hisLaserHeight_%d",i),160,0,16000);
    hisLaserT0[i] = new TH1D(Form("hisLaserT0_%d",i),Form("hisLaserT0_%d",i),300,-20,40);
    hisLaserHT[i] = new TH2D(Form("hisLaserHT_%d",i),Form("hisLaserHT_%d",i),80,0,16000,300,-20,40);
  }

  Int_t nEntries = ch->GetEntries();
  for( int ievt = 0; ievt < nEntries; ievt++){
    ch->GetEntry(ievt);
    for( int iCsi = 0; iCsi < CsiNumber; iCsi++){
      hisLaserHeight[CsiID[iCsi]]->Fill(CsiSignal[iCsi]);
    }
    if( LaserSignal[0] < 200 ){ continue; }
    if( LaserTime[0]   < 50  ){ continue; }
    if( LaserID[0]    != 0   ){ continue; }
    
    for( int iCsi = 0; iCsi < CsiNumber; iCsi++){
      hisLaserHT[CsiID[iCsi]]->Fill(CsiSignal[iCsi],CsiTime[iCsi] - LaserTime[0] );
      if( CsiSignal[CsiID[iCsi]] <400 && CsiSignal[iCsi] > 200 ){ 
	hisLaserT0[CsiID[iCsi]]->Fill( CsiTime[iCsi]-LaserTime[0] );
      }
    }
  }
  for( int i = 0; i< 2716; i++){
    hisLaserHeight[i]->Write();
    hisLaserHT[i]->Write();
    hisLaserT0[i]->Write();
  }
  
  tfOut->Close();
}
