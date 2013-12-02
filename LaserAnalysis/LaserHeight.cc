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
  
  Int_t    CsiNumber;
  Short_t  CsiID[2716];
  Double_t CsiTime[2716];
  Double_t CsiSignal[2716];
  
  Int_t    LaserNumber;
  Short_t  LaserID[5];
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

  TFile* tfOut = new TFile(Form("Data/LaserHeightTiming_%d.root",RunNumber),"RECREATE");
  TTree* trOut = new TTree("trLaser","trLaser");
  //Int_t RunNumber;
  Double_t Entries[2716];
  Double_t ID[2716];
  Double_t Output[2716];
  Double_t Mean[2716];
  Double_t RMS[2716];
  Double_t Error[2716];
  Double_t Chisq[2716];
  
  for( int i = 0; i< 2716; i++){
    ID[i] = -1;
    Entries[i] = 0;
    Output[i] = 0;
    RMS[i] = 0;
    Error[i] = 0;
    Chisq[i] = 0;
    
  }
  
  trOut->Branch("RunNumber",&RunNumber,"RunNumber/I");
  trOut->Branch("Entries",Entries,"Entries[2716]/D");
  trOut->Branch("Output",Output,"Output[2716]/D");
  trOut->Branch("Error",Error,"Error[2716]/D");
  trOut->Branch("Chisq",Chisq,"Chisq[2716]/D");
  trOut->Branch("Mean",Mean,"Mean[2716]/D");
  trOut->Branch("RMS",RMS,"RMS[2716]/D");
  

  TH1D* hisLaserHeight[2716];
  TH1D* hisLaserT0[2716];
  TH2D* hisLaserHT[2716];
  for( int i = 0; i < 2716; i++){
    hisLaserHeight[i] = new TH1D(Form("hisLaserHeight_%d",i),Form("hisLaserHeight_%d",i),1600,0,16000);
    hisLaserT0[i] = new TH1D(Form("hisLaserT0_%d",i),Form("hisLaserT0_%d",i),300,-20,40);
    hisLaserHT[i] = new TH2D(Form("hisLaserHT_%d",i),Form("hisLaserHT_%d",i),80,0,16000,300,-20,40);
  }

  TH1D* hisLaserTime = new TH1D("hisLaserTime","hisLaserTime", 500,0,500);
  TH1D* hisLaserSignal = new TH1D("hisLaserSignal","hisLasersignal",160,0,16000);
  
  Int_t nEntries = ch->GetEntries();
  for( int ievt = 0; ievt < nEntries; ievt++){
  //for( int ievt = 0; ievt < 5000; ievt++){
    ch->GetEntry(ievt);
    if( LaserID[0]     != 0   ){ continue; }
    if( LaserTime[0]   <  140 ){ continue; }
    if( LaserSignal[0] <  100 ){ continue; }
    
    for( int iCsi = 0; iCsi < CsiNumber; iCsi++){
      if( CsiSignal[iCsi] < 10 ){ continue; }
      if( CsiTime[iCsi] - LaserTime[0]  < 0 ||
	  CsiTime[iCsi] - LaserTime[0] > 50 ){
	continue; 
      }
      hisLaserHeight[(int)(CsiID[iCsi])]->Fill(CsiSignal[iCsi]);
    }
    if( LaserSignal[0] < 200 ){ continue; }
    
    for( int iCsi = 0; iCsi < CsiNumber; iCsi++){
      hisLaserHT[CsiID[iCsi]]->Fill(CsiSignal[iCsi],CsiTime[iCsi] - LaserTime[0] );
      if( CsiSignal[CsiID[iCsi]] <400 && CsiSignal[iCsi] > 200 ){ 
	hisLaserT0[CsiID[iCsi]]->Fill( CsiTime[iCsi]-LaserTime[0] );
      }
    }
  }

  for( int i = 0; i<2716; i++){
    Entries[i] = (Double_t)hisLaserHeight[i]->Integral();
    if( Entries[i] < 1 ){ continue; }
    hisLaserHeight[i]->Fit("gaus","","",
			   hisLaserHeight[i]->GetMean()-hisLaserHeight[i]->GetRMS(),
			   hisLaserHeight[i]->GetMean()+hisLaserHeight[i]->GetRMS());

    Mean[i] = hisLaserHeight[i]->GetMean();
    RMS[i]  = hisLaserHeight[i]->GetRMS();
    TF1* func=NULL;
    func = hisLaserHeight[i]->GetFunction("gaus");
    if( func ==NULL ){ continue;}
    Chisq[i] = (func->GetChisquare()/func->GetNDF());
    Output[i] = func->GetParameter(1);
    Error[i]  = func->GetParError(1);
  }
  for( int i = 0; i< 2716; i++){
    hisLaserHeight[i]->Write();
    hisLaserHT[i]->Write();
    hisLaserT0[i]->Write();
  }
  trOut->Fill();
  trOut->Write();
  tfOut->Close();
}
