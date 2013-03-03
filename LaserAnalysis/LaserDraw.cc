#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"



int main( int argc, char** argv ){
  Int_t  iRun  = atoi(argv[1]);
  TFile* tfOut = new TFile(Form("LaserTimeStability_%d.root",iRun),"RECREATE");
  
  TH1D* hisLaserHeight = new TH1D(Form("hisLaserHeight_%d",iRun),Form("hisLaserHeight_%d",iRun),
			    160,0,4800);
  TH1D* hisLaserTime = new TH1D(Form("hisLaserTime_%d",iRun),Form("hisLaserTime_%d",iRun),
			  50,0,500);
  TH1D* hisTimeDeltaLaser[2716];
  TH1D* hisTimeDelta[2716];
  TH1D* hisHeight[2716];
  for( int ich = 0; ich < 2716; ich++){
    hisTimeDeltaLaser[ich] = new TH1D(Form("hisTimeDeltaLaser_%d_%d",ich,iRun),
					    Form("hisTimeDeltaLaser_%d_%d",ich,iRun),
					    500,-50,50);
    hisTimeDelta[ich] = new TH1D(Form("hisTimeDelta_%d_%d",ich,iRun),
				       Form("hisTimeDelta_%d_%d",ich,iRun),
				       200,-20,20);
    hisHeight[ich] = new TH1D(Form("hisHeight_%d_%d",ich,iRun),
				    Form("hisHeight_%d_%d",ich,iRun),
				    160,0,16000);
  }
  std::cout<< "Making Hists" <<std::endl;

  std::cout<< iRun << std::endl;
  std::string  ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string  iFileName    = "%s/laser_output_%d.root";
  TFile* tf = new TFile(Form(iFileName.c_str(),ROOTFILE_WAV.c_str(),iRun));
  TTree* tr = (TTree*)tf->Get("Tree");
  
  Int_t    CsiNumber;
  Short_t  CsiID[2716];
  Double_t CsiSignal[2716];
  Double_t CsiTime[2716];
  Int_t    LaserNumber;
  Short_t  LaserID[5];
  Double_t LaserSignal[5];
  Double_t LaserTime[5];
  
  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiID",CsiID);//CsiNumber
  tr->SetBranchAddress("CsISignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  tr->SetBranchAddress("LaserNumber",&LaserNumber);
  tr->SetBranchAddress("LaserID",LaserID);//LaserNumber
  tr->SetBranchAddress("LaserSignal",LaserSignal);//LaserNumber
  tr->SetBranchAddress("LaserTime",LaserTime);//LaserNumber
  
  Int_t nEntries = tr->GetEntries();
  for( int ievent  = 0; ievent < nEntries; ievent++){
    tr->GetEntry( ievent );
    if( LaserNumber < 5 ){ continue; }
    if( LaserTime[0] < 100 ) {continue; }
    if( LaserSignal[0] < 200 ){ continue; }
    if( CsiID[0] != 0 ){ continue; }
    
    hisLaserTime->Fill(LaserTime[0]);
    hisLaserHeight->Fill(LaserSignal[0]);
    for( int ich = 0; ich < CsiNumber; ich++){
      short tmpCsIID = CsiID[ich];
      hisTimeDeltaLaser[tmpCsIID]->Fill( CsiTime[ich] - LaserTime[0] );
      hisTimeDelta[tmpCsIID]->Fill( CsiTime[ich]- CsiTime[0]);
      hisHeight[tmpCsIID]->Fill(CsiSignal[ich]);
    }
  }
  
  tfOut->cd();

  for( int ich = 0; ich < 2716; ich++){
      hisTimeDeltaLaser[ich]->Write();
      hisTimeDelta[ich]->Write();
      hisHeight[ich]->Write();
  }
}


