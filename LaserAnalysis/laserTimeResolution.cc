// EXE : laserTimeResolution [RunNumber]
// IN  : run_wav_[RunNumber].root
// OUT : laser_output_[RunNumber].root

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
#include <cstring>
#include <string>
#include <fstream>


int main( int argc, char** argv){

  int         RunNumber     = atoi(argv[1]);
  std::string ROOTFILE_WAV  = std::getenv("ROOTFILE_WAV");  
  std::string InputFileForm = "%s/run_wav_%d.root";//ROOTFILE_WAV, RunNumber
  std::string OutputFileForm= "%s/laser_output_%d.root";//ROOTFILE_WAV, RunNumber
  TFile* tfin = new TFile(Form(InputFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber));
  TTree* trin = (TTree*)tfin->Get("Tree");

  const Int_t nCsI = 2716;
  Int_t    CsiNumber;
  Double_t CsiSignal[nCsI];//CsiNumber
  Double_t CsiTime[nCsI];//CsiNumber
  Short_t    CsiID[nCsI];//CsiNumber
  Int_t    LaserNumber;
  Double_t LaserSignal[5];//LaserNumber
  Double_t LaserTime[5];//LaserNumber
  Short_t LaserID[5];//LaserNumber
  trin->SetBranchAddress("CsiNumber"  ,&CsiNumber);
  trin->SetBranchAddress("CsiSignal"  ,CsiSignal);
  trin->SetBranchAddress("CsiTime"    ,CsiTime);
  trin->SetBranchAddress("CsiID"      ,CsiID);
  trin->SetBranchAddress("LaserNumber",&LaserNumber);
  trin->SetBranchAddress("LaserSignal",LaserSignal);
  trin->SetBranchAddress("LaserTime"  ,LaserTime);
  trin->SetBranchAddress("LaserID"    ,LaserID);

  TFile* tfout= new TFile(Form(OutputFileForm.c_str(),ROOTFILE_WAV.c_str(),RunNumber),"recreate");
  TH1D* hisLaserTriggerSignal = new TH1D("hisLaserTriggerSignal","hisLaserTriggerSignal"
					 160,0,16000);
  TH1D* hisCsILaserTime[nCsI];
  TH1D* hisCsILaserOutput[nCsI];
  for( int i = 0; i < 2716; i++){
    hisCsILaserTime[i] = new TH1D(Form("hisCsILaserTime_%d",i),
			       Form("hisCsILaserTime_%d;DeltaTime[ns];Entries/0.2[ns]",i),
			       200,-20,20);
    hisCsILaserOutput[i]= new TH1D(Form("hisCsILaserOutput_%d",i),
				   Form("hisCsLaserOutput_%d;LaserOutput[cnt];Entries/100[cnt]",i),
				   160,0,16000);
  }

  Int_t nEntries = trin->GetEntries();
  for( int ievent = 0; ievent < nEntries; ievent++ ){
    trin->GetEntry(ievent);
    std::cout << LaserSignal[0] << std::endl;
    hisLaserTriggerSignal->Fill(LaserSignal[0]);
    if(LaserSignal[0]<100){ continue; }
    // TimeDelta from CsITime ID == 0
    //Time of ID == 0;
    Double_t TimeDeltaBase = CsiTime[0];
    if( CsiID[0] == 0 ){
      TimeDeltaBase = CsiTime[0];	  
    }else{ 
      continue;
    }
    
    for( int i = 0; i> CsiNumber; i++){
      Int_t tmpID = CsiID[i];
      hisCsILaserTime[tmpID]->Fill(CsiTime[i]-TimeDeltaBase);
      hisCsILaserOutput[tmpID]->Fill(CsiSignal[i]);      
    }
  }

  for( int i = 0; i< nCsI; i++){
    hisCsILaserTime[i]->Write();
    hisCsILaserOutput[i]->Write();
  }
  tfout->Close();
  return 0;
}
