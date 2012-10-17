#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "E14WavReader.h"
#include <cstring>
#include <string>



int
main( int argc ,char ** argv ){
  
  int RunNumber = atoi( argv[1]);
  std::string WAVFILE = std::getenv( "ROOTFILE_WAV");
  TFile* tfin = new TFile(Form("%s/TEMPLATE_FIT_RESULT_1_%d.root",
			       WAVFILE.c_str(),RunNumber));
  /*
  TTree* trin = (TTree*)tfin->Get("WFTree");
  TFile* tfout = new TFile( Form("%s/TEMPLATE_FIT_RESULT_1_%d_TIME.root",WAVFILE.c_str(), RunNumber),"RECREATE");
  */
  TChain* trin = new TChain("WFTree");
  for( int i = 0; i< 12; i++){
    trin->Add(Form("%s/TEMPLATE_FIT_RESULT_1_%d.root",WAVFILE.c_str(),4503+i));
  }
  TFile* tfout = new TFile("Pi0Out.root","RECREATE");
  TH2D* hisTimeDelta = new TH2D("hisTimeDelta","hisTimeDelta",2716, 0, 2716,
				400, -100, 100 );  
  TH2D* hisEnergyTimeDelta[2716];
  for( int i = 0; i< 2716; i++){
    hisEnergyTimeDelta[i] = new TH2D(Form("hisEnergyTimeDelta%d",i),
				     Form("hisEnergyTimeDelta%d",i),
				     40,0,16000,
				     400,-100,100);
  }
    


  E14WavReader* reader = new E14WavReader(trin);
  Long_t entries =  reader->fChain->GetEntries();
  for( int ievent  = 0; ievent < entries ; ievent++){
    reader->GetEntry( ievent  );
    Double_t ScintiTime;
    Int_t    ScintiID;
    Bool_t   ScintiOn = false;
    for( int i = 0; i< reader->EtcNumber ; i++){    
      if( reader->EtcID[i] == 1 ){
	ScintiID   = reader->EtcID[i]; 
	ScintiTime = reader->EtcHHTime[i];
	ScintiOn = true;
	break;
      }
    }
    if( !ScintiOn ){ continue; }
    for( int ich  = 0; ich < reader->CsiNumber; ich++){
      
      if( reader->CsiSignal[ich]      < 50  ||
	  reader->EtcSignal[ScintiID] < 50  ){
	    continue;
      }
      int CsiID  = reader->CsiID[ich];
      double CsiTime = reader->CsiTime[ich];
      double CsiSignal = reader->CsiSignal[ich]; 
      hisTimeDelta->Fill( CsiID, CsiTime - ScintiTime );
      hisEnergyTimeDelta[CsiID]->Fill( CsiSignal, CsiTime- ScintiTime );
    }
  }
  for( int i = 0; i< 2716; i++){
    hisEnergyTimeDelta[i]->Write();
    }
  
  hisTimeDelta->Write();
  tfout->Close();
  return 0;
}
