#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"


int
main( int argc, char** argv ){
  Int_t RunNumber; 
  if( argc != 2 ){
    return -1; 
  }else{
    RunNumber = atoi( argv[1] );
  }
  std::string ROOTFILEWAV = std::getenv("ROOTFILE_WAV");
  std::string InputFilename = Form("%s/TEMPLATE_FIT_RESULT_1_%d.root",ROOTFILEWAV.c_str(), RunNumber);
  std::string OutputFilename = Form("%s/HeightDistribution_%d.root",ROOTFILEWAV.c_str(), RunNumber ); 
  
  TFile* tfin = new TFile( InputFilename.c_str());
  TTree* trin = (TTree*)tfin->Get("WFTree");
  
  Int_t    CsiNumber; 
  Short_t    CsiID[2716];
  Double_t CsiSignal[2716];
  trin->SetBranchAddress( "CsiNumber", &CsiNumber);
  trin->SetBranchAddress( "CsiID", CsiID );
  trin->SetBranchAddress( "CsiSignal",CsiSignal );
  
  TFile* tfout = new TFile( OutputFilename.c_str(),"RECREATE");
  TH1D* hisHeight[2716]; 
  for( int i = 0; i< 2716; i++){
    hisHeight[i] = new TH1D(Form("hisHeight%d",i), Form( "Height%d",i),1600, 0, 16000);
  }
  
  Int_t entries = trin->GetEntries();
  for( int ievent  = 0 ; ievent < entries; ievent++){
    trin->GetEntry( ievent );
    for( int icrystal = 0; icrystal < CsiNumber; icrystal++){      
      Int_t crystalID = CsiID[ icrystal ];
      hisHeight[ crystalID ]-> Fill( CsiSignal[ icrystal ] );
    }
  }
  
  for( int i = 0; i< 2716; i++){
    hisHeight[i]->Write() ;
  }  
  tfout->Close(); 
}

  
    

  
  
  
