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

#include <libgen.h>

int
main( int argc, char** argv ){
  Int_t RunNumber;
  std::string RunNumberList;
  if( argc != 2 ){
    return -1; 
  }else{
    RunNumberList = argv[1]; 
  }
  std::cout << RunNumberList << std::endl;
  
  std::string ROOTFILEWAV    = std::getenv("ROOTFILE_WAV");  
  std::cout<< __LINE__ << std::endl;
  std::string BaseFilename   = basename( argv[1] );
  std::cout<< __LINE__ << std::endl;
  std::string OutputFilename = Form("%s/HeightDistributionMerge_%s.root",
				    ROOTFILEWAV.c_str(),
				    BaseFilename.substr( 0, BaseFilename.size() -4 ).c_str() );
  
  std::cout << "BaseFilename  :" << BaseFilename << std::endl;
  std::cout << "OutputFilename:" << OutputFilename << std::endl; 

  std::ifstream ifsRunNumberList(RunNumberList.c_str());
  TFile* tfout  = new TFile(OutputFilename.c_str(),"RECREATE");

  TH1D*  hisHeightAll[2716]; 
  for( int i = 0; i< 2716; i++){
    hisHeightAll[i] = new TH1D(Form("hisHeightAll%d",i), Form( "HeightAll%d",i),1600, 0, 16000);
  }    
  
  std::cout<< "Start Loop" << std::endl; 
  while ( ifsRunNumberList >> RunNumber ){
    std::cout << RunNumber << std::endl;
    std::string InputFilename = Form("%s/HeightDistribution_%d.root",ROOTFILEWAV.c_str(), RunNumber ); 
    TFile* tfin = new TFile( InputFilename.c_str());    
    TH1D* hisHeight[2716];
    for( int ich = 0; ich < 2716; ich++){
      hisHeight[ich] = (TH1D*)tfin->Get(Form("hisHeight%d",ich));
      hisHeightAll[ich]->Add( hisHeight[ich] );
    }
    tfin->Close();
  }
  
  tfout->cd();
  for( int i = 0; i< 2716; i++){
    hisHeightAll[i]->Write() ;
  }  
  tfout->Close(); 
}

  
    

  
  
  
