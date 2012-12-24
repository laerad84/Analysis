#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"

#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "User_Function.h"


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

  TFile* tfout = new TFile("3Pi0Out.root","RECREATE");
  TTree* trout = new TTree("trCalibration", "Output from Time zero" );
  E14GNAnaDataContainer data; 
  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  data.branchOfKlong( trout );

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
  int nCsIDigi           = 0;
  int CsIID[2716]        = {-1};
  double CsIEnergy[2716] = {0};
  double CsITime[2716]   = {-1};
  double CsIHHTime[2716] = {-1};

  Long_t entries =  reader->fChain->GetEntries();
  
  /*
  for( int ievent  = 0; ievent < entries ; ievent++){
    nCsIDigi = 0;
    for( int ich = 0; ich < 3000; ich++){
      CsIID[ich]     = -1; 
      CsIEnergy[ich] = 0.;
      CsITime[ich]   = -1;
      CsIHHTime[ich] = -1;
    }
    
    reader->GetEntry( ievent  );
    for( int ich  = 0; ich < reader->CsiNumber; ich++){
      
      int CsiID  = reader->CsiID[ich];
      double CsiTime = reader->CsiTime[ich];
      double CsiSignal = reader->CsiSignal[ich]; 
    }
  }
  for( int i = 0; i< 2716; i++){
    hisEnergyTimeDelta[i]->Write();
  }
  */
  hisTimeDelta->Write();
  tfout->Close();
  return 0;
}
