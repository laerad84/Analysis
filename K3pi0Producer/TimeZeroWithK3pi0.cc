#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <string>
#include <list>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

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
#include "TChain.h"

#include "E14WavReader.h"

int
main( int argc ,char ** argv ){
  
  int RunNumber = atoi( argv[1]);
  std::string WAVFILE = std::getenv( "ROOTFILE_WAV");

  TChain* trin = new TChain("Tree"); 
  trin->Add(Form("%s/run_wav_%d.root",WAVFILE.c_str(),RunNumber));

  TFile* tfout = new TFile(Form("%s/run_wav_%d_cl.root",WAVFILE.c_str(),RunNumber),"RECREATE");
  TTree* trout = new TTree("trCalibration", "Output from Time zero" );
  E14GNAnaDataContainer data; 
  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  data.branchOfKlong( trout );

  E14WavReader* reader = new E14WavReader(trin);
  int nCsIDigi           = 0;
  int CsIID[2716]        = {-1};
  double CsIEnergy[2716] = {0};
  double CsITime[2716]   = {-1};
  double CsIHHTime[2716] = {-1};

  std::cout<< __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  Long_t entries =  reader->fChain->GetEntries();
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
      int CsiID        = reader->CsiID[ich];
      double CsiTime   = reader->CsiTime[ich];
      double CsiSignal = reader->CsiSignal[ich]; 
      double CsiEnergy = reader->CsiEnergy[ich];
    }
    
  }

  tfout->Close();
  return 0;
}
