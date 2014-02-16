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
#inlcude "cluster/ClusterFinder.h"
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"
#include "csimap/CsiMap.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TF1.h"


#include "IDHandler.h"
#include "E14WavReader_V1.h"
#include "CosmicTriggerTree.h"
#include "User_Function.h"
#include "User_Functions.h"
#include "GeneralFunctions.h"

int
main( int argc, char** argv ){


  int DataType = atoi(argv[1]);
  //   Id =0 : Data , ID = 2 : MC KL, ID = 3 : MC Beam

  const int nFileType = 3;
  std::vector<int> RunList[nFileType];
  for( int i = 0; i< 100; i++){
    RunList[1].push_back( i );
  }  

  Int_t RunN[24]={4502,4503,4504,4505,4506,4507,4508,4509,4510,4511,4512,4513,
		  4514,4515,4516,4517,4518,4519,4520,4521,4522,4523,4524,4525};
  for( int i = 0; i<24; i++){
    RunList[0].push_back(RunN);
  }

  TChain* ch;
  switch( DataType ){
  case 0:
    break;
  case 1:    
    break;
  default:
    return -1;
  }
  char* eventFilename[nFileType] = {"%s/run_wav_%d.root","%s/Simpi0_1E6_KLBEAM_%d.root",""};
  

  return 0;
}
