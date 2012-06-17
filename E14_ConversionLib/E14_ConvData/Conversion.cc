#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"

#include "TApplication.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TCanvas.h"

#include "TStyle.h"
#include "TSystem.h"

#include "E14Calibrator.h"
#include "E14Fill.h"
#include "E14RawData.h"
#include "E14DataReadVME.h"
#include "E14Mapper.h"
#include "SemiOfflineHist.h"

int
main( int argc, char** argv ){
  std::string inputComp = argv[1];
  std::string runNum    = argv[2];
  std::string crateID   = argv[3];

  E14Fill a;
  a.SetInputCompName(inputComp);
  a.SetRunNumber(atoi(runNum.c_str()));
  a.SetCrateID(atoi(crateID.c_str()));
  a.Fill();
  
  return 0;
}




