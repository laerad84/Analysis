#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <list>
#include <vector>

#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TText.h"
#include "TF1.h"
#include "TStyle.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "E14_CLUSTER_BUILDER/E14ReadSumFile.h"


void SetEnvironment(){
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("neMRIuo");
}


int
main(int argc, char** argv){
  
  
  
}


