#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "DataAnalysis/E14ConvReader.h"
#include "DataAnalysis/E14MultiConvReader.h"
#include "DataAnalysis/E14IDHandler.h"

#include "GeneralMacros.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "WaveformFitter.h"
#include "Environment.h"

const int nCrate = 11;

int  main(int argc,char** argv)
{
  
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);
  int rtn = PrintEnvironment();
  unsigned int rst = GetEnvironment();
  std::cout << rst << std::endl;

  TGraph* gr[E14_NCsI];

  for( int i = 0; i< E14_NCsI ; ++i){
    gr[i]= new TGraph();
    gr[i]->SetMarkerStyle(5);
    gr[i]->SetNameTitle(Form("CH%04d",i),Form("CH%04d",i));
  }

  int runNumber = 4502;
  E14IDHandler* handler = new E14IDHandler();  
  E14ConvReader* conv[nCrate];  
  long nentries;
  
  for( int icrate = 0 ;icrate < nCrate; ++icrate){    
    conv[icrate] = new E14ConvReader();
    conv[icrate]->SetBranchAddress();
    conv[icrate]->AddFile(Form("%s/crate%d/run%d_conv.root",convFileDir.c_str(),icrate,runNumber));
    nentries = conv[icrate]->GetChainEntries();
    std::cout<< nentries << std::endl;
    std::cout<< icrate << std::endl;    
  }

  DEBUGOUT();
  conv[3]->GetChainEntry(0);
  DEBUGOUT();
  return 0;

}
