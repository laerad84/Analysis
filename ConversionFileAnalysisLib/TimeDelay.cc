// Calculation timing difference of each channels 
// Written by jwlee 
// 2012 6 30 
// 
//


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

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "Structs.h"

#include "E14ConvReader.h"
#include "E14IDHandler.h"
#include "E14ConvWriter.h"


#include "TROOT.h"
#include "Math/WrappedFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMethods.h"
#include "Math/WrappedTF1.h"

int main( int argc , char** argv ){
  
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  std::string SUMFILEDIR= std::getenv("ROOTFILE_SUMUP");
  std::cout << ANALIBDIR  << std::endl;
  std::cout << SUMFILEDIR << std::endl; 

  if( argc != 2 ) {
    std::cerr <<" Please Input RUnNumber "  << std::endl;
    return -1;    
  }
  
  Int_t RunNumber  = atoi( argv[1]) ;    
  TFile* tfin = new TFile( Form("run%d_wav.root", RunNumber));
  TTree* trin = (TTree*)tfin->Get("WFTree"); 
  if( tfin ==NULL || trin ==NULL ) { 
    std::cerr << "TFile of Tree is NULL " << std::endl;
    return -1; 
  }

  E14ConvWriter* wConv= new E14ConvWriter( Form("%s/Sum%d.root", SUMFILEDIR.c_str(),RunNumber), trin );
  wConv->AddModule("Csi");
  wConv->AddModule("CC03");
  wConv->AddModule("OEV");
  wConv->AddModule("CV");
  wConv->AddModule("Cosmic");
  wConv->AddModule("Laser");
  wConv->AddModule("Etc");
  
  wConv->Set();
  wConv->SetMap();
  wConv->SetBranchAddress();

  
  TFile* tfout = new TFile(Form("TimeDelta_%d.root", RunNumber),"recreate");
  TTree* trout = new TTree("TimeDelta", "");
  const int nCH = 2716; 
  double deltaTimeFit[nCH][nCH];
  double deltaTimeSpl[nCH][nCH];
  double deltaTimeHH[nCH][nCH];
  

  std::cout << "Check " << std::endl;
  std::cout << "Entries" << trin->GetEntries() << std::endl;
  
  std::cout << "Loop Start" << std::endl; 
  for( int ievent = 0; ievent < trin->GetEntries() ; ievent++ ) {
    trin->GetEntry(ievent);
    
    




  }
  
}
