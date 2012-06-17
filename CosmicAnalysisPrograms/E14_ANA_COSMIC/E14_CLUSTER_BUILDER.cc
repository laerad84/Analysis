#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"

#include "E14ReadSumFile.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "HoughCsI.h"
#include "TThread.h"

TCanvas*   can;
IDHandler* handler;
CsIImage*  image;
TThread*   th;
Bool_t     finished;
void* handle(void *ptr){
  
  can->cd();
  image->Draw();
  can->Modified();
  can->Update();
  
}

int
main(int argc, char** argv){
  finished = kFALSE;
  gStyle->SetPalette(1);
  // ARGV[0]  <InputROOTFile> <OutputROOTFile> [<Vision>]
  std::string InputFile;
  std::string OutputFile;
  
  if( argc == 3){
    InputFile = argv[1];
    OutputFile = argv[2];
  }else{
    std::cerr << "ARGUEMENT ERROR" << std::endl;
    return -1;
  }
  
  //

  TApplication* app = new TApplication("app",&argc, argv);
  E14ReadSumFile* Reader = new E14ReadSumFile();
  Reader->Add(InputFile.c_str());

  TFile* tfout = new TFile(OutputFile.c_str(),"RECREATE");
  TTree* trout = new TTree("tro","output from clustering");
  E14GNAnaDataContainer data;
  data.branchOfClusterList(trout);
  data.branchOfDigi(trout);
  
  handler = new IDHandler("Data/crystal.txt");
  image   = new CsIImage(handler);
  can     = new TCanvas("can","",800,0,800,800);

  
  int nCsIDigi = 0;
  int CsIDigiID[3000]   = {-1};
  int CsIDigiTime[3000] = {0};
  double CsIDigiE[3000] = {0};
  double CalFactor[3000] = {0};
  
  //Set Calibration Factor
  for( int i = 0; i< 2716; i++){
    CalFactor[i] = 14.0/1000;
  }
  
  ClusterFinder cFinder;
  
  long nEntries = Reader->GetEntries();
  std::cout << "TotalEntries::" << nEntries  << std::endl;
  for( int iEntry = 0; iEntry < nEntries; iEntry++){    
    Reader->GetEntry(iEntry);    
    //image->SetTitle(Form("Event:%d",iEntry));
    nCsIDigi = 0;
    for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
      if( Reader->CsiEne[idigi] > 100 ){
	CsIDigiID[nCsIDigi]   = Reader->CsiModID[idigi];
	CsIDigiE[nCsIDigi]    = Reader->CsiEne[idigi]*CalFactor[Reader->CsiModID[idigi]];
	CsIDigiTIme[nCsIDigi] = Reader->CsiTime[idigi]; 	
	nCsIDigi++;
      }
    }  
    std::list<Cluster> clist = clusterFinder.findCluster(nCsIDigi,CsIDigiID,CsIDigiE,CsIDigiTime);
    std::cout << list.size() << std::endl;
  }  
  app->Run();
}
