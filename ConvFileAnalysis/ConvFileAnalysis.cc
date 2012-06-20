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

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "Structs.h"
#include "Environment.h"

#include "ConvFileAnalysis/E14ConvReader.h"
#include "ConvFileAnalysis/E14IDHandler.h"


const int nCrate = 11;

int  main(int argc,char** argv)
{
  if( argc !=2 ){
    std::cerr << "Please Input RunNumber " << std::endl;
    return -1; 
  }
  Int_t RunNumber = atoi( argv[1] );
  // GetEnvironment // 

  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  std::cout << ANALIBDIR << std::endl;
  int envRtn = GetEnvironment();
  PrintEnvironment();

  // Read Channel map  // 

  E14MapReader* map = new E14MapReader(Form("%s/Sum%d.root",sumFileDir,RunNumber));
  struct MapStruct CsiMap;
  struct MapStruct CC03Map;
  struct MapStruct OEVMap;
  struct MapStruct CosmicMap;
  struct MapStruct LaserMap;
  struct MapStruct CVMap;
  struct MapStruct EtcMap;
  
  {
    map->Add("Csi");
    map->Add("CC03");
    map->Add("OEV");
    map->Add("Cosmic");
    map->Add("Laser");
    map->Add("CV");
    map->Add("Etc");
    map->SetMap();
    std::cout<< "Setup Map is done." << std::endl; 
    std::cout<< map->GetNmodule() << std::endl;    
    
    map->CopyMap(0, CsiMap);
    map->CopyMap(1, CC03Map);
    map->CopyMap(2, OEVMap);
    map->CopyMap(3, CosmicMap);
    map->CopyMap(4, LaserMap);
    map->CopyMap(5, CVMap);
    map->CopyMap(6, EtcMap);
  }

  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  
  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    conv[icrate] = new E14ConvReader();
    conv[icrate]->AddFile( Form("%s/crate%d/run%d_conv.root",convFileDir.c_str(), icrate, RunNumber); 
  }

  TFile* tfout = new TFile(Form("run%d_wav.root",RunNumber),"recreate");
  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree"); 
  Int_t    nFit;
  Int_t    const nCsI = 2716;
  Int_t    ID[E14_NCsI];  
  Double_t wfPed[E14_NCsI];
  Double_t wfSignal[E14_NCsI];
  Double_t wfMean[E14_NCsI];
  Double_t wfParameterA[E14_NCsI];
  Double_t wfParameterB[E14_NCsI];
  Double_t wfHalfOfPeak[E14_NCsI];
  Int_t    L1BoardID[E14_NCsI];
  Short_t  TrigNoConv[E14_NCsI];
  Short_t  SpillNoConv[E14_NCsI];
  Short_t  DataArr[E14_NCsI][48];
  Float_t  IntFADC[E14_NCsI];
  Double_t Ped[E14_NCsI];
  Double_t Height[E14_NCsI];
  Double_t Mean[E14_NCsI];

  trout->Branch("nFit"        ,&nFit       ,"nFit/I");
  trout->Branch("ID"          ,ID          ,"ID[nFit]/I");//nFit
  trout->Branch("wfPed"       ,wfPed       ,"wfPed[nFit]/D");//nFit
  trout->Branch("wfSignal"    ,wfSignal    ,"wfSignal[nFit]/D");//nFit
  trout->Branch("wfMean"      ,wfMean      ,"wfMean[nFit]/D");//nFit
  trout->Branch("wfParameterA",wfParameterA,"wfParameterA[nFit]/D");//nFit
  trout->Branch("wfParameterB",wfParameterB,"wfParameterB[nFit]/D");//nFit
  trout->Branch("wfHalfOfPeak",wfHalfOfPeak,"wfHalfOfPeak[nFit]/D");//nFit
  //trout->Branch("nCsI"        ,&nCsI       ,"nCsI/I");
  trout->Branch("L1BoardID"   ,L1BoardID   ,"L1BoardID[nCsI]/I");//nCsI
  trout->Branch("IntFADC"     ,IntFADC     ,"IntFADC[nCsI]/F");//nCsI


  trout->Write();
  tfout->Close();
  app->Run();
  return 0;
}
