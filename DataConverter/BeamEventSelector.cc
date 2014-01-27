#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "csimap/CsiMap.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "L1TrigCounter.h"
#include "CrateIDHandler.h"
#include "CosmicTriggerTree.h"


int main( int argc, char** argv ){

  if( argc != 2 ) { 
    std::cout<< "Wrong argument" << std::endl;
    return -1; 
  }
  Int_t RunNumber = std::atoi( argv[1]);

  std::string ANALYSISLIB = std::getenv("ANALYSISLIB");
  std::string IFDIR = std::getenv("ROOTFILEWAV");
  std::string OFDIR = std::getenv("ROOTFILEBEAM");

  char* InputFilename = Form("%s/run_wav_%d.root",IFDIR.c_str(),RunNumber);
  char* OutputFilename= Form("%s/run_wav_Beam_%d.root",OFDIR.c_str(),RunNumber);
  //char* TimeOffsetFile= Form("%s/Data/TimeOffset/TimeOffset_with_cosmic.dat");
  char* TimeOffsetFile= Form("%s/Data/TimeOffset/testNewWORKCompileOffset.txt",ANALYSISLIB.c_str());

  CsiMap* map = CsiMap::getCsiMap();
  CrateIDHandler*       CIDHandler = new CrateIDHandler();
  L1TrigCounter*        l1         = new L1TrigCounter();
  CosmicTriggerTree*    cosmicTrig = new CosmicTriggerTree();
  E14GNAnaDataContainer data;
  ClusterFinder         cFinder;
  GammaFinder           gFinder;

  l1->ReadMapFile();
  l1->SetThreshold(1800);
  l1->Reset();
  double CSIL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
					1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};

  Double_t TimeOffsetCosmic[2716];
  for( int i = 0; i< 2716; i++){
    TimeOffsetCosmic[i] = 0;
  }
  std::ifstream ifs(TimeOffsetFile);
  Int_t tmpID;
  Double_t tmpOffset, tmpRMS;
  while( ifs >> tmpID >> tmpOffset >> tmpRMS ){
    if( tmpRMS > 100 ){ continue; }
    TimeOffsetCosmic[tmpID] = tmpOffset;
  }

  std::cout<< "TimeOffset" << std::endl;
  TFile* tf = new TFile(InputFilename);
  TTree* tr = (TTree*)tf->Get("Tree");

  Int_t EventNo;

  Int_t CsiNumber;
  Double_t CsiSignal[2716];
  Double_t CsiEne[2716];
  Double_t CsiTime[2716];
  Double_t CsiHHTime[2716];
  Double_t CsiChisq[2716];
  Short_t  CsiID[2716];
  Short_t  CsiNDF[2716];
  Short_t  CsiPosID[2716];
  Short_t  CsiGB[2716];
  Short_t  CsiCrate[2716];
  Short_t  CsiL1[2716];
  
  Int_t    CsiL1nTrig;
  Double_t CsiL1TrigCount[20];

  Int_t LaserNumber;
  Double_t LaserSignal[10];
  Double_t LaserHHTime[10];
  Double_t LaserTime[10];
  Double_t LaserChisq[10];
  Short_t  LaserID[10];
  Short_t  LaserNDF[10];
  
  Int_t OEVNumber;
  Double_t OEVSignal[50];
  Double_t OEVTime[50];
  Double_t OEVChisq[50];
  Short_t  OEVID[50];
  Short_t  OEVNDF[50];

  Int_t CC03Number;
  Double_t CC03Signal[50];
  Double_t CC03Time[50];
  Double_t CC03Chisq[50];
  Short_t  CC03ID[50];
  Short_t  CC03NDF[50];

  Int_t CosmicNumber;
  Double_t CosmicSignal[20];
  Double_t CosmicTime[20];
  Double_t CosmicChisq[20];
  Short_t  CosmicID[20];
  Short_t  CosmicNDF[20];
  
  Int_t EtcNumber;
  Double_t EtcSignal[10];
  Double_t EtcTime[10];
  Double_t EtcHHTime[10];
  Double_t EtcChisq[10];
  Short_t  EtcID[10];
  Short_t  EtcNDF[10];

  tr->SetBranchAddress("EventNo",&EventNo);

  tr->SetBranchAddress("CsiNumber",&CsiNumber);  
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  tr->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  tr->SetBranchAddress("CsiHHTime",CsiHHTime);//CsiNumber
  tr->SetBranchAddress("CsiChisq",CsiChisq);//CsiNumber
  tr->SetBranchAddress("CsiID",CsiID);//CsiNumber
  tr->SetBranchAddress("CsiNDF",CsiNDF);//CsiNumber

  tr->SetBranchAddress("CC03Number",&CC03Number);
  tr->SetBranchAddress("CC03Signal",CC03Signal);//CC03Number
  tr->SetBranchAddress("CC03Time",CC03Time);//CC03Number
  tr->SetBranchAddress("CC03Chisq",CC03Chisq);//CC03Number
  tr->SetBranchAddress("CsiID",CsiID);//CsiNumber
  tr->SetBranchAddress("CC03NDF",CC03NDF);//CC03Number
  tr->SetBranchAddress("CC03ID",CC03ID);//CC03Number

  tr->SetBranchAddress("OEVNumber",&OEVNumber);
  tr->SetBranchAddress("OEVSignal",OEVSignal);//OEVNumber
  tr->SetBranchAddress("OEVTime",OEVTime);//OEVNumber
  tr->SetBranchAddress("OEVChisq",OEVChisq);//OEVNumber
  tr->SetBranchAddress("OEVID",OEVID);//OEVNumber
  tr->SetBranchAddress("OEVNDF",OEVNDF);//OEVNumber
  tr->SetBranchAddress("OEVID",OEVID);//OEVNumber

  tr->SetBranchAddress("LaserNumber",&LaserNumber);
  tr->SetBranchAddress("LaserSignal",LaserSignal);//LaserNumber
  tr->SetBranchAddress("LaserHHTime",LaserHHTime);//LaserNumber
  tr->SetBranchAddress("LaserTime",LaserTime);//LaserNumber
  tr->SetBranchAddress("LaserChisq",LaserChisq);//LaserNumber
  tr->SetBranchAddress("LaserID",LaserID);//LaserNumber
  tr->SetBranchAddress("LaserNDF",LaserNDF);//LaserNumber

  tr->SetBranchAddress("EtcNumber",&EtcNumber);
  tr->SetBranchAddress("EtcSignal",EtcSignal);//EtcNumber
  tr->SetBranchAddress("EtcHHTime",EtcHHTime);//EtcNumber
  tr->SetBranchAddress("EtcTime",EtcTime);//EtcNumber
  tr->SetBranchAddress("EtcChisq",EtcChisq);//EtcNumber
  tr->SetBranchAddress("EtcID",EtcID);//EtcNumber
  tr->SetBranchAddress("EtcNDF",EtcNDF);//EtcNumber

  std::cout<< "Branch" << std::endl;
  TFile* tfOut = new TFile(OutputFilename,"recreate");
  TTree* trOut = new TTree("EventTree","Beam Data");

  trOut->Branch("RunNo",&RunNumber,"RunNumber/I");
  trOut->Branch("EventNo",&EventNo,"EventNo/I");
  data.branchOfClusterList( trOut );
  data.branchOfGammaList( trOut );
  trOut->Branch("CsiNumber",&CsiNumber,"CsiNumber/I");
  trOut->Branch("CsiSignal",CsiSignal,"CsiSignal[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiEne",CsiEne,"CsiEne[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiTime",CsiTime,"CsiTime[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiHHTime",CsiHHTime,"CsiHHTime[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiChisq",CsiChisq,"CsiChisq[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiNDF",CsiNDF,"CsiNDF[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiID",CsiID,"CsiID[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiPosID",CsiPosID,"CsiPosID[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiGB",CsiGB,"CsiGB[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiCrate",CsiCrate,"CsiCrate[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiL1",CsiL1,"CsiL1[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiL1nTrig",&CsiL1nTrig,"CsiL1nTrig/I");
  trOut->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsiL1TrigCount[20]/D");
  int s_arrSize = 120;
  Int_t    GamClusNumbers;
  Int_t    GamClusSizes[120];
  Double_t GamClusCsiSignal[120][120];
  Double_t GamClusCsiChisq[120][120];
  Int_t    GamClusCsiL1[120][120];
  Int_t    GamClusCsiCrate[120][120];

  trOut->Branch("GamClusNumbers",&GamClusNumbers,"GamClusNumbers/I");
  trOut->Branch("GamClusSizes",&GamClusSizes,"GamClusSizes[GamClusNumbers]/I");  
  trOut->Branch("GamClusCsiSignal",GamClusCsiSignal,Form("GamClusCsiSignal[GamClusNumbers][%d]/D",s_arrSize));
  trOut->Branch("GamClusCsiChisq",GamClusCsiChisq,Form("GamClusCsiChisq[GamClusNumbers][%d]/D",s_arrSize));
  trOut->Branch("GamClusCsiL1",GamClusCsiL1,Form("GamClusCsiL1[GamClusNumbers][%d]/I",s_arrSize));
  trOut->Branch("GamClusCsiCrate",GamClusCsiCrate,Form("GamClusCsiCrate[GamClusNumbers][%d]/I",s_arrSize));

  trOut->Branch("OEVNumber",&OEVNumber,"OEVNumber/I");
  trOut->Branch("OEVSignal",OEVSignal,"OEVSignal[OEVNumber]/D");//OEVNumber
  trOut->Branch("OEVTime",OEVTime,"OEVTime[OEVNumber]/D");//OEVNumber
  trOut->Branch("OEVChisq",OEVChisq,"OEVChisq[OEVNumber]/D");//OEVNumber
  trOut->Branch("OEVNDF",OEVNDF,"OEVNDF[OEVNumber]/S");//OEVNumber
  trOut->Branch("OEVID",OEVID,"OEVID[OEVNumber]/S");//OEVNumber

  trOut->Branch("CC03Number",&CC03Number,"CC03Number/I");
  trOut->Branch("CC03Signal",CC03Signal,"CC03Signal[CC03Number]/D");//CC03Number
  trOut->Branch("CC03Time",CC03Time,"CC03Time[CC03Number]/D");//CC03Number
  trOut->Branch("CC03Chisq",CC03Chisq,"CC03Chisq[CC03Number]/D");//CC03Number
  trOut->Branch("CC03NDF",CC03NDF,"CC03NDF[CC03Number]/S");//CC03Number
  trOut->Branch("CC03ID",CC03ID,"CC03ID[CC03Number]/S");//CC03Number


  /*
  trOut->Branch("CosmicNumber",&CosmicNumber,"CosmicNumber/I");
  trOut->Branch("CosmicSignal",CosmicSignal,"CosmicSignal[CosmicNumber]/D");//CosmicNumber
  trOut->Branch("CosmicTime",CosmicTime,"CosmicTime[CosmicNumber]/D");//CosmicNumber
  trOut->Branch("CosmicChisq",CosmicChisq,"CosmicChisq[CosmicNumber]/D");//CosmicNumber
  trOut->Branch("CosmicID",CosmicID,"CosmicID[CosmicNumber]/S");//CosmicNumber
  trOut->Branch("CosmicNDF",CosmicNDF,"CosmicNDF[CosmicNumber]/S");//CosmicNumber  
  */

  trOut->Branch("EtcNumber",&EtcNumber,"EtcNumber/I");
  trOut->Branch("EtcSignal",EtcSignal,"EtcSignal[EtcNumber]/D");//EtcNumber
  trOut->Branch("EtcTime",EtcTime,"EtcTime[EtcNumber]/D");//EtcNumber
  trOut->Branch("EtcHHTime",EtcHHTime,"EtcHHTime[EtcNumber]/D");//EtcNumber
  trOut->Branch("EtcChisq",EtcChisq,"EtcChisq[EtcNumber]/D");//EtcNumber
  trOut->Branch("EtcID",EtcID,"EtcID[EtcNumber]/S");//EtcNumber
  trOut->Branch("EtcNDF",EtcNDF,"EtcNDF[EtcNumber]/S");//EtcNumber  

  trOut->Branch("LaserNumber",&LaserNumber,"LaserNumber/I");
  trOut->Branch("LaserSignal",LaserSignal,"LaserSignal[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserTime",LaserTime,"LaserTime[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserHHTime",LaserHHTime,"LaserHHTime[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserChisq",LaserChisq,"LaserChisq[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserID",LaserID,"LaserID[LaserNumber]/S");//LaserNumber
  trOut->Branch("LaserNDF",LaserNDF,"LaserNDF[LaserNumber]/S");//LaserNumber  

  cosmicTrig->SetBranchAddress(tr);
  cosmicTrig->Branch(trOut);

  for( int i = 0; i< tr->GetEntries(); i++){
    tr->GetEntry(i);
    l1->Reset();

    if( LaserNumber != 0 && LaserSignal[0] > 200 && LaserID[0] == 0){ continue; }
    cosmicTrig->TriggerDecision();
    for( int i = 0; i< 2716; i++){
      CsiPosID[i] = -1;
    }
    Int_t CsiIntID[2716]={0};

    for( int i =0 ; i< CsiNumber; i++){
      l1->Fill(CsiID[i],CsiSignal[i]);
      Double_t x(0),y(0);
      x = map->getX(CsiID[i]);
      y = map->getY(CsiID[i]);
      if( x > 0 ){
	if( y > 0 ){ CsiPosID[i]  = 0; }
	else{ CsiPosID[i] = 1; }
      }else{
	if( y > 0 ){ CsiPosID[i] = 2; }
	else{ CsiPosID[i]  =3; }
      }
      if( y > 0 ){
	if( x < -100 ){
	  CsiGB[i] = -1;
	}else{
	  CsiGB[i] = 1;
	}
      }else{
	if( x < 75 ){
	  CsiGB[i] = -1;
	}else{
	  CsiGB[i] = 1;
	}
      }
      CsiL1nTrig = 0;
      std::vector<double> vecCount = l1->GetCount();
      if( vecCount.size() != 20){
	std::cout<< "VectorSizeError" << std::endl;
      }
      for( int ic = 0; ic< 20; ic++){
	CsiL1TrigCount[ic] = vecCount.at(ic);
	if( CsiL1TrigCount[ic] > CSIL1TrigCountThreshold[ic] ){
	  CsiL1nTrig++;
	}
      }

      CsiIntID[i] = CsiID[i];
      CsiCrate[i] = CIDHandler->GetCrate(CsiID[i]);
      CsiL1[i]    = CIDHandler->GetL1(CsiID[i]);
      CsiTime[i] = CsiTime[i] - TimeOffsetCosmic[CsiID[i]];
    }
    std::list<Cluster> clist;
    std::list<Gamma>   glist; 
    
    clist  = cFinder.findCluster( CsiNumber, CsiIntID, CsiEne, CsiTime );
    gFinder.findGamma( clist, glist );
    data.setData( clist );
    data.setData( glist );
    std::list<Gamma>::iterator git = glist.begin();
    GamClusNumbers = glist.size();
    int clNumber = 0;
    for(; git != glist.end(); git++){
      GamClusSizes[clNumber] = (*git).clusterIdVec().size();
      for( int i = 0; i< (*git).clusterIdVec().size(); i++){
	for( int in = 0; in < CsiNumber; in++){
	  if( CsiID[in] == (*git).clusterIdVec()[i]){
	    GamClusCsiSignal[clNumber][i] = CsiSignal[in];
	    GamClusCsiChisq[clNumber][i]  = CsiChisq[in];
	    GamClusCsiCrate[clNumber][i]  = CsiCrate[in];
	    GamClusCsiL1[clNumber][i]     = CsiL1[in];
	    break;
	  }
	}
      }
      clNumber++;
    }
    trOut->Fill();
  }
  tfOut->Close();
}
