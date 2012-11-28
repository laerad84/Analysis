#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <list>
#include <vector>

#include "EventTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TApplication.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TProfile.h"
#include "TRandom.h"

#include "CsIPoly.h"
#include "IDHandler.h"
#include "PulseGenerator.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include "gamma/Gamma.h"
#include "cluster/Cluster.h"

#include "gamma/GammaFinder.h"
#include "cluster/ClusterFinder.h"

#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
//#include "User_Function.h"
#include "ClusterTimeStructure.h"

int main( int argc, char** argv ){

  Int_t GammaEnergy = atoi( argv[1] );
  Int_t Degree      = atoi( argv[2] );
  Int_t DataIndex   = atoi( argv[3] );

  std::cout<< __LINE__ <<std::endl;

  std::string ROOTFILE_GAMMAHIT  = std::getenv("ROOTFILE_GAMMAHIT");
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMAHIT");

  TChain* ch = new TChain("T");  
  ch->Add(Form("%s/Cluster_%dMeV_%ddeg-1E5-%d.root",ROOTFILE_GAMMAHIT.c_str(),GammaEnergy, Degree, DataIndex ));
  std::cout<< "File is Opened" << std::endl;

  TFile*  tfout = new TFile(Form("%s/Cluster_Time_%dMeV_%ddeg-1E5-%d.root",ROOTFILE_GAMMACLUS.c_str(),GammaEnergy, Degree,DataIndex ),"RECREATE");
  TTree* trout = new TTree("trCluster","");
  const int arrSize = 120;
  Int_t EventNumber;
  Int_t nCluster;
  Int_t nCrystal[arrSize];//[nCluster]
  Int_t ClusterID[arrSize];//[nCluster]

  Double_t ClusterT[arrSize];//[nCluster] 
  Double_t ClusterR[arrSize];//[nCluster]
  Double_t ClusterEnergy[arrSize];//[nCluster]
  Double_t ClusterPhi[arrSize];//[nCluster]
  Double_t ClusterTheta[arrSize];//[nCluster]
  
  Double_t CrystalT[arrSize][arrSize]; 
  Double_t CrystalEnergy[arrSize][arrSize];
  Double_t CrystalR[arrSize][arrSize];
  Double_t CrystalPhi[arrSize][arrSize];

  trout->Branch("EventNumber"  ,&EventNumber ,"EventNumber/I");
  trout->Branch("nCluster"     ,&nCluster    ,"nCluster/I");
  trout->Branch("nCrystal"     ,nCrystal     ,"nCrystal[nCluster]/I");
  trout->Branch("ClusterID"    ,ClusterID    ,"ClusterID[nCluster]/I");
  trout->Branch("ClusterEnergy",ClusterEnergy,"ClusterEnergy[nCluster]/D");
  trout->Branch("ClusterR"     ,ClusterR     ,"ClusterR[nCluster]/D");
  trout->Branch("ClusterT"     ,ClusterT     ,"ClusterT[nCluster]/D");
  trout->Branch("ClusterTheta" ,ClusterTheta ,"ClusterTheta[nCluster]/D");
  trout->Branch("ClusterPhi"   ,ClusterPhi   ,"ClusterPhi[nCluster]/D");
  trout->Branch("CrystalT"     ,CrystalT     ,"CrystalT[nCluster][120]/D");
  trout->Branch("CrystalEnergy",CrystalEnergy,"CrystalEnergy[nCluster][120]/D");
  trout->Branch("CrystalR"     ,CrystalR     ,"CrystalR[nCluster][120]/D");
  trout->Branch("CrystalPhi"   ,CrystalPhi   ,"CrystalPhi[nCluster][120]/D");

  E14GNAnaDataContainer data;
  ClusterTimeAnalyzer* clusterTimeAnalyzer = new ClusterTimeAnalyzer();
  
  data.setBranchAddress( ch );
  std::cout<< ch->GetEntries() << std::endl;
  for( Int_t eventIndex = 0; eventIndex < ch->GetEntries(); eventIndex++){
    ch->GetEntry( eventIndex );
    EventNumber = eventIndex ;
    nCluster = 0;

    for( int iarr = 0; iarr < arrSize; iarr++){
      nCrystal[iarr] = 0;
      ClusterID[iarr] = -1;
      ClusterEnergy[iarr] = 0; 
      ClusterR[iarr] = 0;
      ClusterTheta[iarr] = 0;
      ClusterT[iarr] = 0;      
      for( int jarr  =0; jarr < arrSize; jarr++){
	CrystalR[iarr][jarr] = 0;
	CrystalT[iarr][jarr] = 0;
	CrystalEnergy[iarr][jarr] = 0;
	CrystalPhi[iarr][jarr] = 0;
      }
    }
    
    std::list<Cluster>     clist;
    std::list<ClusterTime> ctlist;
    data.getData( clist );
    //std::cout<< "Number of Cluster:" <<  clist.size() << std::endl;
    clusterTimeAnalyzer->Convert( clist ,ctlist);
    if(( clist.size() != ctlist.size() ) ||
       ( clist.size() == 0 )             ||
       ( ctlist.size() == 0)             ){ continue; }
    nCluster = ctlist.size();    

    std::list<ClusterTime>::iterator  itlist = ctlist.begin();
    std::list<Cluster>::iterator      itCl   = clist.begin();
    Int_t clIndex = 0;
    for( ;
	 itlist != ctlist.end();
	 itlist++ ,itCl++ ,clIndex++ ){
      //std::cout<< (*itList).GetClusterTime() << " : " << (*itList).GetClusterR() << std::endl;
      ClusterID[clIndex]    = (*itCl).id();
      ClusterT[clIndex]     = (*itlist).GetClusterTime();
      ClusterR[clIndex]     = (*itlist).GetClusterR();
      ClusterPhi[clIndex]   = (*itlist).GetClusterPhi();
      ClusterTheta[clIndex] = (double)Degree;
      ClusterEnergy[clIndex]= (*itCl).e();
      nCrystal[clIndex]     = (*itCl).clusterIdVec().size();
      for( Int_t cryIndex = 0; cryIndex < nCrystal[clIndex] ; cryIndex++ ){
	CrystalR[ clIndex ][ cryIndex ]     = (*itlist).clusterRVec()[cryIndex];
	CrystalPhi[ clIndex ][ cryIndex]    = (*itlist).clusterPhiVec()[cryIndex];
	CrystalT[ clIndex ][ cryIndex ]     = (*itlist).clusterTimeDeltaVec()[cryIndex];
	CrystalEnergy[ clIndex ][ cryIndex] = (*itCl).clusterEVec()[cryIndex];	
      }
    }
    trout->Fill();
  }
  
  trout->Write();
  tfout->Close();
}
