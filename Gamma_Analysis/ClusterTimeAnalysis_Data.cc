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
#include "TF1.h"

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
#include "EnergyConverter.h"

double AdjFunc( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double value = p0 + p1*exp(p2*x0) + p3*exp(p4*x0);
  return value;
}


int main( int argc, char** argv ){
  TF1* TimeAdjFunc = new TF1("TimeAdjFunc",AdjFunc, 0, 2000,5);
  Double_t Par[5] = {-0.0905327+0.0112319,1.54915,-0.114423,0.0758477,0.00487457};
  Double_t ParErrors[5] = {0.00245834,0.0263153,0.00188467,0.007594,4.81501e-05};
  TimeAdjFunc->SetParameters(Par);
  TimeAdjFunc->SetParErrors(ParErrors);
  

  std::string ROOTFILE_WAV       = std::getenv("ROOTFILE_WAV");  
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string ANALYSISLIB        = std::getenv("ANALYSISLIB");

  EnergyConverter* Econv = new EnergyConverter();
  Econv->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",ANALYSISLIB.c_str()));

  TChain* ch = new TChain("T");  
  std::ifstream ifs(Form("%s/Data/RunList/KLRunList_2.txt",ANALYSISLIB.c_str()));
  int tmpRunNumber;
  while( ifs >> tmpRunNumber ){
    std::cout << tmpRunNumber << std::endl;
    //ch->Add( Form("%s/run_wav_%d_Cal_CosmicTime.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
    ch->Add( Form("%s/run_wav_%d_Cal_FNL_COS.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
  }
  GammaFinder gFinder;
  TFile*  tfout = new TFile(Form("%s/Data_All_Cal_FNL_COS.root",ROOTFILE_GAMMACLUS.c_str()),"RECREATE");
  TTree*  trout = new TTree("trCluster","");
  const int arrSize = 120;
  Int_t EventNumber;
  Int_t nCluster;
  Int_t nCrystal[arrSize];//[nCluster]
  Int_t ClusterID[arrSize];//[nCluster]

  Double_t ClusterTheta[arrSize];//[nCluster]
  Double_t ClusterT[arrSize];//[nCluster] 
  Double_t ClusterR[arrSize];//[nCluster]
  Double_t ClusterEnergy[arrSize];//[nCluster]
  Double_t ClusterPhi[arrSize];//[nCluster]

  Double_t CrystalT[arrSize][arrSize]; 
  Double_t CrystalEnergy[arrSize][arrSize];
  Double_t CrystalR[arrSize][arrSize];
  Double_t CrystalPhi[arrSize][arrSize];
  Double_t CrystalSignal[arrSize][arrSize];
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
  trout->Branch("CrystalSignal",CrystalSignal,"CrystalSignal[nCluster][120]/D");

  E14GNAnaDataContainer data;
  ClusterTimeAnalyzer* clusterTimeAnalyzer = new ClusterTimeAnalyzer();
  
  data.setBranchAddress( ch );
  std::cout<< ch->GetEntries() << std::endl;
  for( Int_t eventIndex = 0; eventIndex < ch->GetEntries(); eventIndex++){
    //if( eventIndex > 1E5 ){ break ; }
    ch->GetEntry( eventIndex );
    EventNumber = eventIndex ;
    nCluster = 0;

    for( int iarr = 0; iarr < arrSize; iarr++){
      nCrystal[iarr] = 0;
      ClusterID[iarr] = -1;
      ClusterEnergy[iarr] = 0; 
      ClusterR[iarr] = 0;
      ClusterPhi[iarr] = 0;
      ClusterT[iarr] = 0;      
      ClusterTheta[iarr] = 0;
      for( int jarr  =0; jarr < arrSize; jarr++){
	CrystalR[iarr][jarr]      = 0;
	CrystalT[iarr][jarr]      = 0;
	CrystalEnergy[iarr][jarr] = 0;
	CrystalPhi[iarr][jarr]    = 0;
	CrystalSignal[iarr][jarr] = 0;
      }
    }    
    std::list<Cluster>     clist;
    std::list<ClusterTime> ctlist;
    std::list<Gamma>       glist;
    std::vector<Klong>     klVec;
    data.getData( clist );
    data.getData( klVec );
    gFinder.findGamma( clist, glist );
    if( clist.size() < 6 ){ continue; }
    if( glist.size() != 6 ){ continue; }

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
      ClusterEnergy[clIndex]= (*itCl).e();
      ClusterTheta[clIndex] = TMath::ATan2((*itlist).GetClusterR(),(6148. - klVec[0].vz() ));
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
