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

double AdjfuncTH( double*x ,double*par ){  
  double value = par[0]*(-1*TMath::Log(1 + par[1])+TMath::Log( 1 + par[1]*exp( par[2]*x[0] )));
  return value;
}


int main( int argc, char** argv ){
  TF1* TimeAdjFunc = new TF1("TimeAdjFunc",AdjFunc, 0, 2000,5);
  Double_t Par[5] = {-0.0905327+0.0112319,1.54915,-0.114423,0.0758477,0.00487457};
  Double_t ParErrors[5] = {0.00245834,0.0263153,0.00188467,0.007594,4.81501e-05};
  TimeAdjFunc->SetParameters(Par);
  TimeAdjFunc->SetParErrors(ParErrors);

  TF1* TimeAdjFuncHeight = new TF1("TimeAdjFuncHeight",AdjfuncTH,0,16000,3);
  TimeAdjFuncHeight->SetParameters(3.308,0.005644,0.0004385);


  std::string ROOTFILE_WAV       = std::getenv("ROOTFILE_WAV");  
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string ANALYSISLIB        = std::getenv("ANALYSISLIB");

  TChain* ch = new TChain("T");  
  std::ifstream ifs(Form("%s/Data/RunList/KLRunList_2.txt",ANALYSISLIB.c_str()));
  int tmpRunNumber;
  while( ifs >> tmpRunNumber ){
    std::cout << tmpRunNumber << std::endl;
    //ch->Add( Form("%s/run_wav_%d_Cal_CosmicTime.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
    ch->Add( Form("%s/run_wav_%d_Cal_FNL_COS.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
  }
  GammaFinder gFinder;

  TH2Poly* poly = new  TH2Poly("poly","",-1000,1000,-1000,1000);
  Double_t binX0[5] = {-1000,-600,-600,-1000,-1000};
  Double_t binY0[5] = {-600,-600,600,600,-600};
  Double_t binX1[7] = {75,-600,-600,-100,-100,75,75};
  Double_t binY1[7] = {-600,-600,600,600,0,0,-600};
  Double_t binX2[7] = {600,75,75,-100,-100,600,600};
  Double_t binY2[7] = {-600,-600,0,0,600,600,-600};
  Double_t binX3[5] = {600,1000,1000,600,600};
  Double_t binY3[5] = {-600,-600,600,600,-600};
  poly->AddBin( 5, binX0, binY0 );
  poly->AddBin( 7, binX1, binY1 );
  poly->AddBin( 7, binX2, binY2 );
  poly->AddBin( 5, binX3, binY3 );
  static double const Pcor[2]= {6.49003,0.99254};
  static double const CsIX0  = 18.5;//mm



  TFile*  tfout = new TFile(Form("%s/Data_All_GammaTime.root",ROOTFILE_GAMMACLUS.c_str()),"RECREATE");
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
  Double_t ClusterChisq2[arrSize];//[nCluster]
  Double_t ClusterPos[arrSize][3];//[nCluster]
  Double_t GammaPos[arrSize][3];//[nCluster]
  Double_t GammaTOF[arrSize];//[nCluster]
  Double_t GammaDOS[arrSize];//[nCluster]
  Int_t GammaRID[arrSize];//[nCluster]

  Double_t CrystalT[arrSize][arrSize]; 
  Double_t CrystalTAdj[arrSize][arrSize]; 
  Double_t CrystalEnergy[arrSize][arrSize];
  Double_t CrystalR[arrSize][arrSize];
  Double_t CrystalPhi[arrSize][arrSize];
  Double_t CrystalSignal[arrSize][arrSize];
  Double_t CrystalHHT[arrSize][arrSize];
  Int_t    CrystalID[arrSize][arrSize];

  Double_t KlongPos[3];
  Double_t KlongEne;
  Double_t KlongMass;
  Double_t KlongPt;
  Double_t KlongChisqZ;

  trout->Branch("EventNumber"  ,&EventNumber ,"EventNumber/I");
  trout->Branch("nCluster"     ,&nCluster    ,"nCluster/I");
  trout->Branch("nCrystal"     ,nCrystal     ,"nCrystal[nCluster]/I");
  trout->Branch("ClusterID"    ,ClusterID    ,"ClusterID[nCluster]/I");
  trout->Branch("ClusterEnergy",ClusterEnergy,"ClusterEnergy[nCluster]/D");
  trout->Branch("ClusterR"     ,ClusterR     ,"ClusterR[nCluster]/D");
  trout->Branch("ClusterT"     ,ClusterT     ,"ClusterT[nCluster]/D");
  trout->Branch("ClusterTheta" ,ClusterTheta ,"ClusterTheta[nCluster]/D");
  trout->Branch("ClusterChisq2",ClusterChisq2,"ClusterChisq2[nCluster]/D");
  trout->Branch("ClusterPhi"   ,ClusterPhi   ,"ClusterPhi[nCluster]/D");
  trout->Branch("CrystalT"     ,CrystalT     ,"CrystalT[nCluster][120]/D");
  trout->Branch("CrystalTAdj"  ,CrystalTAdj  ,"CrystalTAdj[nCluster][120]/D");
  trout->Branch("CrystalEnergy",CrystalEnergy,"CrystalEnergy[nCluster][120]/D");
  trout->Branch("CrystalR"     ,CrystalR     ,"CrystalR[nCluster][120]/D");
  trout->Branch("CrystalPhi"   ,CrystalPhi   ,"CrystalPhi[nCluster][120]/D");
  trout->Branch("CrystalSignal",CrystalSignal,"CrystalSignal[nCluster][120]/D");
  trout->Branch("CrystalID"    ,CrystalID    ,"CrystalID[nCluster][120]/I");
  trout->Branch("ClusterPos"   ,ClusterPos   ,"ClusterPos[nCluster][3]/D");
  trout->Branch("GammaPos"     ,GammaPos     ,"GammaPos[nCluster][3]/D");  
  trout->Branch("GammaTOF"     ,GammaTOF     ,"GammaTOF[nCluster]/D");
  trout->Branch("GammaDOS"     ,GammaDOS     ,"GammaDOS[nCluster]/D");
  trout->Branch("GammaRID"     ,GammaRID     ,"GammaRID[nCluster]/I");

  trout->Branch("KlongPos"     ,KlongPos     ,"KlongPos[3]/D");
  trout->Branch("KlongEne"     ,&KlongEne    ,"KlongEne/D");
  trout->Branch("KlongMass"    ,&KlongMass   ,"KlongMass/D");
  trout->Branch("KlongPt"      ,&KlongPt     ,"KlongPt/D");
  trout->Branch("KongChisqZ"   ,&KlongChisqZ ,"KlongChisqZ/D");

  E14GNAnaDataContainer data;
  ClusterTimeAnalyzer* clusterTimeAnalyzer = new ClusterTimeAnalyzer();
  
  data.setBranchAddress( ch );
  Int_t    CsiNumber;
  Int_t    CsiDigiID[2716];
  Double_t CsiTime[2716];
  Double_t CsiHHTime[2716];
  Double_t CsiSignal[2716];
  Double_t CsiEne[2716];
  ch->SetBranchAddress("CsiNumber",&CsiNumber);
  ch->SetBranchAddress("CsiModID" ,CsiDigiID);//CsiNumber
  ch->SetBranchAddress("CsiTime"  ,CsiTime   );//CsiNumber
  ch->SetBranchAddress("CsiHHTime",CsiHHTime );//CsiNumber
  ch->SetBranchAddress("CsiSignal",CsiSignal );//CsiNumber
  ch->SetBranchAddress("CsiEne"   ,CsiEne    );//CsiNumber
  //ch->SetBranchAddress

  std::cout<< ch->GetEntries() << std::endl;
  for( Int_t eventIndex = 0; eventIndex < ch->GetEntries(); eventIndex++){
  //for( Int_t eventIndex = 0; eventIndex < 1000; eventIndex++){
    //if( eventIndex > 1E5 ){ break ; }
    ch->GetEntry( eventIndex );
    EventNumber = eventIndex ;
    nCluster = 0;
    KlongPos[0] = 0;
    KlongPos[1] = 0;
    KlongPos[2] = 0;

    for( int iarr = 0; iarr < arrSize; iarr++){
      nCrystal[iarr]      = 0;
      ClusterID[iarr]     = -1;
      ClusterEnergy[iarr] = 0; 
      ClusterR[iarr]      = 0;
      ClusterPhi[iarr]    = 0;
      ClusterT[iarr]      = 0;      
      ClusterTheta[iarr]  = 0;
      ClusterChisq2[iarr] = 0;
      ClusterPos[iarr][0] = 0;
      ClusterPos[iarr][1] = 0;
      ClusterPos[iarr][2] = 0;
      GammaPos[iarr][0]   = 0;
      GammaPos[iarr][1]   = 0;
      GammaPos[iarr][2]   = 0;
      GammaTOF[iarr]      = 0.;
      GammaDOS[iarr]      = 0.;
      GammaRID[iarr]      = 0.;

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

    if( clist.size() != 6 ){ continue; }
    if( glist.size() != 6 ){ continue; }
    if( klVec.size() <= 0 ){ continue; }

    //std::cout<< "Number of Cluster:" <<  clist.size() << std::endl;

    clusterTimeAnalyzer->Convert( clist ,ctlist);    
    if(( clist.size() != ctlist.size() ) ||
       ( clist.size() == 0 )             ||
       ( ctlist.size() == 0)             ){ continue; }
    nCluster = ctlist.size();    

    std::list<ClusterTime>::iterator  itlist = ctlist.begin();
    std::list<Cluster>::iterator      itCl   = clist.begin();
    std::list<Gamma>::iterator        itgamma= glist.begin();
    KlongPos[0] = (klVec[0]).vx();
    KlongPos[1] = (klVec[0]).vy();
    KlongPos[2] = (klVec[0]).vz();
    KlongMass   = (klVec[0]).m();
    KlongChisqZ = (klVec[0]).chisqZ();

    Int_t clIndex = 0;
    for( ;
	 itlist != ctlist.end();
	 itlist++ ,itCl++ ,clIndex++, itgamma++){
      //std::cout<< (*itList).GetClusterTime() << " : " << (*itList).GetClusterR() << std::endl;

      ClusterID[clIndex]    = (*itCl).id();
      ClusterR[clIndex]     = (*itlist).GetClusterR();
      ClusterPhi[clIndex]   = (*itlist).GetClusterPhi();
      ClusterEnergy[clIndex]= (*itCl).e();
      ClusterTheta[clIndex] = TMath::ATan2((*itlist).GetClusterR(),(6148. - klVec[0].vz() ));
      ClusterPos[clIndex][0]= (*itCl).x();
      ClusterPos[clIndex][1]= (*itCl).y();
      ClusterPos[clIndex][2]= (*itCl).z();
      ClusterChisq2[clIndex]= (*itgamma).chisq();
      GammaPos[clIndex][0]  = (*itgamma).x();
      GammaPos[clIndex][1]  = (*itgamma).y();
      GammaPos[clIndex][2]  = (*itgamma).z();
      GammaTOF[clIndex]     = TMath::Sqrt((KlongPos[0]-GammaPos[clIndex][0])*(KlongPos[0]-GammaPos[clIndex][0])
					  +(KlongPos[1]-GammaPos[clIndex][1])*(KlongPos[1]-GammaPos[clIndex][1])
					  +(KlongPos[2]-GammaPos[clIndex][2])*(KlongPos[2]-GammaPos[clIndex][2]))/299.8;
      GammaDOS[clIndex]     = CsIX0*(Pcor[0]+Pcor[1]*log(ClusterEnergy[clIndex]/1000.));//mm

      ClusterT[clIndex]     = (*itlist).GetClusterTime();
      
      if( GammaPos[clIndex][0] < -600 ){
	GammaRID[clIndex] = 1;
      }else if( GammaPos[clIndex][0] > 600 ){
	GammaRID[clIndex] = 4;
      }else{
	if( GammaPos[clIndex][1] < 0 ){
	  if( GammaPos[clIndex][0] < 75 ){
	    GammaRID[clIndex] = 2;
	  }else{
	    GammaRID[clIndex] = 3;
	  }
	}else{
	  if( GammaPos[clIndex][0] < -100 ){
	    GammaRID[clIndex] = 2;
	  }else{ 
	    GammaRID[clIndex] = 3;
	  }
	}
      }
      
      nCrystal[clIndex]     = (*itCl).clusterIdVec().size();

      Int_t    MaxEID      = 0;
      Double_t MaxE        = 0;
      Double_t MaxEHeight  = 0; 
      Double_t MaxETime    = 0;
      for( Int_t cryIndex = 0; cryIndex < nCrystal[clIndex] ; cryIndex++ ){
	CrystalR[clIndex][cryIndex]      = (*itlist).clusterRVec()[cryIndex];
	CrystalPhi[clIndex][cryIndex]    = (*itlist).clusterPhiVec()[cryIndex];
	CrystalT[clIndex][cryIndex]      = (*itlist).clusterTimeDeltaVec()[cryIndex];

	CrystalEnergy[clIndex][cryIndex] = (*itCl).clusterEVec()[cryIndex];
	CrystalID[clIndex][cryIndex]     = (*itlist).clusterIDVec()[cryIndex];
	
	for( int idigi = 0; idigi < CsiNumber; idigi++){
	  if( CrystalID[clIndex][cryIndex]   == CsiDigiID[idigi] ){
	    CrystalSignal[clIndex][cryIndex] =  CsiSignal[idigi];
	    break;
	  }
	}

	CrystalTAdj[clIndex][cryIndex] = (*itlist).clusterTimeDeltaVec()[cryIndex] - TimeAdjFuncHeight->Eval(CrystalSignal[clIndex][cryIndex]);
	/// Set Maximum Energy deposited Crystal Time as ClusterTime;
	if( MaxE < CrystalEnergy[clIndex][cryIndex] ){
	  MaxE     = CrystalEnergy[clIndex][cryIndex];
	  MaxETime = (*itCl).clusterTimeVec()[cryIndex] -TimeAdjFuncHeight->Eval(CrystalSignal[clIndex][cryIndex]);
	}	
      }
      ClusterT[clIndex] = MaxETime;
    }
    trout->Fill();
  }

  trout->Write();
  tfout->Close();
}
