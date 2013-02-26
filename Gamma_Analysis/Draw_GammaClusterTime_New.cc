#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TApplication.h"
#include "ClusterTimeReader.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "CsIPoly.h"
#include "TH2Poly.h"


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
  double value = p0 + 1.2*p1*exp(p2*x0) + 0*p3*exp(p4*x0);
  return value;
}
double AdjFunc1( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double value = p0 + p1*exp(p2*x0);
  return value;
}
double AdjFuncTest( double* x, double* par){
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

int main( int argc, char** argv) {

  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string ROOTFILE_WAV       = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB        = std::getenv("ANALYSISLIB");
  
  TF1* TimeAdjFuncEnergy = new TF1("TimeAdjFuncEnergy",AdjFunc,0,2000,5);
  //TimeAdjFuncEnergy->SetParameters(-7.51860e-01,9.57348e-01,-9.55972e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-9.69218e-01,1.12202e+00,-7.52470e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-6.32987e-01,7.64508e-01,-9.74328e-02,0,0);//E:200-300;
  //TimeAdjFuncEnergy->SetParameters(-6.38508e-01,7.55330e-01,-9.05449e-02,0,0);//E:200-300,T:10-14
  //TimeAdjFuncEnergy->SetParameters(0,0,0,0,0);

  TF1* TimeAdjFuncHeight = new TF1("TimeAdjFuncHeight",AdjfuncTH,0,16000,3);
  TimeAdjFuncHeight->SetParameters(3.308,0.005644,0.0004385);
 
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


  TFile*  tf = new TFile(Form("%s/Data_All_GammaTime.root",ROOTFILE_GAMMACLUS.c_str()));
  TTree*  trIn = (TTree*)tf->Get("trCluster");
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

  trIn->SetBranchAddress("EventNumber"  ,&EventNumber );
  trIn->SetBranchAddress("nCluster"     ,&nCluster    );
  trIn->SetBranchAddress("nCrystal"     ,nCrystal     );
  trIn->SetBranchAddress("ClusterID"    ,ClusterID    );
  trIn->SetBranchAddress("ClusterEnergy",ClusterEnergy);
  trIn->SetBranchAddress("ClusterR"     ,ClusterR     );
  trIn->SetBranchAddress("ClusterT"     ,ClusterT     );
  trIn->SetBranchAddress("ClusterTheta" ,ClusterTheta );
  trIn->SetBranchAddress("ClusterChisq2",ClusterChisq2);
  trIn->SetBranchAddress("ClusterPhi"   ,ClusterPhi   );
  trIn->SetBranchAddress("CrystalT"     ,CrystalT     );
  trIn->SetBranchAddress("CrystalTAdj"  ,CrystalTAdj  );
  trIn->SetBranchAddress("CrystalEnergy",CrystalEnergy);
  trIn->SetBranchAddress("CrystalR"     ,CrystalR     );
  trIn->SetBranchAddress("CrystalPhi"   ,CrystalPhi   );
  trIn->SetBranchAddress("CrystalSignal",CrystalSignal);
  trIn->SetBranchAddress("CrystalID"    ,CrystalID    );
  trIn->SetBranchAddress("ClusterPos"   ,ClusterPos   );
  trIn->SetBranchAddress("GammaPos"     ,GammaPos     );
  trIn->SetBranchAddress("GammaTOF"     ,GammaTOF     );
  trIn->SetBranchAddress("GammaDOS"     ,GammaDOS     );
  trIn->SetBranchAddress("GammaRID"     ,GammaRID     );

  trIn->SetBranchAddress("KlongPos"     ,KlongPos     );
  trIn->SetBranchAddress("KlongEne"     ,&KlongEne    );
  trIn->SetBranchAddress("KlongMass"    ,&KlongMass   );
  trIn->SetBranchAddress("KlongPt"      ,&KlongPt     );
  trIn->SetBranchAddress("KlongChisqZ"  ,&KlongChisqZ );

  Int_t nEntries = trIn->GetEntries();
  
  TFile* tfout= new TFile("TimeOut.root","recreate");
  TProfile2D* profGammaEET  = new TProfile2D("profGammaEET","profGammaEET;GammaE[MeV];GammaE[MeV];MeanTimeDelta[ns]", 100,0,1000,100,0,1000);
  TH2D*       hisGammaTimeResolution[16];
  for( int i = 0; i < 10; i++){
    hisGammaTimeResolution[i] = new TH2D(Form("hisGammaTimeResolution_%d",i),Form("hisGammaTimeResolution_%d",i),
					 100,0,1000,400,-20,20);
  }
  double klmass = 497.648;
  for( int ievt = 0; ievt < nEntries; ievt++){
    trIn->GetEntry(ievt);
    if( KlongMass < klmass -10  || KlongMass > klmass +10 ){ continue; }
    for( int iCluster = 0;iCluster < nCluster-1; iCluster++){
      if( ClusterR[iCluster] < 250 ){ continue; }
      if( abs(GammaPos[iCluster][1])>550){ continue;}
      if( CrystalEnergy[iCluster][0] < 0.2*ClusterEnergy[iCluster] ){ continue;}
      for( int jCluster = iCluster+1; jCluster < nCluster; jCluster++){
	if( ClusterR[jCluster] < 250 ){ continue; }
	if( abs(GammaPos[jCluster][1])>550){ continue;}
	if( CrystalEnergy[jCluster][0] < 0.2*ClusterEnergy[jCluster] ){ continue; }
	Int_t histIndex  = -1;
	if( GammaRID[iCluster] == 1 ){
	  switch( GammaRID[jCluster] ){
	  case 1:
	    histIndex = 1;
	    break;
	  case 2:
	    histIndex = 2;
	    break;
	  case 3:
	    histIndex = 3;
	    break;
	  case 4:
	    histIndex = 4;
	    break;
	  default: 
	    histIndex  = -1;
	  }
	}else if( GammaRID[iCluster] == 2 ){
	  switch( GammaRID[jCluster] ){
	  case 1:
	    histIndex = 2;
	    break;
	  case 2:
	    histIndex = 5;
	    break;
	  case 3:
	    histIndex = 6;
	    break;
	  case 4:
	    histIndex = 7;
	    break;
	  default: 
	    histIndex  = -1;
	  }	  
	}else if( GammaRID[iCluster] == 3 ){
	  switch( GammaRID[jCluster] ){
	  case 1:
	    histIndex = 3;
	    break;
	  case 2:
	    histIndex = 6;
	    break;
	  case 3:
	    histIndex = 8;
	    break;
	  case 4:
	    histIndex = 9;
	    break;
	  default: 
	    histIndex  = -1;
	  }	  
	}else if( GammaRID[iCluster] == 4 ){
	  switch( GammaRID[jCluster] ){
	  case 1:
	    histIndex = 4;
	    break;
	  case 2:
	    histIndex = 7;
	    break;
	  case 3:
	    histIndex = 9;
	    break;
	  case 4:
	    histIndex = 10;
	    break;
	  default: 
	    histIndex  = -1;
	  }	  
	}else{
	  histIndex = -1;
	}

	if( ClusterEnergy[iCluster]*0.9 < ClusterEnergy[jCluster] &&
	    ClusterEnergy[iCluster]*1.1 > ClusterEnergy[jCluster] &&
	    histIndex > 0                                         ){
	  
	  hisGammaTimeResolution[histIndex-1]->Fill(ClusterEnergy[iCluster],
						    ClusterT[iCluster]-GammaTOF[iCluster]- (ClusterT[jCluster]-GammaTOF[jCluster]));
	}

	if( ClusterEnergy[iCluster] >= ClusterEnergy[jCluster] ){
	  profGammaEET->Fill( ClusterEnergy[iCluster],ClusterEnergy[jCluster],
			      ClusterT[iCluster]-GammaTOF[iCluster]- (ClusterT[jCluster]-GammaTOF[jCluster]));
	}else{
	  profGammaEET->Fill( ClusterEnergy[jCluster],ClusterEnergy[jCluster],
			      ClusterT[jCluster]-GammaTOF[jCluster]-(ClusterT[iCluster]-GammaTOF[iCluster]));
	}	    
      }
    }
  }  

  profGammaEET->Write();
  for( int i = 0; i< 10; i++){
    hisGammaTimeResolution[i]->Write();
  }
  tfout->Close();
}
