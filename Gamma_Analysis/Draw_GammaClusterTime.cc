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
  std::string iFileForm = "%s/Data_All_Cal_FNL_COS.root";//ROOTFILE_GAMMACLUS
  std::string oFileForm = "GammaClusterTime.root";
  
  TF1* TimeAdjFuncEnergy = new TF1("TimeAdjFuncEnergy",AdjFunc,0,2000,5);
  //TimeAdjFuncEnergy->SetParameters(-7.51860e-01,9.57348e-01,-9.55972e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-9.69218e-01,1.12202e+00,-7.52470e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-6.32987e-01,7.64508e-01,-9.74328e-02,0,0);//E:200-300;
  //TimeAdjFuncEnergy->SetParameters(-6.38508e-01,7.55330e-01,-9.05449e-02,0,0);//E:200-300,T:10-14
  TimeAdjFuncEnergy->SetParameters(0,0,0,0,0);

  TF1* TimeAdjFuncHeight = new TF1("TimeAdjFuncHeight",AdjfuncTH,0,16000,3);
  TimeAdjFuncHeight->SetParameters(3.308,0.005644,0.0004385);
  
  TFile* tf       = new TFile(Form(iFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr       = (TTree*)tf->Get("trCluster");
  ClusterTimeReader* reader = new ClusterTimeReader(tr);
  Int_t  nEntries = reader->fChain->GetEntries();
  TFile* tfout    = new TFile(oFileForm.c_str(),"recreate");
  
  TH2D* hisRT    = new TH2D("hisRT"   ,"hisRT",80,0,1000,100,-10,10);
  TH2D* hisRTNon = new TH2D("hisRTNon","hisRTNon",80,0,1000,100,-10,10);
  
  for( int ievent = 0; ievent < nEntries; ievent++){
    reader->GetEntry(ievent);
    Int_t    MEGIndex    = -1;
    Double_t MEGEnergy   = 0;
    Double_t MEGTime     = 0;
    Double_t MEGTimeNon  = 0;
    Double_t GTime[6]= {0.};
    Double_t GRad[6] = {0.};    
    Double_t GTimeNon[6]= {0.};

    if( reader->nCluster !=6 ){ continue; }
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      Double_t clusterR      = reader->ClusterR[clusterIndex];
      Double_t clusterTheta  = reader->ClusterTheta[clusterIndex];
      Double_t clusterPhi    = reader->ClusterPhi[clusterIndex];
      Double_t clusterEnergy = reader->ClusterEnergy[clusterIndex];
      Double_t clusterT      = reader->ClusterT[clusterIndex];
      Double_t clusterZ      = clusterR/TMath::Tan(clusterTheta);
      Double_t clusterL      = clusterR/TMath::Sin(clusterTheta);
      Double_t MECHeight     = 0;
      Double_t MECEnergy     = 0;      
      for( int crystalIndex = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	if( reader->ClusterEnergy[crystalIndex] > MECEnergy ){
	  MECEnergy = reader->CrystalEnergy[clusterIndex][crystalIndex];
	  MECHeight = reader->CrystalSignal[clusterIndex][crystalIndex];
	}
      }

      Double_t clusteradjT   = clusterT-TimeAdjFuncHeight->Eval(MECHeight);
      GTimeNon[clusterIndex] = clusterT;
      GTime[clusterIndex] = clusteradjT;
      GRad[clusterIndex]  = clusterR;
      if( clusterEnergy > MEGEnergy ){ 
	MEGEnergy = clusterEnergy;
	MEGTime   = clusteradjT;
	MEGTimeNon= clusterT;
	MEGIndex  = clusterIndex;
      }
    }
    for( int clusterIndex = 0; clusterIndex <reader->nCluster;clusterIndex++){
      if( clusterIndex == MEGIndex){ continue;}
      hisRT->Fill(GRad[clusterIndex],GTime[clusterIndex]- MEGTime);
      hisRTNon->Fill(GRad[clusterIndex],GTimeNon[clusterIndex]- MEGTimeNon);
    }
  }
  hisRT->Write();
  tfout->Write();
}
