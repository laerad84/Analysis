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
  std::string iFileForm = "%s/Data_All_Cal_FNL_COS_HTADJ.root";//ROOTFILE_GAMMACLUS
  //std::string iFileForm = "%s/Data_All_old.root";//ROOTFILE_GAMMACLUS
  std::string oFileForm = "GammaClusterTime.root";
  
  TF1* TimeAdjFuncEnergy = new TF1("TimeAdjFuncEnergy",AdjFunc,0,2000,5);
  //TimeAdjFuncEnergy->SetParameters(-7.51860e-01,9.57348e-01,-9.55972e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-9.69218e-01,1.12202e+00,-7.52470e-02,0,0);
  TimeAdjFuncEnergy->SetParameters(-6.32987e-01,7.64508e-01,-9.74328e-02,0,0);//E:200-300;
  //TimeAdjFuncEnergy->SetParameters(-6.38508e-01,7.55330e-01,-9.05449e-02,0,0);//E:200-300,T:10-14
  //TimeAdjFuncEnergy->SetParameters(0,0,0,0,0);

  TF1* TimeAdjFuncHeight = new TF1("TimeAdjFuncHeight",AdjfuncTH,0,16000,3);
  TimeAdjFuncHeight->SetParameters(3.308,0.005644,0.0004385);
  
  TFile* tf       = new TFile(Form(iFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr       = (TTree*)tf->Get("trCluster");
  ClusterTimeReader* reader = new ClusterTimeReader(tr);
  Int_t  nEntries = reader->fChain->GetEntries();
  TFile* tfout    = new TFile(oFileForm.c_str(),"recreate");
  
  TH2D* hisRTNew = new TH2D("hisRTNew","hisRTNew",80,0,1000,100,-20,20);
  TH2D* hisRT    = new TH2D("hisRT"   ,"hisRT",80,0,1000,100,-50,50);
  TH2D* hisRTNon = new TH2D("hisRTNon","hisRTNon",80,0,1000,100,-20,20);
  TH1D* hisCut   = new TH1D("hisCut","hisCut",20,0,20);
  TH2D* hisWT    = new TH2D("hisWT","hisWT",80,-1000,1000,100,-20,20);
  TH2D* hisHT    = new TH2D("hisHT","hisHT",80,-1000,1000,100,-20,20);
  TProfile2D* profTime =new TProfile2D("profTime","profTime",40,-1000,1000,40,-1000,1000);
  /*
  double ene = gam.e();  
  double sinTheta = sin(gam.p3().theta());
  double newx = gam.coex()-L*sinTheta*cos(gam.p3().phi());
  double newy = gam.coey()-L*sinTheta*sin(gam.p3().phi());
  */
  static double const Pcor[2]={6.49003,0.99254};
  static double const CsIX0=18.5;//mm

  for( int ievent = 0; ievent < nEntries; ievent++){
    reader->GetEntry(ievent);
    hisCut->Fill(0);
    Int_t    MEGIndex    = -1;
    Double_t MEGEnergy   = 0;
    Double_t MEGTime     = 0;
    Double_t MEGTimeNon  = 0;
    Double_t MEGTimeNew  = 0;
    Double_t MEGR        = 1000;
    Double_t GTime[6]    = {0.};
    Double_t GRad[6]     = {0.};    
    Double_t GTimeNon[6] = {0.};
    Double_t GTimeNew[6] = {0.};
    Double_t GPhi[6]    = {0.};
    Double_t GE[6]      = {0.};
    Double_t GCE[6]     = {0.};
    Double_t GCT[6]     = {0.};

    if( reader->nCluster !=6 ){ continue; }
    bool ClusterRegion = true;
    bool ClusterMinRBool = false;//  < 300 mm at least one crystal 

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
      /*
      if( TMath::Abs(clusterR*cos(reader->ClusterPhi[clusterIndex] )) < 100 ||
	  TMath::Abs(clusterR*sin(reader->ClusterPhi[clusterIndex] )) < 100 ){
	ClusterRegion = false;
	continue;
      }
      */
      if( TMath::Abs(clusterR) >850 ){ 
	ClusterRegion = false;  
	continue; 
      }
      /*
      if( clusterR < 400 ){ ClusterMinRBool = true; }
      */

      //ShowerDepth
      double L           = CsIX0*(Pcor[0]+Pcor[1]*log(clusterEnergy/1000.));//mm
      double sinTheta    = sin(clusterTheta);
      double clusterRnew = clusterR-L*sinTheta;

      for( int crystalIndex = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	if( reader->ClusterEnergy[crystalIndex] > MECEnergy ){	  
	  MECEnergy = reader->CrystalEnergy[clusterIndex][crystalIndex];
	  GCE[clusterIndex] = MECEnergy;
	  MECHeight = reader->CrystalSignal[clusterIndex][crystalIndex];
	  GCT[clusterIndex] = clusterT + reader->CrystalT[clusterIndex][clusterIndex];
	}
      }

      Double_t clusteradjT   = clusterT-TimeAdjFuncHeight->Eval(MECHeight)+clusterR/sin(clusterTheta)/300.;
      Double_t clusternewT   = clusteradjT + L/299.8 - L*cos(clusterTheta)/80.;
      GTimeNon[clusterIndex] = clusterT;
      GTime[clusterIndex]    = clusteradjT;
      GTimeNew[clusterIndex] = clusternewT;
      GRad[clusterIndex]     = clusterR;
      GPhi[clusterIndex]    = clusterPhi;
      GE[clusterIndex]      = clusterEnergy;
      //if( clusterEnergy > MEGEnergy ){ 
      if( clusterR < MEGR ){ 
	MEGR      = clusterR;
	MEGEnergy = clusterEnergy;
	MEGTime   = clusteradjT;
	MEGTimeNon= clusterT;
	MEGTimeNew= clusternewT;
	MEGIndex  = clusterIndex;
      }
    }
    
    if( !ClusterRegion ){ continue; }
    hisCut->Fill(1);
    /*
    if( !ClusterMinRBool ){ continue; }
    hisCut->Fill(2);
    */

    for( int clusterIndex = 0; clusterIndex <reader->nCluster;clusterIndex++){
      if( clusterIndex == MEGIndex){ continue;}
      //if( abs(GRad[clusterIndex]*cos(GPhi[clusterIndex])) < 600 &&
      //abs(GRad[clusterIndex]*sin(GPhi[clusterIndex])) < 600 ){
      if( MEGEnergy < 100 || GE[clusterIndex] < 100 && GE[clusterIndex]>500){ continue; }
      hisRT->Fill(GRad[clusterIndex],GTime[clusterIndex]- MEGTime);
      hisRTNon->Fill(GRad[clusterIndex],GTimeNon[clusterIndex]- MEGTimeNon);
      hisRTNew->Fill(GRad[clusterIndex],GTimeNew[clusterIndex]- MEGTimeNew);
      profTime->Fill(GRad[clusterIndex]*cos(GPhi[clusterIndex]),
		     GRad[clusterIndex]*sin(GPhi[clusterIndex]),
		     MEGTime);
      hisWT->Fill(GRad[clusterIndex]*cos(GPhi[clusterIndex]),GTime[clusterIndex]-MEGTime);
      hisHT->Fill(GRad[clusterIndex]*sin(GPhi[clusterIndex]),GTime[clusterIndex]-MEGTime);
      //}
    }
  }
  hisWT->Write();
  hisHT->Write();
  profTime->Write();
  hisCut->Write();
  hisRTNon->Write();
  hisRTNew->Write();
  hisRT->Write();
  tfout->Close();
}
