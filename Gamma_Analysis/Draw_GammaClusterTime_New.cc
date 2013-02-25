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

  TChain* ch = new TChain("T");
  std::ifstream ifs(Form("%s/Data/RunList/KLRunList_2.txt",ANALYSISLIB.c_str()));
  int tmpRunNumber;
  while( ifs >> tmpRunNumber ){
    std::cout<<  tmpRunNumber  << std::endl;
    ch->Add(Form("%s/run_wav_%d_Cal_FNL_COS.root", ROOTFILE_WAV.c_str(),tmpRunNumber));
  }
  E14GNAnaDataContainer data;
  data.setBranchAddress(ch);



  GammaFinder gFinder;
  std::cout<< ch->GetEntries() << std::endl;





  TFile* tfOut = new TFile("GammaTimeOut.root","recreate");
  TTree* trOut = new TTree("trOut","");

  Int_t    GammaRegion[6];
  Int_t    GammaCrystalID[6];
  Double_t GammaEnergy[6];
  Double_t GammaTime[6];
  Double_t GammaTimeV0[6];
  Double_t GammaPos[6][3];
  Double_t GammaTheta[6];
  Double_t LengthOfFlight[6];
  Double_t DepthOfShower[6];

  Double_t KlongPos[3];
  Double_t KlongEne;
  Double_t KlongMass;
  Double_t KlongPt;
  Double_t KlongChisqZ;

}
