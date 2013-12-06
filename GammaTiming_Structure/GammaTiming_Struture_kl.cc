#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"

#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
#include "TVector2.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/Gamma.h"
#include "gamma/GammaFinder.h"
#include "cluster/Cluster.h"
#include "cluster/ClusterFinder.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include <vector>
#include <list>
#include <string>
#include <cstdlib>
#include <cstdio>
#include "TDirectory.h"
#include "TProfile.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "IDHandler.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "csimap/CsiMap.h"


const double KLMass = 497.648;//MeV
double const sol = 299.792458;//[mm/nsec]
double const solc= 80;//[mm/nsec]
double const Pcor[2]={6.49003,0.99254};
double const CsIX0=18.5;//mm

double showerDepth(double x){
  double L = CsIX0*(Pcor[0]+Pcor[1]*log(x/1000.));//mm  
  return L;
}

double showerTimeDelay(Klong kl, Gamma g){
  double depth     = showerDepth(g.e());
  double cosTheta  = abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  double delayTime = 500/solc - depth*(1-cosTheta)/sol;
  return delayTime;
}

double showerTimeDelayAdj(Klong kl, Gamma g){
  double clusterTimingConstant = 0.541812;
  double depth     = showerDepth( g.e() );
  double cosTheta  = TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  //double delayTime = 500/solc - depth*(1-cosTheta)/sol + 0.541812*TMath::Log10(TMath::E())*log(g.e()/1000);
  double delayTime = depth*(1/sol-cosTheta/solc);
  return delayTime;
}

double KLSpectrum(double* x,double*p){
  //double klp   = TMath::Sqrt( x[0]*x[0]-KLMass*KLMass);  
  double klp = x[0];
  double sigma = p[2]*(1-(p[3]+p[4]*klp)*(klp-p[1]));
  double value = p[0]*TMath::Exp(-0.5*(klp-p[1])*(klp-p[1])/(sigma*sigma));
  return value;
}

double gammaLOF( Klong kl, Gamma g){
  double length = 0;
  length = sqrt( pow(g.x()-kl.vx(),2)+ pow(g.y()-kl.vy(),2)+pow(g.z()-kl.vz(),2));
  return length;
}
double HeightDelay( double height ){
  double value = TMath::Log( 1 + 0.03566*TMath::Exp( height/1621));
  return value;
}

//void DistributionTester(){
int main( int argc, char** argv){


  TF1* soltFunc  = new TF1("soltfunc",KLSpectrum,0,12,5);
  TF1* sugarFunc = new TF1("sugarfunc",KLSpectrum,0,12,5);
  const double soltPar[5] = { 1, 1.37,7.48e-1,-2.9e-1,1.68e-2};
  const double sugarPar[5] = {1,1.41991,0.810237,-0.301413,0.017092};
  soltFunc->SetParameters(soltPar);
  sugarFunc->SetParameters(sugarPar);
  IDHandler* handler = new IDHandler();
  CsiMap* map = CsiMap::getCsiMap();

  const int nFile = 1;
  TFile* tf;
  TTree* tr;
  //char* name = "SumUp";
  char* name = "KLCollection_Minimum";

  tf = new TFile(Form("%s.root",name));
  tr = (TTree*)tf->Get(Form("RecTree"));
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  Int_t CsiNumber;
  Int_t CsiModID[3000];
  Double_t CsiSignal[3000];
  Double_t CsiTime[3000];
  Double_t CsiEne[3000];
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  
  TFile* tfOut = new TFile(Form("GammaTimeShape_%s.root",name),"recreate");
  TTree* trOut = new TTree("GammaTimeShape","Time correction");

  Double_t Radius;//Distance from Csi Center
  Double_t ZVtx;//Distance from Csi surface
  Double_t theta;//InjectionAngle
  //Double_t cot;//cotangent of Injection Angle Z/R
  Double_t phi;//radial angle
  Double_t X;
  Double_t Y;
  Int_t    ClusterSize;
  Double_t E[120];//Energy;
  Double_t T[120];//Timing
  Double_t R[120];//Radial distance
  Double_t D[120];//
  Double_t FractionAngle[120];
  Int_t    CutCondition;
  Int_t    EventID;
  Double_t BaseTime;
  Double_t EGamma;
  Double_t ECenter;
  Double_t GammaChisq;
  trOut->Branch("EGamma",&EGamma,"EGamma/D");
  trOut->Branch("ECenter",&ECenter,"ECenter/D");
  trOut->Branch("EventID",&EventID,"EventID/I");
  trOut->Branch("Radius",&Radius,"Radius/D");
  trOut->Branch("X",&X,"X/D");
  trOut->Branch("Y",&Y,"Y/D");
  trOut->Branch("ZVtx",&ZVtx,"ZVtx/D");
  trOut->Branch("theta",&theta,"theta/D");
  trOut->Branch("phi",&phi,"phi/D");
  trOut->Branch("ClusterSize",&ClusterSize,"ClusterSize/I");
  trOut->Branch("BaseTime",&BaseTime,"BaseTime/D");
  trOut->Branch("E",E,"E[ClusterSize]/D");//ClusterSize
  trOut->Branch("T",T,"T[ClusterSize]/D");//ClusterSize
  trOut->Branch("R",R,"R[ClusterSize]/D");//ClusterSize
  trOut->Branch("D",D,"D[ClusterSize]/D");//ClusterSize
  trOut->Branch("FractionAngle",FractionAngle,"FractionAngle[ClusterSize]/D");//ClusterSize
  trOut->Branch("CutCondition",&CutCondition,"CutCondition/I");
  trOut->Branch("GammaChisq",&GammaChisq,"GammaChisq/D");

  TH1D* hisInjectionAngle = new TH1D("hisIndjectionAngle","hisInjectionAngle",100,0,100);
  int nTotal = tr->GetEntries(); 

  int cklmass = 0; 
  int cgammaPosIn =1;
  int cgammaPosOut=2;
  int cgammaPosCenter=3;
  int cgammaPosSmall=4;
  int cgammaPosLarge=5;
  
  for( int ievent  =0; ievent < nTotal; ievent++){
    tr->GetEntry( ievent );
    if( (ievent % 1000) == 0){
      std::cout<< ievent <<"/" << nTotal << std::endl;
    }
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    //std::list<Pi0> plist;

    data.getData(glist);
    //data.getData(plist);    
    data.getData(klVec);
    std::list<Gamma>::iterator git = glist.begin();
    //std::list<Pi0>::iterator   pit = plist.begin();

    CutCondition ==0;

    for( int i  =0; i< 6; i++, git++){
      double x = (*git).x();
      double y = (*git).y();
      if( TMath::Abs(x) < 150 && TMath::Abs(y) < 150 ){
	CutCondition |= 1 << cgammaPosIn;
      }
      if( TMath::Abs(y) > 550  || sqrt(x*x+y*y) > 850 ){
	CutCondition |= 1 << cgammaPosOut;
      }
    }

    for( int i = 0; i< klVec[0].pi0().size(); i++){
      if( TMath::Abs(klVec[0].pi0()[i].g1().x()) < 150 &&
	  TMath::Abs(klVec[0].pi0()[i].g1().y()) < 150 ){
	CutCondition |= 1 << cgammaPosIn;
      }
      if( TMath::Abs(klVec[0].pi0()[i].g2().x()) < 150 &&
	  TMath::Abs(klVec[0].pi0()[i].g2().y()) < 150 ){
	CutCondition |= 1 << cgammaPosIn;
      }
      if( TMath::Sqrt(klVec[0].pi0()[i].g1().x()*klVec[0].pi0()[i].g1().x()+klVec[0].pi0()[i].g1().y()*klVec[0].pi0()[i].g1().y()) > 850){
	CutCondition |= 1 << cgammaPosOut;
      }
      if( TMath::Sqrt(klVec[0].pi0()[i].g2().x()*klVec[0].pi0()[i].g2().x()+klVec[0].pi0()[i].g2().y()*klVec[0].pi0()[i].g2().y()) > 850){
	CutCondition |= 1 << cgammaPosOut;
      }      
    }
    

    int gIndex = 0;
    git = glist.begin();
    //pit = plist.begin();
    int cCutCondition = CutCondition;

    for( int i  =0; i< 3; i++,git++){
      ECenter = 0;
      EGamma = (*git).edep();
      GammaChisq=(*git).chisq();
      X  = (*git).coex();
      Y  = (*git).coey();
      ZVtx   = (*git).z() -klVec[0].pi0()[i].vz();
      CLHEP::Hep3Vector p=CLHEP::Hep3Vector(X,Y,ZVtx);
      Radius = p.perp();
      phi    = p.phi();
      theta  = p.theta();
      hisInjectionAngle->Fill(ZVtx/Radius);
      CutCondition = cCutCondition;
      BaseTime = 0;
      for( int j = 0; j< 120 ; j++){
	E[j] = 0;
	T[j] = 0;
	R[j] = 0;
	D[j] = 0;
	FractionAngle[j] = 0;
      }
      ClusterSize = (*git).clusterIdVec().size();
      for( int j  =0; j< (*git).clusterIdVec().size(); j++){	
	double posx = map->getX((*git).clusterIdVec()[j]);
	double posy = map->getY((*git).clusterIdVec()[j]);
	if( j == 0 ){ 
	  if( abs( posx - X ) > 12.5 ||  abs(posy - Y )> 12.5 ){
	    CutCondition |= 1 << cgammaPosCenter;
	  }
	}
	E[j]= (*git).clusterEVec()[j];
	T[j]= (*git).clusterTimeVec()[j];
	CLHEP::Hep3Vector v = CLHEP::Hep3Vector(posx-X,posy-Y,0);
	FractionAngle[j] = v.phi() - p.phi();
	R[j]=v.mag()*cos(FractionAngle[j]);
	D[j]=v.mag()*sin(FractionAngle[j]);
	if( v.mag() < 20 ){
	  BaseTime += E[j]*T[j];
	  ECenter += E[j];
	}
      }
      BaseTime=BaseTime/ECenter;
      trOut->Fill();
      //if(((i+1)%2) == 0 ){ pit++;}
    }
  }
  hisInjectionAngle->Write();
  trOut->Write();
  tfOut->Close();
}
