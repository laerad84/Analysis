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
  char* name = argv[1];
  tf = new TFile(Form("%s",name));
  tr = (TTree*)tf->Get(Form("RecTree"));
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  
  TFile* tfOut = new TFile(Form("GammaTimingAna_%s",name),"recreate");


  TH1D* hisGammaTimeDelta = new TH1D("hisGammaTimeDelta","hisGammaTimeDelta",100,-10,10);
  TH2D* hisGammaDistTime  = new TH2D("hisGammaDistTime","hisGammaDistTime",1000,0,2000,100,-10,10);
  TH1D* hisGammaTimeChisq[2];
  for( int i = 0; i <2; i++){
    hisGammaTimeChisq[i] = new TH1D(Form("hisGammaTimeChisq_%d",i),Form("hisGammaTimeChisq_%d",i),200,0,50);
  }
  for( int ievent  =0; ievent < tr->GetEntries(); ievent++){
    tr->GetEntry( ievent );
    if( (ievent % 1000) == 0){
      std::cout<< ievent <<"/" << tr->GetEntries()<< std::endl;
    }
    std::list<Gamma>   glist;
    //std::vector<Klong> klVec;
    std::list<Pi0> plist;

    data.getData(glist);
    data.getData(plist);    

    int gIndex = 0;
    std::list<Gamma>::iterator git = glist.begin();
    std::list<Pi0>::iterator   pit = plist.begin();
    for( int i  =0; i< plist.size(); i++,pit++){
      Double_t Dist =TMath::Sqrt(TMath::Power((*pit).g1().x()-(*pit).g2().x(),2)
				 +TMath::Power((*pit).g1().y()-(*pit).g2().y(),2));

      Double_t FL[2];
      Double_t GammaTimeChisq[2];
      Double_t Weight[2] = {0,0};

      FL[0] = TMath::Sqrt( TMath::Power((*pit).g1().x(),2)+TMath::Power((*pit).g1().y(),2)+TMath::Power((*pit).g1().z()-(*pit).vz(),2));
      FL[1] = TMath::Sqrt( TMath::Power((*pit).g2().x(),2)+TMath::Power((*pit).g2().y(),2)+TMath::Power((*pit).g2().z()-(*pit).vz(),2));
      Double_t DeltaTime = (FL[1]-FL[0])/sol;

      for( int isize = 1; isize < (*pit).g1().clusterTimeVec().size(); isize++){
	GammaTimeChisq[0] += TMath::Power((*pit).g1().clusterTimeVec()[isize]-(*pit).g1().clusterTimeVec()[0],2)*(*pit).g1().e();
	Weight[0]    += (*pit).g1().e();
      }
      for( int isize = 1; isize < (*pit).g2().clusterTimeVec().size(); isize++){
	GammaTimeChisq[1] += TMath::Power((*pit).g2().clusterTimeVec()[isize]-(*pit).g2().clusterTimeVec()[0],2)*(*pit).g2().e();
	Weight[1]    += (*pit).g2().e();
      }
      GammaTimeChisq[0] = GammaTimeChisq[0]/Weight[0];
      GammaTimeChisq[1] = GammaTimeChisq[1]/Weight[1];

   
      if( abs((*pit).g2().clusterTimeVec()[0]-(*pit).g1().clusterTimeVec()[0]-DeltaTime) < 1 ){
	hisGammaTimeChisq[0]->Fill(GammaTimeChisq[0]);
	hisGammaTimeChisq[0]->Fill(GammaTimeChisq[1]);
      }else{
	hisGammaTimeChisq[1]->Fill(GammaTimeChisq[0]);
	hisGammaTimeChisq[1]->Fill(GammaTimeChisq[1]);
      }

      hisGammaDistTime->Fill(Dist,(*pit).g2().clusterTimeVec()[0]-(*pit).g1().clusterTimeVec()[0]-DeltaTime);
      hisGammaTimeDelta->Fill( (*pit).g2().clusterTimeVec()[0]-(*pit).g1().clusterTimeVec()[0]-DeltaTime);

    }
  }

  for( int i = 0; i <2; i++){
    hisGammaTimeChisq[i]->Write();
  }
  hisGammaTimeDelta->Write();
  hisGammaDistTime->Write();
  tfOut->Close();
}
