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
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"


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


#include "User_Function.h"

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
  GammaFinder gFinder;
  ClusterFinder clusterFinder;


  const int nFile = 1;
  TFile* tf;
  TTree* tr;
  //char* name = "SumUp";
  char* name = argv[1];
  int CalRes = std::atoi( argv[2]);

  tf = new TFile(Form("%s",name));
  tr = (TTree*)tf->Get(Form("RecTree"));
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  
  TFile* tfOut = new TFile(Form("ShapeChisq_%d_%s",CalRes,name),"recreate");
  TTree* trOut = new TTree("RecTree","");
  E14GNAnaDataContainer dataCopy;
  dataCopy.branchOfPi0List(trOut);

  
  Double_t CalConst[2716];
  std::ifstream ifs(Form("Data/GainRMS_%d.txt",CalRes));
  int tmpId;
  double tmpCal;
  while( ifs >> tmpId >> tmpCal){
    CalConst[tmpId] = tmpCal;
  }

  TH1D* hisShapeChisqGamma[2];
  for( int i = 0; i <2; i++){
    hisShapeChisqGamma[i] = new TH1D(Form("hisShapeChisqGamma_%d",i),Form("hisShapeChisqGamma_%d",i),50,0,2.5);
  }
  for( int ievent  =0; ievent < tr->GetEntries(); ievent++){
    tr->GetEntry( ievent );
    if( (ievent % 1000) == 0){
      std::cout<< ievent <<"/" << tr->GetEntries()<< std::endl;
    }
    dataCopy.reset();
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::list<Gamma>   glistCal;
    //std::vector<Klong> klVec;
    std::list<Pi0> plist;
    std::list<Pi0> plistCal;

    data.getData(glist);
    data.getData(plist);    

    Int_t    nID=0;
    Int_t ID[2716];
    Double_t Energy[2716];
    Double_t Time[2716];
    for( int i = 0; i< 2716; i++){
      ID[i] = -1;
      Energy[i] = 0;
      Time[i] = 0;
    }
    int gIndex = 0;
    std::list<Gamma>::iterator git = glist.begin();
    std::list<Pi0>::iterator   pit = plist.begin();
    for( int i  =0; i< plist.size(); i++,pit++){
      for( int j = 0; j < (*pit).g1().clusterIdVec().size(); j++){
	ID[nID] = (*pit).g1().clusterIdVec()[j];
	Energy[nID] = (*pit).g1().clusterEVec()[j]*CalConst[(*pit).g1().clusterIdVec()[j]];
	Time[nID]   = (*pit).g1().clusterTimeVec()[j];
	nID++;
      }
      for( int j = 0; j < (*pit).g2().clusterIdVec().size(); j++){
	ID[nID] = (*pit).g2().clusterIdVec()[j];
	Energy[nID] = (*pit).g2().clusterEVec()[j]*CalConst[(*pit).g2().clusterIdVec()[j]];
	Time[nID]   = (*pit).g2().clusterTimeVec()[j];
	nID++;
      }
    }
    clist = clusterFinder.findCluster( nID, ID, Energy, Time );
    gFinder.findGamma(clist,glistCal);
    if( clist.size() < 2 ){ continue; }
    if( glistCal.size() != 2 ){ continue; }
    if( user_rec( glistCal, plistCal)){
      user_cut( dataCopy, plistCal);
      dataCopy.setData(plistCal);
      std::list<Pi0>::iterator pit = plistCal.begin();
      if( TMath::Abs( (*pit).g1().clusterTimeVec()[0] -(*pit).g2().clusterTimeVec()[0] ) < 1 ){
	hisShapeChisqGamma[0]->Fill((*pit).g1().chisq());
	hisShapeChisqGamma[0]->Fill((*pit).g2().chisq());
      }else{
	hisShapeChisqGamma[1]->Fill((*pit).g1().chisq());
	hisShapeChisqGamma[1]->Fill((*pit).g2().chisq());
      }
      
      trOut->Fill();
      dataCopy.eventID++;
    }
  }
  
  for( int i = 0; i< 2; i++){
    hisShapeChisqGamma[i]->Write();
  }

  trOut->Write();
  tfOut->Close();
}
