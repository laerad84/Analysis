#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector.h"
#include "TDirectory.h"
#include "TProfile.h"

#include "gnana/E14GNAnaDataContainer.h"
#include "gnana/E14GNAnaFunction.h"
#include "gamma/Gamma.h"
#include "gamma/GammaFinder.h"
#include "cluster/Cluster.h"
#include "cluster/ClusterFinder.h"
#include "pi0/Pi0.h"
#include "klong/Klong.h" 
#include "csimap/CsiMap.h"

#include <vector>
#include <list>
#include <cstdlib>
#include <cstdio>
#include <string>

#include "CLHEP/Vector/ThreeVector.h"

#include "LocalFunction.h"

const double KLMass = 497.648;//MeV
double const sol = 299.792458;//[mm/nsec]
double const solc= 80;//[mm/nsec]
double const Pcor[2]={6.49003,0.99254};
double const CsIX0=18.5;//mm

double gammaLOF( Klong kl, Gamma g){
  double length = 0;
  length = sqrt( pow(g.x()-kl.vx(),2)+ pow(g.y()-kl.vy(),2)+pow(g.z()-kl.vz(),2));
  return length;
}
double THFunction( double *x, double *p){
  double value=0;
  value = p[0]+p[1]*TMath::Log( 1 + p[2]*TMath::Exp(x[0]/2000));
  return value;
}

int main( int argc, char** argv){
  TFile* tf;
  TTree* tr;
  char* name = "WAVNOCV";
  tf = new TFile(Form("Kl_Total_%s.root",name));
  tr = (TTree*)tf->Get("trKL");
  
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  Int_t CsiNumber;
  Int_t CsiModID[3000]    = {-1};
  Double_t CsiSignal[3000]= {0};
  Double_t CsiTime[3000]  = {0};
  Double_t CsiEne[3000]   = {0}; 
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  
  TFile* tfOut = new TFile(Form("kl_Total_%s_GammaTime.root",name),"recreate");
  TTree* trOut = new TTree("trKL","GammaTime Adj");
  Double_t GammaCenterTime[6];
  Double_t GammaWeight[6];
  Double_t GammaTOF[6];
  Double_t GammaTotE[6];
  Double_t GammaTSigma[6];
  Double_t GammaTDeltaMax;
  trOut->Branch("GammaCenterTime",GammaCenterTime,"GammaCenterTime[6]/D");
  trOut->Branch("GammaWeight",GammaWeight,"GammaWeight[6]/D");
  trOut->Branch("GammaTOF",GammaTOF,"GammaTOF[6]/D");
  trOut->Branch("GammaTotE",GammaTotE,"GammaTotE[6]/D");
  trOut->Branch("GammaTSigma",GammaTSigma,"GammaTSigma[6]/D");
  trOut->Branch("GammaTDeltaMax",&GammaTDeltaMax,"GammaTDeltaMax/D");
  data.branchOfKlong( trOut );

  for( int ievt = 0; ievt< tr->GetEntries(); ievt++){
  //for( int ievt = 0; ievt< 100; ievt++){
    tr->GetEntry(ievt);
    std::vector<Klong> klVec;
    data.getData(klVec);
    GammaTDeltaMax = 0;
    for( int ig = 0; ig< 6; ig++){
      GammaCenterTime[ig] = 0;
      GammaWeight[ig] = 0;
      GammaTOF[ig] = 0;
      GammaTotE[ig] = 0;      
      GammaTSigma[ig] = 0;
    }

    int igamma = 0;
    for( int ipi = 0; ipi < klVec[0].pi0().size(); ipi++){
      double tmpTime[2];
      tmpTime[0] = klVec[0].pi0()[ipi].g1().t();
      tmpTime[1] = klVec[0].pi0()[ipi].g2().t();
      SetGammaTime( klVec[0].pi0()[ipi].g1());
      SetGammaTime( klVec[0].pi0()[ipi].g2());
      /*
      std::cout<< tmpTime[0] << "\t" <<   klVec[0].pi0()[ipi].g1().t() << std::endl;
      std::cout<< tmpTime[1] << "\t" <<   klVec[0].pi0()[ipi].g2().t() << std::endl;
      */
      GammaCenterTime[igamma] = GetTiming( klVec[0].pi0()[ipi].g1());
      GammaWeight[igamma]     = GetWeight( klVec[0].pi0()[ipi].g1());
      GammaTOF[igamma]        = CalGammaTOF( klVec[0],klVec[0].pi0()[ipi].g1());
      GammaTotE[igamma]       = klVec[0].pi0()[ipi].g1().e();
      GammaTSigma[igamma]     = GetClusterTSigma( klVec[0].pi0()[ipi].g1());
      igamma++;
      GammaCenterTime[igamma] = klVec[0].pi0()[ipi].g2().t();
      GammaWeight[igamma]     = GetWeight( klVec[0].pi0()[ipi].g2());
      GammaTOF[igamma]        = CalGammaTOF( klVec[0],klVec[0].pi0()[ipi].g2());
      GammaTotE[igamma]       = klVec[0].pi0()[ipi].g2().e();
      GammaTSigma[igamma]     = GetClusterTSigma( klVec[0].pi0()[ipi].g2());

      igamma++;
    }
    for( int ig = 0; ig < igamma-1; ig++){
      for( int jg  =ig+1; jg < igamma; jg++){
	if( TMath::Abs(GammaCenterTime[ig]-GammaTOF[ig] - (GammaCenterTime[jg] - GammaTOF[jg])) > GammaTDeltaMax){
	  GammaTDeltaMax = TMath::Abs(GammaCenterTime[ig]-GammaTOF[ig] - (GammaCenterTime[jg] - GammaTOF[jg]));
	}
      }
    }
    trOut->Fill();
  }
  trOut->Write();
  tfOut->Close();
}
  
