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
#include "TProfile.h"
#include "TProfile2D.h"



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
  
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  TF1* TimeAdjFunc = new TF1("TimeAdjFunc",AdjFunc, 0, 2000,5);

  //Double_t Par[5] = {-0.0905327,1.54915,-0.114423,0.0758477,0.00487457};
  //Double_t ParErrors[5] = {0.00245834,0.0263153,0.00188467,0.007594,4.81501e-05};
  //Double_t Par[5] = {-0.105097,1.52645,-0.10655,0.0620572,0.00910542};
  //Double_t ParErrors[5]={0.00186334,0.0135129,0.000805812,0.000167933,1.18097e-05};
  //Double_t Par[5] = {-0.0644067,1.1759,-0.165316,0.0570758,0.0049958};
  //Double_t ParErrors[5] = {0.00203663,0.0418876,0.00542069,0.00221644,4.08696e-05};
  Double_t Par[5] = {-0.0976609,1.17012,-0.160851,0.0823451,0.00420403};
  Double_t ParErrors[5] = {0.0193245,0.0423249,0.00586115,0.0147254,0.000386874};
  TimeAdjFunc->SetParameters(Par);
  TimeAdjFunc->SetParErrors(ParErrors);

  TFile* tf = new TFile(Form("%s/Data_All.root",ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr = (TTree*)tf->Get("trCluster");
  ClusterTimeReader* reader = new ClusterTimeReader(tr);
  Int_t nEntries = reader->fChain->GetEntries(); 
  
  TFile* tfout = new TFile("CrystalPositionTimeEnergy_0.root","recreate");

  const int nClusterE = 10;
  const int nCrystalE = 10;
  const int nTheta    = 8;
  
  TH1D* hisRESpectrum[10];
  TH1D* hisRTSpectrum[10][10];
  for( int RIndex  = 0; RIndex < 10; RIndex++){
    hisRESpectrum[RIndex]  = new TH1D(Form("hisRESpectrum%d",RIndex),Form("hisRESpectrum%d",RIndex),
				     200,0,500);
    for( int EIndex  = 0; EIndex  < 10; EIndex++){
      hisRTSpectrum[RIndex][EIndex]  = new TH1D(Form("hisRTSpectrum_%d_%d",RIndex,EIndex),
						Form("hisRTSpectrum_%d_%d",RIndex,EIndex),
						200,0,500);
    }   
  }
  
  nEntries = reader->fChain->GetEntries();
  for( int evtIndex = 0; evtIndex < nEntries; evtIndex++){
    reader->GetEntry( evtIndex );
    for( int clusterIndex  =0 ;clusterIndex < reader->nCluster; clusterIndex++){
  
      if( reader->nCrystal[ clusterIndex ] < 6 ){ continue; }
      if( reader->ClusterR[ clusterIndex ] > 500 ){ continue; }
      if( reader->ClusterR[ clusterIndex ] < 200 ){ continue; }

      //if( !blr ){ continue; }
      

      Double_t RCenterCrystal = 1000;
      Double_t ECenterCrystal = 0;
      Double_t TCenterCrystal = 0;
      Double_t EMaximum       = 0;
      for( int crystalIndex = 0; crystalIndex < reader->nCrystal[ clusterIndex]; crystalIndex++){
	Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t RinCluster = reader->CrystalR[clusterIndex][crystalIndex]*TMath::Cos(reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t DinCluster = reader->CrystalR[clusterIndex][crystalIndex]*TMath::Sin(reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t EinCluster = reader->CrystalEnergy[clusterIndex][crystalIndex];
	// Time-Energy relation fixed ? // 
	Double_t TinCluster = reader->CrystalT[clusterIndex][crystalIndex]-TimeAdjFunc->Eval(EinCluster);
	if( EMaximum < EinCluster ){
	  EMaximum = EinCluster;
	}
	if( RadinCluster < RCenterCrystal){ 
	  RCenterCrystal = RadinCluster;
	  ECenterCrystal = EinCluster;
	  TCenterCrystal = TinCluster;
	}
      }

      for( int crystalIndex = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	Double_t PhiinCluster = reader->CrystalPhi[clusterIndex][crystalIndex];
	Double_t CosTheta     = TMath::Cos(PhiinCluster);
	Double_t SinTheta     = TMath::Sin(PhiinCluster);
	Double_t R            = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t RinCluster   = R*CosTheta;
	Double_t DinCluster   = R*SinTheta;
	Double_t EinCluster   = reader->CrystalEnergy[clusterIndex][crystalIndex];
	Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TimeAdjFunc->Eval(EinCluster)-TCenterCrystal;
	Int_t histIndex  = (int)(((int)(R))/10);
	Int_t EIndex     = (int)(((int)(EinCluster))/50);
	if( histIndex < 0  || histIndex >= 10 ){  // WrongIndex // 
	  continue;
	}
	if( EIndex < 0 || EIndex >=10 ){ // WrongIndex // 
	  continue;
	}
	if( evtIndex < 10 ){
	  std::cout<< histIndex << " : " << R << std::endl;
	}
	hisRESpectrum[histIndex]->Fill(EinCluster);
	hisRTSpectrum[histIndex][EIndex]->Fill(TinCluster);
      }
    }
  }
  for( int hisIndex  = 0; hisIndex < 10; hisIndex++){
    hisRESpectrum[hisIndex]->Write();
    for( int EIndex  = 0; EIndex < 10; EIndex++){
      hisRTSpectrum[hisIndex][EIndex]->Write();
    }
  }
  tfout->Close();

  return 0;
}
