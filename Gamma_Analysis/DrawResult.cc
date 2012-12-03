#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TApplication.h"
#include "ClusterTimeReader.h"
#include "TMath.h"
int main( int argc, char** argv ){
  
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");

  TFile* tf = new TFile(Form("%s/Data_All.root",ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr = (TTree*)tf->Get("trCluster");
  ClusterTimeReader* reader = new ClusterTimeReader(tr);
  Int_t nEntries = reader->fChain->GetEntries(); 
  
  TFile* tfout = new TFile("ClusterTimeStructure.root","recreate");  

  const int nE = 6;
  const int nTheta = 8; 
  TH2D* hisRT[nE-1][nTheta-1];
  TH2D* hisRT_l[nE-1][nTheta-1];  
  TH2D* hisRT_r[nE-1][nTheta-1];  
  TH2D* hisDT[nE-1][nTheta-1];
  TH2D* hisDT_r[nE-1][nTheta-1];
  TH2D* hisDT_l[nE-1][nTheta-1];
  TH2D* hisPhiPhi[nE-1][nTheta -1 ];
  Int_t EArr[nE] = {100,200,300,500,800,1300};
  Int_t ThetaArr[nTheta] = {10,15,20,25,30,35,40,45};
  

  for( int iIndex = 0; iIndex < nE-1; iIndex++){
    for( int jIndex = 0; jIndex < nTheta-1; jIndex++){
      hisPhiPhi[iIndex][jIndex] = new TH2D(Form("hisPhiPhi_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  Form("hisPhiPhi_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  60,-1*TMath::Pi(),TMath::Pi(),60,-1*TMath::Pi(),TMath::Pi());
      hisRT[iIndex][jIndex]  = new TH2D(Form("hisRT_E_%d_%d_Theta_%d_%d",
					     EArr[iIndex],EArr[iIndex+1],
					     ThetaArr[jIndex],ThetaArr[jIndex+1]),
					Form("hisRT_E_%d_%d_Theta_%d_%d",
					     EArr[iIndex],EArr[iIndex+1],
					     ThetaArr[jIndex],ThetaArr[jIndex+1]),
					40,-100,100,200,-10,10);
      hisRT_l[iIndex][jIndex]  = new TH2D(Form("hisRT_l_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  Form("hisRT_l_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  40,-100,100,200,-10,10);
      hisRT_r[iIndex][jIndex]  = new TH2D(Form("hisRT_r_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  Form("hisRT_r_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  40,-100,100,200,-10,10);
      hisDT[iIndex][jIndex]  = new TH2D(Form("hisDT_E_%d_%d_Theta_%d_%d",
					     EArr[iIndex],EArr[iIndex+1],
					     ThetaArr[jIndex],ThetaArr[jIndex+1]),
					Form("hisDT_E_%d_%d_Theta_%d_%d",
					     EArr[iIndex],EArr[iIndex+1],
					     ThetaArr[jIndex],ThetaArr[jIndex+1]),
					40,-100,100,200,-10,10);
      hisDT_l[iIndex][jIndex]  = new TH2D(Form("hisDT_l_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  Form("hisDT_l_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  40,-100,100,200,-10,10);
      hisDT_r[iIndex][jIndex]  = new TH2D(Form("hisDT_r_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  Form("hisDT_r_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  40,-100,100,200,-10,10);
    }
  }

  for( int evtIndex = 0; evtIndex < nEntries; evtIndex++){
    reader->GetEntry( evtIndex );
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      if( reader->ClusterR[ clusterIndex] > 500 ){ continue; }
      bool bAbort = true; 
      Int_t EnergyIndex = 0;
      Int_t ThetaIndex  = 0;
      for( int EIndex = 0; EIndex < nE-1; EIndex++){
	if( reader->ClusterEnergy[ clusterIndex ] > EArr[ EIndex] &&
	    reader->ClusterEnergy[ clusterIndex ] <= EArr[ EIndex +1] ){
	  EnergyIndex  = EIndex;
	  bAbort = false;
	  break;
	}
      }
      if( bAbort ) { continue;}
      bAbort = true; 
      for( int TIndex = 0; TIndex < nTheta -1 ; TIndex++){
	if( reader->ClusterTheta[ clusterIndex ]*180./TMath::Pi() >  ThetaArr[ TIndex] &&
	    reader->ClusterTheta[ clusterIndex ]*180./TMath::Pi() <= ThetaArr[ TIndex +1 ]){
	  ThetaIndex = TIndex;
	  bAbort = false;
	  break;
	}
      }
      if( bAbort ){ continue; }

      bool blr = true;// view from downstream , l is false( -1,0 ), r is true( 1 )//
      if( TMath::Cos( reader->ClusterPhi[ clusterIndex ] ) < 0 ){ 
	blr = false;
      }else{ blr = true; }
      
      if( reader->nCrystal[clusterIndex] < 6 ){ continue; }
      if( reader->ClusterR[clusterIndex] > 500 ){ continue; }
      if( reader->ClusterR[clusterIndex] < 250 ){ continue; }
      for( int crystalIndex  = 0; crystalIndex < reader->nCrystal[ clusterIndex ]; crystalIndex++){       
	if( reader->CrystalEnergy[clusterIndex][crystalIndex] <= 3 ) { continue; }
	
	Double_t RinCluster = reader->CrystalR[ clusterIndex][crystalIndex]*TMath::Cos( reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t DinCluster = reader->CrystalR[ clusterIndex][crystalIndex]*TMath::Sin( reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t TinCluster = reader->CrystalT[ clusterIndex][crystalIndex];

	hisPhiPhi[EnergyIndex][ThetaIndex]->Fill( reader->ClusterPhi[clusterIndex], reader->CrystalPhi[clusterIndex][crystalIndex]);
	if( TMath::Abs(DinCluster) < 25*sqrt(2) ){
	  hisRT[EnergyIndex][ThetaIndex]->Fill( RinCluster, TinCluster);
	  if( !blr ){
	    hisRT_l[EnergyIndex][ThetaIndex]->Fill( RinCluster,TinCluster);
	  }else{
	    hisRT_r[EnergyIndex][ThetaIndex]->Fill( RinCluster,TinCluster);
	  }
	}	
	if( TMath::Abs(RinCluster) < 25*sqrt(2) ){
	  hisDT[EnergyIndex][ThetaIndex]->Fill( DinCluster, TinCluster);
	  if( !blr ){
	    hisDT_l[EnergyIndex][ThetaIndex]->Fill( DinCluster,TinCluster);
	  }else{
	    hisDT_r[EnergyIndex][ThetaIndex]->Fill( DinCluster,TinCluster);
	  }
	}
      }
    }
  }

  for( int iIndex = 0; iIndex < nE-1; iIndex++){
    for( int jIndex = 0; jIndex < nTheta-1; jIndex++){
      hisPhiPhi[iIndex][jIndex]->Write();
      hisRT[iIndex][jIndex]->Write();
      hisRT_l[iIndex][jIndex]->Write();
      hisRT_r[iIndex][jIndex]->Write();
      hisDT[iIndex][jIndex]->Write();
      hisDT_l[iIndex][jIndex]->Write();
      hisDT_r[iIndex][jIndex]->Write();
    }
  }

  tfout->Close();
  return 0; 
}
