#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TApplication.h"
#include "ClusterTimeReader.h"
#include "TMath.h"
#include "TChain.h"


int main( int argc, char** argv ){

  Int_t Energy = atoi( argv[1] );
  Int_t Theta  = atoi( argv[2] );
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  
  TChain* ch = new TChain("trCluster");
  for( int iIndex =0; iIndex < 10; iIndex++){
    ch->Add(Form("%s/Cluster_Time_Back_%dMeV_%ddeg-1E5-%d.root",ROOTFILE_GAMMACLUS.c_str(),Energy,Theta,iIndex));
  }

  /*
  TFile* tf = new TFile(Form("%s/Data_All.root",ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr = (TTree*)tf->Get("trCluster");
  */
  ClusterTimeReader* reader = new ClusterTimeReader(ch);  
  Int_t nEntries = reader->fChain->GetEntries(); 
  
  TFile* tfout = new TFile(Form("%s/ClusterTimeStructure_SIM_Back_%dMeV_%ddeg.root",ROOTFILE_GAMMACLUS.c_str(),Energy,Theta),"recreate");  

  const int nE = 6;
  const int nTheta = 8; 
  TH2D* hisRT;
  TH2D* hisRT_l;
  TH2D* hisRT_r;
  TH2D* hisDT;
  TH2D* hisDT_r;
  TH2D* hisDT_l;
  TH2D* hisPhiPhi;
  hisPhiPhi = new TH2D(Form("hisPhiPhi_E_%d_%d_Theta_%d_%d",
			    Energy,Energy,
			    Theta,Theta),
		       Form("hisPhiPhi_E_%d_%d_Theta_%d_%d",
			    Energy,Energy,
			    Theta,Theta),
		       60,-1*TMath::Pi(),TMath::Pi(),60,-1*TMath::Pi(),TMath::Pi());
  hisRT  = new TH2D(Form("hisRT_E_%d_%d_Theta_%d_%d",
			 Energy,Energy,
			 Theta,Theta),
		    Form("hisRT_E_%d_%d_Theta_%d_%d",
			 Energy,Energy,
			 Theta,Theta),
		    40,-100,100,200,-10,10);
  hisRT_l  = new TH2D(Form("hisRT_l_E_%d_%d_Theta_%d_%d",
			   Energy,Energy,
			   Theta,Theta),
		      Form("hisRT_l_E_%d_%d_Theta_%d_%d",
			   Energy,Energy,
			   Theta,Theta),
		      40,-100,100,200,-10,10);
  hisRT_r  = new TH2D(Form("hisRT_r_E_%d_%d_Theta_%d_%d",
			   Energy,Energy,
			   Theta,Theta),
		      Form("hisRT_r_E_%d_%d_Theta_%d_%d",
			   Energy,Energy,
			   Theta,Theta),
		      40,-100,100,200,-10,10);
  hisDT  = new TH2D(Form("hisDT_E_%d_%d_Theta_%d_%d",
			 Energy,Energy,
			 Theta,Theta),
		    Form("hisDT_E_%d_%d_Theta_%d_%d",
			 Energy,Energy,
			 Theta,Theta),
		    40,-100,100,200,-10,10);
  hisDT_l  = new TH2D(Form("hisDT_l_E_%d_%d_Theta_%d_%d",
			   Energy,Energy,
			   Theta,Theta),
		      Form("hisDT_l_E_%d_%d_Theta_%d_%d",
			   Energy,Energy,
			   Theta,Theta),
		      40,-100,100,200,-10,10);
  hisDT_r  = new TH2D(Form("hisDT_r_E_%d_%d_Theta_%d_%d",
			   Energy,Energy,
			   Theta,Theta),
		      Form("hisDT_r_E_%d_%d_Theta_%d_%d",
			   Energy,Energy,
			   Theta,Theta),
		      40,-100,100,200,-10,10);
  for( int evtIndex = 0; evtIndex < nEntries; evtIndex++){
    reader->GetEntry( evtIndex );
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      if( reader->ClusterR[ clusterIndex] > 500 ){ continue; }

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
	
	hisPhiPhi->Fill( reader->ClusterPhi[clusterIndex], reader->CrystalPhi[clusterIndex][crystalIndex]);
	if( TMath::Abs(DinCluster) < 25*sqrt(2) ){
	  hisRT->Fill( RinCluster, TinCluster);
	  if( !blr ){
	    hisRT_l->Fill( RinCluster,TinCluster);
	  }else{
	    hisRT_r->Fill( RinCluster,TinCluster);
	  }
	}	
	if( TMath::Abs(RinCluster) < 25*sqrt(2) ){
	  hisDT->Fill( DinCluster, TinCluster);
	  if( !blr ){
	    hisDT_l->Fill( DinCluster,TinCluster);
	  }else{
	    hisDT_r->Fill( DinCluster,TinCluster);
	  }
	}
      }
    }
  }

  hisPhiPhi->Write();
  hisRT->Write();
  hisRT_l->Write();
  hisRT_r->Write();
  hisDT->Write();
  hisDT_l->Write();
  hisDT_r->Write();


  tfout->Close();
  return 0; 
}
