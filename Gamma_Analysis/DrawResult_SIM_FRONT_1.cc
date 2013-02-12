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
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <limits>
int main( int argc, char** argv ){

  Int_t Energy = atoi( argv[1] );
  Int_t Theta  = atoi( argv[2] );
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string iFileForm = "%s/Cluster_Time_%dMeV_%ddeg-1E5-%d.root";
  std::string oFileForm = "%s/ClusterTimeStructure_SIM_1_%dMeV_%ddeg.root";

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// I/O File Setting ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TChain* ch = new TChain("trCluster");
  for( int iIndex =0; iIndex < 2; iIndex++){
    ch->Add(Form(iFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str(),Energy,Theta,iIndex));
  }
  /*
  TFile* tf = new TFile(Form("%s/Data_All.root",ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr = (TTree*)tf->Get("trCluster");
  */
  ClusterTimeReader* reader = new ClusterTimeReader(ch);
  Int_t nEntries = reader->fChain->GetEntries();   
  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str(),Energy,Theta),"recreate");  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const int nE = 6;
  const int nTheta = 8; 
  const int nBinsR = 41;
  const int nBinsPhi = 41;
  Double_t width = 12.5;
  Int_t    nBinsRD  = 19;
  Double_t RDMax = width*nBinsRD/2.;
  Double_t RDMin = -1*width*nBinsRD/2.; 

  TProfile* profRT = new TProfile(Form("profRT_%d_%d",Energy,Theta),
				  Form("profRT_%d_%d",Energy,Theta),
				  nBinsRD, RDMin,RDMax);
  TH2D* hisRT = new TH2D(Form("hisRT_%d_%d",Energy,Theta),
			 Form("hisRT_%d_%d",Energy,Theta),
			 nBinsRD,RDMin,RDMax,200,-10,10);
  TProfile* profDT = new TProfile(Form("profDT_%d_%d",Energy,Theta),
				  Form("profDT_%d_%d",Energy,Theta),
				  nBinsRD, RDMin,RDMax);
  TH2D* hisDT = new TH2D(Form("hisDT_%d_%d",Energy,Theta),
			 Form("hisDT_%d_%d",Energy,Theta),
			 nBinsRD,RDMin,RDMax,200,-10,10);
  
  for( int eventIndex = 0; eventIndex < nEntries; eventIndex++){
    reader->GetEntry( eventIndex );
     if( eventIndex %1000 == 0 && eventIndex != 1000){ std::cout<< eventIndex <<"/" << nEntries<< std::endl;}
    //if( reader->nCluster > 1 ){ continue; }
    
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){      
      if( reader->ClusterR[clusterIndex] > 500 ){continue;}
      if( reader->ClusterR[clusterIndex] < 200 ){continue;}
      if( reader->nCrystal[clusterIndex] < 4 ){ continue; }

      Double_t RCenterCrystal = 1000;
      Double_t ECenterCrystal = 0;
      Double_t TCenterCrystal = 0;
      Double_t RCenter = 0;
      Double_t DCenter = 0;
      for( int crystalIndex  = 0; crystalIndex < reader->nCrystal[ clusterIndex ]; crystalIndex++){       
	if( reader->CrystalEnergy[clusterIndex][crystalIndex] <= 100){ continue; }	
	
	// Get CenterCrystal Data
	for( int crystalIndex = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	  Double_t PhiinCluster = reader->CrystalPhi[clusterIndex][crystalIndex];
	  Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	  Double_t RinCluster   = RadinCluster*TMath::Cos(PhiinCluster);
	  Double_t DinCluster   = RadinCluster*TMath::Sin(PhiinCluster);
	  Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex];
	  Double_t EinCluster   = reader->CrystalEnergy[clusterIndex][crystalIndex];
	  if( RadinCluster < TMath::Abs(RCenterCrystal)){
	    RCenterCrystal = RadinCluster;
	    ECenterCrystal = EinCluster;
	    TCenterCrystal = TinCluster;
	    RCenter = RinCluster;
	    DCenter = DinCluster;
	  }
	}
      }    

      for( int crystalIndex = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t PhiinCluster = reader->CrystalPhi[clusterIndex][crystalIndex];
	Double_t RinCluster   = RadinCluster*TMath::Cos(PhiinCluster);
	Double_t DinCluster   = RadinCluster*TMath::Sin(PhiinCluster);
	Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TCenterCrystal;
	Double_t EinCluster   = reader->CrystalEnergy[clusterIndex][crystalIndex];
	if( TinCluster > 20  || TinCluster < -20  ){ continue;}
	if( RinCluster > 200 || RinCluster < -200 ){ continue;}
	if( DinCluster > 200 || DinCluster < -200 ){ continue;}
	if( EinCluster == 0 ){ continue; }
	if( TinCluster == 0 ){ continue; }
	if( RinCluster == std::numeric_limits<double>::quiet_NaN()){continue;}
	if( DinCluster == std::numeric_limits<double>::quiet_NaN()){continue;}
	if( TinCluster == std::numeric_limits<double>::quiet_NaN()){continue;}
	if( TinCluster != TinCluster){ continue; }
	if( TMath::Abs(DinCluster < 25 )){ 
	  profRT->Fill(RinCluster,TinCluster);
	  hisRT->Fill(RinCluster,TinCluster);
	}
	if( TMath::Abs(RinCluster < 25 )){ 
	  profDT->Fill(DinCluster,1.0);//,TinCluster);
	  hisDT->Fill(DinCluster,TinCluster);
	}
	
      }
    }    
  }//
  std::cout<< profRT->GetEntries() << std::endl;
  profRT->Write();
  hisRT->Write();
  profDT->Write();
  hisDT->Write();
  tfout->Close();

}
