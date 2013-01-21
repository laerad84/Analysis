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
#include "TProfile2D.h"

int main( int argc, char** argv ){

  Int_t Energy = atoi( argv[1] );
  Int_t Theta  = atoi( argv[2] );
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string iFileForm = "%s/Cluster_Time_%dMeV_%ddeg-1E5-%d.root";
  std::string oFileForm = "%s/ClusterTimeStructure_SIM_%dMeV_%ddeg.root";


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// I/O File Setting ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TChain* ch = new TChain("trCluster");
  for( int iIndex =0; iIndex < 10; iIndex++){
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

  TH1D* hisRPhiTime[ nBinsR ][ nBinsPhi ];
  TH2D* hisRPhi;
  TH2D* hisRPhiChisqN;

  TH2D* hisRT;
  TH2D* hisRT_l;
  TH2D* hisRT_r;
  TH2D* hisDT;
  TH2D* hisDT_r;
  TH2D* hisDT_l;
  TH2D* hisPhiPhi;
  for( int kIndex = 0; kIndex < nBinsR; kIndex++){
    for( int lIndex  = 0; lIndex <nBinsPhi; lIndex++){
      hisRPhiTime[ kIndex ][ lIndex ] = new TH1D(Form("hisRPhiTime_E_%d_Theta_%d_%d_%d",
						      Energy,Theta,kIndex,lIndex),
						 Form("hisRPhiTime_E_%d_Theta_%d_%d_%d",
						      Energy,Theta,kIndex,lIndex),
						 200, -10, 10);
    }
  }
  hisRPhi   = new TH2D(Form("hisRPhi_E_%d_%d_Theta_%d_%d",
			    Energy,Energy,Theta,Theta),
		       Form("hisRPhi_E_%d_%d_Theta_%d_%d",
			    Energy,Energy,Theta,Theta),
		       41,-105,105,41,-105,105);
  hisRPhiChisqN = new TH2D(Form("hisRPhiChisqN_E_%d_%d_Theta_%d_%d",
						Energy,Energy,Theta,Theta),
					   Form("hisRPhiChisqN_E_%d_%d_Theta_%d_%d",
						Energy,Energy,Theta,Theta),
					   41,-105,105,41,-105,105);
  
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
	if( TinCluster == 0 ){ continue; }
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
	int iBinR = (int)((RinCluster+105)/5) -1;
	int iBinPhi = (int)((DinCluster+105)/5) -1;
	if( iBinR <0 || iBinPhi < 0) { continue; }
	else if( iBinR >= nBinsR || iBinPhi >= nBinsPhi ){continue; }
	else{
	  hisRPhiTime[iBinR][iBinPhi]->Fill(TinCluster);
	}			  
      }
    }
  }

  Double_t TimeMean[ nBinsR ][ nBinsPhi ];
  Double_t TimeSigma[ nBinsR ][ nBinsPhi ];
  Double_t TimeChisqN[ nBinsR ][ nBinsPhi ];
  for( int kIndex  = 0; kIndex < nBinsR; kIndex++){
    for( int lIndex  = 0; lIndex < nBinsPhi; lIndex++){
      TimeMean[ kIndex][ lIndex]   = -10;
      TimeSigma[ kIndex][ lIndex]  = -1;
      TimeChisqN[ kIndex][ lIndex] = -1;
      Double_t Mean  =  hisRPhiTime[kIndex][lIndex]->GetMean();
      Double_t RMS   =  hisRPhiTime[kIndex][lIndex]->GetRMS();	  
      Bool_t Flag = kTRUE;
      if( hisRPhiTime[kIndex][lIndex]->GetEntries() < 25 ){ 
	Flag = kFALSE;
      }else{
	hisRPhiTime[kIndex][lIndex]->Fit("gaus","Q","",Mean-RMS,Mean+RMS);
	TF1* func = hisRPhiTime[kIndex][lIndex]->GetFunction("gaus");
	if( func->GetParameter(0 ) < 10 ){  
	  hisRPhiTime[kIndex][lIndex]->GetListOfFunctions()->Delete();	      
	  Flag = kFALSE; 
	}else{
	  Mean = func->GetParameter(1);
	  RMS  = func->GetParameter(2);	 	  
	  if( RMS < 0.5 ){ RMS  = 0.5;} 	      
	  hisRPhiTime[kIndex][lIndex]->Fit("gaus","Q","",Mean-2*RMS,Mean+2*RMS);
	  func = hisRPhiTime[kIndex][lIndex]->GetFunction("gaus");
	  TimeMean[ kIndex][ lIndex] = func->GetParameter( 1 );
	  if( TimeMean[ kIndex][ lIndex] > 2 || 
	      TimeMean[ kIndex][ lIndex] <-2 ){
	    hisRPhiTime[kIndex][lIndex]->GetListOfFunctions()->Delete();	      
	    Flag = kFALSE;
	  }else{
	    TimeSigma[ kIndex][ lIndex] = func->GetParameter( 2 );
	    TimeChisqN[ kIndex][ lIndex] = func->GetChisquare()/func->GetNDF();	  
	    hisRPhiTime[kIndex][lIndex]->GetListOfFunctions()->Delete();
	  }
	}
      }
      hisRPhi->SetBinContent( lIndex+1, kIndex+1, TimeMean[kIndex][lIndex]);
      hisRPhi->SetBinError( lIndex+1, kIndex+1, TimeSigma[kIndex][lIndex]);	
      hisRPhiChisqN->SetBinContent( lIndex+1, kIndex+1, TimeChisqN[kIndex][lIndex]);
    }
  }

  for( int kIndex = 0; kIndex < nBinsR; kIndex++){
    for( int lIndex = 0; lIndex < nBinsPhi; lIndex++){
      hisRPhiTime[kIndex][lIndex] ->Write();
    }
  }

  hisRPhi->Write();
  hisRPhiChisqN->Write();
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
