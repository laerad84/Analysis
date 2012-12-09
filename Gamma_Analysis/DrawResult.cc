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
  TH2D* hisRPhi[nE-1][nTheta -1 ];
  TH2D* hisRPhiChisqN[nE-1][nTheta -1];
  TH2D* hisRPhiE[nE-1][nTheta -1];

  
  const int nBinsR   = 41;
  const int nBinsPhi = 41;
  
  TH1D* hisRPhiTime[nE-1][nTheta-1][41][41];
  TH1D* hisRPhiEnergy[nE-1][nTheta-1][41][41];

  Int_t EArr[nE] = {100,200,300,500,800,1300};
  Int_t ThetaArr[nTheta] = {10,15,20,25,30,35,40,45};
  
  Int_t RArr[ nBinsR ];
  Int_t PhiArr[ nBinsPhi ];
  for( int kIndex  = 0; kIndex < nBinsR; kIndex++){
    RArr[ kIndex ] = -105+ 5*kIndex;
  }
  for( int lIndex = 0; lIndex < nBinsPhi; lIndex++){
    PhiArr[ lIndex ] = -105+5*lIndex;
  }
  
  for( int iIndex = 0; iIndex < nE-1; iIndex++){
    for( int jIndex = 0; jIndex < nTheta-1; jIndex++){
      for( int kIndex  = 0; kIndex  < nBinsR; kIndex++){
	for( int lIndex = 0; lIndex < nBinsPhi; lIndex++){
	  hisRPhiTime[ iIndex][jIndex][kIndex][lIndex] 
	    = new TH1D(Form("hisRPhiTime_%d_%d_%d_%d",iIndex,jIndex,kIndex,lIndex),
		       Form("hisRPhiTime_E_%d_%d_Theta_%d_%d_xy_%d_%d",
			    EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1],
			    RArr[kIndex],PhiArr[lIndex]),
		       200,-10,10);
	  hisRPhiEnergy[iIndex][jIndex][kIndex][lIndex]
	    = new TH1D(Form("hisRPhiEnergy_%d_%d_%d_%d",iIndex,jIndex,kIndex,lIndex),
		       Form("hisRPhiEnergy_E_%d_%d_Theta_%d_%d_xy_%d_%d",
			    EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1],
			    RArr[kIndex],PhiArr[lIndex]),
		       100,0,1000);					       
	}
      }
      hisRPhi[iIndex][jIndex] = new TH2D(Form("hisRPhi_E_%d_Theta_%d",
					      iIndex,jIndex),
					 Form("hisRPhi_E_%d_%d_Theta_%d_%d",
					      EArr[iIndex], EArr[iIndex+1],
					      ThetaArr[jIndex], ThetaArr[jIndex+1]),
					 41,-105,105,41,-105,105);
      hisRPhiE[iIndex][jIndex] = new TH2D(Form("hisRPhiE_E_%d_Theta_%d",
					      iIndex,jIndex),
					 Form("hisRPhiE_E_%d_%d_Theta_%d_%d",
					      EArr[iIndex], EArr[iIndex+1],
					      ThetaArr[jIndex], ThetaArr[jIndex+1]),
					 41,-105,105,41,-105,105);
      
      hisRPhiChisqN[iIndex][jIndex] = new TH2D(Form("hisRPhiChisqN_E_%d_%d_Theta_%d_%d",
						    iIndex,jIndex),
					       Form("hisRThetaChisqN_E_%d_%d_Theta_%d_%d",
						    EArr[iIndex], EArr[iIndex+1],
						    ThetaArr[jIndex], ThetaArr[jIndex+1]),
					       41,-105,105,41,-105,105);
      hisPhiPhi[iIndex][jIndex] = new TH2D(Form("hisPhiPhi_E_%d_Theta_%d",
						iIndex,jIndex),
					   Form("hisPhiPhi_E_%d_%d_Theta_%d_%d",
						EArr[iIndex],EArr[iIndex+1],
						ThetaArr[jIndex],ThetaArr[jIndex+1]),
					   60,-1*TMath::Pi(),TMath::Pi(),60,-1*TMath::Pi(),TMath::Pi());
      hisRT[iIndex][jIndex]  = new TH2D(Form("hisRT_E_%d_Theta_%d",
					     iIndex,jIndex), 
					Form("hisRT_E_%d_%d_Theta_%d_%d",
					     EArr[iIndex],EArr[iIndex+1],
					     ThetaArr[jIndex],ThetaArr[jIndex+1]),
					40,-100,100,200,-10,10);
      hisRT_l[iIndex][jIndex]  = new TH2D(Form("hisRT_l_E_%d_Theta_%d",
					       iIndex,jIndex),
					  Form("hisRT_l_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  40,-100,100,200,-10,10);
      hisRT_r[iIndex][jIndex]  = new TH2D(Form("hisRT_r_E_%d_Theta_%d",
					       iIndex,jIndex),
					  Form("hisRT_r_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  40,-100,100,200,-10,10);
      hisDT[iIndex][jIndex]  = new TH2D(Form("hisDT_E_%d_Theta_%d",
					     iIndex,jIndex),
					Form("hisDT_E_%d_%d_Theta_%d_%d",
					     EArr[iIndex],EArr[iIndex+1],
					     ThetaArr[jIndex],ThetaArr[jIndex+1]),
					40,-100,100,200,-10,10);
      hisDT_l[iIndex][jIndex]  = new TH2D(Form("hisDT_l_E_%d_Theta_%d",
					       iIndex,jIndex),
					  Form("hisDT_l_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  40,-100,100,200,-10,10);
      hisDT_r[iIndex][jIndex]  = new TH2D(Form("hisDT_r_E_%d_Theta_%d",
					       iIndex,jIndex),
					  Form("hisDT_r_E_%d_%d_Theta_%d_%d",
					       EArr[iIndex],EArr[iIndex+1],
					       ThetaArr[jIndex],ThetaArr[jIndex+1]),
					  40,-100,100,200,-10,10);
    }
  }

  for( int evtIndex = 0; evtIndex < nEntries; evtIndex++){
    reader->GetEntry( evtIndex );
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      if( reader->ClusterR[clusterIndex] > 500 ){ continue; }
      if( reader->ClusterR[clusterIndex] < 200 ){ continue; }

      bool bAbort = true; 
      Int_t EnergyIndex = 0;
      Int_t ThetaIndex  = 0;
      for( int EIndex = 0; EIndex < nE-1; EIndex++){
	if( reader->ClusterEnergy[clusterIndex] >  EArr[EIndex] &&
	    reader->ClusterEnergy[clusterIndex] <= EArr[EIndex +1] ){
	  EnergyIndex  = EIndex;
	  bAbort = false;
	  break;
	}
      }
      if( bAbort ) { continue;}
      bAbort = true; 
      for( int TIndex = 0; TIndex < nTheta -1 ; TIndex++){
	if( reader->ClusterTheta[clusterIndex]*180./TMath::Pi() >  ThetaArr[TIndex] &&
	    reader->ClusterTheta[clusterIndex]*180./TMath::Pi() <= ThetaArr[TIndex+1]){
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
	//if( reader->CrystalEnergy[clusterIndex][crystalIndex] > 400 ){ continue; }
	Double_t RinCluster = reader->CrystalR[clusterIndex][crystalIndex]*TMath::Cos(reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t DinCluster = reader->CrystalR[clusterIndex][crystalIndex]*TMath::Sin(reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t TinCluster = reader->CrystalT[clusterIndex][crystalIndex];
	Double_t EinCluster = reader->CrystalEnergy[clusterIndex][crystalIndex];

	if( EinCluster == 0){ continue; }
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
	  hisDT[EnergyIndex][ThetaIndex]->Fill( DinCluster,TinCluster);
	  if( !blr ){
	    hisDT_l[EnergyIndex][ThetaIndex]->Fill( DinCluster,TinCluster);
	  }else{
	    hisDT_r[EnergyIndex][ThetaIndex]->Fill( DinCluster,TinCluster);
	  }
	}
	int iBinR = (int)((RinCluster+105)/5) - 1;
	int iBinTheta = (int)((DinCluster+105)/5) -1;
	if( iBinR < 0 || iBinTheta < 0 ){
	  continue;
	}else if( iBinR >= nBinsR || iBinTheta >= nBinsPhi ){
	  continue;
	}else{
	  if( TinCluster != 0){
	    hisRPhiTime[EnergyIndex][ThetaIndex][iBinR][iBinTheta]->Fill(TinCluster);
	  }
	  hisRPhiEnergy[EnergyIndex][ThetaIndex][iBinR][iBinTheta]->Fill(EinCluster);
	}
      }
    }
  }
  

  Double_t TimeMean[ nE-1 ][ nTheta-1 ][ nBinsR ][ nBinsPhi ];
  Double_t TimeSigma[ nE-1 ][ nTheta-1 ][ nBinsR ][ nBinsPhi ];
  Double_t TimeChisqN[ nE-1 ][ nTheta-1 ][ nBinsR ][ nBinsPhi ];

  Double_t EnergyMean[ nE-1 ][ nTheta-1 ][ nBinsR ][ nBinsPhi ];
  Double_t EnergySigma[ nE-1 ][ nTheta-1 ][ nBinsR ][ nBinsPhi ];
  Double_t EnergyChisqN[ nE-1 ][ nTheta-1 ][ nBinsR ][ nBinsPhi ];
  for( int iIndex  = 0; iIndex <nE-1; iIndex++){
    for( int jIndex  = 0; jIndex < nTheta-1; jIndex++){
      for( int kIndex = 0; kIndex < nBinsR; kIndex++){
	for( int lIndex  = 0; lIndex < nBinsPhi; lIndex++){
	  TimeMean[ iIndex][ jIndex][ kIndex][ lIndex]   = -10;
	  TimeSigma[ iIndex][ jIndex][ kIndex][ lIndex]  = -1;
	  TimeChisqN[ iIndex][ jIndex][ kIndex][ lIndex] = -1;
	  Double_t Mean  =  hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->GetMean();
	  Double_t RMS   =  hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->GetRMS();	  
	  Bool_t Flag = kTRUE;
	  if( hisRPhiTime[ iIndex][jIndex][ kIndex][lIndex]->GetEntries() < 25 ){ 
	    Flag = kFALSE;
	  }else{
	    hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->Fit("gaus","Q","",Mean-RMS,Mean+RMS);
	    TF1* func = hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->GetFunction("gaus");
	    if( func->GetParameter(0 ) < 10 ){  
	      hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->GetListOfFunctions()->Delete();	      
	      Flag = kFALSE; 
	    }else{
	      Mean = func->GetParameter(1);
	      RMS  = func->GetParameter(2);	 	  
	      if( RMS < 0.5 ){ RMS  = 0.5;} 	      
	      hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->Fit("gaus","Q","",Mean-2*RMS,Mean+2*RMS);
	      func = hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->GetFunction("gaus");
	      TimeMean[ iIndex][ jIndex][ kIndex][ lIndex] = func->GetParameter( 1 );
	      if( TimeMean[ iIndex][ jIndex][ kIndex][ lIndex] > 2 || 
		  TimeMean[ iIndex][ jIndex][ kIndex][ lIndex] <-2 ){
		hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->GetListOfFunctions()->Delete();	      
		Flag = kFALSE;
	      }else{
		TimeSigma[ iIndex][ jIndex][ kIndex][ lIndex] = func->GetParameter( 2 );
		TimeChisqN[ iIndex][ jIndex][ kIndex][ lIndex] = func->GetChisquare()/func->GetNDF();	  
		hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->GetListOfFunctions()->Delete();
	      }
	    }
	  }
	  hisRPhi[iIndex][jIndex]->SetBinContent( lIndex+1, kIndex+1, TimeMean[iIndex][jIndex][kIndex][lIndex]);
	  hisRPhi[iIndex][jIndex]->SetBinError( lIndex+1, kIndex+1, TimeSigma[iIndex][jIndex][kIndex][lIndex]);	
	  hisRPhiChisqN[iIndex][jIndex]->SetBinContent( lIndex+1, kIndex+1, TimeChisqN[iIndex][jIndex][kIndex][lIndex]);
	}
      }
    }
  }
  for( int iIndex  = 0; iIndex <nE-1; iIndex++){
    for( int jIndex  = 0; jIndex < nTheta-1; jIndex++){
      for( int kIndex = 0; kIndex < nBinsR; kIndex++){
	for( int lIndex  = 0; lIndex < nBinsPhi; lIndex++){
	  EnergyMean[iIndex][jIndex][kIndex][lIndex] = hisRPhiEnergy[iIndex][jIndex][kIndex][lIndex]->GetMean();
	  EnergySigma[iIndex][jIndex][kIndex][lIndex] = hisRPhiEnergy[iIndex][jIndex][kIndex][lIndex]->GetRMS();
	  hisRPhiE[iIndex][jIndex]->SetBinContent( lIndex+1, kIndex+1, EnergyMean[iIndex][jIndex][kIndex][lIndex]);
	  hisRPhiE[iIndex][jIndex]->SetBinError( lIndex+1, kIndex+1, EnergySigma[iIndex][jIndex][kIndex][lIndex]);	
	  //hisRPhiChisqN[iIndex][jIndex]->SetBinContent( lIndex+1, kIndex+1, EnergyChisqN[iIndex][jIndex][kIndex][lIndex]);
	  std::cout << iIndex << ":" << jIndex << ":" << kIndex << ":" << lIndex << " : " << EnergyMean[iIndex][jIndex][kIndex][lIndex] << std::endl;
	}
      }
    }
  }
  
  for( int iIndex = 0; iIndex < nE-1; iIndex++){
    for( int jIndex = 0; jIndex < nTheta-1; jIndex++){
      for( int kIndex = 0; kIndex < nBinsR; kIndex++){
	for( int lIndex = 0; lIndex < nBinsPhi; lIndex++){
	  hisRPhiTime[iIndex][jIndex][kIndex][lIndex]->Write();
	  hisRPhiEnergy[iIndex][jIndex][kIndex][lIndex]->Write();
	}
      }
      hisRPhi[iIndex][jIndex]->Write();
      hisRPhiChisqN[iIndex][jIndex]->Write();      
      hisRPhiE[iIndex][jIndex]->Write();
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
