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

int main( int argc, char** argv ){
  
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");

  TFile* tf = new TFile(Form("%s/Data_All.root",ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr = (TTree*)tf->Get("trCluster");
  ClusterTimeReader* reader = new ClusterTimeReader(tr);
  Int_t nEntries = reader->fChain->GetEntries(); 
  
  TFile* tfout = new TFile("CrystalTimeEnergy.root","recreate");  

  const int nClusterE  = 10;
  const int nCrystalE  = 18;
  const int nTheta     = 8;
  const int nBinsR     = 41;
  const int nBinsPhi   = 41;
  
  TH2D* hisTimeEnergyClusterE[nClusterE-1][nTheta-1];
  TH2D* hisTimeEnergyCrystalE[nCrystalE-1][nTheta-1];

  TProfile* profTimeEnergyClusterE[nClusterE-1][nTheta-1];
  TProfile* profTimeEnergyCrystalE[nCrystalE-1][nTheta-1];

  TProfile2D* profTimeEnergyClusterCrystalE[nTheta-1];
  TProfile2D* profEnergyPosition[nTheta-1];
  TProfile2D* profTimePosition[nTheta-1];

  TProfile* profTimeEnergyLowEnergy = new TProfile("profTimeEnergyLowEnergy",
						   "profTimeEnergyLowEnergy",
						   125,0,500);

  TProfile* profTimeEnergyCrystalE_Wide[nCrystalE-1][nTheta-1];
  TProfile* profTimeEnergyCrystalEMerge[nCrystalE-1];
  Int_t EClusterArr[nClusterE] = {100,200,300,400,500,600,700,800,900,1000};
  Int_t ECrystalArr[nCrystalE] = {0,12.5,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400};
  Int_t ThetaArr[nTheta]       = {10,15,20,25,30,35,40,45};
  Int_t RArr[nBinsR];
  Int_t PhiArr[nBinsPhi];
  for( int kIndex  = 0; kIndex < nBinsR; kIndex++){
    RArr[kIndex] = -105+5*kIndex;
  }
  for( int lIndex = 0; lIndex < nBinsPhi; lIndex++){
    PhiArr[lIndex] = -105+5*lIndex;
  }

  for( int iIndex = 0 ; iIndex < nCrystalE-1; iIndex++){
    profTimeEnergyCrystalEMerge[iIndex] = new TProfile(Form("profTimeEnergyCrystalMerge_%d",iIndex),
						       Form("profTimeEnergyCrystalMerge_%d_%d",
							    ECrystalArr[iIndex],ECrystalArr[iIndex+1]),
						       20,0,500);						       
  }

  for( int jIndex = 0; jIndex < nTheta-1; jIndex++){
    profTimeEnergyClusterCrystalE[jIndex] = new TProfile2D(Form("profTimeEnergyClusterCrystalE_%d",jIndex),
						     Form("profTimeEnergyClusterCrystalE_%d_%d;ClusterE;CrystalE",ThetaArr[jIndex],ThetaArr[jIndex+1]),
						     40,0,1000,125,0,500);
    profEnergyPosition[jIndex] = new TProfile2D(Form("profEnergyPosition_%d",jIndex),
						Form("profEnergyPosition_%d_%d;X;Y",ThetaArr[jIndex],ThetaArr[jIndex+1]),
						80,-100,100,80,-100,100);
    profTimePosition[jIndex] = new TProfile2D(Form("profTimePosition_%d",jIndex),
					      Form("profTimePosition_%d_%d;X;Y",ThetaArr[jIndex],ThetaArr[jIndex+1]),
					      80,-100,100,80,-100,100);
    
    for( int iIndex = 0; iIndex < nClusterE-1; iIndex++){
      
      hisTimeEnergyClusterE[iIndex][jIndex] = new TH2D(Form("hisTimeEnergyCluster_%d_%d",iIndex,jIndex),
						       Form("hisTimeEnergyCluster_%d_%d_%d_%d",
							    EClusterArr[iIndex],EClusterArr[iIndex+1],
							    ThetaArr[jIndex],ThetaArr[jIndex+1]),
						       125,0,500,250,-5,5);      
      profTimeEnergyClusterE[iIndex][jIndex] = new TProfile(Form("profTimeEnergyCluster_%d_%d",iIndex,jIndex),
							    Form("profTimeEnergyCluster_%d_%d_%d_%d",
							    EClusterArr[iIndex],EClusterArr[iIndex+1],
							    ThetaArr[jIndex],ThetaArr[jIndex+1]),
							    125,0,525);      
    }
    for( int iIndex = 0; iIndex < nCrystalE-1; iIndex++){
      hisTimeEnergyCrystalE[iIndex][jIndex] = new TH2D(Form("hisTimeEnergyCrystal_%d_%d",iIndex,jIndex),
						       Form("hisTimeEnergyCrystal_%d_%d_%d_%d",
							    ECrystalArr[iIndex],ECrystalArr[iIndex+1],
							    ThetaArr[jIndex],ThetaArr[jIndex+1]),
						       125,0,500,250,-5,5);      
      
      profTimeEnergyCrystalE[iIndex][jIndex] = new TProfile(Form("profTimeEnergyCrystal_%d_%d",iIndex,jIndex),
							    Form("profTimeEnergyCrystal_%d_%d_%d_%d",
								 ECrystalArr[iIndex],ECrystalArr[iIndex+1],
								 ThetaArr[jIndex],ThetaArr[jIndex+1]),
							    125,0,500);      
      
      profTimeEnergyCrystalE_Wide[iIndex][jIndex] = new TProfile(Form("profTimeEnergyCrystal_Wide_%d_%d",iIndex,jIndex),
								 Form("profTimeEnergyCrystal_Wide_%d_%d_%d_%d",
								 ECrystalArr[iIndex],ECrystalArr[iIndex+1],
								 ThetaArr[jIndex],ThetaArr[jIndex+1]),
								 20,0,500);      
      
    }
  }
  

  for( int evtIndex  = 0; evtIndex < nEntries; evtIndex++){
    reader->GetEntry( evtIndex ); 
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      if( reader->ClusterR[clusterIndex] > 500 ){ continue; }
      if( reader->ClusterR[clusterIndex] < 200 ){ continue; }
      bool bAbort  =true;
      Int_t EnergyIndex;
      Int_t ThetaIndex;
      Int_t ECrystalIndex; 
      for( int EIndex  = 0 ; EIndex < nClusterE-1; EIndex++){
	if( reader->ClusterEnergy[clusterIndex] > EClusterArr[EIndex] && 
	    reader->ClusterEnergy[clusterIndex] <=EClusterArr[EIndex+1] ){
	  EnergyIndex = EIndex;
	  bAbort = false;
	  break;
	}		
      }
      if( bAbort ){ continue; }
      bAbort = true; 
      for( int TIndex  = 0; TIndex < nTheta-1; TIndex ++){
	if( reader->ClusterTheta[clusterIndex]*180./TMath::Pi() > ThetaArr[TIndex] &&
	    reader->ClusterTheta[clusterIndex]*180./TMath::Pi() <=ThetaArr[TIndex+1] ){
	  ThetaIndex = TIndex;
	  bAbort  = false;
	  break;
	}
      }
      if( bAbort ){ continue; }
      bool blr = true;
      if( TMath::Cos( reader->ClusterPhi[clusterIndex] ) < 0 ){
	blr = false;
      }else{
	blr = true; 
      }

      bAbort = true; 
      Double_t CrystalMaxE = 0;      
      Double_t MinimumR    = 1000;
      for( int crystalIndex  = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	if( MinimumR > reader->CrystalR[clusterIndex][crystalIndex] ){
	  MinimumR = reader->CrystalR[clusterIndex][crystalIndex];
	  CrystalMaxE = reader->CrystalEnergy[clusterIndex][crystalIndex];
	}
      }
      for( int ECryIndex  = 0; ECryIndex < nCrystalE -1 ; ECryIndex++){
	if( CrystalMaxE > ECrystalArr[ECryIndex] &&
	    CrystalMaxE <=ECrystalArr[ECryIndex+1] ){
	  ECrystalIndex = ECryIndex;
	  bAbort = false;
	  break;
	}
      }
      if( bAbort ){ continue; }

      if( reader->nCrystal[ clusterIndex ] < 6 ){ continue; }
      if( reader->ClusterR[ clusterIndex ] > 500 ){ continue; }
      if( reader->ClusterR[ clusterIndex ] < 200 ){ continue; }

      if( !blr ){ continue; }
      if( CrystalMaxE < reader->ClusterEnergy[clusterIndex]/3. ){ continue; }
      for( int crystalIndex  = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	//if( reader->CrystalEnergy[clusterIndex][crystalIndex] > reader->ClusterEnergy[clusterIndex]/3. ){ continue; }
	if( reader->CrystalEnergy[clusterIndex][crystalIndex] < 3 ){ continue; }
	Double_t CosTheta = TMath::Cos( reader->CrystalPhi[clusterIndex][crystalIndex] );
	Double_t SinTheta = TMath::Sin( reader->CrystalPhi[clusterIndex][crystalIndex] );
	Double_t R        = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t RinCluster = R*CosTheta;
	Double_t DinCluster = R*SinTheta;
	Double_t TinCluster = reader->CrystalT[clusterIndex][crystalIndex];
	Double_t EinCluster = reader->CrystalEnergy[clusterIndex][crystalIndex];
	profEnergyPosition[ThetaIndex]->Fill(RinCluster,DinCluster,EinCluster);
	profTimePosition[ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);
	if( EinCluster == 0){ continue; }
	if( TinCluster == 0){ continue; }
	if( R < 30 && RinCluster< 25){
	  hisTimeEnergyClusterE[EnergyIndex][ThetaIndex]   ->Fill(EinCluster, TinCluster);
	  hisTimeEnergyCrystalE[ECrystalIndex][ThetaIndex] ->Fill(EinCluster, TinCluster);
	  profTimeEnergyClusterE[EnergyIndex][ThetaIndex]  ->Fill(EinCluster, TinCluster);
	  profTimeEnergyCrystalE[ECrystalIndex][ThetaIndex]->Fill(EinCluster, TinCluster);
	  profTimeEnergyClusterCrystalE[ThetaIndex]        ->Fill(reader->ClusterEnergy[clusterIndex],EinCluster,TinCluster);
	  profTimeEnergyCrystalE_Wide[ECrystalIndex][ThetaIndex]->Fill(EinCluster,TinCluster);
	  if( ThetaIndex < 3 ){
	    profTimeEnergyCrystalEMerge[ECrystalIndex]->Fill( EinCluster,TinCluster);
	    if( EnergyIndex == 3 ){
	      profTimeEnergyLowEnergy->Fill( EinCluster,TinCluster);
	    }
	  }
	}

      }      
    }
  }


  for( int jIndex=0; jIndex < nTheta-1; jIndex++){
    for( int iIndex  =0; iIndex < nClusterE-1 ; iIndex++){
      hisTimeEnergyClusterE[iIndex][jIndex]->Write();
      profTimeEnergyClusterE[iIndex][jIndex]->SetLineColor(iIndex+1);
      profTimeEnergyClusterE[iIndex][jIndex]->Write();
    }
    for( int iIndex  =0; iIndex < nCrystalE-1; iIndex++){
      hisTimeEnergyCrystalE[iIndex][jIndex]->Write();
      profTimeEnergyCrystalE[iIndex][jIndex]->SetLineColor(iIndex+1);
      profTimeEnergyCrystalE[iIndex][jIndex]->Write();
      profTimeEnergyCrystalE_Wide[iIndex][jIndex]->Write();
    }
    profTimeEnergyClusterCrystalE[jIndex]->Write();
    profEnergyPosition[jIndex]->Write();
    profTimePosition[jIndex]->Write();    
  }
  for( int iIndex = 0; iIndex < nCrystalE-1; iIndex++){
    profTimeEnergyCrystalEMerge[iIndex]->Write();
  }
  profTimeEnergyLowEnergy->Write();
  tfout->Close();
  return 0;
}
