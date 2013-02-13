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
#include "TProfile2D.h"
#include "TProfile.h"
double AdjFunc( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double value = p0 + 1.2*p1*exp(p2*x0) + 0*p3*exp(p4*x0);
  return value;
}
double AdjFunc1( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double value = p0 + p1*exp(p2*x0);
  return value;
}
double AdjFuncTest( double* x, double* par){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double value = p0 + p1*exp(p2*x0) + p3*exp(p4*x0);
  return value;

}
double AdjfuncTH( double*x ,double*par ){  
  double value = par[0]*(-1*TMath::Log(1 + par[1])+TMath::Log( 1 + par[1]*exp( par[2]*x[0] )));
  return value;
}

int main( int argc, char** argv ){
  int Energy = atoi(argv[1]);
  int Theta  = atoi(argv[2]);

  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string iFileForm = "%s/Cluster_Time_%dMeV_%ddeg-1E5-%d.root";
  std::string oFileForm = "ClusterTimeStructure_%dMeV_%ddeg.root";
  

  TF1* TimeAdjFuncEnergy = new TF1("TimeAdjFuncEnergy",AdjFunc,0,2000,5);
  //TimeAdjFuncEnergy->SetParameters(-7.51860e-01,9.57348e-01,-9.55972e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-9.69218e-01,1.12202e+00,-7.52470e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(3.10325e-01,-4.18253e-01,-1.65285e-02,0,0);
  TF1* TimeAdjFuncHeight = new TF1("TimeAdjFuncHeight",AdjfuncTH,0,16000,3);
  TimeAdjFuncHeight->SetParameters(3.308,0.005644,0.0004385);
  
  TFile* tf       = new TFile(Form(iFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str(),Energy,Theta,5));
  TTree* tr       = (TTree*)tf->Get("trCluster");
  ClusterTimeReader* reader = new ClusterTimeReader(tr);
  Int_t  nEntries = reader->fChain->GetEntries();
  TFile* tfout    = new TFile(Form(oFileForm.c_str(),Energy,Theta),"recreate");
  

  TProfile2D* profRDT = new TProfile2D("profRDT","profRDT",
				   21,-105,105,21,-105,105);
  TProfile*   profRT  = new TProfile("profRT","profRT",
				     21,-105,105);
  TProfile*   profDT  = new TProfile("profDT","profDT",
				     21,-105,105);
  
  
  TH2D*  hisThetaRDistribution = new TH2D("hisThetaRDistribution",
					  "hisThetaRDistribution;R[mm];Theta[degree]",
					  300,200,500,60,0,60);
  TH2D*  hisERDistribution     = new TH2D("hisREDistribution",
					  "hisREDistributionl;R[mm];E[MeV]",
					  300,200,500,300,0,900);
  TH2D*  hisEThetaDistribution = new TH2D("hisRThetaDistribution",
					  "hisRThetaDistribution;Theta[deg];E[MeV]",
					  60,0,60,300,0,900);
  TH1D*  hisThetaDistribution  = new TH1D("hisThetaDistribution",
					  "hisThetaDistribution;Theta[deg];N/1[deg]",
					  60,0,60);
  TH1D*  hisRDistribution      = new TH1D("hisRDistribution",
					  "hisRDistribution;R[mm];N/1[mm]",
					  300,200,500);
  TH1D*  hisEDistribution      = new TH1D("hisEDistribution",
					  "hisEDistribution;E[MeV];N/3Mev[mm]",
					  300,0,900);
  TH1D*  hisClusterLength       = new TH1D("hisClusterLength",
					  "hisClusterLength;Length[mm];N/5[mm]",
					  24,0,300);
  TH1D*  hisClusterWidth       = new TH1D("hisClusterWidth",
					  "hisClusterWidth;Width[mm];N/5[mm]",
					  24,0,300);
  TH2D*  hisClusterEnergyLength= new TH2D("hisClusterEnergyLength",
					  "hisClusterEnergyLength;Energy[MeV];Length[mm]",
					  300,0,900,24,0,300);
  TH2D*  hisClusterEnergyWidth = new TH2D("hisClusterEnergyWidth",
					  "hisClusterEnergyWidth;Energy[MeV];Length[mm]",
					  300,0,900,24,0,300);

  TProfile2D* profRDE = new TProfile2D("profRDE","profRDE",
				       21,-105,105,21,-105,105);
  TProfile2D* profEETDependency = new TProfile2D("profEETDependency",
						"profEETDependency;E[MeV];E[MeV]",
						200,0,600,200,0,600);
  TProfile*   profETDependency = new TProfile("profETDependency",
					      "profETDependency;E[MeV];T[ns]",
					      300,0,900);
  
  const Double_t InnerRadCut  = 200;
  const Double_t OuterRadCut  = 500;
  const Double_t LowEnergyCut = 3;
  std::cout<< nEntries << std::endl;
  for( int ievent = 0; ievent < nEntries; ievent++){
    reader->GetEntry(ievent);
    
    if( reader->nCluster >= 2 ) { continue; }
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      //if( reader->ClusterChisq2[clusterIndex] > 2.5 ){ continue; }
      if( reader->ClusterR[clusterIndex] > OuterRadCut ){ continue; }
      if( reader->ClusterR[clusterIndex] < InnerRadCut ){ continue; }
      if( reader->ClusterEnergy[clusterIndex] < 100 ){ continue; }
      hisThetaDistribution->Fill(reader->ClusterTheta[clusterIndex]*180/TMath::Pi());
      hisRDistribution    ->Fill(reader->ClusterR[clusterIndex]);
      hisEDistribution    ->Fill(reader->ClusterEnergy[clusterIndex]);
      hisThetaRDistribution->Fill(reader->ClusterR[clusterIndex]    ,reader->ClusterTheta[clusterIndex]*180/TMath::Pi());
      hisEThetaDistribution->Fill(reader->ClusterTheta[clusterIndex]*180./TMath::Pi(),reader->ClusterEnergy[clusterIndex]);
      hisERDistribution    ->Fill(reader->ClusterR[clusterIndex]    ,reader->ClusterEnergy[clusterIndex] );
      
      //if(TMath::Tan(reader->ClusterTheta[clusterIndex]) > 0.3 ){ continue; }
      //if( reader->ClusterEnergy[clusterIndex] > 500 ){ continue; }
      //if( reader->ClusterEnergy[clusterIndex] < 300 ){ continue; }
      Double_t RadCenterCrystal = 1000;
      Double_t RCenterCrystal = 0;
      Double_t DCenterCrystal = 0;
      Double_t TCenterCrystal = 0;
      Double_t ECenterCrystal = 0;
      Double_t EMaxCrystal    = 0;
      Double_t ClusterLength  = 0;
      Double_t ClusterWidth   = 0;
      Double_t DMinCrystal    = 1000;
      Double_t DMaxCrystal    = -1000;
      Double_t RMinCrystal    = 1000;
      Double_t RMaxCrystal    = -1000;

      for(int crystalIndex = 0; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	Double_t EinCluster = reader->CrystalEnergy[clusterIndex][crystalIndex];
	if( EinCluster > EMaxCrystal ){
	  EMaxCrystal = EinCluster; 
	}
      }
      
      
      for(int crystalIndex = 0 ;crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t PhiinCluster = reader->CrystalPhi[clusterIndex][crystalIndex];	
	Double_t RinCluster   = RadinCluster*TMath::Cos(PhiinCluster);
	Double_t DinCluster   = RadinCluster*TMath::Sin(PhiinCluster);
	Double_t SiginCluster = reader->CrystalSignal[clusterIndex][crystalIndex];
	Double_t EinCluster   = reader->CrystalEnergy[clusterIndex][crystalIndex];
	//Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TimeAdjFuncHeight->Eval(SiginCluster))-TimeAdjFuncEnergy->Eval(EinCluster);
	Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex];//-TimeAdjFuncEnergy->Eval(EinCluster);//-TimeAdjFuncHeight->Eval(SiginCluster);
	
	if( RinCluster > RMaxCrystal ){ RMaxCrystal = RinCluster;}
	if( RinCluster < RMinCrystal ){ RMinCrystal = RinCluster;}
	if( DinCluster > DMaxCrystal ){ DMaxCrystal = DinCluster;}
	if( DinCluster < DMinCrystal ){ DMinCrystal = DinCluster;}

	if( RadinCluster < RadCenterCrystal ){
	  RadCenterCrystal= RadinCluster;
	  RCenterCrystal = RinCluster;
	  DCenterCrystal = DinCluster;
	  ECenterCrystal = EinCluster;
	  TCenterCrystal = TinCluster;
	}
      }
      //if( ECenterCrystal < 0.3*reader->ClusterEnergy[clusterIndex] ){ continue; }

      hisClusterWidth->Fill(DMaxCrystal-DMinCrystal);
      hisClusterLength->Fill(RMaxCrystal-RMinCrystal);
      hisClusterEnergyWidth->Fill(reader->ClusterEnergy[clusterIndex],DMaxCrystal-DMinCrystal);
      hisClusterEnergyLength->Fill(reader->ClusterEnergy[clusterIndex],RMaxCrystal-RMinCrystal);

      std::vector<double> RadVec;
      std::vector<double> PhiVec;
      std::vector<double> RVec;
      std::vector<double> DVec;
      std::vector<double> TVec;
      std::vector<double> EVec;
      std::vector<double> HVec;
      std::vector<double>::iterator itRadVec;
      std::vector<double>::iterator itPhiVec;
      std::vector<double>::iterator itRVec;
      std::vector<double>::iterator itDVec;
      std::vector<double>::iterator itEVec;
      std::vector<double>::iterator itHVec;

      for(int crystalIndex = 0 ; crystalIndex < reader->nCrystal[clusterIndex]; crystalIndex++){
	if( reader->nCrystal[clusterIndex] < 5 ){ continue; }

	Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t PhiinCluster = reader->CrystalPhi[clusterIndex][crystalIndex];	
	Double_t RinCluster   = RadinCluster*TMath::Cos(PhiinCluster);
	Double_t DinCluster   = RadinCluster*TMath::Sin(PhiinCluster);
	Double_t SiginCluster = reader->CrystalSignal[clusterIndex][crystalIndex];
	Double_t EinCluster   = reader->CrystalEnergy[clusterIndex][crystalIndex];
	//Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TimeAdjFuncHeight->Eval(SiginCluster)-TimeAdjFuncEnergy->Eval(EinCluster)-TCenterCrystal;
	Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TCenterCrystal;//-TimeAdjFuncEnergy->Eval(EinCluster);//-TimeAdjFuncHeight->Eval(SiginCluster);

	if( TMath::IsNaN(TinCluster)) { continue; }
	if( TMath::IsNaN(RinCluster)) { continue; }
	if( TMath::IsNaN(DinCluster)) { continue; }
	if( TMath::IsNaN(RadinCluster)){continue;}
	if( TMath::IsNaN(PhiinCluster)){continue;}
	//if( DinCluster < 0.5*25 ){ continue; }
	profRDE->Fill(RinCluster,DinCluster,EinCluster);
	profRDT->Fill(RinCluster,DinCluster,TinCluster);
	if( TMath::Abs(DinCluster) < 25 ){
	  profRT->Fill(RinCluster,TinCluster);
	}
	if( TMath::Abs(RinCluster) < 25 ){
	  profDT->Fill(DinCluster,TinCluster);
	}
	
	if( abs(RinCluster) > 25 ){ continue; }

	RadVec.push_back(RadinCluster);
	PhiVec.push_back(PhiinCluster);
	RVec.push_back(RinCluster);
	DVec.push_back(DinCluster);
	HVec.push_back(SiginCluster);
	TVec.push_back(TinCluster);
	EVec.push_back(EinCluster);
      }

      if( RadVec.size() < 2 ){ continue; }
      for( int iIndex = 0; iIndex < RadVec.size(); iIndex++){
	for( int jIndex = iIndex+1; jIndex < RadVec.size(); jIndex++){
	  if(abs(abs(DVec.at(jIndex))-abs(DVec.at(iIndex))) < 0.5*25 ){
	    if(EVec.at(jIndex) > EVec.at(iIndex)){
	      profEETDependency->Fill( EVec.at(jIndex), EVec.at(iIndex), TVec.at(jIndex)-TVec.at(iIndex));
	      if( EVec.at(iIndex) < 9 ){ 
		profETDependency->Fill( EVec.at(jIndex),TVec.at(jIndex)-TVec.at(iIndex));
	      }
	    }else{
	      profEETDependency->Fill( EVec.at(iIndex), EVec.at(jIndex), TVec.at(iIndex)-TVec.at(jIndex));
	      if( EVec.at(jIndex) < 9 ){ 
		profETDependency->Fill( EVec.at(iIndex),TVec.at(iIndex)-TVec.at(jIndex));
	      }
	    }
	  }
	}
      }
    }
  }

  profRDE->Write();
  profRDT->Write();
  profRT->Write();
  profDT->Write();
  hisThetaDistribution ->Write();
  hisRDistribution     ->Write();
  hisEDistribution     ->Write();
  hisThetaRDistribution->Write();
  hisEThetaDistribution->Write();
  hisERDistribution    ->Write();
  hisClusterWidth->Write();
  hisClusterLength->Write();
  hisClusterEnergyWidth->Write();
  hisClusterEnergyLength->Write();
  profEETDependency->Write();
  profETDependency->Write();
  tfout->Write();
}

