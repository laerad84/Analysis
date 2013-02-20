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
#include "TChain.h"

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

int main( int argc, char** argv) {
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");

  std::string iFileForm = "%s/Cluster_Time_%dMeV_%ddeg-1E5-%d.root";//ROOTFILE_GAMMACLUE
  std::string oFileForm = "ClusterTimeStructure_SIM_All.root";

  //std::string iFileForm = "%s/Data_All_Cal_FNL_COS_HTADJ.root";//ROOTFILE_GAMMACLUS
  //std::string oFileForm = "ClusterTimeStructure_Data_FNL_COS_TimeOffsetSearch.root";
  
  TF1* TimeAdjFuncEnergy = new TF1("TimeAdjFuncEnergy",AdjFunc,0,2000,5);
  //TimeAdjFuncEnergy->SetParameters(-7.51860e-01,9.57348e-01,-9.55972e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-9.69218e-01,1.12202e+00,-7.52470e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-6.32987e-01,7.64508e-01,-9.74328e-02,0,0);//E:200-300;
  //TimeAdjFuncEnergy->SetParameters(-6.38508e-01,7.55330e-01,-9.05449e-02,0,0);//E:200-300,T:10-14
  TimeAdjFuncEnergy->SetParameters(0,0,0,0,0);

  TF1* TimeAdjFuncHeight = new TF1("TimeAdjFuncHeight",AdjfuncTH,0,16000,3);
  //TimeAdjFuncHeight->SetParameters(3.308,0.005644,0.0004385);
  TimeAdjFuncHeight->SetParameters(0,0,0);

  std::cout<< "Chain" << std::endl;
  TChain* tr = new TChain("trCluster");
  std::cout<< __LINE__ <<std::endl;

  /*
  double EArr[7]  = {100, 145, 210, 310, 450, 650, 950};
  double TArr[14] = {10, 12, 14, 16, 18, 20, 22, 25, 28, 30, 32, 35, 40, 50};
  for( int iarr  =0; iarr < 7; iarr++){
    for( int jarr = 0; jarr < 15; jarr++){
      for( int i = 0; i < 4; i++){
	tr->Add(Form(iFileForm.c_str(),EArr[iarr],TArr[jarr],i));
      }
    }
  }
  */
  Int_t SIMTheta = 40; 
  for( int i = 0; i< 4; i++){
    tr->Add(Form(iFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str(),450,SIMTheta,i));
  }
  //TFile* tf       = new TFile(Form(iFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str()));
  //TTree* tr       = (TTree*)tf->Get("trCluster");
  ClusterTimeReader* reader = new ClusterTimeReader(tr);
  Int_t  nEntries = reader->fChain->GetEntries();
  TFile* tfout    = new TFile(oFileForm.c_str(),"recreate");
  
  std::cout<< __LINE__ <<std::endl;

  TH2D*  hisThetaRDistribution = new TH2D("hisThetaRDistribution" ,"hisThetaRDistribution;R[mm];Theta[degree]"    ,300,200,500,60,0,60);
  TH2D*  hisERDistribution     = new TH2D("hisREDistribution"     ,"hisREDistributionl;R[mm];E[MeV]"              ,300,200,500,300,0,900);
  TH2D*  hisEThetaDistribution = new TH2D("hisRThetaDistribution" ,"hisRThetaDistribution;Theta[deg];E[MeV]"      ,60,0,60,300,0,900);
  TH1D*  hisThetaDistribution  = new TH1D("hisThetaDistribution"  ,"hisThetaDistribution;Theta[deg];N/1[deg]"     ,60,0,60);
  TH1D*  hisRDistribution      = new TH1D("hisRDistribution"      ,"hisRDistribution;R[mm];N/1[mm]"               ,300,200,500);
  TH1D*  hisEDistribution      = new TH1D("hisEDistribution"      ,"hisEDistribution;E[MeV];N/3Mev[mm]"           ,300,0,900);
  TH1D*  hisClusterLength      = new TH1D("hisClusterLength"      ,"hisClusterLength;Length[mm];N/5[mm]"          ,24,0,300);
  TH1D*  hisClusterWidth       = new TH1D("hisClusterWidth"       ,"hisClusterWidth;Width[mm];N/5[mm]"            ,24,0,300);
  TH2D*  hisClusterEnergyLength= new TH2D("hisClusterEnergyLength","hisClusterEnergyLength;Energy[MeV];Length[mm]",300,0,900,24,0,300);
  TH2D*  hisClusterEnergyWidth = new TH2D("hisClusterEnergyWidth" ,"hisClusterEnergyWidth;Energy[MeV];Length[mm]" ,300,0,900,24,0,300);
  TProfile2D* profEETDependency= new TProfile2D("profEETDependency","profEETDependency;E[MeV];E[MeV]",200,0,600,200,0,600);
  TProfile*   profETDependency = new TProfile("profETDependency"  ,"profETDependency;E[MeV];T[ns]",300,0,900);


  const int nRegion= 6;  
  TProfile*   profET[nRegion];
  TProfile*   profETR[3][nRegion];
  TProfile*   profCET[nRegion];
  TProfile*   profCETR[3][nRegion];
  TProfile*   profCRT[nRegion];
  TProfile*   profCER[nRegion];
  TProfile*   profCRE[nRegion];

  TProfile*   profTR_E_restrict[nRegion];
  TH2D*       hisTR_E_restrict[nRegion];
  TProfile*   profTD_E_restrict[nRegion];
  TH2D*       hisTD_E_restrict[nRegion];
  TProfile2D* prof2ndET[nRegion];
  TProfile* prof2ndETR[nRegion];
  TProfile* prof2ndETD[nRegion];
  TH2D*     his2ndETR[nRegion];
  TH2D*     his2ndETD[nRegion];

  Int_t EnergyBoundary[nRegion+1] = {100,200,300,500,800,1300,2000};
  Int_t CEnergyBoundary[nRegion+1]= {50,100,150,250,400,650,1000};
  const int nTRegion  = 6;
  Int_t ThetaBoundary[nTRegion+1] = {0,10,14,18,26,40,90};
  
  std::cout<< "Setting" << std::endl;
  for( int Index = 0; Index < nRegion; Index++){
    
    profET[Index] = new TProfile(Form("profET_%d",Index),Form("profET_%d_%d",EnergyBoundary[Index],EnergyBoundary[Index+1]),300,0,300);
    profCET[Index] = new TProfile(Form("profCET_%d",Index),Form("profCET_%d_%d",CEnergyBoundary[Index],CEnergyBoundary[Index+1]),300,0,600);
    profCER[Index] = new TProfile(Form("profCER_%d",Index),Form("profCER_%d_%d",CEnergyBoundary[Index],CEnergyBoundary[Index+1]),300,0,600);
    profCRE[Index] = new TProfile(Form("profCRE_%d",Index),Form("profCRE_%d_%d",CEnergyBoundary[Index],CEnergyBoundary[Index+1]),25,0,25);
    profCRT[Index] = new TProfile(Form("profCRT_%d",Index),Form("profCRT_%d_%d",CEnergyBoundary[Index],CEnergyBoundary[Index+1]),25,0,25);
    
    profTR_E_restrict[Index] = new TProfile(Form("profTR_E_restrict_%d",Index),Form("profTR_E_restrict_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105);
    hisTR_E_restrict[Index]  = new TH2D(Form("hisTR_E_restrict_%d",Index),Form("hisTR_E_restrict_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105,100,-10,10);
    profTD_E_restrict[Index] = new TProfile(Form("profTD_E_restrict_%d",Index),Form("profTD_E_restrict_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105);
    hisTD_E_restrict[Index]  = new TH2D(Form("hisTD_E_restrict_%d",Index),Form("hisTD_E_restrict_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105,100,-10,10);
    prof2ndET[Index]         = new TProfile2D(Form("prof2ndET_%d",Index),Form("prof2ndET_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105,21,-105,105);
    prof2ndETR[Index]        = new TProfile(Form("prof2ndETR_%d",Index),Form("prof2ndETR_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105);
    prof2ndETD[Index]        = new TProfile(Form("prof2ndETD_%d",Index),Form("prof2ndETD_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105);
    his2ndETR[Index]         = new TH2D(Form("his2ndETR_%d",Index),Form("his2ndETR_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105,100,-10,10);
    his2ndETD[Index]         = new TH2D(Form("his2ndETD_%d",Index),Form("his2ndETD_%d_%d",ThetaBoundary[Index],ThetaBoundary[Index+1]),21,-105,105,100,-10,10);

    for( int IIndex = 0; IIndex < 3; IIndex++){
      profETR[IIndex][Index] = new TProfile(Form("profETR_%d_%d",IIndex,Index),Form("profETR_%d_%d_%d",IIndex,EnergyBoundary[Index],EnergyBoundary[Index+1]),50,0,300);
      profCETR[IIndex][Index] = new TProfile(Form("profCETR_%d_%d",IIndex,Index),Form("proCfCETR_%d_%d_%d",IIndex,CEnergyBoundary[Index],CEnergyBoundary[Index+1]),50,0,300);
    }
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Cut Condition 
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  const Double_t InnerRadCut  = 200;
  const Double_t OuterRadCut  = 500;
  const Double_t LowEnergyCut = 3;
  
  std::cout<< "Loop" << std::endl;
  for( int ievent = 0; ievent < nEntries; ievent++){
    reader->GetEntry(ievent);
    
    if( reader->nCluster > 1 ){ continue; }
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      Double_t clusterR = reader->ClusterR[clusterIndex];
      Double_t clusterTheta = reader->ClusterTheta[clusterIndex];
      Double_t clusterPhi = reader->ClusterPhi[clusterIndex];

      if( reader->ClusterChisq2[clusterIndex] > 2.5 ){ continue; }
      if( reader->ClusterR[clusterIndex] > OuterRadCut ){ continue; }
      //if( reader->ClusterR[clusterIndex] < InnerRadCut ){ continue; }
      if( TMath::Abs(clusterR*TMath::Cos(clusterPhi)) < 150 && TMath::Abs( clusterR*TMath::Sin(clusterPhi)) < 150 ){ continue; }
      if( reader->ClusterEnergy[clusterIndex] < 100 ){ continue; }

      hisThetaDistribution->Fill(reader->ClusterTheta[clusterIndex]*180/TMath::Pi());
      hisRDistribution    ->Fill(reader->ClusterR[clusterIndex]);
      hisEDistribution    ->Fill(reader->ClusterEnergy[clusterIndex]);
      hisThetaRDistribution->Fill(reader->ClusterR[clusterIndex]    ,reader->ClusterTheta[clusterIndex]*180/TMath::Pi());
      hisEThetaDistribution->Fill(reader->ClusterTheta[clusterIndex]*180./TMath::Pi(),reader->ClusterEnergy[clusterIndex]);
      hisERDistribution    ->Fill(reader->ClusterR[clusterIndex]    ,reader->ClusterEnergy[clusterIndex] );
      
      Int_t EnergyIndex = -1;
      for( int Index  = 0; Index < nRegion+1; Index++){
	if( reader->ClusterEnergy[clusterIndex] > EnergyBoundary[Index] && reader->ClusterEnergy[clusterIndex] <= EnergyBoundary[Index+1] ){
	  EnergyIndex  = Index;
	  break;
	}
      }

      Int_t ThetaIndex = -1;
      /*
      for( int Index = 0; Index <nRegion+1; Index++){
	if( reader->ClusterTheta[clusterIndex]*180 > ThetaBoundary[Index]*TMath::Pi() && reader->ClusterTheta[clusterIndex]*180 <= ThetaBoundary[Index+1]*TMath::Pi()){
	  ThetaIndex = Index;
	  break;
	}
      }
      */
      for( int Index = 0; Index <nRegion+1; Index++){
	if( SIMTheta > ThetaBoundary[Index]*TMath::Pi() && SIMTheta <= ThetaBoundary[Index+1]*TMath::Pi()){
	  ThetaIndex = Index;
	  break;
	}
      }
      EnergyIndex = 0;
      ThetaIndex  = 0;

      //std::cout<< ThetaIndex << std::endl;
      //if(TMath::Tan(reader->ClusterTheta[clusterIndex]) > 0.3 ){ continue; }
      //if( reader->ClusterEnergy[clusterIndex] > 300 ){ continue; }
      //if( reader->ClusterEnergy[clusterIndex] < 200 ){ continue; }
      //if( reader->ClusterTheta[clusterIndex]*180/TMath::Pi() < 28  ){ continue; }
      //if( reader->ClusterTheta[clusterIndex]*180/TMath::Pi() >= 28 ){ continue; }

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
	Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TimeAdjFuncHeight->Eval(SiginCluster)-TimeAdjFuncEnergy->Eval(EinCluster);
	//Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TimeAdjFuncHeight->Eval(SiginCluster);


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
      Int_t CEnergyIndex = -1;
      for( int Index  = 0; Index < nRegion+1; Index++){
	if( ECenterCrystal > CEnergyBoundary[Index] && ECenterCrystal <= CEnergyBoundary[Index+1] ){
	  CEnergyIndex  = Index;
	}
      }
      

      if( ECenterCrystal < 0.4*reader->ClusterEnergy[clusterIndex] ){ continue; }
      if( ECenterCrystal > 500 ){ continue; }
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

	Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t PhiinCluster = reader->CrystalPhi[clusterIndex][crystalIndex];	
	Double_t RinCluster   = RadinCluster*TMath::Cos(PhiinCluster);
	Double_t DinCluster   = RadinCluster*TMath::Sin(PhiinCluster);
	Double_t SiginCluster = reader->CrystalSignal[clusterIndex][crystalIndex];
	Double_t EinCluster   = reader->CrystalEnergy[clusterIndex][crystalIndex];
	Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TimeAdjFuncHeight->Eval(SiginCluster)-TimeAdjFuncEnergy->Eval(EinCluster)-TCenterCrystal;//-0.05*RadinCluster;
	if( TinCluster > 10  || TinCluster < -10 ){ continue;}
	//Double_t TinCluster   = reader->CrystalT[clusterIndex][crystalIndex]-TCenterCrystal-TimeAdjFuncHeight->Eval(SiginCluster);
	if( EnergyIndex < 0 || EnergyIndex >= nRegion ){ continue; }
	if( RadinCluster <= RadCenterCrystal ){ continue; }
	profET[EnergyIndex]->Fill(EinCluster,TinCluster);

	if(ECenterCrystal > 100 && ECenterCrystal < 500){
	  if( EinCluster > 50 ){
	    prof2ndET[ThetaIndex]->Fill( RinCluster,DinCluster,TinCluster);
	    if( TMath::Abs(DinCluster) < 25 ){ 
	      prof2ndETR[ThetaIndex]->Fill( RinCluster,TinCluster);
	      his2ndETR[ThetaIndex]->Fill( RinCluster,TinCluster);
	    }
	    if( TMath::Abs(RinCluster) < 25 ){ 
	      prof2ndETD[ThetaIndex]->Fill( TMath::Abs(DinCluster),TinCluster);
	      his2ndETD[ThetaIndex]->Fill( TMath::Abs( DinCluster ), TinCluster);
	    }
	  }
	}
	
	if( ThetaIndex >= 0 && ThetaIndex < nRegion ){
	  if( ECenterCrystal > 200 && ECenterCrystal < 500 ){
	    if( EinCluster > 3 && EinCluster < 24 ){
	      if( TMath::Abs(DinCluster) < 62.5){
		profTR_E_restrict[ThetaIndex]->Fill(RinCluster,TinCluster);
		hisTR_E_restrict[ThetaIndex]->Fill(RinCluster,TinCluster);
	      }
	      if( TMath::Abs(RinCluster) < 62.5){
		profTD_E_restrict[ThetaIndex]->Fill(DinCluster,TinCluster);
		hisTD_E_restrict[ThetaIndex]->Fill(DinCluster,TinCluster);
	      }
	    }
	  }
	}
	if( RadinCluster < 25&& RadCenterCrystal<6.3 ){
	  if( TMath::Abs( PhiinCluster ) < TMath::Pi()/4. ){
	    profETR[0][EnergyIndex]->Fill(EinCluster,TinCluster);
	    if( CEnergyIndex >= 0 )
	      profCETR[0][CEnergyIndex]->Fill(EinCluster,TinCluster);
	  }else if( TMath::Abs(PhiinCluster) < TMath::Pi()*3./4.){
	    profETR[1][EnergyIndex]->Fill(EinCluster,TinCluster);
	    if( CEnergyIndex >= 0 )
	      profCETR[1][CEnergyIndex]->Fill(EinCluster,TinCluster);
	    if( CEnergyIndex >= 0 ){
	      profCET[CEnergyIndex]->Fill(EinCluster,TinCluster);
	      profCER[CEnergyIndex]->Fill(EinCluster,RadinCluster);
	      profCRE[CEnergyIndex]->Fill(RadinCluster,EinCluster);
	      profCRT[CEnergyIndex]->Fill(RadinCluster,TinCluster);
	    }
	  }else{
	    profETR[2][EnergyIndex]->Fill(EinCluster,TinCluster);
	    if( CEnergyIndex >= 0 )
	      profCETR[2][CEnergyIndex]->Fill(EinCluster,TinCluster);

	  }
	}

	//if( DinCluster < 0.5*25 ){ continue; }
	if( abs(RinCluster) > 25 ){ continue; }
	RadVec.push_back(RadinCluster);
	PhiVec.push_back(PhiinCluster);
	RVec.push_back(RinCluster);
	DVec.push_back(DinCluster);
	HVec.push_back(SiginCluster);
	TVec.push_back(TinCluster);
	EVec.push_back(EinCluster);
      }
      //std::cout<< "Test" << std::endl;
      if( RadVec.size() < 2 ){ continue; }
      for( int iIndex = 0; iIndex < RadVec.size(); iIndex++){
	for( int jIndex = iIndex+1; jIndex < RadVec.size(); jIndex++){
	  if(abs(abs(DVec.at(jIndex))-abs(DVec.at(iIndex))) < 0.5*25 ){
	    if(EVec.at(jIndex) > EVec.at(iIndex)){
	      profEETDependency->Fill( EVec.at(jIndex), EVec.at(iIndex), TVec.at(jIndex)-TVec.at(iIndex));
	      if( EVec.at(iIndex) < 9){
	      profETDependency->Fill( EVec.at(jIndex),TVec.at(jIndex)-TVec.at(iIndex));
	      }
	    }else{
	      profEETDependency->Fill( EVec.at(iIndex), EVec.at(jIndex), TVec.at(iIndex)-TVec.at(jIndex));
	      if( EVec.at(jIndex) <9){
	      profETDependency->Fill( EVec.at(iIndex),TVec.at(iIndex)-TVec.at(jIndex));
	      }
	    }
	  }
	}
      }
    }
  }

  for( int Index = 0; Index < nRegion; Index++){
    profET[Index]->Write();
    profCET[Index]->Write();
    profCER[Index]->Write();
    profCRE[Index]->Write();
    profCRT[Index]->Write();
    
    profTR_E_restrict[Index]->Write();
    hisTR_E_restrict[Index]->Write();
    profTD_E_restrict[Index]->Write();
    hisTD_E_restrict[Index]->Write();
    prof2ndET[Index]->Write();
    prof2ndETR[Index]->SetLineColor(Index+1);
    prof2ndETD[Index]->SetLineColor(Index+1);

    prof2ndETR[Index]->Write();
    prof2ndETD[Index]->Write();
    his2ndETR[Index]->Write();
    his2ndETD[Index]->Write();

    for( int IIndex  = 0; IIndex < 3; IIndex++){
      profETR[IIndex][Index]->SetLineColor(IIndex+1);
      profETR[IIndex][Index]->Write();
      profCETR[IIndex][Index]->SetLineColor(IIndex+1);
      profCETR[IIndex][Index]->Write();
    }
  }

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
