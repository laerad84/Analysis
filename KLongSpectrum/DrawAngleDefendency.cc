#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include <iostream>
#include "GeneralFunctions.h"
int main( int argc, char** argv){


  TFile* tf = new TFile("kl_GammaTimeShape_DATA_NONTIMECAL.root");
  TTree* tr = (TTree*)tf->Get("GammaTimeShape");

  Double_t Radius;//Distance from Csi Center
  Double_t ZVtx;//Distance from Csi surface
  Double_t theta;//InjectionAngle
  //Double_t cot;//cotangent of Injection Angle Z/R
  Double_t phi;//radial angle
  Double_t X;
  Double_t Y;
  Int_t    ClusterSize;
  Double_t E[120];//Energy;
  Double_t T[120];//Timing
  Double_t R[120];//Radial distance
  Double_t D[120];//
  Double_t FractionAngle[120];
  Int_t    CutCondition;
  Int_t    EventID;
  Double_t BaseTime;
  Double_t EGamma;
  Double_t ECenter;
  Double_t GammaChisq;
  tr->SetBranchAddress("EGamma",&EGamma);
  tr->SetBranchAddress("ECenter",&ECenter);
  tr->SetBranchAddress("EventID",&EventID);
  tr->SetBranchAddress("Radius",&Radius);
  tr->SetBranchAddress("X",&X);
  tr->SetBranchAddress("Y",&Y);
  tr->SetBranchAddress("ZVtx",&ZVtx);
  tr->SetBranchAddress("theta",&theta);
  tr->SetBranchAddress("phi",&phi);
  tr->SetBranchAddress("ClusterSize",&ClusterSize);
  tr->SetBranchAddress("BaseTime",&BaseTime);
  tr->SetBranchAddress("E",E);//ClusterSize
  tr->SetBranchAddress("T",T);//ClusterSize
  tr->SetBranchAddress("R",R);//ClusterSize
  tr->SetBranchAddress("D",D);//ClusterSize
  tr->SetBranchAddress("FractionAngle",FractionAngle);//ClusterSize
  tr->SetBranchAddress("CutCondition",&CutCondition);
  tr->SetBranchAddress("GammaChisq",&GammaChisq);

  TFile* tfOut = new TFile("GammaTimeShape.root","recreate");

  int nAngle = 9;
  Double_t *AngleDist = GenLogArray(nAngle+1,1,40);
  const int      nDiv      = 4;
  const int      nDivTime  = 200;
  TH1D* hisAngle = new TH1D("hisAngle","hisAngle",9,AngleDist);
  TH1D* hisDist_D = new TH1D("hisDist_D","hisDist_D",nDiv+1,-12.5,112.5);

  TH2D* hisRT[nAngle];
  TH2D* hisDT[nAngle];
  TH2D* hisRDT[nAngle][nDiv+1];
  Double_t DivLength = 25;
  for( int i = 0; i< nAngle; i++){
    for( int j = 0; j < nDiv+1;j++){ 
      hisRDT[i][j] = new TH2D(Form("hisRDT_%d_%d",i,j),Form("hisRDT_%d_%d",i,j),
			  2*nDiv+1,-(nDiv +0.5)*DivLength,(nDiv+0.5)*DivLength,
			  200,-10,10);			      
    }
    hisRT[i] = new TH2D(Form("hisRT_%d",i),Form("hisRT_%d",i),
			2*nDiv+1,-(nDiv +0.5)*DivLength,(nDiv+0.5)*DivLength,
			200,-10,10);
    hisDT[i] = new TH2D(Form("hisDT_%d",i),Form("hisDT_%d",i),
			2*nDiv+1,-(nDiv +0.5)*DivLength,(nDiv+0.5)*DivLength,
			200,-10,10);
  }


  for( int ievent = 0; ievent< tr->GetEntries(); ievent++){
    tr->GetEntry(ievent);    
    if( (ievent %1000) == 0 ){ std::cout<< ievent << "/" << tr->GetEntries() << std::endl;}
    Int_t AngleIndex = hisAngle->Fill(ZVtx/Radius)-1;
    if( AngleIndex >= nAngle ){ continue; }
    if( AngleIndex <0 ){ continue; }

    if( CutCondition != 0 ){ continue; }
    if( EGamma < 100 ){ continue; }
    for( int i = 0; i< ClusterSize; i++){
      if( E[i] < 10 ){ continue;}
      if( TMath::Abs(D[i]) < 12.5 ){
	hisRT[AngleIndex]->Fill(R[i],T[i]-BaseTime);
      }
      if( TMath::Abs(R[i]) < 12.5 ){
	hisDT[AngleIndex]->Fill(D[i],T[i]-BaseTime);
      }
      int DIndex = hisDist_D->Fill(TMath::Abs(D[i]))-1;
      hisRDT[AngleIndex][DIndex]->Fill(R[i],T[i]-BaseTime);
    }
  }
  
  for( int i = 0; i< nAngle; i++){
    hisRT[i]->Write();
    hisDT[i]->Write();
    for( int j = 0; j< nDiv+1; j++){
      hisRDT[i][j]->Write();
    }
  }

  tfOut->Close();
  return 0; 
}
