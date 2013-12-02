#include <iostream>
#include <string>
#include <fstream>
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TChain.h"
#include "TMath.h"
int main( int argc, char** argv){
  TFile* tfTemp = new TFile("~/local/Analysis/AnalysisLib/Data/Temperature_Factor/TemperatureCorrectionFactor.root");
  TTree* trTemp = (TTree*)tfTemp->Get("TemperatureCorrectionCsI");
  int  RN;
  double temp;  
  trTemp->SetBranchAddress("RunNumber",&RN);
  trTemp->SetBranchAddress("Temperature",&temp);

  TFile* tfLaserLin = new TFile("~/local/Analysis/AnalysisLib/Data/HeightLinearity_Laser.root");
  TGraph* grLin[3];
  for( int i = 0; i< 3; i++){
    grLin[i] = (TGraph*)tfLaserLin->Get(Form("heightLinearity_%d",i));
  }


  const int nAll = 495;
  TFile* tf[nAll];
  std::vector<int> RunNumberVec;


  TChain* ch = new TChain("trLaser");
  
  std::ifstream ifs("../RunList/RunAll.csv");
  int tmpRunNumber;
  while( ifs >> tmpRunNumber ){
    RunNumberVec.push_back( tmpRunNumber );
    ch->Add(Form("Data/LaserHeightTiming_%d.root",tmpRunNumber));
  }
  Int_t RunNumber;
  Double_t Entries[2716];
  Double_t Output[2716];
  Double_t Mean[2716];
  Double_t RMS[2716];
  Double_t Error[2716];
  Double_t Chisq[2716];

  for( int i = 0; i< 2716; i++){

    Entries[i] = 0;
    Output[i] = 0;
    RMS[i] = 0;
    Error[i] = 0;
    Chisq[i] = 0;
  }

  ch->SetBranchAddress("RunNumber",&RunNumber);
  ch->SetBranchAddress("Entries",Entries);
  ch->SetBranchAddress("Output",Output);
  ch->SetBranchAddress("Error",Error);
  ch->SetBranchAddress("Chisq",Chisq);
  ch->SetBranchAddress("Mean",Mean);
  ch->SetBranchAddress("RMS",RMS);


  TFile* tfOut = new TFile("Output.root","recreate");

  TGraphErrors* gr[2716];
  TGraphErrors* grAdj[2716];
  TGraphErrors* grBase[2716];
  TGraphErrors* grBaseAdj[2716];
  TGraphErrors* grBaseHeight[2716];
  TGraphErrors* grBaseTemp[2716];
  for( int i = 0; i< 2716; i++){
    gr[i] = new TGraphErrors();
    gr[i]->SetNameTitle(Form("LaserOutput_%d",i),Form("LaserOutput_%d",i));
    grAdj[i] = new TGraphErrors();
    grAdj[i]->SetNameTitle(Form("LaserOutputAdj_%d",i),Form("LaserOutputAdj_%d",i));
    grBase[i] = new TGraphErrors();
    grBase[i]->SetNameTitle(Form("LaserOutputBase_%d",i),Form("LaserOutputBase_%d",i));
    grBaseAdj[i] = new TGraphErrors();
    grBaseAdj[i]->SetNameTitle(Form("LaserOutputBaseAdj_%d",i),Form("LaserOutputBaseAdj_%d",i));
    grBaseHeight[i] = new TGraphErrors();
    grBaseHeight[i]->SetNameTitle(Form("LaserOutputBaseHeight_%d",i),Form("LaserOutputBaseHeight_%d",i));
    grBaseTemp[i] = new TGraphErrors();
    grBaseTemp[i]->SetNameTitle(Form("LaserOutputBaseTemp_%d",i),Form("LaserOutputBaseTemp_%d",i));

  }
  TH2D* hisRatio = new TH2D("hisRatio","hisRatio",160,0,16000,200,0,10);
  TH2D* hisRatioAdj = new TH2D("hisRatioAdj","hisRatioAdj",160,0,16000,200,0,10);
  for( int ientries = 0 ; ientries < ch->GetEntries(); ientries++){
    ch->GetEntry( ientries );
    trTemp->GetEntry( RunNumber );
    Double_t BaseOut = Output[22];
    for( int i = 0; i< 2716; i++){
      int CorrectionID = 0;
      if( i < 2240 ){ CorrectionID = 0;}
      else{
	CorrectionID = 2; 
      }
      if( Chisq < 0 || TMath::IsNaN( Chisq[i] )){ 
	gr[i]->SetPoint( gr[i]->GetN(), RunNumber,0);
	gr[i]->SetPointError( gr[i]->GetN()-1, 0,0);
	grAdj[i]->SetPoint( grBaseAdj[i]->GetN(), RunNumber, 0);
	grBase[i]->SetPoint( grBaseAdj[i]->GetN(), RunNumber, 0);
	grBaseAdj[i]->SetPoint( grBaseAdj[i]->GetN(), RunNumber, 0);
	grAdj[i]->SetPointError( grBaseAdj[i]->GetN()-1,0,0);
	grBase[i]->SetPointError( grBaseAdj[i]->GetN()-1,0,0);
	grBaseAdj[i]->SetPointError( grBaseAdj[i]->GetN()-1,0,0);
      }else{
	gr[i]->SetPoint( gr[i]->GetN(), RunNumber,Output[i]);
	gr[i]->SetPointError( gr[i]->GetN()-1, 0, Error[i]);
	grAdj[i]    ->SetPoint( grAdj[i]->GetN(),     RunNumber, Output[i]/grLin[CorrectionID]->Eval(Output[i],0,"S"));
	grBase[i]   ->SetPoint( grBase[i]->GetN(),    RunNumber, Output[i]/BaseOut);
	grBaseAdj[i]->SetPoint( grBaseAdj[i]->GetN(), RunNumber, Output[i]/grLin[CorrectionID]->Eval(Output[i],0,"S")/(BaseOut/grLin[CorrectionID]->Eval(BaseOut,0,"S")));
	grAdj[i]    ->SetPointError( grAdj[i]->GetN()-1,0,Error[i]/BaseOut);
	grBase[i]   ->SetPointError( grBase[i]->GetN()-1,0,Error[i]/BaseOut);
	grBaseAdj[i]->SetPointError( grBaseAdj[i]->GetN()-1,0,Error[i]/BaseOut);
	grBaseHeight[i]->SetPoint(grBaseHeight[i]->GetN(),Output[i],Output[i]/BaseOut);
	grBaseHeight[i]->SetPointError( grBaseHeight[i]->GetN()-1,0,Error[i]/BaseOut);
	grBaseTemp[i]->SetPoint(grBaseTemp[i]->GetN(),temp,Output[i]/BaseOut);
      if( i <2240){
	hisRatioAdj->Fill( Output[i],Output[i]/grLin[CorrectionID]->Eval(Output[i],0,"S")/(BaseOut/grLin[CorrectionID]->Eval(BaseOut,0,"S")));
	hisRatio->Fill( Output[i],Output[i]/BaseOut);
	}
      }
    }
  }

  for( int i = 0; i< 2716;i++ ){
    gr[i]->Write();
    grBase[i]->Write();
    grAdj[i]->Write();
    grBaseAdj[i]->Write();
    grBaseHeight[i]->Write();
    grBaseTemp[i]->Write();
  }
  hisRatio->Write();
  hisRatioAdj->Write();
  tfOut->Close();
}
