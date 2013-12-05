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
#include "IDHandler.h"

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
  IDHandler* handler =new IDHandler();

  const int nAll = 495;
  TFile* tf[nAll];
  std::vector<int> RunNumberVec;


  TChain* ch = new TChain("trLaser");
  
  //std::ifstream ifs("../RunList/RunAll.csv");
  std::ifstream ifs("../RunList/RunAllAdj.txt");
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

  Double_t Initial[2716];

  for( int i = 0; i< 2716; i++){
    Initial[i] = 0;
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
  TGraphErrors* grBaseHeightAdj[2716];
  TGraphErrors* grBaseDrift[2716];
  TGraphErrors* grBaseDriftAdj[2716];
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
    grBaseHeightAdj[i] = new TGraphErrors();
    grBaseHeightAdj[i]->SetNameTitle(Form("LaserOutputBaseHeightAdj_%d",i),Form("LaserOutputBaseHeightAdj_%d",i));
    grBaseDrift[i] = new TGraphErrors();
    grBaseDrift[i]->SetNameTitle(Form("LaserOutputDrift_%d",i),Form("LaserOutputDrift_%d",i));
    grBaseDriftAdj[i] = new TGraphErrors();
    grBaseDriftAdj[i]->SetNameTitle(Form("LaserOutputDriftAdj_%d",i),Form("LaserOutputDriftAdj_%d",i));
    grBaseTemp[i] = new TGraphErrors();
    grBaseTemp[i]->SetNameTitle(Form("LaserOutputBaseTemp_%d",i),Form("LaserOutputBaseTemp_%d",i));

  }




  int baseID[4] = {22, 35,2212, 2227};
  /*
     2|3
    -----
     0|1
  */



  TH2D* hisRatio = new TH2D("hisRatio","hisRatio",160,0,16000,200,0,10);
  TH2D* hisRatioAdj = new TH2D("hisRatioAdj","hisRatioAdj",160,0,16000,200,0,10);
  for( int ientries = 0 ; ientries < ch->GetEntries(); ientries++){
    ch->GetEntry( ientries );
    if( Entries[0] < 100 ){ continue; }
    trTemp->GetEntry( RunNumber );
    if( ientries == 0){
      for( int i = 0; i< 2716; i++){
	Initial[i] = Output[i];
      }
    }

    Double_t BaseOut[4];
    for( int i = 0; i< 4; i++){
      BaseOut[i] = Output[baseID[i]];
    }

    Double_t BaseInit = Initial[22];
    for( int i = 0; i< 2716; i++){
      int RegionID;
      double x,y;
      handler->GetMetricPosition(i,x,y);
      if( y < 0 ){
	if( x < 0 ){
	  RegionID = 0;
	}else{
	  RegionID = 1;
	}
      }else{
	if( x < 0 ){
	  RegionID = 2;
	}else{
	  RegionID = 3;
	}
      }

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
	double cconst[2];
	cconst[0] = grLin[CorrectionID]->Eval(Output[i],0,"S");
	cconst[1] = grLin[CorrectionID]->Eval(BaseOut[RegionID],0,"S");
	
	if( Initial[i] != 0 ){
	  /*
	  grBaseDrift[i]->SetPoint(grBaseDrift[i]->GetN(),Output[i],Output[i]/BaseOut[RegionID]/Initial[i]*BaseInit);
	  grBaseDrift[i]->SetPointError( grBaseDrift[i]->GetN()-1,0,Error[i]/BaseOut[RegionID]/Initial[i]*BaseInit);
	  grBaseDriftAdj[i]->SetPoint(grBaseDriftAdj[i]->GetN(),Output[i],(Output[i]/cconst[0])/(BaseOut[RegionID]/cconst[1])/Initial[i]*BaseInit);
	  grBaseDriftAdj[i]->SetPointError( grBaseDriftAdj[i]->GetN()-1,0,Error[i]/(BaseOut[RegionID]/cconst[1])/Initial[i]*BaseInit);
	  */
	  if( Error[i]/Output[i]  > 0.01 ){ continue; }
	  gr[i]->SetPoint( gr[i]->GetN(), RunNumber,Output[i]);
	  gr[i]->SetPointError( gr[i]->GetN()-1, 0, Error[i]);
	  grAdj[i]    ->SetPoint( grAdj[i]->GetN(),     RunNumber, Output[i]/grLin[CorrectionID]->Eval(Output[i],0,"S"));
	  grBase[i]   ->SetPoint( grBase[i]->GetN(),    RunNumber, Output[i]/BaseOut[RegionID]);
	  grBaseAdj[i]->SetPoint( grBaseAdj[i]->GetN(), RunNumber, Output[i]/grLin[CorrectionID]->Eval(Output[i],0,"S")/(BaseOut[RegionID]/grLin[CorrectionID]->Eval(BaseOut[RegionID],0,"S")));
	  grAdj[i]    ->SetPointError( grAdj[i]->GetN()-1,0,Error[i]/BaseOut[RegionID]);
	  grBase[i]   ->SetPointError( grBase[i]->GetN()-1,0,Error[i]/BaseOut[RegionID]);
	  grBaseAdj[i]->SetPointError( grBaseAdj[i]->GetN()-1,0,Error[i]/BaseOut[RegionID]);
	  grBaseHeight[i]->SetPoint(grBaseHeight[i]->GetN(),Output[i],Output[i]/BaseOut[RegionID]);
	  grBaseHeight[i]->SetPointError( grBaseHeight[i]->GetN()-1,0,Error[i]/BaseOut[RegionID]);
	  grBaseHeightAdj[i]->SetPoint(grBaseHeightAdj[i]->GetN(),Output[i],(Output[i]/cconst[0])/(BaseOut[RegionID]/cconst[1]));
	  grBaseHeightAdj[i]->SetPointError( grBaseHeightAdj[i]->GetN()-1,0,Error[i]/(BaseOut[RegionID]/cconst[1]));
	  grBaseTemp[i]->SetPoint(grBaseTemp[i]->GetN(),temp,Output[i]/BaseOut[RegionID]);
	  
	  grBaseDrift[i]->SetPoint(grBaseDrift[i]->GetN(),RunNumber,Output[i]/BaseOut[RegionID]/Initial[i]*BaseInit);
	  grBaseDrift[i]->SetPointError( grBaseDrift[i]->GetN()-1,0,Error[i]/BaseOut[RegionID]/Initial[i]*BaseInit);
	  grBaseDriftAdj[i]->SetPoint(grBaseDriftAdj[i]->GetN(),RunNumber,(Output[i]/cconst[0])/(BaseOut[RegionID]/cconst[1])/Initial[i]*BaseInit);
	  grBaseDriftAdj[i]->SetPointError( grBaseDriftAdj[i]->GetN()-1,0,Error[i]/(BaseOut[RegionID]/cconst[1])/Initial[i]*BaseInit);
	}
      if( i <2240){
	hisRatioAdj->Fill( Output[i],Output[i]/grLin[CorrectionID]->Eval(Output[i],0,"S")/(BaseOut[RegionID]/grLin[CorrectionID]->Eval(BaseOut[RegionID],0,"S")));
	hisRatio->Fill( Output[i],Output[i]/BaseOut[RegionID]);
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
    grBaseHeightAdj[i]->Write();
    grBaseDrift[i]->Write();
    grBaseDriftAdj[i]->Write();
    grBaseTemp[i]->Write();
  }
  hisRatio->Write();
  hisRatioAdj->Write();
  tfOut->Close();
}
