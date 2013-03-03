//void DrawTimeStablity(){
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TApplication.h"

int main( int argc ,char** argv ){
  TApplication* app = new TApplication("app",&argc,argv);
  const Int_t StrRunNumber = 4158;
  const Int_t EndRunNumber = 4738;

  TGraphErrors* grTimeDelta   = new TGraphErrors();
  TGraphErrors* grHeight      = new TGraphErrors();
  TGraphErrors* grTimeHeight  = new TGraphErrors();
  TGraphErrors* grHeightLaser = new TGraphErrors();
  TH1D* hisTimeDelta[1000];
  TH1D* hisHeight[1000];
  
  TCanvas* can = new TCanvas("can","can",800,800);
  Int_t ChannelNumber = 1206;
  for( Int_t RunNumber = 0; RunNumber <= EndRunNumber-StrRunNumber; RunNumber++){
    std::cout<< RunNumber + StrRunNumber << std::endl;
    TFile* tf = new TFile(Form("Data/LaserTimeStability_%d.root",RunNumber+StrRunNumber));
    if( tf == NULL ) { continue; }
    if( tf->GetSize() < 1000 ){ continue; }
    hisTimeDelta[RunNumber] = (TH1D*)tf->Get(Form("hisTimeDeltaLaser_%d_%d",ChannelNumber,RunNumber+StrRunNumber));    
    hisHeight[RunNumber]    = (TH1D*)tf->Get(Form("hisHeight_%d_%d",ChannelNumber,RunNumber+StrRunNumber));
    if(hisTimeDelta[RunNumber] == NULL ){ continue; }
    if(hisTimeDelta[RunNumber]->Integral() < 100   ){ continue; }
    if(hisTimeDelta[RunNumber]->GetMean()  > 14000 ){ continue; }
    
    std::cout<< hisTimeDelta[RunNumber]->GetMean() << std::endl;
    TF1* fTime = new TF1("funcGaus","[0]*TMath::Gaus(x,[1],[2])",-20,20);
    fTime->SetParameter(0, hisTimeDelta[RunNumber]->Integral());
    fTime->SetParLimits(0, hisTimeDelta[RunNumber]->Integral()*0.5,hisTimeDelta[RunNumber]->Integral()*1.5);
    fTime->SetParameter(1, hisTimeDelta[RunNumber]->GetMean());
    fTime->SetParLimits(1, hisTimeDelta[RunNumber]->GetMean()-3,hisTimeDelta[RunNumber]->GetMean()+3);
    hisTimeDelta[RunNumber]->Fit(fTime,"Q","",hisTimeDelta[RunNumber]->GetBinCenter(hisTimeDelta[RunNumber]->GetMaximumBin())-1,hisTimeDelta[RunNumber]->GetBinCenter(hisTimeDelta[RunNumber]->GetMaximumBin())+1);
    grTimeDelta->SetPoint(grTimeDelta->GetN()  , RunNumber+StrRunNumber, fTime->GetParameter(1));
    grTimeDelta->SetPointError(grTimeDelta->GetN()-1,0,fTime->GetParameter(2));
    grHeight->SetPoint(grHeight->GetN(), RunNumber+StrRunNumber, hisHeight[RunNumber]->GetMean());
    grHeight->SetPointError(grHeight->GetN()-1, 0, hisHeight[RunNumber]->GetRMS());			
    grTimeHeight->SetPoint(grTimeHeight->GetN(),hisHeight[RunNumber]->GetMean(),fTime->GetParameter(1));
    //grTimeHeight->SetPointError(grTimeHeight->GetN()-1,hisHeight[RunNumber]->GetRMS(),gaus->GetParameter(2));
    tf->Close();
  }
  can->Divide(2,2);
  can->cd(1);  
  grTimeDelta->Draw("AP");  
  can->cd(2);
  grTimeHeight->SetMarkerStyle(21);
  grTimeHeight->Draw("AP");
  can->cd(3);
  grHeight->Draw("AP");
  app->Run();
}
