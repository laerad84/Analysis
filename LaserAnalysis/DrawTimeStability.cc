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
#include "TMath.h"
#include "TF1.h"

double laserTimeDelay( double* x , double* par ){
  double value = 0;
  double xpar = x[0];
  value = par[0]*TMath::Log(1-TMath::Exp(-par[1]*par[2])+TMath::Exp(par[1]*(xpar-par[2])))+par[3];
  return value;
}

double gausfunc( double *x, double* p ){
  double value = 0;
  value = p[0]*TMath::Gaus(x[0],p[1],p[2]);
  return value;
}
int main( int argc ,char** argv ){
  TApplication* app = new TApplication("app",&argc,argv);
  const Int_t StrRunNumber = 4158;
  const Int_t EndRunNumber = 4738;

  TGraphErrors* grTimeDelta   = new TGraphErrors();
  TGraphErrors* grTimeDelta1   = new TGraphErrors();
  TGraphErrors* grHeight      = new TGraphErrors();
  TGraphErrors* grTimeHeight  = new TGraphErrors();
  TGraphErrors* grTimeHeight1 = new TGraphErrors();
  TGraphErrors* grHeightLaser = new TGraphErrors();
  TGraphErrors* grHTLaser0    = new TGraphErrors();
  TGraphErrors* grHTLaser1    = new TGraphErrors();

  Int_t ChannelNumber = 1207;
  TGraphErrors* grChDelay = new TGraphErrors();
  grChDelay->SetNameTitle("grChDelay","grChDelay");

  TH1D* hisTimeDelta[1000];
  TH1D* hisHeight[1000];
  TH1D* hisTimeDeltaLaser[1000];
  TF1* funct              = new TF1("funct",laserTimeDelay,0,20000,4);
  TFile* tfDelay          = new TFile("LaserTimeTotal_m100p50.root");
  TGraphErrors * grDel = (TGraphErrors*)tfDelay->Get(Form("grHeightDeltaT_%d",ChannelNumber));
  for( int ip = 0; ip < grDel->GetN(); ip++){
    grChDelay->SetPoint( ip, grDel->GetX()[ip],grDel->GetY()[ip]);
    grChDelay->SetPointError( ip, 0,grDel->GetEY()[ip]);
  }

  funct->SetParameter(0,1);
  funct->SetParameter(1,0.001);
  funct->SetParameter(2,6000);
  funct->SetParameter(3,0);

  funct->SetParLimits(0,1,1);
  funct->SetParLimits(1,0.00001,0.1);
  funct->SetParLimits(2,3000,16000);
  funct->SetParLimits(3,-0.5,0.5);
  grChDelay->Fit(funct,"","",0,16000); 
  //getchar();

  TCanvas* can            = new TCanvas("can","can",1200,800);   

  for( Int_t RunNumber = 0; RunNumber <= EndRunNumber-StrRunNumber; RunNumber++){
    std::cout<< RunNumber + StrRunNumber << std::endl;
    TFile* tf = new TFile(Form("Data/LaserTimeStability_%d.root",RunNumber+StrRunNumber));
    if( tf == NULL ) { continue; }
    if( tf->GetSize() < 1000 ){ continue; }
    hisTimeDelta[RunNumber] = (TH1D*)tf->Get(Form("hisTimeDeltaLaser_%d_%d",ChannelNumber,RunNumber+StrRunNumber));    
    hisHeight[RunNumber]    = (TH1D*)tf->Get(Form("hisHeight_%d_%d",ChannelNumber,RunNumber+StrRunNumber));
    if(hisTimeDelta[RunNumber] == NULL ){ continue; }
    if(hisTimeDelta[RunNumber]->Integral() < 100   ){ continue; }
    if(hisHeight[RunNumber]->GetMean()  > 9000 ){ continue; }

    std::cout<< hisTimeDelta[RunNumber]->GetMean() << std::endl;
    TF1* fTime = new TF1("funcGaussian",gausfunc,-20,20,3);

    fTime->SetParameter(0, hisTimeDelta[RunNumber]->Integral());
    fTime->SetParLimits(0, 0,hisTimeDelta[RunNumber]->Integral()*1.5);
    fTime->SetParameter(1, hisTimeDelta[RunNumber]->GetMean());
    fTime->SetParLimits(1, hisTimeDelta[RunNumber]->GetMean()-3,hisTimeDelta[RunNumber]->GetMean()+3);
    fTime->SetParameter(2, hisTimeDelta[RunNumber]->GetRMS());
    fTime->SetParLimits(2, 0, hisTimeDelta[RunNumber]->GetRMS()*2);

    double minboundary = hisTimeDelta[RunNumber]->GetBinCenter(hisTimeDelta[RunNumber]->GetMaximumBin())-1;
    double maxboundary = hisTimeDelta[RunNumber]->GetBinCenter(hisTimeDelta[RunNumber]->GetMaximumBin())+1;

    hisTimeDelta[RunNumber]->Fit(fTime,"","",minboundary,maxboundary);
    std::cout<< "FitResult" << std::endl;
    std::cout<< fTime->GetParameter(0) << std::endl;
    std::cout<< fTime->GetParameter(1) << std::endl;
    std::cout<< fTime->GetParameter(2) << std::endl;
    if( fTime->GetParError(2) > 0.2){ continue; }

    /*
    hisTimeDelta[RunNumber]->Draw();
    can->Update();
    can->Modified();
    can->Update();
    getchar();
    */

    grTimeDelta->SetPoint(grTimeDelta->GetN(), 
			  RunNumber+StrRunNumber,
			  fTime->GetParameter(1));//-funct->Eval(hisHeight[RunNumber]->GetMean()));
    grTimeDelta->SetPointError(grTimeDelta->GetN()-1,0,fTime->GetParError(2));
    grTimeDelta1->SetPoint(grTimeDelta->GetN(),
			   RunNumber+StrRunNumber,
			   fTime->GetParameter(1)-funct->Eval(hisHeight[RunNumber]->GetMean()));
    grTimeDelta1->SetPointError(grTimeDelta->GetN()-1,0,fTime->GetParError(2));

    grHeight->SetPoint(grHeight->GetN(), RunNumber+StrRunNumber, hisHeight[RunNumber]->GetMean());
    grHeight->SetPointError(grHeight->GetN()-1, 0, hisHeight[RunNumber]->GetRMS());			
    
    grTimeHeight->SetPoint(grTimeHeight->GetN(),hisHeight[RunNumber]->GetMean(),fTime->GetParameter(1));
    grTimeHeight1->SetPoint(grTimeHeight1->GetN(),hisHeight[RunNumber]->GetMean(),fTime->GetParameter(1)-funct->Eval(hisHeight[RunNumber]->GetMean()));
    //grTimeHeight->SetPointError(grTimeHeight->GetN()-1,hisHeight[RunNumber]->GetRMS(),gaus->GetParameter(2));
    tf->Close();
  }

  grTimeDelta->SetNameTitle("grTimeDelta","grTimeDelta;RunNumber;DelayedTime[ns]");
  grHeight->SetNameTitle("grHeight","grHeight;RunNumber;Height of Signal[cnt]");
  grTimeHeight->SetNameTitle("grTimeHeight","grTimeHeight;Height[cnt];Delayed Time[ns]");

  funct->SetParLimits( 3, -20,20);
  grTimeHeight->Fit(funct,"","",0,16000);

  can->Divide(3,2);
  can->Draw();
  can->cd(1);

  funct->Draw();
  can->cd(2);
  grChDelay->Draw("AP");
  can->Update();
  can->Modified();
  can->Update();
  // getchar();

  can->cd(3);  
  grTimeDelta->SetMarkerStyle(20);
  grTimeDelta1->SetMarkerStyle(21);
  grTimeDelta->SetMarkerSize(0.5);
  grTimeDelta1->SetMarkerSize(0.5);
  grTimeDelta1->SetMarkerColor(2);
  grTimeDelta->Draw("AP");  
  grTimeDelta1->Draw("P");
  can->cd(4);
  grTimeHeight->SetMarkerStyle(21);
  grTimeHeight1->SetMarkerStyle(22);
  grTimeHeight1->SetMarkerColor(2);
  grTimeHeight->Draw("AP");

  grTimeHeight1->Draw("P");
  can->cd(5);
  grHeight->Draw("AP");
  TCanvas* can1  = new TCanvas("can1","can1",800,400);
  can1->Draw();
  TPad* pad1 = new TPad("pad1","pad1",0,0.45,1,1);
  TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.45);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  grTimeDelta->Draw("AP");  
  grTimeDelta1->Draw("P");
  pad2->cd();  
  grHeight->Draw("AP");
  app->Run();
}
