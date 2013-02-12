#include <iostream>
#include "TFile.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "E14WaveFitter.h"
#include "E14WaveformAnalyzer.h"
#include <string>
#include <vector>
#include "TApplication.h"
#include "TCanvas.h"
#include "TTree.h"
#include "Waveform_Reader.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStyle.h"

double func( double*x ,double*par ){  
  double value = par[0]*(-1*TMath::Log(1 + par[1])+TMath::Log( 1 + par[1]*exp( par[2]*x[0] )));
  return value;
}
int main( int argc, char** argv ){
  gStyle->SetOptFit(111111111);
  std::string iFilename = "Waveform_Height.root";
  std::string iDataname = "run_conv_mode5.root";

  TFile* tf = new TFile(iFilename.c_str());

  TH2D* hisTemplate[120];
  TProfile* prof[120];
  TGraph*   gr[120];
  for( int i = 0; i< 120; i++){
    hisTemplate[i] = (TH2D*)tf->Get(Form("hisTemplate_%d",i));
    prof[i] = hisTemplate[i]->ProfileX();
    gr[i] = new TGraph();
    gr[i]->SetNameTitle(Form("Template_%d",i),
			Form("Template_%d",i*5));   
  }
  TFile* tfData = new TFile(iDataname.c_str());
  TTree* trData = (TTree*)tfData->Get("EventTree");
  Waveform_Reader* r = new Waveform_Reader(trData);

  TApplication* app = new TApplication("app",&argc,argv);
  TCanvas*can = new TCanvas("can","",800,800);
  std::vector<int> voltList;
  for( int i = 0; i< 120; i++){
    if( hisTemplate[i]->GetEntries() != 0){
      voltList.push_back( i*5 );
      Double_t x,y;
      for( int ipoint = 0; ipoint < prof[i]->GetNbinsX();ipoint++){
	x = prof[i]->GetBinCenter(ipoint+1);
	y = prof[i]->GetBinContent(ipoint+1);	
	gr[i]->SetPoint(ipoint,x,y);
      }
    }
  }

  std::cout << __LINE__ << std::endl;
  /// Set initial waveform /// 
  TSpline3* wav            = new TSpline3("wav",gr[10]);  
  TGraph*   grVoltTime     = new TGraph();
  E14WaveFitter* wavFitter = new E14WaveFitter();
  TH1D* hisHeight[120];
  TH1D* hisTime[120];
  for( int i = 0; i< 120; i++){
    hisHeight[i] = new TH1D(Form("hisHeight_%d",i),Form("hisHeight_%d",i),160,0,16000);
    hisTime[i]   = new TH1D(Form("hisTime_%d",i)  ,Form("hisTime_%d",i)  ,500,0,500);
  }
  TGraph* grWave = new TGraph();
  std::cout << __LINE__ << std::endl;
  for( int ievent  = 0; ievent < trData->GetEntries(); ievent++){
    r->GetEntry(ievent);
    Int_t voltIndex = ((r->volt)*(-1))/5;
    if( voltIndex < 0 || voltIndex > 119 ){ std::cout<< "Out of region" << std::endl; }
    grWave->Set(0);
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      grWave->SetPoint( ipoint, r->xaxis[ipoint]*8, r->data[ipoint]);
    }
    /*
    wavFitter->InitPar();
    wavFitter->SetWaveform(wav);
    wavFitter->SetParameter(0,r->param[0]);
    wavFitter->SetParLimits(0,r->param[0]*0.8,r->param[0]*1.2);
    wavFitter->SetParameter(1,r->param[1]);
    wavFitter->SetParLimits(1,r->param[1]-8,r->param[1]+8);
    wavFitter->SetParameter(2,r->param[4]);
    wavFitter->SetParLimits(2,r->param[4]-8,r->param[4]+8);

    wavFitter->Fit(grWave);
    std::cout<< wavFitter->GetParameter(0) << "\t" 
	     << wavFitter->GetParameter(1) << "\t"
	     << wavFitter->GetParameter(2) << std::endl;
    std::cout<< r->param[0] << "\t"
	     << r->param[1] << "\t"
	     << r->param[2] << "\t"
	     << r->param[3] << "\t"
	     << r->param[4] << std::endl;
    grWave->Draw("AP");
    can->Update();
    can->Modified();
    getchar();
    hisHeight[voltIndex]->Fill( wavFitter->GetParameter(0));
    hisTime[voltIndex]->Fill( wavFitter->GetParameter(1));				
    */    
    hisHeight[voltIndex]->Fill(r->param[0]);
    hisTime[voltIndex]->Fill(r->param[1]);				
    grWave->GetListOfFunctions()->Delete();
    wavFitter->Clear();
  }

  TGraphErrors* grTimeDelta =new TGraphErrors();
  for( int i = 0; i< 120; i++){
    if( hisHeight[i] ->GetEntries() == 0 ){ continue; }
    grTimeDelta->SetPoint(grTimeDelta->GetN(),hisHeight[i]->GetMean(),hisTime[i]->GetMean()-148.25);
    //grTimeDelta->SetPointError(grTimeDelta->GetN()-1,hisHeight[i]->GetRMS(),hisTime[i]->GetRMS());
  }
  TFile* tfLaser = new TFile("LaserTimeTotal_m100p50.root");
  TH2D*  hisHeightTimeLaser = (TH2D*)tfLaser->Get("hisDelta2D_325");
  TH1D*  hisDeltaStd        = (TH1D*)tfLaser->Get("hisDeltaStd_325");
  TProfile* profTimeLaser = hisHeightTimeLaser->ProfileX();
  TGraph* grTimeDeltaLaser = new TGraph();
  for( int ibin = 0; ibin < profTimeLaser->GetNbinsX(); ibin++){
    if( profTimeLaser->GetBinEntries(ibin+1) < 10 ){ continue; }
    grTimeDeltaLaser->SetPoint( grTimeDeltaLaser->GetN(), profTimeLaser->GetBinCenter(ibin+1),profTimeLaser->GetBinContent(ibin+1)-hisDeltaStd->GetMean());
  }

  //wav->Draw("C");
  //grVoltTime->SetMarkerStyle(4);
  //grVoltTime->Draw("AP");
  grTimeDelta->SetMarkerStyle(4);
  grTimeDeltaLaser->SetMarkerStyle(5);
  grTimeDeltaLaser->Draw("AP");
  grTimeDelta->SetMarkerColor(2);
  grTimeDelta->Draw("P");
  TF1* Fitfunc = new TF1("func",func,0,16000,3);
  Fitfunc->SetParameter(0, 1);
  Fitfunc->SetParameter(1, 0.001);
  Fitfunc->SetParameter(2, 0.001);
  grTimeDelta->Fit(Fitfunc,"","",0,16000);

  app->Run();
  return 0;
}
