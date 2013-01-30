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
int main( int argc, char** argv ){
  
  std::string iFilename = "Waveform_Height.root";
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
  /// Set initial waveform /// 
  TSpline3* wav = new TSpline3("wav",gr[20]);  
  TGraph* grVoltTime = new TGraph();
  E14WaveFitter* wavFitter = new E14WaveFitter();
  wavFitter->SetWaveform(wav);
  for( std::vector<int>::iterator it = voltList.begin();
       it != voltList.end();
       it++){
    wavFitter->InitPar();
    wavFitter->SetParameter(0,1);
    wavFitter->SetParameter(1,0);
    wavFitter->SetParameter(2,0);
    wavFitter->SetParLimits(0,0.9,1.1);
    wavFitter->SetParLimits(1,-10,10);
    wavFitter->SetParLimits(2,-0.01,0.01);
    wavFitter->Fit(gr[(*it)/5]);
    gr[(*it)/5]->Draw("AP");
    std::cout<< wavFitter->GetParameter(0) << "\t" 
	     << wavFitter->GetParameter(1) << "\t" 
	     << wavFitter->GetParameter(2) << std::endl;
    can->Update();
    can->Modified();
    grVoltTime->SetPoint( grVoltTime->GetN(),*it, wavFitter->GetParameter(1));
    //getchar();
  }
  

  for( std::vector<int>::iterator it = voltList.begin();
       it != voltList.end();
       it++){
    
  }
  


  //wav->Draw("C");
  grVoltTime->SetMarkerStyle(4);
  grVoltTime->Draw("AP");
  app->Run();
  return 0;
}
