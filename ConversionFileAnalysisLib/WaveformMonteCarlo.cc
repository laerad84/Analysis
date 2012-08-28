#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TTree.h"

#include "TGraph.h"
#include "TSpline.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

#include "E14WaveFitter.h"

int
main( int argc, char** argv ){
  const Int_t nChannel = 2716;
  const Double_t FluctPedestal = 1;
  TApplication* app = new TApplication("App", &argc , argv );

  E14WaveFitter* Fitter = new E14WaveFitter();

  TCanvas* can = new TCanvas("can", "can" ,0,0,800,800);
  can->Divide(2,2);

  TFile* tfin = new TFile("TEMPLATE_SPLINE_250_500.root");
  
  std::cout << "Read Data File " << std::endl;

  TGraph* gr[nChannel];
  TSpline3* spl[nChannel];
  for( Int_t ch = 0; ch < 2716; ch++){
    std::cout<< ch << "\t";
    gr[ch] = (TGraph*)tfin->Get(Form("Template_graph_%d",ch));
    if( gr[ch]->GetN() == 0 ){ continue; }
    spl[ch] = new TSpline3(Form("Spline_%d",ch),gr[ch]);
  }
  std::cout<< std::endl;

  std::cout << "End Read File " << std::endl;

  Double_t TimeOffset = 150;
  Double_t TimeDelta;
  Double_t Height;
  Double_t Noise; 
  Double_t Signal;

  TGraph* grWaveform = new TGraph();
  grWaveform->SetMarkerStyle(6);
  
  std::cout << "Draw " << std::endl;
  can->cd(1);
  gr[0]->Draw("AP");
  can->cd(2);
  spl[0]->Draw("AP");
  can->Update();
  can->Modified();

  
  getchar();
  
  Fitter->SetWaveform( spl[0] );
  TH1D* hisTime = new TH1D( "hisTimeDelta" ,"", 800,-8,8);
  TH1D* hisHeight = new TH1D("hisHeight","Height",800,-0.2,0.2);
  TH1D* hisChisquare = new TH1D("hisChisquare","Chisquare",1000,0,1000);
  Int_t nLoop = 0; 
  


  while(1){

    TimeDelta = 8*gRandom->Rndm();
    Height    = 1000; 
    grWaveform->Set(0);
    for( int i = 0; i< 48; i++){

      Noise = gRandom->Gaus(0,FluctPedestal);
      Signal = (int)(Height*spl[0]->Eval( i*8 -TimeOffset + TimeDelta ) + Noise); 
      std::cout<< i << " : " << Signal << std::endl;
      grWaveform->SetPoint(i,i*8, Signal);      
    }

    Fitter->InitPar();
    Fitter->Fit(grWaveform);
    grWaveform->GetListOfFunctions()->Delete();
    std::cout << Fitter->GetChisquare()/ Fitter->GetNDF() << std::endl;
    hisTime->Fill(Fitter->GetParameter(1) -(TimeOffset - TimeDelta));
    hisHeight->Fill((Fitter->GetParameter(0)- Height)/Height);
    hisChisquare->Fill(Fitter->GetChisquare()/Fitter->GetNDF());
    if( nLoop % 1000 == 0 ){
      can->cd(1);
      gPad->SetGridx();
      gPad->SetGridy();
      hisChisquare->Draw();
      can->cd(2);
      gPad->SetGridx();
      gPad->SetGridy();
      hisHeight->Draw();
      can->cd(3);
      gPad->SetGridx();
      gPad->SetGridy();
      grWaveform->Draw("AP");
      Fitter->m_FitFunc->SetLineWidth(1);
      Fitter->m_FitFunc->Draw("same");
      can->cd(4);
      hisTime->Draw();
      can->Update();
      can->Modified();
      getchar();
    }

    Fitter->Clear();
    nLoop++;
  }
  
  app->Run();
  return 0;
}
