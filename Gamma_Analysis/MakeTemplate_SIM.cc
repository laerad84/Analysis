#include "PulseGenerator.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1.h"
#include "TApplication.h"
#include "TCanvas.h"
#include <vector>
#include "WaveformFitter.h"
#include "TFile.h"
#include "TH2D.h"
#include "TSpline.h"
#include "TProfile.h"
#include "TSpline.h"
#include "E14WaveFitter.h"
#include "E14WaveformAnalyzer.h"
int main( int argc, char** argv ){
  
  TApplication* app = new TApplication("app",&argc,argv);
  
  WaveformFitter* wavFitter = new WaveformFitter(48,kFALSE);
  PulseGenerator* gen = new PulseGenerator();
  TCanvas* can = new TCanvas("can","",800,800);
  std::vector<double> EVector;
  std::vector<double> TVector;
  Double_t meanTime = 0; 
  Double_t Energy   = 0; 
  TFile* tfout = new TFile("WaveTemplate_Sim.root","RECREATE");
  
  TH2D* hisTemplate_Waveform = new TH2D(Form("hisTemplate_Waveform"),
					 Form("hisTemplate_Waveform"),
					450,-150,300,150,-0.25,1.25);  
  TH1D* hisTimeDelta = new TH1D("hisTimeDelta","hisTimeDelta",
				500,150,200);
  E14WaveFitter* wavSplFitter = new E14WaveFitter();
  E14WaveformAnalyzer* wavAnalyzer = new E14WaveformAnalyzer();
  for( int i = 0; i< 10000 ; i++){
    //std::cout<< "Generator" << std::endl;
    can->cd(1);
    //gen->Reset();
    Energy = (gRandom->Rndm()*200 +200)/17.52;
    meanTime = 8*(gRandom->Rndm()-0.5);
    TF1* func = gen->GetWaveform( Energy,meanTime);
    const Int_t nPoint= 48;
    TGraph* gr = new TGraph(nPoint);
    for( int ipoint = 0 ; ipoint < nPoint; ipoint++){
      double offset = 8*(gRandom->Rndm()-0.5);
      double tsum   = 0.;
      double d_tmp  = func->Eval(ipoint*8.-offset)+gRandom->Gaus(0.,2.4);
      gr->SetPoint(ipoint,ipoint*8,d_tmp);
    }
    
    bool fit = wavFitter->Fit( gr);
    if( !fit ){ continue; }

    TF1* fitFunc = wavFitter->GetFunction();
    if( fitFunc->GetChisquare()/fitFunc->GetNDF() > 100 ){ continue; }
    Double_t m_Signal = fitFunc->GetParameter(0);
    Double_t m_Time   = fitFunc->GetParameter(1);
    Double_t m_Pedestal = fitFunc->GetParameter(4);
    for( int i = 0; i< 48; i++){
      hisTemplate_Waveform->Fill(gr->GetX()[i] - m_Time, (gr->GetY()[i]-m_Pedestal)/m_Signal);
    }
    hisTimeDelta->Fill(m_Time-meanTime);
    delete gr;
  }

  TGraph* gre = new TGraph();
  TProfile* profTemplate_Waveform = hisTemplate_Waveform->ProfileX();
  for( int xIndex = 0; xIndex< hisTimeDelta->GetNbinsX();xIndex++){
    Double_t xCenter = profTemplate_Waveform->GetBinCenter(xIndex+1);
    Double_t yValue  = profTemplate_Waveform->GetBinContent(xIndex+1);
    gre->SetPoint( xIndex, xCenter,yValue );
  }
  TSpline3* spl = new TSpline3("Spline",gre);

  // Renormalize gre // 
  TGraph* grn = new TGraph();
  gre->Fit("pol0","","",-120,-100);
  TF1* func = gre->GetFunction("pol0");
  Double_t gnd  = func->GetParameter(0);
  Double_t peak = spl->Eval(0);
  for( int i = 0; i< gre->GetN(); i++){
    grn->SetPoint( i, gre->GetX()[i], (gre->GetY()[i]-gnd)/peak);
  }
  TSpline3* splTemplate = new TSpline3("splTemplate",grn);
  wavSplFitter->SetWaveform(splTemplate);
  TH2D* hisTemplate_Waveform_Delta_min = new TH2D("hisTemplate_Waveform_Delta_min",
						  "hisTemplate_Waveform_Delta_min",
						  450,-150,300,150,-0.25,1.25);
  TH1D* hisTimeDelta_Waveform_min = new TH1D("hisTimeDelta_Waveform_min",
					     "hisTimeDelta_Waveform_min",
					     500,150,200);
  for( int i = 0; i< 10000 ; i++){
    wavSplFitter->InitPar();
    wavAnalyzer->_Clear();
    //std::cout<< "Generator" << std::endl;

    can->cd(1);
    //gen->Reset();
    Energy = (gRandom->Rndm()*200 +200)/17.52*0.1;
    meanTime = 8*(gRandom->Rndm()-0.5);
    TF1* func = gen->GetWaveform( Energy,meanTime);
    const Int_t nPoint= 48;
    TGraph* gr = new TGraph(nPoint);
    for( int ipoint = 0 ; ipoint < nPoint; ipoint++){
      double offset = 8*(gRandom->Rndm()-0.5);
      double tsum   = 0.;
      double d_tmp  = func->Eval(ipoint*8.-offset)+gRandom->Gaus(0.,2.4);
      gr->SetPoint(ipoint,ipoint*8,d_tmp);
    }
    wavAnalyzer->AnalyzeWaveform( gr->GetY() );
    wavSplFitter->SetParameter(wavAnalyzer);
    bool fit = wavSplFitter->Fit( gr);
    if( !fit ){ continue; }
    TF1* fitFunc = wavSplFitter->GetFunction();
    Double_t m_Signal = fitFunc->GetParameter(0);
    Double_t m_Time   = fitFunc->GetParameter(1);
    Double_t m_Pedestal = fitFunc->GetParameter(4);
    for( int i = 0; i< 48; i++){
      hisTemplate_Waveform_Delta_min->Fill(gr->GetX()[i] - m_Time, (gr->GetY()[i]-m_Pedestal-fitFunc->Eval(gr->GetX()[i]))/m_Signal);
    }    
    hisTimeDelta_Waveform_min->Fill(m_Time-meanTime);
    delete gr;
  }
  hisTimeDelta->Write();
  hisTimeDelta_Waveform_min->Write();
  hisTemplate_Waveform->Write();
  hisTemplate_Waveform_Delta_min->Write();
  tfout->Close();
  //app->Run();

}
