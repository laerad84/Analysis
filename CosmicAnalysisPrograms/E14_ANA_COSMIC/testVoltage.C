#include <fstream>
#include <iostream>
#include "TTree.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"

void testVoltage(){
  gSystem->AddIncludePath("include");
  gSystem->Load("lib/libtest.so");
  gStyle->SetPalette(1);
  
  IDHandler* handler = new IDHandler("Data/crystal.txt");
  CsIImage*  image   = new CsIImage(handler);
  
  
  
  TF1* GainCurve_Small = new TF1("GainCurve_Small","TMath::Power(x-657.8,0.9635)",0,5000);
  TF1* GainCurve_Large = new TF1("GainCurve_Large","TMath::Power(x-627.9,1.369)",0,5000);
  
  TFile* tf = new TFile("test_peak3508.root");
  TGraphErrors* grGain = (TGraphErrors*)tf->Get("gain");
  
  grGain->Draw("AP");
  gPad->Modified();
  gPad->Update();
  const Double_t corFactor1500V_Small = GainCurve_Small->Eval(1500)/GainCurve_Small->Eval(1750);
  const Double_t corFactor1500V_Large = GainCurve_Large->Eval(1500)/GainCurve_Large->Eval(1750);    
  Double_t DestinationGain      = 1400.0;
  std::cout << corFactor1500V_Small << std::endl;
  std::cout << corFactor1500V_Large << std::endl;
  
  std::ofstream ofs("GainAdjustVoltage.txt");
  
  TH1D* hisGainRaw = new TH1D("hisGainRaw","hisGainRaw",400,0,16000);
  TH1D* hisGainCal = new TH1D("hisGainCal","hisGainCal",400,0,16000);
  TH1D* hisVoltage = new TH1D("hisVoltage","hisVoltage",180,0,1800);
  
  int nBadCH = 0;
  for( int i = 0; i< 2716; i++){
    if( i >=2240     && i< 2240+120 ){continue;}
    if( i >=2716-120 && i< 2716     ){continue;}
    
    Double_t GainOutput;
    Double_t CorrectedGain;
    Double_t Voltage;
    Double_t EndAdjGain;
    Double_t CorrectionFactor;
    Double_t Ratio;
    
    
    TF1* GainCurve;
    
    if( i < 2240 ){//in case of Small
      GainOutput = grGain->GetY()[i];
      GainCurve = GainCurve_Small;
      CorrectionFactor = corFactor1500V_Small;
    }else{//in case of Large
      GainOutput = grGain->GetY()[i]/2;
      GainCurve  = GainCurve_Large; 
      CorrectionFactor = corFactor1500V_Large;
    }
    
    CorrectedGain = GainOutput*CorrectionFactor;
    
    if( CorrectedGain > 0){
      std::cout << "TEST" << std::endl;
      Raitio = DestinationGain/CorrectedGain;
    }else{
      Ratio = 0; 
    }
    
    std::cout <<CorrectedGain << " " << DestinationGain << " \t" << DestinationGain/CorrectedGain<< "\t" << Ratio << std::endl;
    
    if( GainOutput < 500 ){
      Voltage = -1;
      //std::cout << " ID: " << i << "\tPeak: " << GainOutput << std::endl;
      std::ofs  << i << "\t" << Voltage << "\t" << CorrectedGain << std::endl;      
      image->Fill(i, GainOutput);
      nBadCH++;
    }else{
      Voltage    = GainCurve->GetX(DestinationGain/CorrectedGain,600,5000);
      std::cout<< Voltage << std::endl;
      EndAdjGain = DestinationGain;
      //std::cout << Voltage << "\t" << EndAdjGain << std::endl;
      if( Voltage < 1000 ){
	Voltage = 1000;
	EndAdjGain = GainOutput * CorrectionFactor * GainCurve->Eval(1000);  
      }else if( Voltage >1750 ){
	Voltage = 1750;
	EndAdjGain = GainOutput * CorrectionFactor * GainCurve->Eval(1750);
      }
      std::ofs  << i << "\t" << Voltage << "\t" << EndAdjGain << std::endl;     
    }
    
    hisGainRaw->Fill(GainOutput);
    hisGainCal->Fill(EndAdjGain);
    hisVoltage->Fill(Voltage);
    
  }
  
  //image->Draw();
  hisGainRaw->Draw();
  hisGainCal->Draw("same");
}



