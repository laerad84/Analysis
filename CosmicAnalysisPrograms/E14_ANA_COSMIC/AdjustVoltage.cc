#include <fstream>
#include <iostream>
#include "TTree.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"

//void testVoltage(){
int
main( int argc, char** argv){

  gStyle->SetPalette(1);
  TApplication* app = new TApplication("App",&argc,argv);
  
  IDHandler* handler = new IDHandler("Data/crystal.txt");
  CsIImage*  image   = new CsIImage(handler);
  TF1* GainCurve_Small_raw = new TF1("GainCurve_Small_raw","TMath::Power(x-657.8,0.9635)",0,5000);
  TF1* GainCurve_Large_raw = new TF1("GainCurve_Large_raw","TMath::Power(x-627.9,1.369)",0,5000);

  TF1* GainCurve_Small = new TF1("GainCurve_Small",
				 Form("%lf*TMath::Power(x-657.8,0.9635)",1/GainCurve_Small_raw->Eval(1750)),
				 0,5000);
  TF1* GainCurve_Large = new TF1("GainCurve_Large",
				 Form("%lf*TMath::Power(x-627.9,1.369)",1/GainCurve_Large_raw->Eval(1750)),
				 0,5000);
  
  TFile* tf = new TFile("test_peak3508.root");
  TGraphErrors* grGain = (TGraphErrors*)tf->Get("gain");
  
  grGain->Draw("AP");
  gPad->Modified();
  gPad->Update();
  
  Double_t DestinationGain      = 1400.0;
  
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
    }else{//in case of Large
      GainOutput = grGain->GetY()[i]/2;
      GainCurve  = GainCurve_Large; 
    }
    
    if( GainOutput > 0){
      std::cout << "TEST" << std::endl;
      Ratio = DestinationGain/GainOutput;
    }else{
      Ratio = 0; 
    }
    
    std::cout << GainOutput      << "\t" 
	      << DestinationGain << " \t"
	      << Ratio           << std::endl;
    if( GainOutput < 500 ){
      
      Voltage = 1750;
      EndAdjGain = GainOutput;

      //std::cout << " ID: " << i << "\tPeak: " << GainOutput << std::endl;

      ofs  << i << "\t" <<GainOutput <<"\t" <<   Voltage << "\t" << EndAdjGain << std::endl;      
      image->Fill(i, GainOutput);
      nBadCH++;
      
    }else{      
      
      Voltage    = GainCurve->GetX(Ratio,700,5000);
      std::cout<< Voltage << std::endl;
      EndAdjGain = DestinationGain;
      //std::cout << Voltage << "\t" << EndAdjGain << std::endl;
      if( Voltage < 1250 ){
	Voltage = 1250;
	EndAdjGain = GainOutput * GainCurve->Eval(1250);  
      }else if( Voltage >1750 ){
	Voltage = 1750;
	EndAdjGain = GainOutput * GainCurve->Eval(1750);
      }
      ofs<<  i << "\t" << GainOutput << "\t" << Voltage << "\t" << EndAdjGain << std::endl;
      std::cout<< GainCurve->Eval(1750) << std::endl;

    }
    
    hisGainRaw->Fill(GainOutput);
    hisGainCal->Fill(EndAdjGain);
    hisVoltage->Fill(Voltage);    

  }
  
  TCanvas* can = new TCanvas("can","",1600,800);
  can->Divide(2,1);
  can->cd(1);
  //image->Draw();
  hisGainRaw->Draw();
  hisGainCal->SetLineColor(2);
  hisGainCal->Draw("same");
  can->cd(2);
  hisVoltage->Draw();

  app->Run();
}



