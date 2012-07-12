#include <cstdlib>
#include <cstdio>

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
  const Double_t HVLowLimit  = 1200;//V
  const Double_t HVHighLimit = 1750;//V    
  Double_t DestinationGain   = 1500.0;
  if( argc ==2){
    DestinationGain = atof(argv[1]);
  }

  std::cout << " //////////////////////////////////////////////////////////////////////////////\n"
	    << " Gain Destination : " << DestinationGain
	    << " //////////////////////////////////////////////////////////////////////////////\n"
	    << std::endl;
  
  TApplication* app = new TApplication("App",&argc,argv);  

  gStyle->SetPalette(1);
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("NeMRI");
  
  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////

  std::cout << "Initialize IMAGE" << std::endl;
  IDHandler* handler    = new IDHandler();
  CsIImage*  imageGain[3];
  CsIImage*  imageVolt[3];  
  for( int i = 0; i< 3; i++){
    imageGain[i] = new CsIImage(handler);
    imageVolt[i] = new CsIImage(handler);
  }  
  imageGain[0]->SetTitle("Before Gain Tune");
  imageGain[1]->SetTitle("After  Gain Tune");
  imageGain[2]->SetTitle("Ratio After/Before Tune");
  imageVolt[0]->SetTitle("Before Voltage");
  imageVolt[1]->SetTitle("After Voltage");
  imageVolt[2]->SetTitle("Ratio Volatage");  
  std::cout << "Initialized Image" << std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////
  
  std::cout << "Gain Curve" << std::endl;
  
  TF1* GainCurve_Small = new TF1("GainCurve_Small_raw","TMath::Power(x-657.8,0.9635)",0,5000);
  TF1* GainCurve_Large = new TF1("GainCurve_Large_raw","TMath::Power(x-627.9,1.369)",0,5000);

    
  std::cout <<"ReadFile" << std::endl;
  std::ifstream ifs("Data/CH_HV_LIST.txt");
  Int_t    CH;
  Double_t HVSet;
  Double_t CH_HV_MAP[2716]      = {0};
  Double_t VoltageInitial[2716] = {0};
  Double_t VoltageFinetune[2716]= {0};
  Double_t GainBefore[2716]     = {0};
  Double_t GainAfter[2716]      = {0};
  
  while( !ifs.eof() ){
    ifs >> CH >> HVSet;
    CH_HV_MAP[CH] = HVSet;
    VoltageInitial[CH] = HVSet;
  }
  
  //TFile* tf = new TFile("cosmic_3599_3600.root");
  TFile* tf = new TFile("ANA_Cosmic_FEB_2.root");
  
  std::ofstream ofs("GainAdjustVoltage_3599_3600.txt");
  
  TH1D* hisGainRaw = new TH1D("hisGainRaw","hisGainRaw",120,0,6000);
  TH1D* hisGainCal = new TH1D("hisGainCal","hisGainCal",120,0,6000);
  TH1D* hisVoltageInit = new TH1D("hisVoltageInit","hisVoltageInit",180,0,1800);
  TH1D* hisVoltageFine = new TH1D("hisVoltageFine","hisVoltageFine",180,0,1800);
  
  
  int nBadCH = 0;
  for( int i = 0; i< 2716; i++){
    if( i >=2240     && i< 2240+120 ){continue;}
    if( i >=2716-120 && i< 2716     ){continue;}

    std::cout<< "CHANNEL:" <<  i << std::endl;
    
    TH1D* his = (TH1D*)tf->Get(Form("his_CH%04d",i));        
    Double_t GainOutput;
    Double_t CorrectedGain    = 0;
    Double_t Voltage          = 0;
    Double_t EndAdjGain       = 0;
    Double_t CorrectionFactor = 0;
    Double_t Ratio            = 0;
    
    TF1* GainCurve;
    TF1* gainlandau=NULL;

    if( i < 2240 ){//in case of Small
      //GainOutput = grGain->GetY()[i];
      //GainCurve = GainCurve_Small;
      gainlandau = his->GetFunction("landau");
      if(gainlandau ==NULL){continue;}
      GainOutput = gainlandau->GetParameter(1);
      GainCurve = GainCurve_Small;
    }else{//in case of Large
      //GainOutput = grGain->GetY()[i]/2;
      //GainCurve = GainCurve_Small;
      gainlandau = his->GetFunction("landau");
      if(gainlandau ==NULL){continue;}
      GainOutput = gainlandau->GetParameter(1)/2;
      GainCurve = GainCurve_Large;
    }
    
    if( GainOutput > 0){
      std::cout << "TEST" << std::endl;
      Ratio = DestinationGain/GainOutput;
    }else{
      Ratio = 0; 
    }
    
    std::cout << GainOutput      << "\t" 
	      << DestinationGain << "\t"
	      << Ratio           << std::endl;
    

    
    VoltageFinetune[i] = GainCurve->GetX(Ratio * (GainCurve->Eval(CH_HV_MAP[i])),700,5000);
    EndAdjGain         = GainOutput*GainCurve->Eval(VoltageFinetune[i])/GainCurve->Eval(CH_HV_MAP[i]);
    if( VoltageFinetune[i] < HVLowLimit ){
      VoltageFinetune[i] = HVLowLimit;
      EndAdjGain         = GainOutput*GainCurve->Eval(HVLowLimit)/GainCurve->Eval(CH_HV_MAP[i]);
    }else if( VoltageFinetune[i] > HVHighLimit ){
      VoltageFinetune[i] = HVHighLimit;
      EndAdjGain         = GainOutput*GainCurve->Eval(HVHighLimit)/GainCurve->Eval(CH_HV_MAP[i]);
    }

    ofs<< i                  << "\t" 
       << (int)VoltageFinetune[i] << "\t" 
       << GainOutput         << "\t" 
       << EndAdjGain         << std::endl;
    
    imageGain[0]->Fill(i, GainOutput);
    imageGain[1]->Fill(i, EndAdjGain);
    if( GainOutput > 0){
      imageGain[2]->Fill(i, EndAdjGain/GainOutput);
    }
    imageVolt[0]->Fill(i,CH_HV_MAP[i]);
    imageVolt[1]->Fill(i,VoltageFinetune[i]);
    if( VoltageInitial[i] >0){
      imageVolt[2]->Fill(i,VoltageFinetune[i]/CH_HV_MAP[i]);
    }
    
    std::cout<< GainCurve->Eval(HVHighLimit) << std::endl;
        
    hisGainRaw->Fill(GainOutput);
    hisGainCal->Fill(EndAdjGain);
    hisVoltageFine->Fill(VoltageFinetune[i]);    
    hisVoltageInit->Fill(CH_HV_MAP[i]);
  }
  
  TCanvas* can = new TCanvas("can","",1600,800);
  can->Divide(4,2);
  
  can->cd(1);
  imageVolt[0]->DrawWithRange("colz",1000,1750);
  can->cd(2);
  imageVolt[1]->DrawWithRange("colz",1000,1750);
  can->cd(3);
  imageVolt[2]->DrawWithRange("colz",0.4,1.6);
  can->cd(4);

  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();

  hisVoltageInit->GetYaxis()->SetRangeUser(0.01,2000);
  hisVoltageInit->Draw();
  hisVoltageFine->SetLineColor(2);
  hisVoltageFine->Draw("same");
  can->cd(5);
  imageGain[0]->DrawWithRange("colz",0,4000);
  can->cd(6);
  imageGain[1]->DrawWithRange("colz",0,4000);
  can->cd(7);
  imageGain[2]->DrawWithRange("colz",0,3);
  can->cd(8);

  gPad->SetLogy();
  gPad->SetGridx();
  gPad->SetGridy();

  hisGainRaw->GetYaxis()->SetRangeUser(0.1,1000);  
  hisGainRaw->Draw();
  hisGainCal->SetLineColor(2);
  hisGainCal->Draw("same");

  app->Run();
}



