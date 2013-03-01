#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>

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
#include "IDHandler.h"

//void testVoltage(){
void PrintUsage(char **argv ){
  std::cout << "Usage: "
	    << argv[0] 
	    << " [Destination Gain] [Initial GainMap] [Initital HVMap] \n" 
	    << "Destination Gain : Integer \n"
	    << "Initial GainMap input File Format : [ID] [Gain] [Sigma]\n"
	    << "Initial HVMap   input File Format : [ID] [HV]  \n"
	    << "Adjusted HVMAP output File Format : [ID] [HV]  \n"
	    << std::endl;
}

IDHandler* handler  = new IDHandler("Data/crystal.txt");
Double_t DestinationGain[2]   = {1500,1500};
const Double_t HVLowLimit     = 1200;//V
const Double_t HVHighLimit    = 1750;//V    

void SearchVoltage(int ID,  double  initialGain, double initialVoltage, double& FinalGain, double&FinalVoltage);

TF1* GainCurve_Small = new TF1("GainCurve_Small_raw","TMath::Power(x-657.8,0.9635)",0,5000);
TF1* GainCurve_Large = new TF1("GainCurve_Large_raw","TMath::Power(x-627.9,1.369)",0,5000);

int
main( int argc, char** argv){
  //Double_t DestinationGain   = 1500.0;

  std::string Initial_HVMAP;
  std::string Initial_Gain_MAP;
  std::string Final_HVMAP;
  std::string Final_Gain_MAP;
  
  if( argc == 4){
    DestinationGain = atof(argv[1]);    
  }else{ 
    PrintUsage(argv);
    return -1;
  }
  
  Initial_Gain_MAP = argv[2];
  Initial_HVMAP    = argv[3];
  Final_Gain_MAP   = Form("%s_AfterAdjusted.txt",
			  Initial_Gain_MAP.substr(0,Initial_Gain_MAP.size()-4).c_str());
  Final_HVMAP      = Form("%s_AfterAdjusted.txt",
			  Initial_HVMAP.substr(0,Initial_HVMAP.size()-4).c_str());

  std::cout << Final_Gain_MAP << std::endl;
  std::cout << Final_HVMAP    << std::endl;

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
  IDHandler* handler    = new IDHandler("Data/crystal.txt");
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
  
  std::cout <<"ReadFile" << std::endl;
  
  Int_t    CH;
  Double_t HVSet;
  Double_t Gain;
  Double_t Sigma;

  Double_t InitialVoltage[2716] = {0};
  Double_t FinalVoltage[2716]   = {0};
  Double_t InitialGain[2716]     = {0};
  Double_t FinalGain[2716]      = {0};
  
  std::ifstream ifsGain(Initial_Gain_MAP.c_str());
  std::ifstream ifsHV(Initial_HVMAP.c_str());
  std::ofstream ofsGain(Final_Gain_MAP.c_str());
  std::ofstream ofsHV(Final_HVMAP.c_str());
  
  while( !ifsGain.eof() ){
    ifsGain >> CH >> Gain >> Sigma;
    if( ifsGain.eof() ){ break;}
    InitialGain[CH] = Gain;
  }
  while( !ifsHV.eof() ) {
    ifsHV >> CH >> HVSet;
    if( ifsHV.eof() ){ break; }
    InitialVoltage[CH] = HVSet;
  }
  
  
  TH1D* hisGainRaw = new TH1D("hisGainRaw","hisGainRaw",120,0,6000);
  TH1D* hisGainCal = new TH1D("hisGainCal","hisGainCal",120,0,6000);
  TH1D* hisVoltageInit = new TH1D("hisVoltageInit","hisVoltageInit",180,0,1800);
  TH1D* hisVoltageFine = new TH1D("hisVoltageFine","hisVoltageFine",180,0,1800);
  
  
  int nBadCH = 0;
  for( int i = 0; i< 2716; i++){
    if( i >=2240     && i< 2240+120 ){continue;}
    if( i >=2716-120 && i< 2716     ){continue;}
    
    SearchVoltage(i, InitialGain[i], InitialVoltage[i], FinalGain[i], FinalVoltage[i]);
    std::cout<< "CHANNEL:" <<  i  << "\t"
	     << InitialVoltage[i] << "\t"
	     << InitialGain[i]    << "\t"
	     << FinalVoltage[i]   << "\t" 
	     << FinalGain[i]      << "\t"
	     << std::endl;
    ofsHV   << i << "\t" << (int)(FinalVoltage[i]) << std::endl;
    ofsGain << i << "\t" << FinalGain[i]           << std::endl;   
    
    imageGain[0]->Fill(i, InitialGain[i]);
    imageGain[1]->Fill(i, FinalGain[i]);
    if( InitialGain[i] > 0){
      imageGain[2]->Fill(i, FinalGain[i]/InitialGain[i]);
    }
    imageVolt[0]->Fill(i,InitialVoltage[i]);
    imageVolt[1]->Fill(i,FinalVoltage[i]);
    if( InitialVoltage[i] >0){
      imageVolt[2]->Fill(i,FinalVoltage[i]/InitialVoltage[i]);
    }
    
    hisGainRaw->Fill(InitialGain[i]);
    hisGainCal->Fill(FinalGain[i]);
    hisVoltageFine->Fill(FinalGain[i]);    
    hisVoltageInit->Fill(InitialGain[i]);
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
  ofsGain.close();
  ofsHV.close();
  app->Run();
}

void SearchVoltage( int ID, double  initialGain, double initialVoltage, double& finalGain, double&finalVoltage){
  TF1*   GainCurve;
  double Ratio;
  if( ID <2240 ){
    GainCurve = GainCurve_Small;
  }else{
    GainCurve = GainCurve_Large;
  }
  if( initialGain == 0){
    finalGain    = 0;
    finalVoltage = 0;
    return;
  } 
  if( initialVoltage == 0){
    finalGain    = 0;
    finalVoltage = 0;
    return;    
  }
  
  if( initialGain > 0){
    std::cout << "TEST" << std::endl;
    Ratio = DestinationGain/initialGain;
  }else{
    Ratio = 0; 
  } 
 
  finalVoltage = GainCurve->GetX(Ratio * (GainCurve->Eval(initialVoltage)),700,5000);
  finalGain    = initialGain*GainCurve->Eval(finalVoltage)/GainCurve->Eval(initialVoltage);

  if( finalVoltage < HVLowLimit ){
    finalVoltage = HVLowLimit;
    finalGain    = initialGain*GainCurve->Eval(HVLowLimit)/GainCurve->Eval(initialVoltage);
  }else if( finalVoltage > HVHighLimit ){
    finalVoltage = HVHighLimit;
    finalGain    = initialGain*GainCurve->Eval(HVHighLimit)/GainCurve->Eval(initialVoltage);
  }
  return;
}




