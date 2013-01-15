#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

void GetRMS(){
  std::string str = std::getenv("ANALYSISLIB");
  str += "/lib/libAnalysisLib.so";
  std::cout<< str << std::endl;
  gSystem->Load(str.c_str());
  IDHandler* handler = new IDHandler();
  CsIImage* image= new CsIImage(handler); 
  
  
  TFile* tf = new TFile("Calibration_Data/CalhistListADV_KLRunList_2.txt_0.root");
  //TFile* tf = new TFile("/home/had/jwlee/workdir/RootFiles/3pi0Run/SimCalibration_Data/CalhistList_0.root");
  
  TH1D* his[2716];
  TH1D* hisFactorSigma = new TH1D("hisFactorSigma","hisFactorSigma",100,0,0.1);
  for( int i = 0; i< 2716; ++i ){
    his[i] = (TH1D*)tf->Get(Form("hisCalibrationFactor_%d",i));
    if( his[i]->GetEntries() != 0){
      hisFactorSigma->Fill(his[i]->GetRMS());
      image->Fill(i, his[i]->GetRMS());
    }
  }
  
  TCanvas* can = new TCanvas("can","",800,800);
  //image->DrawWithRange("colz",0.04,0.1);
  hisFactorSigma->Draw();
}
