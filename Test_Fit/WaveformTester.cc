#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "E14WaveformAnalyzer.h"

int main( int argc ,char** argv ){  
  if( argc != 2 ){
    std::cerr << "Arguement Error " << std::endl;  
    return -1; 
  } 

  //Int_t RunNumber = 4205;
  Int_t RunNumber  = atoi(argv[1]);
  std::cout<< RunNumber << std::endl;
  std::string ROOTFILEWAV= std::getenv("ROOTFILE_WAV");
  std::string Filename   = Form("%s/TEMPLATE_FIT_RESULT_DUMP_%d.root",ROOTFILEWAV.c_str(),RunNumber);

  std::cout<< __LINE__ << std::endl;
  
  TApplication* app = new TApplication( "app" , &argc , argv );
  TCanvas* can = new TCanvas("can","",800,400);
  
  
#ifdef DEBUG
  std::cout<< "DEFINE INPUTFILE" << std::endl; 
#endif
  std::cout<< __LINE__ << std::endl;

  TFile* tf = new TFile(Filename.c_str());
  TTree* tr = (TTree*)tf->Get("Waveform");
  
  Int_t    EventNumber;
  Int_t    ModuleNumber;
  Double_t Waveform[48];
  Short_t  TimeInfo[48];
  Double_t ChisqNDF;
  Double_t PeakTime;
  Double_t HHTime;
  Double_t Height;
  Double_t Pedestal;
  

  tr->SetBranchAddress("EventNumber",&EventNumber);
  tr->SetBranchAddress("ModuleNumber",&ModuleNumber);
  tr->SetBranchAddress("Waveform",Waveform);
  tr->SetBranchAddress("TimeInfo",TimeInfo);
  tr->SetBranchAddress("PeakTime",&PeakTime);
  tr->SetBranchAddress("HHTime",&HHTime);
  tr->SetBranchAddress("Height",&Height);
  tr->SetBranchAddress("Pedestal",&Pedestal);
  
  std::cout<< __LINE__ << std::endl;

  TFile* tfOut           = new TFile(Form("OutputDistrib_Tester_%d.root",RunNumber),"RECREATE");
  std::cout<< __LINE__ << std::endl;
  E14WaveformAnalyzer*  WaveAnalyzer = new E14WaveformAnalyzer(48);
  std::cout<< __LINE__ << std::endl;

  Int_t TotalEvent = tr->GetEntries(); 
  TGraph* grWave = new TGraph();

  can->Divide( 2,1);
  for( int ievent  =0 ; ievent < TotalEvent ; ievent++){
  //for( int ievent  =0 ; ievent < 100 ; ievent++){
    tr->GetEntry(ievent);    
    grWave->Set(0);
    std::cout<< __LINE__ << std::endl;
    can->cd(1);
    WaveAnalyzer->_Draw( Waveform );
    can->cd(2);
    WaveAnalyzer->m_peakGraph->SetMarkerStyle(6);
    WaveAnalyzer->m_peakGraph->Draw("AP");
    WaveAnalyzer->m_peakFunc->Draw("same");
    std::cout<< __LINE__ << std::endl;
    can->Update();
    can->Modified();
    getchar();
    WaveAnalyzer->_Clear();
  }

  tfOut->Close();
  app->Run();
  return 0;
}

