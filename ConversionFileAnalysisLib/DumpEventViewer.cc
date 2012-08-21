#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TTreePlayer.h"
#include "TTreePerfStats.h"

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "Structs.h"

#include "E14ConvReader.h"
#include "E14IDHandler.h"
#include "E14ConvWriter.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TText.h"
#include "TROOT.h"
#include "Math/WrappedFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMethods.h"
#include "Math/WrappedTF1.h"

#include "E14WaveFitter.h"

int
main( int argc , char* argv[]){
  //TApplication* app = new TApplication("App",&argc, argv);
  TCanvas* can = new TCanvas("can","", 800,800);
  TFile* tfTemplate = new TFile("TEMPLETE_OUT_HEIGHT_0_OK.root");

  TGraphErrors* tempGr[2716];
  TSpline3*     tempSpl[2716]; 

  for( int ispline = 0; ispline< 2716; ispline++){
    tempGr[ispline]  = NULL;
    tempSpl[ispline] = NULL;
    tempGr[ispline]  = (TGraphErrors*)tfTemplate->Get(Form("Waveform_Height_%d_0",ispline));
    if( tempGr[ispline]->GetN()!= 0){
      tempSpl[ispline] = new TSpline3(Form("waveformSpl_%d",ispline),(TGraph*)(tempGr[ispline]));
    }else{
      //std::cout<< "Non Exist channel:" << ispline << std::endl;
    }
  }
  tfTemplate->Close();

  
  TGraph* gr = new TGraph(); 
  TGraph* gr1 = new TGraph();
  gr->SetMarkerStyle(6);
  gr1->SetMarkerStyle(6);
  
  std::string Env = std::getenv("ROOTFILE_WAV");

  TFile* tf= new TFile(Form("%s/TEMPLATE_FIT_RESULT_DUMP_4503.root",Env.c_str()));  
  TTree* tr = (TTree*)tf->Get("Waveform");
  
  Int_t RunNumber;
  Int_t EventNumber;
  Int_t ModuleNumber;
  Double_t Waveform[48];
  Short_t  TimeInfo[48];
  
  Double_t ChisqNDF;
  Double_t PeakTime;
  Double_t HHTime;
  Double_t Height;
  Double_t Pedestal;
  
  tr->SetBranchAddress("RunNumber",&RunNumber);
  tr->SetBranchAddress("EventNumber",&EventNumber);
  tr->SetBranchAddress("ModuleNumber",&ModuleNumber);
  tr->SetBranchAddress("Waveform",Waveform);
  tr->SetBranchAddress("TimeInfo",TimeInfo);
  tr->SetBranchAddress("ChisqNDF",&ChisqNDF);
  tr->SetBranchAddress("PeakTime",&PeakTime);
  tr->SetBranchAddress("HHTime",&HHTime);
  tr->SetBranchAddress("Height",&Height);
  tr->SetBranchAddress("Pedestal",&Pedestal);

  E14WaveFitter* Fitter     = new E14WaveFitter();

  Int_t nTotal = tr->GetEntries();


  TFile* tempFile = new TFile("tempFile.root","recreate");
  TH2D* hisTestWaveform[5];
  for( int i = 0; i< 5; i++){
    hisTestWaveform[i] = new TH2D(Form("hisTestWaveform%d",i),Form("hisTestWaveform%d",i),
				  500,-250,250,200,-0.5,1.5);
  }


  for( int ievent  = 0 ; ievent < nTotal ; ievent++){
    tr->GetEntry(ievent);    
    if( ChisqNDF < 10 || PeakTime> 325 || PeakTime < 100){
      continue ;
    }
    gr->Set(0);
    gr1->Set(0);
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      if( Waveform[ipoint] >16000 || Waveform[ipoint] ==0){ continue; }
      gr->SetPoint( gr->GetN() , TimeInfo[ipoint], Waveform[ipoint]);
      gr1->SetPoint( gr1->GetN() , TimeInfo[ipoint], Waveform[ipoint]);
    }

    Fitter->InitPar();
    Fitter->SetWaveform( tempSpl[ModuleNumber] );
    //bool fit = Fitter->Fit(gr);

    if( ModuleNumber != 1200 ){ continue; }

    std::cout<< Height << std::endl;
    for( int ipoint  =0; ipoint < 48; ipoint++){

      Double_t time  = TimeInfo[ipoint] - PeakTime;
      Double_t peak  = (Waveform[ipoint] - Pedestal)/Height;
      if( Height> 10 && Height <= 100 ){
	hisTestWaveform[0]->Fill( time , peak );
      }
      if( Height >100 && Height <= 500 ){
	hisTestWaveform[1]->Fill( time , peak );
      }
      if( Height >500 && Height <= 1000 ){
	hisTestWaveform[2]->Fill( time , peak );
      }
      if( Height >1000 && Height <=2000 ){
	hisTestWaveform[3]->Fill( time , peak );
      }
      if( Height > 2000 && Height <=6000 ){
	hisTestWaveform[4]->Fill( time , peak );
      }
    }
    std::cout<< hisTestWaveform[4]->GetEntries() << std::endl;
      /*
    if( fit ){
      std::cout << fit << std::endl ;
      
      //gr->Draw("AP");
      gr1->Draw("AP");
      Fitter->m_FitFunc->Draw("same");
      std::cout << gr->GetN() << std::endl;

      can->Update();
      can->Modified();    
      

      std::cout<< ChisqNDF << "\t" << ChisqNDF/Height << std::endl; 
      
      getchar(); 
      Fitter->Clear();
    }
      */
  }
  for( int i = 0; i< 5; i++){
    hisTestWaveform[i]->Write();
  }
  tempFile->Close();
  //app->Run();
  return 0; 
}
