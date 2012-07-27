#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>

#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TPDF.h"
#include "TPostScript.h"
#include "TApplication.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
int main( int argc, char** argv ){
  
  
  std::string InputFilename = argv[1];
  TApplication* app = new TApplication( "app", &argc, argv );
  
  TFile * inputFile = new TFile(InputFilename.c_str());
  TTree * inputTree = (TTree*)inputFile->Get("WFTree");
  Int_t Data[48]; 
  Int_t Time[48]; 
  Int_t ID;
  Int_t EventNumber;
  inputTree->SetBranchAddress("Data",Data);
  inputTree->SetBranchAddress("Time",Time);
  inputTree->SetBranchAddress("ID",&ID);
  inputTree->SetBranchAddress("EventNumber",&EventNumber);
  
  TCanvas* can = new TCanvas("can","",800,800);
  TGraph* gr  = new TGraph();
  gr->SetMarkerStyle(6);
  TGraph* grHeightRMS = new TGraph();
  for( int i = 0; i< inputTree->GetEntries(); i++){
    //for( int i = 0; i< 5000; i++){
    gr->Set(0);
    inputTree->GetEntry(i);
    Double_t max = 0;
    Double_t min = 18000;
    Double_t Sigsq = 0;
    Double_t PeakTime = 0; 
    for( int index = 0; index < 48; index++ ){
      gr->SetPoint( index , Time[index], Data[index] );
      if( Data[index] < min ){ min = Data[index] ;}
      if( Data[index] > max ){ max = Data[index] ; PeakTime = Time[index]; }      
    }
    for( int index = 0; index < 48; index++){
      Sigsq += ( Data[index] - min )*( Data[index] - min );
    }
    
    grHeightRMS->SetPoint( grHeightRMS->GetN(), PeakTime, TMath::Sqrt( Sigsq )/( max-min ));
    gr->Draw("AP");
    //std::cout << gr->GetRMS(2) << std::endl;
    //gr->GetYaxis()->SetRangeUser(0.1, 16000);    
    //gPad->SetLogy();
    if( TMath::Sqrt( Sigsq ) / ( max -min )  > 2.8 ){ continue; } 
    if( max - min < 20 ){continue; }
    can->Update();
    can->Modified();
    getchar();         
   
  }
  app->Run();
}
