#include <fstream>
#include <iostream>
#include <cstring>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"

int
main( int argc, char** argv){
  if( argc  != 2 ) {
    std::cerr << "input Filename " << std::endl;
    return -1;
  }
  std::string InputFilename  = argv[1];
  std::string OutFilename    = Form("%s.txt",InputFilename.substr(0,InputFilename.size()-5).c_str());
  TApplication* app = new TApplication("app",&argc, argv);
  //gSystem->Load("lib/libtest.so");
  
  IDHandler* handler      = new IDHandler();
  CsIImage*  imagePeak    = new CsIImage(handler);
  CsIImage*  imageHit     = new CsIImage(handler);
  CsIImage*  imageNorm    = new CsIImage(handler);
  CsIImage*  imageLowCH   = new CsIImage(handler);
  imagePeak->SetTitle("Distribution On image; x[mm]; y[mm]");
  imageHit->SetTitle("Number Entries; x[mm];y[mm]");
  imageNorm->SetTitle("Fit Function Height/Entries; x[mm]; y[mm]");
  
  TH1D* hisWidth      = new TH1D("Sigma_gain","Distribution of Sigma/Gain",
				 100,0,0.2);
  TH1D* gainHist      = new TH1D("gain","Distribution of gain;Integrated ADC;N/10cnt",500,0,5000);				 
  TH1D* hisChi2Ndf   = new TH1D("hisChi2ndf","",100,0,20);
  TH1D* hisSigmaPeak = new TH1D("hisSigmaPeak","",100,0,0.5);

  std::ofstream ofs(OutFilename.c_str());  
  
  Double_t gain;
  Double_t Norm;
  Double_t HistEntries;
  std::cout<< InputFilename << std::endl;
  TFile* tf = new TFile(InputFilename.c_str());
  TTree* tr = (TTree*)tf->Get("GainFitPar");
  Double_t Peak;
  Double_t Sigma;
  Int_t ID;
  Double_t Chi2;
  Int_t Ndf;
  tr->SetBranchAddress("Chi2" ,&Chi2);
  tr->SetBranchAddress("Ndf"  ,&Ndf);
  tr->SetBranchAddress("Peak" ,&Peak);
  tr->SetBranchAddress("Sigma",&Sigma);
  tr->SetBranchAddress("ID"   ,&ID);


  for( int i = 0; i< tr->GetEntries(); i++){    
    tr->GetEntry(i);
    if( Chi2/Ndf > 4.5  || Sigma/Peak > 0.12){
      std::cout<< ID         << "\t"
	       << Peak       << "\t" 
	       << Sigma      << "\t"
	       << Sigma/Peak << std::endl;
    }
    

    ofs << ID << "\t" << Peak << "\t" << Sigma << std::endl;
    hisSigmaPeak->Fill(Sigma/Peak);
    imagePeak->Fill(ID,Sigma/Peak);
    gainHist->Fill(Peak);
    hisChi2Ndf->Fill(Chi2/Ndf);
  }
  
  ofs.close();
  TCanvas* can = new TCanvas("can","can",800,800);
  can->Divide(2,2);
  can->cd(1);
  gainHist->Draw();
  can->cd(2);
  imagePeak->Draw();
  can->cd(3);
  hisSigmaPeak->Draw();
  can->cd(4);
  hisChi2Ndf->Draw();
  app->Run();
  
}
