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
  
  IDHandler* handler      = new IDHandler("Data/crystal.txt");
  CsIImage*  imagePeak    = new CsIImage(handler);
  CsIImage*  imageHit     = new CsIImage(handler);
  CsIImage*  imageNorm    = new CsIImage(handler);
  CsIImage*  imageLowCH   = new CsIImage(handler);
  imagePeak->SetTitle("Distribution On image; x[mm]; y[mm]");
  imageHit->SetTitle("Number Entries; x[mm];y[mm]");
  imageNorm->SetTitle("Fit Function Height/Entries; x[mm]; y[mm]");

  TH1D* gainHist      = new TH1D("gain","Distribution of gain;Integrated ADC;N/10cnt",500,0,5000);				 
  
  std::ofstream ofs(OutFilename.c_str());  
  
  Double_t gain;
  Double_t Sigma;
  Double_t Norm;
  Double_t HistEntries;
  Double_t Chi2;
  Double_t Ndf;

  std::cout<< InputFilename << std::endl;
  TFile* tf = new TFile(InputFilename.c_str());

  for( int i = 0; i< 2716; i++){    
    TH1D* hist  = (TH1D*)tf->Get(Form("his_CH%04d",i));
    TF1*  func  = NULL;
    func        = hist->GetFunction("landau");
    HistEntries = hist->GetEntries();
    
    if( func == NULL ){
      continue;
    }
    if( i < 2240){      
      Norm = func->GetParameter(0);
      gain = func->GetParameter(1);      
      Sigma= func->GetParameter(2);
      Chi2 = func->GetChisquare();
      Ndf  = func->GetNDF();
    }else{
      Norm = func->GetParameter(0);
      gain = func->GetParameter(1)/2;      
      Sigma= func->GetParameter(2)/2;
      Chi2 = func->GetChisquare();
      Ndf  = func->GetNDF();
    }
    
    if( gain>300 )
      imagePeak->Fill(i, gain);
    if( hist->GetEntries() > 0)
      imageHit->Fill(i, hist->GetEntries());
    if( Norm >0)
      imageNorm->Fill(i, Norm);
    if( gain< 700 && gain>0 ){
      std::cout<< i << " \t" << gain << "\t" << Chi2/Ndf << std::endl;
      imageLowCH->Fill(i, gain);
    }
    ofs << i     << "\t" 
	<< gain  << "\t"
	<< Sigma << std::endl;             
    gainHist->Fill(gain);
    
  }
  
  ofs.close();
  TCanvas* can = new TCanvas("can","",1500,1000);
  can->Divide(3,2);  
  can->cd(1);
  imagePeak->Draw();
  can->cd(2);
  imageHit->Draw();
  can->cd(3);
  imageNorm->Draw();
  can->cd(4);
  gainHist->Draw();
  can->cd(5);
  imageLowCH->Draw();
  app->Run();
}
