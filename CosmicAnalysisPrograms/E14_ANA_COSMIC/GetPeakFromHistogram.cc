#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSpectrum.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"


//void GetPeakFromHistogram(int runNumber){
int main( int argc, char** argv){
  if( argc != 3){
    return -1;
  }

  std::string ifilename = argv[1];
  std::string ofilename = argv[2];

  TFile* tf= new TFile( ifilename.c_str());

  

  TCanvas* can = new TCanvas("can","",800,800);
  can->Draw();
  TFile* tf = new TFile(filename);
  if( ! tf->IsOpen() ){
    return -1;
  }
  
  TFile* otf = new TFile(Form("%s/CosmicOutput.root",path),"update");
  TGraphErrors* gr = new TGraphErrors();
  gr->SetNameTitle(Form("grCosmicOutput%04d",runNumber),Form("CosmicOutput%04d",runNumber));
  
  for( int i = 0 ; i< 2716; i++){
    std::cout << "ID:" << i << std::endl;
    TH1D*  hisCosmic = (TH1D*)tf->Get(Form("%s%04d",hisname,i));
    TF1*   fitFunction = NULL;
    Double_t Peak;
    Double_t Sigma; 
    if( hisCosmic->GetEntries() <= 10 ){
      gr->SetPoint(gr->GetN(), i, 0);
      gr->SetPointError(gr->GetN()-1 ,0 ,0);
      continue;
    }
    
    TSpectrum* spec  = new TSpectrum(10,10);
    Int_t nPeak      = spec->Search(hisCosmic,10,"",0.2);
    if( nPeak != 0 ){
      Float_t* xpeaks  = spec->GetPositionX();
      Float_t* ypeaks  = spec->GetPositionY();
      
      Float_t PeakX   =0;
      Float_t lowerX  =0;
      Float_t CurrentX=0;
      Float_t PeakY   =0;
      Float_t lowerY  =0;
      Float_t CurrentY=0;
      
      for( int ipeak = 0; ipeak <nPeak; ipeak++){
	CurrentY = ypeaks[ipeak];
	CurrentX = xpeaks[ipeak];
	if( CurrentY > lowerY ){
	  lowerY = CurrentY;
	  lowerX = CurrentX;
	  if( CurrentX > 500 ){
	    PeakX = CurrentX;
	    PeakY = CurrentY;
	  }
	}    
	
      }
      if(PeakX == 0){
	PeakX = lowerX;
	PeakY = lowerY;
      }
      
      std::cout << "PEAK" << nPeak << std::endl;    
      /*
	hisCosmic->Fit("gaus","","",
	300,
	hisCosmic->GetMean()+4*hisCosmic->GetRMS());
      */
      hisCosmic->Fit("gaus","","",
		     PeakX-2*hisCosmic->GetRMS(),
		     PeakX+2*hisCosmic->GetRMS());
      fitFunction = hisCosmic->GetFunction("gaus");
      std::cout <<fitFunction->GetParameter(1) << std::endl;
      Peak = fitFunction->GetParameter(1);
      Sigma= fitFunction->GetParameter(2);
      
      gr->SetPoint(gr->GetN(), i, fitFunction->GetParameter(1));
      gr->SetPointError(gr->GetN()-1,0, fitFunction->GetParError(1));
    }else{
      gr->SetPoint(gr->GetN(), i, 0);
      gr->SetPointError(gr->GetN()-1,0,0);
      
    }
    //can->Update();
    //can->Modified();
    //can->Update();
    //getchar();
  }						
  gr->GetYaxis()->SetRangeUser(0,5000);
  gr->SetMarkerStyle(3);
  gr->Draw("AP");
  std::cout << gr->GetN() << std::endl;
  gr->Write();
  otf->Close();
  //gr->Draw("AP");
}


  
