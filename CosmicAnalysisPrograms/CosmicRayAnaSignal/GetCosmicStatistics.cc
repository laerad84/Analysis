#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <libgen.h>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"


bool searchPeak(TH1D* hisCosmic, Double_t& Norm, Double_t& Peak, Double_t& Sigma, Double_t& chi2, Int_t& ndf){
  TSpectrum* spec    = new TSpectrum(10,10);
  Int_t nPeak        = spec->Search(hisCosmic,10,"",0.2);
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
      if( CurrentX > 4 && CurrentX > hisCosmic->GetRMS()){
	if( CurrentY > lowerY ){
	  lowerY = CurrentY;
	  lowerX = CurrentX;
	}
      }          
    }
    
    PeakX = lowerX;
    PeakY = lowerY;
    /*
    if(PeakX < 300){    
    Peak  = -1;
    Sigma = -1;
    Norm  = -1;
      chi2  = -1;
      ndf   = -1;
      return false;
    }
    */
    //TF1* fitFunction = new TF1("fitFunc","landau");
    Double_t LowLimit = PeakX*0.8;
    Double_t HighLimit= PeakX*1.4;
    if( 2*hisCosmic->GetRMS() > PeakX*0.4){
      HighLimit = PeakX + 2.0*hisCosmic->GetRMS();
    }

    TFitResultPtr result = hisCosmic->Fit("landau","Q","",
					  LowLimit,
					  HighLimit);
    int fitresult = result;
    if( fitresult == -1 ){
      Peak  = -1;
      Sigma = -1;
      Norm  = -1;
      chi2  = -1;
      ndf   = -1;
      return false;
    }else{      
      TF1* fitFunction = hisCosmic->GetFunction("landau");
      Norm = fitFunction->GetParameter(0);
      Peak = fitFunction->GetParameter(1);
      Sigma= fitFunction->GetParameter(2);
      chi2 = fitFunction->GetChisquare();
      ndf  = fitFunction->GetNDF();
      if( Peak <  PeakX-0.5*hisCosmic->GetRMS()){
	
	Peak  = -1;
	Sigma = -1;
	Norm  = -1;
	chi2  = -1;
	ndf   = -1;
	
	return false;
      }
      return true;
    }
  }else{
    Peak  = -1; 
    Sigma = -1;
    Norm  = -1;
    chi2  = -1;
    ndf   = -1;
    return false;
  }
}


int 
main( int argc, char** argv){

  std::string inputFileList;
  std::string ROOTFILE_COSMIC = std::getenv("ROOTFILE_COSMIC");
  TChain* trCosmic = new TChain("CosmicOut");  

  Int_t RunIDStart;
  Int_t RunIDEnd;
  std::vector<int> RunList;
  if( argc ==2 ){
    inputFileList = argv[1];
    int tmpRunList;
    std::ifstream ifs( inputFileList.c_str());
    while( ifs >> tmpRunList ){
      RunList.push_back(tmpRunList);
      trCosmic->Add(Form("%s/run%d_cosmic.root",ROOTFILE_COSMIC.c_str(),tmpRunList));
    }

  }else if( argc ==3){
    RunIDStart= atoi( argv[1]);
    RunIDEnd = atoi( argv[2]);
    for( int i = RunIDStart; i <= RunIDEnd; i++){
      trCosmic->Add(Form("%s/run%d_cosmic.root",ROOTFILE_COSMIC.c_str(),i));
    }
  }else{
    return -1; 
  } 
  std::cout<< trCosmic->GetEntries() << std::endl;
}
 
