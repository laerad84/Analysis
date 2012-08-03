#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TSpline.h"

#include "TCanvas.h"
#include "TApplication.h"

int
main( int argc, char** argv ){
  std::string  ROOTFILE_COSMIC = std::getenv("ROOTFILE_COSMIC");    
  TFile* tf_CalFactor = new TFile(Form("%s/CosmicResult_20120209.root",
				       ROOTFILE_COSMIC.c_str()));
  TTree* tr_CalFactor = (TTree*)tf_CalFactor->Get("GainFitPar");
  Double_t Norm;
  Double_t Peak;
  Double_t Sigma;
  Double_t Chi2;
  Int_t    Ndf;
  Int_t    ID;

  tr_CalFactor->SetBranchAddress("Norm",&Norm);
  tr_CalFactor->SetBranchAddress("Peak",&Peak);
  tr_CalFactor->SetBranchAddress("Sigma",&Sigma);
  tr_CalFactor->SetBranchAddress("Ndf",&Ndf);
  tr_CalFactor->SetBranchAddress("ID",&ID);
  
  Double_t CalFactor[2716] = {0};
  for( int iID = 0; iID< tr_CalFactor->GetEntries(); iID++){
    tr_CalFactor->GetEntry(iID);
    if( Peak > 0 ){
      CalFactor[ID] = 14./Peak; // MeV / Signal Height;; 
    }
  }
  tf_CalFactor->Close();

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  TGraphErrors* gr_waveform_Height[2716];
  TGraphErrors* gr_waveform_Cosmic[2716];
  TFile* tfwf_height = new TFile("TEMPLETE_OUT_HEIGHT_0.root");
  TFile* tfwf_cosmic = new TFile("TEMPLETE_OUT_COSMIC.root");
  for( int ich = 0; ich < 2716; ich ++){
    gr_waveform_Height[ich] = (TGraphErrors*)tfwf_height->Get(Form("Waveform_Height_%d_0",ich));
    gr_waveform_Cosmic[ich] = (TGraphErrors*)tfwf_cosmic->Get(Form("Waveform_Cosmic_%d",ich));
  }
  return 0;
}



