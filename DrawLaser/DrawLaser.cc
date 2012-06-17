#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

#include "E14ReadSumFile.h"
#include "IDHandler.h"

int
main( int argc, char **argv){
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("neRMIuo");
  if( argc  != 2 ){
    std::cerr << "Please Input RunNumber" << std::endl;
    return -1; 
  }
    
  int RunNumber = atoi(argv[1]);  

  TFile* tfile = new TFile(Form("LaserOutRootData/LaserOutMeanData%d.root",RunNumber),"recreate");
  TTree* trLaserOut = new TTree("LaserTree","Laser Out");
			   
  double  CsISigma[N_TOTAL_CSI];
  double  CsIOut[N_TOTAL_CSI];
  double  RatioSigma[N_TOTAL_CSI];
  double  RatioOut[N_TOTAL_CSI];
  double  LaserOut[N_TOTAL_CSI];
  double  LaserSigma[N_TOTAL_CSI];
  double  CsIErr[N_TOTAL_CSI];
  double  RatioErr[N_TOTAL_CSI];
  double  LaserErr[N_TOTAL_CSI];

  trLaserOut->Branch("RunNumber" ,&RunNumber ,"RunNumber/I");
  trLaserOut->Branch("CsIOut"    ,CsIOut     ,Form("CsIOut[%d]/D",N_TOTAL_CSI));
  trLaserOut->Branch("CsISigma"  ,CsISigma   ,Form("CsISigma[%d]/D",N_TOTAL_CSI));
  trLaserOut->Branch("CsIErr"    ,CsIErr     ,Form("CsIErr[%d]/D",N_TOTAL_CSI));
  trLaserOut->Branch("RatioOut"  ,RatioOut   ,Form("RatioOut[%d]/D",N_TOTAL_CSI));
  trLaserOut->Branch("RatioSigma",RatioSigma ,Form("RatioSigma[%d]/D",N_TOTAL_CSI));
  trLaserOut->Branch("RatioErr"  ,RatioErr   ,Form("RatioErr[%d]/D",N_TOTAL_CSI));
  trLaserOut->Branch("LaserOut"  ,LaserOut   ,Form("LaserOut[%d]/D",N_TOTAL_CSI));
  trLaserOut->Branch("LaserSigma",LaserSigma ,Form("LaserSigma[%d]/D",N_TOTAL_CSI));
  trLaserOut->Branch("LaserErr"  ,LaserErr   ,Form("LaserErr[%d]/D",N_TOTAL_CSI));

  TFile* tf = new TFile(Form("LaserMergeData/LaserMerged%d.root",RunNumber));
  TTree* tr = (TTree*)tf->Get("LaserOut");
  int nDigi;
  double SumOut[N_TOTAL_CSI];
  double PINOut[N_TOTAL_CSI];
  int ID[N_TOTAL_CSI];
  tr->SetBranchAddress("nDigi",&nDigi);
  tr->SetBranchAddress("SumOut",SumOut);//nDigi
  tr->SetBranchAddress("PINOut",PINOut);//nDigi
  tr->SetBranchAddress("ID",ID);//nDigi


  TH1D* hisCsI[N_TOTAL_CSI];
  TH1D* hisRatio[N_TOTAL_CSI];
  TH1D* hisLaser[N_TOTAL_CSI];
  for( int ich = 0; ich < N_TOTAL_CSI; ++ich){
    hisCsI[ich]   = new TH1D(Form("hisCsI%d",ich),Form("CsI%d",ich),2000,0,200000);
    hisRatio[ich] = new TH1D(Form("hisRatio%d",ich),Form("Ratio%d",ich),2000,0,100);
    hisLaser[ich] = new TH1D(Form("hisLaser%d",ich),Form("Laser%d",ich),2000,0,500000);
  }  

  
  for( int ientry = 0; ientry < tr->GetEntries(); ientry++){
    if( ientry %100 ==0 ){
      std::cout << ientry <<std::endl;
    }
    tr->GetEntry(ientry);
    for( int idigi = 0; idigi < nDigi; ++idigi){
      int idofCsI = ID[idigi];
      hisCsI[ idofCsI ]->Fill(SumOut[idigi]);
      hisLaser[ idofCsI ]->Fill(PINOut[idigi]);
      hisRatio[ idofCsI ]->Fill(SumOut[idigi]/PINOut[idigi]);
    }
  }


  for( int ich = 0; ich < N_TOTAL_CSI; ++ich){

    /*
    tr->Project(hisCsI[ich]->GetName(),"SumOut",Form("ID==%d"));
    tr->Project(hisRatio[ich]->GetName(),"SumOut/PINOut",Form("ID==%d"));
    tr->Project(hisLaser[ich]->GetName(),"PINOut",Form("ID==%d"));
    */

    CsIOut[ich]   = hisCsI[ich]->GetMean();
    CsIErr[ich]   = hisCsI[ich]->GetMeanError();
    CsISigma[ich] = hisCsI[ich]->GetRMS();
    LaserOut[ich] = hisLaser[ich]->GetMean();
    LaserErr[ich] = hisLaser[ich]->GetMeanError();
    LaserSigma[ich] = hisLaser[ich]->GetRMS();
    RatioOut[ich] = hisRatio[ich]->GetMean();
    RatioErr[ich] = hisRatio[ich]->GetMeanError();
    RatioSigma[ich] = hisRatio[ich]->GetRMS();

  }	     	       
  trLaserOut->Fill();
  tfile->cd();
  trLaserOut->Write();
  tfile->Close();
}

  

  



  


