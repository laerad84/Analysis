#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <vector>

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPDF.h"
#include "TText.h"
#include "TProfile.h"
#include "TMath.h"

#include "IDHandler.h"
#include "CsIImage.h"
#include "E14ReadSumFile.h"
#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "gamma/GammaFinder.h"



int
main( int argc, char** argv){
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptStat("neMRIuo");
  gStyle->SetOptFit(111111111);
  
  std::cout << "KLong Monitor" << std::endl;  
  std::string RunListFilename = "Calibration_Data/KLRunList_2.txt"; 
  
  int RunNumber;
  int RunNumberFinal;

  std::string RunName="Klong";

  std::vector<int> RunList;
  std::ifstream ifs(RunListFilename.c_str());
  
  if( !ifs.is_open()){return -1;}
  
  while( ifs >> RunNumber){
    RunList.push_back(RunNumber);
  }
  
  RunNumberFinal = RunList[RunList.size()-1];  
  std::cout<< "Number Of Run: " << RunList.size() << std::endl;
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  
  IDHandler* handler = new IDHandler();    
  
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////
  // Define TChain
  TChain* ch;
  ch = new TChain("trCalibration");
  std::vector<int>::iterator it; 
  for( it = RunList.begin();
       it != RunList.end();
       it++){
    ch->Add(Form("Calibration_Data/Calibration_new_GammaEneFixing_10_6_Normalize/CalibrationADV_%04d_15.root",*it));
  }
  
  Int_t eventID;
  TFile* tfOut = new TFile("KLHists.root","RECREATE");
  TH1D* his_COEX = new TH1D("his_COEX","COEX of Klong", 400,-200,200);
  TH1D* his_COEY = new TH1D("his_COEY","COEY of Klong", 400,-200,200);  
  TH1D* his_COEX_ZMass = new TH1D("his_COEX_ZMass","COEX of Klong",400,-200,200);
  TH1D* his_COEY_ZMass = new TH1D("his_COEY_ZMass","COEY of Klong",400,-200,200);
  TH1D* his_COEX_KLCUT = new TH1D("his_COEX_KLCUT","COEX of Klong",400,-200,200);
  TH1D* his_COEY_KLCUT = new TH1D("his_COEY_KLCUT","COEY of Klong",400,-200,200);

 ///////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////////
 //Define data;
  
  E14GNAnaDataContainer data;
  data.setBranchAddress(ch);    
  long nentries = ch->GetEntries();
  
  for( int ievent = 0; ievent< nentries; ievent++){      
    ch->GetEntry(ievent);      
    Klong kl;
    Double_t klcoex=0;
    Double_t klcoey=0;
    
    // Data Analysis 
    
    data.getData(kl);
    for( int pi0Index  = 0; pi0Index < kl.pi0().size(); pi0Index++){
      klcoex += kl.pi0()[pi0Index].g1().x() * kl.pi0()[pi0Index].g1().e();
      klcoex += kl.pi0()[pi0Index].g2().x() * kl.pi0()[pi0Index].g2().e();
      klcoey += kl.pi0()[pi0Index].g1().y() * kl.pi0()[pi0Index].g1().e();
      klcoey += kl.pi0()[pi0Index].g2().y() * kl.pi0()[pi0Index].g2().e();
    }
    klcoex = klcoex/kl.e();
    klcoey = klcoey/kl.e();
    
    if( data.CutCondition == 0){
      his_COEX->Fill( klcoex );
      his_COEY->Fill( klcoey );
    }
    if( (data.CutCondition& (1+8)) == 0){ 
      his_COEX_ZMass->Fill( klcoex );
      his_COEY_ZMass->Fill( klcoey );      
    }
    if( (data.CutCondition& (1+2+4+8)) == 0){ 
      his_COEX_KLCUT->Fill( klcoex );
      his_COEY_KLCUT->Fill( klcoey );      
    }
    
  }  
  his_COEX->Write();
  his_COEY->Write();
  his_COEX_ZMass->Write();
  his_COEY_ZMass->Write();
  his_COEX_KLCUT->Write();
  his_COEY_KLCUT->Write();

  tfOut->Close();

}

