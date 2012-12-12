#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "E14WavReader.h"
#include <cstring>
#include <string>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "User_Function.h"
#include "EnergyConverter.h"

#include "ClusterFinder_EDIT.h"
#include "IDHandler.h"
#include "CsIImage.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "HoughCsI.h"
#include "Chisq_cosmic.h"
#include "E14CosmicAnalyzer.h"

#include "TChain.h"


int main( int argc, char** argv){

    std::string RunNumberList = argv[1];
    gStyle->SetOptFit(11111111);
    std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
    std::string ANALIBDIR    = std::getenv("ANALYSISLIB");
    
    TChain* trin = new TChain("trOut");	
    std::ifstream ifs(RunNumberList.c_str());	
    int tmpRunNumber;
    while( ifs >>tmpRunNumber ){
      trin->Add(Form("%s/CosmicOut_Converted_%d.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
    }
    const int nCSI = 2716;
    Int_t    RunNumber;
    Int_t    EventNumber;
    Double_t ScintiSignal;
    Double_t ScintiHHTime;
    Double_t ScintiTime;
    Int_t    nCSIDigi;
    Double_t CSIDigiE[nCSI];//nCSIDigi
    Double_t CSIDigiTime[nCSI];//nCSIDigi
    Double_t CSIDigiHHTime[nCSI];//nCSIDigi
    Int_t    CSIDigiID[nCSI];//nCSIDigi
    Double_t CSIDigiSignal[nCSI];//nCSIDigi
    Double_t FitP0[2];
    Double_t FitP1[2];
    Double_t FitChisq[2];
    Double_t CSIDigiDeltaT0[nCSI];//nCSIDigi
    Double_t CSIDigiDeltaT1[nCSI];//nCSIDigi
    Int_t    CosmicTrigUp;
    Int_t    CosmicTrigDn;
    Double_t Roh;
    Double_t Theta;
    
    trin->SetBranchAddress( "RunNumber"     , &RunNumber      );
    trin->SetBranchAddress( "EventNumber"   , &EventNumber    );
    trin->SetBranchAddress( "ScintiSignal"  , &ScintiSignal   );
    trin->SetBranchAddress( "ScintiHHTimne" , &ScintiHHTime   );
    trin->SetBranchAddress( "ScintiTime"    , &ScintiTime     );
    trin->SetBranchAddress( "nCSIDigi"      , &nCSIDigi       );
    trin->SetBranchAddress( "CSIDigiE"      , CSIDigiE        );
    trin->SetBranchAddress( "CSIDigiTime"   , CSIDigiTime     );
    trin->SetBranchAddress( "CSIDigiHHTime" , CSIDigiHHTime   );
    trin->SetBranchAddress( "CSIDigiID"     , CSIDigiID       );
    trin->SetBranchAddress( "CSIDigiSignal" , CSIDigiSignal   );
    trin->SetBranchAddress( "CSIDigiDeltaT0", CSIDigiDeltaT0  );
    trin->SetBranchAddress( "CSIDigiDeltaT1", CSIDigiDeltaT1  );
    trin->SetBranchAddress( "FitP0"         , FitP0           );
    trin->SetBranchAddress( "FitP1"         , FitP1           );
    trin->SetBranchAddress( "FitChisq"      , FitChisq        );
    trin->SetBranchAddress( "CosmicTrigUp"  , &CosmicTrigUp   );
    trin->SetBranchAddress( "CosmicTrigDn"  , &CosmicTrigDn   );
    trin->SetBranchAddress( "Roh"           , &Roh            );
    trin->SetBranchAddress( "Theta"         , &Theta          );
    


    TFile* tfout = new TFile("TimeEnergy_CosmicEnergyRange.root","RECREATE");
    TH2D* hisTimeEnergy[2716];
    TProfile* profTimeEnergy[2716];
    TH2D* hisTimeSignal[2716];
    TProfile* profTimeSignal[2716];
    
    for( int i = 0; i< 2716; i++){
      hisTimeEnergy[i] = new TH2D(Form("hisTimeEnergy_%d",i),
				  Form("hisTimeEnergy_%d;Energy[MeV];Time[ns]",i),
				  40,0,40,150,-15,15);
      profTimeEnergy[i] = new TProfile(Form("profTimeEnergy_%d",i),
				   Form("profTimeEnergy_%d;Energy[MeV];Time[ns]",i),
				   40,0,40);
      
      hisTimeSignal[i] = new TH2D(Form("hisTimeSignal_%d",i),
				  Form("hisTimeSignal_%d;Signal[MeV];Time[ns]",i),
				  40,0,400,150,-15,15);
      profTimeSignal[i] = new TProfile(Form("profTimeSignal_%d",i),
				   Form("profTimeSignal_%d;Signal[MeV];Time[ns]",i),
				   40,0,400);
      
    }
    
    for( int ievent = 0; ievent < trin->GetEntries() ; ievent++){
      trin->GetEntry(ievent);
      if( nCSIDigi > 200 ){ continue; }
      if( FitChisq[2] > 5){ continue; }
      for( int idigi = 0; idigi < nCSIDigi; idigi++){

	hisTimeEnergy[CSIDigiID[idigi]]->Fill(CSIDigiE[idigi],CSIDigiDeltaT1[idigi]);
	profTimeEnergy[CSIDigiID[idigi]]->Fill(CSIDigiE[idigi],CSIDigiDeltaT1[idigi]);
	hisTimeSignal[CSIDigiID[idigi]]->Fill(CSIDigiSignal[idigi],CSIDigiDeltaT1[idigi]);
	profTimeSignal[CSIDigiID[idigi]]->Fill(CSIDigiSignal[idigi],CSIDigiDeltaT1[idigi]);
      }
    }
    
    for( int i = 0; i< 2716; i++){
      hisTimeEnergy[i]->Write();
      hisTimeSignal[i]->Write();
      profTimeEnergy[i]->Write();
      profTimeSignal[i]->Write();
    }
    
    tfout->Close();
}
