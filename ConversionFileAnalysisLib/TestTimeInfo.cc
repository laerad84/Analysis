#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "TROOT.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

#include "TCanvas.h"
#include "TApplication.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
int main(int argc, char** argv ){

  TApplication* app = new TApplication("app",&argc, argv );
  
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////


  std::string TimeOffsetFile = "/home/jwlee/local/Analysis/CosmicAnalysisPrograms/CosmicRayAnaSignal/testNewWORKCompileOffset.txt";

  std::ifstream ifs(TimeOffsetFile.c_str());
  
  Double_t TimeOffset[2716];
  
  Int_t ID;
  Double_t Offset;
  Double_t RMS;

  if( !ifs.is_open() ){
    std::cerr << "File is not exist." << std::endl; 
    return -1; 
  }
  while ( ifs >>  ID >> Offset >> RMS ){
    TimeOffset[ID] = Offset;
  }
  
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  std::string wavDir = std::getenv("ROOTFILE_WAV");
  std::cout<< wavDir << std::endl;
  wavDir+="/TEMPLATE_FIT_RESULT_4503.root";
  TFile* tf =new TFile(wavDir.c_str());
  TTree* tr = (TTree*)tf->Get("WFTree");
  
  Int_t    CsinDigi;
  Double_t CsiSignal[2716];
  Double_t CsiTime[2716];
  Double_t CsiHHTime[2716];
  Double_t CsiChisq[2716];
  Short_t  CsiNDF[2716];
  Short_t  CsiID[2716];

  Int_t    ScinDigi;
  Double_t SciSignal[5];
  Double_t SciTime[5];
  Double_t SciChisq[5];
  Short_t  SciNDF[5];
  Short_t  SciID[5];

  tr->SetBranchAddress("CsiNumber", &CsinDigi);
  tr->SetBranchAddress("CsiSignal", CsiSignal);
  tr->SetBranchAddress("CsiTime"  , CsiTime);
  tr->SetBranchAddress("CsiHHTime", CsiHHTime);
  tr->SetBranchAddress("CsiChisq" , CsiChisq);
  tr->SetBranchAddress("CsiNDF"   , CsiNDF);
  tr->SetBranchAddress("CsiID"    , CsiID);

  tr->SetBranchAddress("EtcNumber", &ScinDigi);
  tr->SetBranchAddress("EtcSignal", SciSignal);
  tr->SetBranchAddress("EtcTime"  , SciTime);
  tr->SetBranchAddress("EtcChisq" , SciChisq);
  tr->SetBranchAddress("EtcNDF"   , SciNDF);
  tr->SetBranchAddress("EtcID"    , SciID);

  std::cout<< tr->GetEntries() << std::endl;
  
  TFile* tfout = new TFile("CrystalTime.root","recreate");

  TH2D* his2D = new TH2D("his2D","",2716,0,2716,280,-50,90);
  TH1D* his1D[40];
  for( int i = 0; i< 40; i++){    
    his1D[i] = new TH1D(Form("hisTest%d",i),"",280,-50,90); 
  }
  TH1D* hisLarge = new TH1D("hisLarge","",280,-50,90);
  TH1D* hisSmall = new TH1D("hisSmall","",280,-50,90);

  for( int ievent = 0; ievent < tr->GetEntries() ; ievent++){
    tr->GetEntry(ievent);

    for( int icsi = 0 ; icsi < CsinDigi; icsi++){
      if( CsiTime[icsi] > 50.  &&
	  CsiTime[icsi] < 325. &&
	  CsiSignal[icsi] > 500 &&
	  CsiChisq[icsi]/CsiNDF[icsi] < 40 &&
	  SciSignal[0] < 10000 &&
	  SciSignal[0]> 100      &&
	  SciTime[0] > 150     &&
	  SciTime[0] < 200  
	  
	  ){
	//std::cout << CsiID[icsi] << std::endl;
	his2D->Fill( CsiID[icsi],
		     SciTime[0] - (CsiHHTime[icsi] - TimeOffset[CsiID[icsi]]));
	bool Fill = false;
	if( CsiID[icsi] < 2240 ){
	  hisSmall->Fill(SciTime[0] - (CsiHHTime[icsi] - TimeOffset[CsiID[icsi]]));
	}else{
	  hisLarge->Fill(SciTime[0] - (CsiHHTime[icsi] - TimeOffset[CsiID[icsi]]));
	}	  
	for( int i  =0; i< 40; i++){
	  if( CsiID[icsi] == 1140 + i ){
	    his1D[i]->Fill(SciTime[0] - (CsiHHTime[icsi] - TimeOffset[CsiID[icsi]]));
	  }
	}
      }
    }
  }
    
  his2D->Draw("col");
  TProfile* pro = his2D->ProfileX();
  pro->Draw("same");
  his2D->Write();
  for( int i = 0; i< 40; i++){
    his1D[i]->Write();
  }
  hisLarge->Write();
  hisSmall->Write();
  tfout->Close();
  //app->Run();

}


  
