#include "cluster/ClusterFinder.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "Structs.h"

#include "E14ConvReader.h"
#include "E14IDHandler.h"
#include "E14ConvWriter.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TChain.h"
#include "TROOT.h"
#include "Math/WrappedFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMethods.h"
#include "Math/WrappedTF1.h"
#include "TProfile.h"

int
main( int argc, char** argv ){
  
  std::ifstream ifs("/home/jwlee/local/Analysis/RunList/EtRunList.txt");
  
  TChain* ch  = new TChain("Waveform");
  std::string Env = std::getenv("ROOTFILE_WAV");
  Int_t RunNumber; 
  Int_t ModuleNumber;
  Double_t Waveform[48];
  Short_t  TimeInfo[48];
  Double_t ChisqNDF;
  Double_t PeakTime;
  Double_t HHTime;
  Double_t Height;
  Double_t Pedestal;
  
  ch->SetBranchAddress("ModuleNumber", (void*)&ModuleNumber);
  ch->SetBranchAddress("Waveform",(void*)Waveform);
  ch->SetBranchAddress("TimeInfo",(void*)TimeInfo);
  ch->SetBranchAddress("PeakTime",(void*)(&PeakTime));
  ch->SetBranchAddress("HHTime",(void*)(&HHTime));
  ch->SetBranchAddress("Height",(void*)(&Height));
  ch->SetBranchAddress("ChisqNDF",(void*)(&ChisqNDF));
  while( ifs >> RunNumber ){
    ch->Add(Form("%s/TEMPLATE_FIT_RESULT_DUMP_%d.root",
		 Env.c_str(), RunNumber)); 
  }
      
  TFile* tfout = new TFile("TEMPLATE_HEIGHT_500_1000.root","recreate");
  TH2D* hisTemp_CsI[2716];
  TProfile* prof[2716];
  TH1D* hisInsertedChannel = new TH1D("hisInsertedChannel","hisInsertedChannel",
				      2716,0,2716);
  for( int i = 0; i< 2716; i++){
    hisTemp_CsI[i] = new TH2D(Form("hisTempCsI_Height_500_1000_%d", i),
			      Form("hisTempCsI_Height_500_1000_%d", i),
			      400,-150,250,150,-0.25,1.25);
  }
  
  TGraph* gr = new TGraph();  
  for( int ievent  = 0; ievent < ch->GetEntries(); ievent++){
    ch->GetEntry(ievent);
    //std::cout <<ChisqNDF << std::endl;
    if( ChisqNDF >30 ) {continue; }    
    if(  Height < 500 || Height > 1000 ){ continue; }
    //gr->Set(0);
    bool frontSignal = false;
    bool backSignal  = false;

    for( int ipoint  = 0; ipoint < 48; ipoint++){
      Double_t Timing = TimeInfo[ipoint] - PeakTime;
      Double_t Signal = (Waveform[ipoint] - Pedestal)/Height;
      
      //gr->Set( ipoint , Timing, Signal);
      if( Timing < -80 ){ 
	if( Signal > 0.1 ){
	  frontSignal  = true;
	}
      }else if( Timing > 100 ){
	if( Signal > 0.2 ){
	  backSignal = true;
	}
      }      
    }

    if( !frontSignal && !backSignal ){
      for( int ipoint = 0; ipoint < 48; ipoint++ ){
	Double_t Timing = TimeInfo[ipoint] - PeakTime;
	Double_t Signal = (Waveform[ipoint] - Pedestal)/Height;
	hisTemp_CsI[ModuleNumber]->Fill( Timing, Signal );
      }
      hisInsertedChannel->Fill(ModuleNumber);
    }
  }
  
  for( Int_t ich  =0; ich < 2716; ich++){
    prof[ich] = hisTemp_CsI[ich]->ProfileX();
    prof[ich]->Write();
    hisTemp_CsI[ich]->Write();
  }
  hisInsertedChannel->Write();
  tfout->Close();

}  

  
  
  
