#include <iostream>
#include <iomanip>
#include <cstring>
#include <string>

#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TROOT.h"

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "Structs.h"
#include "Environment.h"
#include "WaveformFitter.h"

#include "TApplication.h"


int
main(int argc, char** argv ){
  TApplication* app = new TApplication("app",&argc,argv);
  std::string TestFileDir="/Volume0/ExpData/2012_Feb_Beam/wave_data/run4502_wavDump_CV.root";  
  TFile* tf = new TFile(TestFileDir.c_str());  
  TTree* tr = (TTree*)tf->Get("WFTree");
  
  Int_t Data[48];
  Int_t Time[48];
  Int_t ID;
  Int_t EventNumber;
  
  tr->SetBranchAddress("Data",Data);
  tr->SetBranchAddress("Time",Time);
  tr->SetBranchAddress("ID"  ,&ID);
  tr->SetBranchAddress("EventNumber",&EventNumber);

  WaveformFitter* wavFitter = new WaveformFitter( 48, kFALSE );
  TGraph* wf[16];
  for( int i = 0; i< 16; i++){
    wf[i] = new TGraph();
  }
  int igr = 0;
  TCanvas* can = new TCanvas( "can","", 1200,1000);
  can->Divide(4,4);
  std::cout<< tr->GetEntries() << std::endl;  
  long nentries = tr->GetEntries();   
  for( int ievent = 0; ievent < nentries; ievent++){
    tr->GetEntry( ievent );
    igr = ievent%16;
    wf[igr]->Set(0);
    wf[igr]->SetNameTitle(Form("Gr%d",igr),Form("%d:CH%d",EventNumber, ID));
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      wf[igr]->SetPoint( ipoint, Time[ipoint], Data[ipoint] );
    }
    bool Fit = wavFitter->Fit(wf[igr]);
    can->cd( igr +1 );
    wf[igr]->Draw("AP");
    wavFitter->Clear();    
    if(( ievent+1 )%16 == 0 ){
      can->Update();
      can->Modified();
      getchar();
    }    
  }
  app->Run();
}

