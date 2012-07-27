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
#include "WaveformFitter.h"

#include "TApplication.h"
#include "TPostScript.h"

int
main(int argc, char** argv ){
  TApplication* app = new TApplication("app",&argc,argv);
  std::string TestFileDir="/Volume0/ExpData/2012_Feb_Beam/RootFile_wav/run3949_wavDump_Cosmic.root";
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
  TGraph* wf[20];
  for( int i = 0; i< 20; i++){
    wf[i] = new TGraph();
  }
  int igr = 0;

  double sum;
  double sqsum;
  double rms;
  double min;
  double max;

  std::cout<< tr->GetEntries() << std::endl;  
  long nentries = tr->GetEntries();   

  TH1D* hisRMS[20];
  TH1D* hisSpc[20];
  for( int i = 0; i< 20; i++){
    hisRMS[i] = new TH1D(Form("hisRMS%d",i),"",200,0,10000);
    hisSpc[i] = new TH1D(Form("hisSpc%d",i),"",200,0,100000);
  }

  for( int ievent = 0; ievent < nentries; ievent++){
    tr->GetEntry( ievent );
    igr = ievent%20;
    wf[igr]->Set(0);
    wf[igr]->SetNameTitle(Form("Gr%d",igr),Form("%d:CH%d",EventNumber, ID));

    sum   = 0;
    sqsum = 0;
    rms   = 0;
    min   = 16000;
    max   = 0; 
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      sum += Data[ipoint];
      sqsum += Data[ipoint]*Data[ipoint];
      if( Data[ipoint] < min ){ min = Data[ipoint]; }
      if( Data[ipoint] > max ){ max = Data[ipoint]; }

      wf[igr]->SetPoint( ipoint, Time[ipoint], Data[ipoint] );
      
    }
    rms = TMath::Sqrt( sqsum - sum*sum );
    hisRMS[ID]->Fill( rms );
    hisSpc[ID]->Fill( sum - min*48 );
    bool Fit = wavFitter->Fit(wf[igr]);
  }

  TPostScript* ps = new TPostScript("Cosmic.ps",111);
  TCanvas* can = new TCanvas( "can", "", 1250,1000);
  can->Divide(5,4);
  ps->NewPage();
  for( int i = 0; i< 20; i++){
    can->cd(i+1);
    hisRMS[i]->Draw();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
  }
  can->Update();
  can->Modified();
  ps->NewPage();
  for( int i = 0; i< 20; i++){
    can->cd(i+1);
    hisSpc[i]->Draw();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
  }
  can->Update();
  can->Modified();
   
  ps->Close();

  app->Run();
}

