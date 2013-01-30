#include "Waveform_Reader.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include <iostream>
#include "WaveformFitter.h"

int main( int argc , char** argv ){
  
  std::string InputFilename = "run_conv_mode5.root";
  std::string OutputFilename= "Waveform_Height.root";
  TApplication* app = new TApplication("app",&argc, argv);
  
  TFile* tf = new TFile(InputFilename.c_str());
  TTree* tr = (TTree*)tf->Get("EventTree");
  TFile* tfout= new TFile( OutputFilename.c_str(), "RECREATE");
  
  WaveformFitter* wavFitter =new WaveformFitter( 48, kFALSE);
  Waveform_Reader* r = new Waveform_Reader(tr);
  
  Int_t nEntries = tr->GetEntries();
  TGraph* gr = new TGraph();
  gr->SetMarkerStyle(4);
  TCanvas* can = new TCanvas("can","",800,800);
  TH2D* hisTemplate[120];
  for( int i = 0; i< 120; i++){
    hisTemplate[i] = new TH2D(Form("hisTemplate_%d",i),
			      Form("hisTemplate_%d;Time[ns];Normalized Pulse",i*5),500,-250,300,150,-0.25,1.25);
  }
  for( int ievent = 0; ievent < nEntries; ievent++){
    r->GetEntry(ievent);
    gr->Set(0);
    int Index = (int)(TMath::Abs(r->volt)/5); 
    for( int iPoint = 0; iPoint < 48; iPoint++){
      hisTemplate[Index]->Fill( r->xaxis[iPoint]*8-r->param[1],
				(r->data[iPoint]-r->param[4])/r->param[0]);
    }
    


    /*
    gr->Draw("AP");
    can->Update();
    can->Modified();
    getchar();
    */
    /*
    bool fit = wavFitter->Fit(gr);
    if( fit ){
      TF1* func = wavFitter->GetFunction();
      wavFitter->Clear();
    }
    */
  }

  for( int i = 0; i< 120; i++){
    hisTemplate[i]->Write();
  }
  tfout->Close();
  //app->Run();
}
