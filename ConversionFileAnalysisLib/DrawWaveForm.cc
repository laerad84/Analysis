#include <iostream>
#include <fstream>
#include <cstdio>
#include <iomanip>

#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPDF.h"
#include "TLine.h"
int
main( int argc , char** argv ){
  TFile* tf = new TFile("TEMPLETE_OUT_COSMIC.root");
  TGraphErrors* waveform[2716];
  for( int i = 0; i< 2716; i++){
    waveform[i] = (TGraphErrors*)tf->Get(Form("Waveform_Cosmic_%d",i));
  }
  TPDF* pdf = new TPDF("waveform_Cosmic.pdf",111);
  TCanvas* can = new TCanvas("can","",0,0,800,800);
  gPad->SetGridx();
  gPad->SetGridy();
  TLine* line = new TLine(-200,0,200,0);
  for( int i = 0; i< 2716;i++){
    waveform[i]->Draw("AP");
    if( waveform[i]->GetN() ==  0 ){ 
      if(i< 2240 || 
	 (i> 2240+119 && i < 2716 -119 )){
	std::cout<< "Entries0 :" << i << std::endl;
      }
    }
    can->Update();
    can->Modified();
  }
  pdf->Close();
}
