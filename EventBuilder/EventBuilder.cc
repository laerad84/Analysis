////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// EventBuilder - convert convfile to wavfile 
// EventBuilder [RunNumber]
// Please Edit include/E14EventBuilder_V0.h && src/E14EventBuilder_V0.cc for fitting your system 
// before Convert data.
// Environment 
// ROOTFILE_CONV - ConvFile Directory
// ROOTFILE_WAV  - WAVFile Directory
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include "E14EventBuilder_V0.h"
#include "TApplication.h"
#include "TPostScript.h"
#include "sys/stat.h"
#include <string>
int
main( int argc, char** argv ){
  
  int RunNumber;
  RunNumber = atoi( argv[1] );
  
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");

  TFile* tfOut = new TFile(Form("%s/run_wav_%d.root",ROOTFILE_WAV.c_str(),RunNumber),"RECREATE");
  TTree* trOut = new TTree("Tree","");
  TCanvas* can  = new TCanvas("can","",800,1200);
  E14EventBuilder_V0* Converter = new E14EventBuilder_V0(trOut,RunNumber);
  tfOut->cd();
  Converter->LoopAll();
  trOut->Write();
  tfOut->Close();

}
