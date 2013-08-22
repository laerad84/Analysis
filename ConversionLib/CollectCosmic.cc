#include "TFile.h"
#include "TTree.h"

#include <iostream>
#include <fstream>

int main(int argc, char** argv){

  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  Int_t RunNumber = atoi(argv[1]);
  TFile* tfin = new TFile(Form("%s/run_wav_%d.root", ROOTFILE_WAV.c_str(),RunNumber));

  TTree* tr = (TTree*)tfin->Get("Tree");
  TFile* tfout = new TFile(Form("%s/CosmicRootFile/run_wav_%d.root",ROOTFILE_WAV.c_str(),RunNumber),"RECREATE");
  TTree* ntr = tr->CopyTree("CosmicTrig==1");
  ntr->Write();
  tfout->Close();
}
