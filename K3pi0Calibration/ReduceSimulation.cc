#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cstdio>
int main( int argc, char** argv ){
  //int Number = std::atoi( argv[1] );
  //std::cout<< "RunNumber :" << Number << std::endl;  
  //std::string ROOTFILE_GSIM = "";  
  //TFile* tfReduced = new TFile(Form("%s/Out_KL3pi0Fast_Reduced_1E7_%d.root",ROOTFILE_GSIM.c_str(),Number),"recreate");
  std::ofstream ofs(argv[3]);
  
  TFile* tfReduced = new TFile(argv[2],"recreate");
  TTree* trC;
  TFile* tf = new TFile(argv[1]);
  TTree* tr = (TTree*)tf->Get("eventTree00");
  tfReduced->cd();
  trC = tr->CopyTree("GenParticle.briefTracks.end_v.fZ[2]==6148&&GenParticle.briefTracks.end_v.fZ[3]==6148&&GenParticle.briefTracks.end_v.fZ[5]==6148&&GenParticle.briefTracks.end_v.fZ[6]==6148&&GenParticle.briefTracks.end_v.fZ[8]==6148&&GenParticle.briefTracks.end_v.fZ[9]==6148");
  Int_t nEntries = trC->GetEntries();
  ofs << nEntries << std::endl;
  ofs.close();
  trC->Write();  
  tfReduced->Close();
}
