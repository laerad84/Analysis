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
  int Number = std::atoi( argv[1] );
  std::cout<< "RunNumber :" << Number << std::endl;
  
  std::string ROOTFILE_GSIM = "";
  
  TFile* tfReduced = new TFile(Form("%s/Out_KL3pi0Fast_Reduced_1E7_%d.root",ROOTFILE_GSIM.c_str(),Number),"recreate");
  TTree* trC;
  TFile* tf = new TFile(Form("%s/Out_KL3pi0Fast.mac_10000000_%d.root",ROOTFILE_GSIM.c_str(),Number));
  TTree* tr = (TTree*)tf->Get("eventTree00");
  tfReduced->cd();
  trC = tr->CopyTree("GenParticle.briefTracks.end_v.fZ[2]==6148&&GenParticle.briefTracks.end_v.fZ[3]==6148&&GenParticle.briefTracks.end_v.fZ[5]==6148&&GenParticle.briefTracks.end_v.fZ[6]==6148&&GenParticle.briefTracks.end_v.fZ[8]==6148&&GenParticle.briefTracks.end_v.fZ[9]==6148");
  trC->Write();  
  tfReduced->Close();
}
