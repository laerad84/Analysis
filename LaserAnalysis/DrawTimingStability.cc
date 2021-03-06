
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <iostream>
#include <fstream>
#include "TMath.h"

Double_t LaserDelayFunc(double* x, double* p ){
  double value = p[0] + p[1]*TMath::Log( 1 + p[2]*TMath::Exp(x[0]/2000));
  return value;
}


int main( int argc, char** argv){
//void DrawTimingStability(){

  Double_t par0[3] = {0,0,0};
  Double_t par1[3] = {1.933,2.851,1.234};
  Double_t par2[3] = {0.01316,0.01282,0.04478};

  std::ifstream ifs("/home/had/jwlee/local/Analysis/AnalysisLib/Data/ch_map_CsI_L1.txt");
  if( !ifs.is_open()){std::cout<< "Error" << std::endl;return -1; }
  Int_t Crate[2716]={-1};
  Int_t FADC[2716]={-1};
  Int_t Ch[2716]={-1};
  Int_t L1[2716]={-1};
  Int_t tmpID,tmpC,tmpF,tmpH,tmpL;
  while( ifs >> tmpID >> tmpC >> tmpF >> tmpH >> tmpL ){
    std::cout<< tmpID << "\t" << tmpC << "\t" << tmpF << "\t" << tmpH << "\t" << tmpL << std::endl;
    Crate[tmpID] = tmpC;
    FADC[tmpID]  = tmpF;
    Ch[tmpID] = tmpH;
    L1[tmpID] = tmpL;
  }

  TFile* tf  =new TFile("LaserTiming.root","recreate");
  TTree *tr  = new TTree("LaserTimingTree","");
  Int_t RunNumber;
  Int_t ID[2716];
  Double_t Timing[2716];
  Double_t TimingRMS[2716];
  Double_t Height[2716];
  Double_t HeightRMS[2716];
  Double_t Entries[2716];

  tr->Branch("RunNumber",&RunNumber,"RunNumber/I");
  tr->Branch("ID",ID,"ID[2716]/I");
  tr->Branch("Timing",Timing,"Timing[2716]/D");
  tr->Branch("TimingRMS",TimingRMS,"TimingRMS[2716]/D");
  tr->Branch("Height",Height,"Height[2716]/D");
  tr->Branch("HeightRMS",HeightRMS,"HeightRMS[2716]/D");
  tr->Branch("Entries",Entries,"Entries[2716]/D");

  for( int i = 4158; i < 4738; i++){
  //for( int i = 4158; i < 4200; i++){
    if( i== 4225 ){ continue; }
    if( i== 4354 ){ continue;}
    TFile* tfin  = new TFile(Form("Data/LaserTimeStability_%d.root",i));
    if( !tfin->IsOpen() ){ continue; }
    if( tfin->GetSize() < 1000 ){ continue; }
    if( tfin == NULL ){ continue; }

    TH1D*  hisTimeDelta[2716];
    TH1D*  hisHeight[2716];
    RunNumber = i;
    std::cout<< RunNumber << std::endl;
    for( int j  =0; j< 2716; j++){
      
      hisTimeDelta[j] = (TH1D*)tfin->Get(Form("hisTimeDelta_%d_%d",j,i));
      hisHeight[j] = (TH1D*)tfin->Get(Form("hisHeight_%d_%d",j,i));
      Timing[j] = 0;
      Height[j] = 0;
      TimingRMS[j] = 0;
      HeightRMS[j] = 0;
      Entries[j] = 0;      
      ID[j] = j;
      if( hisTimeDelta[j] == NULL || hisTimeDelta[j]->GetEntries() < 10 ){ continue; }
      Timing[j]= hisTimeDelta[j]->GetMean();
      Height[j]= hisHeight[j]->GetMean();
      TimingRMS[j] = hisTimeDelta[j]->GetRMS();
      HeightRMS[j] = hisHeight[j]->GetRMS();
      Entries[j]   = hisHeight[j]->GetEntries();
      delete hisTimeDelta[j];
      delete hisHeight[j];
    }
    tfin->Close();
    tr->Fill();
    
  }
  tf->cd();
  tr->Write();
  tf->Close();


}
