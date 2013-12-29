
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <iostream>


int main( int argc, char** argv){
//void DrawTimingStability(){

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
