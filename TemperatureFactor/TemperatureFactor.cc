#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"

int 
main( int argc ,char** argv){

  std::string DataFileDirStr = "/nfs/g/had/koto/ps/klea/work/jwlee/local/Analysis/Data/2012_FEB/RunInfo_TEMP_ROOTFILE";
  std::string DataNameStr    = "/RunInfo_Temp_%04d.root";
  std::string HisNameStr     = "hisTemp_CH%02d";
  std::string FileName;
  FileName += DataFileDirStr.c_str();
  FileName += DataNameStr.c_str();
  int ChannelID[4] = {44,53,46,55};

  TFile* tfout = new TFile("TemperatureCorrectionFactor.root","recreate");
  TTree* trTemp= new TTree("TemperatureCorrectionCsI",
			   "Temperature Correction Factor Tree for CsI Calorimeter");  
  int    RunNumber;
  double Temperature;
  double CorrectionFactor;
  int    CorrectionFlag;
  trTemp->Branch("CorrectionFlag",&CorrectionFlag,"CorrectionFlag/I");
  trTemp->Branch("RunNumber",&RunNumber,"RunNumber/I");
  trTemp->Branch("Temperature",&Temperature,"Temperature/D");
  trTemp->Branch("CorrectionFactor",&CorrectionFactor,"CorrectionFactor/D");
  
  TFile* tf = NULL;
  TH1D*  hisTemp[4];
  for( int i = 0; i< 4; ++i){
    hisTemp[i] = NULL;
  }

  double TempBCK=0;
  double CorrBCK=0;

  for( int iRun = 0; iRun < 5000; ++iRun ){
    tf = TFile::Open(Form( FileName.c_str(),iRun));
    RunNumber = iRun;
    CorrectionFlag = 0; 
    if( tf == NULL ){      
      Temperature = -273;
      CorrectionFactor = 1;
      trTemp->Fill();
      CorrectionFlag = -1;
      continue; 
    }

    if( iRun %100 == 0 ){ 
      std::cout << iRun << std::endl;
    }

    Temperature = 0;     
    int tempDataNumber = 0; 
    for( int iCH = 0; iCH<4 ; ++iCH){
      hisTemp[iCH] = (TH1D*)tf->Get(Form(HisNameStr.c_str(), ChannelID[iCH]));
      double localTemp =  hisTemp[iCH]->GetMean();      
      Temperature +=localTemp;
      if(localTemp > 0 && localTemp < 100 && hisTemp[iCH]->Integral() > 0){
	++tempDataNumber;
      }
      delete hisTemp[iCH];
      hisTemp[iCH] = NULL;      
    }
    Temperature = Temperature/4; 
    CorrectionFactor = 1.38286 -0.01418*Temperature;

    if(tempDataNumber != 4){
      Temperature = TempBCK;
      CorrectionFactor = CorrBCK;
      trTemp->Fill();
      tf->Close();
      tf = NULL;
      continue;
    }
    
    CorrectionFlag = 1;
    TempBCK  = Temperature;
    CorrBCK  = CorrectionFactor;
    trTemp->Fill();
    tf->Close();
    tf = NULL;
  }
  tfout->cd();
  trTemp->Write();
  tfout->Close();
  
}



    
    

    

    



  
  
  

  


