#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "csimap/CsiMap.h"

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>

int main( int argc, char** argv ){

  if( argc != 2 ) { 
    std::cout<< "Wrong argument" << std::endl;
    return -1; 
  }
  Int_t RunNumber = std::atoi( argv[1]);

  std::string IFDIR = std::getenv("ROOTFILEWAV");
  std::string OFDIR = std::getenv("ROOTFILELASER");
  char* InputFilename = Form("%s/run_wav_%d.root",IFDIR.c_str(),RunNumber);
  char* OutputFilename= Form("%s/run_wav_Laser_%d.root",OFDIR.c_str(),RunNumber);

  CsiMap* map = CsiMap::getCsiMap();

  TFile* tf = new TFile(InputFilename);
  TTree* tr = (TTree*)tf->Get("Tree");

  Int_t EventNo;

  Int_t CsiNumber;
  Double_t CsiSignal[2716];
  Double_t CsiTime[2716];
  Double_t CsiHHTime[2716];
  Double_t CsiChisq[2716];
  Short_t  CsiID[2716];
  Short_t  CsiNDF[2716];
  Short_t  CsiPosID[2716];

  Int_t LaserNumber;
  Double_t LaserSignal[10];
  Double_t LaserHHTime[10];
  Double_t LaserTime[10];
  Double_t LaserChisq[10];
  Short_t  LaserID[10];
  Short_t  LaserNDF[10];
  
  Int_t OEVNumber;
  Double_t OEVSignal[50];
  Double_t OEVTime[50];
  Double_t OEVChisq[50];
  Short_t  OEVID[50];
  Short_t  OEVNDF[50];

  Int_t CC03Number;
  Double_t CC03Signal[50];
  Double_t CC03Time[50];
  Double_t CC03Chisq[50];
  Short_t  CC03ID[50];
  Short_t  CC03NDF[50];

  tr->SetBranchAddress("EventNo",&EventNo);

  tr->SetBranchAddress("CsiNumber",&CsiNumber);  
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  tr->SetBranchAddress("CsiHHTime",CsiHHTime);//CsiNumber
  tr->SetBranchAddress("CsiChisq",CsiChisq);//CsiNumber
  tr->SetBranchAddress("CsiID",CsiID);//CsiNumber
  tr->SetBranchAddress("CsiNDF",CsiNDF);//CsiNumber

  tr->SetBranchAddress("CC03Number",&CC03Number);
  tr->SetBranchAddress("CC03Signal",CC03Signal);//CC03Number
  tr->SetBranchAddress("CC03Time",CC03Time);//CC03Number
  tr->SetBranchAddress("CC03Chisq",CC03Chisq);//CC03Number
  tr->SetBranchAddress("CsiID",CsiID);//CsiNumber
  tr->SetBranchAddress("CC03NDF",CC03NDF);//CC03Number
  tr->SetBranchAddress("CC03ID",CC03ID);//CC03Number

  tr->SetBranchAddress("OEVNumber",&OEVNumber);
  tr->SetBranchAddress("OEVSignal",OEVSignal);//OEVNumber
  tr->SetBranchAddress("OEVTime",OEVTime);//OEVNumber
  tr->SetBranchAddress("OEVChisq",OEVChisq);//OEVNumber
  tr->SetBranchAddress("OEVID",OEVID);//OEVNumber
  tr->SetBranchAddress("OEVNDF",OEVNDF);//OEVNumber
  tr->SetBranchAddress("OEVID",OEVID);//OEVNumber

  tr->SetBranchAddress("LaserNumber",&LaserNumber);
  tr->SetBranchAddress("LaserSignal",LaserSignal);//LaserNumber
  tr->SetBranchAddress("LaserHHTime",LaserHHTime);//LaserNumber
  tr->SetBranchAddress("LaserTime",LaserTime);//LaserNumber
  tr->SetBranchAddress("LaserChisq",LaserChisq);//LaserNumber
  tr->SetBranchAddress("LaserID",LaserID);//LaserNumber
  tr->SetBranchAddress("LaserNDF",LaserNDF);//LaserNumber

  TFile* tfOut = new TFile(OutputFilename,"recreate");
  TTree* trOut = new TTree("LaserTree","LaserOut");

  trOut->Branch("RunNo",&RunNumber,"RunNumber/I");
  trOut->Branch("EventNo",&EventNo,"EventNo/I");

  trOut->Branch("CsiNumber",&CsiNumber,"CsiNumber/I");
  trOut->Branch("CsiSignal",CsiSignal,"CsiSignal[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiTime",CsiTime,"CsiTime[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiHHTime",CsiHHTime,"CsiHHTime[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiChisq",CsiChisq,"CsiChisq[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiNDF",CsiNDF,"CsiNDF[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiID",CsiID,"CsiID[CsiNumber]/S");//CsiNumber
  trOut->Branch("CsiPosID",CsiPosID,"CsiPosID[CsiNumber]/S");//CsiNumber


  trOut->Branch("OEVNumber",&OEVNumber,"OEVNumber/I");
  trOut->Branch("OEVSignal",OEVSignal,"OEVSignal[OEVNumber]/D");//OEVNumber
  trOut->Branch("OEVTime",OEVTime,"OEVTime[OEVNumber]/D");//OEVNumber
  trOut->Branch("OEVChisq",OEVChisq,"OEVChisq[OEVNumber]/D");//OEVNumber
  trOut->Branch("OEVNDF",OEVNDF,"OEVNDF[OEVNumber]/S");//OEVNumber
  trOut->Branch("OEVID",OEVID,"OEVID[OEVNumber]/S");//OEVNumber

  trOut->Branch("CC03Number",&CC03Number,"CC03Number/I");
  trOut->Branch("CC03Signal",CC03Signal,"CC03Signal[CC03Number]/D");//CC03Number
  trOut->Branch("CC03Time",CC03Time,"CC03Time[CC03Number]/D");//CC03Number
  trOut->Branch("CC03Chisq",CC03Chisq,"CC03Chisq[CC03Number]/D");//CC03Number
  trOut->Branch("CC03NDF",CC03NDF,"CC03NDF[CC03Number]/S");//CC03Number
  trOut->Branch("CC03ID",CC03ID,"CC03ID[CC03Number]/S");//CC03Number

  trOut->Branch("LaserNumber",&LaserNumber,"LaserNumber/I");
  trOut->Branch("LaserSignal",LaserSignal,"LaserSignal[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserTime",LaserTime,"LaserTime[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserHHTime",LaserHHTime,"LaserHHTime[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserChisq",LaserChisq,"LaserChisq[LaserNumber]/D");//LaserNumber
  trOut->Branch("LaserID",LaserID,"LaserID[LaserNumber]/S");//LaserNumber
  trOut->Branch("LaserNDF",LaserNDF,"LaserNDF[LaserNumber]/S");//LaserNumber  

  TH1D*  hisLaserSignal = new TH1D("hisLaserSignal","hisLaserSignal",100,0,5000);
  TH2D*  hisL1 = new TH2D("hisLaserSignalCsiNumber",";CsiNumber;LaserSignal[cnt]",2716,0,2716,500,0,5000);
  TH2D*  hisL2 = new TH2D("hisLaserSignalLaserTime",";LaserTime[ns];LaserSignal[cnt]",500,0,500,500,0,5000);
  for( int i = 0; i< tr->GetEntries(); i++){
    tr->GetEntry(i);
    for( int i = 0; i< 2716; i++){
      CsiPosID[i] = -1;
    }
    for( int i =0 ; i< CsiNumber; i++){
      Double_t x(0),y(0);
      x = map->getX(CsiID[i]);
      y = map->getY(CsiID[i]);
      if( x > 0 ){
	if( y > 0 ){ CsiPosID[i]  = 0; }
	else{ CsiPosID[i] = 1; }
      }else{
	if( y > 0 ){ CsiPosID[i] = 2; }
	else{ CsiPosID[i]  =3; }
      }
    }

    if( LaserNumber != 0){
      hisLaserSignal->Fill(LaserSignal[0]);
      hisL1->Fill(CsiNumber,LaserSignal[0]);
      hisL2->Fill(LaserTime[0],LaserSignal[0]);
    }
    if( LaserSignal[0]> 800 ){
      trOut->Fill();
    }
  }

  hisLaserSignal->Write();
  hisL1->Write();
  hisL2->Write();
  tfOut->Close();
}
