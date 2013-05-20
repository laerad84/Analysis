#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

int main( int argc, char** argv){
  
  std::string ANALYSISLIB = std::getenv("ANALYSISLIB");
  TChain* runTree = new TChain("runTree");
  runTree->Add(Form("%s/Data/RunDataAll.root",ANALYSISLIB.c_str()));  
  /*
  TFile* tf = new TFile(Form("%s/Data/RunDataAll.root",ANALYSISLIB.c_str()));  
  TTree* runTree = (TTree*)tf->Get("runTree");
  */
  Int_t runID;
  UInt_t StartDate;
  UInt_t EndDate;
  Int_t  bufSize;
  Int_t TrigType;
  Int_t nTrigMAX;
  Int_t trigdelay;
  Int_t Threshold;
  Int_t Threshold_offset;
  Int_t SpillEtRecordCycle;
  Int_t SECTMonRecordCycle;
  UInt_t TotalnSpill;
  UInt_t TotalnSEC;
  UInt_t TotalnTMon;
  UInt_t TotalnTrigRequested;
  UInt_t TotalnTrigAccepted;
  UInt_t TotalnEt;
  UInt_t TotalnEt_withVeto;
  UInt_t TotalnTrigOR;
  UInt_t TotalnCalibTrig;
  UInt_t TotalnChamberTrig;
  UInt_t TotalnVetoCount;

  TBranch        *b_runID;   //!
  TBranch        *b_StartDate;   //!
  TBranch        *b_EndDate;   //!
  TBranch        *b_bufSize;   //!
  TBranch        *b_TrigType;   //!
  TBranch        *b_nTrigMAX;   //!
  TBranch        *b_trigdelay;   //!
  TBranch        *b_Threshold;   //!
  TBranch        *b_Threshold_offset;   //!
  TBranch        *b_SpillEtRecordCycle;   //!
  TBranch        *b_SECTMonRecordCycle;   //!
  TBranch        *b_TotalnSpill;   //!
  TBranch        *b_TotalnSEC;   //!
  TBranch        *b_TotalnTMon;   //!
  TBranch        *b_TotalnTrigRequested;   //!
  TBranch        *b_TotalnTrigAccepted;   //!
  TBranch        *b_TotalnEt;   //!
  TBranch        *b_TotalnEt_withVeto;   //!
  TBranch        *b_TotalnTrigOR;   //!
  TBranch        *b_TotalnCalibTrig;   //!
  TBranch        *b_TotalnChamberTrig;   //!
  TBranch        *b_TotalnVetoCount;   //!
  
  runTree->SetBranchAddress("runID", &runID, &b_runID);
  runTree->SetBranchAddress("StartDate", &StartDate, &b_StartDate);
  runTree->SetBranchAddress("EndDate", &EndDate, &b_EndDate);
  runTree->SetBranchAddress("bufSize", &bufSize, &b_bufSize);
  runTree->SetBranchAddress("TrigType", &TrigType, &b_TrigType);
  runTree->SetBranchAddress("nTrigMAX", &nTrigMAX, &b_nTrigMAX);
  runTree->SetBranchAddress("trigdelay", &trigdelay, &b_trigdelay);
  runTree->SetBranchAddress("Threshold", &Threshold, &b_Threshold);
  runTree->SetBranchAddress("Threshold_offset", &Threshold_offset, &b_Threshold_offset);
  runTree->SetBranchAddress("SpillEtRecordCycle", &SpillEtRecordCycle, &b_SpillEtRecordCycle);
  runTree->SetBranchAddress("SECTMonRecordCycle", &SECTMonRecordCycle, &b_SECTMonRecordCycle);
  runTree->SetBranchAddress("TotalnSpill", &TotalnSpill, &b_TotalnSpill);
  runTree->SetBranchAddress("TotalnSEC", &TotalnSEC, &b_TotalnSEC);
  runTree->SetBranchAddress("TotalnTMon", &TotalnTMon, &b_TotalnTMon);
  runTree->SetBranchAddress("TotalnTrigRequested", &TotalnTrigRequested, &b_TotalnTrigRequested);
  runTree->SetBranchAddress("TotalnTrigAccepted", &TotalnTrigAccepted, &b_TotalnTrigAccepted);
  runTree->SetBranchAddress("TotalnEt", &TotalnEt, &b_TotalnEt);
  runTree->SetBranchAddress("TotalnEt_withVeto", &TotalnEt_withVeto, &b_TotalnEt_withVeto);
  runTree->SetBranchAddress("TotalnTrigOR", &TotalnTrigOR, &b_TotalnTrigOR);
  runTree->SetBranchAddress("TotalnCalibTrig", &TotalnCalibTrig, &b_TotalnCalibTrig);
  runTree->SetBranchAddress("TotalnChamberTrig", &TotalnChamberTrig, &b_TotalnChamberTrig);
  runTree->SetBranchAddress("TotalnVetoCount", &TotalnVetoCount, &b_TotalnVetoCount);


	
  std::vector<int> RunList;
  
  std::ifstream ifs("/home/jwlee/local/Analysis/RunList/KLRunList_2.txt");
  int tmpRunNo;
  while( ifs >> tmpRunNo ){
    RunList.push_back(tmpRunNo);
    std::cout<< tmpRunNo << std::endl;
  }

  Int_t SumSEC = 0;
  Int_t SumTMON = 0; 
  Int_t SumTrigReq = 0;
  Int_t SumTrigAcc = 0; 
  Int_t SumTrigOR  = 0;

  std::vector<int>::iterator it;

  for( int i = 0; i< runTree->GetEntries(); i++){
    runTree->GetEntry(i);
    //std::cout<< runID << std::endl;
    //std::cout<< TotalnSEC <<"\t" << TotalnTMon << "\t" 
    //<< TotalnTrigRequested << "\t" << TotalnTrigAccepted << std::endl; 
    if( runID <  4249 ){ continue;}
    if( runID >= 4625 ){ break;}
    for( it = RunList.begin(); it != RunList.end(); it++){
      //std::cout << runID << "\t" << (*it) << std::endl;
      if( runID == (*it) ){
	SumSEC += TotalnSEC;
	SumTMON += TotalnTMon;
	SumTrigReq += TotalnTrigRequested;
	SumTrigAcc += TotalnTrigAccepted;
	SumTrigOR  += TotalnTrigOR;
	break;
      }   
    } 
  }
  std::cout << "SEC TOTAL       : " << SumSEC << std::endl;
  std::cout << "TMon TOTAL      : " << SumTMON << std::endl;
  std::cout << "TrigReq TOTAL   : " << SumTrigReq << std::endl;
  std::cout << "TrigAcc TOTAL   : " << SumTrigAcc << std::endl;
  std::cout << "Avg Accept ratio: " << (double)SumTrigAcc/(double)SumTrigOR  << std::endl;
  std::cout << "Avg Accept ratio: " << (double)SumTrigAcc/(double)SumTrigReq << std::endl;
}

