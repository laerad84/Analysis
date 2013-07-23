#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector> 

#include "RunTree.h"
#include "TFile.h"
#include "TTree.h"
int main( int argc, char** argv){
  
  std::string RUNListFile=argv[1];
  std::string RUNTRIGFILE="/Volume0/ExpData/2012_Feb_Beam/TrigCrateData/RunDataAll.root";
  int tmpRunNumber;
  
  std::cout<< "RunNumber List File name" << std::endl;
  std::cout<< RUNListFile << std::endl;
  std::ifstream ifs( RUNListFile.c_str());
  TFile* tf = new TFile(Form(RUNTRIGFILE.c_str()));
  TTree* tr = (TTree*)tf->Get("runTree");
  RunTree* runTree =new RunTree(tr);
  std::vector<int> runNumberVec;
  while( ifs >> tmpRunNumber ){
    //std::cout<< tmpRunNumber  << std::endl;
    runNumberVec.push_back( tmpRunNumber ); 
  }
  std::vector<int>::iterator it = runNumberVec.begin();
  Int_t  secNumber= 0; 
  Int_t  tmonNumber=0;
  Int_t TotalnTrigRequested=0;
  Int_t TotalnTrigAccepted=0;
  Int_t TotalnCalibTrig=0;
  
  std::cout<<"LOOP"<< std::endl;
  std::cout<< tr->GetEntries() << std::endl;
  for( int ievent  = 0; ievent < tr->GetEntries(); ievent++){
    tr->GetEntry(ievent);
    if( (*it) > runTree->runID){
      continue;
    }else if( (*it) < runTree->runID ){
      while( (*it) < runTree->runID ){
	std::cout<< "RunNumber Missed:" << (*it) << std::endl;
	it++;
      }
    }
    if( (*it) == runTree->runID){
      secNumber += runTree->TotalnSEC;
      TotalnTrigAccepted  += runTree->TotalnTrigAccepted;
      TotalnTrigRequested += runTree->TotalnTrigRequested;
      TotalnCalibTrig     += runTree->TotalnCalibTrig;
      std::cout<< (*it) << "\t" << runTree->TotalnSEC << "\t" 
	       << runTree->TotalnTrigAccepted << "\t"
	       << runTree->TotalnTrigRequested << "\t"
	       << (double)(runTree->TotalnTrigAccepted)/(double)(runTree->TotalnTrigRequested) << std::endl;
      //std::cout<< (*it) << "\t" << secNumber << std::endl;
      it++;
    }

    if( it == runNumberVec.end() ) {break; }
  }  
  std::cout<< "TotalnSEC          : " << secNumber           << std::endl;
  std::cout<< "TotalnTrigAccepted : " << TotalnTrigAccepted  << std::endl;
  std::cout<< "TotalnTrigRequested: " << TotalnTrigRequested << std::endl;
  std::cout<< "Avg. Accept Ratio  : " << (double)TotalnTrigAccepted/(double)TotalnTrigRequested << std::endl;
}
