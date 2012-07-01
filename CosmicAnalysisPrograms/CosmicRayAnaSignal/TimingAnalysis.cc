#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>

#include "TROOT.h"
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TClonesArray.h"


#include "IDHandler.h"
#include "E14CosmicOut.h"

int
main( int argc , char** argv){
  const int nChannel = 2716;
  IDHandler* handler = new IDHandler();

  TFile* tf = new TFile("test_All_3897.root");
  TTree* tr = (TTree*)tf->Get("CosmicOut");

  E14CosmicOut* CosmicOut = new E14CosmicOut( tr ); 
  std::cout << CosmicOut->fChain << std::endl;
  std::cout << CosmicOut->fChain->GetEntries() << std::endl;
  std::cout<< __LINE__ << std::endl;
  
  TFile* tfout = new TFile("TestTimeOutLogic.root","recreate");
  std::cout<< __LINE__ << std::endl;  
  std::cout<< "Entries::" << CosmicOut->GetEntries() << std::endl;
  TH1D* hisTimeDiff[nChannel];
  for( int i = 0; i< nChannel; i++){
    hisTimeDiff[i] =  new TH1D(Form("his__%d", i),Form("__%d;ns;",i), 50,-25,25);
  }
  double x[2]; y[2];
  for( int iID = 0; iID < nChannel ; iID++){
    handler->GetMetricPosition( iID, x[0], y[0] );

    for( int i = 0; i< nChannel; i++){
      hisTimeDiff[i]->Reset();
      hisTimeDiff[i]->SetNameTitle(Form("his_%d_%d", iID, i),Form("%d_%d;ns;",iID,i));
      //hisTimeDiff[i]  = new TH1D(Form("his_%d_%d", iID, i),Form("%d_%d;ns;",iID,i),50,-25,25);
    }
    
    for( int ievent = 0; ievent < CosmicOut->GetEntries(); ievent++){
      if( (ievent % 1000) == 0  && ievent ){ std::cout << "EVENT:" << ievent << std::endl; }        
      CosmicOut->GetEntry(ievent);
      // Except Laser Event // 
      if( CosmicOut->CosmicFit != 1  && CosmicOut->nDigi!=0 && CosmicOut->nDigi < 300 ){ continue; }
      if( CosmicOut->CsIHHTiming[iID] < -500 ){ continue; }
      for( int jID = 0; jID < nChannel; jID++){
	if( CosmicOut->CsIHHTiming[jID] < -500 ){ continue; }
	handler->GetMetricPosition( jID, x[1], y[1] );
	hisTimeDiff[jID]->Fill(CosmicOut->CsIHHTiming[iID] - CosmicOut->CsIHHTiming[jID] );
      }
    }

    for( int i = 0; i< nChannel; i++){
      hisTimeDiff[i]->Write();
    }    
  }

  for( int i = 0; i< nChannel; i++){
    delete hisTimeDiff[i];
    hisTimeDiff[i] = NULL;
  }    
  tfout->Close();
}
