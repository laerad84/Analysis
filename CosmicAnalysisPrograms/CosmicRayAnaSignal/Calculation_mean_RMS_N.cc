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
#include "TMath.h"

#include "IDHandler.h"
#include "E14CosmicOut.h"

int
main( int argc , char** argv){
  const int nChannel = 2716;
  IDHandler* handler = new IDHandler();

  TFile* tf = new TFile("Cosmic_3897.root");
  TTree* tr = (TTree*)tf->Get("CosmicOut");

  E14CosmicOut* CosmicOut = new E14CosmicOut( tr ); 
  std::cout << CosmicOut->fChain << std::endl;
  std::cout << CosmicOut->fChain->GetEntries() << std::endl;
  std::cout<< __LINE__ << std::endl;
  
  TFile* tfout = new TFile("TestTimeOutLogic_3897.root","recreate");
  std::cout<< __LINE__ << std::endl;  


  Double_t *Delta_Mean[2716];
  Double_t *Delta_SquareMean[2716];
  Double_t *Delta_RMS[2716];
  Double_t *Delta_Error[2716];
  Int_t    *Delta_Entries[2716];

  for( int i = 0; i< 2716; i++){
    Delta_RMS[i]        = new Double_t [2716];
    Delta_Error[i]      = new Double_t [2716];
    Delta_Mean[i]       = new Double_t [2716];
    Delta_SquareMean[i] = new Double_t [2716];
    Delta_Entries[i]    = new Int_t    [2716];
  }

  for( int i = 0; i< 2716; i++){
    for( int j  =0; j< 2716; j++){
      Delta_Mean[i][j]       = 0;
      Delta_SquareMean[i][j] = 0;
      Delta_Entries[i][j]    = 0; 
      Delta_Error[i][j]      = 9999;
      Delta_RMS[i][j]        = 0;
    }
  }
  
  Double_t x[2];
  Double_t y[2];
  std::cout << "Entries" << CosmicOut->GetEntries() << std::endl;
  for( int ievent = 0; ievent < CosmicOut->GetEntries(); ievent++){
    if( (ievent % 1000) == 0  && ievent ){ std::cout << "EVENT:" << ievent << std::endl; }        
    CosmicOut->GetEntry(ievent);

    // Except Laser Event // 
    if( CosmicOut->CosmicFit != 1  && CosmicOut->nDigi!=0 && CosmicOut->nDigi < 300 ){ continue; }        
    if( CosmicOut->CalFactor < 0.8 ) { continue; }
    Double_t cos       = TMath::Cos(CosmicOut->theta*TMath::Pi()/180.);
    Double_t sin       = TMath::Sin(CosmicOut->theta*TMath::Pi()/180.);	
    
    for( int idigi = 0 ;idigi < CosmicOut->nDigi; idigi++){            
      if(CosmicOut->CsIdepE[idigi]< 30 ){ continue; }
      if( CosmicOut->CsIFitTiming[idigi] < -500 ){ continue; }
      handler->GetMetricPosition(CosmicOut->CsIID[idigi], x[0], y[0] );
      double idist = TMath::Abs((x[0] -(1./cos*(CosmicOut->roh - y[0]*sin ))) * cos );
      if( idist >50){ continue; }
      for( int jdigi = idigi+1; jdigi < CosmicOut->nDigi; jdigi++){
	
	if(CosmicOut->CsIdepE[jdigi]< 30 ){ continue; }
	if( CosmicOut->CsIFitTiming[jdigi] < -500 ){ continue; }	
	handler->GetMetricPosition(CosmicOut->CsIID[jdigi], x[1], y[1] );
	double jdist = TMath::Abs((x[1] -(1./cos*(CosmicOut->roh - y[1]*sin ))) * cos );
	if( jdist > 50 ){ continue; }

	Double_t Distance  = (x[1]-x[0])*-1*sin + (y[1]-y[0])*cos;
	Double_t TimeDelay = Distance/299.792458;
	Double_t DeltaTime = CosmicOut->CsIFitTiming[idigi] -CosmicOut->CsIFitTiming[jdigi]-TimeDelay;
	//std::cout<< (x[1]-x[0]) << " : " << (y[1]-y[0]) << " : " << TimeDelay << std::endl;
	//mm/sec
	
	Delta_Mean[CosmicOut->CsIID[idigi]][CosmicOut->CsIID[jdigi]]       += DeltaTime;
	//Delta_SquareMean[CosmicOut->CsIID[idigi]][CosmicOut->CsIID[jdigi]] += TMath::Power(DeltaTime,2);
	Delta_Entries[CosmicOut->CsIID[idigi]][CosmicOut->CsIID[jdigi]]++;

      }
    }
  }

  for( int ievent = 0; ievent < CosmicOut->GetEntries(); ievent++){
    if( (ievent % 1000) == 0  && ievent ){ std::cout << "EVENT:" << ievent << std::endl; }        
    CosmicOut->GetEntry(ievent);

    // Except Laser Event // 
    if( CosmicOut->CosmicFit != 1  && CosmicOut->nDigi!=0 && CosmicOut->nDigi < 300 ){ continue; }        
    if( CosmicOut->CalFactor < 0.8 ) { continue; }
    Double_t cos       = TMath::Cos(CosmicOut->theta*TMath::Pi()/180.);
    Double_t sin       = TMath::Sin(CosmicOut->theta*TMath::Pi()/180.);	
    
    for( int idigi = 0 ;idigi < CosmicOut->nDigi; idigi++){            
      if(CosmicOut->CsIdepE[idigi]< 30 ){ continue; }
      if( CosmicOut->CsIFitTiming[idigi] < -500 ){ continue; }
      handler->GetMetricPosition(CosmicOut->CsIID[idigi], x[0], y[0] );
      double idist = TMath::Abs((x[0] -(1./cos*(CosmicOut->roh - y[0]*sin ))) * cos );
      if( idist >50){ continue; }
      for( int jdigi = idigi+1; jdigi < CosmicOut->nDigi; jdigi++){
	
	if(CosmicOut->CsIdepE[jdigi]< 30 ){ continue; }
	if( CosmicOut->CsIFitTiming[jdigi] < -500 ){ continue; }	
	handler->GetMetricPosition(CosmicOut->CsIID[jdigi], x[1], y[1] );
	double jdist = TMath::Abs((x[1] -(1./cos*(CosmicOut->roh - y[1]*sin ))) * cos );
	if( jdist > 50 ){ continue; }

	Double_t Distance  = (x[1]-x[0])*-1*sin + (y[1]-y[0])*cos;
	Double_t TimeDelay = Distance/299.792458;
	Double_t DeltaTime = CosmicOut->CsIFitTiming[idigi] -CosmicOut->CsIFitTiming[jdigi]-TimeDelay;
	//std::cout<< (x[1]-x[0]) << " : " << (y[1]-y[0]) << " : " << TimeDelay << std::endl;
	//mm/sec
	
	//Delta_Mean[CosmicOut->CsIID[idigi]][CosmicOut->CsIID[jdigi]]       += DeltaTime;
	DeltaTime = DeltaTime-Delta_Mean[CosmicOut->CsIID[idigi]][CosmicOut->CsIID[jdigi]]/Delta_Entries[CosmicOut->CsIID[idigi]][CosmicOut->CsIID[jdigi]];
	Delta_SquareMean[CosmicOut->CsIID[idigi]][CosmicOut->CsIID[jdigi]] += TMath::Power(DeltaTime,2);
	//Delta_Entries[CosmicOut->CsIID[idigi]][CosmicOut->CsIID[jdigi]]++;

      }
    }
  }














  TTree*   trout = new TTree("TimeDeltaCosmic","");
  Int_t    IDFirst;
  Int_t    IDSecond;
  Int_t    Entries; 
  Double_t Mean;
  Double_t RMS;
  Double_t Error;
  trout->Branch("IDFirst" , &IDFirst ,"IDFirst/I");
  trout->Branch("IDSecond", &IDSecond,"IDSecond/I");
  trout->Branch("Mean"    , &Mean    ,"Mean/D");
  trout->Branch("RMS"     , &RMS     ,"RMS/D");
  trout->Branch("Error"   , &Error   ,"Error/D");
  trout->Branch("Entries" , &Entries ,"Entries/I");

  for( int iID = 0; iID < 2716; iID++){
    for( int jID = 0; jID < 2716; jID++ ){
      IDFirst  = iID;
      IDSecond = jID;
      Entries = (Delta_Entries[iID][jID]);
      if( Entries != 0){
	Mean     = Delta_Mean[iID][jID]/Entries; 
	//RMS      = TMath::Sqrt( Delta_SquareMean[iID][jID] -TMath::Power( Delta_Mean[iID][jID], 2)) )/Entries;
	RMS      = TMath::Sqrt( Delta_SquareMean[iID][jID] )/Entries;
	Error=  RMS / TMath::Sqrt( (double)(Entries) );
      }else{
	Mean = 0;
	RMS  = 10000;
	Error= 1e8;
      }
      trout->Fill();
    }
  }    

  trout->Write();
  tfout->Close();

}
