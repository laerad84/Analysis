#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>

#include "TROOT.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

#include "TClonesArray.h"
#include "IDHandler.h"
#include "E14CosmicOut.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TProfile.h"
#include "TMath.h"
int
main ( int argc, char** argv ){

  TApplication* app = new TApplication("app", & argc ,argv );
  const int nChannel = 2716;
  IDHandler* handler = new IDHandler();

  /*  
      {}
      {

      TFile* tf = new TFile("TestTimeOutLogic_3897.root");
      TTree* tr = (TTree*)tf->Get("TimeDeltaCosmic");
      Int_t IDFirst;
      Int_t IDSecond;
      Double_t Mean;
      Double_t RMS;
      Double_t Error;
      Int_t Entries;
      tr->SetBranchAddress("IDFirst",&IDFirst);
      tr->SetBranchAddress("IDSecond",&IDSecond);
      tr->SetBranchAddress("Mean",&Mean);
      tr->SetBranchAddress("RMS",&RMS);
      tr->SetBranchAddress("Error",&Error);
      tr->SetBranchAddress("Entries",&Entries);

      TTree* trRead = new TTree("ReadFile","");
      trRead->ReadFile("testNewWORKCompileOffset.txt","CHID/I:Offset/D:OffsetError/D");
  
      int CHID;
      double Offset;
      trRead->SetBranchAddress("CHID",&CHID);
      trRead->SetBranchAddress("Offset",&Offset);

      const int nCH = 2716;
      int CHIDList[nCH]={-1};
      double OffsetList[nCH]={0};
      int CHFlag[nCH] = {0};
      if( trRead->GetEntries() != nCH ){
      std::cout<< "Error" << std::endl; 
      return 0 ;
      }

      for( int i = 0; i< 2716; i++){    
      trRead->GetEntry(i);
      CHIDList[CHID]   = CHID;
      OffsetList[CHID] = Offset;    
      }
  
      TH2D* hisDistrib = new TH2D( "hisDistrib" , "" , 2716,0,2716, 200, -10, 10);

      for( int i = 0; i< tr->GetEntries() ; i++){

      tr->GetEntry(i); 
      //if( Entries <= 100 || Error >= 5 ){ continue; }
      if( Entries == 0 ){ continue; }
      if( IDFirst != 0 ) { continue; }
      hisDistrib->Fill( IDSecond, Mean - OffsetList[IDFirst] + OffsetList[IDSecond]);
      }
      hisDistrib->Draw("colz");
      }
  */

  
    
  TFile* tf = new TFile("Cosmic_3897.root");
  TTree* tr = (TTree*)tf->Get("CosmicOut");
    
  E14CosmicOut* CosmicOut = new E14CosmicOut( tr ); 
  std::cout << CosmicOut->fChain << std::endl;
  std::cout << CosmicOut->fChain->GetEntries() << std::endl;
  std::cout<< __LINE__ << std::endl;
  TTree* trRead = new TTree("ReadFile","");
  trRead->ReadFile("testNewWORKCompileOffset.txt","CHID/I:Offset/D:OffsetError/D");

  int CHID;
  double Offset;
  trRead->SetBranchAddress("CHID",&CHID);
  trRead->SetBranchAddress("Offset",&Offset);
  
  const int nCH = 2716;
  int CHIDList[nCH]={-1};
  double OffsetList[nCH]={0};
  int CHFlag[nCH] = {0};
  if( trRead->GetEntries() != nCH ){
    std::cout<< "Error" << std::endl; 
    return 0 ;
  }
  TH1D* hisDis = new TH1D("hisDis","",200,-20,20);
  TH1D* hisDist= new TH1D("hisDist","",200,-20,20);
  TH2D* hisDistrib = new TH2D( "hisDistrib" , "" , 2716,0,2716, 200, -20, 20);
  {}
  {    
    for( int i = 0; i< 2716; i++){    
      trRead->GetEntry(i);
      CHIDList[CHID]   = CHID;
      OffsetList[CHID] = Offset;    
    }
    
    std::cout<< __LINE__ << std::endl;  
    double x[2],y[2];
    int TestID = 1000;
    if( argc ==2 ){
      TestID = atoi( argv[1] );
    }
    
    
    for( int ievent  =0; ievent < CosmicOut->GetEntries(); ievent++){
      if( ievent %1000 ==0  && ievent ) {
	std::cout<< "Event : " << ievent << std::endl; 
      }
      CosmicOut->GetEntry(ievent);
      if( CosmicOut->CalFactor < 0.8){ continue; }
      if( CosmicOut->CosmicFit != 1   ||
	  CosmicOut->nDigi     == 0   ||
	  CosmicOut->nDigi     >  300 ){ continue; }
      
      Double_t cos       = TMath::Cos(CosmicOut->theta*TMath::Pi()/180.);
      Double_t sin       = TMath::Sin(CosmicOut->theta*TMath::Pi()/180.);	
      
      int IDFirst = -1;
      int IDSecond = -1; 
      double time0 = -500;
      for( int idigi = 0; idigi < CosmicOut->nDigi; idigi++){
	if( CosmicOut->CsIID[idigi] == TestID ){ 
	  IDFirst = TestID;
	  time0 = CosmicOut->CsIFitTiming[idigi];
	  continue; 
	}
      }
      
      handler->GetMetricPosition( TestID ,x[0],y[0] );
      double idist = TMath::Abs((x[0] -(1./cos*(CosmicOut->roh - y[0]*sin ))) * cos );
      if( idist >50){ continue; }
      
      for( int jdigi = 0; jdigi <CosmicOut->nDigi; jdigi++){
	IDSecond = CosmicOut->CsIID[jdigi];
	if( IDSecond == IDFirst ) { continue ; }	
	handler->GetMetricPosition(CosmicOut->CsIID[jdigi], x[1], y[1] );	
	double jdist = TMath::Abs((x[1] -(1./cos*(CosmicOut->roh - y[1]*sin ))) * cos );
	if( jdist > 50 ){ continue; }
	//Double_t DeltaTime = time0 - CosmicOut->CsIFitTiming[jdigi] - OffsetList[IDFirst] + OffsetList[IDSecond];
	Double_t DeltaTime = time0 - CosmicOut->CsIFitTiming[jdigi] ;//- OffsetList[IDFirst] + OffsetList[IDSecond];
	hisDistrib->Fill(IDSecond, DeltaTime);	
	hisDis->Fill(DeltaTime+OffsetList[IDFirst]-OffsetList[IDSecond]);
	hisDist->Fill(DeltaTime);
      }
    }
  }
  
  TProfile* pro = hisDistrib->ProfileX();  
  pro->SetLineColor( 2 );
  pro->SetLineWidth( 2 );
  pro->SetMarkerStyle(7);

  TCanvas* can = new TCanvas("can","",800,800);
  //can->Divide(2,1);
  //can->cd(1);
  hisDistrib->Draw("colz");
  pro->Draw("same");
  //can->cd(2);
  /*
  hisDist->Draw();
  hisDis->SetLineColor(2);
  hisDis->Draw("same");
  pro->Draw("same");
  */
  app->Run();
      
}    
