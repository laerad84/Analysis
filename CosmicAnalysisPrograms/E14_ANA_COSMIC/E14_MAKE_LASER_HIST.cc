#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TFile.h"

#include "E14ReadSumFile.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "HoughCsI.h"

const Double_t COSMIC_THRESHOLD[20] = {10000,10000,10000,10000,10000,
				       10000,10000,10000,10000,10000,
				       10000,10000,10000,10000,10000,
				       10000,10000,10000,10000,10000};


int
main(int argc, char** argv){
  gStyle->SetPalette(1);
  // ARGV[0]  <InputROOTFile> <OutputROOTFile> [<Vision>]
  std::string InputFile;
  std::string OutputFile;


  if( argc == 3){
    InputFile = argv[1];
    OutputFile = argv[2];
  }else{
    std::cerr << "ARGUEMENT ERROR" << std::endl;
    return -1;
  }
  //TApplication* app = new TApplication("app",&argc,argv);
  //
  E14ReadSumFile* Reader = new E14ReadSumFile(0);
  Reader->Add(InputFile.c_str());

  IDHandler* handler = new IDHandler();
  HoughCsI*  hough   = new HoughCsI();
  CsIImage*  image   = new CsIImage(handler);

  TFile* tfOut  = new TFile(OutputFile.c_str(),"Recreate");
  TTree* trOut  = new TTree("CosmicOut","");
  Int_t CosmicTrigger[20]={0};
  Int_t UpperID[5];
  Int_t DownID[5];
  Int_t nHitUp  =  0;
  Int_t nHitDn  =  0;
  Bool_t bCosmicTrigger =  false;
  Int_t TriggerIndex = -1;
  Double_t CsIdepE[N_TOTAL_CSI];
  Double_t CsIID[N_TOTAL_CSI];
  Double_t CalFactor;
  Int_t nDigi;
  Int_t Trigger;

  TH1D* hisLaser[2716];
  for( int i = 0; i< 2716; i++){
    hisLaser[i] = new TH1D(Form("laser%d",i), Form("laser%d",i),600,0,600000);
  }


  trOut->Branch("nDigi",&nDigi,"nDigi/I");
  trOut->Branch("CsIdepE",CsIdepE,"CsIDepE[nDigi]/D");//nDigi;
  trOut->Branch("CsIID",CsIID,"CsIID[nDigi]/D");//nDigi
  trOut->Branch("UpperID",&UpperID,"UpperID[5]/I");
  trOut->Branch("DownID",&DownID,"DownID[5]/I");
  trOut->Branch("nHitUp",&nHitUp,"nHitUp/I");
  trOut->Branch("nHitDn",&nHitDn,"nHitDn/I");
  trOut->Branch("Trigger",&Trigger,"Trigger/I");
  trOut->Branch("CalFactor",&CalFactor,"CalFactor/D");

  /*
  TH1D*  hisCosmic[25][N_TOTAL_CSI];
  TH1D*  hisCosmicCal[25][N_TOTAL_CSI];
  TH1D*  hisCosmicAll[N_TOTAL_CSI];
 
  std::cout << "MakeHistogram" << std::endl;
  for( int j = 0; j< N_TOTAL_CSI; j++){
    hisCosmicAll[j] = new TH1D(Form("hisCosmicAll_%04d",j),
			       Form("hisCosmic_With_All_Trigger_%04d",j),
			       16000,
			       -1600,
			       14400);
    for( int i = 0; i< 25; i++){
      
      hisCosmic[i][j] = new TH1D(Form("hisCosmic_%04d_U%d_D%d",j, i/5, i%5),
				 Form("hisCosmic_%04d_UpsideScinti_%d_DnsideScinti%d",j, i/5, i%5),
				 16000,
				 -1600,
				 14400);
      hisCosmicCal[i][j] = new TH1D(Form("hisCosmicCal_%04d_U%d_D%d",j, i/5, i%5),
				    Form("hisCosmicCal_%04d_UpsideScinti_%d_DnsideScinti%d",j, i/5, i%5),
				    16000,
				    -1600,
				    14400);      
    }
  }
  std::cout <<"End Make Histogram\n" << std::endl;
  */

  
  long nEntries = Reader->GetEntries();
  std::cout << "TotalEntries::" << nEntries  << std::endl;
  
  Double_t roh;
  Double_t theta;
  TH1D* OutputDistribution = new TH1D("Dist","",150,0,600000);
  TH1D* outputForth[4];
  for( int i =0 ; i< 4; i++){
    outputForth[i] = new TH1D(Form("Forth%d",i),"",150,0,600000);
  }
    

  for( int iEntry = 0; iEntry < nEntries; iEntry++){
    Reader->GetEntry(iEntry);    
    if( iEntry%100 == 0 ){
      std::cout << "\r" << iEntry << "/" << nEntries << std::endl;
      std::cout << std::flush;
    }
    /// Init
    nDigi = 0;
    for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
      if( Reader->CsiEne[idigi] > 300 ){
	nDigi++;
      }
    }
    if( nDigi > 1000 ){
      for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
	if( Reader->CsiEne[idigi] > 300 ){
	  hisLaser[Reader->CsiModID[idigi]]->Fill(Reader->CsiEne[idigi]);
	}
      }
    }
  }

  for( int i = 0; i< 2716; i++){
    if( hisLaser[i]->GetMean() != 0){
      /*
      std::cout<<i<< "\t" 
	       <<  hisLaser[i]->GetMean() << "\t" 
	       <<  hisLaser[i]->GetRMS() <<"\t"
	       <<  hisLaser[i]->GetRMS()/ hisLaser[i]->GetMean() << std::endl;
      */
      OutputDistribution->Fill(hisLaser[i]->GetMean());
      image->Fill(i,hisLaser[i]->GetMean());
      double x,y;
      
      handler->GetMetricPosition(i, x,y);
      int indexHist; 
      if( x < 0 && y >=0){
	indexHist = 0;
      }else if( x< 0 && y < 0 ){
	indexHist = 1;
      }else if( x> 0 && y >= 0){
	indexHist  =2;
      }else if( x> 0 && y < 0){
	indexHist = 3;
      }
      outputForth[indexHist]->Fill(hisLaser[i]->GetMean());
      if(hisLaser[i]->GetMean() < 5000 ){
	std::cout<< i << "\t" << hisLaser[i]->GetMean() << std::endl;
      }
    }

    hisLaser[i]->Write();
  }
  gStyle->SetPalette(1);
  TCanvas* can = new TCanvas("can","",1600,800);
  can->Divide(2,1);
  can->cd(1);  
  OutputDistribution->Draw();
  for( int i = 0; i< 4; i++){
    outputForth[i]->SetLineColor(i+1);
    outputForth[i]->Draw("same");

    std::cout <<outputForth[i]->GetMean() << std::endl;
  }
  can->cd(2);
  image->Draw();
  /*
  for( int i = 0; i< N_TOTAL_CSI; i++){
    hisCosmicAll[i]->Write();
    for( int j = 0; j<25; j++){
      hisCosmic[j][i]->Write();
      hisCosmicCal[j][i]->Write();
    }
  }
  */
  //trOut->Write();
  tfOut->Close();

  //app->Run();


}
