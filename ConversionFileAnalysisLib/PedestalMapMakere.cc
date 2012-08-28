#include <iostream>
#include <iomanip>
#include <fstream>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include <list>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TApplication.h"
#include "GeneralMacros.h"
#include "GeneralTypes.h"
#include "E14ConvReader.h"
#include "E14MapReader.h"
#include "E14IDHandler.h"
#include "E14ConvWriter.h"

const int nCrate = 11;
int
main( int argc, char* argv[] ){

  if( argc != 2 ){
    std::cerr << "Number of argument must be 2 " << std::endl;
    std::cerr << Form("%s [RunNUmber]",argv[0])  << std::endl;
  }

  Int_t RunNumber = atoi( argv[1] );

  TApplication* app = new TApplication("app",&argc, argv);
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV");
  
  std::cout<< ANALIBDIR   << std::endl;
  std::cout<< CONVFILEDIR << std::endl;
  std::cout<< SUMFILEDIR  << std::endl;
  std::cout<< WAVEFILEDIR << std::endl;


  TFile* tf[nCrate];
  E14ConvReader* conv[nCrate]; 
  for( int iCrate = 0; iCrate < nCrate; iCrate++){
    std::cout << Form( "%s/crate%d/run%d_conv.root", CONVFILEDIR.c_str() ,iCrate, RunNumber ) << std::endl;
    tf[iCrate] = new TFile(Form( "%s/crate%d/run%d_conv.root", CONVFILEDIR.c_str() ,iCrate, RunNumber ));
    conv[iCrate] = new E14ConvReader((TTree*)tf[iCrate]->Get("EventTree"));
  }
  
  /*
  TFile tfTemplate = new TFile("TEMPLATE_SPLINE_500_1000.root");
  TGraphErrors* tempGr[2716];
  TSpline3*     tempSpl[2716];
  for( int ispline = 0; ispline < 2716; ispline++){
    tempGr[ispline] = NULL;
    tempSpl[ispline] = NULL;
    tempGr[ispline] = (TGraphErrors*)tfTemplate->Get(Form("Template_Graph_%d",ispline));
  }
  */

  TFile * tfout = new TFile(Form("%s/Pedestal_RMS_%d.root", WAVEFILEDIR.c_str(), RunNumber), "recreate");
  TTree * trout = new TTree("Dummy","Dummy");
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(), RunNumber),
					    trout);

  std::cout<< "Set Map" << std::endl;
  tfout->cd();
  wConv->AddModule("Csi");
  wConv->AddModule("CC03");
  wConv->AddModule("OEV");
  wConv->AddModule("CV");
  wConv->AddModule("Cosmic");
  wConv->AddModule("Laser");
  wConv->AddModule("Etc");
  
  wConv->Set();
  wConv->SetMap();
  wConv->Branch();
  
  std::cout<< "Check Entries " << std::endl;
  int nEntries = conv[0]->GetEntries(); 
  for( int iCrate = 0; iCrate < nCrate; iCrate++){
    std::cout<< conv[iCrate]->GetEntries() << std::endl;
  }

  for( int iCrate = 0; iCrate < nCrate; iCrate++){
    if( nEntries != conv[iCrate]->GetEntries()){
      std::cout<< "Entries is Different" << std::endl;
    }
  }

  TGraph* grPedestal = new TGraph();
  TGraph* grPedestalRMS = new TGraph(); 
  grPedestal->SetName("grPedestal");
  grPedestal->SetTitle("Pedestal;ID;Pedestal");
  grPedestalRMS->SetName("grPedestalRMS");
  grPedestalRMS->SetTitle("PedestalRMS;ID;Pedestal");
  
  TH1D* hisPedestalDistrib[2716];
  TH1D* hisPedestalRMSDistrib[2716];
  for( int ich = 0; ich < 2716; ich ++){
    hisPedestalDistrib[ich] = new TH1D(Form("hisPedestal%d",ich), Form("hisPedstal%d",ich),100,0,1000);
    hisPedestalRMSDistrib[ich] = new TH1D(Form("hisPedestalRMS%d",ich), Form("hisPedestalRMS%d",ich),1000,0,20);
  }
  
  TGraph* grTemp = new TGraph();
  TCanvas* can = new TCanvas("can","can",800,800);
  for( int ievent  = 0; ievent < nEntries; ievent++){

    wConv->InitData();
    for( int iCrate = 0;iCrate < nCrate; iCrate++){
      conv[iCrate]->GetEntry(ievent);
    }
    wConv->m_RunNo = RunNumber;
    wConv->m_EventNo = ievent; 
    
    int iCsiMod = 0; 
    int nSubCsiModule = wConv->GetNsubmodule(iCsiMod);
    if( nSubCsiModule  <= 0 ){ continue; }
    for( int iSubMod = 0; iSubMod< nSubCsiModule; iSubMod++){
      grTemp->Set(0);
      const Int_t nSampPoints = 7;      
      Double_t MinSigma    = 10000;
      Double_t MinPedestal = -999999;
      Double_t SumUp       = 0.; 
      if( wConv->SetGraph( iCsiMod, iSubMod, conv, grTemp) == 0 ){ 
	continue;
      }
      Int_t nSample  = grTemp->GetN();
      /*
	can->cd();
	grTemp->Draw("AP");
	can->Update();
	can->Modified();
	getchar();
      */

      

      //for( int ipoint  = 0; ipoint < nSample - nSampPoints +1; ipoint++){
      for( int ipoint = 0; ipoint < 1; ipoint++){
	Double_t Sigma = 0.;
	Double_t Pedestal = 0.;
	SumUp += grTemp->GetY()[ipoint];
	for( int isubpoint  = 0; isubpoint < nSampPoints; isubpoint++){
	  Pedestal += grTemp->GetY()[ipoint+isubpoint];
	  Sigma    += TMath::Power(grTemp->GetY()[ipoint+isubpoint],2);
	  
	}
	Pedestal = Pedestal / nSampPoints;
	Sigma   = TMath::Sqrt( Sigma/nSampPoints - (Pedestal*Pedestal) );
	//std::cout<< grTemp->GetY()[ipoint] << " : "  << Pedestal << " : " << Sigma << std::endl;
	if( Sigma < MinSigma ){ 
	  MinSigma = Sigma;
	  MinPedestal = Pedestal; 
	}
      }
      



      SumUp = SumUp- MinPedestal*(grTemp->GetN());
      
      if( SumUp > 300 ) { continue; }
      //std::cout << MinSigma << " : " << MinPedestal << std::endl;
      hisPedestalDistrib[iSubMod]->Fill(MinPedestal);
      hisPedestalRMSDistrib[iSubMod]->Fill(MinSigma);
    }
  }
  for( int ich = 0; ich< 2716; ich++){
    grPedestal->SetPoint( ich, ich , hisPedestalDistrib[ich]->GetMean());
    grPedestalRMS->SetPoint( ich, ich, hisPedestalRMSDistrib[ich]->GetMean());
  }
  for( int ich  = 0; ich < 2716; ich++){
    hisPedestalDistrib[ich]->Write();
    hisPedestalRMSDistrib[ich]->Write();
  }
  grPedestal->Write();
  grPedestalRMS->Write();
  tfout->Close();
  //app->Run();
  return 0;
}

