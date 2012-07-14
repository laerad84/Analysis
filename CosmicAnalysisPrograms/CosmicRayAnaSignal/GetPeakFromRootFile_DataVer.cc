#include <iostream>
#include <fstream>

#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TSpectrum.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TApplication.h"
#include "TChain.h"

#include "CsI_Module.h"

bool searchPeak(TH1D* hisCosmic, Double_t& Norm, Double_t& Peak, Double_t& Sigma, Double_t& chi2, Int_t& ndf){
  TSpectrum* spec    = new TSpectrum(10,10);
  Int_t nPeak        = spec->Search(hisCosmic,10,"",0.2);
  if( nPeak != 0 ){
    Float_t* xpeaks  = spec->GetPositionX();
    Float_t* ypeaks  = spec->GetPositionY();
    
    Float_t PeakX   =0;
    Float_t lowerX  =0;
    Float_t CurrentX=0;
    Float_t PeakY   =0;
    Float_t lowerY  =0;
    Float_t CurrentY=0;
    
    for( int ipeak = 0; ipeak <nPeak; ipeak++){
      CurrentY = ypeaks[ipeak];
      CurrentX = xpeaks[ipeak];
      if( CurrentX > 4 && CurrentX > hisCosmic->GetRMS()){
	if( CurrentY > lowerY ){
	  lowerY = CurrentY;
	  lowerX = CurrentX;
	}
      }          
    }
    
    PeakX = lowerX;
    PeakY = lowerY;
    /*
    if(PeakX < 300){    
      Peak  = -1;
      Sigma = -1;
      Norm  = -1;
      chi2  = -1;
      ndf   = -1;
      return false;
    }
    */
    //TF1* fitFunction = new TF1("fitFunc","landau");
    Double_t LowLimit = PeakX*0.8;
    Double_t HighLimit= PeakX*1.4;
    if( 2*hisCosmic->GetRMS() > PeakX*0.4){
      HighLimit = PeakX + 2.0*hisCosmic->GetRMS();
    }

    TFitResultPtr result = hisCosmic->Fit("landau","Q","",
					  LowLimit,
					  HighLimit);
    int fitresult = result;
    if( fitresult == -1 ){
      Peak  = -1;
      Sigma = -1;
      Norm  = -1;
      chi2  = -1;
      ndf   = -1;
      return false;
    }else{      
      TF1* fitFunction = hisCosmic->GetFunction("landau");
      Norm = fitFunction->GetParameter(0);
      Peak = fitFunction->GetParameter(1);
      Sigma= fitFunction->GetParameter(2);
      chi2 = fitFunction->GetChisquare();
      ndf  = fitFunction->GetNDF();
      if( Peak <  PeakX-0.5*hisCosmic->GetRMS()){
	
	Peak  = -1;
	Sigma = -1;
	Norm  = -1;
	chi2  = -1;
	ndf   = -1;
	
	return false;
      }
      return true;
    }
  }else{
    Peak  = -1; 
    Sigma = -1;
    Norm  = -1;
    chi2  = -1;
    ndf   = -1;
    return false;
  }
}


int 
main( int argc, char** argv){

  std::string inputFile;
  std::string outputFile;
  
  if( argc != 3){
    std::cerr << "Argument Error" << std::endl;
    return -1; 
  }else{
    inputFile  = argv[1];
    outputFile = argv[2];
  }

  TApplication* app = new TApplication("app",&argc, argv);
  TCanvas* can  =  new TCanvas("can","can",0,0,800,800);
  can->Draw();


  CsI_Module* csi = new CsI_Module("csi");

  //TFile* tf =new TFile(inputFile.c_str());
  //TTree* trCosmic = (TTree*)tf->Get("CosmicOut");
  TChain* trCosmic = new TChain("CosmicOut");
  Int_t RunID;
  std::ifstream ifs( inputFile.c_str() );
  while ( ifs >> RunID ){
    trCosmic->Add(Form("Cosmic_%d.root",RunID));
  }

  TFile* tfout = new TFile(outputFile.c_str(),"recreate");
  TGraph* gr = new TGraph();
  TGraphErrors* grGain = new TGraphErrors();
  TH1D* hisGain = new TH1D("hisGain","",200,0,40);
  TTree* tr = new TTree("GainFitPar","Cosmic Ray Peak");
  
  Double_t LandauNorm;
  Double_t LandauPeak;
  Double_t LandauSigma;
  Double_t Chi2;
  Int_t    ID;
  Int_t    Ndf;
  Int_t    Nentries;
  tr->Branch("Norm" ,&LandauNorm,"Norm/D");
  tr->Branch("Peak" ,&LandauPeak,"Peak/D");
  tr->Branch("Sigma",&LandauSigma,"Sigma/D");
  tr->Branch("ID"   ,&ID,"ID/I");
  tr->Branch("Chi2" ,&Chi2,"Chi2/D");
  tr->Branch("Ndf"  ,&Ndf,"Ndf/I");

  gr->SetNameTitle("cosmic","cosmic");
  grGain->SetNameTitle("gain","gain");
  TH1D* CosmicHist = NULL;
  TH1D* CosmicTempHist =NULL;
  for( int i = 0; i < 2716 ; i++){    
    //for( int i = 0; i < 100 ; i++){    
    std::cout << i << std::endl;
    
    CosmicTempHist = new TH1D( Form("hisTEMPCH%04d",i),Form("CosmicTemp_%04d",i),
			       200,-200,800);
    trCosmic->Project(CosmicTempHist->GetName(),
		      "CsidepE*CalFactor",
		      Form("CsIID==%d && CalFactor>0.9",i));
      
    ID = i;
    Double_t Peak;
    Double_t Sigma;
    Double_t Norm;    
    if( searchPeak( CosmicTempHist,Norm, Peak, Sigma ,Chi2,Ndf) ){
      std::cout<< i << " : " <<  Peak << " : " << Sigma << std::endl;
      CosmicHist = new TH1D(Form("his_CH%04d",i),Form("Cosmic_%04d",i),
			    160,0,8*Peak);
      trCosmic->Project(CosmicHist->GetName(),
			"CsidepE*CalFactor",
			Form("CsIID==%d && CalFactor>0.9",i));
      if( searchPeak( CosmicHist,Norm, Peak, Sigma ,Chi2,Ndf) ){      	
	if( i >= 2240){
	  Peak  = Peak /2;
	  Sigma = Sigma/2;
	}
	
	gr->SetPoint(gr->GetN(), Peak, Sigma);
	grGain->SetPoint(grGain->GetN(), i, Peak);
	grGain->SetPointError(grGain->GetN()-1, 0, Sigma);
	hisGain->Fill(Peak);      
	std::cout<< Norm << std::endl;
	LandauNorm  = Norm;
	LandauPeak  = Peak;
	LandauSigma = Sigma; 

	csi->Fill( i, LandauPeak );

      }else{
	gr->SetPoint(gr->GetN(), 0, 0);
	grGain->SetPoint(grGain->GetN(), i, 0);
	grGain->SetPointError(grGain->GetN()-1, 0, 0);
	LandauNorm = -1;
	LandauPeak = -1;
	LandauSigma = -1;
	Chi2 = -1;
	Ndf  = -1;
      }
    }else{
      gr->SetPoint(gr->GetN(), 0, 0);
      grGain->SetPoint(grGain->GetN(), i, 0);
      grGain->SetPointError(grGain->GetN()-1, 0, 0);
      LandauNorm = -1;
      LandauPeak = -1;
      LandauSigma = -1;
      Chi2 = -1;
      Ndf  = -1;
    }
    tr->Fill();
    
    //can->Modified();
    //can->Update();
    //getchar();
    
    if( CosmicTempHist ){
      CosmicTempHist->Write();
      CosmicTempHist->Clear();
      CosmicTempHist = NULL;
    }
    if( CosmicHist ){
      CosmicHist->Write();
      CosmicHist->Clear();
      CosmicHist = NULL;
    }
  }

  std::cout<< "Analysis is Done" << std::endl;

  TCanvas* canNew = new TCanvas("canNew","",1600,800);    
  canNew->Divide(2,1);
  canNew->cd(1);
  gr->Draw("AP");
  canNew->cd(2);

  csi->Draw("colz");

  can->Modified();
  can->Update();
  gr->Write();

  grGain->Write();
  hisGain->Write();

  tr->Write();
  tfout->Close();
  
  app->Run();

}
 
