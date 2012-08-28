
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TTreePlayer.h"
#include "TTreePerfStats.h"


#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "Structs.h"

#include "E14ConvReader.h"
#include "E14IDHandler.h"
#include "E14ConvWriter.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TText.h"
#include "TROOT.h"
#include "Math/WrappedFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMethods.h"
#include "Math/WrappedTF1.h"

#include "E14WaveFitter.h"


int
main( int argc, char** argv){
  if( argc != 2){
    std::cerr << "Please Input RunNumber" << std::endl;
    return -1;
  }
  Int_t RunNumber  = atoi( argv[1]);
  
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB");
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV");
  std::string SUBFILEDIR  = std::getenv("ROOTFILE_SUMUP");


  WaveformFitter* wavFitter = new WaveformFitter( 48, kFALSE );
  E14WaveFitter*  Fitter    = new E14WaveFitter();
  TApplication*   app       = new TApplication("app", &argc, argv);
  
  std::cout << "Setting IO File" << std::endl;

  TFile* tfTemplate = new TFile("TEMPLATE_SPLINE_250_500.root");
  TGraph*    TempGr[2716];
  TSpline3*  TempSpl[2716]; 
  for( int ich = 0; ich < 2716; ich++){
    TempGr[ich] = NULL;
    TempSpl[ich] = NULL;
    TempGr[ich] = (TGraph*)tfTemplate->Get(Form("Template_graph_%d",ich));
    if( TempGr[ich]->GetN() != 0 ){ 
      TempSpl[ich]= new TSpline3(Form("TempSpl%d",ich),TempGr[ich]);
    }else{      
      std::cout<< "Non Exist Channel: " << ich << std::endl;
    }
  }
  tfTemplate->Close();

  TFile* tf  = new TFile(Form("%s/TEMPLATE_FIT_RESULT_DUMP_%d.root",WAVEFILEDIR.c_str(),RunNumber));
  TTree* tr  = (TTree*)tf->Get("Waveform");
  
  Int_t EventNumber;
  Int_t ModuleNumber;
  Double_t Waveform[48];
  Short_t  TimeInfo[48];
  Double_t ChisqNDF;
  Double_t PeakTime;
  Double_t HHTime;
  Double_t Height;
  Double_t Pedestal;

  tr->SetBranchAddress("EventNumber",&EventNumber);
  tr->SetBranchAddress("ModuleNumber",&ModuleNumber);
  tr->SetBranchAddress("Waveform",Waveform);
  tr->SetBranchAddress("TimeInfo",TimeInfo);
  tr->SetBranchAddress("ChisqNDF",&ChisqNDF);
  tr->SetBranchAddress("PeakTime",&PeakTime);
  tr->SetBranchAddress("HHTime",&HHTime);
  tr->SetBranchAddress("Pedestal",&Pedestal);

  TFile* tfout = new TFile(Form("%s/Template_Fit_DUMP_REANALYSIS%d.root",WAVEFILEDIR.c_str(), RunNumber),"recreate");
  

  TH2D* hisChisqHeightTotal = new TH2D("hisChisqTotal","ChisqTotal;Height;ChisqNDF",500,0,10000,200,0,200);
  TH2D* hisChisqHeight[2716];
  for( int ich = 0; ich < 2716; ich++) {
    hisChisqHeight[ich] = new TH2D(Form("hisChisqHeight%d",ich),Form("hisChisqHeight%d",ich),
				   500,0,10000,200,0,200);
  }

 
  Int_t nEntries = tr->GetEntries();
  
  TGraph* grTemp = new TGraph();
  TGraph* grAspect = new TGraph();
  grTemp->SetMarkerStyle(7);
  TCanvas* can = new TCanvas("can","",1000,1000);
  can->Divide( 3,3);
  TGraph* grChisqCorr      = new TGraph();
  TGraph* grChisqDeltaped  = new TGraph();
  TGraph* grChisqDeltaTail = new TGraph();
  TGraph* grChisqAspect    = new TGraph();
  TGraph* grChisqHeight    = new TGraph();
  TGraph* grChisqTime      = new TGraph();
  grChisqTime     ->SetMarkerStyle(6);
  grChisqHeight   ->SetMarkerStyle(6);
  grChisqDeltaped ->SetMarkerStyle(6);
  grChisqDeltaTail->SetMarkerStyle(6);
  grChisqAspect   ->SetMarkerStyle(6);
  grChisqCorr     ->SetMarkerStyle(6);
  for( int ievent = 0; ievent < nEntries; ievent++){
    tr->GetEntry(ievent);

    grTemp->Set(0);
    grAspect->Set(0);

    Double_t SumUp = 0.;
    for( int ipoint = 0; ipoint < 48; ipoint++){
      SumUp += Waveform[ipoint] - Pedestal;
      grTemp->SetPoint(ipoint,TimeInfo[ipoint], Waveform[ipoint]);
    }    
    Fitter->SetWaveform( TempSpl[ModuleNumber] );
    Fitter->InitPar();
    bool Fit = Fitter->Fit(grTemp);
    if( Fit ){
      if( Fitter->GetChisquare()/ Fitter->GetNDF() < 50 ){ continue; }
      Double_t h = Fitter->GetParameter(0);
      if( h > 500 || h < 100 ){ continue; }
      Double_t SumAllSignal   = 0.;
      Double_t SumAllFit      = 0.;
      Double_t Delta_Pedestal = 0.;
      Double_t Delta_tail     = 0.;
      Double_t Chisq          = 0.;
      Double_t Dynamic_PED_SIG= 0.;
      Double_t Dynamic_PED    = 0.;
      Int_t    nPoint_Ped  = 0;
      Int_t    nPoint_Sig  = 0;
      Int_t    nPoint_Tail = 0; 
      TF1* func = Fitter->GetFunction();
      Chisq = Fitter->GetChisquare() / Fitter->GetNDF();
      for( int ipoint  = 0; ipoint < grTemp->GetN(); ipoint++){
	double xPos = grTemp->GetX()[ipoint];
	double yPos = grTemp->GetY()[ipoint]; 

	if( xPos > Fitter->GetParameter(1)-100 &&
	    xPos < Fitter->GetParameter(1)-70  ){
	  Delta_Pedestal += yPos - func->Eval(xPos);
	    Dynamic_PED += yPos;
	    Dynamic_PED_SIG += yPos*yPos;
	    nPoint_Ped++;
	  }	
	  if( xPos > Fitter->GetParameter(1)-100 &&
	      xPos < Fitter->GetParameter(1)+50  ){
	    SumAllSignal += grTemp->GetY()[ipoint]-Fitter->GetParameter(2);
	    SumAllFit += func->Eval(xPos)-Fitter->GetParameter(2);
	    nPoint_Sig++;
	  }	
	  if( xPos > Fitter->GetParameter(1)+50  &&
	      xPos < Fitter->GetParameter(1)+200 ){
	    Delta_tail += yPos-func->Eval(xPos);	  
	    nPoint_Tail++;
	  }
	  if( xPos > Fitter->GetParameter(1) - 100 &&
	      xPos < Fitter->GetParameter(1) + 50  ){
	    grAspect->SetPoint( grAspect->GetN(), yPos, func->Eval(xPos));
	  }


	}
	Dynamic_PED = Dynamic_PED/ nPoint_Ped;
	Dynamic_PED_SIG = TMath::Sqrt(( Dynamic_PED_SIG / nPoint_Ped) - Dynamic_PED*Dynamic_PED );
	std::cout<<"PED:" <<  Dynamic_PED_SIG << std::endl;
	//Chisq = Chisq/Dynamic_PED_SIG;
	std::cout<< Chisq << std::endl ;
	if( TMath::Abs(SumAllSignal/SumAllFit - 1 ) < 0.05 &&
	    TMath::Abs(Delta_Pedestal/nPoint_Ped)   < 3    &&
	    grAspect->GetCorrelationFactor()        > 0.998 ){

	  grChisqTime     ->SetPoint( grChisqTime->GetN()     , Fitter->GetParameter(1)  , Chisq);
	  grChisqHeight   ->SetPoint( grChisqHeight->GetN()   , Fitter->GetParameter(0)  , Chisq);
	  grChisqDeltaped ->SetPoint( grChisqDeltaped->GetN() , Delta_Pedestal/nPoint_Ped, Chisq);
	  grChisqAspect   ->SetPoint( grChisqAspect->GetN()   , SumAllSignal/SumAllFit   , Chisq);
	  grChisqDeltaTail->SetPoint( grChisqDeltaTail->GetN(), Delta_tail/nPoint_Tail   , Chisq);
	  grChisqCorr     ->SetPoint( grChisqCorr->GetN(), grAspect->GetCorrelationFactor(), Chisq);	  
	  
	  std::cout<< Fitter->GetChisquare() / Fitter->GetNDF() << std::endl;
	  if( ievent % 100 == 99 ){
	    can->cd(1);
	    grTemp->GetListOfFunctions()->Delete();
	    grTemp->Draw("AP");
	    func->Draw("same");
	    can->cd(2);
	    grChisqDeltaped->Draw("AP");
	    can->cd(3);
	    grChisqDeltaTail->Draw("AP");
	    can->cd(4);
	    grChisqAspect->Draw("AP");
	    can->cd(5);
	    grChisqHeight->Draw("AP");
	    can->cd(6);
	    grChisqCorr->Draw("AP");
	    can->cd(7);
	    grChisqTime->Draw("AP");
	    can->Update();
	    can->Modified();
	    getchar();
	  }
	}
    }
  }

  app->Run();
}
