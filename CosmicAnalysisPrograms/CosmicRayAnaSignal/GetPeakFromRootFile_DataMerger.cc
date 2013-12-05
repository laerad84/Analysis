#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include <libgen.h>

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"


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

  std::string inputFileList;
  //std::string outputFile;
  /*
  if( argc != 2){
    std::cerr << "Argument Error" << std::endl;
    return -1; 
  }else{
    inputFileList = argv[1];
  }

  char* basenameStr = basename( argv[1] );
  Int_t YMD = 0; 
  sscanf( basenameStr , "CosmicList_%d.txt", &YMD);
  std::cout<< basenameStr << " : " << YMD << std::endl ;

  std::string ROOTFILE_COSMIC = std::getenv("ROOTFILE_COSMIC");

  //TFile* tf =new TFile(inputFile.c_str());
  //TTree* trCosmic = (TTree*)tf->Get("CosmicOut");
  std::cout << "Input TChain Settig" << std::endl; 
  TChain* trCosmic = new TChain("CosmicOut");  
  Int_t RunID;
  std::ifstream ifs( inputFileList.c_str() );
  if( !ifs.is_open() ){ std::cout<< "confirm File path" << std::endl; }
  while ( ifs >> RunID ){
    std::cout<< "Read:" << Form("%s/run%d_cosmic.root",ROOTFILE_COSMIC.c_str(),RunID) << std::endl;
    trCosmic->Add(Form("%s/run%d_cosmic.root",ROOTFILE_COSMIC.c_str(),RunID));
  }
  */
  std::string ROOTFILE_COSMIC = std::getenv("ROOTFILE_COSMIC");
  TChain* trCosmic = new TChain("CosmicOut");  


  Int_t RunIDStart;
  Int_t RunIDEnd;
  std::vector<int> RunList;
  if( argc ==2 ){
    inputFileList = argv[1];
    int tmpRunList;
    std::ifstream ifs( inputFileList.c_str());
    while( ifs >> tmpRunList ){
      RunList.push_back(tmpRunList);
      trCosmic->Add(Form("%s/run%d_cosmic.root",ROOTFILE_COSMIC.c_str(),tmpRunList));
    }

  }else if( argc ==3){
    RunIDStart= atoi( argv[1]);
    RunIDEnd = atoi( argv[2]);
    for( int i = RunIDStart; i <= RunIDEnd; i++){
      trCosmic->Add(Form("%s/run%d_cosmic.root",ROOTFILE_COSMIC.c_str(),i));
    }
  }else{
    return -1; 
  }
  


  Int_t    nDigi;
  Int_t    CsIID[2716];//nDigi
  Double_t CsIdepE[2716];//nDigi
  Double_t CalFactor;
  Int_t    HitUp;
  Int_t    HitDn;
  trCosmic->SetBranchAddress("nDigi",&nDigi);
  trCosmic->SetBranchAddress("CsIID",CsIID);//nDigi
  trCosmic->SetBranchAddress("CsIdepE",CsIdepE);//nDigi
  trCosmic->SetBranchAddress("CalFactor",&CalFactor);
  trCosmic->SetBranchAddress("HitUp",&HitUp);
  trCosmic->SetBranchAddress("HitDn",&HitDn);

  std::cout <<  "Output confirmation" << std::endl; 
  TFile* tfout;
  if( argc ==2 ){    
    const char* base = basename(argv[1]);
    tfout= new TFile(Form("%s/CosmicResult_%s.root",
			  ROOTFILE_COSMIC.c_str(), 
			  base),
		     "recreate");
  }else if( argc ==3 ){
    tfout= new TFile(Form("%s/CosmicResult_%d_%d.root",
			  ROOTFILE_COSMIC.c_str(), 
			  RunIDStart,RunIDEnd),
		     "recreate");
  }

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

  std::cout << "Define Histogram" << std::endl; 
  TH1D* CosmicHist[2716];
  TH1D* CosmicTempHist[2716];
  for( int i = 0; i< 2716; i++){
    CosmicTempHist[i] = new TH1D( Form("hisTempCH%04d",i), Form("CosmicTemp_%04d",i),
				  200,0,100);
  }

  std::cout << "///////////////////////////////////////////////////////" << std::endl;
  std::cout << "Total Entries : " << trCosmic->GetEntries() << std::endl; 
  std::cout << "///////////////////////////////////////////////////////" << std::endl;
  std::cout << "First Loop Start" << std::endl;
  for( int ievent  =0 ;ievent < trCosmic->GetEntries(); ievent++){
    trCosmic->GetEntry(ievent);
    if( CalFactor <= 0.9 ) continue; 
    for( int idigi = 0; idigi < nDigi; idigi++){
      CosmicTempHist[ CsIID[ idigi ] ]->Fill( CsIdepE[ idigi ]*CalFactor );
    }
  }

  Double_t Peak;
  Double_t Sigma; 
  Double_t Norm; 
  std::cout<< " Get First  Result and Set Output Histogram " << std::endl; 
  for( int ich = 0 ; ich< 2716; ich++){
    if( searchPeak ( CosmicTempHist[ich], Norm ,Peak, Sigma, Chi2, Ndf ) ){
      //std::cout<< ich << " : " << Peak << " : " << Sigma << std::endl;
      CosmicHist[ich] = new TH1D( Form( "his_CH%04d", ich ),
				  Form( "his_CH%04d", ich ),
				  200,0,100);
    }else{
      CosmicHist[ich] = new TH1D( Form( "his_CH%04d", ich ),
				  Form( "his_CH%04d", ich ),
				  200,0, 100);
    }
  }

  std::cout<< "Second Loop Start" << std::endl; 
  for( int ievent  =0 ;ievent < trCosmic->GetEntries(); ievent++){
    //std::cout<< ievent << std::endl;
    trCosmic->GetEntry(ievent);
    //if( CalFactor <= 0.9 ) continue; 
    if( (HitUp & 7) == 0 ){ continue; }
    if( (HitDn & 7) == 0 ){ continue; }

    for( int idigi = 0; idigi < nDigi; idigi++){
      CosmicHist[ CsIID[ idigi ] ]->Fill( CsIdepE[ idigi ]);//*CalFactor );
    }
  }

  
  std::cout<< "Get Result " << std::endl;

  for( int ich = 0; ich < 2716; ich ++){ 
    ID= ich;
    if( CosmicHist[ich]->GetEntries() > 100){
      if( searchPeak( CosmicHist[ich], Norm, Peak, Sigma, Chi2, Ndf ) ){
	if( ich >= 2240 ){
	  Peak  = Peak / 2 ;
	  Sigma = Sigma / 2; 
	}

	gr->SetPoint( gr->GetN(), Peak, Sigma );
	grGain->SetPoint( grGain->GetN() ,ich, Peak );
	grGain->SetPointError( grGain->GetN() -1 , 0, Sigma );
	hisGain->Fill(Peak);
	LandauNorm = Norm;
	LandauPeak = Peak;
	LandauSigma= Sigma;
	
      }else{
	gr->SetPoint( gr->GetN(), 0, 0);
	grGain->SetPoint( grGain->GetN(), ich, 0);
	grGain->SetPointError( grGain->GetN()-1, 0, 0 );
	LandauNorm = -1;
	LandauPeak = -1;
	LandauSigma = -1; 
	Chi2 = -1;
	Ndf = -1;
      }      
    }else{
      gr->SetPoint( gr->GetN(), 0, 0);
      grGain->SetPoint( grGain->GetN(), ich, 0);
      grGain->SetPointError( grGain->GetN()-1, 0, 0 );
      LandauNorm = -1;
      LandauPeak = -1;
      LandauSigma = -1; 
      Chi2 = -1;
      Ndf = -1;      
    }
    tr->Fill();
  }

  std::cout<< "Analysis is Done" << std::endl;

  for( int ich  =0 ;ich < 2716; ich++){
    CosmicTempHist[ich]->Write();
    CosmicHist[ich]->Write();
  }

  tr->Write();
  tfout->Close();

}
 
