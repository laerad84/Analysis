#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
int
main( int argc, char** argv ){
  
  if( argc != 3){
    std::cout<< "Number Of argument is must be 3"      << std::endl;
    std::cout<< argv[0] << " [RunListFile] [RegionID]" << std::endl; 
    return 0;
  }
  std::string InputFileName  = argv[1];  
  int RegionID               = atoi(argv[2]);
  std::ifstream ifs(InputFileName.c_str());
  if( !ifs.is_open() ){ std::cout<< "File is not exist" << std::endl;return -1;}
  std::vector<int> RunNumberList;
  int RunNumber; 
  while( ifs >> RunNumber){
    RunNumberList.push_back( RunNumber );
  }
  std::string ROOTFILE_WAV=std::getenv("ROOTFILE_WAV");
  TFile* tfout = new TFile(Form("TEMPLETE_OUT_HEIGHT_%d.root",RegionID),
			   "recreate");
  TGraphErrors* grTemp = new TGraphErrors();
  TGraphErrors* grWave[2716];
  for( int i = 0; i< 2716; i++ ){
    grWave[i] = new TGraphErrors();
    grWave[i]->SetNameTitle(Form("Waveform_Height_%d_%d",i,RegionID),
			    Form("Waveform_Height_%d_%d",i,RegionID));
  }

  TH2D* hisTempCosmic[2716];
  for( int ich = 0; ich < 2716; ich++){
    hisTempCosmic[ich] = new TH2D( Form("TempHist_Height_%d_%d",ich,RegionID),
				   Form("Templete_Height_%d_%d",ich,RegionID),
				   400,-200,200,
				   200,-0.5,1.5);
  }

  std::cout << "SumUp" << std::endl;
  std::cout << "RunNumberListSize:" <<RunNumberList.size() << std::endl ;
  for( std::vector<int>::iterator it = RunNumberList.begin();
       it != RunNumberList.end();
       it++){
    std::cout<< *it << std::endl;
    TFile* tf = new TFile(Form("%s/TEMPLETE_HEIGHT_%d_%d.root",
			       ROOTFILE_WAV.c_str(),*it,RegionID));    

    for( int ich = 0; ich< 2716; ich++){
      hisTempCosmic[ich]->Add((TH2D*)tf->Get(Form("hisTempCsI_Height_%d_%d",
						  ich,RegionID)));
    }
    tf->Close();
  }
 
  tfout->cd();
  for( int ich  = 0; ich < 2716; ich++){
    TH1D* hisSlice[400];
    grTemp->Set(0);
    if( hisTempCosmic[ich]->GetEntries() != 0 ){
      int ipoint = 0;
      for( int ibin = 0; ibin < 400; ibin++){
	hisSlice[ibin] = hisTempCosmic[ich]->ProjectionY(Form("bin%d",ibin),ibin,ibin+1);
	if( hisSlice[ibin]->GetEntries() < 100 ){ continue; }
	if( hisTempCosmic[ich]->GetBinCenter(ibin) < -150 ){ continue; }
	hisSlice[ibin]->Fit("gaus","Q","",
			    hisSlice[ibin]->GetMean()-2*hisSlice[ibin]->GetRMS(),
			    hisSlice[ibin]->GetMean()+2*hisSlice[ibin]->GetRMS());
	TF1* func = hisSlice[ibin]->GetFunction("gaus");
	grTemp->SetPoint( ipoint ,hisTempCosmic[ich]->GetBinCenter( ibin ),
			  func->GetParameter(1));
	grTemp->SetPointError(ipoint,0,func->GetParameter(2));      
	ipoint++;
      }
    double minRMS   = 10000;
    int    minpoint = 0;
    double pedestal = 0;
    double height   = 0;
    for( int ipoint = 0; ipoint < grTemp->GetN()-7; ipoint++){
      if( grTemp->GetX()[ipoint+6] > 0 ){ continue; }
      double rms        = 0; 
      double partialMax = -1;
      double point      = 0; 
      for( int isubpoint  = 0; isubpoint < 7; isubpoint++){
	if( partialMax < grTemp->GetY()[ipoint+isubpoint]){
	  partialMax = grTemp->GetY()[ipoint+isubpoint];
	}
      }
      for( int isubpoint  = 0; isubpoint < 7; isubpoint++){
	rms += TMath::Power(partialMax-grTemp->GetY()[ipoint +isubpoint],2);
      }
      if( rms < minRMS ){
	minRMS    = rms;
	minpoint  = ipoint;
      }
    }
    for( int ipoint = minpoint; ipoint < minpoint +7; ipoint++){
      pedestal += grTemp->GetY()[ipoint];
    }
    pedestal = pedestal / 7;
    TSpline3* spl = new TSpline3("spl",(TGraph*)grTemp);
    height = spl->Eval(0) - pedestal; 
    for( int ipoint = 0; ipoint < grTemp->GetN(); ipoint++){
      grWave[ich]->SetPoint( ipoint , grTemp->GetX()[ipoint], 
			     (grTemp->GetY()[ipoint]-pedestal)/height);
      grWave[ich]->SetPointError( ipoint, 0 , grTemp->GetEY()[ipoint]/height ); 
    }
      spl->Delete();
    }
    /*
    TSpline3* splWave = new TSpline3(Form("WaveForm_Cosmic_spl_%d",ich),
				    (TGraph*)grWave[ich]);
    splWave->Write();
    */
    grWave[ich]->Write();
    hisTempCosmic[ich]->Write();
  }
  tfout->Close();
}
