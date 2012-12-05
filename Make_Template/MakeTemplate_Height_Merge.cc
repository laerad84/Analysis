////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// MakeTemplate_Height_Merge - Merge TemplateFile(ouput of MakeTemplate_Height) with listFile.
// MakeTemplate_Height_Merge [RunList].txt [LowLimit] [HighLimit] 
// You must define environment argument 
// $ROOTFILE_WAV
// NEED FILE LIST
// ListFile  : [RUNList].txt <-be careful...
// InputFile : TEMPLETE_HEIGHT_[RunNumber]_[LowHeightLimit]_[HighHeightLimit].root
// OutputFile: TEMPLETE_HEIGHT_[RUNList]_[LowHeightLimit]_[HighHeightLimit].root
// Format of listFile
// [runNumber]\n
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <libgen.h>

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
  
  if( argc != 4){
    std::cerr<< "Number Of argument is must be 3"      << std::endl;
    std::cerr<< argv[0] << " [RunListFile] [MinHeight] [MaxHeight]" << std::endl; 
    return -1;
  }

  std::string InputFileName  = argv[1];  
  Int_t MinHeight            = atoi(argv[2]);
  Int_t MaxHeight            = atoi(argv[3]);
  if( MinHeight > MaxHeight ){ 
    std::cerr<< "Error" << std::endl; 
    return -1; 
  }

  std::ifstream ifs(InputFileName.c_str());
  if( !ifs.is_open() ){ std::cout<< "File is not exist" << std::endl;return -1;}
  std::vector<int> RunNumberList;
  int RunNumber; 
  while( ifs >> RunNumber){
    RunNumberList.push_back( RunNumber );
  }
  
  std::string Basename = basename( argv[1] );
  std::string ROOTFILE_WAV=std::getenv("ROOTFILE_WAV");
  TFile* tfout = new TFile(Form("TEMPLETE_OUT_HEIGHT_%s_%d_%d.root",
				Basename.substr(0, Basename.size()-4).c_str(),
				MinHeight,
				MaxHeight),
			   "recreate");

  TGraphErrors* grTemp = new TGraphErrors();
  TGraphErrors* grWave[2716];
  TH2D* hisTempCosmic[2716];

  for( int i = 0; i< 2716; i++ ){
    grWave[i] = new TGraphErrors();
    grWave[i]->SetNameTitle(Form("Waveform_Height_%d_%d_%d",MinHeight,MaxHeight,i),
			    Form("Waveform_Height_%d_%d_%d",MinHeight,MaxHeight,i));
    hisTempCosmic[i] = new TH2D( Form("TempHist_Height_%d_%d_%d",MinHeight,MaxHeight,i),
				 Form("Template_Height_%d_%d_%d",MinHeight,MaxHeight,i),
				 450,-150,300,
				 150,-0.25,1.25);
  }  

  std::cout << "SumUp" << std::endl;
  std::cout << "RunNumberListSize:" << RunNumberList.size() << std::endl ;
  for( std::vector<int>::iterator it = RunNumberList.begin();
       it != RunNumberList.end();
       it++){
    std::cout<< *it << std::endl;
    TFile* tf = new TFile(Form("%s/TEMPLETE_HEIGHT_%d_%d_%d.root",
			       ROOTFILE_WAV.c_str(),*it,MinHeight,MaxHeight));    
    if( tf == NULL ){ return -1; }
    
    for( int ich = 0; ich< 2716; ich++){
      hisTempCosmic[ich]->Add((TH2D*)tf->Get(Form("hisTemplateCsI_Height_%d_%d_%d",
						  MinHeight,MaxHeight,ich)));
    }
    tf->Close();
  }
 
  tfout->cd();
  for( int ich  = 0; ich < 2716; ich++){
    std::cout<< ich << std::endl;
    TH1D* hisSlice[400];
    grTemp->Set(0);
    if( hisTempCosmic[ich]->GetEntries() != 0 ){
      int ipoint = 0;
      for( int ibin = 0; ibin < 450; ibin++){
	hisSlice[ibin] = hisTempCosmic[ich]->ProjectionY(Form("bin%d",ibin),ibin,ibin+1);
	if( hisSlice[ibin]->GetEntries() < 100 ){ continue; }
	if( hisTempCosmic[ich]->GetBinCenter(ibin) < -150 ){ continue; }
	
	hisSlice[ibin]->Fit("gaus","Q","",
			    hisSlice[ibin]->GetMean()-1.5*hisSlice[ibin]->GetRMS(),
			    hisSlice[ibin]->GetMean()+1.5*hisSlice[ibin]->GetRMS());
	TF1* func = hisSlice[ibin]->GetFunction("gaus");
	grTemp->SetPoint( ipoint ,hisTempCosmic[ich]->GetBinCenter( ibin ),
			  func->GetParameter(1));
	
	/*
	grTemp->SetPoint( ipoint , hisTempCosmic[ich]->GetBinCenter( ibin ),
			  hisSlice[ibin]->GetMean());
	*/
	//grTemp->SetPointError(ipoint,0,func->GetParameter(2));      
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
      Int_t Wavepoint = 0;
      for( int ipoint = 0; ipoint < grTemp->GetN(); ipoint++){
	if( grTemp->GetY()[ipoint] > 1.25 || grTemp->GetY()[ipoint] < -0.25 ){
	  continue;
	}
	grWave[ich]->SetPoint( Wavepoint , grTemp->GetX()[ipoint], 
			       (grTemp->GetY()[ipoint]-pedestal)/height);
	grWave[ich]->SetPointError( Wavepoint, 0 , grTemp->GetEY()[ipoint]/height ); 
	Wavepoint++;
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
