#include "cluster/ClusterFinder.h"

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
#include "TChain.h"
#include "TROOT.h"
#include "Math/WrappedFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMethods.h"
#include "Math/WrappedTF1.h"
#include "TProfile.h"

int
main( int argc, char** argv ){
  if( argc != 4 ){ 
    std::cerr << "Number of argumemt must be 3" << std::endl;
    std::cerr << Form( "%s [list of RunNumber] [Height Low Limit] [Height High Limit]", argv[0])
	      << std::endl;
    return -1;
  }

  std::string Filename  = argv[1];
  Double_t LowLimit  = atof( argv[2] );
  Double_t HighLimit = atof( argv[3] );
  if( LowLimit > HighLimit || LowLimit < 0 ){
    std::cerr << "Please Confirm Limit\n"
	      << "Invalid Number for Limit" << std::endl;
    return -1; 
  }



  //std::ifstream ifs("/home/jwlee/local/Analysis/RunList/EtRunList.txt");
  std::ifstream ifs(Filename.c_str());
  TChain* ch  = new TChain("Waveform");
  std::string Env = std::getenv("ROOTFILE_WAV");
  Int_t RunNumber; 
  Int_t ModuleNumber;
  Double_t Waveform[48];
  Short_t  TimeInfo[48];
  Double_t ChisqNDF;
  Double_t PeakTime;
  Double_t HHTime;
  Double_t Height;
  Double_t Pedestal;
  
  ch->SetBranchAddress("ModuleNumber",(void*)(&ModuleNumber));
  ch->SetBranchAddress("Waveform"    ,(void*)Waveform       );
  ch->SetBranchAddress("TimeInfo"    ,(void*)TimeInfo       );
  ch->SetBranchAddress("PeakTime"    ,(void*)(&PeakTime)    );
  ch->SetBranchAddress("HHTime"      ,(void*)(&HHTime)      );
  ch->SetBranchAddress("Height"      ,(void*)(&Height)      );
  ch->SetBranchAddress("ChisqNDF"    ,(void*)(&ChisqNDF)    );
  ch->SetBranchAddress("Pedestal"    ,(void*)(&Pedestal)    );
  while( ifs >> RunNumber ){
    ch->Add(Form("%s/TEMPLATE_FIT_RESULT_DUMP_%d.root",
		 Env.c_str(), RunNumber)); 
  }
      
  TFile* tfout = new TFile(Form("TEMPLATE_HEIGHT_%d_%d.root",(int)LowLimit,(int)HighLimit),
			   "recreate");
  TH2D* hisTemp_CsI[2716];
  TProfile* prof[2716];
  TGraph*   grTempRaw[2716];
  TGraph*   grTemp[2716];
  TH1D* hisInsertedChannel = new TH1D("hisInsertedChannel","hisInsertedChannel",
				      2716,0,2716);
  for( int i = 0; i< 2716; i++){
    hisTemp_CsI[i] = new TH2D(Form("hisTempCsI_Height_%d_%d_%d",(int)LowLimit,(int)HighLimit, i),
			      Form("hisTempCsI_Height_%d_%d_%d",(int)LowLimit,(int)HighLimit, i),
			      400,-150,250,150,-0.25,1.25);
    grTempRaw[i] = new TGraph();
    grTempRaw[i]->SetName(Form("Template_Raw_%d",i ));
    grTempRaw[i]->SetTitle(Form("Template_Raw_%d",i ));
    grTemp[i] = new TGraph();
    grTemp[i]->SetName(Form("Template_%d",i));
    grTemp[i]->SetTitle(Form("Template_%d",i));
      
  }
  
  TGraph* gr = new TGraph();  
  for( int ievent  = 0; ievent < ch->GetEntries(); ievent++){
    ch->GetEntry(ievent);
    //std::cout <<ChisqNDF << std::endl;
    if( ChisqNDF >30 ) {continue; }    
    if(  Height < LowLimit || Height > HighLimit ){ continue; }
    gr->Set(0);
    bool frontSignal = false;
    bool backSignal  = false;
    
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      Double_t Timing = TimeInfo[ipoint] - PeakTime;
      Double_t Signal = (Waveform[ipoint] - Pedestal)/Height;

      /*
      std::cout << ModuleNumber << " : " 
		<< Timing       << " : " 
		<< Signal       << " : " 
		<< Waveform[ipoint] << " : " 
		<< Pedestal     << " : " 
		<< Height       << std::endl; 
      */

      gr->SetPoint( ipoint , Timing, Signal);
       if( Timing < -80 ){ 
	if( Signal > 0.1 ){
	  frontSignal  = true;
	}
      }else if( Timing > 100 ){
	if( Signal > 0.2 ){
	  backSignal = true;
	}
      } 
    }

    if( !frontSignal && !backSignal ){
      for( int ipoint = 0; ipoint < 48; ipoint++ ){
	Double_t Timing = TimeInfo[ipoint] - PeakTime;
	Double_t Signal = (Waveform[ipoint] - Pedestal)/Height;
	hisTemp_CsI[ModuleNumber]->Fill( Timing, Signal );
      }
      hisInsertedChannel->Fill(ModuleNumber);
    }
  }
  std::cout << "Event Analysis is end" << std::endl;
  /*
  for( int ich = 0; ich < 2716; ich++){
    TH1D* hisSlice[400];
    grTempRaw[ich]->Set(0);
    if( hisTemp_CsI[ich]->GetEntries() != 0 ){
      int ipoint  = 0;
      for( int ibin  =0; ibin < 400; ibin++){
	hisSlice[ibin] = hisTemp_CsI[ich]->ProjectionY( Form("bin%d",ibin), 
							ibin, ibin+1);
	if( hisSlice[ibin]->GetEntries() < 16){ continue; }
	if( hisTemp_CsI[ich]->GetBinCenter(ibin)-150){ continue;}
	grTempRaw[ich]->SetPoint( ipoint,
				  hisTemp_CsI[ich]->GetBinCenter(ibin ),
				  hisSlice[ich]->GetMean());
	ipoint++;
      }
    }
    double minRMS = 10000; 
    int    minpoint = 0; 
    double pedestal = 0;
    double height = 0; 
    int    nsample  = 7;
    for( int ipoint = 0; ipoint < grTempRaw[ich]->GetN()-nsample; ipoint++){
      if( grTempRaw[ich]->GetX()[ipoint+nsample-1] > 0 ){ continue ;}
      double rms = 0.;
      double point  = 0; 
      double mean = 0;
      double sqmean = 0;
      
      for( int isubpoint = 0; isubpoint < nsample; isubpoint++){
	mean += grTempRaw[ich]->GetY()[ipoint + isubpoint ];
	sqmean += grTempRaw[ich]->GetY()[ipoint + isubpoint];
      }

      mean   /= nsample;
      sqmean /= nsample;
      
      rms = TMath::Sqrt( sqmean - mean*mean );
      if( rms < minRMS ){
	minpoint  = ipoint;
	minRMS    = rms; 
	pedestal  = mean;
      }

    }

    TSpline3* spl =  new TSpline3("spl",grTempRaw[ich]);
    height = spl->Eval(0) - pedestal;
    for( int ipoint  =0; ipoint < grTempRaw[ich]->GetN(); ipoint++){
      grTemp[ich]->SetPoint( ipoint , grTempRaw[ich]->GetX()[ipoint],
			     (grTempRaw[ich]->GetY()[ipoint] - pedestal)/height);
    }
    delete spl;
    }
  */
  for( Int_t ich  =0; ich < 2716; ich++){
    prof[ich] = hisTemp_CsI[ich]->ProfileX();
    prof[ich]->Write();
    hisTemp_CsI[ich]->Write();
    //grTemp[ich]->Write();    
  }
  hisInsertedChannel->Write();
  tfout->Close();
  
}  

  
  
  
