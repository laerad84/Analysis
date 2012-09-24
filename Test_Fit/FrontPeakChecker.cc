#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TApplication.h"


Double_t CheckWave_FrontMean( Double_t *Wave );
Bool_t   GetMinimum_FrontPeak( Double_t *Wave ,
			       Int_t& MaxPoint, Double_t& Maximum,
			       Int_t& MinPoint, Double_t& Minimum);
Bool_t  CheckWidth( Double_t* Wave,
		    Double_t& LowBoundary,
		    Double_t& HighBoundary);

int main( int argc ,char** argv ){  
  if( argc != 2 ){
    std::cerr << "Arguement Error " << std::endl;  
    return -1; 
  } 
  //Int_t RunNumber = 4205;
  Int_t RunNumber  = atoi(argv[1]);
  std::string ROOTFILEWAV= std::getenv("ROOTFILE_WAV");
  std::string Filename   = Form("%s/TEMPLATE_FIT_RESULT_DUMP_%d.root",ROOTFILEWAV.c_str(),RunNumber);

  TFile* tf = new TFile(Filename.c_str());
  TTree* tr = (TTree*)tf->Get("Waveform");

  Int_t    EventNumber;
  Int_t    ModuleNumber;
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
  tr->SetBranchAddress("PeakTime",&PeakTime);
  tr->SetBranchAddress("HHTime",&HHTime);
  tr->SetBranchAddress("Height",&Height);
  tr->SetBranchAddress("Pedestal",&Pedestal);

  TFile* tfOut           = new TFile(Form("OutputDistrib_%d.root",RunNumber),"RECREATE");
  TH2D*  hisDelta        = new TH2D("DeltaHis"    ,"",1600,0,16000,100,-50,50);
  TH2D*  hisMinMax       = new TH2D("hisMinMax"   ,"hisMinMax;Minimum;Maximum",48,0,48,48,0,48);
  TH2D*  hisMaxPoint     = new TH2D("hisMaxPeak"  ,"hisMaxPeak;Maximum;Height",48,0,48,1600,0,16000);
  TH2D*  hisWaveform     = new TH2D("hisWaveform" ,"",450,-150,300,150,-0.25,1.25);
  TH2D*  hisWaveWidth    = new TH2D("hisWaveWidth","",500,0,500,1600,0,16000);
  TH2D*  hisWaveWidthLow = new TH2D("hisWaveWidthLow","",500,0,500,160,0,160);

  Int_t TotalEvent = tr->GetEntries(); 
  TGraph* grWave = new TGraph();
  for( int ievent  =0 ; ievent < TotalEvent ; ievent++){
  //for( int ievent  =0 ; ievent < 100 ; ievent++){
    tr->GetEntry(ievent);
    grWave->Set(0);
    Double_t fMean =CheckWave_FrontMean( Waveform );
    Double_t Maximum;
    Double_t Minimum;
    Int_t    MaxPoint;
    Int_t    MinPoint;
    Bool_t   Rst = GetMinimum_FrontPeak( Waveform, MaxPoint, Maximum, MinPoint, Minimum);
    if( Rst &&
	fMean - Minimum  <  9 && 
	fMean - Minimum  > -5 &&
	MaxPoint - MinPoint >= 5 &&
	MaxPoint            < 30 ){


      hisDelta   ->Fill( Maximum - Minimum , fMean - Minimum );
      hisMinMax  ->Fill( MinPoint, MaxPoint);
      hisMaxPoint->Fill( MaxPoint, Maximum - Minimum);
      Double_t LowBoundary;
      Double_t HighBoundary;
      Bool_t RstWidth = CheckWidth(Waveform, LowBoundary, HighBoundary);
      hisWaveWidth->Fill( HighBoundary - LowBoundary,Maximum - Minimum);     
      hisWaveWidthLow->Fill( HighBoundary - LowBoundary,Maximum - Minimum);

      if( HighBoundary - LowBoundary < 100 &&
	  HighBoundary - LowBoundary > 10  &&
	  Maximum - Minimum          > 300 ){
	
	for( int i = 0; i< 48; i++){
	  hisWaveform->Fill(TimeInfo[i] - PeakTime, (Waveform[i]-Pedestal)/Height );
	}
      }

      if( HighBoundary - LowBoundary > 48*8 ){
	std::cout<< HighBoundary << "\t" << LowBoundary <<  "\t" << HighBoundary- LowBoundary << std::endl;
      }

    }
#ifdef DEBUG
    std::cout << "--------------------------------------------------------\n"
	      << "Waveform : " << ievent << "\n"
	      << "MaxPoint : " << MaxPoint<< "\t" << Maximum << "\n"
	      << "MinPoint : " << MinPoint<< "\t" << Minimum << "\n"
	      << "FrontMean: " << fMean   << "\t" << fMean -Minimum << "\n"
	      << "--------------------------------------------------------"
	      << std::endl;
#endif
  }  
  hisMinMax->Write();
  hisDelta->Write();
  hisMaxPoint->Write();
  hisWaveform->Write();
  hisWaveWidth->Write();
  hisWaveWidthLow->Write();
  tfOut->Close();
  return 0; 
}

Double_t CheckWave_FrontMean( Double_t *Wave ){  
  Double_t FrontMean = 0;
  for( int i = 0; i < 4; i++){
    FrontMean += Wave[i];
  }
  FrontMean = FrontMean / 4; 
  return FrontMean; 
}

Bool_t  GetMinimum_FrontPeak( Double_t *Wave ,
			      Int_t& MaxPoint, Double_t& Maximum,
			      Int_t& MinPoint, Double_t& Minimum){
  Minimum = (Double_t)(0xFFFF);
  Maximum = 0;
  MaxPoint = 0;
  MinPoint = 0;
  for( int i = 4; i< 48-5; i++){
    if( Maximum < Wave[i] ){
      Maximum = Wave[i];
      MaxPoint = i;
    }
  }
  for( int i = 4; i< MaxPoint; i++){
    if( Minimum > Wave[i] ){
      Minimum  = Wave[i];
      MinPoint = i;
    }    
  }
  if( Maximum - Minimum < 10 ){
    return false;
  }else{
    return true;
  }
}

Bool_t CheckWidth( Double_t* Wave,
		   Double_t& LowBoundary,
		   Double_t& HighBoundary){

  Int_t    MaxPoint;
  Int_t    MinPoint;
  Double_t Maximum;
  Double_t Minimum;
  Double_t fMean = CheckWave_FrontMean( Wave );
  Bool_t   Rst   = GetMinimum_FrontPeak( Wave, MaxPoint, Maximum, MinPoint, Minimum);
  
  LowBoundary  = 0;
  HighBoundary = 376;
  
  // Find LowBoundary //
  for( int ipoint = MaxPoint; ipoint > MinPoint; ipoint--){
    if(( Wave[ ipoint - 1 ] - Minimum ) < (Maximum - Minimum) /2 &&
       ( Wave[ ipoint ]     - Minimum ) > (Maximum - Minimum)  /2 ){
      LowBoundary = 8*(ipoint-1) + 8*(( Maximum - Minimum )/2  - (Wave[ipoint-1] - Minimum))/(Wave[ipoint]- Wave[ipoint-1]);
      break;
    }else if( Wave[ipoint] == (Maximum - Minimum)/2 ){
      LowBoundary = ipoint*8;
      break;
    }
}
  
  // Find HighBoundary // 
  for( int ipoint = MaxPoint; ipoint < 48; ipoint++){
    if( (Wave[ ipoint +1 ] - Minimum ) < (Maximum - Minimum ) /2 && 
	(Wave[ ipoint ]    - Minimum ) > (Maximum - Minimum ) /2 ){      
      HighBoundary = 8*(ipoint+1) - 8*(( Maximum - Minimum )/2 - (Wave[ipoint+1] - Minimum))/( Wave[ ipoint ] - Wave[ ipoint + 1 ]);      
      
      break;
      
	}else if (( (Wave[ipoint] - Minimum) == (Maximum - Minimum )/2 )){
      HighBoundary = ipoint*8;
      break;
    }
  }
    
  if( LowBoundary == 0 || HighBoundary == 48*8 ){
    return false;
  }else{
    return true;
  }
}
  
			      
