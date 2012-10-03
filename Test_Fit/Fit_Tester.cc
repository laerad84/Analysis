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
#include "TSpline.h"
#include "E14WaveformAnalyzer.h"
#include "E14WaveFitter.h"

int main( int argc ,char** argv ){  
  Int_t Number = -1;
  if( argc != 2 && argc != 3){
    std::cerr << "Arguement Error " << std::endl;  
    return -1; 
  }else if ( argc == 3 ){
    Number = atoi( argv[2] );
  }
  

  //Int_t RunNumber = 4205;
  Int_t RunNumber  = atoi(argv[1]);
  std::cout<< RunNumber << std::endl;
  std::string ROOTFILEWAV= std::getenv("ROOTFILE_WAV");
  std::string Filename   = Form("%s/TEMPLATE_FIT_RESULT_DUMP_%d.root",ROOTFILEWAV.c_str(),RunNumber);

#ifdef DEBUG
  std::cout<< __LINE__ << std::endl;
#endif  
  TApplication* app = new TApplication( "app" , &argc , argv );
  TCanvas* can = new TCanvas("can","",500,1000);
  
  
#ifdef DEBUG
  std::cout<< "DEFINE INPUTFILE" << std::endl; 
#endif
  std::cout<< __LINE__ << std::endl;

  TFile* tf = new TFile(Filename.c_str());
  TTree* tr = (TTree*)tf->Get("Waveform");
  
  Int_t    EventNumber;
  Int_t    ModuleNumber;
  Double_t Waveform[48];
  Short_t  TimeInfo[48];
  /*
  Double_t ChisqNDF;
  Double_t PeakTime;
  Double_t HHTime;
  Double_t Height;
  Double_t Pedestal;
  */
  {
    tr->SetBranchAddress("EventNumber",&EventNumber);
    tr->SetBranchAddress("ModuleNumber",&ModuleNumber);
    tr->SetBranchAddress("Waveform",Waveform);
    tr->SetBranchAddress("TimeInfo",TimeInfo);
    /*
    tr->SetBranchAddress("PeakTime",&PeakTime);
    tr->SetBranchAddress("HHTime",&HHTime);
    tr->SetBranchAddress("Height",&Height);
    tr->SetBranchAddress("Pedestal",&Pedestal);
    */
  }

  TFile* tfTemplate  =new TFile("TEMPLATE_COMPLETE_GRAPH_3pi0RunList_200_400.root");
  TGraph*   grWaveform[2716];
  TSpline3* splWaveform[2716];
  for( int i = 0; i< 2716; i++){
    grWaveform[i] = NULL;
    splWaveform[i] = NULL;
    grWaveform[i] = (TGraph*)tfTemplate->Get(Form("Template_%d",i));
    if( grWaveform[i] != NULL ){
      splWaveform[i] = new TSpline3(Form("spl%d",i), grWaveform[i]);
    }
  }


  #ifdef DEBUG
  std::cout<< __LINE__ << std::endl;
  #endif
  TFile* tfOut           = new TFile(Form("Output_Fit_Tester_%d.root",RunNumber),"RECREATE");
  const int nPointInGraph = 48;
  E14WaveformAnalyzer*  WaveAnalyzer = new E14WaveformAnalyzer(48);
  E14WaveFitter*        WaveFitter   = new E14WaveFitter(48);

  Int_t TotalEvent = tr->GetEntries(); 
  TGraph* grWave = new TGraph();
  
  TTree* trWaveAna = new TTree("WaveAna","");
  Double_t Height;
  Double_t PeakTime;
  Double_t MinTime;
  Double_t Maximum;
  Double_t Minimum;
  Double_t PeakPointMaximum;
  
  Double_t Pedestal;
  Double_t PedestalSigma;
  Double_t HeadMean;
  Double_t TailMean;
  Double_t HeadSigma;
  Double_t TailSigma;
  Double_t BoundaryHead;
  Double_t BoundaryTail;
  Double_t Width;
  Double_t ADC;
  Double_t SlopeDelta;
  Int_t    StartPoint; 
  Double_t FitHeight;
  Double_t FitTime;
  Double_t FitPedestal;
  Double_t ChisqNDF;
  Double_t FitResult;
  Int_t    EntryNumber;
  trWaveAna->Branch("EntryNumber"  ,&EntryNumber  , "EntryNumber/I"  );
  trWaveAna->Branch("RunNumber"    ,&RunNumber    , "RunNumber/I"    );
  trWaveAna->Branch("EventNumber"  ,&EventNumber  , "EventNumber/I"  );
  trWaveAna->Branch("ModuleNumber" ,&ModuleNumber , "ModuleNumber/I" );
  trWaveAna->Branch("Height"       ,&Height       , "Height/D"       );
  trWaveAna->Branch("PeakTime"     ,&PeakTime     , "PeakTime/D"     );
  trWaveAna->Branch("Maximum"      ,&Maximum      , "Maximum/D"      );
  trWaveAna->Branch("Minimum"      ,&Minimum      , "Minimum/D"      );
  trWaveAna->Branch("Pedestal"     ,&Pedestal     , "Pedestal/D"     );
  trWaveAna->Branch("PedestalSigma",&PedestalSigma, "PedestalSigma/D");
  trWaveAna->Branch("HeadMean"     ,&HeadMean     , "HeadMean/D"     );
  trWaveAna->Branch("HeadSigma"    ,&HeadSigma    , "HeadSigma/D"    );
  trWaveAna->Branch("TailMean"     ,&TailMean     , "TailMean/D"     );
  trWaveAna->Branch("TailSigma"    ,&TailSigma    , "TailSigma/D"    );
  trWaveAna->Branch("BoundaryHead" ,&BoundaryHead , "BoundaryHead/D" );
  trWaveAna->Branch("BoundaryTail" ,&BoundaryTail , "BoundaryTail/D" );
  trWaveAna->Branch("Width"        ,&Width        , "Width/D"        );
  trWaveAna->Branch("ADC"          ,&ADC          , "ADC/D"          );
  trWaveAna->Branch("SlopeDelta"   ,&SlopeDelta   , "SlopeDelta/D"   );
  trWaveAna->Branch("StartPoint"   ,&StartPoint   , "StartPoint/I"   );
  trWaveAna->Branch("PeakPointMaximum",&PeakPointMaximum , "PeakPointMaximum/D");
  trWaveAna->Branch("FitTime"      ,&FitTime      , "FitTime/D");
  trWaveAna->Branch("FitPedestal"  ,&FitPedestal  , "FitPedestal/D");
  trWaveAna->Branch("FitHeight"    ,&FitHeight    , "FitHeight/D");
  trWaveAna->Branch("FitResult"    ,&FitResult    , "FitResult/D");
		    


  can->Divide( 1, 2);
  //for( int ievent  =0 ; ievent < TotalEvent; ievent++){

  for( int ievent  =0 ; ievent < 2716*1000 ; ievent++){    
    if( Number >= 0 ){
      tr->GetEntry( Number );
    }else{
      tr->GetEntry(ievent);    
    }
    EntryNumber= ievent;
    FitTime     = 0;
    FitPedestal = 0;
    FitHeight   = 0;
    FitResult   = 0;
      
    grWave->Set(0);
    for( int i = 0; i< nPointInGraph; i++){
      grWave->SetPoint( grWave->GetN(), i*8, Waveform[i]); 
    }
      
    WaveAnalyzer->_GetMeanHead( Waveform, HeadMean    , HeadSigma     );
    WaveAnalyzer->_GetMeanTail( Waveform, TailMean    , TailSigma     );
    WaveAnalyzer->_GetMaximum ( Waveform, Maximum     , PeakTime      );
    WaveAnalyzer->_GetMinimum ( Waveform, Minimum     , MinTime       );
    WaveAnalyzer->_GetPedestal( Waveform, Pedestal    , PedestalSigma );
    WaveAnalyzer->_GetWidth   ( Waveform, BoundaryHead, BoundaryTail  );
    WaveAnalyzer->_GetSumSlope( Waveform, StartPoint  , SlopeDelta    );
    
    Int_t rst = WaveAnalyzer->AnalyzeWaveform( Waveform ); 
    Bool_t Fitrst = false;
    if( (rst & 3) != 0){ 
      std::cerr << "Waveform is Over/Under flow " << std::endl;
    }else{
      Double_t TestHeight;
      Double_t TestPeakTime;
      Double_t TestPedestal; 
      WaveAnalyzer->GetParameters( TestPedestal, TestHeight, TestPeakTime);
      WaveFitter->SetParameter( WaveAnalyzer );
      if( splWaveform[ModuleNumber] != NULL ){ 
	WaveFitter->SetWaveform( splWaveform[ModuleNumber] );
	Fitrst      = WaveFitter->Fit(grWave);	
	FitHeight   = WaveFitter->GetParameter(0);
	FitTime     = WaveFitter->GetParameter(1);
	FitPedestal = WaveFitter->GetParameter(2);
	FitResult   = WaveFitter->GetFitResult();
      }else{
	Fitrst = false; 
      }
    }
    
    Width  = BoundaryTail - BoundaryHead;
    Height = Maximum      - Pedestal;
    WaveAnalyzer->_GetADC( Waveform , ADC);
    Double_t SumUp = 0;
    TGraph* grWaveSum = new TGraph();
    for( int ipoint  = 0; ipoint < 48; ipoint++){
      SumUp += Waveform[ ipoint ] - Pedestal; 
      grWaveSum->SetPoint( ipoint , ipoint*8, SumUp);
    }
    
    
    //if( Height  < 10 && SlopeDelta > 100&& SlopeDelta/Height>16 && PeakTime>(StartPoint+4)*8 && PeakTime< (StartPoint+12)*8 ){
    /*
    if(PeakTime>(StartPoint+4)*8&&
       PeakTime<(StartPoint+12)*8&&
       Height<1000&&
       SlopeDelta<8000&&
       BoundaryHead>0&&
       BoundaryTail<384&&
       PeakTime<250&&
       Height>200 &&
       EventNumber == 3417
       ){
    */

    

      can->cd(1);
      gPad->SetGridx();
      gPad->SetGridy();
      WaveAnalyzer->Draw( Waveform );

      can->cd(2);
      gPad->SetGridx();
      gPad->SetGridy();
      WaveFitter->Draw(grWave);
      can->Update();
      can->Modified();
      getchar();
      getchar();
      

    trWaveAna->Fill();
    WaveAnalyzer->_Clear();
    if( Number >= 0 ){ break ; }

  }
  trWaveAna->Write();
  tfOut->Close();
  //app->Run();
  return 0;
}

