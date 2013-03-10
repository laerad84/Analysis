////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// DON'T CONCERN ABOUT SPELL MISS!!!                                                              //
// MakeTemplate_Height - Making Template histogram with function-waveform-fitter                  //
// ./MakeTemplate_Height [RunNumber] [MinimumHeight] [MaximumHeight]                              //
// You must define environment arguement ...                                                      //
// $ROOTFILE_CONV  : convfile( which has waveform data ) directory                                //
// $ROOTFILE_WAV   : output folder                                                                //
// $ROOTFILE_SUMUP : for read channel map                                                         //
// NEED FILE LIST                                                                                 //
// InputFile : convFile( which has waveform information                                           //
// OutputFile: TEMPLETE_HEIGHT_[RunNumber]_[MinimumHeight]_[MaximumHeight].root                   //
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

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

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "E14WaveformAnalyzer.h"
#include "E14WaveFitter.h"
#include "Structs.h"

#include "E14ConvReader.h"
#include "E14IDHandler.h"
#include "E14ConvWriter.h"
#include "sys/stat.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "Style.h"
struct stat st;

const int nCrate = 11;

const Double_t COSMIC_THRESHOLD[20] = {100,100,100,100,100,
				       100,100,100,100,100,
				       100,100,100,100,100,
				       100,100,100,100,100};

double PeakRegion_1pct[3][16]={{6313 ,7947 ,8946 ,9935 ,10509,10919,11313,11698,
				12016,12329,12576,12799,12995,13171,13411,14232},
			       {6942 ,8102 ,9214 ,10085,10616,11030,11437,11835,
				12189,12483,12739,12950,13118,13304,13708,14000},
			       {6359 ,8078 ,8952 ,10020,10531,10892,11306,11715,
				12100,12402,12636,12888,13066,13216,13406,13701}};
double PeakRegion_2pct[3][8]={{7947 ,9935 ,10919,11698,12329,12799,13171,14232},
			      {8102 ,10085,11030,11835,12483,12950,13304,14000},
			      {8078 ,10020,10892,11715,12402,12888,13216,13701}};

int  main(int argc,char** argv)
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Preparation 
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  SetStyle();
  if( argc != 4 ){ 
    std::cout << argv[0] << "\t"
	      << "[RunNumber]\t" << "[EventNumber]\t" << "[ChannelNumber]" << std::endl;
    return -1;
  }

  // argv[1]:channelListName
  // argv[2]:RunNumber=
  
  Int_t RunNumber             = atoi( argv[1] );
  Int_t EventNumber           = atoi( argv[2] );
  Int_t ChannelNumber         = atoi( argv[3] );

  // GetEnvironment // 
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::string SUMFILEDIR_1= "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/2012_FEB/Sumup";
  
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;

  std::string iFileForm = "%s/crate%d/run%d_conv.root";
  std::string mapFileDir;
  if( stat(Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),&st) == 0 ){
    mapFileDir = SUMFILEDIR;
  }else{
    mapFileDir = SUMFILEDIR_1;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Setting IO & Analysis Class 
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  // Setting  Classes //
  WaveformFitter*      wavFitter    = new WaveformFitter(48, kFALSE);  
  E14WaveFitter*       Fitter       = new E14WaveFitter();
  E14WaveformAnalyzer* wavAnalyzeer = new E14WaveformAnalyzer();
  TFile* tf[nCrate];
  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
  }

  //TFile* tfout = new TFile(Form("run%d_wav.root",RunNumber),"recreate");
  //TFile* tfout = new TFile(Form("%s/WaveformExtract_%d.root",
  //WAVEFILEDIR.c_str(),RunNumber),
  //"recreate");

  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree");   
  
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",mapFileDir.c_str(),RunNumber),
					    trout);
  //tfout->cd();
  {
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
    std::cout<< "Check Entries" << std::endl;    
    for( int icrate = 0; icrate < nCrate; icrate++){
      std::cout<< conv[icrate]->GetEntries() << std::endl;
    }
    int nentries = conv[0]->GetEntries();  
    for( int icrate = 1; icrate < nCrate; icrate++){
      if( nentries != conv[icrate]->GetEntries() ){
	std::cout << "Entries is Different" << std::endl;
      }
    }
  }
  // Additional Branch
  Double_t LaserTrigWF[48];
  Double_t ChannelOutWF[2716][48];
  Double_t SignalTime[48];
  trout->Branch("LaserTrigWF",LaserTrigWF,"LaserTrigWF[48]/D");
  trout->Branch("ChannelOutWF",ChannelOutWF,"ChannelOutWF[2716][48]/D");
  trout->Branch("SignalTime",SignalTime,"SignalTime[48]/D");
  for( int iPoint  =0; iPoint < 48; iPoint ++){
    SignalTime[iPoint] = 8*iPoint;
  }
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////// 
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  TH2D* hisCosmicTemplete[20];
  TH2D* hisLaserTemplete[5];
  TH2D* hisGammaTemplete[20][8];
  TApplication* app = new TApplication("app", &argc , argv );  
  TCanvas* can      = new TCanvas( "can ", "Canvas" ,800,800);
  TGraph* gr        = new TGraph();
  gr->SetMarkerStyle(21);
  
  int CsiModuleID    = wConv->GetModuleID("Csi");
  int CosmicModuleID = wConv->GetModuleID("Cosmic");
  int CVModuleID     = wConv->GetModuleID("CV");
  int LaserModuleID  = wConv->GetModuleID("Laser");
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set Trigger Map Cosmic Laser CV 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  const int nCVModule     = 10;
  const int nCosmicModule = 20; 
  const int CosmicArr[20] = {4 ,5 ,2 ,3 ,6 ,7 ,0 ,1 ,12,13,10,11,14,15,8 ,9 ,16,17,18,19};
  
  if((wConv->ModMap[CVModuleID]).nMod     != 10){ std::cout<< "CV nModule is not equal" << std::endl;}
  if((wConv->ModMap[CosmicModuleID]).nMod != 20){ std::cout<< "Cosmic nModule is not equal" << std::endl;}  
  
  int LaserCFC[3];
  int CVCFC[10][3];  
  int CosmicCFC[20][3];
  for( int subID = 0; subID < 3; subID++ ){    
    LaserCFC[subID] = (wConv->ModMap[LaserModuleID]).Map[0][subID];
  }
  for( int icv = 0; icv < (wConv->ModMap[CVModuleID]).nMod; icv++){
    for( int subID = 0; subID < 3; subID++){
      CVCFC[icv][subID] = (wConv->ModMap[LaserModuleID]).Map[0][subID];
    }
  }
  for( int icosmic = 0; icosmic < (wConv->ModMap[CosmicModuleID]).nMod; icosmic++){
    for( int subID = 0; subID < 3; subID++){
      CosmicCFC[icosmic][subID] = (wConv->ModMap[LaserModuleID]).Map[0][subID];
    }
  }
  
  double CosmicSignal[20];
  double CosmicTime[20];
  double CVSignal[10];
  double CVTime[10];
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  
  std::cout <<"Analysis Started" <<std::endl;
  
  //for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
  for( int ievent  = EventNumber; ievent < EventNumber+1; ievent++ ){
    wConv->InitData();
    for( int icosmic = 0; icosmic < 20; icosmic++){
      CosmicSignal[icosmic] = 0;
      CosmicTime[icosmic]   = 0;
    }
    for( int icv = 0; icv < 10; icv++){
      CVSignal[icv] = 0;
      CVTime[icv]   = 0; 
    }
    
    for( int icrate = 0; icrate < nCrate; icrate++){
      conv[icrate]->GetEntry(ievent);
    } 
    
    wConv->m_RunNo   = RunNumber;
    wConv->m_EventNo = ievent;
    if( ievent %100 == 0 && ievent ){ std::cout<< ievent << "/" << conv[0]->GetEntries() << std::endl;} 
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// 1.Fitting Wave form /// 
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // CsI data : only taking Data out of #1200;
    for( int iMod = 0 ; iMod <1; iMod++){
      int nSubModule = (wConv->ModMap[iMod]).nMod;
      for( int iSubMod = ChannelNumber ; iSubMod < ChannelNumber+1 ; iSubMod++){

	const double HeightMinimum     = 5.;
	const double slopeDeltaMinimum = 50.;
	const double widthThreshold    = 32.;
	const double heightThreshold   = 30.;

	//if( iSubMod %100 != 0 || iSubMod > 1000){ continue;}
	//std::cout<< iSubMod << std::endl;
	int iCrate = 9999;
	int iSlot  = 9999;
	int iCh    = 9999;	
	wConv->GetCFC( iMod, iSubMod, iCrate, iSlot, iCh);	
	if( iCrate == 9999 || iSlot == 9999 || iCh == 9999 ) continue;	
	wConv->SetGraph( iMod, iSubMod, conv, gr);
	for( int iPoint = 0; iPoint < 48 ; iPoint++){
	  ChannelOutWF[iSubMod][iPoint] = gr->GetY()[iPoint];
	}	
	bool fit     = wavFitter->Fit( gr ); 
	std::cout<< wavFitter->GetFunction()->GetChisquare() << std::endl;
	gr->Draw("AP");
	can->Update();
	can->Modified();
	can->Update();
	getchar();
	int chIndex  = (wConv->mod[iMod])->m_nDigi;	 
	if( fit ){ 
	  TF1* fitFunc      = wavFitter->GetFunction();	    
	  double halfHeight = fitFunc->GetParameter(0)/2 + fitFunc->GetParameter(4);
	  double halfTiming = fitFunc->GetX( halfHeight,
					     fitFunc->GetParameter(1)-48, fitFunc->GetParameter(1));
	  TF1* linearFunction = new TF1("func","pol1",halfTiming - 12, halfTiming + 12);
	  gr->Fit( linearFunction, "Q", "", halfTiming -12, halfTiming +12 );
	  double halfFitTiming = linearFunction->GetX( halfHeight, halfTiming -12, halfTiming +12);
	  gr->GetListOfFunctions()->Delete();
	  delete linearFunction;
	  wavFitter->Clear();	  
	}
      }
    }
    //trout->Fill();      
  }
  std::cout<< "end Loop" <<std::endl;
  std::cout<< "All Convert is Done" << std::endl;
  std::cout<< "Close" << std::endl;
  return 0;
}
