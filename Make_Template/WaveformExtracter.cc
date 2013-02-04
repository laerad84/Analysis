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
  
  if( argc != 2 ){ 
    std::cout<< argv[0] << "\t" << "[ChannelListName]" << "\t" << "[RunNumber]" << std::endl;
    return -1;
  }

  // argv[1]:channelListName
  // argv[2]:RunNumber=
  Int_t RunNumber             = atoi( argv[1] );

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
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  

  TFile* tf[nCrate];
  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
  }

  //TFile* tfout = new TFile(Form("run%d_wav.root",RunNumber),"recreate");
  TFile* tfout = new TFile(Form("%s/WaveformExtract_%d.root",
				WAVEFILEDIR.c_str(),RunNumber),
			   "recreate");
  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree");   
  
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",mapFileDir.c_str(),RunNumber),
					    trout);
  tfout->cd();
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
  gr->SetMarkerStyle(6);

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

  std::cout <<"Loop " <<std::endl;
  for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
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
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fit waveform with function for all modules, all channels;
    // loop for Modules( CsI, CC03, CV, .... )
    // loop for channels;
    // 1.Fitting waveform for all channels,
    // 2.TriggerDicision ( check Laser event, cosmic event ),
    // 3.No Laser && No Cosmic Event -> Maybe Gamma event,,,,-> Using in Template making..
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////




    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// 1.Fitting Wave form /// 
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // CsI data : only taking Data out of #1200;
    for( int iMod = 0 ; iMod <1; iMod++){
      int nSubModule = (wConv->ModMap[iMod]).nMod;
      for( int iSubMod = 0 ; iSubMod <wConv->GetNsubmodule(0) ; iSubMod++){
	//if( iSubMod %100 != 0 || iSubMod > 1000){ continue;}
	//std::cout<< iSubMod << std::endl;
	int iCrate = 9999;
	int iSlot  = 9999;
	int iCh    = 9999;	
	wConv->GetCFC( iMod, iSubMod, iCrate, iSlot, iCh);	
	// Ignore unmapped channel // 
	if( iCrate == 9999 || iSlot == 9999 || iCh == 9999 ) continue;       
	
	// SetData to Graph // 	
	wConv->SetGraph( iMod, iSubMod, conv, gr);
	for( int iPoint = 0; iPoint < 48 ; iPoint++){
	  ChannelOutWF[iSubMod][iPoint] = gr->GetY()[iPoint];
	}
	
	// Fit graph with waveform Fitter // 
	bool fit     = wavFitter->Fit( gr ); 
	int chIndex  = (wConv->mod[iMod])->m_nDigi;	 
	if( fit ){ 
	  TF1* fitFunc      = wavFitter->GetFunction();	    
	  double halfHeight = fitFunc->GetParameter(0)/2 + fitFunc->GetParameter(4);
	  double halfTiming = fitFunc->GetX( halfHeight,
					     fitFunc->GetParameter(1)-48, fitFunc->GetParameter(1));
	  TF1* linearFunction = new TF1("func","pol1",halfTiming - 12, halfTiming + 12);
	  gr->Fit( linearFunction, "Q", "", halfTiming -12, halfTiming +12 );
	  double halfFitTiming = linearFunction->GetX( halfHeight, halfTiming -12, halfTiming +12);
	  
	  /// Setting Data for write : in this program wConv->Write() had been comment out
	  wConv->mod[iMod]->m_FitHeight[chIndex] = 1;
	  wConv->mod[iMod]->m_ID[chIndex]        = iSubMod;
	  wConv->mod[iMod]->m_Pedestal[chIndex]  = fitFunc->GetParameter(4);
	  wConv->mod[iMod]->m_Signal[chIndex]    = fitFunc->GetParameter(0);
	  wConv->mod[iMod]->m_Time[chIndex]      = fitFunc->GetParameter(1);
	  wConv->mod[iMod]->m_HHTime[chIndex]    = halfTiming;
	  wConv->mod[iMod]->m_ParA[chIndex]      = fitFunc->GetParameter(3);
	  wConv->mod[iMod]->m_ParB[chIndex]      = fitFunc->GetParameter(2);
	  wConv->mod[iMod]->m_FitTime[chIndex]   = halfFitTiming;
	  wConv->mod[iMod]->m_nDigi++;	      	    	  	      
	  
	  //std::cout << iMod << ":" << iSubMod << ":" << gr->GetMean(0) << std::endl; 	      
	  //gr->Draw("AP");
	  //can->Update();
	  //can->Modified();
	  //getchar();	  

	  gr->GetListOfFunctions()->Delete();
	  delete linearFunction;
	  wavFitter->Clear();
	  
	}
      }
    }

    // Except CsIData   
    for( int iMod = 1; iMod < wConv->GetNmodule(); iMod++ ){
      if( iMod != LaserModuleID ){ continue; }
#ifdef DEBUG
      std::cout<< __LINE__ << std::endl;
      std::cout<< iMod << std::endl;
#endif      
      int nSubModule = (wConv->ModMap[iMod]).nMod;      
      for( int iSubMod = 0; iSubMod < nSubModule; iSubMod++){	

	// Channel Mapping
	// Set C(rate), F(ADC), C(hannel) event
	int iCrate = 9999;
	int iSlot  = 9999;
	int iCh    = 9999;	
	wConv->GetCFC( iMod, iSubMod, iCrate, iSlot, iCh);	
	// Ignore unmapped channel // 
	if( iCrate == 9999 || iSlot == 9999 || iCh == 9999 ) continue;       
	
	// SetData to Graph // 	
	wConv->SetGraph( iMod, iSubMod, conv, gr);
	if( iMod == LaserModuleID && iSubMod == 0 ){
	  for( int iPoint = 0; iPoint < 48; iPoint++){
	    LaserTrigWF[iPoint] = gr->GetY()[iPoint];
	  }
	}
	// Fit graph with waveform Fitter // 
	bool fit     = wavFitter->Fit( gr ); 
	int chIndex  = (wConv->mod[iMod])->m_nDigi;	 
	if( fit ){ 
	  TF1* fitFunc      = wavFitter->GetFunction();	    
	  double halfHeight = fitFunc->GetParameter(0)/2 + fitFunc->GetParameter(4);
	  double halfTiming = fitFunc->GetX( halfHeight,
					     fitFunc->GetParameter(1)-48, fitFunc->GetParameter(1));
	  TF1* linearFunction = new TF1("func","pol1",halfTiming - 12, halfTiming + 12);
	  gr->Fit( linearFunction, "Q", "", halfTiming -12, halfTiming +12 );
	  double halfFitTiming = linearFunction->GetX( halfHeight, halfTiming -12, halfTiming +12);
	  
	  /// Setting Data for write : in this program wConv->Write() had been comment out
	  wConv->mod[iMod]->m_FitHeight[chIndex] = 1;
	  wConv->mod[iMod]->m_ID[chIndex]        = iSubMod;
	  wConv->mod[iMod]->m_Pedestal[chIndex]  = fitFunc->GetParameter(4);
	  wConv->mod[iMod]->m_Signal[chIndex]    = fitFunc->GetParameter(0);
	  wConv->mod[iMod]->m_Time[chIndex]      = fitFunc->GetParameter(1);
	  wConv->mod[iMod]->m_HHTime[chIndex]    = halfTiming;
	  wConv->mod[iMod]->m_ParA[chIndex]      = fitFunc->GetParameter(3);
	  wConv->mod[iMod]->m_ParB[chIndex]      = fitFunc->GetParameter(2);
	  wConv->mod[iMod]->m_FitTime[chIndex]   = halfFitTiming;
	  wConv->mod[iMod]->m_nDigi++;	      	    	  	      
	  
	  //std::cout << iMod << ":" << iSubMod << ":" << gr->GetMean(0) << std::endl; 	      
	  //gr->Draw("AP");
	  //can->Update();
	  //can->Modified();
	  //getchar();	  

	  gr->GetListOfFunctions()->Delete();
	  delete linearFunction;
	  wavFitter->Clear();
	}else{
	  ;
	}
      }      
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // End Waveform Fitting 
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef DEBUG
    std::cout << "TriggerDicision" << std::endl;
#endif    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////    
    /// 2Trigger dicision ///
    /// You must rewrite for fit your setting //
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////    
    //std::cout<< "Trigger Processing " << std::endl;

    for( int iMod = 0; iMod < wConv->GetNmodule(); iMod++ ){      
      if( iMod == LaserModuleID ){
	if( wConv->mod[iMod]->m_nDigi >0 ){
	  if( wConv->mod[iMod]->m_Signal[0] > 100 ){
	    wConv->m_LaserTrig = 1;
	    wConv->m_TrigFlag |= 1;
	  }
	}
      }else if( iMod == CosmicModuleID ){
	//int nSubMod = (wConv->ModMap[iMod]).nMod;
	int nSubMod = wConv->mod[iMod]->m_nDigi;
	for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++ ){	    
	  CosmicSignal[CosmicArr[wConv->mod[iMod]->m_ID[iSubMod]]] = wConv->mod[iMod]->m_Signal[iSubMod];
	  CosmicTime[CosmicArr[wConv->mod[iMod]->m_ID[iSubMod]]]   = wConv->mod[iMod]->m_Time[iSubMod];
	}
	for( int iCosmic = 0; iCosmic < 5; iCosmic++){
	  if( CosmicSignal[ iCosmic    ] > COSMIC_THRESHOLD[ iCosmic     ] ||
	      CosmicSignal[ iCosmic+10 ] > COSMIC_THRESHOLD[ iCosmic +10 ]){
	    wConv->m_CosmicTrigFlagUp |= 1 << iCosmic;
	  }
	  if( CosmicSignal[ iCosmic+5  ] > COSMIC_THRESHOLD[ iCosmic +5  ] ||
	      CosmicSignal[ iCosmic+15 ] > COSMIC_THRESHOLD[ iCosmic +15 ]){
	    wConv->m_CosmicTrigFlagDn |= 1 << iCosmic;
	  }	  
	}
	if( wConv->m_CosmicTrigFlagUp && 
	    wConv->m_CosmicTrigFlagDn ){
	  wConv->m_CosmicTrig = 1; 
	  wConv->m_TrigFlag  |= 2;
	}
      }else if( iMod == CVModuleID ){
	//int nSubMod = (wConv->ModMap[iMod]).nMod;
	int nSubMod = wConv->mod[iMod]->m_nDigi;
	for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++ ){	    
	  CVSignal[wConv->mod[iMod]->m_ID[iSubMod]] = wConv->mod[iMod]->m_Signal[iSubMod];
	  CVTime[wConv->mod[iMod]->m_ID[iSubMod]]   = wConv->mod[iMod]->m_Time[iSubMod];
	}
	for( int iSubMod = 0; iSubMod < nSubMod; iSubMod++){
	  if( wConv->mod[iMod]->m_Signal[iSubMod] > 500 ){
	    wConv->m_CVTrig    = 1;
	    wConv->m_TrigFlag |= 4;
	  }
	}
      }else{
	continue;
      } 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // End trigger dicision 
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////


    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// Fill data to Template histogram ///
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
    //if( (wConv->m_TrigFlag & 1)  != 0  || (wConv->m_TrigFlag & 2) !=0 ){
    
    //std::cout<< wConv->m_TrigFlag << std::endl;

    trout->Fill();
    
  }
    
  std::cout<< "end Loop" <<std::endl;
  std::cout<< "All Convert is Done" << std::endl;
  std::cout<< "Close" << std::endl;
  trout->Write();
  tfout->Close();
  return 0;
}
