// The Propose of this program is making Templetes of Trigger Signal
// For example, Cosmic, Laser and CV
// Make CsI templete by 2%
//
//

//#include "gnana/DigiReader.h"
//#include "gnana/E14GNAnaDataContainer.h"
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

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "Structs.h"

#include "E14ConvReader.h"
#include "E14IDHandler.h"
#include "E14ConvWriter.h"

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
  if( argc !=3 ){
    std::cerr << "Please Input RunNumber and section Number" << std::endl;
    std::cerr << argv[0] << " [RunNumber] [MinimumHeight] [MaximumHeight] " << std::endl; 
    return -1; 
  }

  Int_t RunNumber = atoi( argv[1] );
  Int_t MinHeight = atoi( argv[2] );
  Int_t MaxHeight = atoi( argv[3] );
  if( MinHeight >= MaxHeight ){
    std::cerr << "MinHeight >= MaxHeight" << std::endl;
    return -1;
  }

  // GetEnvironment // 
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;

  // Setting  Classes //
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Setting IO 
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TFile* tf[nCrate];

  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
  }
  //TFile* tfout = new TFile(Form("run%d_wav.root",RunNumber),"recreate");
  TFile* tfout = new TFile(Form("%s/TEMPLETE_HEIGHT_%d_%d_%d.root",
				WAVEFILEDIR.c_str(),RunNumber,MinHeight,MaxHeight),
			   "recreate");
  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree");   
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),
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
  Double_t HeightArr[91];
  HeightArr[0] = 1.;
  HeightArr[1] = 10.;
  HeightArr[2] = 50.;
  for( int i = 3; i < 91; i++){    
    HeightArr[i] = 50*TMath::Power(2,(i-2)*0.1);
  }

  TH2D* hisTempCsI_Height[2716];
  TH1D* hisTempCsI_Spectrum[2716];
  for( int i = 0; i< 2716; i++){
    hisTempCsI_Height[i] = new TH2D(Form("hisTemplateCsI_Height_%d_%d_%d",MinHeight,MaxHeight,i),
				    Form("hisTemplateCsI_Height_%d_%d_%d",MinHeight,MaxHeight,i),
				    400,-200, 200,200, -0.5, 1.5);
    hisTempCsI_Spectrum[i] = new TH1D(Form("hisTemplateSpectrum_%d_%d_%d",MinHeight,MaxHeight,i),
				      Form("hisTemplateSpectrum_%d_%d_%d",MinHeight,MaxHeight,i),
				      90,HeightArr);
  }  

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop Start  /// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout <<"Loop " <<std::endl;
  for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
  //for( int ievent  = 0; ievent < 10000; ievent++){
    //std::cout<< ievent << std::endl;

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
      //std::cout << icrate << std::endl;
      conv[icrate]->GetEntry(ievent);
    } 

    wConv->m_RunNo   = RunNumber;
    wConv->m_EventNo = ievent;

    if( ievent %100 == 0 && ievent ){ std::cout<< ievent << "/" << conv[0]->GetEntries() << std::endl;} 

    // Conversion Convfile to Wav File 
    //for( int iMod = 1; iMod < 2; iMod++){
    //std::cout<< "Event Processing " << std::endl ;
    for( int iMod = 0; iMod < wConv->GetNmodule(); iMod++ ){      
      int nSubModule = (wConv->ModMap[iMod]).nMod;      
      for( int iSubMod = 0; iSubMod < nSubModule; iSubMod++){	
	int iCrate = 9999;
	int iSlot  = 9999;
	int iCh    = 9999;
	iCrate = (wConv->ModMap[iMod]).Map[iSubMod][0];
	iSlot  = (wConv->ModMap[iMod]).Map[iSubMod][1];
	iCh    = (wConv->ModMap[iMod]).Map[iSubMod][2];
	  
	// Ignore unmapped channel // 
	if( iCrate == 9999 || iSlot == 9999 || iCh == 9999 ) continue;       
	  
	gr->Set(0);
	for( int ipoint = 0; ipoint< 48; ipoint++){
	  if(conv[iCrate]->Data[iSlot][iCh][ipoint]>16000){ continue; }
	  gr->SetPoint( gr->GetN(), ipoint*8, conv[iCrate]->Data[iSlot][iCh][ipoint]);
	}
	bool fit = wavFitter->Fit( gr ); 
	int chIndex  = (wConv->mod[iMod])->m_nDigi;	 
	if( fit ){ 
	  TF1* fitFunc      = wavFitter->GetFunction();	    
	  double halfHeight = fitFunc->GetParameter(0)/2 + fitFunc->GetParameter(4);
	  double halfTiming = fitFunc->GetX( halfHeight,
					     fitFunc->GetParameter(1)-48, fitFunc->GetParameter(1));
	  wConv->mod[iMod]->m_FitHeight[chIndex]      = 1;
	  wConv->mod[iMod]->m_ID[chIndex]       = iSubMod;
	  wConv->mod[iMod]->m_Pedestal[chIndex] = fitFunc->GetParameter(4);
	  wConv->mod[iMod]->m_Signal[chIndex]   = fitFunc->GetParameter(0);
	  wConv->mod[iMod]->m_Time[chIndex]   = fitFunc->GetParameter(1);
	  wConv->mod[iMod]->m_HHTime[chIndex] = halfTiming;
	  wConv->mod[iMod]->m_ParA[chIndex]     = fitFunc->GetParameter(3);
	  wConv->mod[iMod]->m_ParB[chIndex]     = fitFunc->GetParameter(2);
	  wConv->mod[iMod]->m_nDigi++;	      	    
	      
	  TF1* linearFunction = new TF1("func","pol1",halfTiming - 12, halfTiming + 12);
	  gr->Fit( linearFunction, "Q", "", halfTiming -12, halfTiming +12 );
	  double halfFitTiming = linearFunction->GetX( halfHeight, halfTiming -12, halfTiming +12);
	  wConv->mod[iMod]->m_FitTime[chIndex]= halfFitTiming;
	  //delete linearFunction;
	  //std::cout << iMod << ":" << iSubMod << ":" << gr->GetMean(0) << std::endl; 	      
	  //gr->Draw("AP");
	  //can->Update();
	  //can->Modified();
	  //getchar();
	  gr->GetListOfFunctions()->Delete();
	  wavFitter->Clear();
	}
      }      
    }
    
    /// Trigger dicision ///
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

    if( wConv->m_TrigFlag  == 1  || wConv->m_TrigFlag ==2 ){
      continue;
    }else{ //Other Beam Event// 
      for( int idigi = 0; idigi < wConv->mod[CsiModuleID]->m_nDigi; idigi++ ){
	if( wConv->mod[CsiModuleID]->m_Signal[idigi]  <  MinHeight ||
	    wConv->mod[CsiModuleID]->m_Signal[idigi]  >= MaxHeight ){
	  continue;
	}
	if( wConv->mod[CsiModuleID]->m_Chisq[idigi]/wConv->mod[CsiModuleID]->m_NDF[idigi] > 250 ){
	  // 250 is just set value ...
	  continue;
	}
	
	int iSubMod = wConv->mod[CsiModuleID]->m_ID[idigi]; 
	int CompensateID;
	if( iSubMod  < 2240 ){
	  CompensateID = 0;
	}else{ 
	  CompensateID = 2;
	}
	int iCrate = 9999;
	int iSlot  = 9999;
	int iCh    = 9999;
	wConv->GetCFC( CsiModuleID, iSubMod, iCrate, iSlot, iCh );

	/*
	iCrate = (wConv->ModMap[CsiModuleID]).Map[iSubMod][0];
	iSlot  = (wConv->ModMap[CsiModuleID]).Map[iSubMod][1];
	iCh    = (wConv->ModMap[CsiModuleID]).Map[iSubMod][2];	
	*/

	hisTempCsI_Spectrum[iSubMod]->Fill( wConv->mod[CsiModuleID]->m_Signal[idigi] );
	for( int ipoint = 0; ipoint < 48; ipoint++){
	  hisTempCsI_Height[iSubMod]->Fill( ipoint*8 - wConv->mod[CsiModuleID]->m_Time[idigi],
					    (conv[iCrate]->Data[iSlot][iCh][ipoint]- wConv->mod[CsiModuleID]->m_Pedestal[idigi])/wConv->mod[CsiModuleID]->m_Signal[idigi]);
	}
      }
    }
    
    /// All Convert is done ///
    /// Trigger Setting     ///
    ////////////////////////////////

    //trout->Fill();

  }
    
  for( int ch = 0; ch < 2716; ch++){
    hisTempCsI_Height[ch]->Write();    
    hisTempCsI_Spectrum[ch]->Write();
  }
  std::cout<< "end Loop" <<std::endl;
  std::cout<< "All Convert is Done" << std::endl;
  std::cout<< "Close" << std::endl;
  //trout->Write();
  tfout->Close();
  return 0;
}
