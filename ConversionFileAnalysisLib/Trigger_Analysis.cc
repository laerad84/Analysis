
// The Propose of this program is making Templetes of Trigger Signal
// For example, Cosmic, Laser and CV


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


const int nCrate = 11;

const Double_t COSMIC_THRESHOLD[20] = {100,100,100,100,100,
				       100,100,100,100,100,
				       100,100,100,100,100,
				       100,100,100,100,100};


int  main(int argc,char** argv)
{
  if( argc !=2 ){
    std::cerr << "Please Input RunNumber " << std::endl;
    return -1; 
  }
  Int_t RunNumber = atoi( argv[1] );

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // GetEnvironment // 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////  
  // Setting  Classes //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  
  E14WaveFitter* Fitter  = new E14WaveFitter();
  //TApplication* app = new TApplication("app", &argc , argv );  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Setting IO File
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Setting IO File" << std::endl;
  TFile* tf[nCrate];
  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    std::cout << Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber) << std::endl;
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
  }

  //TFile* tfout = new TFile(Form("run%d_wav.root",RunNumber),"recreate");
  std::cout << Form("%s/TEMPLETE_COSMIC_%d.root",WAVEFILEDIR.c_str(),RunNumber) << std::endl;


  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Set Template  /// 
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  
  TFile* tfTemplate = new TFile("TEMPLETE_OUT_HEIGHT_0_OK.root");
  TGraphErrors* tempGr[2716];
  TSpline3*     tempSpl[2716]; 
  for( int ispline = 0; ispline< 2716; ispline++){
    tempGr[ispline]  = NULL;
    tempSpl[ispline] = NULL;
    tempGr[ispline]  = (TGraphErrors*)tfTemplate->Get(Form("Waveform_Height_%d_0",ispline));
    if( tempGr[ispline]->GetN()!= 0){
      tempSpl[ispline] = new TSpline3(Form("waveformSpl_%d",ispline),(TGraph*)(tempGr[ispline]));
    }else{
      std::cout<< "Non Exist channel:" << ispline << std::endl;
    }
  }
  tfTemplate->Close();

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Set input / output File 
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  //  TCanvas* test = new TCanvas("test","",400,400);
  TFile* tfout  = new TFile(Form("%s/Trigger_Analysis_%d.root",WAVEFILEDIR.c_str(),RunNumber),
			    "recreate");
  TTree* trout  = new TTree("WFTree","Waveform Analyzed Tree");   
  std::cout << Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber) << std::endl;  
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),
					    trout);
  std::cout<< "Setting Map" << std::endl;
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

  std::cout << "Setting IO File End" << std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout<< "Setting Hist" << std::endl;
  
  TCanvas* can      = new TCanvas( "can ", "Canvas" ,800,800);
  TGraph* gr        = new TGraph();
  gr->SetMarkerStyle(6);

  int CsiModuleID    = wConv->GetModuleID("Csi");
  int CosmicModuleID = wConv->GetModuleID("Cosmic");
  int CVModuleID     = wConv->GetModuleID("CV");
  int LaserModuleID  = wConv->GetModuleID("Laser");

  ///////////////////////////////////////////////////////////////////////////////////////////////////////  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set Trigger Map Cosmic Laser CV 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout<< "Setting Trigger" << std::endl;
  const int nCVModule     = 10;
  const int nCosmicModule = 20; 
  const int CosmicArr[20] = {4 ,5 ,2 ,3 ,6 ,7 ,0 ,1 ,12,13,10,11,14,15,8 ,9 ,16,17,18,19};
  
  if((wConv->ModMap[CVModuleID]).nMod     != 10)
    { 
      std::cout<< "CV nModule is not equal"     << std::endl;
    }
  if((wConv->ModMap[CosmicModuleID]).nMod != 20)
    { 
      std::cout<< "Cosmic nModule is not equal" << std::endl;
    }
  
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
  TH2D* hisCosmic[20][2];
  TH2D* hisLaser[2];
  TH2D* hisCV[10][2];
  hisLaser[0] = new TH2D("hisLaser_NoHit","hisLaser_NoHit",48,0,48*8,320,0,16000);
  hisLaser[1] = new TH2D("hisLaser_Hit"  ,"hisLaser_Hit"  ,48,0,48*8,320,0,16000);
  for( int i = 0; i< 20; i++){
    hisCosmic[i][0] = new TH2D(Form("hisCosmicNoHit_%d",i),Form("hisCosmicNoHit_%d",i),48,0,48*8,320,0,16000);
    hisCosmic[i][1] = new TH2D(Form("hisCosmicHit_%d",i)  ,Form("hisCosmicHit_%d",i)  ,48,0,48*8,320,0,16000);
  }
  for( int i = 0; i< 10; i++){
    hisCV[i][0] = new TH2D(Form("hisCVNoHit_%d",i),Form("hisCVNoHit_%d"), 48, 0,48*8, 320,0,16000);
    hisCV[i][1] = new TH2D(Form("hisCVHit_%d",i)  ,Form("hisCVHit_%d")  , 48, 0,48*8, 320,0,16000);
  }
  
  TH1D* hisTrig_Raw = new TH1D("hisTrig_Raw","hisTrig_Raw",8,0,8);
  char* binIndex[8] = {"Laser",
		       "Cosmic",
		       "CV",
		       "Laser&Cosmic",
		       "Laser&CV",
		       "Cosmic&CV",
		       "Laser&Cosmic&CV",
		       "Other Triggers"};
  for( int ibin  = 1; ibin < 8+1; ibin++){
    hisTrig_Raw->GetXaxis()->SetBinLabel(ibin,binIndex[ibin-1]);
  }

  //TText* text = new TText(0,0,"");  

  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop Start  /// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout <<"Loop " <<std::endl;
  for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
  //for( int ievent  = 0; ievent < 2000 ; ievent++){
  //for( int ievent  = 0; ievent < 5 ; ievent++){
    //std::cout<< ievent << std::endl;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// Init Data Component;
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
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
    
    //// Init Data Component; 

    wConv->m_RunNo   = RunNumber;
    wConv->m_EventNo = ievent;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Trigger dicision  && Other detector...
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    for( int iMod = 0; iMod < wConv->GetNmodule(); iMod++ ){      
      //std::cout << iMod << std::endl;
      if( iMod == wConv->GetModuleID("Csi") ){ continue; }
      int nSubModule = wConv->GetNsubmodule( iMod );
      if( nSubModule <= 0 ){ continue ;}      
      for( int iSubMod = 0; iSubMod < nSubModule; iSubMod++){	
	if( wConv->SetGraph( iMod, iSubMod ,conv , gr ) == 0 ){ continue; }	       
	
	bool fit = wavFitter->Fit( gr ); 
	int chIndex  = (wConv->mod[iMod])->m_nDigi;	 
	if( fit ){ 
	  TF1* fitFunc      = wavFitter->GetFunction();	    
	  double halfHeight = fitFunc->GetParameter(0)/2 + fitFunc->GetParameter(4);
	  double halfTiming = fitFunc->GetX( halfHeight,
					     fitFunc->GetParameter(1)-48, fitFunc->GetParameter(1));
	  wConv->mod[iMod]->m_FitHeight[chIndex]= 1;
	  wConv->mod[iMod]->m_ID[chIndex]       = iSubMod;
	  wConv->mod[iMod]->m_Pedestal[chIndex] = fitFunc->GetParameter(4);
	  wConv->mod[iMod]->m_Signal[chIndex]   = fitFunc->GetParameter(0);
	  wConv->mod[iMod]->m_Time[chIndex]     = fitFunc->GetParameter(1);
	  wConv->mod[iMod]->m_HHTime[chIndex]   = halfTiming;
	  wConv->mod[iMod]->m_ParA[chIndex]     = fitFunc->GetParameter(3);
	  wConv->mod[iMod]->m_ParB[chIndex]     = fitFunc->GetParameter(2);
	  wConv->mod[iMod]->m_nDigi++;	      	    
	  
	  TF1* linearFunction = new TF1("func","pol1",halfTiming - 12, halfTiming + 12);
	  gr->Fit( linearFunction, "Q", "", halfTiming -12, halfTiming +12 );
	  double halfFitTiming = linearFunction->GetX( halfHeight, halfTiming -12, halfTiming +12);
	  wConv->mod[iMod]->m_FitTime[chIndex]= halfFitTiming;
	  //fitFunc->Clear();
	  delete linearFunction;
	  wavFitter->Clear();	  
	}
      }	
      
      if( iMod == LaserModuleID ){
	if( wConv->mod[iMod]->m_nDigi > 0){
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
	//std::cout<< __LINE__ << std::endl;
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
	//std::cout<< __LINE__ << std::endl;
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

    if( wConv->m_CosmicTrig ){
      for( int iSubMod = 0; iSubMod < wConv->mod[CosmicModuleID]->m_nDigi; iSubMod++){
	int crate;
	int FADC;
	int channel;
	wConv->GetCFC( CosmicModuleID , wConv->mod[CosmicModuleID]->m_ID[iSubMod], crate, FADC, channel );	
	for( int ipoint  = 0; ipoint < 48; ipoint++){
	  hisCosmic[wConv->mod[CosmicModuleID]->m_ID[iSubMod]][1]->Fill( ipoint*8,conv[crate]->Data[ FADC ][ channel ][ ipoint ]);
	}
      }
    }else{
      for( int iSubMod = 0; iSubMod < (wConv->ModMap[CosmicModuleID]).nMod; iSubMod++){
	int crate;
	int FADC;
	int channel;
	wConv->GetCFC( CosmicModuleID, iSubMod, crate, FADC, channel );
	for( int ipoint  = 0;ipoint < 48; ipoint++){
	  hisCosmic[iSubMod][0]->Fill( ipoint* 8 , conv[crate]->Data[ FADC ][ channel ][ ipoint ] );
	}
      }
    }

    if( wConv->m_LaserTrig ){
      int crate;
      int FADC;
      int channel; 
      wConv->GetCFC( LaserModuleID, 0, crate, FADC , channel );
      for( int ipoint = 0; ipoint < 48; ipoint++ ){
	hisLaser[1]->Fill( ipoint*8, conv[crate]->Data[ FADC ][ channel ][ ipoint ] );
      }
    }else{
      int crate;
      int FADC; 
      int channel;
      wConv->GetCFC( LaserModuleID, 0 , crate, FADC, channel );
      for( int ipoint = 0; ipoint < 48; ipoint++){
	hisLaser[0]->Fill( ipoint*8, conv[crate]->Data[ FADC ][ channel ][ ipoint ] );
      }
    }
    if( wConv->m_LaserTrig ){
      hisTrig_Raw->Fill(0);
      if( wConv->m_CosmicTrig ){
	hisTrig_Raw->Fill(3);
	if( wConv->m_CVTrig ){
	  hisTrig_Raw->Fill(6);
	}
      }else if( wConv->m_CVTrig ){
	hisTrig_Raw->Fill(4);
      }
    }
    
    if( wConv->m_CosmicTrig ){
      hisTrig_Raw->Fill(1);
      if( wConv->m_CVTrig ){
	hisTrig_Raw->Fill(5);
      }
    }
    
    if( wConv->m_CVTrig){
      hisTrig_Raw->Fill(2);
    }
    
    if( !(wConv->m_CVTrig)     &&
	!(wConv->m_CosmicTrig) &&
	!(wConv->m_LaserTrig)  ){    
      hisTrig_Raw->Fill(7);
    }
    
  }
  
  for( int i = 0; i< 20; i++){
    hisCosmic[i][0]->Write();
    hisCosmic[i][1]->Write();
  }
  

  hisLaser[0]->Write();
  hisLaser[1]->Write();

  hisTrig_Raw->Write();

  trout->Write();
  tfout->Close();
}
