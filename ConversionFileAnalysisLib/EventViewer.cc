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
  if( argc !=5 ){
    std::cerr << "Please Input RunNumber " << std::endl;
    return -1; 
  }
  Int_t RunNumber   = atoi( argv[1] );
  Int_t EventNumber = atoi( argv[2] );
  Int_t ModNumber   = atoi( argv[3] );
  Int_t SubModNumber= atoi( argv[4] );

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
  E14WaveFitter* Fitter  = new E14WaveFitter();
  
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

  TFile* tfTemplete = new TFile("TEMPLETE_OUT_HEIGHT_0_OK.root");
  TGraphErrors* tempGr[2716];
  TSpline3* tempSpl[2716]; 
  for( int ispline = 0; ispline< 2716; ispline++){
    tempGr[ispline]  = (TGraphErrors*)tfTemplete->Get(Form("Waveform_Height_%d_0",ispline));    

    tempSpl[ispline] = new TSpline3(Form("waveform_Spl_%d",ispline),(TGraph*)(tempGr[ispline]));
    
  }
  //  tfTemplete->Close();

  //  TFile* tfout = new TFile(Form("%s/TEMPLETE_FIT_RESULT_%d.root",WAVEFILEDIR.c_str(),RunNumber),
  //"recreate");
  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree");   

  std::cout << Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber) << std::endl;

  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),
					    trout);
  std::cout<< "Setting Map" << std::endl;
  //  tfout->cd();
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
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<< "Setting Hist" << std::endl;
  
  TApplication* app = new TApplication("app", &argc , argv );  
  TCanvas* can      = new TCanvas( "can ", "Canvas" ,800,800);
  //  can->Divide( 2, 2);
  TGraph* gr        = new TGraph();
  gr->SetMarkerStyle(6);

  int CsiModuleID    = wConv->GetModuleID("Csi");
  int CosmicModuleID = wConv->GetModuleID("Cosmic");
  int CVModuleID     = wConv->GetModuleID("CV");
  int LaserModuleID  = wConv->GetModuleID("Laser");
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Set Trigger Map Cosmic Laser CV 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
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
  
  TText* text = new TText(0,0,"");  
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  // Loop Start  /// 
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout <<"Loop " <<std::endl;
  //for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
  int ievent = EventNumber;
  for( int icrate = 0; icrate < nCrate; icrate++){
    //std::cout << icrate << std::endl;
    conv[icrate]->GetEntry(ievent);
  } 
    
  int iMod   = ModNumber;
  int iSubMod= SubModNumber;
  std::cout<< __LINE__ << std::endl;
  if( iMod >= wConv->GetNmodule()){ return -1 ;}
  std::cout<< __LINE__ << std::endl;
  //if( iMod == wConv->GetModuleID("Csi") ){ return -1;}
  std::cout<< __LINE__ << std::endl;
  int nSubModule = wConv->GetNsubmodule( iMod );
  if( nSubModule <= 0 ){ return -1 ;}      
  std::cout<< __LINE__ << std::endl;
  std::cout<< iSubMod << " : " << nSubModule << std::endl;
  if( iSubMod >= nSubModule ){ return -1 ;}
  std::cout<< __LINE__ << std::endl;
  TGraph* grother = new TGraph();
  grother->SetMarkerStyle(7);
  if( wConv->SetGraph( iMod, iSubMod ,conv , gr ) == 0 ){ return -1; }	       
  if( wConv->SetGraph( iMod, iSubMod ,conv , grother ) == 0 ){ return -1; }	       
  int crate, fadc, channel; 
  wConv->GetCFC(iMod, iSubMod, crate, fadc ,channel );
  std::cout<< crate <<  " : " << fadc << " : " << channel << std::endl;

  TGraph*deltagr = new TGraph();
  TGraph*deltagrother= new TGraph();
  deltagr->SetMarkerStyle(7);
  deltagrother->SetMarkerStyle(7);
  
  //can->cd(1);
  can->Divide(2,2);
  can->cd(1);
  wavFitter->Fit(gr);
  gr->Draw("AP");
  TF1* func = wavFitter->GetFunction();
  for( int ipoint = 0; ipoint < gr->GetN(); ipoint++){
    deltagr->SetPoint( ipoint, ipoint*8, gr->GetY()[ipoint] - func->Eval( gr->GetX()[ipoint] ));
  }
  

  can->cd(2);
  //Fitter->SetWaveform( tempSpl[iSubMod]);
  //Fitter->Fit(grother);
  //std::cout<< Fitter->m_FitFunc->Eval( grother->GetX()[2]) << std::endl;
  Fitter->SetWaveform( tempSpl[iSubMod] );

  std::cout << E14WaveFitter::m_spl->Eval(10)   << std::endl;
  std::cout << Fitter->CheckWaveform( grother ) << std::endl; 
  std::cout << tempSpl[iSubMod]->Eval(10)       << std::endl;
  
  tempSpl[iSubMod]->Draw();
  for( int ipoint = 0; ipoint < grother->GetN(); ipoint++){
    deltagrother->SetPoint( ipoint, ipoint*8, grother->GetY()[ipoint] - Fitter->m_FitFunc->Eval( grother->GetX()[ipoint] ));
  }

  std::cout << Fitter->m_FitFunc->GetChisquare()/ Fitter->m_FitFunc->GetNDF() << std::endl;
  std::cout << func->GetChisquare()/func->GetNDF() << std::endl;
  
  grother->Draw("AP");
  can->Update();
  can->Modified();
  can->cd(3);
  deltagr->Draw("AP");
  can->cd(4);
  deltagrother->Draw("AP");


  
  /*
  getchar();
  if( iMod != CsiModuleID ){
    bool fit = wavFitter->Fit( gr ); 
    int chIndex  = (wConv->mod[iMod])->m_nDigi;	 
    if( fit ){ 
      TF1* fitFunc      = wavFitter->GetFunction();	    
      double halfHeight = fitFunc->GetParameter(0)/2 + fitFunc->GetParameter(4);
      double halfTiming = fitFunc->GetX( halfHeight,
					 fitFunc->GetParameter(1)-48, fitFunc->GetParameter(1));
      wConv->mod[iMod]->m_Fit[chIndex]      = 1;
      wConv->mod[iMod]->m_ID[chIndex]       = iSubMod;
      wConv->mod[iMod]->m_Pedestal[chIndex] = fitFunc->GetParameter(4);
      wConv->mod[iMod]->m_Signal[chIndex]   = fitFunc->GetParameter(0);
      wConv->mod[iMod]->m_Timing[chIndex]   = fitFunc->GetParameter(1);
      wConv->mod[iMod]->m_HHTiming[chIndex] = halfTiming;
      wConv->mod[iMod]->m_ParA[chIndex]     = fitFunc->GetParameter(3);
      wConv->mod[iMod]->m_ParB[chIndex]     = fitFunc->GetParameter(2);
      wConv->mod[iMod]->m_nDigi++;	      	    
      
      TF1* linearFunction = new TF1("func","pol1",halfTiming - 12, halfTiming + 12);
      gr->Fit( linearFunction, "Q", "", halfTiming -12, halfTiming +12 );
      double halfFitTiming = linearFunction->GetX( halfHeight, halfTiming -12, halfTiming +12);
      wConv->mod[iMod]->m_FitTiming[chIndex]= halfFitTiming;
      delete linearFunction;
      wavFitter->Clear();	  
    }
  }else{
    
    can->cd(2);
    tempGr[ iSubMod ]->Draw("AP");
    for( int i  =0; i< tempGr[iSubMod]->GetN() ; i++){
      std::cout<< tempGr[iSubMod]->GetX()[i] << std::endl;
    }
    TSpline3* spltest = new TSpline3("spltest",(TGraph*)tempGr[iSubMod]);
    spltest->Draw();
    tempSpl[iSubMod]->SetLineColor(2);
    //tempSpl[iSubMod]->Draw("same");
    can->Update();
    can->Modified();
    getchar();
    std::cerr << tempSpl[ iSubMod ] ->GetXmin() << " : " << tempSpl[iSubMod]->GetXmax() << std::endl;
    std::cout<< iSubMod << std::endl;
    Fitter->SetWaveform( tempSpl[iSubMod]);
    //tempSpl[iSubMod]->Eval(-145.0);
    //Fitter->SetWaveform( spltest );
    Fitter->InitPar();
    std::cout<< "Fit" << std::endl;
    bool fit = Fitter->Fit(gr, 0, 400);
    int chIndex = (wConv->mod[iMod])->m_nDigi;
    std::cout<< RunNumber << " : " << ievent << " : " << iSubMod << std::endl ;
    if( fit ){ 
      can->cd();
      gr->SetNameTitle(Form("gr_%d_%d",iMod,iSubMod),Form("gr_%d_%d",iMod,iSubMod));
      gr->Draw("AP");
      
      text->DrawTextNDC(0.2,0.8,Form("CHISQ/NDF:%lf",Fitter->GetFitResult()));
      Fitter->m_FitFunc->Draw("same");
      can->Update();
      can->Modified();
      
      std::cout<< Fitter->GetFitResult() << std::endl;
      
      wConv->mod[iMod]->m_Fit[chIndex]      = 1;
      wConv->mod[iMod]->m_ID[chIndex]       = iSubMod;
      wConv->mod[iMod]->m_Pedestal[chIndex] = Fitter->GetParameter(2);
      wConv->mod[iMod]->m_Signal[chIndex]   = Fitter->GetParameter(0);
      wConv->mod[iMod]->m_Timing[chIndex]   = Fitter->GetParameter(1);
      wConv->mod[iMod]->m_HHTiming[chIndex] = Fitter->GetConstantFraction();
      wConv->mod[iMod]->m_nDigi++;
    }
      
    Fitter->Clear();
  }
  */
  std::cout<< "end Loop" <<std::endl;
  app->Run();
  std::cout<< "Close" << std::endl;
  //trout->Write();
  //tfout->Close();
  return 0;
}
