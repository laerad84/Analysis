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
#include "TSpline.h"

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


TSpline3* splCsi[2716];
TSpline3* splLaser;
TSpline3* splCsiTemp[2716];
TSpline3* splLaserTemp;
TSpline3* spl;

double FuncSplCsi( double *x, double *p){
  double value =0;
  if( x[0] > -150 && x[0] < 250 ){
    value = spl->Eval(x[0]-p[1])*p[0]+p[2];
  }else{
    value = 0;
  }
  return value;
}
double FuncSplLaser( double* x, double *p ){
  double value = 0;
  if( x[0] -150 && x[0] < 250 ){
    value = splLaser->Eval(x[0]-p[1])*p[0]+p[2];
  }else{
    value = 0;
  }
  return 0; 
}




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
  std::string ANALYSISLIB = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string ROOTFILE_CONV= std::getenv("ROOTFILE_CONV");
  std::string ROOTFILE_WAV= std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::string SUMFILEDIR_1= "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/2012_FEB/Sumup";
  
  std::cout << ANALYSISLIB << std::endl;  
  std::cout << ROOTFILE_CONV << std::endl;
  std::cout << ROOTFILE_WAV << std::endl;
  std::cout << SUMFILEDIR  << std::endl;

  std::string iFileForm = "%s/crate%d/run%d_conv.root";
  std::string mapFileDir;
  if( stat(Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),&st) == 0 ){
    mapFileDir = SUMFILEDIR;
  }else{
    mapFileDir = SUMFILEDIR_1;
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  /////////// Read Wave form && Set Spline //////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////


  TH2D*     tempLaser;
  TProfile* profLaser;
  TGraph*   grLaser;
  TGraph*   grLaserTemp;

  TH2D*     tempCsi[2716]={NULL};
  TProfile* profCsi[2716]={NULL};
  TGraph*   grCsi[2716]={NULL};
  TGraph*   grCsiTemp[2716]={NULL};
  

  std::string TemplateFile = "%s/Data/laserTemplate.root";
  std::string LaserTemplateFile = "%s/Data/LaserWaveTemplate.root";

  // Set STemplate to hist // 
  TFile* tfLaserTemplate = new TFile(Form(LaserTemplateFile.c_str(),ANALYSISLIB.c_str()));
  tempLaser = (TH2D*)tfLaserTemplate->Get("hisLaserTemplate");
  profLaser = tempLaser->ProfileX();
  profLaser->SetLineColor(1);
  grLaser     = new TGraph();
  grLaserTemp = new TGraph();
  for( int iBin = 0; iBin < profLaser->GetNbinsX();iBin++){
    grLaserTemp->SetPoint( grLaserTemp->GetN(),profLaser->GetBinCenter(iBin+1),profLaser->GetBinContent(iBin+1));
  }
  splLaserTemp = new TSpline3("splLaserTemp",grLaserTemp);
  for( int iBin = 0; iBin < profLaser->GetNbinsX();iBin++){
    grLaser->SetPoint( grLaser->GetN(),profLaser->GetBinCenter(iBin+1),profLaser->GetBinContent(iBin+1)/splLaserTemp->Eval(0));
  }
  splLaser = new TSpline3("splLaser",grLaser);
  tfLaserTemplate->Close();
  
  TFile* tfTemplate = new TFile(Form(TemplateFile.c_str(),ANALYSISLIB.c_str()));
  for( int i = 0; i< 2716;i++){
    tempCsi[i]=(TH2D*)tfTemplate->Get(Form("ChannelTemplate_%d",i));
    if( tempCsi[i]->GetEntries() < 4800 ){ continue; }
    profCsi[i]=tempCsi[i]->ProfileX();
    grCsi[i] = new TGraph();
    grCsiTemp[i] = new TGraph();
    for( int iBin = 0; iBin < profCsi[i]->GetNbinsX();iBin++){
      grCsiTemp[i]->SetPoint( grCsiTemp[i]->GetN(),profCsi[i]->GetBinCenter(iBin+1),profCsi[i]->GetBinContent(iBin+1));
    }
    splCsiTemp[i] = new TSpline3(Form("splCsiTemp_%d",i),grCsiTemp[i]);
     for( int iBin = 0; iBin < profCsi[i]->GetNbinsX();iBin++){
       grCsi[i]->SetPoint( grCsi[i]->GetN(),profCsi[i]->GetBinCenter(iBin+1),profCsi[i]->GetBinContent(iBin+1)/splCsiTemp[i]->Eval(0));
     }
     splCsi[i] = new TSpline3(Form("splCsi_%d",i),grCsi[i]);
  }
  tfTemplate->Close();
  
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

  TFile* tfout = new TFile(Form("%s/LaserWaveform_%d.root",ROOTFILE_WAV.c_str(),RunNumber),"recreate");
  TTree* trWave = new TTree("trWave","Selected Waveform");
  Int_t    EventNumber;
  UShort_t LaserTrigWF[48];
  Int_t    CsiNumber;
  Int_t    CsiID[2716];
  UShort_t CsiOutWF[2716][48];
  Double_t LaserTime;
  Double_t LaserHeight;
  Double_t LaserPedestal;
  Double_t CsiTime[2716];
  Double_t CsiHeight[2716];
  Double_t CsiPedestal[2716];
  
  trWave->Branch("RunNumber",&RunNumber,"RunNumber/I");
  trWave->Branch("EventNumber",&EventNumber,"EventNumber/I");

  //trWave->Branch("LaserTrigWF",LaserTrigWF,"LaserTrigWF[48]/s");
  trWave->Branch("LaserHeight",&LaserHeight,"LaserHeight/D");
  trWave->Branch("LaserTime",&LaserTime,"LaserTime/D");
  trWave->Branch("LaserPedestal",&LaserPedestal,"LaserPedestal/D");
  trWave->Branch("CsiNumber",&CsiNumber,"CsiNumber/I");
  trWave->Branch("CsiID",CsiID,"CsiID[CsiNumber]/I");     //CsiNumber
  //trWave->Branch("CsiOutWF",CsiOutWF,"CsiOutWF[CsiNumber][48]/s"); //CsiNumber
  trWave->Branch("CsiTime",CsiTime,"CsiTime[CsiNumber]/D");        //CsiNumber
  trWave->Branch("CsiHeight",CsiHeight,"CsiHeight[CsiNumber]/D");  //CsiNumber
  trWave->Branch("CsiPedestal",CsiPedestal,"CsiPedestal[CsiNumber]/D");//CsiNumber
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

  TF1* funcLaser = new TF1("funcLaser",FuncSplLaser, -100,200,3);
  TF1* funcCsi   = new TF1("funcCsi"  ,FuncSplCsi  , -100,200,3);

  std::cout <<"Loop " <<std::endl;
  for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
  //for( int ievent  = 0; ievent < 100; ievent++){
    EventNumber = ievent;
    wConv->InitData();
    CsiNumber = 0;
    for( int i = 0; i< 2716; i++){
      CsiID[i] = -1;
      for( int j = 0; j< 48; j++){
	CsiOutWF[i][j] = 0;
      }
    }

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
    
    //// Check LaserWaveform //// 
    {
      int iMod       = LaserModuleID;
      int nSubModule = (wConv->ModMap[iMod]).nMod;      
      int iSubMod    = 0;
      // Csi Mapping
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
	  LaserTrigWF[iPoint] = (UShort_t)(gr->GetY()[iPoint]);
	}
      }
      // Fit graph with waveform Fitter // 
      bool fit  = wavFitter->Fit(gr);
      if( fit ){
	TF1* fitFunc = wavFitter->GetFunction();
	Double_t tmpHeight = fitFunc->GetParameter(0);
	Double_t tmpTime   = fitFunc->GetParameter(1);
	Double_t tmpPed    = fitFunc->GetParameter(2);

	funcLaser->SetParameter(0,tmpHeight);
	funcLaser->SetParameter(1,tmpTime);
	funcLaser->SetParameter(2,tmpPed);
	funcLaser->SetParLimits(0,tmpHeight*0.5,tmpHeight*1.5);
	funcLaser->SetParLimits(1,tmpTime-24,tmpTime+24);
	funcLaser->SetParLimits(2,tmpPed-2,tmpPed+2);
	gr->Fit(funcLaser, "Q","",tmpTime-100,tmpTime+50);
	LaserHeight = funcLaser->GetParameter(0);
	LaserTime   = funcLaser->GetParameter(1);
	LaserPedestal= funcLaser->GetParameter(2);
	gr->GetListOfFunctions()->Delete();
      }else{
	LaserHeight = 0;
	LaserTime   = 0; 
      }
    }

    if( LaserHeight < 100 ){ continue; }
    if( LaserTime   < 50 || LaserTime > 300 ){ continue; }

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
	if( splCsi[iSubMod] == NULL ){continue; }
	spl = splCsi[iSubMod];
	// SetData to Graph // 	
	wConv->SetGraph( iMod, iSubMod, conv, gr);
	for( int iPoint = 0; iPoint < 48 ; iPoint++){
	  CsiOutWF[CsiNumber][iPoint] = (UShort_t)(gr->GetY()[iPoint]);
	}
	CsiID[CsiNumber] = iSubMod;
	/// further Analysis ///
	bool fit  = wavFitter->Fit(gr);
	if( fit ){
	  TF1* fitFunc = wavFitter->GetFunction();
	  Double_t tmpHeight = fitFunc->GetParameter(0);
	  Double_t tmpTime   = fitFunc->GetParameter(1);
	  Double_t tmpPed    = fitFunc->GetParameter(2);
	  
	  funcCsi->SetParameter(0,tmpHeight);
	  funcCsi->SetParameter(1,tmpTime);
	  funcCsi->SetParameter(2,tmpPed);
	  funcCsi->SetParLimits(0,tmpHeight*0.5,tmpHeight*1.5);
	  funcCsi->SetParLimits(1,tmpTime-24,tmpTime+24);
	  funcCsi->SetParLimits(2,tmpPed-2,tmpPed+2);
	  
	  gr->Fit(funcCsi,"Q","",tmpTime-100,tmpTime+50);
	  CsiHeight[CsiNumber] = funcCsi->GetParameter(0);
	  CsiTime[CsiNumber]   = funcCsi->GetParameter(1);
	  CsiPedestal[CsiNumber]=funcCsi->GetParameter(2);
	  CsiNumber++;
	  gr->GetListOfFunctions()->Delete();
	}
      }
    }
    trWave->Fill();
  }
  
  std::cout<< "end Loop" <<std::endl;
  std::cout<< "All Convert is Done" << std::endl;
  std::cout<< "Close" << std::endl;
  //trout->Write();
  trWave->Write();
  tfout->Close();
  return 0;
}
