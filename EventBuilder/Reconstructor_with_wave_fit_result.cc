
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
#include "TChain.h"

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
#include "TSpline.h"
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
#include "CsI_Module.h"
#include "IDHandler.h"


#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/GammaFinder.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "user_func.h"




const int nCrate  = 11;
const int nPoints = 48; 
const Double_t COSMIC_THRESHOLD[20] = {100,100,100,100,100,
				       100,100,100,100,100,
				       100,100,100,100,100,
				       100,100,100,100,100};


int  main(int argc,char** argv)
{
  
  if( argc !=2  && argc != 3  ){
    std::cerr << "Please Input RunNumber " << std::endl;
    return -1; 
  }
  Int_t RunNumber = atoi( argv[1] );

  Int_t iDiv;
  if( argc ==3 ){
    iDiv = atoi( argv[2] );
  }

  // GetEnvironment // 
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::string COSMICFILEDIR = std::getenv("ROOTFILE_COSMIC");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;
  std::cout << COSMICFILEDIR << std::endl;

  // Setting  Classes //
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  
  E14WaveFitter* Fitter     = new E14WaveFitter();
  TApplication* app         = new TApplication("app", &argc , argv );  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Setting Template 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  std::string TemplateRootFile = Form("%s/Data/Template/TEMPLATE_COMPLETE_GRAPH_3pi0RunList_200_400.root",
				      ANALIBDIR.c_str());
  TFile* tfTemplate = new TFile( TemplateRootFile.c_str());
  TGraph*   grTemplate[2716];
  TSpline3* splTemplate[2716];
  for( int i = 0; i< 2716; i++){
    grTemplate[i] = NULL; 
    grTemplate[i] = (TGraph*)tfTemplate->Get(Form("Template_%d", i )); 
    if( grTemplate[i] != NULL ){
      splTemplate[i] = new TSpline3(Form("Spl_Template_%d",i), grTemplate[i] );
    }else{
      splTemplate[i] = NULL;
    }
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Setting Calibration Factor
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  std::string CalibrationRootFile = Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",
					 ANALIBDIR.c_str());
  EnergyConverter* EConverter = new EnergyConverter();
  EConverter->ReadCalibrationRootFile( CalibrationRootFile.c_str() );

  std::cout<< "//////////////////////////" << std::endl;
  std::cout<< "calibration File loaded."   << std::endl;
  std::cout<< "//////////////////////////" << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read TimeOffset // 
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  std::ifstream  ifsTimeOffset("testNewWORKCompileOffset.txt");
  Double_t TimeOffset[2716]={0};
  Int_t    tempID;
  Double_t tempOffset;
  Double_t tempChisq;
  while(  ifsTimeOffset >> tempID >> tempOffset >> tempChisq ){
    TimeOffset[tempID] = tempOffset;
  }

  std::cout<< "/////////////////////////////\n";
  std::cout<< "Time Offset loaded.\n";           
  std::cout<< "/////////////////////////////\n";

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read ID map // 
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  IDHandler* idHandler = new IDHandler();

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Set input / output File 
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "/////////////////////////////////////////////\n";
  std::cout << "Set I/O File" <<std::endl;  
  std::cout << "/////////////////////////////////////////////\n";

  TFile* tfConv[nCrate];
  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate] = new TFile( Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber ));
    conv[icrate] = new E14ConvReader((TTree*)tfConv[icrate]->Get("EventTree"));
  }
  TFile* tfout = new TFile(Form("%s/run_Wave_Analyzed_%d.root",WAVEFILEDIR.c_str(), RunNumber));
  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree");
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),
					    trout);
  std::cout<< "Setting Map" << std::endl;
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
  }
  const int iCsiMod    = 0; 
  const int iCC03Mod   = 1;
  const int iOEVMod    = 2;
  const int iCVMod     = 3;
  const int iCosmicMod = 4; 
  const int iLaserMod  = 5;
  const int iEtcMod    = 6;

  std::cout<< " Prepare OutputFile " << std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Set I/O Class 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  E14GNAnaDataContainer data;
  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  data.branchOfPi0List( trout );
  ClusterFinder clusterFinder;
  GammaFinder   gFinder;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////  
  int nCSIDigi = 0;
  int CSIDigiID[3000] = {0};
  double CSIDigiE[3000]={0}, CSIDigiTime[3000]={0};
  double CSICalFactor[3000]={0};
  for( int ich = 0; ich < 3000; ++ich){
    CSICalFactor[ich] = 1;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
  std::cout<< "Pi0 Calibration" << std::endl;
  // loop analysis
  int nloop = conv->GetEntries();
  Double_t TimeLimit[2]={0};
  Double_t TimeDistrib[nPoints]; // nPoints = 48 for this program 
  Double_t TimeWeightDistrib[ nPoints ];
  
  for( int ievet = 0 ; ievet < nloop; ievt++){
    wConv->InitData();
    for( int icrate  =0; icrate < nCrate; icrate++){
      conv->GetEntry( ievent );
    }
    if( ievent %100 == 0 && ievent ){ std::cout << ievent << "/" << nloop << "\n"; }
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
	  delete linearFunction;
	  wavFitter->Clear();	  
	}
      }	
      
      if( iMod == LaserModuleID ){
	if( wConv->mod[iMod]->m_nDigi > 0){
	  if( wConv->mod[iMod]->m_Signal[0] > 100 ){
	    wConv->m_LaserTrig = 1;
	    wConv->m_TrigFlag |= 1;
	    TotalTriggerFlag  |= 1;
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
	  TotalTriggerFlag   |= 2;
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
	    TotalTriggerFlag  |= 4;
	  }
	}
      }else{
	continue;
      } 
    }
    std::cout<< "End Trigger Dicision " << std::endl;







  }



  TH1D* hist = new TH1D("his","",1000,0,2000);  
  std::cout<<"start loop analysis"<<std::endl;
  std::cout<<"# of entry : "<<nloop<<std::endl;
  for( int ievt=0; ievt<nloop; ievt++ ){

    if(ievt%(nloop/10)==0&&nloop>100)
      std::cout<<ievt/(nloop/100)<<"%"<<std::endl;
    
    for( int ichannel = 0; ichannel < 3000; ichannel++){
      CSIDigiID[ichannel]  = 0; 
      CSIDigiE[ichannel]   = 0; 
      CSIDigiTime[ichannel]= 0;
    }
    trin->GetEntry( ievt );
    timePeak = wConv->m_TimePeak;
    timeSigma= wConv->m_TimeSigma;
    nCSIDigi = 0;
    for( int ich = 0; ich < wConv->mod[iCsiMod]->m_nDigi;ich++){
      int CsiChannelID = wConv->mod[iCsiMod]->m_ID[ich];
      double Energy  = EConverter->ConvertToEnergy( CsiChannelID, wConv->mod[iCsiMod]->m_Signal[ich]);
      double TimeData= wConv->mod[iCsiMod]->m_Time[ich]-TimeOffset[ich];
      
      if( Energy > 3 ){
	CSIDigiID[nCSIDigi]   = CsiChannelID;
	CSIDigiE[nCSIDigi]    = Energy;
	CSIDigiTime[nCSIDigi] = TimeData;
	nCSIDigi++;
      }
    }

    std::cout<< nCSIDigi << std::endl;
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    clist = clusterFinder.findCluster(nCSIDigi, CSIDigiID, CSIDigiE, CSIDigiTime);
    for( std::list<Cluster>::iterator it  = clist.begin();
	 it != clist.end();
	 ++it){
      std::cout<< (*it).e() << std::endl;
    }
    gFinder.findGamma(clist,glist);

    // pi0 reconstrunction and position correction for angle dependency
    if( glist.size() != 2 ) continue; 
    std::list<Pi0> piList; 
    double mass=0;

    //double position =3484.; //shiomisan:3534.;
    double position = 3526;
    if(!user_rec(glist,piList,mass,position)) continue;
    std::cout<< mass << std::endl;
    hist->Fill(mass);
    // cuts
    //user_cut(data,piList);    
    
    // fill data to TTree
    data.setData( piList );
    trout->Fill();
    data.eventID++;

  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  trout->Write();
  tfout->Close();
  
  std::cout<< "Close" << std::endl;
  return 0;

}
