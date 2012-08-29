//////////////////////////////////////////////////////////////////////
// Aug. 29. 2012
// J.W.Lee @ Osaka Univ.
// Minimum Set to read conv file.
// Need to Define
// ANALYSISLIB    : Folder of AnalysisLib ~/../Analysis/AnalysisLib
// ROOTFILE_CONV  : Folder of convfile
// ROOTFILE_SUMUP : Folder of sumupfile
// and add ${ANALYSISLIB}/lib to LD_LIBRARY_PATH for compile
//////////////////////////////////////////////////////////////////////

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
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


#include "TROOT.h"
#include "Math/WrappedFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMethods.h"
#include "Math/WrappedTF1.h"


const int nCrate = 11;
int  main(int argc,char** argv)
{
  if( argc !=2 ){
    std::cerr << "Please Input RunNumber " << std::endl;
    return -1; 
  }
  Int_t RunNumber = atoi( argv[1] );
  
  // Read Environment // 
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;


  // Open convfile // 
  TFile* tf[nCrate];
  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
  }

  // Open output file // 
  TFile* tfout = new TFile(Form("%s/run%d_wav.root",WAVEFILEDIR.c_str(),RunNumber),
			   "recreate");
  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree");   
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),
					    trout);
  tfout->cd();
  { // just for folding // 
    ////////////////////////////////////////
    // Set Modules to E14ConvWriter class //
    ////////////////////////////////////////
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
    ////////////////////////////////////////////////////////////////////////
    // Check Entries of each convfile if different each other, Send Alert // 
    ////////////////////////////////////////////////////////////////////////
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

  ////////////////////////////////
  // Set graph to dump waveform //
  ////////////////////////////////
  TGraph* gr = new TGraph();
  gr->SetMarkerStyle(6);

  TApplication* app  = new TApplication("app",&argc, argv );
  TCanvas* can = new TCanvas("can","can", 800,800);

  /////////////////
  // Start  Loop //
  /////////////////
  std::cout <<" Loop " <<std::endl;  
  for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
    // Initialize data members of E14ConvWriter class // 
    wConv->InitData();
    
    // Set Data // 
    for( int icrate = 0; icrate < nCrate; icrate++){
      conv[icrate]->GetEntry(ievent);
    } 
    
    // Confirm process // 
    if( ievent %100 == 0 && ievent ){ std::cout<< ievent << "/" << conv[0]->GetEntries() << std::endl;} 
    
    /////////////////////////////////////////////////////////////////////
    // Analysis for All registed channel //
    // Module : Each Detector ex: Csi, CC03, OEV...
    // SubMod : Channel of each Detector
    // If you want to see only Csi, then Use 
    // for( int iMod = CsiModuleNumber; iMod == CsiModuleNumber; iMod++){
    /////////////////////////////////////////////////////////////////////

    for( int iMod = 0; iMod < wConv->GetNmodule(); iMod++ ){      
      int nSubModule = (wConv->ModMap[iMod]).nMod;
      for( int iSubMod = 0; iSubMod < nSubModule; iSubMod++){		

	//////////////////////////////////////////////////////////////////////////////////
	// Set waveform of ModuleID = iMod && SubModuleID( channel ID ) = iSubMod to graph
	//////////////////////////////////////////////////////////////////////////////////
	gr->Set(0);
	if( wConv->SetGraph( iMod, iSubMod ,conv , gr ) == 0 ){ continue; }			

	///////////////////////////////////
	///////////////////////////////////
	// Add code to Analyze waveform  //  
	///////////////////////////////////
	///////////////////////////////////
	// example 
	/*
	gr->Draw("AP");
	can->Update();
	can->Modified();
	getchar();
	*/
      }
    }	
    //trout->Fill();
  }

  std::cout<< "end Loop" <<std::endl;
  //app->Run();
  std::cout<< "Close" << std::endl;
  trout->Write();
  tfout->Close();
  app->Run();
  return 0;

}
