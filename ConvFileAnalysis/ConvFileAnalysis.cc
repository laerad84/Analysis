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

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "Structs.h"
#include "Environment.h"

#include "ConvFileAnalysis/E14ConvReader.h"
#include "ConvFileAnalysis/E14IDHandler.h"
#include "ConvFileAnalysis/E14ConvWriter.h"

const int nCrate = 11;

int  main(int argc,char** argv)
{
  if( argc !=2 ){
    std::cerr << "Please Input RunNumber " << std::endl;
    return -1; 
  }
  Int_t RunNumber = atoi( argv[1] );

  // GetEnvironment // 
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  std::cout << ANALIBDIR << std::endl;
  int envRtn = GetEnvironment();
  PrintEnvironment();

  // Setting  Classes //
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  
  TFile* tf[nCrate];

  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",convFileDir, icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
    //conv[icrate]->SetBranchAddress();
    //conv[icrate]->AddFile( Form("%s/crate%d/run%d_conv.root",convFileDir, icrate, RunNumber)); 
  }  
  TFile* tfout = new TFile(Form("run%d_wav.root",RunNumber),"recreate");
  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree");   
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",sumFileDir,RunNumber),
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
  
  TApplication* app = new TApplication("app", &argc , argv );
  
  TCanvas* can = new TCanvas( "can ", "Canvas" ,800,800);
  std::cout <<"Loop " <<std::endl;
  TGraph* gr = new TGraph();
  gr->SetMarkerStyle(6);
  
  for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
    wConv->InitData();
    for( int icrate = 0; icrate < nCrate; icrate++){
      //std::cout << icrate << std::endl;
      conv[icrate]->GetEntry(ievent);
    } 
    if( ievent %100 == 0 && ievent ){ std::cout<< ievent << "/" << conv[0]->GetEntries() << std::endl;} 
    
    // Analysis
    //std::cout<< "Total Number of Module :" << wConv->GetNmodule() << std::endl;
    

    //for( int iMod = 1; iMod < 2; iMod++){
    for( int iMod = 0; iMod < wConv->GetNmodule(); iMod++ ){
      
      int nSubModule = (wConv->ModMap[iMod]).nMod;
      //std::cout<< iMod << std::endl;
      //std::cout<< iMod << " : " << nSubModule << std::endl;
      
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
	      TF1* fitFunc = wavFitter->GetFunction();	    
	      wConv->mod[iMod]->m_Fit[chIndex]      = 1;
	      wConv->mod[iMod]->m_ID[chIndex]       = iSubMod;
	      wConv->mod[iMod]->m_Pedestal[chIndex] = fitFunc->GetParameter(4);
	      wConv->mod[iMod]->m_Signal[chIndex]   = fitFunc->GetParameter(0);
	      wConv->mod[iMod]->m_Timing[chIndex]   = fitFunc->GetParameter(1);
	      wConv->mod[iMod]->m_HHTiming[chIndex] 
		= fitFunc->GetX(fitFunc->GetParameter(0)/2+fitFunc->GetParameter(4),fitFunc->GetParameter(1)-56, fitFunc->GetParameter(1));
	      wConv->mod[iMod]->m_ParA[chIndex]     = fitFunc->GetParameter(3);
	      wConv->mod[iMod]->m_ParB[chIndex]     = fitFunc->GetParameter(2);
	      wConv->mod[iMod]->m_nDigi++;	      	    
	      
	    //std::cout << iMod << ":" << iSubMod << ":" << gr->GetMean(0) << std::endl; 
	    
	    //gr->Draw("AP");
	    //can->Update();
	    //can->Modified();
	    //getchar();
	    wavFitter->Clear();	  
	  }
      }
    }	
    trout->Fill();
  }

  std::cout<< "end Loop" <<std::endl;
  //app->Run();
  std::cout<< "Close" << std::endl;
  trout->Write();
  tfout->Close();
  return 0;
}
