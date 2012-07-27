#include <iostream>
#include <fstream>
#include <string>

#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TGraph.h"

#include "TPDF.h"
#include "TCanvas.h"

#include "WaveformFitter.h"
#include "E14ConvReader.h"

int
main( int argc , char** argv){
  // ./bin/MakeCosmicTemplete [ runNumber ] [ CosmicListFile ] [outputFile]
  if( argc != 4){
    std::cerr << "Wrong number of Arguments" << std::endl;
    std::cerr << "./bin/MakeCosmicTemplete [ runNumber ] [ CosmicResultFile ] [outputFile]" << std::endl;
    return -1;
  }

  gStyle->SetOptFit ( 11111111  );
  gStyle->SetOptStat( "neMRiuo" );  
  gStyle->SetPalette(1);

  const int CosmicArr[20]= {4 ,5 ,2 ,3 ,6 ,7 ,0 ,1 ,12,13,
			    10,11,14,15,8 ,9 ,16,17,18,19};

  int RunNumber = atoi(argv[1]);
  std::string CosmicResultFile = argv[2];
  std::string OutputFile       = argv[3];
  
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"   );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV" );
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV"  );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  
  TFile* tf[nCrate];
  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",
				  CONVFILEDIR.c_str(), icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
  }
  TFile* tfout = new TFile(Form("%s/run%d_wav.root",WAVEFILEDIR.c_str(),RunNumber),
			   "recreate");
  TTree* trout = new TTree("WFTree","Waveform Analyzed Tree");   
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),trout);
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

  TH2D* hisTemp_Csi[(wConv->ModMap[0]).nMod];
  TH2D* hisTemp_CC03[(wConv->ModMap[1]).nMod];
  TH2D* hisTemp_OEV[(wConv->ModMap[2]).nMod];
  TH2D* hisTemp_CV[(wConv->ModMap[3]).nMod];
  TH2D* hisTemp_Cosmic[(wConv->ModMap[4]).nMod];
  TH2D* hisTemp_Laser[(wConv->ModMap[5]).nMod];
  for( int iSubMod = 0; iSubMod < (wConv->ModMap[0]).nMod; iSubMod++){
    hisTemp_Csi[iSubMod] = new TH2D(Form("Temp_Csi_%d",iSubMod),"",300,-100,200,150,0,1.5);
  }
  for( int iSubMod = 0; iSubMod < (wConv->ModMap[1]).nMod; iSubMod++){
    hisTemp_CC03[iSubMod] = new TH2D(Form("Temp_CC03_%d",iSubMod),"",300,-100,200,150,0,1.5);
  }
  for( int iSubMod = 0; iSubMod < (wConv->ModMap[2]).nMod; iSubMod++){
    hisTemp_OEV[iSubMod] = new TH2D(Form("Temp_OEV_%d",iSubMod),"",300,-100,200,150,0,1.5);
  }
  for( int iSubMod = 0; iSubMod < (wConv->ModMap[3]).nMod; iSubMod++){
    hisTemp_CV[iSubMod] = new TH2D(Form("Temp_CV_%d",iSubMod),"",300,-100,200,150,0,1.5);
  }

  TGraph* gr = new TGraph();
  gr->SetMarkerStyle( 6 );
  for( int  ievent  = 0; ievent < conv[0]->GetEntries(); ievent++ ){
    wConv->InitData();
    for( int icrate = 0; icrate < nCrate; icrate++){
      conv[icrate]->GetEntry(ievent);
    }
    if( ievent % 100 = 0 && ievent ){ 
      std::cout<< ievent << "/" << conv[0]->GetEntries() << std::endl; 
    }
    for( int iMod  = 0; iMod <2; iMod++){
      int nSubModule = (wConv->ModMap[iMod]).nMod;
      for( int iSubMod  = 0 ; iSubMod < nSubModule; iSubMod++){
	int iCrate  = 9999;
	int iSlot   = 9999;
	int iCh     = 9999;
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
	  double halfHeight = fitFunc->GetParameter(0)/2 + fitFunc->GetParameter(4);
	  double halfTiming = fitFunc->GetX( halfHeight,
					     fitFunc->GetParameter(1)-48, 
					     fitFunc->GetParameter(1));
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
	  TSpline3* spl    = new TSpline3( "spl", gr);
	  double splTiming = TSplineGetX(spl, halfHeight, halfTiming-12,  halfTiming +12 );
	  wConv->mod[iMod]->m_SplTiming[chIndex]=splTiming;

	  delete spl;	    
	  delete linearFunction;

	  wavFitter->Clear();	  
	}
      }
    }

    // Event Process .. for Trigger 
    Int_t LaserModID = 5;
    Int_t LaserID    = 5;
    Double_t LaserSignal[5]= {0};
    Double_t LaserTime[5]  = {0};

    for( int iSubMod = 0; iSubMod < wConv->mod[LaserModID]->m_nDigi; iSubMod++){
      LaserID[iSubMod]     = wConv->mod[iMod]->m_ID[iSubMod];
      LaserSignal[iSubMod] = wConv->mod[iMod]->m_Signal[iSubMod];
      LaserTime[iSubMod]   = wConv->mod[iMod]->m_Timint[iSubMod];
    }
    if(LaserSignal[0] > 200 ){
      continue; 
    }

    trout->Fill();
  }
  std::cout<< "End Loop" << std::endl;
  std::cout<< "Close"    << std::endl;
  trout->Write();
  tfout->Close();

  return 0;
}
