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

#include "CsI_Module.h"


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

  // GetEnvironment // 
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;

 //  TCanvas* test = new TCanvas("test","",400,400);
  TFile* tfin  = new TFile(Form("%s/TEMPLATE_FIT_RESULT_%d.root",WAVEFILEDIR.c_str(),RunNumber));
  TTree* trin  = (TTree*)tfin->Get("WFTree");
  std::cout << Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber) << std::endl;  
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),
					    trin);
  std::cout<< "Setting Map" << std::endl;
  tfin->cd();
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
    wConv->SetBranchAddress();
  }

  TApplication* app = new TApplication("App",&argc, argv );

  TCanvas* can = new TCanvas("can","",1600,800);  
  can->Divide(2,1);
  CsI_Module* csi[2];
  csi[0] = new CsI_Module("csi0");
  csi[1] = new CsI_Module("csi1");

  std::cout<<"Entries:" << trin->GetEntries() << std::endl; 
  for( int ievent = 0; ievent < trin->GetEntries() ; ievent++){
    trin->GetEntry(ievent);
    csi[0]->Reset();
    csi[1]->Reset();
    std::cout << ievent << " : " 
	      << wConv->mod[0]->m_nDigi << std::endl;
    for( int ich =0 ;ich < wConv->mod[0]->m_nDigi; ich++ ){
      csi[0]->Fill(wConv->mod[0]->m_ID[ich], wConv->mod[0]->m_Signal[ich]);
      csi[1]->Fill(wConv->mod[0]->m_ID[ich], wConv->mod[0]->m_ADC[ich]);
    }
    can->cd(1);
    csi[0]->Draw("colz");
    gPad->SetLogz();
    can->cd(2);
    csi[1]->Draw("colz");
    gPad->SetLogz();
    can->Update();
    can->Modified();
    getchar();
  }
  app->Run();
}
