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


class GSpline { 
private:
  const TSpline3* fspline;
  const double fY0;
public:
  GSpline ( const TSpline3* spline, double y):fspline(spline), fY0(y) {}
  double operator() (double x) const {
    return fspline->Eval(x) - fY0; 
  }
};
double TSplineGetX( TSpline3* spl,  double y, double xmin, double xmax, int npx = 100){

  GSpline g(spl, y);
  ROOT::Math::WrappedFunction<GSpline> wf1(g);
  ROOT::Math::BrentRootFinder brf;
  brf.SetFunction( wf1,xmin, xmax);
  brf.SetNpx(npx);
  brf.SetLogScan(false);
  brf.Solve(npx, 1.E-10,1.E-10);
  if( brf.Root() == xmax || brf.Root() == xmin ){
    return 1.E-10;
  }
  return brf.Root();  
}


int  main(int argc,char** argv)
{
  /// Usage /// 
  // DumpWaveForm [RunNumber] [ Name of Module ] [ Threshold  ] 
  if( argc != 4 ){
    std::cout<< "Usage:\n"
	     << Form(" %s [ RunNumber ] [ Name of Module ] [ Threshold ] ", argv[0] )
	     << std::endl;
    return -1;
  }

  // Get Arguement 
  Int_t    RunNumber  = atoi( argv[1] );
  char*    ModName    = argv[2];
  Double_t Thereshold = atof( argv[2] );


  // GetEnvironment // 

  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEIDR  = std::getenv("ROOTFILE_SUMUP");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;

  // Setting  Classes //
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);

  TFile* tf[nCrate];
  E14ConvReader* conv[nCrate];

  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
  }

  std::cout<< "SetIO" <<std::endl ;
  //TFile* tfout = new TFile(Form("run%d_wav.root",RunNumber),"recreate");
  TFile* tfout = new TFile(Form("%s/run%d_wavDump_%s.root",WAVEFILEDIR.c_str(),RunNumber,ModName),
			   "recreate");
  TTree* trout = new TTree("WFTree",Form("Waveform of %s",ModName));

  Int_t Data[48];
  Int_t Time[48];
  Int_t ID;
  Int_t EventNumber;
  trout->Branch("Data"       ,Data        ,"Data[48]/I");
  trout->Branch("Time"       ,Time        ,"Time[48]/I");
  trout->Branch("ID"         ,&ID         ,"ID/I");
  trout->Branch("EventNumber",&EventNumber,"EventNumber/I");
  TTree* trdummy = new TTree("dummy","");
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),
					    trdummy);
  tfout->cd();
  {
    wConv->AddModule(ModName);
    std::cout << __LINE__ << std::endl;
    wConv->Set();
    std::cout << __LINE__ << std::endl;
    wConv->SetMap();
    std::cout << __LINE__ << std::endl;
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

  std::cout << "Loop" << std::endl;
  for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
    for( int icrate = 0; icrate < nCrate; icrate++){
      conv[icrate]->GetEntry(ievent);
    } 
    
    if( ievent %100 == 0 && ievent ){ std::cout<< ievent << "/" << conv[0]->GetEntries() << std::endl;} 
    
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
	
	for( int ipoint = 0; ipoint< 48; ipoint++){
	  if(conv[iCrate]->Data[iSlot][iCh][ipoint]>16000){ continue; }
	  Data[ipoint] = conv[iCrate]->Data[iSlot][iCh][ipoint];
	  Time[ipoint] = ipoint*8;
	}
	ID = iSubMod;
	EventNumber = ievent;
	trout->Fill();
      }
    }	
  }
  
  std::cout<< "end Loop" <<std::endl;
  //app->Run();
  std::cout<< "Close" << std::endl;
  trout->Write();
  tfout->Close();
  return 0;
}
