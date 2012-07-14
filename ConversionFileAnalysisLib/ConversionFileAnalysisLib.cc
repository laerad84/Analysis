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

  // Setting  Classes //
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  
  TFile* tf[nCrate];

  E14ConvReader* conv[nCrate];
  for( int icrate = 0; icrate < nCrate; icrate++){
    tf[icrate]   = new TFile(Form("%s/crate%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber)); 
    conv[icrate] = new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
  }
  //TFile* tfout = new TFile(Form("run%d_wav.root",RunNumber),"recreate");
  TFile* tfout = new TFile(Form("%s/run%d_wav.root",WAVEFILEDIR.c_str(),RunNumber),
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
  
  TApplication* app = new TApplication("app", &argc , argv );  
  TCanvas* can = new TCanvas( "can ", "Canvas" ,800,800);
  std::cout <<"Loop " <<std::endl;
  TGraph* gr = new TGraph();
  gr->SetMarkerStyle(6);
  
  for( int ievent  = 0; ievent < conv[0]->GetEntries(); ievent++){
  //for( int ievent  = 0; ievent < 100; ievent++){
    wConv->InitData();
    for( int icrate = 0; icrate < nCrate; icrate++){
      //std::cout << icrate << std::endl;
      conv[icrate]->GetEntry(ievent);
    } 
    if( ievent %100 == 0 && ievent ){ std::cout<< ievent << "/" << conv[0]->GetEntries() << std::endl;} 
    
    // Analysis
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
	  TSpline3* spl    = new TSpline3( "spl", gr);
	  double splTiming = TSplineGetX(spl, halfHeight, halfTiming-12,  halfTiming +12 );
	  wConv->mod[iMod]->m_SplTiming[chIndex]=splTiming;

	  delete spl;	    
	  delete linearFunction;
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
