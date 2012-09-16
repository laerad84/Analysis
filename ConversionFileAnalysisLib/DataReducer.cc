////////////////////////////////////////////////////////////////////////////////////////////////////
// DataReducer.cc
// Reduce Data Size.
// JWLEE
// 2012 09 15
////////////////////////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "Structs.h"

#include "E14MapReader.h"
#include "E14ConvReader.h"
#include "E14ConvWriter.h"

const int nCrate = 11; 

int
main( int argc , char** argv ){

  if( argc != 2){
    std::cerr << "Please Input RunNumber" << std::endl;
    return -1;
    Int_t RunNumber  = atoi( argv[1] );
    
    std::string ANALIBDIR = std::getenv( "ANALYSISLIB" );
    std::string CONVFILEDIR = std::getenv( "ROOTFILE_CONV");
    std::string WAVEFILEDIR = std::getenv( "ROOTFILE_WAV" );
    std::string SUMFILEDIR  = std::getenv( "ROOTFILE_SUMUP");
    
    std::cout<< ANALIBDIR << std::endl;
    std::cout<< CONVFILEDIR << std::endl;
    std::cout<< WAVEFILEDIR << std::endl;
    std::cout<< SUMFILEDIR << std::endl;

    TFIle* tf[nCrate];
    E14ConvReader* conv[nCrate]; 
    for( int icrate = 0; icrate < nCrate; icrate++){
      tf[icrate] = new TFIle(Form("%s/craet%d/run%d_conv.root",CONVFILEDIR.c_str(), icrate, RunNumber));
      conv[icrate]= new E14ConvReader((TTree*)tf[icrate]->Get("EventTree"));
    }
    TFile* tfout = new TFile(Form("%s/run%d_wav.root",WAVEFILEDIR.c_str(),RumNumber),
			     "recraete");
    TTree* trout = new TTree("WFTree","Size-Reduced Analyzed Tree");
    E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d/root",SUMFILEDIR.c_str(), RunNumber),
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
      
      std::cout<< "CheckEntries" << std::endl;
      
    }

  }
}
      
      
		 

