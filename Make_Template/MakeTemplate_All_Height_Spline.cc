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
#include "TChain.h"
#include "TROOT.h"
#include "Math/WrappedFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMethods.h"
#include "Math/WrappedTF1.h"
#include "TProfile.h"

int
main( int argc, char** argv){
  if( argc != 3 ){ 
    std::cerr << "Number of argumemt must be 3" << std::endl;
    std::cerr << Form( "%s [Height Low Limit] [Height High Limit]", argv[0])
	      << std::endl;
    return -1;
  }

  Double_t LowLimit  = atof( argv[1] );
  Double_t HighLimit = atof( argv[2] );
  if( LowLimit > HighLimit || LowLimit < 0 ){
    std::cerr << "Please Confirm Limit\n"
	      << "Invalid Number for Limit" << std::endl;
    return -1; 
  }

  TFile* tf = new TFile(Form("TEMPLATE_HEIGHT_%d_%d.root",(int)LowLimit,(int)HighLimit));
  TFile* tfout = new TFile(Form("TEMPLATE_SPLINE_%d_%d.root",(int)LowLimit, (int)HighLimit) ,"recreate");
  
  TGraph* grTemp[2716];
  TGraph* grTempRaw[2716];
  TH2D*   hisTemplate[2716];
  TSpline3* splTemp[2716];
  for( int ich = 0; ich < 2716; ich++){
    hisTemplate[ich] = (TH2D*)tf->Get(Form("hisTempCsI_Height_%d_%d_%d",(int)LowLimit, (int)HighLimit,ich));
    grTemp[ich] = new TGraph();
    grTempRaw[ich] = new TGraph();
    grTemp[ich]->SetName(Form("Template_graph_%d", ich));
    grTemp[ich]->SetTitle(Form("Template_graph_%d", ich));
    /*
    grTempRaw[ich]->SetName(Form("Template_graphRaw_%d", ich));
    grTempRaw[ich]->SetTitle(Form("Template_graphRaw_%d", ich));
    */
  }
  std::cout << "Set grTempRaw " << std::endl;
  for( int ich = 0; ich < 2716; ich++){
    //std::cout << ich << std::endl;
    TH1D* hisSlice[400];
    grTempRaw[ich]->Set(0);
    if( hisTemplate[ich]->GetEntries()==0){ continue ;}
    int ipoint  = 0; 

    for( int  ibin  =0; ibin < 400 ; ibin++) {
      hisSlice[ibin] = hisTemplate[ich]->ProjectionY(Form("bin%d", ibin),
						      ibin,ibin+1);
      if( hisSlice[ibin]->GetEntries() < 4 ) { continue; }
      if( hisTemplate[ich]->GetBinCenter( ibin ) < -150 ){ continue; }      
      grTempRaw[ich]->SetPoint( ipoint,
				hisTemplate[ich]->GetBinCenter( ibin ),
				hisSlice[ibin]->GetMean());
      ipoint++;
    }
  }

  std::cout<< "Set grTemp" << std::endl;  
  for( int ich = 0; ich  < 2716; ich++){
    //std::cout<< ich << std::endl;
    double minRMS = 10000000;
    double minpoint = 0.;
    double ped;
    double sig;
    
    int nSample  = 7; 
    for( int ipoint = 0; ipoint < grTempRaw[ich]->GetN()-nSample; ipoint++){
      if( grTempRaw[ich]->GetX()[ipoint+nSample-1] > 0 ){ break; }
      double rms    = 0; 
      double mean   = 0;
      double sqmean = 0;
      
      for( int isubpoint  = 0; isubpoint < nSample; isubpoint++){ 
	mean   += grTempRaw[ich]->GetY()[ipoint + isubpoint];
	sqmean += grTempRaw[ich]->GetY()[ipoint + isubpoint];
	
	mean   /= nSample;
	sqmean /= nSample;
	
	rms = TMath::Sqrt( sqmean - mean*mean );
	if( rms  < minRMS ){
	  minRMS = rms;
	  ped    = mean;
	}
      }
    }

    if( grTempRaw[ich]->GetN() == 0){ continue; }
    TSpline3* spl = new TSpline3("spl", grTempRaw[ich] );
    sig = spl->Eval(0)-ped;
    for( int ipoint = 0; ipoint < grTempRaw[ich]->GetN(); ipoint++){
      grTemp[ich]->SetPoint( ipoint , grTempRaw[ich]->GetX()[ipoint],
			     (grTempRaw[ich]->GetY()[ipoint]-ped)/sig);
    }
    delete spl;

    std::cout<< grTemp[ich]->GetX()[ grTemp[ich]->GetN() -1 ] << std::endl;
    if( grTemp[ich]->GetX()[grTemp[ich]->GetN()-1 ] >0 ){
      std::cout<<"Fit:" << ich << std::endl;
      grTemp[ich]->Fit("pol1","Q","",180,300);
      TF1* func = grTemp[ich]->GetFunction("pol1");
      if( func->GetParameter(1) >= 0 ){ std::cout << "Unexpected Fit result : " << ich << std::endl;}
      else{
	while( 1 ){	  
 	  Double_t xpos = grTemp[ich]->GetX()[ grTemp[ich]->GetN() -1 ]+1;
	  Double_t ypos = func->Eval(xpos);
	  if( xpos <= 300 ){
	    if( xpos == 300 ){
	      if( ypos > 0 ){ std::cout<< "Check Channel:" <<  ich << std::endl; }
	    }
	    
	    if( ypos >= 0 ){
	      grTemp[ich]->SetPoint( grTemp[ich]->GetN(), xpos, ypos );
	    }else{
	      grTemp[ich]->SetPoint( grTemp[ich]->GetN(), xpos, 0 );
	    }
	  }else{
	    break;
	  }	    
	}
	grTemp[ich]->GetListOfFunctions()->Delete();
      }
    } 
  }				
  std::cout<< "End Fill" << std::endl;

  tfout->cd();


  for( int ich = 0; ich < 2716; ich++){
    std::cout << ich << std::endl;
    grTemp[ich]->Write();
    std::cout << ich << std::endl;
  }
	   
  tfout->Close(); 

  return 0;
}
