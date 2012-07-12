#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"

#include "E14ReadSumFile.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "HoughCsI.h"
#include "TThread.h"
#include "TROOT.h"

#include "CsI_Module.h"
#include "OEV_Module.h"
#include "CC03_Module.h"

TCanvas*   can;
IDHandler* handler;
CsIImage*  image;
TThread*   th;
Bool_t     finished;
CsI_Module* csi = new CsI_Module("csi");
OEV_Module* oev = new OEV_Module("oev");
CC03_Module* cc03 = new CC03_Module("cc03");


void* handle(void *ptr){  
  can->cd();
  image->Draw();
  /*
  csi->DrawWithRange(0,1000,"    ");
  oev->DrawWithRange(0,10000,"same ");
  cc03->DrawWithRange(0,10000,"same ");
  */
  can->Update();
  can->Modified();
  can->Update();  
}

int
main(int argc, char** argv){
  gROOT->SetStyle("plain");
  std::cout<< "Start" << std::endl;
  Double_t threshold=5;
  finished = kFALSE;
  gStyle->SetPalette(1);
  // ARGV[0]  <InputROOTFile> <OutputROOTFile> [<Vision>]

  std::string InputFile;
  std::string OutputFile;
  
  if( argc == 3){
    InputFile = argv[1];
    threshold = atof(argv[2]);
  }else{
    std::cerr << "ARGUEMENT ERROR" << std::endl;
    return -1;
  }
  
  //

  TApplication* app = new TApplication("app",&argc, argv);
  std::cout<< "Reader"<< std::endl;
  E14ReadSumFile* Reader = new E14ReadSumFile(0);
  Reader->Add(InputFile.c_str());
  
  handler = new IDHandler();
  image   = new CsIImage(handler);
  can     = new TCanvas("can","",800,0,800,800);
  
  long nEntries = Reader->GetEntries();
  std::cout << "TotalEntries::" << nEntries  << std::endl;
  int iEntry = 0; 

  while(1){
    char c = getchar();
    if( c == 'n'){
      std::cout<< "EventNumber:"<< std::flush;
      std::cin >> iEntry;
    }
    image->Reset();

    if( th != 0 ){
      th->Join();
    }

    Reader->GetEntry(iEntry);    
    image->SetTitle(Form("Event:%d",iEntry));
    csi->Reset();
    cc03->Reset();
    oev->Reset();
    for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
      if( Reader->CsiEne[idigi] > threshold){
	double x,y;
	image->Fill(Reader->CsiModID[idigi],Reader->CsiEne[idigi]);

	std::cout <<  Reader->CsiModID[idigi] << " : " 
		  <<  Reader->CsiEne[idigi]   << std::endl;
	csi->Fill(Reader->CsiModID[idigi], Reader->CsiEne[idigi]);	
      }	
    }  
    for( int idigi = 0; idigi < Reader->CC03Number; idigi++){
      if( Reader->CC03Ene[idigi] > 1000){
	cc03->Fill( Reader->CC03ModID[idigi] , Reader->CC03Ene[idigi] );
      }
    }
    for( int idigi = 0; idigi < Reader->OEVNumber; idigi++){
      if( Reader->OEVEne[idigi] > 1000 ){
	//std::cout<<  Reader->OEVModID[idigi]<< " : "
	//	 <<  Reader->OEVEne[idigi]  << std::endl;
	oev->Fill( Reader->OEVModID[idigi] , Reader->OEVEne[idigi] );
      }
    }

    iEntry++;

    //th = new TThread("t0",handle,(void*) 0);
    //th->Run();    


    image->Draw();
    can->Update();
    can->Modified();
    getchar();

    
  }
  
  app->Run();
}
