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

TCanvas*   can;
IDHandler* handler;
CsIImage*  image;
TThread*   th;
Bool_t     finished;

void* handle(void *ptr){  
  can->cd();
  image->Draw();
  can->Modified();
  can->Update();  
}

int
main(int argc, char** argv){
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
  E14ReadSumFile* Reader = new E14ReadSumFile();
  Reader->Add(InputFile.c_str());
  
  handler = new IDHandler("Data/crystal.txt");
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
    for( int idigi = 0; idigi< Reader->CsiNumber; idigi++){
      if( Reader->CsiEne[idigi] > threshold){
	double x,y;
	image->Fill(Reader->CsiModID[idigi],Reader->CsiEne[idigi]);
	std::cout << Reader->CsiEne[idigi] << std::endl;
      }
    }  
    iEntry++;
    th = new TThread("t0",handle,(void*) 0);
    th->Run();    
  }
  
  app->Run();
}
