#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include "E14EventBuilder_V0.h"
#include "TApplication.h"
#include "TPostScript.h"
int
main( int argc, char** argv ){  

  int RunNumber;
  int EventNumber;
  RunNumber = atoi( argv[1] );
  EventNumber = atoi( argv[2] );

  TFile* tfOut = new TFile("test.root","RECREATE");
  TTree* trOut = new TTree("trOut","");
  TApplication* app  = new TApplication("app", &argc, argv);
  TCanvas* can  = new TCanvas("can","",800,1200);
  E14EventBuilder_V0* test = new E14EventBuilder_V0(trOut,RunNumber);
  
  tfOut->cd();
  //test->LoopAll();
  
  for( int ievent = EventNumber; ievent < EventNumber+100; ievent++){
    std::cout << "Analysis started : " << ievent  << std::endl;
    test->EventProcess(ievent);
    std::cout<< "Draw result" << std::endl;
    can->Clear();
    test->DrawEvent(can);
    can->Update();
    can->Modified();
    getchar();
  }

  //trOut->Write();
  //tfOut->Close();


  app->Run();


}
