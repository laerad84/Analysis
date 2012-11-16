#include <iostream>
#include <fstream>
#include "EDepositAnalysis.h"
#include "TApplication.h"
#include "TCanvas.h"

int main( int argc, char** argv){
  //TApplication* app = new TApplication("app",&argc, argv);
  //TCanvas* can =  new TCanvas("can","",800,800);
  EDepositAnalysis* EDep = new EDepositAnalysis();	

  /*
  for( int i = 0; i< 100; i++){
    EDep->CsIEne->Reset();
    EDep->EventProcess(i);
    can->cd();
    EDep->DrawEvent();
    EDep->Export();
    can->Update();
    can->Modified();
    getchar();
  }
  */
  EDep->Loop();  
  EDep->Close();
  //app->Run();
}
