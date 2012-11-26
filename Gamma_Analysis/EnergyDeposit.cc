#include <iostream>
#include <fstream>
#include "EDepositAnalysis.h"
#include "TApplication.h"
#include "TCanvas.h"
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <string>

int main( int argc, char** argv){
  //TApplication* app = new TApplication("app",&argc, argv);
  //TCanvas* can =  new TCanvas("can","",800,800);
  std::cout<< "Start" << std::endl;
  std::string ROOTFILE_GAMMAHIT = std::getenv("ROOTFILE_GAMMAHIT");
  std::string ROOTFILE_GAMMACLUS= std::getenv("ROOTFILE_GAMMACLUS");
  std::cout << ROOTFILE_GAMMAHIT << std::endl;
  std::cout << ROOTFILE_GAMMACLUS << std::endl;
  Int_t Energy = atoi( argv[1] );
  Int_t Degree = atoi( argv[2] );
  Int_t Index  = atoi( argv[3] );
  std::cout<< Energy  << " : " << Degree << " : " << Index << std::endl;
  std::string InputFilename = Form("%s/template_gamma_%dMeV_%ddeg-1E5-%d.root",
				   ROOTFILE_GAMMAHIT.c_str(),Energy, Degree, Index);
  std::string OutputFilename = Form("%s/Cluster_%dMeV_%ddeg-1E5-%d.root",
				    ROOTFILE_GAMMACLUS.c_str(),Energy, Degree, Index);
  EDepositAnalysis* EDep = new EDepositAnalysis(InputFilename.c_str(),
						OutputFilename.c_str());	

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
