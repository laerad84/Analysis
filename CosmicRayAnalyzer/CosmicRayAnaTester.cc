#include "E14CosmicHough.h"
#include "E14CosmicChisq.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TApplication.h"
#include "TRandom.h"

#include <iostream>
#include <fstream>

int main( int argc, char** argv){

  TGraph* gr = new TGraph();
  for( int i = 0; i< 20; i++){
    gr->SetPoint( i, i*25+gRandom->Gaus(0,12.5), gRandom->Gaus(0,12.5)+i*25 );  
  }

  TApplication* app = new TApplication("app",&argc, argv );
  E14CosmicHough* hough = new E14CosmicHough();
  E14CosmicChisq* chisq = new E14CosmicChisq(); 
  bool HoughTester = hough->CosmicJudgment( gr );
  bool ChisqTester = chisq->SetFunction( gr ); 
  double Chisq_Chisq = chisq->CalChisq();
  TH1D* his = new TH1D("hisDistance","hisDistance",50, 0, 50 );
  

  if( !HoughTester ){ 
    std::cout << "This is not Cosmic ray event" << std::endl; 
    return -1;
  }
  if( !ChisqTester ){ 
    std::cout << "The function was not set" << std::endl;
    return -1;
  }

  Double_t hough_roh, hough_theta;
  Double_t chisq_roh, chisq_theta;
  hough_roh = hough->GetRoh();
  hough_theta = hough->GetTheta();
  chisq_roh = chisq->GetRoh();
  chisq_theta = chisq->GetTheta();
  for( int i = 0; i< gr->GetN(); i++){
    double distance = chisq->CalDistance( gr->GetX()[i], gr->GetY()[i] );
    his->Fill( distance );
  }

  std::cout<< hough_roh << "\t" << hough_theta << std::endl;
  std::cout<< chisq_roh << "\t" << chisq_theta << std::endl;

  TCanvas* can = new TCanvas("can","can",800,800);
  TLine* l1 = chisq->GetLine();
  TLine* l2 = hough->GetLine();
  his->Draw();

  app->Run();

}
