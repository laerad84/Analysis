#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "TGraph.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGraph.h"
#include "Chisq_cosmic.h"
#include "TLine.h"

int
main(int argc, char** argv){
  
  TApplication* app = new TApplication("app",&argc, argv);
  TGraph* gr = new TGraph();
  gr->SetPoint(0,1.1,0.9);
  gr->SetPoint(1,1.9,2.2);
  gr->SetPoint(2,3.4,3.5);
  Chisq_cosmic* chi2_cosmic = new Chisq_cosmic();
  chi2_cosmic->SetFunction(gr);
  chi2_cosmic->SetRange(1,-45);
  chi2_cosmic->CalChisq();
  std::cout<< chi2_cosmic->GetDistance(0,0) << std::endl;
  TCanvas* can = new TCanvas("can","",800,800);
  TLine* line = chi2_cosmic->GetLine();
  gr->Draw("AP");
  std::cout <<line << std::endl;
  line->Draw("same");


  app->Run();
}
