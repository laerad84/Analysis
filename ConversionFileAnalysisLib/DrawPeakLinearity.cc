#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TSpline.h"
#include "TApplication.h"
#include "PeakCompensater.h"
int main( int argc, char** argv ){

  TApplication* app = new TApplication("app", &argc, argv);
  PeakCompensater* Peak = new PeakCompensater();
  std::cout<< Peak->Compensate(0,12000) << std::endl;

  //app->Run();
  return 0; 
}
