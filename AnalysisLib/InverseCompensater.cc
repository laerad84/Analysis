#include "PeakCompensater.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "TCanvas.h"
#include "TF1.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "TSpline.h"
#include "EnergyConverter.h"

int main( int argc, char** argv ){
  std::string ANALYSISLIB = std::getenv("ANALYSISLIB");
  EnergyConverter* conv = new EnergyConverter();
  conv->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",ANALYSISLIB.c_str()));
 
  std::cout<< conv->GetCalibrationConstant(0) << std::endl;
  for( int i = 1; i< 160; i++){
    std::cout<< i*10 <<"\t" << conv->ConvertToHeight(0,i*10)<< std::endl;
  }
}


