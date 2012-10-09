#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cstdio>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "EnergyConverter.h"
#include "PeakCompensater.h"
int
main( int argc , char** argv ){

  std::string ANALYSISLIB = std::getenv("ANALYSISLIB");
  std::string CalibrationRootFile = Form( "%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root", ANALYSISLIB.c_str()); 
  EnergyConverter* converter = new EnergyConverter();
  converter->ReadCalibrationRootFile( CalibrationRootFile.c_str());
  for( int i = 0; i< 2716; i++){
    std::cout << i << " \t:\t " << converter->GetCalibrationConstant( i ) << "\n";
  }

}
