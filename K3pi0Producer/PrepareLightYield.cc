#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "EnergyConverter.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"


int main( int argc, char** argv ){ 
  std::string ANALYSISLIB = std::getenv("ANALYSISLIB");
  std::ifstream  ifs("GainAdjustVoltage.txt");
  double Out0[2716]  ={0};
  double Out1[2716]  ={0};
  double HVAll[2716] ={0};

  int tmpID;
  double tmpOut0;
  double tmpOut1;
  double HV;
  while(ifs >> tmpID >> tmpOut0 >> HV >> tmpOut1){
    if(tmpOut0 != 0  ){      
      if( (tmpID < 2716) && ( tmpID >= 0 ) ){
	Out0[tmpID] = tmpOut0;
	Out1[tmpID] = tmpOut1;
	HVAll[tmpID]= HV;
      }
    }
  }
  EnergyConverter* Converter = new EnergyConverter();
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",ANALYSISLIB.c_str()));
  double CalFactor0 = 0.08485;// 14MeV/165Cnt;  

  std::ofstream ofs("RelativeLY.txt");
  double LYRatio[2716]={0};  
  Double_t StdOut = 2578;// 12p.e/MeV
  for( int i = 0; i< 2716; i++){
    if( Out0[i]> 0){
      LYRatio[i] = Out0[i]/StdOut;
    }else{
      LYRatio[i] = Converter->GetCalibrationConstant(i)/CalFactor0;
    }
    ofs << i << "\t" << LYRatio[i] << "\n";
  }
  ofs.close();
  return 0; 
}
