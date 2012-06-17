#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"

#include "klong/Klong.h"

#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"

#include <cstdlib>
#include <cstdio>
#include <list>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <list>

#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"


int main( int argc ,char** argv){
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("neRMI");
  const int nCSI = 2716;  
  
  Int_t CalibrationNumber  = -1; //iterationNumber 
  std::string runListFile;
  if( argc == 3){
    CalibrationNumber = atoi(argv[1]);
    runListFile       = argv[2];
  }else{
    std::cerr << "<<<>>> Argument Error <<<>>> " << "\n"
	      << "Usage:" << argv[0] 
	      << " [CalibrationNumber] [Filename of List] "
	      << std::endl;
    return -1;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // Set initial CalibrationFactor and Files 
  //////////////////////////////////////////////////////////////////////////////

  std::vector<int> VecRunNum;
  std::ifstream ifsList(runListFile.c_str());
  if( !ifsList.is_open() ){
    return -1; 
  }
  while( !ifsList.eof() ){
    Int_t runNumber;
    if( ifsList >>  runNumber){
      VecRunNum.push_back(runNumber );
    }
  }
  
  int    nCal[nCSI];
  double CalFactor[nCSI];
  double CalRMS[nCSI];
  double preCalFactor[nCSI];  
  
  for( int i = 0; i< nCSI; i++){
    preCalFactor[i] = 1;
  }
  if( CalibrationNumber!=0 ){
    std::ifstream ifs(Form("Calibration_Data/CalibrationFactor_%d.dat",
			   CalibrationNumber));
    if(!ifs.is_open()){ return -1;}
    Int_t id;
    Double_t precal;
    while( !ifs.eof() ){
      if( ifs >> id >> precal){
	preCalFactor[id] = precal;
      }
    }
  }
  int initialRunNumber;
  int finalRunNumber;
  Int_t nextCalNum = CalibrationNumber +1;
  std::ofstream ofs(Form("Calibration_Data/CalibrationFactor_%d.dat",
			 nextCalNum));
  std::ofstream ofs1(Form("Calibration_Data/CalibrationStatics_%d.dat",
			  nextCalNum));
  std::string listFilename=runListFile.substr(runListFile.find_last_of('/')+1);
  std::cout<< listFilename << std::endl;
  TFile* tfOut 
    = new TFile(Form("Calibration_Data/CalhistList_%s_%d.root",
		     listFilename.substr(0,-4).c_str(),
		     CalibrationNumber),
		"RECREATE");
  

  ////////////////////////////////////////////////////////////////////////////  
  // Calibration Factor Histogram
  ////////////////////////////////////////////////////////////////////////////  
  
  TH1D* hisCalibrationFactor[nCSI];
  for( int i =0; i < nCSI; i++){
    hisCalibrationFactor[i] = new TH1D(Form("hisCalibrationFactor_%d",i),
				       Form("CalibrationFactor:%d",i),
				       500,0,2);
  }  
  
  ////////////////////////////////////////////////////////////////////////////
  // Loop
  ////////////////////////////////////////////////////////////////////////////
  //- Add Histogram for all Run;
  
  for( std::vector<int>::iterator iRun = VecRunNum.begin();
       iRun!=VecRunNum.end();
       ++iRun){
    TH1D*  hisCalFactor[nCSI];  
    TFile* tf = new TFile(Form("Calibration_Data/Calibration_%d_%d.root",
			       *iRun,CalibrationNumber));        
    std::cout<< "Channel Loop" << std::endl;
    for( int i = 0; i< nCSI; i++){
      hisCalFactor[i] = (TH1D*)tf->Get(Form("his_Calibration_%04d",i));
      nCal[i]         = hisCalFactor[i]->Integral();      
      hisCalibrationFactor[i]->Add(hisCalFactor[i]);
    }
    tf->Close();
  }
  
  tfOut->cd();
  //- Aquire Calibration Factor
  for( int i = 0; i< nCSI; i++){
    hisCalibrationFactor[i]->Write();
    ofs1 << i << "\t" << hisCalibrationFactor[i]->Integral() << std::endl;
    if( hisCalibrationFactor[i]->Integral() < 20 ){
      continue;
    }           
    hisCalibrationFactor[i]->Fit("gaus","Q","");
    TF1* calFunction= hisCalibrationFactor[i]->GetFunction("gaus");
    CalFactor[i]    = hisCalibrationFactor[i]->GetMean();
    CalRMS[i]       = hisCalibrationFactor[i]->GetRMS();
    //CalFactor[i] = calFunction->GetParameter(1);
    //CalRMS[i] = calFunction->GetParameter(2);
    ofs << i << "\t" << CalFactor[i]*preCalFactor[i] << std::endl;    
  }

  tfOut->Close();
  ofs.close();      
  ofs1.close();
}
