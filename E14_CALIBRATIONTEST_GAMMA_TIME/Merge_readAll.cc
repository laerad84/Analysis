#include <iostream>
#include <string>
#include <cstring>

#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"

#include "TPostScript.h"

#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#inclucd "CalibrationTree.h"

//void Merge(int iterationNumber){

int 
main( int argc, char** argv ){
  int iterationNumber  = 0;

  TFile* tf = new TFile("SimCalibration_Data/MergeOutput.root","recreate");
  TH2D*  his_klpos_delta_All=new TH2D("his_klpos_delta_All","recklz-simklz:recklz",
				      75,0,7500,100,-100,900);
  TH2D*  his_klpos_delta_Fit=new TH2D("his_klpos_delta_Fit","recklz-simklz:recklz",
				      75,0,7500,100,-100,900);
  TH2D*  his_klpos_delta_Mis=new TH2D("his_klpos_delta_Mis","recklz-simklz:recklz",
				      75,0,7500,100,-100,900);
  TH1D*  his_klpos_All = new TH1D("his_klpos_All","recklz",
				  75,0,7500);
  TH1D*  his_klpos_Fit = new TH1D("his_klpos_Fit","recklz",
				  75,0,7500);  
  TH1D*  his_klpos_Mis = new TH1D("his_klpos_Mis","recklz",
				  75,0,7500);
  TH1D*  his_delta_All = new TH1D("his_delta_All","recklz-simklz",100,-100,900);
  TH1D*  his_delta_Fit = new TH1D("his_delta_Fit","recklz-simklz",100,-100,900);
  TH1D*  his_delta_Mis = new TH1D("his_delta_Mis","recklz-simklz",100,-100,900);

  TH2D*  his_klpos_klChisqZ_All = new TH2D("his_klpos_klChisqZ_All","klChisqZ:recklz",
					   75,0,7500,100,0,100);
  TH2D*  his_klpos_klChisqZ_Fit = new TH2D("his_klpos_klChisqZ_Fit","klChisqZ:recklz",
					   75,0,7500,100,0,100);
  TH2D*  his_klpos_klChisqZ_Mis = new TH2D("his_klpos_klChisqZ_Mis","klChisqZ:recklz",
					   75,0,7500,100,0,100);

  TH1D*  his_klChisqZ_All  = new TH1D("his_klChisqZ_All","klChisqZ",100,0,100);
  TH1D*  his_klChisqZ_Fit  = new TH1D("his_klChisqZ_Fit","klChisqZ",100,0,100);
  TH1D*  his_klChisqZ_Mis  = new TH1D("his_klChisqZ_Mis","klChisqZ",100,0,100);
  TH1D*  his_SecklChisqZ_All  = new TH1D("his_SecklChisqZ_All","Second klChisqZ",100,0,100);
  TH1D*  his_SecklChisqZ_Fit  = new TH1D("his_SecklChisqZ_Fit","Second klChisqZ",100,0,100);
  TH1D*  his_SecklChisqZ_Mis  = new TH1D("his_SecklChisqZ_Mis","Second klChisqZ",100,0,100);

  TH1D* his_GammaE_All = new TH1D("his_GammaE_All","GammaE",100,0,2000);
  TH1D* his_GammaE_Fit = new TH1D("his_GammaE_Fit","GammaE",100,0,2000);
  TH1D* his_GammaE_Mis = new TH1D("his_GammaE_Mis","GammaE",100,0,2000);
  
  TH1D* his_klE_All = new TH1D("his_klE_All","klEnergy",80,0,8000);
  TH1D* his_klE_Fit = new TH1D("his_klE_Fit","klEnergy",80,0,8000);
  TH1D* his_klE_Mis = new TH1D("his_klE_Mis","klEnergy",80,0,8000);
  
  TH1D* his_klMass_All = new TH1D("his_klMass_All","klMass",200,400,600);
  TH1D* his_klMass_Fit = new TH1D("his_klMass_Fit","klMass",200,400,600);
  TH1D* his_klMass_Mis = new TH1D("his_klMass_Mis","klMass",200,400,600);

  TH1D* his_klDeltaChisqZ_All = new TH1D("his_klDeltaChisqZ_All","klDeltaChisqZ",100,0,100);
  TH1D* his_klDeltaChisqZ_Fit = new TH1D("his_klDeltaChisqZ_Fit","klDeltaChisqZ",100,0,100);
  TH1D* his_klDeltaChisqZ_Mis = new TH1D("his_klDeltaChisqZ_Mis","klDeltaChisqZ",100,0,100);
  
  /*
  TH1D* his_GammaR_All = new TH1D("his_GammaR_All","R of Gamma Position",100,0,1000);
  TH1D* his_GammaR_Fit = new TH1D("his_GammaR_Fit","R of Gamma Position",100,0,1000);
  TH1D* his_GammaR_Mis = new TH1D("his_GammaR_Mis","R of Gamma Position",100,0,1000);
  */

  std::string Filename = "SimCalibration_Data/Calibration_%04d_%d.root";
  TChain* ch = new TChain("trCalibration");
  for( int i = 0; i < 200; ++i){
    ch->Add(Form(Filename.c_str(),i*10,iterationNumber));
  }
  
  E14GNAnaDataContaioner data;
  CalibrationTree        caldata;  
  data.setBranchAddress( ch );
  caldata.setBranchAdderess( ch );

  std::cout <<ch->GetEntries() << std::endl;
  std::cout<< __LINE__<< std::endl;



  return 0;
}
    
  
