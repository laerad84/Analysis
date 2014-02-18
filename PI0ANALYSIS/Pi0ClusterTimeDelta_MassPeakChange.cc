/// Convert simulation data with trigger simulation.


#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <string>
#include <list>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
//#include "ClusterFinder_EDIT.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "User_Functions.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "IDHandler.h"

//#include "E14WavReader.h"
#include "E14WavReader_V1.h"
//#include "E14WaveReader_V2.h"
//#include "L1TrigCounter.h"
#include "EnergyConverter.h"
#include "TRandom.h"
#include "TF1.h"

double funcResolutionInvSq( double* x, double* p){
  double value = 0;
  if( x[0] >  0  && p[0] > 0){
    //value = 10000./(1.26*1.26*1e3/(x[0]*p[0])+16900/(x[0]*x[0]) + 0.76*0.76);
    value = 12.7*x[0]*p[0];
  }
  return value;
}

int
main( int argc ,char ** argv ){
  Int_t Types  = atoi( argv[1] );
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");
  
  std::string ROOTFILE_SIMCONV  = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/ConvFile";
  std::string ROOTFILE_SIMPI0   = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/SIMPI0";//kekcc
  //std::string ROOTFILE_SIMPI0   = "/Volume0/Simulation/Pi0Run/NewPi0Data_2013";//local
  //std::string iFileForm          = "%s/SimPi0_1E6_LYRES_%d.root";     // ROOTFILE_SIMCONV
  //std::string oFileForm          = "%s/SimPi0_1E6_LYRES_Merged.root"; // ROOTFILE_SIM3PI0

  std::string iFileForm;
  std::string oFileForm;
  switch( Types ){
  case 0:
    iFileForm = "%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0_nocal_nopi0peak_notempcorr.root";
    oFileForm = "%s/Pi0Mass_nocal_nopi0peak_notemp.root";
    break;
  case 1:
    iFileForm = "%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0_nopi0peak_notempcorr.root";
    oFileForm = "%s/Pi0Mass_nopi0peak_notempcorr.root";
    break;
  case 2:
    iFileForm = "%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0.root";
    oFileForm = "%s/Pi0Mass_All_Corr.root";
    break;
  case 3:
    iFileForm = "%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0_nocal.root";
    oFileForm = "%s/Pi0Mass_nocal.root";
    break;
  case 4:
    iFileForm = "%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0_nopi0peak.root";
    oFileForm = "%s/Pi0Mass_nopi0peak.root";
    break;
  default :
    return -1;
  }


  Int_t RunN[24]={4502,4503,4504,4505,4506,4507,4508,4509,4510,4511,4512,4513,
		  4514,4515,4516,4517,4518,4519,4520,4521,4522,4523,4524,4525};
  //  TChain* trin = new TChain("T");
  TFile* tf[24]; 
  TTree* tr[24];
  for( int i = 0; i < 24; i++){
    tf[i] = new TFile(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(), RunN[i]));
    tr[i] = (TTree*)tf[i]->Get("T");
    //trin->Add(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(),RunN[i]));
  }

  int    RunNo;
  int    EventNumber;
  int    CsiNumber;
  int    CsiModID[2716];
  double CsiEne[2716];
  double CsiTime[2716];
  double CsiHHTime[2716];
  double CsiSignal[2716];
  int    CsiL1nTrig;
  double CsiL1TrigCount[20];
  for( int i = 0; i< 24; i++){
    
  tr[i]->SetBranchAddress("RunNumber",&RunNo);
  tr[i]->SetBranchAddress("EventNumber",&EventNumber);
  tr[i]->SetBranchAddress("CsiNumber",&CsiNumber);
  tr[i]->SetBranchAddress("CsiModID",CsiModID);
  tr[i]->SetBranchAddress("CsiEne",CsiEne);
  tr[i]->SetBranchAddress("CsiTime",CsiTime);
  tr[i]->SetBranchAddress("CsiHHTime",CsiHHTime);
  tr[i]->SetBranchAddress("CsiSignal",CsiSignal);
  tr[i]->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  tr[i]->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  }
  /*
  int    nTrack;
  UShort_t  track[200];
  int    pid[200];
  float mass[200];
  float ek[200];
  float end_ek[200];
  double p[200][3];
  double end_p[200][3];
  double v[200][3];
  double end_v[200][3];

  trin->SetBranchAddress("nTrack",&nTrack);
  trin->SetBranchAddress("track",track);
  trin->SetBranchAddress("pid",pid);
  trin->SetBranchAddress("mass",mass);
  trin->SetBranchAddress("ek",ek);
  trin->SetBranchAddress("end_ek",end_ek);
  trin->SetBranchAddress("p",p);
  trin->SetBranchAddress("end_p",end_p);
  trin->SetBranchAddress("v",v);
  trin->SetBranchAddress("end_v",end_v);
  */

  //double CsiL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
  //					1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};
  double CsiL1TrigCountThreshold[20] = {1000,3000,3000,3400,3400,3400,2200,2200,2400,2400,
					2400,1000,1000,1000,1000,1000,1000,1000,1000,1000};
  double CsiL1TrigHighThreshold = 50000;
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_SIMPI0.c_str()),"recreate");
  TH1D* hisPi0Mass[24];
  for( int i = 0; i< 24; i++){
    hisPi0Mass[i] = new TH1D(Form("hisPi0Mass_%d",RunN[i]),Form("hisPi0Mass_%d;Rec.Mass[MeV]",RunN[i]),200,0,400);
  }

  E14GNAnaDataContainer data; 
  const int nHist  =1;
  char* Name[nHist] = {"Data"};



  std::cout<< "LOOP" << std::endl;
  for( int iFile = 0; iFile < 24; iFile++){
    
    data.setBranchAddress(tr[iFile]);
    long entries = tr[iFile]->GetEntries();
    std::cout<< entries << std::endl;
    for( int ievent  = 0; ievent < entries ; ievent++){
      //for( int ievent  = 0; ievent < 100 ; ievent++){
      tr[iFile]->GetEntry(ievent);
      std::list<Cluster> clist;
      std::list<Gamma>   glist;
      std::list<Pi0>     plist;
      
      data.getData(plist);
      
      Int_t hisID = -1;
      
      std::list<Pi0>::iterator pit = plist.begin();
      hisID = 0;
      int nTrig = 0; 
      for( int i = 1; i< 11; i++){
	if( CsiL1TrigCount[i] > CsiL1TrigCountThreshold[i]  &&
	    CsiL1TrigCount[i] < CsiL1TrigHighThreshold ){
	  nTrig++;
	}
      }
      if( nTrig >=2 ){
	if( hisID >= 0){      
	  double x[2]; 
	  x[0] = (*pit).g1().x();
	  x[1] = (*pit).g2().x();
	  double y[2]; 
	  y[0] = (*pit).g1().y();
	  y[1] = (*pit).g2().y();
	  // Gamma Position Cut 
	  double R[2];
	
	  bool bPosition = true;
	  if( TMath::Abs(y[0])> 550  ||
	      TMath::Abs(y[1])> 550  ){ 
	    bPosition = false;
	  }
	  for( int ig = 0; ig < 2; ig++){
	    R[ig] = TMath::Sqrt( x[ig]*x[ig] + y[ig]*y[ig]);
	    if( R[ig] > 850){ bPosition = false; }
	    if( TMath::Abs(y[ig]) < 150 && TMath::Abs(x[ig])< 150 ){
	      bPosition = false;
	    }	  
	  }
	  if( !bPosition ){ continue; }
	  double cosTheta = TMath::Abs( x[0]*x[1]+y[0]*y[1] )/TMath::Sqrt((x[0]*x[0]+y[0]*y[1])*(x[1]*x[1]+y[1]*y[1]));
	  double Eg1 = (*pit).g1().e();
	  double Eg2 = (*pit).g2().e();
	  double gchisq_1 = (*pit).g1().chisq();
	  double gchisq_2 = (*pit).g2().chisq();
	  double pi0pt    = TMath::Sqrt((*pit).p3()[0]*(*pit).p3()[0]+ (*pit).p3()[1]*(*pit).p3()[1]);
	  if( Eg1 > 350 &&
	      Eg2 > 150 &&
	      gchisq_1 < 5 && 
	    gchisq_2 < 5 &&
	      pi0pt  > 50  &&
	      cosTheta < 0.9 ){
	    hisPi0Mass[iFile]->Fill((*pit).m());
	  }
	}  
      }    
    }
  }
  for( int i = 0; i< 24; i++ ){
    hisPi0Mass[i]->Write();
  }
  tfout->Close();
  return 0;
}
