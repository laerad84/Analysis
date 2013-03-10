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
#include "klong/Klong.h"
#include "cluster/ClusterFinder.h"
//#include "ClusterFinder_EDIT.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "User_Function.h"

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
#include "L1TrigCounter.h"
#include "EnergyConverter.h"

int
main( int argc ,char ** argv ){
  
  int RunNumber = atoi( argv[1]);
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");
  
  //std::string iFileForm="%s/run_wav_%d.root";
  //std::string oFileForm="%s/run_wav_%d_Cal_FNL_COS_newTimeOffset.root";
  std::string ROOTFILE_SIMCONV = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/3pi0Run/ConvFile";
  std::string ROOTFILE_SIM3PI0 = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/3pi0Run/SIM3PI0";

  //std::string iFileForm        = "%s/Conv_KL3pi0.mac_1000000_%d.root"; //ROOTFILE_SIMCONV
  //std::string oFileForm        = "%s/Sim3pi0_wav_%d.root";     //ROOTFILE_SIM3PI0
  std::string iFileForm          = "%s/Sim3pi0_wav_%d.root";    //ROOTFILE_SIM3PI0
  std::string oFileForm          = "%s/Sim3pi0_Wav_KL_%d.root"; // ROOTFILE_SIM3PI0

  TFile* tfin = new TFile(Form(iFileForm.c_str(),ROOTFILE_SIMCONV.c_str(),RunNumber));
  TTree* trin = (TTree*)tfin->Get("Tree");

  int    EventNumber;
  int    CsiNumber;
  int    CsiModID[2716];
  double CsiEne[2716];
  double CsiTime[2716];
  double CsiHHTime[2716];
  double CsiSignal[2716];
  int    CsiL1nTrig;
  double CsiL1TrigCount[20];
  
  trin->SetBranchAddress("EventNum",&EventNumber);
  trin->SetBranchAddress("CsiNumber",&CsiNumber);
  trin->SetBranchAddress("CsiModID",CsiModID);
  trin->SetBranchAddress("CsiEne",CsiEne);
  trin->SetBranchAddress("CsiTime",CsiTime);
  trin->SetBranchAddress("CsiHHTime",CsiHHTime);
  trin->SetBranchAddress("CsiSignal",CsiSignal);
  trin->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  trin->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);

  double CsiL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
					 1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_SIM3PI0.c_str(),RunNumber),"recreate");
  TTree* trout = new TTree("T", "Output from Time zero" );  
  
  int nCSIDigi = 0;
  int CSIDigiID[2716]={-1};
  double CSIDigiE[2716] = {0.};
  double CSIDigiTime[2716] ;
  double CSIDigiHHTime[2716];
  double CSIDigiSignal[2716]={0};
  double CSIL1TrigCount[20];
  int    CSIL1nTrig;

  trout->Branch("RunNumber"     ,&RunNumber    ,"RunNumber/I");
  trout->Branch("EventNumber"   ,&EventNumber  ,"EventNumber/I");
  trout->Branch("CsiNumber"     ,&nCSIDigi     ,"CsiNumber/I");
  trout->Branch("CsiModID"      ,CSIDigiID     ,"CsiModID[CsiNumber]/I");//nCSIDigi
  trout->Branch("CsiEne"        ,CSIDigiE      ,"CsiEne[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiTime"       ,CSIDigiTime   ,"CsiTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiHHTime"     ,CSIDigiHHTime ,"CsiHHTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiSignal"     ,CSIDigiSignal ,"CsiSignal[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiL1nTrig"    ,&CSIL1nTrig   ,"CsiL1nTrig/I");
  trout->Branch("CsiL1TrigCount",CSIL1TrigCount,"CsiL1TrigCount[20]/D");

  E14GNAnaDataContainer data; 
  data.branchOfClusterList(trout);
  data.branchOfDigi(trout);
  data.branchOfKlong(trout);
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  GammaFinder gFinder;
  ClusterFinder clusterFinder;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int nCsI               = 0;
  int CsIID[2716]        = {-1};
  double CsIEnergy[2716] = {0.};
  double CsISignal[2716] = {0.};
  double CsITime[2716]   = {-1};
  double CsIHHTime[2716] = {-1};
  
  std::cout<< __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  std::cout<< "Start Loop" << std::endl;

  //Long_t entries =  reader->fChain->GetEntries();
  long entries = trin->GetEntries();
  std::cout<< entries << std::endl;
  for( int ievent  = 0; ievent < entries ; ievent++){
    //for( int ievent  = 0; ievent < 100 ; ievent++){
    trin->GetEntry(ievent);
    data.reset();

    if( (ievent%1000) ==0 && ievent ){ std::cout<< ievent << std::endl;}
    //std::cout<< ievent << std::endl;
    /////// Initialize data /////////
    for( int ich = 0; ich < 2716; ich++){
      CsIID[ich]     = -1; 
      CsIEnergy[ich] = 0.;
      CsITime[ich]   = -1;
      CsIHHTime[ich] = -1;
      CsISignal[ich] = 0.;

      CSIDigiID[ich] = -1;
      CSIDigiE[ich] = 0;
      CSIDigiTime[ich]   = 0;
      CSIDigiHHTime[ich] = 0;
      CSIDigiSignal[ich] = 0.;
    }       
    nCsI = 0; 
    nCSIDigi = 0;
    CSIL1nTrig = 0; 
    nCSIDigi = 0;
    //std::cout<< __PRETTY_FUNCTION__ << std::endl;
    for( int ich = 0; ich < CsiNumber; ich++){
      int tmpID        = CsiModID[ich];
      double tmpTime   = CsiTime[ich];
      double tmpHHTime = CsiHHTime[ich]-50;// Just Set //       
      double tmpSignal = CsiSignal[ich];
      double tmpEne    = CsiEne[ich];
      if( tmpSignal > 5 && tmpEne > 0.5 ){
	CSIDigiID[nCSIDigi]     = tmpID;
	CSIDigiE[nCSIDigi]      = tmpEne;
	CSIDigiSignal[nCSIDigi] = tmpSignal;
	CSIDigiTime[nCSIDigi]   = tmpTime;
	CSIDigiHHTime[nCSIDigi] = tmpHHTime;
	nCSIDigi++;
      }
    }
    for( int i = 0; i< 20; i++){
      CSIL1TrigCount[i] = CsiL1TrigCount[i]; 
      if( CSIL1TrigCount[i] > CsiL1TrigCountThreshold[i] ){
	CSIL1nTrig++;
      }
    }
    if( nCSIDigi < 5 ){ continue;}
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    clist = clusterFinder.findCluster( nCSIDigi, CSIDigiID, CSIDigiE,CSIDigiTime);
    gFinder.findGamma(clist,glist);
    if( glist.size() == 6 ){
      if( user_rec(glist,klVec)){
	data.setData(clist);
	data.setData(glist);
	user_cut( data, klVec);
	data.setData(klVec);
      }
    }
    trout->Fill();
  }

  std::cout<< "End" << std::endl;
  trout->Write();
  std::cout<< "Write" << std::endl;
  tfout->Close();
  return 0;
}
