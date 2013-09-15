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
  //std::string iFileForm          = "%s/Sim3pi0_wav_fast_5E6_%d_Calibration.root";    //ROOTFILE_SIM3PI0
  //std::string oFileForm          = "%s/Sim3pi0_wav_fast_KL_RES_LY_pe_5E6_%d_Calibration.root"; // ROOTFILE_SIM3PI0
  //std::string iFileForm          = "%s/Sim3pi0_wav_fast_5E6_%d_Calibration_mis_1.root";    //ROOTFILE_SIM3PI0
  //std::string oFileForm          = "%s/Sim3pi0_wav_fast_KL_RES_LY_pe_5E6_%d_Calibration_mis_1.root"; // ROOTFILE_SIM3PI0

  std::string iFileForm          = "%s/Sim3pi0_wav_ALCV_5E8_%d.root";    //ROOTFILE_SIM3PI0
  std::string oFileForm          = "%s/Sim3pi0_wav_ALCV_KL_RES_LY_pe_5E8KL_%d.root"; // ROOTFILE_SIM3PI0


  //std::string iFileForm = "%s/Sim3pi0_wav_fast_5E6_%d_Calibration.root";
  //std::string oFileForm = "%s/Sim3pi0_wav_fast_KL_RES_LY_pe_5E6_%d_Calibration_%s.root";


  TF1* func = new TF1("ResFunc", funcResolutionInvSq, 0, 10000,1);
  /*
  EnergyConverter* Converter = new EnergyConverter();
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",ANALYSISLIB.c_str()));
  double CalFactor0 = 0.08485;// 14MeV/165Cnt;  
  */
  double RelativeLY[2716] = {0};
  std::ifstream ifs(Form("%s/local/Analysis/K3pi0Producer/Data/RelativeLY.txt",HOME.c_str()));
  int listID;
  double listLY;
  while( ifs >> listID >> listLY ){
    RelativeLY[listID] = listLY;
    //std::cout<< listLY << std::endl;
  }

  /*
  TFile* tfin = new TFile(Form(iFileForm.c_str(),ROOTFILE_SIM3PI0.c_str(),RunNumber));
  TTree* trin = (TTree*)tfin->Get("Tree");
  */

  TChain* trin = new TChain("Tree");
  for( int i = RunNumber*100; i < (RunNumber+1)*100; i++){
    trin->Add(Form(iFileForm.c_str(),ROOTFILE_SIM3PI0.c_str(),i));
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
  
  trin->SetBranchAddress("RunNumber",&RunNo);
  trin->SetBranchAddress("EventNumber",&EventNumber);
  trin->SetBranchAddress("CsiNumber",&CsiNumber);
  trin->SetBranchAddress("CsiModID",CsiModID);
  trin->SetBranchAddress("CsiEne",CsiEne);
  trin->SetBranchAddress("CsiTime",CsiTime);
  trin->SetBranchAddress("CsiHHTime",CsiHHTime);
  trin->SetBranchAddress("CsiSignal",CsiSignal);
  trin->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  trin->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);


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

  trout->Branch("RunNumber"     ,&RunNo        ,"RunNumber/I");
  trout->Branch("EventNumber"   ,&EventNumber  ,"EventNumber/I");
  trout->Branch("CsiNumber"     ,&nCSIDigi     ,"CsiNumber/I");
  trout->Branch("CsiModID"      ,CSIDigiID     ,"CsiModID[CsiNumber]/I");//nCSIDigi
  trout->Branch("CsiEne"        ,CSIDigiE      ,"CsiEne[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiTime"       ,CSIDigiTime   ,"CsiTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiHHTime"     ,CSIDigiHHTime ,"CsiHHTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiSignal"     ,CSIDigiSignal ,"CsiSignal[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiL1nTrig"    ,&CSIL1nTrig   ,"CsiL1nTrig/I");
  trout->Branch("CsiL1TrigCount",CSIL1TrigCount,"CsiL1TrigCount[20]/D");


  trout->Branch("nTrack",&nTrack,"nTrack/I");
  trout->Branch("track" ,track  ,"track[nTrack]/S"   );//nTrack
  trout->Branch("pid"   ,pid    ,"pid[nTrack]/I"     );//nTrack
  trout->Branch("mass"  ,mass   ,"mass[nTrack]/F"    );//nTrack
  trout->Branch("ek"    ,ek     ,"ek[nTrack]/F"      );//nTrack
  trout->Branch("end_ek",end_ek ,"end_ek[nTrack]/F"  );//nTrack
  trout->Branch("p"     ,p      ,"p[nTrack][3]/D"    );//nTrack
  trout->Branch("end_p" ,end_p  ,"end_p[nTrack][3]/D");//nTrack
  trout->Branch("v"     ,v      ,"v[nTrack][3]/D"    );//nTrack
  trout->Branch("end_v" ,end_v  ,"end_v[nTrack][3]/D");//nTrack


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
	func->SetParameter(0,RelativeLY[tmpID]);
	double value = func->Eval(tmpEne);
	double mont  = gRandom->PoissonD(value)/value;
	
	CSIDigiID[nCSIDigi]     = tmpID;
	CSIDigiE[nCSIDigi]      = tmpEne*mont;
	CSIDigiSignal[nCSIDigi] = tmpSignal*mont;
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
    if( clist.size() < 6 ){ continue; }
    if( glist.size() == 6 ){
      if( user_rec(glist,klVec)){
	data.setData(clist);
	data.setData(glist);
	user_cut( data, klVec);
	data.setData(klVec);
      }
      trout->Fill();
    }
  }

  std::cout<< "End" << std::endl;
  trout->Write();
  std::cout<< "Write" << std::endl;
  tfout->Close();
  return 0;
}
