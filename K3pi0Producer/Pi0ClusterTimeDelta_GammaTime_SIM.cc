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
#include "TF1.h"
#include "IDHandler.h"
#include "TRandom.h"

#include "E14WavReader_V1.h"
//#include "L1TrigCounter.h"
#include "EnergyConverter.h"

#include "User_Functions.h"
#include "GeneralFunctions.h"





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
  
  std::string ROOTFILE_SIMCONV  = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/ConvFile";
  std::string ROOTFILE_SIMPI0   = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/SIMPI0";

  //std::string iFileForm          = "%s/SimPi0_1E6_Al_BS_%d.root";        // ROOTFILE_SIMCONV
  //std::string oFileForm          = "%s/SimPi0_1E6_LYRES_Al_BS_%d.root"; // ROOTFILE_SIM3PI0

  std::string iFileForm          = "%s/SimPi0_1E6_KLBEAM_%d.root";        // ROOTFILE_SIMCONV
  std::string oFileForm          = "%s/SimPi0_1E6_LYRES_KL_%d.root"; // ROOTFILE_SIM3PI0

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
    trin->Add(Form(iFileForm.c_str(),ROOTFILE_SIMPI0.c_str(),i));
  }

  int    RunNo;
  int    EventNumber;
  int    CsiNumber;
  int    CsiModID[2716];
  double CsiEne[2716];
  double CsiTime[2716];
  double CsiHHTime[2716];
  double CsiSignal[2716];
  double CsiChisq[2716];
  short  CsiNDF[2716];
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
  trin->SetBranchAddress("CsiChisq",CsiChisq);
  trin->SetBranchAddress("CsiNDF",CsiNDF);
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


  int     CVNumber      = 0;
  Short_t CVModID[100]  = {-1};
  double  CVEne[100] = {0};
  double  CVTime[100]   = {0};
  int     SciNumber = 0;
  double  SciEne[1] ={0};

  trin->SetBranchAddress("CVNumber" ,&CVNumber);
  trin->SetBranchAddress("CVModID"  ,CVModID  );//CVNumber
  trin->SetBranchAddress("CVEne"    ,CVEne );//CVNumber
  trin->SetBranchAddress("CVTime"   ,CVTime   );//CVNumber
  trin->SetBranchAddress("SciNumber",&SciNumber);
  trin->SetBranchAddress("SciEne"   ,SciEne);//SciNumber

  
  double CsiL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
					1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  GammaCut* gammaCut = new GammaCut();
  CsiCut*   csiCut   = new CsiCut();

  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_SIMPI0.c_str(),RunNumber),"recreate");
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
  csiCut->Branch(trout);
  gammaCut->Branch(trout);

  /*
  trout->Branch("CsiNumber"     ,&nCSIDigi     ,"CsiNumber/I");
  trout->Branch("CsiModID"      ,CSIDigiID     ,"CsiModID[CsiNumber]/I");//nCSIDigi
  trout->Branch("CsiEne"        ,CSIDigiE      ,"CsiEne[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiTime"       ,CSIDigiTime   ,"CsiTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiHHTime"     ,CSIDigiHHTime ,"CsiHHTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiSignal"     ,CSIDigiSignal ,"CsiSignal[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiL1nTrig"    ,&CSIL1nTrig   ,"CsiL1nTrig/I");
  trout->Branch("CsiL1TrigCount",CSIL1TrigCount,"CsiL1TrigCount[20]/D");
  */


  trout->Branch("nTrack",&nTrack,"nTrack/I"          );
  trout->Branch("track" ,track  ,"track[nTrack]/S"   );//nTrack
  trout->Branch("pid"   ,pid    ,"pid[nTrack]/I"     );//nTrack
  trout->Branch("mass"  ,mass   ,"mass[nTrack]/F"    );//nTrack
  trout->Branch("ek"    ,ek     ,"ek[nTrack]/F"      );//nTrack
  trout->Branch("end_ek",end_ek ,"end_ek[nTrack]/F"  );//nTrack
  trout->Branch("p"     ,p      ,"p[nTrack][3]/D"    );//nTrack
  trout->Branch("end_p" ,end_p  ,"end_p[nTrack][3]/D");//nTrack
  trout->Branch("v"     ,v      ,"v[nTrack][3]/D"    );//nTrack
  trout->Branch("end_v" ,end_v  ,"end_v[nTrack][3]/D");//nTrack

  trout->Branch("CVNumber" ,&CVNumber ,"CVNumber/I");
  trout->Branch("CVModID"  ,CVModID   ,"CVModID[CVNumber]/D");//CVNumber
  trout->Branch("CVEne"    ,CVEne     ,"CVEne[CVNumber]/D");//CVNumber
  trout->Branch("CVTime"   ,CVTime    ,"CVTime[CVNumber]/D");//CVNumber
  trout->Branch("SciNumber",&SciNumber,"SciNumber/I");
  trout->Branch("SciEne"   ,SciEne    ,"SciEne[SciNumber]/D");//SciNumber
  

  E14GNAnaDataContainer data; 
  data.branchOfClusterList(trout);
  data.branchOfPi0List(trout);
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  GammaFinder gFinder;
  ClusterFinder cFinder;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int nCsI               = 0;
  int CsIID[2716]        = {-1};
  double CsIEnergy[2716] = {0.};
  double CsISignal[2716] = {0.};
  double CsITime[2716]   = {-1};
  double CsIHHTime[2716] = {-1};
  

  // Alz : 2624 from CsI
  double AlzPosition = 3526;//mm

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
    /////// Initialize data /////////

    if( CsiNumber<4){ continue;}

    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::list<Gamma>   glistTCut;
    std::list<Pi0>     plist;

    csiCut->Decision( CsiNumber, CsiModID, CsiEne,CsiTime, CsiSignal, CsiChisq,CsiNDF);
    //clist = cFinder.findCluster( CsiNumber,CSIDigiID,CSIDigiE,CSIDigiTime);
    clist = cFinder.findCluster( csiCut->CsiNumber,csiCut->CsiID,csiCut->CsiEne,csiCut->CsiTime);
    gFinder.findGamma( clist, glist );
    gammaCut->Decision( glist );
    SetGammaTime( glist );


    //////////////////////    //////////////////////
    //////////////////////    //////////////////////
    //////////////////////    //////////////////////
    //////////////////////    //////////////////////



    if( clist.size() < 2 ){ continue; }
    if( glist.size() < 2 ){ continue; }
    
    std::list<Gamma>::iterator git = glist.begin();
    for( ; git != glist.end(); git++){
      SetGammaTime( (*git));
    }
    GammaTimeDeltaCut( glist, glistTCut,2);
    data.setData( clist );
    data.setData( glistTCut );
    std::list<Gamma>::iterator gitT = glistTCut.begin();
    if( glistTCut.size() ==2 ){
      if( User_RecG2(glist,plist)){
	std::list<Pi0>::iterator it = plist.begin();
	(*it).setVtx( 0,0,AlzPosition );
	(*it).updateVars();
	data.setData( clist );
	data.setData( glist );
	data.setData(plist);    
	trout->Fill();
      }
    }
  }

  std::cout<< "End" << std::endl;
  trout->Write();
  std::cout<< "Write" << std::endl;
  tfout->Close();
  return 0;
}
