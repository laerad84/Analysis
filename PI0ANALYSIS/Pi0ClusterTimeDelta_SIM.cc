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
//#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"

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

#include "L1TrigCounter.h"
#include "EnergyConverter.h"

int
main( int argc ,char ** argv ){
  int RunNumber = atoi( argv[1]);

  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");

  std::string ROOTFILE_DATACONV = std::getenv("ROOTFILE_WAV");
  std::string ROOTFILE_DATAPI0  = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/2012_FEB/Pi0_data";
  std::string ROOTFILE_SIMCONV  = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/ConvFile";
  std::string ROOTFILE_SIMPI0   = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/SIMPI0";

  //std::string iFileForm = "%s/Conv_e14_Al_Target_BS.mac_1000000_%d.root"; //ROOTFILE_SIMCONV
  //std::string oFileForm = "%s/SimPi0_1E6_Al_BS_%d.root";              //ROOTFILE_SIMPI0
  std::string iFileForm = "%s/Conv_e14_AL_Target_KL.mac_1000000_%d.root"; //ROOTFILE_SIMCONV
  std::string oFileForm = "%s/SimPi0_1E6_KLBEAM_%d.root";              //ROOTFILE_SIMPI0
  


  // In case of Data, Appling Temperature & Calibration is applied //
  // No Need to apply Temperature & Calibration constant //
  //std::string TCalFile = Form("%s/Data/TimeOffset/TimeOffset_with_cosmic.dat",ANALYSISLIB.c_str());  
  std::string TCalFile = Form("%s/Data/TimeOffset/testNewWORKCompileOffset.txt",ANALYSISLIB.c_str());  
  std::string ECalFile = Form("%s/local/Analysis/K3pi0Producer/Data/CalibrationFactorADV_15.dat",HOME.c_str());
  std::string TempCalibrationFilename = Form("%s/Data/Temperature_Factor/TemperatureCorrectionFactor.root",ANALYSISLIB.c_str());

  EnergyConverter* Converter = new EnergyConverter();
  Converter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",ANALYSISLIB.c_str()));
  Double_t Pi0PeakCorFactor = 0.9937;

  //TChain* trin = new TChain("T"); 
  //trin->Add(Form(iFileForm.c_str(),ROOTFILE_SIMCONV.c_str(),RunNumber));
  TFile* tfin = new TFile(Form(iFileForm.c_str(),ROOTFILE_SIMCONV.c_str(),RunNumber));
  TTree* trin = (TTree*)tfin->Get("T");

  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_SIMPI0.c_str(),RunNumber),"recreate");
  TTree* trout = new TTree("Tree", "Output from Time zero" );

  int EventNumber;
  int CsiNumber;
  int CsiModID[2716];
  double CsiEne[2716];
  double CsiTime[2716];
  trin->SetBranchAddress("EventNum",&EventNumber);
  trin->SetBranchAddress("CsiNumber",&CsiNumber);
  trin->SetBranchAddress("CsiModID",CsiModID);
  trin->SetBranchAddress("CsiEne",CsiEne);
  trin->SetBranchAddress("CsiTime",CsiTime);
  
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

  trout->Branch("CVNumber" ,&CVNumber ,"CVNumber/I");
  trout->Branch("CVModID"  ,CVModID   ,"CVModID[CVNumber]/D");//CVNumber
  trout->Branch("CVEne"    ,CVEne     ,"CVEne[CVNumber]/D");//CVNumber
  trout->Branch("CVTime"   ,CVTime    ,"CVTime[CVNumber]/D");//CVNumber
  trout->Branch("SciNumber",&SciNumber,"SciNumber/I");
  trout->Branch("SciEne"   ,SciEne    ,"SciEne[SciNumber]/D");//SciNumber
  

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

  /*
  int    NTrack;
  short  Track[200];
  int    Pid[200];
  float  Mass[200];
  float  Ek[200];
  float  end_Ek[200];
  double P[200];
  double end_P[200];
  double V[200];
  double end_V[200];
  
  trin->Branch("nTrack",&NTrack,"nTrack/I");
  trin->Branch("track",Track,"track[nTrack]/S");//NTrack
  trin->Branch("pid",Pid,"pid[nTrack]/F");//NTrack
  trin->Branch("mass",Mass,"mass[nTrack]/F");//NTrack
  trin->Branch("ek",Ek,"ek[nTrack]/F");//NTrack
  trin->Branch("end_ek",end_Ek,"end_ek[nTrack]/F");//NTrack
  trin->Branch("p",P,"p[nTrack]/D");//NTrack
  trin->Branch("end_p",end_P,"end_p[nTrack]/D");//NTrack
  trin->Branch("v",V,"v[nTrack]/D");//NTrack
  trin->Branch("end_v",end_V,"end_v[nTrack]/D");//NTrack
  */

  trout->Branch("nTrack",&nTrack,"nTrack/I");
  trout->Branch("track",track,"track[nTrack]/S");//nTrack
  trout->Branch("pid",pid,"pid[nTrack]/I");//nTrack
  trout->Branch("mass",mass,"mass[nTrack]/F");//nTrack
  trout->Branch("ek",ek,"ek[nTrack]/F");//nTrack
  trout->Branch("end_ek",end_ek,"end_ek[nTrack]/F");//nTrack
  trout->Branch("p",p,"p[nTrack][3]/D");//nTrack
  trout->Branch("end_p",end_p,"end_p[nTrack][3]/D");//nTrack
  trout->Branch("v",v,"v[nTrack][3]/D");//nTrack
  trout->Branch("end_v",end_v,"end_v[nTrack][3]/D");//nTrack

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  
  double TimeDelta[2716];
  double TimeDeltaSig[2716];
  double CalibrationFactor[2716] = {1.};
  double TimeDeltaLength[2716]={0};
  int tmpID;
  double tmpDelta;
  double tmpDeltaSig;
  double tmpCalFactor; 
  for( int i = 0; i< 2716; i++){
    CalibrationFactor[i] = 1;
  }
  
  //std::ifstream ifs(Form("%s/local/Analysis/K3pi0Producer/Data/Pi0Peak.dat",HOME.c_str()));
  std::ifstream ifsTCal(Form(TCalFile.c_str(),ANALYSISLIB.c_str()));
  if( !ifsTCal.is_open() ) { std::cerr <<"File does not exist."<< Form(TCalFile.c_str(),ANALYSISLIB.c_str())  << std::endl; return -1;}
  while( ifsTCal >> tmpID >> tmpDelta >> tmpDeltaSig ){
    TimeDelta[ tmpID ]    = tmpDelta;
    TimeDeltaSig[ tmpID ] = tmpDeltaSig; 
  }

  std::ifstream ifsECal(Form(ECalFile.c_str(),HOME.c_str()));
  if( !ifsECal.is_open() ){ std::cerr << "File does not exist." << Form(ECalFile.c_str(),HOME.c_str()) << std::endl; return -1; }
  while( ifsECal >> tmpID >> tmpCalFactor ){
    CalibrationFactor[ tmpID ] = tmpCalFactor;
  }
  
  L1TrigCounter* l1 = new L1TrigCounter();
  l1->ReadMapFile();
  l1->SetThreshold(1000);
  l1->Reset();
  
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
    l1->Reset();    
    nCSIDigi = 0;
    //std::cout<< __PRETTY_FUNCTION__ << std::endl;
    for( int ich = 0; ich < CsiNumber; ich++){
      int tmpID = CsiModID[ich];
      double tmpTime   = CsiTime[ich];
      double tmpHHTime = CsiTime[ich]-50;// Just Set //       
      double tmpSignal = Converter->ConvertToHeight(CsiModID[ich],CsiEne[ich])/Pi0PeakCorFactor/CalibrationFactor[CsiModID[ich]];
      double tmpEne = CsiEne[ich];
      if( tmpSignal > 5 && tmpEne > 0.5 ){
	CSIDigiID[nCSIDigi]     = tmpID;
	//CSIDigiE[nCSIDigi]      = tmpEne/CalibrationFactor[CsiModID[ich]];
	//CSIDigiE[nCSIDigi]      = tmpEne/CalibrationFactor[CsiModID[ich]]*(1+0.02*(tmpSignal/10000)*(tmpSignal/10000));// For Calibration Reconstruction Test.
	CSIDigiE[nCSIDigi]      = tmpEne;// /CalibrationFactor[CsiModID[ich]]*(1+0.02*(tmpSignal/10000)*(tmpSignal/10000))*(1.013-9.953e-3*exp(tmpSignal/1.5314e4));// First Adjustment
	CSIDigiSignal[nCSIDigi] = tmpSignal;
	CSIDigiTime[nCSIDigi]   = tmpTime;
	CSIDigiHHTime[nCSIDigi] = tmpHHTime;
	nCSIDigi++;
	if(tmpSignal == TMath::Infinity()){ 
	  std::cout<< tmpID << "\t" 
		   << tmpEne << "\t"
		   << tmpSignal << "\t"
		   << Converter->GetCalibrationConstant(tmpID) << "\t"	    
		   << CalibrationFactor[CsiModID[ich]] << "\t"
		   << std::endl;
	}
	l1->Fill(tmpID,tmpSignal);
      }
    }
    

    CSIL1nTrig = 0;
    std::vector<double> vecCount    = l1->GetCount(); 
    for( int i = 0; i< vecCount.size(); i++){
      CSIL1TrigCount[i] = vecCount.at(i);
      if( vecCount.at(i) > 1500 ){ 
	CSIL1nTrig++;
      }
    }

    trout->Fill();
  }
  std::cout<< "End" << std::endl;
  tfout->cd();
  trout->Write();
  std::cout<< "Write" << std::endl;
  tfout->Close();
  return 0;
}
