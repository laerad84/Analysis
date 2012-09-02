
// The Propose of this program is making Templetes of Trigger Signal
// For example, Cosmic, Laser and CV


//#include "gnana/DigiReader.h"
//#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TTreePlayer.h"
#include "TTreePerfStats.h"
#include "TChain.h"

#include "GeneralTypes.h"
#include "GeneralMacros.h"
#include "WaveformFitter.h"
#include "E14MapReader.h"
#include "Structs.h"

#include "E14ConvReader.h"
#include "E14IDHandler.h"
#include "E14ConvWriter.h"
#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TText.h"
#include "TROOT.h"
#include "Math/WrappedFunction.h"
#include "Math/BrentRootFinder.h"
#include "Math/BrentMethods.h"
#include "Math/WrappedTF1.h"

#include "E14WaveFitter.h"
#include "CsI_Module.h"
#include "IDHandler.h"


#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/GammaFinder.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "user_func.h"



const int nCrate = 11;

const Double_t COSMIC_THRESHOLD[20] = {100,100,100,100,100,
				       100,100,100,100,100,
				       100,100,100,100,100,
				       100,100,100,100,100};


int  main(int argc,char** argv)
{
  if( argc !=2  && argc != 3  ){
    std::cerr << "Please Input RunNumber " << std::endl;
    return -1; 
  }
  Int_t RunNumber = atoi( argv[1] );

  Int_t iDiv;
  if( argc ==3 ){
    iDiv = atoi( argv[2] );
  }

  // GetEnvironment // 
  std::string ANALIBDIR   = std::getenv("ANALYSISLIB"  );
  std::string CONVFILEDIR = std::getenv("ROOTFILE_CONV");
  std::string WAVEFILEDIR = std::getenv("ROOTFILE_WAV" );
  std::string SUMFILEDIR  = std::getenv("ROOTFILE_SUMUP");
  std::string COSMICFILEDIR = std::getenv("ROOTFILE_COSMIC");
  std::cout << ANALIBDIR   << std::endl;  
  std::cout << CONVFILEDIR << std::endl;
  std::cout << WAVEFILEDIR << std::endl;
  std::cout << SUMFILEDIR  << std::endl;
  std::cout << COSMICFILEDIR << std::endl;
  // Setting  Classes //
  WaveformFitter* wavFitter = new WaveformFitter(48, kFALSE);  
  E14WaveFitter* Fitter     = new E14WaveFitter();
  TApplication* app = new TApplication("app", &argc , argv );  
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  //  Setting Calibration Factor
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  Double_t CosmicPeak[3000]={0};  
  Double_t CosmicCalibrationFactor[3000]={0};
  Double_t CosmicOut;
  Int_t    CosmicID;
  TFile*   tfCosmic = new TFile(Form("%s/CosmicResult_20120209.root",COSMICFILEDIR.c_str()));
  TTree*   trCosmic = (TTree*)tfCosmic->Get("GainFitPar");
  trCosmic->SetBranchAddress("ID",&CosmicID);
  trCosmic->SetBranchAddress("Peak",&CosmicOut);
  for( int i = 0; i< trCosmic->GetEntries(); i++ ){
    trCosmic->GetEntry(i);    
    CosmicPeak[CosmicID] = CosmicOut;
    if(CosmicOut != 0){
      CosmicCalibrationFactor[CosmicID] = 14.014/CosmicOut;
    }else{
      CosmicCalibrationFactor[CosmicID] = 0;
    }
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read TimeOffset // 
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  std::ifstream  ifsTimeOffset("testNewWORKCompileOffset.txt");
  Double_t TimeOffset[2716]={0};
  Int_t    tempID;
  Double_t tempOffset;
  Double_t tempChisq;
  while(  ifsTimeOffset >> tempID >> tempOffset >> tempChisq ){
    TimeOffset[tempID] = tempOffset;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read ID map // 
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  IDHandler* idHandler = new IDHandler();

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Set input / output File 
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout<< "Set I/O File" <<std::endl;  

  TFile* tfin;
  if( argc == 2){
    tfin = new TFile(Form("%s/TEMPLATE_FIT_RESULT_1_%d.root",WAVEFILEDIR.c_str(),RunNumber));
  }else if( argc == 3){
    tfin = new TFile(Form("%s/TEMPLATE_FIT_RESULT_1_%d_%d.root",WAVEFILEDIR.c_str(),RunNumber,iDiv));
  }

  TTree* trin  = (TTree*)tfin->Get("WFTree");

  std::cout << Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber) << std::endl;  
  E14ConvWriter* wConv = new E14ConvWriter( Form("%s/Sum%d.root",SUMFILEDIR.c_str(),RunNumber),
					    trin);

  std::cout<< "Setting Map" << std::endl;
  {
    wConv->AddModule("Csi");
    wConv->AddModule("CC03");
    wConv->AddModule("OEV");
    wConv->AddModule("CV");
    wConv->AddModule("Cosmic");
    wConv->AddModule("Laser");
    wConv->AddModule("Etc");
    wConv->Set();
    wConv->SetMap();
    wConv->SetBranchAddress();
  }
  std::cout<< " Prepare OutputFile " << std::endl;

  TFile* tfout;
  if( argc == 2){
    tfout = new TFile(Form("%s/RECONSTRUCT_RESULT_%d.root",SUMFILEDIR.c_str(),RunNumber),"recreate");
  }else if( argc == 3){
    tfout = new TFile(Form("%s/RECONSTRUCT_RESULT_%d_%d.root",SUMFILEDIR.c_str(),RunNumber,iDiv),"recreate");
  }
  TTree* trout  =new TTree("Tree","Reconstruct CsiEvent");
  std::cout << "Setting IO File End" << std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  /// Set I/O Class 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  E14GNAnaDataContainer data;
  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  data.branchOfPi0List( trout );
  ClusterFinder clusterFinder;
  GammaFinder   gFinder;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////  
  int nCSIDigi = 0;
  int CSIDigiID[3000] = {0};
  double CSIDigiE[3000]={0}, CSIDigiTime[3000]={0};
  double CSICalFactor[3000]={0};
  for( int ich = 0; ich < 3000; ++ich){
    CSICalFactor[ich] = 1;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
  std::cout<< "Pi0 Calibration" << std::endl;
  // loop analysis
  int nloop = trin->GetEntries();
  int iCsiMod = 0;
  TH1D* hist = new TH1D("his","",1000,0,2000);  
  std::cout<<"start loop analysis"<<std::endl;
  std::cout<<"# of entry : "<<nloop<<std::endl;
  for( int ievt=0; ievt<nloop; ievt++ ){

    if(ievt%(nloop/10)==0&&nloop>100)
      std::cout<<ievt/(nloop/100)<<"%"<<std::endl;
    
    for( int ichannel = 0; ichannel < 3000; ichannel++){
      CSIDigiID[ichannel]  = 0; 
      CSIDigiE[ichannel]   = 0; 
      CSIDigiTime[ichannel]= 0;
    }
    trin->GetEntry( ievt );
    double timePeak = wConv->m_TimePeak;
    double timeSigma= wConv->m_TimeSigma;
    if( timeSigma >=3.5 ) { continue; }
    nCSIDigi = 0;
    for( int ich = 0; ich < wConv->mod[iCsiMod]->m_nDigi;ich++){
      int CsiChannelID = wConv->mod[iCsiMod]->m_ID[ich];
      double Energy  = wConv->mod[iCsiMod]->m_Signal[ich] * CosmicCalibrationFactor[CsiChannelID];
      double TimeData= wConv->mod[iCsiMod]->m_Time[ich]-TimeOffset[ich];
      
      if( Energy > 3 && abs(TimeData - timePeak) < 8 ){
	CSIDigiID[nCSIDigi]   = CsiChannelID;
	CSIDigiE[nCSIDigi]    = Energy;
	CSIDigiTime[nCSIDigi] = TimeData;
	nCSIDigi++;
      }
    }

    std::cout<< nCSIDigi << std::endl;
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    clist = clusterFinder.findCluster(nCSIDigi, CSIDigiID, CSIDigiE, CSIDigiTime);
    for( std::list<Cluster>::iterator it  = clist.begin();
	 it != clist.end();
	 ++it){
      std::cout<< (*it).e() << std::endl;
    }
    gFinder.findGamma(clist,glist);

    // pi0 reconstrunction and position correction for angle dependency
    if( glist.size() != 2 ) continue; 
    std::list<Pi0> piList; 
    double mass=0;

    //double position =3484.; //shiomisan:3534.;
    double position = 3526;
    if(!user_rec(glist,piList,mass,position)) continue;
    std::cout<< mass << std::endl;
    hist->Fill(mass);
    // cuts
    //user_cut(data,piList);    
    
    // fill data to TTree
    data.setData( piList );
    trout->Fill();
    data.eventID++;

  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  trout->Write();
  tfout->Close();
  
  std::cout<< "Close" << std::endl;
  return 0;

}
