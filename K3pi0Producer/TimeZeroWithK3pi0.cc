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
#include "ClusterFinder_EDIT.h"
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

#include "E14WavReader.h"
#include "E14WavReader_V1.h"
#include "E14WaveReader_V2.h"


int
main( int argc ,char ** argv ){
  
  int RunNumber = atoi( argv[1]);
  std::string WAVFILE = std::getenv( "ROOTFILE_WAV");

  TChain* trin = new TChain("Tree"); 
  trin->Add(Form("%s/run_wav_%d.root",WAVFILE.c_str(),RunNumber));

  TFile* tfout = new TFile(Form("%s/run_wav_%d_cl.root",WAVFILE.c_str(),RunNumber),"recreate");
  TTree* trout = new TTree("TreeCalibration", "Output from Time zero" );
  int nCSIDigi = 0;
  int CSIDigiID[2716]={-1};
  double CSIDigiE[2716] = {0.};
  double CSIDigiTime[2716]  = {-1};
  double CSIDigiHHTime[2716]= {-1};

  trout->Branch("CsiNumber",&nCSIDigi,"CsiNumber/I");
  trout->Branch("CsiID",CSIDigiID,"CsiID[CsiNumber]/I");//nCSIDigi
  trout->Branch("CsiEne" ,CSIDigiE,"CsiEne[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiTime",CSIDigiTime,"CsiTime[CsiNumber]/D");//nCSIDigi
  trout->Branch("CsiHHTime",CSIDigiHHTime,"CsiHHTime[CsiNumber]/D");//nCSIDigi


  E14GNAnaDataContainer data; 

  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  data.branchOfKlong( trout );

  double TimeDelta[2716];
  double TimeDeltaSig[2716];
  int tmpID;
  double tmpDelta;
  double tmpDeltaSig;
  
  std::string ANAFILEDIR = std::getenv("HOME");
  std::ifstream ifs(Form("%s/local/Analysis/K3pi0Producer/Data/Pi0Peak.dat",ANAFILEDIR.c_str()));
  while( ifs >> tmpID >> tmpDelta >> tmpDeltaSig ){
    TimeDelta[ tmpID ]    = tmpDelta;
    TimeDeltaSig[ tmpID ] = tmpDeltaSig; 
  }

  E14WaveReader_V2* reader = new E14WaveReader_V2(trin);
  GammaFinder   gFinder;
  ClusterFinder_EDIT clusterFinder;

  int nCsI               = 0;
  int CsIID[2716]        = {-1};
  double CsIEnergy[2716] = {0.};
  double CsISignal[2716] = {0.};
  double CsITime[2716]   = {-1};
  double CsIHHTime[2716] = {-1};

  
  std::cout<< __PRETTY_FUNCTION__ << " : " << __LINE__ << std::endl;
  Long_t entries =  reader->fChain->GetEntries();
  for( int ievent  = 0; ievent < entries ; ievent++){
    nCsI = 0;
    nCSIDigi = 0;
    for( int ich = 0; ich < 3000; ich++){
      CsIID[ich]     = -1; 
      CsIEnergy[ich] = 0.;
      CsITime[ich]   = -1;
      CsIHHTime[ich] = -1;
      CsISignal[ich] = 0.;
      CSIDigiID[ich] = 0;

      CSIDigiE[ich] = 0;
      CSIDigiTime[ich]   = 0;
      CSIDigiHHTime[ich] = 0;

    }    
    if( (ievent%100) ==0 && ievent ){ std::cout<< ievent << std::endl;}
    reader->GetEntry( ievent  );

    if(( reader->TrigFlag & 3 ) != 0 ){ continue; }
    if( reader->CsinTimeCluster == 0 ){ continue; }
    if( reader->CsinTimeCluster  > 2  ){ continue ;}    
    if( reader->CsinTimeCluster == 2 && reader->CsiTimeClusterHead[0] < 50 ){ continue; }
    //// assumption :: All REAL csi Event if TimeCluster #0 
    std::cout<< "Analysis" << std::endl;

    nCsI = 0;
    for( int ich  = 0; ich < reader->CsiNumber; ich++){      
      int CsiID        = reader->CsiID[ich];
      double CsiTime   = reader->CsiTime[ich];
      double CsiSignal = reader->CsiSignal[ich]; 
      double CsiEnergy = reader->CsiEne[ich];      
      double CsiHHTime = reader->CsiHHTime[ich];
      int CsiTimeClusterID = reader->CsiTimeClusterID[ich];
      if( CsiTimeClusterID == 0){
	CsIID[nCsI]     =  CsiID;
	CsISignal[nCsI] =  CsiSignal;
	CsIEnergy[nCsI] =  CsiEnergy;
	CsITime[nCsI]   =  CsiTime;
	CsIHHTime[nCsI] =  CsiHHTime;
	nCsI++;
      }
    }
    /// Adjustment ///

    nCSIDigi=0;
    for( int idigi = 0; idigi< nCsI; idigi++){
      if( CsISignal[idigi] > 5 && CsIEnergy[idigi]>0.5){
	CSIDigiID[ nCSIDigi ]     = CsIID[idigi];
	CSIDigiE[ nCSIDigi ]      = CsIEnergy[idigi];
	CSIDigiTime[ nCSIDigi ]   = CsITime[idigi]-TimeDelta[CsIID[idigi]];
	CSIDigiHHTime[ nCSIDigi ] = CsIHHTime[idigi];
	nCSIDigi++;
      }
    }
    
    if( nCSIDigi<5){ continue;}

    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    
    
    std::cout<< "Cluster" << std::endl;
    std::cout << nCSIDigi << std::endl;
    for( int i = 0; i< nCSIDigi; i++ ){

      std::cout<< i            << " / "
	       << nCSIDigi     << " : " 
	       << CSIDigiID[i] << " : " 
	       << CSIDigiE[i]  << " : "
	       << CSIDigiTime[i] << std::endl;
    }

    clist = clusterFinder.findCluster( nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);
    std::cout<< "ClusterFinder" << std::endl;
    /*

    gFinder.findGamma( clist, glist );
    std::cout<< clist.size() << "\t" << glist.size() << std::endl;
    if( clist.size() < 6 ){ continue; }
    if( glist.size() !=6 ){ continue; }
    if( glist.size() == 6 ){

      data.setData( clist );
      data.setData( glist );
      user_cut( data, klVec );
      data.setData(klVec);

    
    }
  */
    trout->Fill();
  }
  std::cout<< "End" << std::endl;
  trout->Write();
  std::cout<< "Write" << std::endl;

  tfout->Close();
  return 0;
}
