#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <list>

#include "EventTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TH1.h"
#include "TApplication.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMarker.h"
#include "TProfile.h"
#include "TRandom.h"

#include "CsIPoly.h"
#include "IDHandler.h"
#include "PulseGenerator.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include "gamma/Gamma.h"
#include "cluster/Cluster.h"

#include "gamma/GammaFinder.h"
#include "cluster/ClusterFinder.h"

#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
//#include "User_Function.h"
#include "ClusterTimeStructure.h"

int main( int argc, char** argv ){

  Int_t GammaEnergy = atoi( argv[1] );
  Int_t Degree      = atoi( argv[2] );
  Int_t DataIndex   = atoi( argv[3] );

  std::cout<< __LINE__ <<std::endl;
  TChain* ch = new TChain("T");  
  ch->Add(Form("/Volume0/gamma/Cluster_%dMeV_%ddeg-1E5-%d.root",GammaEnergy, Degree, DataIndex ));
  std::cout<< "File is Opened" << std::endl;

  E14GNAnaDataContainer data;
  ClusterTimeAnalyzer* clusterTimeAnalyzer = new ClusterTimeAnalyzer();

  data.setBranchAddress( ch );
  std::cout<< ch->GetEntries() << std::endl;
  for( Int_t eventIndex = 0; eventIndex < ch->GetEntries(); eventIndex++){
    ch->GetEntry( eventIndex );
    std::list<Cluster> clist;
    std::list<ClusterTime> ctlist;
    data.getData( clist );
    std::cout<< "Number of Cluster:" <<  clist.size() << std::endl;
    if( clist.size() == 0 ){ continue; }
    clusterTimeAnalyzer->Convert( clist ,ctlist);
    std::cout<< __PRETTY_FUNCTION__ << ":" << __LINE__ <<  std::endl; 
    for( std::list<ClusterTime>::iterator itList = ctlist.begin();
	 itList != ctlist.end();
	 itList++ ){
      std::cout<< (*itList).GetClusterTime() << " : " << (*itList).GetClusterR() << std::endl;
    }    
  }
}
