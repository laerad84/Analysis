#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TStyle.h"
#include "TROOT.h"

#include <gnana/DigiReader.h>
#include <gnana/E14GNAnaDataContainer.h>
#include <cluster/ClusterFinder.h>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>


int 
main( int  argc,char** argv)
{
  //read argument 
  // ./E14_ANA_COSMIC threshold_filename inputFilename outputFilename nevent
  
  int nloop = -1;
  std::string thFilename;
  std::string inFilename;
  std::string outFilename;

  int nloop = -1;
  if( argc == 5 ){
    nloop = atoi(argv[4]);
  }else if ( argc == 3){
    thFilename  = argv[1];
    inFilename  = argv[2];
    outFilename = argv[3];
  }else{
    std::cerr << "The number of arguement is 4 or 5" << std::endl;
    std::cerr << Form("%s [thresholdFilename] [inputFilename] [outputFilename];",argv[0]) << std::endl;
    return -1; 
  }

  TChain *trInput = new TChain("eventTree00");
  trInput->Add(ifname.c_str());

  DigiReader digiReader( trin );

  int outerDetID = digiReader.addDetector("os.");
  int CosmicDet  = digiReader.addDetector("COSMIC.");

  TFile tfOut  = new TFile(outFilename.c_str(), "RECREATE");
  TTree* trOut = new TTree("tro","Output from E14_ANA_COSMIC");
  
  
  E14GNAnaDataContainer data;
  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );


  // declare  ClusterFinder and variables
  int nCSIDigi=0;
  int CSIDigiID[3000]={0};
  double CSIDigiE[3000]={0},CSIDigiTime[3000]={0};
  ClusterFinder clusterFinder;
  
  // loop analysis
  int nentry = trin->GetEntries();
  std::cout<<"# of entry in input tree =="<<nentry<<std::endl;
  if( nloop<0 || nloop>nentry ) nloop = nentry;
  
  std::cout<<"\n start loop analysis for "<<nloop<<" events..."<<std::endl;
  for(int ievt=0;ievt<nloop;ievt++){
    if(nloop>100&&ievt%(nloop/100)==0)
      std::cout<<ievt/(nloop/100)<<"%"<<std::endl;
    
    trin->GetEntry(ievt);
    
    //Clustering
    digiReader.getCsiDigi(nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);
    std::list<Cluster> clist = clusterFinder.findCluster(nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);

    //filling digi-data and Cluster infomation in TTree
    data.setData( digiReader );
    data.setData( clist );
    trout->Fill();
    data.eventID++;
  }

  
  // end of analysis
  trout->Write();
  fout->Close();
  std::cout<<"finish!"<<std::endl;
  
  return 0;
}
