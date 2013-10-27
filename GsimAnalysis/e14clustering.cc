#include <TFile.h>
#include <TTree.h>
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


Int_t main(Int_t argc,char** argv)
{
  //read argument 
  int nloop = -1;
  if(argc==4){
    nloop = atoi(argv[3]);
    std::cout<<"nRequest== "<<nloop<<"events."<<std::endl;
  }else if(argc!=3){
    std::cerr << "Argument error."<<std::endl
	      <<"usege:  bin/e14clustering input output [nRequest]" <<std::endl;
    return 1;
  }

  std::string ifname = argv[1];
  std::string ofname = argv[2];
  std::cout<<"input file: "<<ifname<<std::endl;
  std::cout<<"output file: "<<ofname<<std::endl;

  // set input file
  TChain *trin = new TChain("eventTree00");
  trin->Add(ifname.c_str());
  DigiReader digiReader( trin );
  int outerDetID = digiReader.addDetector("os.");
  
  // set output file
  TFile *fout = new TFile(ofname.c_str(),"RECREATE");
  TTree *trout = new TTree("tro","output from e14clustering");  
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
