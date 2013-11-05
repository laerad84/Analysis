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

  if(argc==5){
    nloop = std::atoi(argv[3]);
    std::cout<<"nRequest== "<<nloop<<"events."<<std::endl;
  }else if(argc!=4 && argc!= 3){
    std::cerr << "Argument error."<<std::endl
	      <<"usege:  bin/e14clustering input output CalibrationPrecision [nRequest]" <<std::endl;
    return 1;
  }

  std::string ifname = argv[1];
  std::string ofname = argv[2];
  std::cout<<"input file: "<<ifname<<std::endl;
  std::cout<<"output file: "<<ofname<<std::endl;
  Int_t CalibrationPrecision = 0;
  if( argc >=4  ){
    CalibrationPrecision = std::atoi(argv[3]);
  }
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
  double CSICal[3000]={1};
  for( int i = 0; i< 3000; i++){
    CSICal[i] = 1;
  }
  std::string HOMEDIR = std::getenv("HOME");
  if( CalibrationPrecision != 0 ){
    std::ifstream ifs(Form("%s/local/Analysis/GsimAnalysis/CalFactor/CalFactor_%d.txt",HOMEDIR.c_str(),CalibrationPrecision));
    if( !ifs.is_open() ){ 
      std::cout<< "Failed to open" << std::endl;
      return -1; 
    }
    int tmpID; 
    double tmpCal;
    while( ifs >> tmpID >> tmpCal ){
      CSICal[tmpID] = tmpCal;
    }
  }

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
    double AdjCSIDigiE[3000]={0}; 
    for( int i = 0; i < nCSIDigi;i++){
      AdjCSIDigiE[i]=CSIDigiE[i]*CSICal[CSIDigiID[i]];
    }
    std::list<Cluster> clist = clusterFinder.findCluster(nCSIDigi,CSIDigiID,AdjCSIDigiE,CSIDigiTime);

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
