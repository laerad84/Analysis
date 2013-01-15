////////////////////////////////////////////////////////////////////////
// 
// E14_CLUSTER_BUILDER.cc
// 
////////////////////////////////////////////////////////////////////////
// Make ClusterFile+KlongFile From SumFile.
// Include VETO Information.
//

#include <TFile.h>
#include <TTree.h>
#include "TChain.h"
#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
#include "gnana/E14GNAnaFunction.h"
#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "gamma/GammaFinder.h"
#include "E14_CLUSTER_BUILDER/E14ReadSumFile.h"
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

bool user_rec(std::list<Gamma> const &glist, std::vector<Klong> &klVec);
void user_cut(E14GNAnaDataContainer &data, std::vector<Klong> const &klVec);

int main(int argc, char* argv[])
{  
  std::cout << " Start Cluster " << std::endl; 
  std::string ANALIBstr = std::getenv("ANALYSISLIB");
  std::cout << ANALIBstr << std::endl;
  //read argument 
  int nloop = -1;
  std::string calibrationFile;

  if(argc==4){
    //nloop = atoi(argv[3]);
    //std::cout<<"nRequest== "<<nloop<<"events."<<std::endl;

    calibrationFile = argv[3];
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

  /*
    TChain *trin = new TChain("eventTree00");
    trin->Add(ifname.c_str());
    DigiReader digiReader( trin );
    int outerDetID = digiReader.addDetector("os.");
  */
  
  E14ReadSumFile* ReadSum = new  E14ReadSumFile();
  ReadSum->Add(ifname.c_str());    
  ReadSum->SetOutputFile(ofname.c_str());
  ReadSum->BranchtoTree(ReadSum->trOut);
  // set output file
  //TFile *fout = new TFile(ofname.c_str(),"RECREATE");
  //TTree *trout = new TTree("tro","output from e14clustering");  

  E14GNAnaDataContainer data;
  data.branchOfClusterList( ReadSum->trOut );
  data.branchOfDigi( ReadSum->trOut );  
  data.branchOfKlong( ReadSum->trOut );

  // declare  ClusterFinder and variables
  int nCSIDigi=0;
  int CSIDigiID[3000]      = {0};
  double CSIDigiE[3000]    = {0};
  double CSIDigiTime[3000] = {0};
  double CSICalFactor[3000]= {0};
  ClusterFinder clusterFinder;
  GammaFinder   gFinder;
  
  if( argc  == 3){
    for( int i =0 ;i< 2716; i++){
      CSICalFactor[i] = 1;
    }
  }else{ // argc == 4;
    std::ifstream ifs(calibrationFile.c_str());
    int id;
    double gain;
    double rms;
    while( !ifs.eof() ){
      ifs >> id >> gain >> rms ;
      if( !(gain==0) ){
	CSICalFactor[id] = 14./gain;
      }else{
	CSICalFactor[id] = 0;
      }
    }
  }  
  
  // loop analysis
  //int nentry = trin->GetEntries();
  int nentry = ReadSum->GetEntries();  
  
  std::cout<<"# of entry in input tree =="<<nentry<<std::endl;
  if( nloop<0 || nloop>nentry ) nloop = nentry;
  
  std::cout<<"\n start loop analysis for "<<nloop<<" events..."<<std::endl;  

  for(int ievt=0;ievt<nloop;++ievt){
    if(nloop>100&&ievt%(nloop/100)==0)
      std::cout<<ievt/(nloop/100)<<"%"<<std::endl;
    ReadSum->GetEntry(ievt);
    ReadSum->SumData();
    //trin->GetEntry(ievt);    
    //Clustering
    //digiReader.getCsiDigi(nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);
    
    nCSIDigi = 0; 
    for( int i = 0; i< ReadSum->CsiNumber; ++i){
      Double_t Energy = ReadSum->CsiEne[i]*CSICalFactor[ReadSum->CsiModID[i]];
      if( Energy > 3 ){
	CSIDigiID[nCSIDigi]  = ReadSum->CsiModID[i];
	CSIDigiE[nCSIDigi]   = Energy;
	CSIDigiTime[nCSIDigi]= ReadSum->CsiTime[i];
	nCSIDigi++;
      }
    }
    
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    if( ReadSum->SetTriggerBit()){
      clist = clusterFinder.findCluster(nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);
      gFinder.findGamma(clist,glist);
      if( glist.size() == 6){
	if(user_rec(glist,klVec)){
	  user_cut(data,klVec);
	}
      }	 

      //filling digi-data and Cluster infomation in TTree
      //data.setData( digiReader );
      //trout->Fill();

    }
    data.setData( clist );
    data.setData(klVec);

    ReadSum->Fill();
    data.eventID++;
  }
  // end of analysis
  //trout->Write();
  //fout->Close();

  ReadSum->Write();
  std::cout << "finish!" << std::endl;
  
  return 0;
}
