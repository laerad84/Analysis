#include "TFile.h"
#include "TChain.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "cluster/ClusterFinder.h"
#include "cluster/Cluster.h"
#include <fstream>
#include <cstdlib>
#include <cstdio>


bool user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList);
void user_cut(E14GNAnaDataContainer &data,std::list<Pi0> const &piList);

int main(int argc,char** argv){
  // read argument  
  if(argc!=3 && argc !=4){
    std::cout<<"arg err:<< \n usage : bin/e14g2ana input output"<<std::endl
	     <<"input file should be the output of e14gnana/e14clustering."<<std::endl;
    return 0;
  }
  std::string ifname = argv[1];
  std::string ofname = argv[2];
  Int_t       CalPrecision=0; 
  if( argc == 4 ){
    CalPrecision= atoi( argv[3]);
  }
  E14GNAnaDataContainer data;

  Double_t Output[2716];
  for( int i = 0; i< 2716; i++){
    Output[i] = 1;
  }
  int tmpID;
  double tmpOutput;
  std::string ANALYSISLIB=std::getenv("ANALYSISLIB");
  if( CalPrecision != 0){
    std::ifstream ifs(Form("%s/Data/SIM_Gain/Gain_%d.txt",ANALYSISLIB.c_str(),CalPrecision));
    if(!ifs.is_open()){ std::cout<< "File Not Opened" << std::endl;return -1;}
    while( ifs >> tmpID >> tmpOutput ){
      Output[tmpID] = tmpOutput;
      std::cout<< Output[tmpID] << std::endl;
    }

  }

  // set input file
  TChain *inputTree = new TChain("tro");
  inputTree->Add(ifname.c_str());
  std::cout<<"input file: "<<ifname<<std::endl;
  data.setBranchAddress( inputTree );

  // set output file
  TFile *outputFile = new TFile(ofname.c_str(),"RECREATE");
  TTree *outputTree = new TTree("Tree","output from e14g2ana");  
  std::cout<<"output file: "<<ofname<<std::endl;
  data.branchOfPi0List( outputTree );
  //  data.branchOfDigi( outputTree );
  //
  ClusterFinder cFinder;
  GammaFinder gFinder;

  // loop analysis
  int nloop = inputTree->GetEntries();
  std::cout<<"start loop analysis"<<std::endl;
  std::cout<<"# of entry : "<<nloop<<std::endl;
  for( int ievt=0; ievt<nloop; ievt++ ){
    if(ievt%(nloop/10)==0&&nloop>100)
      std::cout<<ievt/(nloop/100)<<"%"<<std::endl;

    // read data
    inputTree->GetEntry( ievt );
    std::list<Cluster> clist;
    data.getData(clist);
    std::list<Cluster>::iterator cit;

    // Apply Output to cluster //
    int nDigi=0;
    int CsiID[2716];
    double CsiE[2716];
    double CsiT[2716];
    nDigi = 0;
    for( int i = 0; i< 2716; i++){
      CsiID[i] = 0;
      CsiE[i] = 0;
      CsiT[i] = 0;
    }

    for( cit = clist.begin(); cit != clist.end(); cit++){
      for( int i = 0; i< (*cit).clusterIdVec().size(); i++){
	CsiID[i] = (*cit).clusterIdVec().at(i);
	CsiE[i]  = (*cit).clusterEVec().at(i) * Output[i];
	CsiT[i]  = (*cit).clusterTimeVec().at(i);
	nDigi++;
      }
    }
    std::list<Cluster> clistAdj=cFinder.findCluster(nDigi,CsiID,CsiE,CsiT);

    // find gammacluster 
    std::list<Gamma> glist;
    gFinder.findGamma(clistAdj,glist);

    // pi0 reconstrunction and position correction for angle dependency
    if( glist.size() != 2 ) continue; 
    std::list<Pi0> piList; 
    if(!user_rec(glist,piList)) continue;


    // cuts
    user_cut(data,piList);
    
  
    // fill data to TTree
    data.setData( piList );
    outputTree->Fill();
    data.eventID++;
  }
  
  outputTree->Write();
  outputFile->Close();
}

