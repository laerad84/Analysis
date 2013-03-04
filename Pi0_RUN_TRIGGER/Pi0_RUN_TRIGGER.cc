#include "TApplication.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/GammaFinder.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "user_func.h"

#include "E14ReadSumFile.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <list>


int 
main(int argc,char** argv){
  // read argument  
  // 
  if( argc != 3 && argc != 4){
    std::cout<< "<<<<>>>Arguement Error<<<>>>"
	     <<  "Usage:" << argv[0] << " [inputFile] [outputFile] :[CalibrationFile]:"
	     << std::endl;
    return -1;
  }

  std::string ifname = argv[1];
  std::string ofname = argv[2];
  std::string cfname;
  if( argc ==4){
    cfname = argv[3];
  }

  //TApplication* app = new TApplication("app",&argc, argv);

  // set input file 
  E14ReadSumFile* ReadSum = new E14ReadSumFile();
  ReadSum->Add(ifname.c_str());
  
  // set output file
  TFile *outputFile = new TFile(ofname.c_str(),"RECREATE");
  TTree *outputTree = new TTree("Tree","output from e14g2ana");  

  // set Data file 
  E14GNAnaDataContainer data;
  data.branchOfClusterList( outputTree );
  data.branchOfDigi( outputTree );
  data.branchOfPi0List( outputTree );

  ClusterFinder clusterFinder;
  GammaFinder gFinder;

  // set Variables
  int nCSIDigi = 0;
  int CSIDigiID[3000] = {0};
  double CSIDigiE[3000]={0}, CSIDigiTime[3000]={0};
  double CSICalFactor[3000]={0};
  for( int ich = 0; ich < 3000; ++ich){
    CSICalFactor[ich] = 1;
  }

  // if Calibration Factor File is defined //

  if( argc ==4){
    std::cout << "OpenFile:" << cfname << std::endl;
    std::ifstream ifs(cfname.c_str());
    if( !ifs.is_open() ){ return -1;}
    int id;
    double CalibrationFactor;
    while( ifs >> id >> CalibrationFactor){
      std::cout<< CalibrationFactor << std::endl;
      CSICalFactor[id] = CalibrationFactor;
    }
  }

  /////////////////////////////////////////////////////////////////////////
  // loop
  /////////////////////////////////////////////////////////////////////////

  std::cout<< "Pi0 Calibration" << std::endl;
  // loop analysis
  int nloop = ReadSum->GetEntries();
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

    // read data
    ReadSum->GetEntry( ievt );
    nCSIDigi = 0; 
    for( int ich = 0; ich < ReadSum->CsiNumber; ++ich){
      double Energy = ReadSum->CsiEne[ich];//*CSICalFactor[ReadSum->CsiModID[ich]]/0.99;
      //std::cout<< ReadSum->CsiModID[ich]  << " : " <<  Energy << std::endl;
      //double Energy = ReadSum->CsiEne[ich];
      if( Energy > 3 ){
	CSIDigiID[nCSIDigi]   = ReadSum->CsiModID[ich];
	CSIDigiE[nCSIDigi]    = Energy;
	CSIDigiTime[nCSIDigi] = ReadSum->CsiTime[ich];
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
    
    // Al target position 2622 from CsI

    //double position =3484.;
    //double position = 3534.;//shiomisan    
    double position = 3526;//20120906

    if(!user_rec(glist,piList,mass,position)) continue;
    std::cout<< mass << std::endl;
    hist->Fill(mass);
    // cuts
    //user_cut(data,piList);    
    
    // fill data to TTree
    data.setData( piList );
    outputTree->Fill();
    data.eventID++;
  }

  //TCanvas* can = new TCanvas("can","",800,800);
  //hist->Draw();
  hist->Write();
  outputTree->Write();
  outputFile->Close();  
  //app->Run();

}

