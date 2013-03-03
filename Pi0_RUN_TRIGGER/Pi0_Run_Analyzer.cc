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
#include "Pi0_RUN_TRIGGER/user_func.h"

//#include "E14ReadSumFile.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <list>

// Pi0_Run_Analyzer [Type:0(SIM),1(SUM),2(WAV)] RunNumber
int 
main(int argc, char** argv){
  Int_t FileType = atoi(argv[0]);
  std::string InputFilename;
  std::string OutputFilename;
  switch ( FileType ){
  case 0:
    InputFilename = "%s/%d.root";//ConvFileDir, RunNumber
    OutputFilename= "%s/%d.root";
    break;
  case 1:
    InputFilename = "%s/%d.root";
    OutputFilename= "%s/%d.root";
    break;
  case 2:
    InputFilename = "%s/%d.root";
    OutputFilename= "%s/%d.root";
    break;
  default:
    return -1;
  }

  TFile* tfIn = new TFile(Form(InputFilename.c_str(),ROOTFILE_PI0CONV.c_str(),RunNumber));
  TTree* trIn = (TTree*)trIn->Get("T");
  const Int_t nCsI = 2716;
  Int_t CsiNumber;
  Short_t CsiID[nCsI];
  Double_t CsiEne[nCsI];
  Double_t CsiTime[nCsI];

  switch ( FileType ){
  case 0:
    trIn->SetBranchAddress("CsiNumber",&CsiNumber);
    trIn->SetBranchAddress("CsiModID",CsiID);//CsiNumber
    trIn->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
    trin->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
    break;
  case 1:
    trIn->SetBranchAddress("CsiNumber",&CsiNumber);
    trIn->SetBranchAddress("CsiID",CsiID);//CsiNumber
    trIn->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
    trin->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
    break;
  case 2:
    trIn->SetBranchAddress("CsiNumber",&CsiNumber);
    trIn->SetBranchAddress("CSIDigiID",CsiID);//CsiNumber
    trIn->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
    trin->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
    break;
  }
  

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
    for( int ich = 0; ich < CsiNumber;ich++){
	double Energy = CsiEne[ich];
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

