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
  std::string ROOTFILE_SUM   = std::getenv("ROOTFILE_SUMUP");
  std::string ROOTFILE_WAV   = std::getenv("ROOTFILE_WAV");
  std::string ROOTFILE_SIM   = std::getenv("ROOTFILE_PI0SIMCONV");
  std::string ROOTFILE_CONV;
  switch ( FileType ){
  case 0:
    ROOTFILE_CONV = ROOTFILE_SUM.c_str();
    InputFilename = "%s/Conv_e14_AL_Target.mac_1000000_%d.root";//ConvFileDir, RunNumber
    OutputFilename= "Pi0SIM.root";    
    break;
  case 1:
    ROOTFILE_CONV = ROOTFILE_SUM.c_str();
    InputFilename = "%s/Sum%d.root";
    OutputFilename= "Pi0SUM.root";
    break;
  case 2:
    ROOTFILE_CONV = ROOTFILE_WAV.c_str();
    InputFilename = "%s/run_wav_%d_Cal_FNL_COS.root";
    OutputFilename= "Pi0WAV.root";
    break;
  default:
    return -1;
  }

  const Int_t nCsI = 2716;
  Int_t CsiNumber;
  Short_t CsiID[nCsI];
  Double_t CsiEne[nCsI];
  Double_t CsiTime[nCsI];

  TChain* trIn = new TChain("T");
  switch(FileType){
  case 0:
    for( int i = 0; i< 200; i++){
      trIn->Add(Form(InputFilename.c_str(),ROOTFILE_CONV.c_str(),i));
    }
    break;
  case 1:
    for( int i = 4502; i < 4526; i++){
      trIn->Add(Form(InputFilename.c_str(),ROOTFILE_CONV.c_str(),i));
    }
    break;
  case 2:
    for( int i = 4502; i< 4526; i++){
      trIn->Add(Form(InputFilename.c_str(),ROOTFILE_CONV.c_str(),i));
    }
    break;
  }
  // Set Branch Address // 
  switch ( FileType ){
  case 0:
    trIn->SetBranchAddress("CsiNumber",&CsiNumber);
    trIn->SetBranchAddress("CsiModID" ,CsiID);//CsiNumber
    trIn->SetBranchAddress("CsiEne"   ,CsiEne);//CsiNumber
    trIn->SetBranchAddress("CsiTime"  ,CsiTime);//CsiNumber
    break;
  case 1:
    trIn->SetBranchAddress("CsiNumber",&CsiNumber);
    trIn->SetBranchAddress("CsiModID" ,CsiID);//CsiNumber
    trIn->SetBranchAddress("CsiEne"   ,CsiEne);//CsiNumber
    trIn->SetBranchAddress("CsiTime"  ,CsiTime);//CsiNumber
    break;
  case 2:
    trIn->SetBranchAddress("CsiNumber",&CsiNumber);
    trIn->SetBranchAddress("CSIDigiID",CsiID);//CsiNumber
    trIn->SetBranchAddress("CsiEne"   ,CsiEne);//CsiNumber
    trIn->SetBranchAddress("CsiTime"  ,CsiTime);//CsiNumber
    break;
  }

  // set output file
  TFile *outputFile = new TFile(Form(OutputFilename.c_str(),""),"RECREATE");
  TTree *outputTree = new TTree("Tree","output from e14g2ana");  

  // set Data file 
  E14GNAnaDataContainer data;
  data.branchOfClusterList( outputTree );
  data.branchOfDigi( outputTree );
  data.branchOfPi0List( outputTree );

  ClusterFinder clusterFinder;
  GammaFinder gFinder;

  // set Variables
  int    nCSIDigi          = 0;
  int    CSIDigiID[3000]   = {0};
  double CSIDigiE[3000]    = {0};
  double CSIDigiTime[3000] = {0};

  /////////////////////////////////////////////////////////////////////////
  // loop
  /////////////////////////////////////////////////////////////////////////

  std::cout<< "Pi0 Calibration" << std::endl;
  // loop analysis
  int nloop = trIn->GetEntries();
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
    trIn->GetEntry( ievt );
    nCSIDigi = 0; 
    for( int ich = 0; ich < CsiNumber;ich++){
	if( CsiEne[ich] > 3 ){
	  CSIDigiID[nCSIDigi]   = CsiID[ich];
	  CSIDigiE[nCSIDigi]    = CsiEne[ich];
	  CSIDigiTime[nCSIDigi] = CsiTime[ich];
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

