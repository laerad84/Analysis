#include "TFile.h"
#include "TChain.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TApplication.h"
#include "TH1D.h"
#include "TCanvas.h"

bool user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList,double& mass);
void user_cut(E14GNAnaDataContainer &data,std::list<Pi0> const &piList);

int main(int argc,char** argv){
  // read argument  
  if(argc!=3){
    std::cout<<"arg err:<< \n usage : bin/e14g2ana input output"<<std::endl
	     <<"input file should be the output of e14gnana/e14clustering."<<std::endl;
    return 0;
  }
  std::string ifname = argv[1];
  std::string ofname = argv[2];
  
  TApplication* app = new TApplication("app",&argc, argv);
  E14GNAnaDataContainer data;
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
  GammaFinder gFinder;
  TH1D* hist = new TH1D("his","",1000,0,2000);

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

    // find gammacluster 
    std::list<Gamma> glist;
    gFinder.findGamma(clist,glist);

    // pi0 reconstrunction and position correction for angle dependency
    if( glist.size() != 2 ) continue; 
    std::list<Pi0> piList; 
    double mass;
    if(!user_rec(glist,piList,mass)) continue;
    std::cout<< mass << std::endl;
    hist->Fill(mass);
    
    // cuts
    user_cut(data,piList);    
    
    // fill data to TTree
    data.setData( piList );
    outputTree->Fill();
    data.eventID++;
  }

  TCanvas* can = new TCanvas("can","",800,800);
  hist->Draw();
  can->Update();
  hist->Write();
  outputTree->Write();
  //outputFile->Close();  
  app->Run();
}

