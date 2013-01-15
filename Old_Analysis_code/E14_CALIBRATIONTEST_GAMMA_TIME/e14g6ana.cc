#include "TFile.h"
#include "TChain.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"

bool user_rec(std::list<Gamma> const &glist, std::vector<Klong> &klVec);
void user_cut(E14GNAnaDataContainer &data, std::vector<Klong> const &klVec);

int main(int argc,char** argv){
  // read argument  
  if(argc!=3){
    std::cout<<"arg err:<< \n usage : bin/e14g6ana input output"<<std::endl
	     <<"input file should be the output of e14gnana/e14clustering."<<std::endl;
    return 0;
  }
  std::cout<<"\n** an example of KL->3pi0 reconstruction**"<<std::endl;
  for(int i=0;i<argc;i++)     std::cout<<" "<<argv[i];
  std::cout<<"\n"<<std::endl;

  std::string ifname = argv[1];
  std::string ofname = argv[2];

  E14GNAnaDataContainer data;

  // set input file
  TChain *inputTree = new TChain("tro");
  inputTree->Add(ifname.c_str());
  std::cout<<"input file: "<<ifname<<std::endl;
  data.setBranchAddress( inputTree );

  // set output file
  TFile *outputFile = new TFile(ofname.c_str(),"RECREATE");
  TTree *outputTree = new TTree("Tree","output from e14g6ana");  
  std::cout<<"output file: "<<ofname<<std::endl;
  data.branchOfKlong( outputTree );

  //
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

    // gamma finding
    std::list<Gamma> glist;
    gFinder.findGamma(clist,glist);

    // kl reconstrunction
    if( glist.size() != 6 ) continue; 
    std::vector<Klong> klVec;
    if(!user_rec(glist,klVec)) continue;
    
    //cuts
    user_cut(data,klVec);
  
    data.setData( klVec );
    outputTree->Fill();
    data.eventID++;
  }
  
  outputTree->Write();
  outputFile->Close();
}



