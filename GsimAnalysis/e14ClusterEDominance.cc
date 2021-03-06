#include "TFile.h"
#include "TChain.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TH1.h"
bool user_rec(std::list<Gamma> const &glist, std::list<Pi0>& piList);
void user_cut(E14GNAnaDataContainer &data,std::list<Pi0> const &piList);
const double pi0mass = 134.9766;
int main(int argc,char** argv){
  // read argument  
  if(argc!=3){
    std::cout<<"arg err:<< \n usage : bin/e14g2ana input output"<<std::endl
	     <<"input file should be the output of e14gnana/e14clustering."<<std::endl;
    return 0;
  }
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
  TTree *outputTree = new TTree("Tree","output from e14g2ana");  
  std::cout<<"output file: "<<ofname<<std::endl;
  data.branchOfPi0List( outputTree );
  //  data.branchOfDigi( outputTree );

  //
  GammaFinder gFinder;


  TH1D* hisSumEnergy = new TH1D("hisSumEnergy","hisSumEnergy",100,0,1);  
  TH1D* hisGammaDepRaw = new TH1D("hisGammaDepRaw","hisGammaDepRaw",100,0,1);
  TH1D* hisGammaDep[10];

  for( int i = 0; i< 10; i++){
    hisGammaDep[i] = new TH1D(Form("hisGammaDep_%d",i),Form("hisGammaDep_%d",i),100,0,1);
  }

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
    if(!user_rec(glist,piList)) continue;
    // cuts
    user_cut(data,piList);


    std::list<Pi0>::iterator pit = piList.begin();
    double EnergySum = 0;
    double FractionSqSum = 0;
    if( TMath::Abs((*pit).m() - pi0mass) < 2 ){
      for( int i = 0; i< (*pit).g1().clusterEVec().size(); i++){
	if( i >= 10 ){ continue; }
	double fractionE = (*pit).g1().clusterEVec()[i]/(*pit).g1().edep();
	hisGammaDep[i]->Fill(fractionE);
	EnergySum+=(fractionE);
	FractionSqSum+=fractionE*fractionE;
      }  
    }
    hisGammaDepRaw->Fill(sqrt(FractionSqSum));

    hisSumEnergy->Fill(EnergySum);
    // fill data to TTree
    data.setData( piList );
    outputTree->Fill();
    data.eventID++;
  }
  hisGammaDepRaw->Write();
  for( int i = 0; i< 10; i++){
    hisGammaDep[i]->Write();
  }
  hisSumEnergy->Write();
  outputTree->Write();
  outputFile->Close();
}

