#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"

#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
#include "TVector2.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/Gamma.h"
#include "cluster/Cluster.h"
#include "cluster/ClusterFinder.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include <vector>
#include <list>
#include <string>
#include <cstdlib>
#include <cstdio>
#include "T0Manager.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "IDHandler.h"

const double KLMass = 497.648;//MeV
double const sol = 299.792458;//[mm/nsec]
double const solc= 80;//[mm/nsec]
double const Pcor[2]={6.49003,0.99254};
double const CsIX0=18.5;//mm

double showerDepth(double x){
  double L = CsIX0*(Pcor[0]+Pcor[1]*log(x/1000.));//mm  
  return L;
}

double showerTimeDelay(Klong kl, Gamma g){
  double depth     = showerDepth(g.e());
  double cosTheta  = abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  double delayTime = 500/solc - depth*(1-cosTheta)/sol;
  return delayTime;
}

double showerTimeDelayAdj(Klong kl, Gamma g){
  double clusterTimingConstant = 0.541812;
  double depth     = showerDepth( g.e() );
  double cosTheta  = TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  //double delayTime = 500/solc - depth*(1-cosTheta)/sol + 0.541812*TMath::Log10(TMath::E())*log(g.e()/1000);
  double delayTime = depth*(1/sol-cosTheta/solc);
  return delayTime;
}

double KLSpectrum(double* x,double*p){
  //double klp   = TMath::Sqrt( x[0]*x[0]-KLMass*KLMass);  
  double klp = x[0];
  double sigma = p[2]*(1-(p[3]+p[4]*klp)*(klp-p[1]));
  double value = p[0]*TMath::Exp(-0.5*(klp-p[1])*(klp-p[1])/(sigma*sigma));
  return value;
}

double gammaLOF( Klong kl, Gamma g){
  double length = 0;
  length = sqrt( pow(g.x()-kl.vx(),2)+ pow(g.y()-kl.vy(),2)+pow(g.z()-kl.vz(),2));
  return length;
}
double HeightDelay( double height ){
  double value = TMath::Log( 1 + 0.03566*TMath::Exp( height/1621));
  return value;
}

//void DistributionTester(){
int main( int argc, char** argv){


  TF1* soltFunc  = new TF1("soltfunc",KLSpectrum,0,12,5);
  TF1* sugarFunc = new TF1("sugarfunc",KLSpectrum,0,12,5);
  const double soltPar[5] = { 1, 1.37,7.48e-1,-2.9e-1,1.68e-2};
  const double sugarPar[5] = {1,1.41991,0.810237,-0.301413,0.017092};
  soltFunc->SetParameters(soltPar);
  sugarFunc->SetParameters(sugarPar);
  IDHandler* handler = new IDHandler();
  const int nFile = 1;
  TFile* tf;
  TTree* tr;
  //char* name = "DATA_NONTIMECAL";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"
  char* name = "DATA_NONTIMECALNOCV";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"

  tf = new TFile(Form("Kl_Total_%s.root",name));
  tr = (TTree*)tf->Get(Form("trKL"));
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  Int_t CsiNumber;
  Int_t CsiModID[3000]={0};
  Double_t CsiSignal[3000]={0};
  Double_t CsiTime[3000]={0};
  Double_t CsiEne[3000]={0};
  E14GNAnaDataContainer data;
  //data.setBranchAddress( tr );
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiModID",CsiModID);//CsiNumber
  tr->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  tr->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  TFile* tfOut = new TFile(Form("kl_Cluster_%s.root",name),"recreate");
  TTree* trOut = new TTree("trCluster","Time Offset Adj");
  data.branchOfClusterList(trOut);
  //data.branchOfDigi(trOut);
  trOut->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsiL1TrigCount[20]/D");
  trOut->Branch("CsiL1nTrig",&CsiL1nTrig,"CsiL1nTrig/I");
  /*
  trOut->Branch("CsiNumber",&CsiNumber,"CsiNumber/I");
  trOut->Branch("CsiSignal",CsiSignal,"CsiSignal[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiModID",CsiModID,"CsiModID[CsiNumber]/I");//CsiNumber
  trOut->Branch("CsiEne",CsiEne,"CsiEne[CsiNumber]/D");//CsiNumber
  trOut->Branch("CsiTime",CsiTime,"CsiTime[CsiNumber]/D");//CsiNumber
  */
  ClusterFinder cFinder;

  Double_t TimeOffset[2716];
  for( int i = 0; i < 2716; i++){
    TimeOffset[i] = 0;
  }
  std::ifstream ifs("TimeOffset_ShowerHeight_15.dat");
  if( !ifs.is_open() ){
    std::cerr << "Timing Calibration File is not opened" << std::endl;
    return -1; 
  }
  int tmpID;
  double tmpOffset;
  while( ifs >> tmpID >> tmpOffset ){
    TimeOffset[tmpID] = tmpOffset;     
  }

  int nTotal = tr->GetEntries(); 
  for( int ievent  =0; ievent < nTotal; ievent++){
      tr->GetEntry( ievent );
      if( CsiNumber > 500 ){ continue; }
      if( (ievent % 1000) == 0){
	std::cout<< ievent <<"/" << nTotal << std::endl;
      }
      for( int i = 0; i< CsiNumber; i++){
	CsiTime[i]=CsiTime[i]-TimeOffset[CsiModID[i]]-HeightDelay(CsiSignal[i]);	
      }
      std::list<Cluster> clist = cFinder.findCluster( CsiNumber, CsiModID, CsiEne, CsiTime );
      data.setData( clist );
      trOut->Fill();
      data.eventID++;
  }
  trOut->Write();
  tfOut->Close();
}
