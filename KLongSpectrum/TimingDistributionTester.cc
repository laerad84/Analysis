#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"

#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
#include "TVector2.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/Gamma.h"
#include "cluster/Cluster.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include <vector>
#include <list>
#include "TH2.h"
#include <string>
#include <cstdlib>
#include <cstdio>
#include "T0Manager.h"
#include "TDirectory.h"
#include "TProfile.h"
const double KLMass = 497.648;//MeV

double KLSpectrum(double* x,double*p){
  //double klp   = TMath::Sqrt( x[0]*x[0]-KLMass*KLMass);  
  double klp = x[0];
  double sigma = p[2]*(1-(p[3]+p[4]*klp)*(klp-p[1]));
  double value = p[0]*TMath::Exp(-0.5*(klp-p[1])*(klp-p[1])/(sigma*sigma));
  return value;
}


//void DistributionTester(){
int main( int argc, char** argv){

  int nTimeIteration = atoi(argv[1]);

  TF1* soltFunc  = new TF1("soltfunc",KLSpectrum,0,12,5);
  TF1* sugarFunc = new TF1("sugarfunc",KLSpectrum,0,12,5);
  const double soltPar[5] = { 1, 1.37,7.48e-1,-2.9e-1,1.68e-2};
  const double sugarPar[5] = {1,1.41991,0.810237,-0.301413,0.017092};
  soltFunc->SetParameters(soltPar);
  sugarFunc->SetParameters(sugarPar);

  const int nFile = 1;
  TFile* tf;
  TTree* tr;
  char* name = "WAV";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"

  tf = new TFile(Form("Kl_Total_%s.root",name));
  tr = (TTree*)tf->Get(Form("trKL"));
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];

  TFile* tfOut = new TFile(Form("DistributionTimeData_%s_%d.root",name,nTimeIteration),"recreate");
  TH2D* hisTimeID;
  TH2D* hisAdjTimeID;
  hisTimeID    = new TH2D(Form("hisTimeID_%d",nTimeIteration),Form("hisTimeID_%d",nTimeIteration),
			  2716,0, 2716,500,-40,60);
  hisAdjTimeID = new TH2D(Form("hisAdjTimeID_%d",nTimeIteration),Form("hisAdjTimeID_%d",nTimeIteration),
			  2716,0, 2716,500,-40,60);
  double sol = 299.792458;//[mm/nsec]
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  
  /*T0Manager* man = new T0Manager();
  if( nTimeIteration > 0 ){
    if( !(man->ReadFile(Form("TimeOffset_%d.dat",nTimeIteration)))){
      std::cout<< "No file exist" << std::endl;
      return -1;
    }
  }
  */
  //man->PrintOffset();

  Double_t TimeOffset[2716]={0};
  if( nTimeIteration > 0 ){
    std::ifstream ifs(Form("TimeOffset_%d.dat",nTimeIteration));
    if( !ifs.is_open() ){
      std::cout<< "No CalibrationFile" << std::endl;
      return -1;
    }
    int tmpID;
    double tmpOffset;
    while( ifs >> tmpID >> tmpOffset ){
      TimeOffset[tmpID] = tmpOffset;
      std::cout<<tmpID << "\t" <<  TimeOffset[tmpID] << std::endl;
    }
  }



  for( int ievent = 0; ievent < tr->GetEntries(); ievent++){      
    tr->GetEntry(ievent);
    //if( ievent  >= 100000 ){ break ; } 
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    data.getData(clist);
    data.getData(glist);
    data.getData(klVec);
    
    
    double klpos[3]={0};
    klpos[0] = klVec[0].vx();
    klpos[1] = klVec[1].vy();
    klpos[2] = klVec[2].vz();
    if( klVec[0].chisqZ() > 10 ){ continue; }
    std::list<Gamma>::iterator git = glist.begin();
    int    g0crystalID = (*git).clusterIdVec()[0];
    double g0time   = (*git).clusterTimeVec()[0];
    double g0length = TMath::Sqrt( TMath::Power(((*git).x() - klVec[0].vx()),2) 
				 + TMath::Power(((*git).y() - klVec[0].vy()),2) 
				 + TMath::Power(((*git).z() - klVec[0].vz()),2));
    double g0Offset = TimeOffset[g0crystalID];//man->GetT0Offset(g0crystalID);
    double g0Delta  = g0Offset-g0length/sol;//man->GetT0Offset(g0crystalID);
    git++;
    if( g0crystalID >= 2240 ){continue; }
    
    double g0Ene = (*git).clusterEVec()[0];
    if(g0Ene > 400 ){ continue; }
    if(g0Ene < 100 ){ continue; }

    for( int igamma = 1; igamma < 6; igamma++,git++){
      int crystalID = (*git).clusterIdVec()[0];
      double Ene = (*git).clusterEVec()[0];
      if(Ene > g0Ene){ continue; }
      if(Ene < 100 ){ continue; }
      double Offset = TimeOffset[crystalID];//man->GetT0Offset(crystalID);
      double length = TMath::Sqrt( TMath::Power(((*git).x() - klVec[0].vx()),2) 
				   + TMath::Power(((*git).y() - klVec[0].vy()),2) 
				   + TMath::Power(((*git).z() - klVec[0].vz()),2));
      double Delta = Offset-length/sol;//man->GetT0Offset(crystalID);
      hisTimeID->Fill(crystalID,((*git).clusterTimeVec()[0]-Offset)-(g0time-g0Offset));
      hisAdjTimeID->Fill(crystalID,((*git).clusterTimeVec()[0]-Delta)-(g0time-g0Delta));
    }
  }
    /// Evaluation T0 /// 

  TF1* func = new TF1("func","gaus",-20,20);
  hisTimeID->FitSlicesY(func);
  hisAdjTimeID->FitSlicesY(func);
  TProfile* profTimeID = hisTimeID->ProfileX();
  TProfile* profAdjTimeID = hisAdjTimeID->ProfileX();

  TH1D *hisTimeID_0 = (TH1D*)gDirectory->Get("hisTimeID_0");
  TH1D *hisTimeID_1 = (TH1D*)gDirectory->Get("hisTimeID_1");
  TH1D *hisTimeID_2 = (TH1D*)gDirectory->Get("hisTimeID_2");
  TH1D *hisAdjTimeID_0 = (TH1D*)gDirectory->Get("hisAdjTimeID_0");
  TH1D *hisAdjTimeID_1 = (TH1D*)gDirectory->Get("hisAdjTimeID_1");
  TH1D *hisAdjTimeID_2 = (TH1D*)gDirectory->Get("hisAdjTimeID_2");

  
  std::cout << "Slice" << std::endl;
  int calibrated[2716] = {-1};
  int nCalibrated = 0;
  double mean=0;
  std::ofstream ofs(Form("TimeOffset_%d.dat",nTimeIteration+1));
  for( int i = 1; i<= profAdjTimeID->GetNbinsX(); i++){
    if( profAdjTimeID->GetBinEntries(i) < 1 ){ continue; }
    calibrated[i-1] = 1;
    mean += profAdjTimeID->GetBinContent(i);
    nCalibrated++;
    std::cout<< i-1 << "\t" << profAdjTimeID->GetBinContent(i) << std::endl;
  }
  mean /= nCalibrated;
  double timeCalFactor[2716] = {0};
  TH1D* TimeOffsetDistribution = new TH1D(Form("TimeOffsetDistribution_%d",nTimeIteration),
					  Form("TimeOffsetDistribution_%d",nTimeIteration),
					  200,-20,20);
  for( int i = 0; i< 2716; i++){
    if( calibrated[i] >0){
      timeCalFactor[i] = profAdjTimeID->GetBinContent(i+1) - mean;
      TimeOffsetDistribution->Fill( profAdjTimeID->GetBinContent(i+1) - mean);
    }else{
      timeCalFactor[i] = 0;
    }
    ofs << i  << "\t" << timeCalFactor[i]+TimeOffset[i] << "\n";
  }
  ofs.close();
  

  hisTimeID->Write();
  hisAdjTimeID->Write();
  //TimeOffsetDistribution->Write();
  /*
  hisTimeID_0->Write();
  hisTimeID_1->Write();
  hisTimeID_2->Write();
  hisAdjTimeID_0->Write();
  hisAdjTimeID_1->Write();
  hisAdjTimeID_2->Write();
  */
  tfOut->Close();
}
