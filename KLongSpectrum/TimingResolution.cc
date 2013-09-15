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

  TFile* tfOut = new TFile(Form("TimeResolution_%s.root",name),"recreate");
  TH2D* hisResolution = new TH2D("hisResolution","hisResolution",100,0,400,100,-20,20);
  TH2D* hisResolutionAdj = new TH2D("hisResolutionAdj","hisResolutionAdj",100,0,400,100,-20,20);
  TH2D* hisResolutionLY[4];//0: good/good 1:good/bad 2:bad/good 3:bad/bad
  TH2D* hisX = new TH2D("hisX","hisX",2716,0,2716,80,-1000,1000);
  TH2D* hisY = new TH2D("hisY","hisY",2716,0,2716,80,-1000,1000);
  for( int i = 0; i< 4; i++){
    hisResolutionLY[i] = new TH2D(Form("hisResolutionLY_%d",i),Form("hisResolutionLY_%d",i),100,0,400,100,-20,20);
  }
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
  std::ifstream ifs("TimeOffset_10.dat");
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
    if( g0crystalID >= 2240 ){continue; }
    double g0time   = (*git).clusterTimeVec()[0];
    double g0length = TMath::Sqrt( TMath::Power(((*git).x() - klVec[0].vx()),2) 
				 + TMath::Power(((*git).y() - klVec[0].vy()),2) 
				 + TMath::Power(((*git).z() - klVec[0].vz()),2));
    double g0Offset = TimeOffset[g0crystalID];//man->GetT0Offset(g0crystalID);
    double g0Delta  = g0Offset-g0length/sol;//man->GetT0Offset(g0crystalID);
    bool bg0good=false;
    double g0xy[2];

    g0xy[0] = (*git).x();
    g0xy[1] = (*git).y();
    if( g0xy[1] > 0 ){ 
      if( g0xy[0] < -100 ){
	bg0good = true;
      }else{
	bg0good = false;
      }      
    }else{
      if( g0xy[0] < 75 ){
	bg0good = true;
      }else{
	bg0good = false;
      }
    }
    hisX->Fill(g0crystalID, g0xy[0] );
    hisY->Fill(g0crystalID, g0xy[1] );
    git++;
    
    //std::cout << g0crystalID  << "\t" << g0xy[0] << "\t" << g0xy[1] << std::endl;

    
    double g0Ene = (*git).clusterEVec()[0];
    git = glist.begin();
    
    for( int igamma  = 0; igamma < 6; igamma++,git++){
      if((*git).clusterEVec()[0] > 300 && (*git).clusterEVec()[0] < 400 ){
	int id0= (*git).clusterIdVec()[0];
	int id1= (*git).clusterIdVec()[1];
	hisResolution->Fill( (*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-((*git).clusterTimeVec()[0])); 
	hisResolutionAdj->Fill( (*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0])); 
	bool bggood =false;
	double gxy[2];
	gxy[0] = (*git).x();
	gxy[1] = (*git).y();
	if( gxy[1] > 0 ){ 
	  if( gxy[0] < -100 ){
	    bggood = true;
	  }else{
	    bggood = false;
	  }      
	}else{
	  if( gxy[0] < 75 ){
	    bggood = true;
	  }else{
	    bggood = false;
	  }
	}
	if( bg0good ){
	  if( bggood ){
	    hisResolutionLY[0]->Fill((*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0]));
	  }else{
	    hisResolutionLY[1]->Fill((*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0]));
	  }
	}else{
	  if( bggood ){
	    hisResolutionLY[2]->Fill((*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0]));
	  }else{
	    hisResolutionLY[3]->Fill((*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0]));
	  }
	}
      }
    }
  }

  hisResolution->Write();
  hisResolutionAdj->Write();
  hisResolutionLY[0]->Write();
  hisResolutionLY[1]->Write();
  hisResolutionLY[2]->Write();
  hisResolutionLY[3]->Write();
  hisX->Write();
  hisY->Write();
  tfOut->Close();
}
