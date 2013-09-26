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
#include "IDHandler.h"

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

  IDHandler* handler = new IDHandler();
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
  TH2D* hisResolutionLY[2];//0: good/good 1:good/bad 2:bad/good 3:bad/bad
  TH2D* hisResolutionLY_Neighbor[2]; 
  TH2D* hisX = new TH2D("hisX","hisX",2716,0,2716,80,-1000,1000);
  TH2D* hisY = new TH2D("hisY","hisY",2716,0,2716,80,-1000,1000);
  for( int i = 0; i< 2; i++){
    hisResolutionLY[i] = new TH2D(Form("hisResolutionLY_%d",i),Form("hisResolutionLY_%d",i),100,0,400,100,-20,20);
    hisResolutionLY_Neighbor[i] = new TH2D(Form("hisResolutionLY_Neighbor_%d",i),Form("hisResolutionLY_Neighbor_%d",i),100,0,400,100,-20,20);
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
      bool bggood =false;
      int id0= (*git).clusterIdVec()[0];
      int id1= (*git).clusterIdVec()[1];
      if((*git).clusterEVec()[0] > 300 && (*git).clusterEVec()[0] < 400 ){
	hisResolution->Fill( (*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-((*git).clusterTimeVec()[0])); 
	hisResolutionAdj->Fill( (*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0]));

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
	if( bggood ){
	  hisResolutionLY[0]->Fill((*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0]));
	}else{
	  hisResolutionLY[1]->Fill((*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0]));
	}
      }

      double x0 = (*git).coex();
      double y0 = (*git).coey();
      for( int iid = 0; iid < (*git).clusterIdVec().size()-1; iid++){
	//for( int iid = 0; iid < 1; iid++){
	int TestID[2]={-1};
	double x[2]={0};
	double y[2]={0}; 
	double R[2]={0};
	double D[2]={0};
	double t[2]={0};
	double L[2]={0};
	double E[2]={0};
	double T[2]={0};
	TestID[0]=(*git).clusterIdVec()[iid];

	handler->GetMetricPosition((*git).clusterIdVec()[iid],x[0],y[0]);
	x[0] = -x[0];
	L[0] = sqrt( (x0-x[0])*(x0-x[0])+(y0-y[0])*(y0-y[0]));
	R[0] = ((x[0]-x0)*x0+(y[0]-y0)*y0)/sqrt( x0*x0+y0*y0 );
	D[0] = sqrt(L[0]*L[0]-R[0]*R[0]);
	E[0] = (*git).clusterEVec()[iid];	  
	T[0] = (*git).clusterTimeVec()[iid];
	for( int jid = iid+1; jid < (*git).clusterIdVec().size(); jid++){

	  TestID[1]=(*git).clusterIdVec()[jid];
	  handler->GetMetricPosition((*git).clusterIdVec()[iid],x[1],y[1]);
	  x[1] = -x[1];
	  L[1] = sqrt( (x0-x[1])*(x0-x[1])+(y0-y[1])*(y0-y[1]));
	  R[1] = ((x[1]-x0)*x0+(y[1]-y0)*y0)/sqrt( x0*x0+y0*y0 );
	  D[1] = sqrt(L[1]*L[1]-R[1]*R[1]);
	  E[1] = (*git).clusterEVec()[jid];	  
	  T[1] = (*git).clusterTimeVec()[jid];	  
	  if( D[0] > 20 || D[1] > 20 ){ continue; }
	  if( abs(R[0]) > 7 || abs(R[1]) > 7 ){ continue; }
	  if( abs(R[0]-R[1]) > 12.5 ){ continue; }
	  if( E[1] < E[0]*0.8 ){ continue; }
	  if( bggood ){
	    hisResolutionLY_Neighbor[0]->Fill(E[1],T[1]-TimeOffset[TestID[1]]-(T[0]-TimeOffset[TestID[0]]));
	  }else{
	    hisResolutionLY_Neighbor[1]->Fill(E[1],T[1]-TimeOffset[TestID[1]]-(T[0]-TimeOffset[TestID[0]]));
	  }
	}
      }
    }
  }


  hisResolution->Write();
  hisResolutionAdj->Write();
  for( int i = 0; i< 2; i++){
    hisResolutionLY[i]->Write();
    hisResolutionLY_Neighbor[i]->Write();
  }

  hisX->Write();
  hisY->Write();
  tfOut->Close();
}
