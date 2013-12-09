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
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include <vector>
#include <list>
#include <string>
#include <cstdlib>
#include <cstdio>
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
  char* name = "KLCollection_Minimum";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"

  tf = new TFile(Form("%s.root",name));
  tr = (TTree*)tf->Get(Form("RecTree"));
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];

  TFile* tfOut = new TFile(Form("TimeResolution_ShowerTimeShape_%s.root",name),"recreate");

  TH2D* hisTimeD = new TH2D("hisTimeD","hisTimeD",12,0,150,100,-15,15);
  TH2D* hisED    = new TH2D("hisED"   ,"hisED"   ,12,0,150,100,0,400);
  TH2D* hisTimeR = new TH2D("hisTimeR","hisTimeR",13,-150-12.5,150+12.5,150,-15,15);
  TH2D* hisER    = new TH2D("hisER"   ,"hisER"   ,13,-150-12.5,150+12.5,100,0,400);
  TH2D* hisTimeL = new TH2D("hisTimeL","hisTimeL",12,0,150,150,-15,15);
  TH2D* hisEL    = new TH2D("hisEL"   ,"hisEL"   ,12,0,150,100,0,400);
  TH2D* hisTimeZitter= new TH2D("hisTimeZitter","hisTimeZitter",2716,0,2716,300,-30,30);
  TH3D* hisRDE   = new TH3D("hisRDE","hisRDE",13,-150-12.5,150+12.5,12,0,150,100,  0,400);
  TH3D* hisRDT   = new TH3D("hisRDT","hisRDT",13,-150-12.5,150+12.5,12,0,150,150,-15,15);

  const int nAngleCut = 7;
  Double_t AngleCut[nAngleCut]={0.02,0.03,0.04,0.05,0.06,0.07,0.4};
  TH2D* hisTimeRTheta[nAngleCut];
  TH2D* hisERTheta[nAngleCut];
  TH2D* hisTimeDTheta[nAngleCut];
  TH2D* hisEDTheta[nAngleCut];
  TH3D* hisRDT_Theta[nAngleCut];
  TH3D* hisRDE_Theta[nAngleCut];
  for( int i = 0; i< nAngleCut; i++){
    hisTimeRTheta[i] = new TH2D(Form("hisTimeR_%d",i),Form("hisTimeR_%d",i),13,-150-12.5,150+12.5,150,-15,15);
    hisERTheta[i]    = new TH2D(Form("hisER_%d",i)   ,Form("hisER_%d",i)   ,13,-150-12.5,150+12.5,100,0,400);
    hisTimeDTheta[i] = new TH2D(Form("hisTimeD_%d",i) ,Form("hisTimeD_%d",i),12,0,150,150,-15,15);
    hisEDTheta[i]    = new TH2D(Form("hisED_%d",i)   ,Form("hisED_%d",i)   ,12,0,150,100,0,400);
    hisRDT_Theta[i]  = new TH3D(Form("hisRDT_Theta_%d",i),Form("hisRDT_Theta_%d",i),13,-150-12.5,150+12.5,12,0,150,150,-15,15);
    hisRDE_Theta[i]  = new TH3D(Form("hisRDE_Theta_%d",i),Form("hisRDE_Theta_%d",i),13,-150-12.5,150+12.5,12,0,150,100,  0,400);
  }

  TH1D* hisInjectionAngle = new TH1D("hisInjectionAngle","hisInjectionAngle",100,0,1);


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

  /*
  Double_t TimeOffset[2716]={0};
  std::ifstream ifs("TimeOffset_Shower_10.dat");
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
  */

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
    if( klVec[0].chisqZ() > 50 ){ continue; }
    if( klpos[2] > 5500 ){continue; }
    Double_t x0(0),y0(0);
    Double_t x1(0),y1(0);
    Double_t t0(0),t1(0); 
    Double_t L(0);
    Double_t E0(0),E1(0);
    Double_t R(0),D(0);
    Double_t InjectionAngle={0};
    std::list<Gamma>::iterator git = glist.begin();
    for( int igamma  = 0; igamma < 6; igamma++,git++){

      //if((*git).clusterEVec()[0] > 300 && (*git).clusterEVec()[0] < 400 ){

	if( abs((*git).x()) > 600 || abs((*git).y() > 600 )){ continue; } 
	Double_t TimeMean=0;
	Int_t nCenter=0;

	for( int iid = 0; iid <(*git).clusterIdVec().size(); iid++){
	  if( (*git).clusterEVec()[iid] > 20 ){
	    TimeMean+= (*git).clusterTimeVec()[iid];//-TimeOffset[(*git).clusterIdVec()[iid]];
	    nCenter++;
	  }
	  if( iid > 2 ){ continue; }
	}
	t0 = TimeMean/nCenter;	
	//t0 = (*git).clusterTimeVec()[0];
	x0 = (*git).coex();
	y0 = (*git).coey();
	double length = TMath::Sqrt( TMath::Power(((*git).x() - klVec[0].vx()),2) 
				     + TMath::Power(((*git).y() - klVec[0].vy()),2) 
				     + TMath::Power(((*git).z() - klVec[0].vz()),2));	

	InjectionAngle = acos(abs((*git).z() - klVec[0].vz())/length)/TMath::Pi();
	hisInjectionAngle->Fill(InjectionAngle);
	//handler->GetMetricPosition((*git).clusterIdVec()[0],x0,y0);
	//t0 = (*git).clusterTimeVec()[0];
	E0 = (*git).clusterEVec()[0];
	int ThetaIndex=0;
	for( int index = 0; index < nAngleCut; index++){
	  if( InjectionAngle < AngleCut[index] ){ ThetaIndex = index; break;}
	}
	if( InjectionAngle > AngleCut[nAngleCut-1] ){ThetaIndex  = nAngleCut;} 

	for( int iid = 0; iid < (*git).clusterIdVec().size() ; iid++){
	  handler->GetMetricPosition((*git).clusterIdVec()[iid],x1,y1); 	  
	  x1 = -x1;
	  L = sqrt( (x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));
	  t1 = (*git).clusterTimeVec()[iid];//-TimeOffset[(*git).clusterIdVec()[iid]];
	  E1  =(*git).clusterEVec()[iid];
	  R = ((x1-x0)*x0+(y1-y0)*y0)/sqrt( x0*x0+y0*y0 );
	  D = sqrt(L*L-R*R);
	  if(E1 < 20 ){ continue; }
	  if( iid < 3 && abs(R) < 25 ){
	    hisTimeZitter->Fill((*git).clusterIdVec()[iid],t1-t0);	  
	  }
	  hisTimeL->Fill( L,t1-t0 );
	  hisEL->Fill( L, E1 );
	  hisRDT->Fill(R,D,t1-t0);
	  hisRDE->Fill(R,D,E1);	    
	  
	  if( ThetaIndex < nAngleCut ){
	    hisRDT_Theta[ThetaIndex]->Fill(R,D,t1-t0);
	    hisRDE_Theta[ThetaIndex]->Fill(R,D,E1);	    
	  }

	  if( abs(D) <= 25 ){
	    hisTimeR->Fill( R,t1-t0 );
	    hisER->Fill( R, E1 );
	    if( ThetaIndex < nAngleCut ){
	    hisTimeRTheta[ThetaIndex]->Fill( R,t1-t0 );
	    hisERTheta[ThetaIndex]->Fill( R, E1 );
	    }
	  }
	  if( abs(R) <= 25 ){
	    hisTimeD->Fill( D,t1-t0 );
	    hisED->Fill( D, E1 );
	    if( ThetaIndex < nAngleCut){
	      if( D <25 ){ 
		if( R > 12.5 ){ continue;}
	      } 
	      hisTimeDTheta[ThetaIndex]->Fill( D,t1-t0 );
	      hisEDTheta[ThetaIndex]->Fill( D, E1 );	      
	    }
	  }

	  //}
      }
    }
  }


  hisTimeR->Write();
  hisER->Write();
  hisTimeD->Write();
  hisED->Write();
  hisTimeL->Write();
  hisEL->Write();
  hisTimeZitter->Write();
  hisInjectionAngle->Write();
  hisRDT->Write();
  hisRDE->Write();
  for( int i = 0; i< nAngleCut; i++){
    hisTimeRTheta[i]->Write();
    hisTimeDTheta[i]->Write();
    hisERTheta[i]->Write();
    hisEDTheta[i]->Write();
    hisRDT_Theta[i]->Write();
    hisRDE_Theta[i]->Write();
  }
  tfOut->Close();
}
