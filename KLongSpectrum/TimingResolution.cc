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
#include "TSpline.h"
#include "TGraphErrors.h"
const double KLMass = 497.648;//MeV

double KLSpectrum(double* x,double*p){
  //double klp   = TMath::Sqrt( x[0]*x[0]-KLMass*KLMass);  
  double klp = x[0];
  double sigma = p[2]*(1-(p[3]+p[4]*klp)*(klp-p[1]));
  double value = p[0]*TMath::Exp(-0.5*(klp-p[1])*(klp-p[1])/(sigma*sigma));
  return value;
}

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
  double cosTheta  = TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  double delayTime = 500/solc - depth*(1-cosTheta)/sol;
  return delayTime;
}
/*
double showerTimeDelayAdj(Klong kl, Gamma g){
  double clusterTimingConstant = 0.541812;
  double depth     = CsIX0*(Pcor[0]);
  double cosTheta  = TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  double delayTime = 500/solc - depth*(1-cosTheta)/sol + 0.541812*TMath::Log10(TMath::E())*log(g.e()/1000);
  return delayTime;
}
*/
double showerTimeDelayAdj(Klong kl, Gamma g){
  double clusterTimingConstant = 0.541812;
  double depth     = showerDepth( g.e() );
  double cosTheta  = TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  //double delayTime = 500/solc - depth*(1-cosTheta)/sol + 0.541812*TMath::Log10(TMath::E())*log(g.e()/1000);
  double delayTime = depth*(1/sol-cosTheta/solc);
  return delayTime;
}

double ThetaConstant(Klong kl, Gamma g){
  double cosTheta  =1-sol/solc*TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  return cosTheta;
}

double gammaLOF( Klong kl, Gamma g){
  double length = 0;
  length = sqrt( pow(g.x()-kl.vx(),2)+ pow(g.y()-kl.vy(),2)+pow(g.z()-kl.vz(),2));
  return length;
}

bool LYRegion( double x, double y ){
  bool bRegion = false;
  if( y > 0 ){
    if( x > -75 ){
      bRegion = true; 
    }else{
      bRegion = false; 
    }
  }else{
    if( x > 100 ){
      bRegion = true; 
    }else{
      bRegion = false;
    }
  }
  return bRegion;
}
/*
TFile* tflin = new TFile("TimingLinearityFuncLaser.root");
TGraphErrors* grTimingLinearity = (TGraphErrors*)tflin->Get("grTimingLinearityLaser");
TSpline3* spl = new TSpline3("spl",grTimingLinearity);
double HeightAdjFunc(double* x, double* par){
  return spl->Eval(x[0]);
}
*/
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
  char* name = "DATA_NONTIMECAL";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"

  tf = new TFile(Form("Kl_Total_%s.root",name));
  tr = (TTree*)tf->Get(Form("trKL"));
  Int_t    CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  int      CsiNumber;
  double   CsiSignal[3000];
  int      CsiModID[3000];
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiModID",CsiModID);//CsiNumber



  TFile* tfOut = new TFile(Form("TimeResolution_%s_Height.root",name),"recreate");
  TH2D* hisResolution = new TH2D("hisResolution","hisResolution",100,0,400,100,-20,20);
  TH2D* hisResolutionAdj = new TH2D("hisResolutionAdj","hisResolutionAdj",100,0,400,100,-20,20);
  TH2D* hisResolutionLY[2];//0: good/good 1:good/bad 2:bad/good 3:bad/bad
  TH2D* hisResolutionHeight[2];
  TH2D* hisResolutionLY_Neighbor[2]; 
  TH2D* hisX = new TH2D("hisX","hisX",2716,0,2716,80,-1000,1000);
  TH2D* hisY = new TH2D("hisY","hisY",2716,0,2716,80,-1000,1000);
  TH2D* hisEnergyDep = new TH2D("hisEnergyDep","hisEnergyDep",100,0,4000,100,-20,20);
  TH2D* hisEnergyDepAdj = new TH2D("hisEnergyDepAdj","hisEnergyDepAdj",100,0,4000,100,-20,20);
  TH2D* hisCEnergyDep = new TH2D("hisCEnergyDep","hisCEnergyDep",100,0,4000,100,-20,20);
  TH2D* hisCEnergyDepAdj = new TH2D("hisCEnergyDepAdj","hisCEnergyDepAdj",100,0,4000,100,-20,20);

  TH2D* hisEnergyDepLY[2];
  TH2D* hisEnergyDepLYAdj[2];
  TH2D* hisCEnergyDepLY[2];
  TH2D* hisCEnergyDepLYAdj[2];
  TH2D* hisEnergyRatioLY[2];
  TH2D* hisEnergyRatioLYLog[2];
  TH2D* hisEnergyRatioLYAdj[2];
  TH2D* hisEnergyRatioEnergyLY[2];
  TH2D* hisEnergyRatioEnergy1LY[2];
  TH2D* hisCEnergyRatioLY[2];
  TH2D* hisCEnergyRatioEnergyLY[2];
  TH2D* hisCEnergyRatioEnergy1LY[2];
  TH2D* hisDeltaXLY[2];
  TH2D* hisDeltaTime[2];
  TH2D* hisThetaDep[2];
  TH2D* hisThetaDepAdj[2];
  TH2D* hisThetaCDep[2];
  TH2D* hisThetaCDepAdj[2];
  TH2D* hisThetaCDep1[2];
  TH2D* hisThetaCDep1Adj[2];
  TH2D* hisThetaDep1[2];
  TH2D* hisThetaDep1Adj[2];
  for( int i = 0; i< 2; i++){
    hisResolutionLY[i] = new TH2D(Form("hisResolutionLY_%d",i),Form("hisResolutionLY_%d",i),100,0,400,100,-20,20);
    hisResolutionHeight[i] = new TH2D(Form("hisResolutionHeight_%d",i),Form("hisResolutionHeight_%d",i),160,0,16000,100,-20,20);

    hisResolutionLY_Neighbor[i] = new TH2D(Form("hisResolutionLY_Neighbor_%d",i),Form("hisResolutionLY_Neighbor_%d",i),100,0,400,100,-20,20);
    hisEnergyDepLY[i] = new TH2D(Form("hisEnergyDepLY_%d",i),Form("hisEnergyDepLY_%d",i),100,0,400,100,-20,20);
    hisCEnergyDepLY[i] = new TH2D(Form("hisCEnergyDepLY_%d",i),Form("hisCEnergyDepLY_%d",i),100,0,400,100,-20,20);
    hisEnergyDepLYAdj[i] = new TH2D(Form("hisEnergyDepLYAdj_%d",i),Form("hisEnergyDepAdjLY_%d",i),100,0,400,100,-20,20);
    hisCEnergyDepLYAdj[i] = new TH2D(Form("hisCEnergyDepLYAdj_%d",i),Form("hisCEnergyDepAdjLY_%d",i),100,0,400,100,-20,20);
    hisEnergyRatioLY[i] = new TH2D(Form("hisEnergyRatioLY_%d",i),Form("hisEnergyRatioLY_%d",i),100,0,1,100,-20,20);
    hisEnergyRatioLYAdj[i] = new TH2D(Form("hisEnergyRatioLYAdj_%d",i),Form("hisEnergyRatioLYAdj_%d",i),100,0,1,100,-20,20);
    hisEnergyRatioLYLog[i] = new TH2D(Form("hisEnergyRatioLYLog_%d",i),Form("hisEnergyRatioLYLog_%d",i),100,-1,0,100,-20,20);
    hisEnergyRatioEnergyLY[i] = new TH2D(Form("hisEnergyRatioEnergyLY_%d",i),Form("hisEnergyRatioEnergyLY_%d",i),100,0,1,100,0,1000);
    hisEnergyRatioEnergy1LY[i] = new TH2D(Form("hisEnergyRatioEnergy1LY_%d",i),Form("hisEnergyRatioEnergy1LY_%d",i),100,0,1,100,0,1000);					 
    hisCEnergyRatioLY[i] = new TH2D(Form("hisCEnergyRatioLY_%d",i),Form("hisCEnergyRatioLY_%d",i),100,0,1,100,-20,20);
    hisCEnergyRatioEnergyLY[i] = new TH2D(Form("hisCEnergyRatioEnergyLY_%d",i),Form("hisCEnergyRatioEnergyLY_%d",i),100,0,1,100,-20,20);
    hisCEnergyRatioEnergy1LY[i] = new TH2D(Form("hisCEnergyRatioEnergy1LY_%d",i),Form("hisCEnergyRatioEnergy1LY_%d",i),100,0,1,100,-20,20);					 
    hisDeltaXLY[i] = new TH2D(Form("hisDeltaXLY_%d",i),Form("hisDeltaXLY_%d",i),100,0,1,200,-200,200);
    hisDeltaTime[i] = new TH2D(Form("hisDeltaTime_%d",i),Form("hisDeltaTime_%d",i),200,-200,200,100,-20,20);
    hisThetaDep[i] = new TH2D(Form("hisThetaDep_%d",i),Form("hisThetaDep_%d",i),800,-20,20,100,-20,20);
    hisThetaDepAdj[i] = new TH2D(Form("hisThetaDepAdj_%d",i),Form("hisThetaDepAdj_%d",i),800,-20,20,100,-20,20);
    hisThetaDep1[i] = new TH2D(Form("hisThetaDep1_%d",i),Form("hisThetaDep1_%d",i),800,-20,20,100,-20,20);
    hisThetaDep1Adj[i] = new TH2D(Form("hisThetaDep1Adj_%d",i),Form("hisThetaDep1Adj_%d",i),800,-20,20,100,-20,20);
    hisThetaCDep[i] = new TH2D(Form("hisThetaCDep_%d",i),Form("hisThetaCDep_%d",i),800,-2,2,100,-20,20);
    hisThetaCDepAdj[i] = new TH2D(Form("hisThetaCDepAdj_%d",i),Form("hisThetaCDepAdj_%d",i),800,-2,2,100,-20,20);
    hisThetaCDep1[i] = new TH2D(Form("hisThetaCDep1_%d",i),Form("hisThetaCDep1_%d",i),800,-2,2,100,-20,20);
    hisThetaCDep1Adj[i] = new TH2D(Form("hisThetaCDep1Adj_%d",i),Form("hisThetaCDep1Adj_%d",i),800,-2,2,100,-20,20);
  }

  Double_t TimeOffset[2716]={0};
  std::ifstream ifs("TimeOffset_ShowerHeight_15.dat");
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
    if(CsiNumber > 500 ){ continue; }
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
    if( klpos[2] > 5000 ){ continue; }

    std::list<Gamma>::iterator git = glist.begin();
    for( int igamma  = 0; igamma < 6; igamma++,git++){
      bool bggood =false;
      int id0= (*git).clusterIdVec()[0];
      int id1= (*git).clusterIdVec()[1];
      int nFit = 0;

      double h0=0;
      double h1=0;
      double t0=0;
      double t1=0;
      for( int i = 0; i< CsiNumber; i++){	
	if( id0 == CsiModID[i] ){
	  h0= CsiSignal[i];
	  if( id0 > id1 ){
	    break;
	  }else{
	    continue;
	  }
	}
	if( id1 == CsiModID[i]){
	  h1= CsiSignal[i];
	  if( id1 > id0 ){ 
	    break;
	  }else{
	    continue;
	  }
	}
      }
      double x0 = (*git).coex();
      double y0 = (*git).coey();
      //if( TMath::Abs( x0 ) < 20 && TMath::Abs(y0) < 20 ){ continue; }
      //if( TMath::Abs( x0 ) > 50 && TMath::Abs(y0) > 50 ){ continue; }

      double gxy[2];
      gxy[0] = (*git).x();
      gxy[1] = (*git).y();
      bggood = LYRegion( gxy[0], gxy[1] );

      if( h1/h0 > 0.9 && h0/h1 > 0.9 ){
	if( bggood ){
	  hisResolutionHeight[0]->Fill(h0,(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0])-t0+t1);
	}else{
	  hisResolutionHeight[1]->Fill(h0,(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0])-t0+t1);
	}
      }
      
      
      if((*git).clusterEVec()[0] > 300 && (*git).clusterEVec()[0] < 400 ){
	hisResolution->Fill( (*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-((*git).clusterTimeVec()[0])); 
	hisResolutionAdj->Fill( (*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0])-t0+t1);
	if( bggood ){
	  hisResolutionLY[0]->Fill((*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0])+t0-t1);
	  
	}else{
	  hisResolutionLY[1]->Fill((*git).clusterEVec()[1],(*git).clusterTimeVec()[1]-TimeOffset[id1]-((*git).clusterTimeVec()[0]-TimeOffset[id0])+t0-t1);
	}
      }
      
      for( int iid = 0; iid < (*git).clusterIdVec().size()-1; iid++){
	if( iid > 3 ){ continue; } 
	int TestID[2]={-1};
	double x[2]={0};
	double y[2]={0}; 
	double R[2]={0};
	double D[2]={0};
	double t[2]={0};
	double L[2]={0};
	double E[2]={0};
	double T[2]={0};
	double h[2]={0};	
	TestID[0]=(*git).clusterIdVec()[iid];
	for( int i = 0; i< CsiNumber; i++){	
	  if( TestID[0] == CsiModID[i] ){
	    h[0]= CsiSignal[i];
	    continue;
	  }
	}
	handler->GetMetricPosition((*git).clusterIdVec()[iid],x[0],y[0]);	
	x[0] = -x[0];
	L[0] = sqrt( (x0-x[0])*(x0-x[0])+(y0-y[0])*(y0-y[0]));
	R[0] = ((x[0]-x0)*x0+(y[0]-y0)*y0)/sqrt( x0*x0+y0*y0 );
	D[0] = sqrt(L[0]*L[0]-R[0]*R[0]);
	E[0] = (*git).clusterEVec()[iid];	  
	T[0] = (*git).clusterTimeVec()[iid];
	//for( int jid = iid+1; jid < iid+3; jid++){
	for( int jid = iid+1; jid < (*git).clusterIdVec().size(); jid++){
	  if(jid>4){ break;}
	  TestID[1]=(*git).clusterIdVec()[jid];
	  for( int i = 0; i< CsiNumber; i++){	
	    if( TestID[1] == CsiModID[i] ){
	      h[1]= CsiSignal[i];
	      continue;
	    }
	  }
	  if( h[1] < 10 || h[0] < 10 ){ continue; }

	  handler->GetMetricPosition((*git).clusterIdVec()[jid],x[1],y[1]);
	  x[1] = -x[1];
	  L[1] = sqrt( (x0-x[1])*(x0-x[1])+(y0-y[1])*(y0-y[1]));
	  R[1] = ((x[1]-x0)*x0+(y[1]-y0)*y0)/sqrt( x0*x0+y0*y0 );
	  //R[1] = sqrt( pow( x[0]-x[1] ,2 )+ pow( y[0] - y[1] ,2) );
	  D[1] = sqrt(L[1]*L[1]-R[1]*R[1]);
	  E[1] = (*git).clusterEVec()[jid];	  
	  T[1] = (*git).clusterTimeVec()[jid];	  
	  //if( R[1] > 26 ){ continue; }


	  if( TMath::Abs( R[0]-R[1] ) >10 ){ continue;}
	  if( TMath::Abs(D[0]) > 20 || TMath::Abs(D[1]) > 20 ){ continue; }
	  
	  if( TMath::Abs(R[0]) > 7 || TMath::Abs(R[1]) > 7 ){ continue; }
	  //if( TMath::Abs(R[0]-R[1]) > 12.5 ){ continue; }

	  if( h[0]/h[1] > 0.9 || h[1]/h[0] > 0.9 ){
	    if( bggood ){
	      hisResolutionHeight[0]->Fill(h[0],T[1]-TimeOffset[TestID[1]]-(T[0]-TimeOffset[TestID[0]]));
	    }else{
	      hisResolutionHeight[1]->Fill(h[0],T[1]-TimeOffset[TestID[1]]-(T[0]-TimeOffset[TestID[0]]));
	    }
	  }

	  if( E[1] < E[0]*0.85 || E[1] > E[0]*1.15 ){ continue; }
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
    hisResolutionHeight[i]->Write();
  }

  hisEnergyDep->Write();
  hisCEnergyDep->Write();
  hisEnergyDepAdj->Write();
  hisCEnergyDepAdj->Write();
  for( int i = 0; i < 2; i++){
    hisEnergyDepLY[i]->Write();
    hisEnergyDepLYAdj[i]->Write();
    hisCEnergyDepLY[i]->Write();
    hisCEnergyDepLYAdj[i]->Write();
    hisEnergyRatioLY[i]->Write();
    hisEnergyRatioLYLog[i]->Write();
    hisEnergyRatioLYAdj[i]->Write();
    hisEnergyRatioEnergyLY[i]->Write();
    hisEnergyRatioEnergy1LY[i]->Write();
    hisCEnergyRatioLY[i]->Write();
    hisCEnergyRatioEnergyLY[i]->Write();
    hisCEnergyRatioEnergy1LY[i]->Write();
    hisDeltaXLY[i]->Write();
    hisDeltaTime[i]->Write();


    hisThetaDep[i]->Write();
    hisThetaDepAdj[i]->Write();
    hisThetaDep1[i]->Write();
    hisThetaDep1Adj[i]->Write();
    hisThetaCDep[i]->Write();
    hisThetaCDepAdj[i]->Write();
    hisThetaCDep1[i]->Write();
    hisThetaCDep1Adj[i]->Write();
  }
  hisX->Write();
  hisY->Write();

  tfOut->Close();
}
