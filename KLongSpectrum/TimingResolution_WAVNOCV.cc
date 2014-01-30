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

double TimingHeightDelay( double* x, double *p ){
  double value  = p[0] + TMath::Log(1+p[1]*TMath::Exp(x[0]/p[2]));
  return value; 
}

double clusterTimeDelay( double * x , double* p ){
  double value  = p[0] + p[1]*TMath::Log(1+p[2]*TMath::Exp(x[0]/2000));
  return value; 
}


int main( int argc, char** argv){
  TF1* cDelayFunc = new TF1("cDelayFunc",clusterTimeDelay,0,20000,3);
  cDelayFunc->SetParameters(0,1.491,0.04457);
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
  char* name = "WAVNOCV";//"DATA_NONTIMECALNOCV";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"
  
  TF1* TimingDelayFunc = new TF1("TimingDelayFunc",TimingHeightDelay,0,20000,3);
  TimingDelayFunc->SetParameters(0,0.03566,1621);

  tf = new TFile(Form("Kl_Total_%s.root",name));
  tr = (TTree*)tf->Get(Form("trKL"));
  Int_t    CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  int      CsiNumber;
  double   CsiSignal[3000];
  int      CsiModID[3000];
  int s_arrsize = 120;
  int   GamClusNumbers;
  double   GamClusSizes[120];
  double   GamClusCsiSignal[120][120];
  
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiModID",CsiModID);//CsiNumber
  tr->SetBranchAddress("GamClusNumbers",&GamClusNumbers);
  tr->SetBranchAddress("GamClusSizes",GamClusSizes);//GamClusNumbers
  tr->SetBranchAddress("GamClusCsiSignal",GamClusCsiSignal);//GamClusNumbers

  TFile* tfOut = new TFile(Form("TimeDistribution_%s.root",name),"recreate");

  Double_t TimeOffset[2716]={0};
  for( int i = 0; i < 2716; i++){
    TimeOffset[i] = 0;
  }
  /*
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
  */

  //0:  Offset+ TOF , 1: ShowerTime
  TH2D* hisGammaDeltaTime[2];
  TH2D* hisGamClusDeltaTime[2];
  for( int i = 0; i< 2; i++){
    hisGammaDeltaTime[i]  = new TH2D(Form("hisGammaDeltaTime_%d",i),Form("hisGammaDeltaTime_%d",i),40,0,20000,400,-20,20);
    hisGamClusDeltaTime[i]= new TH2D(Form("hisGamClusDeltaTime_%d",i),Form("hisGamClusDeltaTime_%d",i),40,0,20000,400,-20,20);
  }
  
  TTree* trOut = new TTree("TimeTree","Timetree");
  int    GammaID[6]; 
  double GammaTime[6];
  double GammaHeight[6];
  double TOFOffset[6];
  double ShowerOffset[6];
  double HeightOffset[6];
  double AllOffset[6];
  double GammaX[6];
  double GammaY[6];
  double GammaEnergy[6];
  double MeanTimeDelta[6];
  trOut->Branch("GammaID"     ,GammaID,"GammaID[6]/I");
  trOut->Branch("GammaTime"   ,GammaTime,"GammaTime[6]/D");
  trOut->Branch("GammaHeight" ,GammaHeight,"GammaHeight[6]/D");
  trOut->Branch("TOFOffset"   ,TOFOffset,"TOFOffset[6]/D");
  trOut->Branch("ShowerOffset",ShowerOffset,"ShowerOffset[6]/D");
  trOut->Branch("GammaEnergy" ,GammaEnergy,"GammaEnergy[6]/D");
  trOut->Branch("HeightOffset",HeightOffset,"HeightOffset[6]/D");
  trOut->Branch("MeanTimeDelta",MeanTimeDelta,"MeanTimeDelta[6]/D");


  for( int ievent = 0; ievent < tr->GetEntries(); ievent++){      
  //for( int ievent = 0; ievent < 1500000; ievent++){      
    tr->GetEntry(ievent);
    if(CsiNumber > 500 ){ continue; }
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    data.getData(clist);
    data.getData(glist);
    data.getData(klVec);
    double klpos[3]={0,0,0};
    klpos[0] = klVec[0].vx();
    klpos[1] = klVec[0].vy();
    klpos[2] = klVec[0].vz();

    if( glist.size() < 6 ){ continue;}
    if( klVec[0].chisqZ() > 10 ){ continue; }
    if( klpos[2] > 5000 ){ continue; }

    for( int i = 0; i< 6; i++){
      GammaID[i] = -1;
      GammaHeight[i] = 0;
      TOFOffset[i] = 0;      
      GammaTime[i] = 0; 
      TOFOffset[i] = 0;
      ShowerOffset[i] = 0;
      HeightOffset[i] = 0;
      AllOffset[i] = 0;
      GammaX[i] = 0;
      GammaY[i] = 0;
      GammaEnergy[i] = 0; 
      MeanTimeDelta[i] = 0;
    }
   
    std::list<Gamma>::iterator git = glist.begin();
    for( int igamma  = 0; igamma < 6; igamma++,git++){
      if( git == glist.end()){ break; }
      GammaTime[igamma]         = (*git).clusterTimeVec()[0];
      GammaID[igamma]           = (*git).clusterIdVec()[0];
      GammaX[igamma]            = (*git).coex();
      GammaY[igamma]            = (*git).coey();
      GammaEnergy[igamma]       = (*git).e();
      TOFOffset[igamma]         = gammaLOF(klVec[0],(*git))/sol;
      ShowerOffset[igamma]      = showerTimeDelayAdj( klVec[0],(*git));      
    }

    //GammaPositionCut //
    bool bGammaPosition=true;
    for( int i = 0; i <6; i++){
      if( TMath::Abs(GammaX[i]) < 150 && TMath::Abs(GammaY[i]) < 150 ){
	bGammaPosition = false;
      }
      if(GammaX[i]*GammaX[i] +GammaY[i]*GammaY[i] > 850*850 ){
	bGammaPosition = false;
      }
      if( TMath::Abs(GammaY[i]) > 550 ){ 
	bGammaPosition = false;
      }      
    }
    if( !bGammaPosition){ continue; }


    int gIndex = 0;
    for( int i = 0; i< klVec[0].pi0().size(); i++){
      for( int j = 1; j< klVec[0].pi0()[i].g1().clusterIdVec().size(); j++){
	if( j >= 120 ){ break; }
	if( GamClusCsiSignal[gIndex][j] < 2000 && GamClusCsiSignal[gIndex][j] > 1000 ){
	  hisGamClusDeltaTime[0]->Fill( GamClusCsiSignal[gIndex][j], klVec[0].pi0()[i].g1().clusterTimeVec()[j] - klVec[0].pi0()[i].g1().clusterTimeVec()[0]);
	  hisGamClusDeltaTime[1]->Fill( GamClusCsiSignal[gIndex][j], klVec[0].pi0()[i].g1().clusterTimeVec()[i] - klVec[0].pi0()[i].g1().clusterTimeVec()[0] - cDelayFunc->Eval(GamClusCsiSignal[gIndex][0]) +cDelayFunc->Eval(GamClusCsiSignal[gIndex][0]));
	}
      }
      gIndex++;
      for( int j = 1; j< klVec[0].pi0()[i].g2().clusterIdVec().size(); j++){
	if( j >= 120 ){ break; }
	if( GamClusCsiSignal[gIndex][j] < 2000 && GamClusCsiSignal[gIndex][j] > 1000 ){
	  hisGamClusDeltaTime[0]->Fill( GamClusCsiSignal[gIndex][j], klVec[0].pi0()[i].g2().clusterTimeVec()[j] - klVec[0].pi0()[i].g2().clusterTimeVec()[0]);
	  hisGamClusDeltaTime[1]->Fill( GamClusCsiSignal[gIndex][j], klVec[0].pi0()[i].g2().clusterTimeVec()[i] - klVec[0].pi0()[i].g2().clusterTimeVec()[0] - cDelayFunc->Eval(GamClusCsiSignal[gIndex][0]) +cDelayFunc->Eval(GamClusCsiSignal[gIndex][0]));
	}
      }
      gIndex++;
    }

    int nsg = 0; 
    for( int i = 0; i< CsiNumber; i++){
      for( int j = 0; j< 6; j++){
	if( GammaID[j] == CsiModID[i] ){
	  GammaHeight[j] = CsiSignal[i]; 
	  nsg++;
	}
      }
      if( nsg == 6 ){ break; }
    }
    
    if(nsg < 6 ){ std::cout<<ievent << ":??" << std::endl;}
    for( int i = 0; i< 6; i++){
      HeightOffset[i] = TimingDelayFunc->Eval(GammaHeight[i]);
      AllOffset[i]    = TimeOffset[GammaID[i]] + TOFOffset[i];
    }

    for( int i = 0; i < 6; i++){
      for( int j = 0; j<6; j++){
	if( i == j ){ continue; }
	MeanTimeDelta[i] += GammaTime[j]-TOFOffset[j];
      }
      MeanTimeDelta[i] = MeanTimeDelta[i]/5;
    }
    trOut->Fill();
  }
  trOut->Write();
  tfOut->Close();
}
