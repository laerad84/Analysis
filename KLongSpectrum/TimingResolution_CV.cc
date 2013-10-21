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
  
  TF1* TimingDelayFunc = new TF1("TimingDelayFunc",TimingHeightDelay,0,20000,3);
  TimingDelayFunc->SetParameters(0,0.03566,1621);

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

  TFile* tfOut = new TFile(Form("TimeDistribution_%s.root",name),"recreate");

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

  TH1D* hisTimingDeltaDist[3];
  TH2D* hisEnergyTimingDist[3];
  TH2D* hisTimingCutDist[2];
  TH2D* hisTimingCutAllDist[2];
  TH1D* hisTimingRMS   = new TH1D("hisTimingRMS","hisTimingRMS",400,0,50);
  TH1D* hisTimingMaxDelta = new TH1D("hisTimingMaxDelta","hisTimingMaxDelta",400,-50,50);
  TH2D* hisKLMassvsRMS = new TH2D("hisKLMassvsRMS","hisKLMassvsRMS",200,400,600,200,0,25);
  TH1D* hisKLMassRMS[2];
  TH1D* hisKLMassPeak[2];
  TH2D* hisKLMassDelta = new TH2D("hisKLMassDelta","hisKLMassDelta",200,400,600,400,-40,40);
  for( int i = 0; i< 2; i++){
    hisKLMassRMS[i]= new TH1D(Form("hisKLMassRMS_%d",i),Form("hisKLMassRMS_%d",i),200,400,600);
    hisKLMassPeak[i] = new TH1D(Form("hisKLMassPeak_%d",i),Form("hisKLMassPeak_%d",i),200,0,50);
  }
  
  for( int i = 0; i< 3; i++){
    hisTimingDeltaDist[i] = new TH1D(Form("hisTimingDeltaDist_%d",i),Form("hisTimingDeltaDist_%d",i),400,-50,50);
    hisEnergyTimingDist[i] = new TH2D(Form("hisEnergyTimingDist_%d",i),Form("hisEnergyTimingDist_%d",i),200,0,2000,400,-50,50);
  }
  for( int i = 0; i< 2; i++){
    hisTimingCutDist[i] = new TH2D(Form("hisTimingCutDist_%d",i),Form("hisTimingCutDist_%d",i),200,0,2000,400,-50,50);
    hisTimingCutAllDist[i] = new TH2D(Form("hisTimingCutAllDist_%d",i),Form("hisTimingCutAllDist_%d",i),200,0,2000,400,-50,50);
  }

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
      ShowerOffset[igamma] = showerTimeDelayAdj( klVec[0],(*git));      
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
      AllOffset[i]    = TimeOffset[GammaID[i]] + TOFOffset[i]+ShowerOffset[i]+HeightOffset[i]; 
    }

    double MinTimeRMS = 9000000;
    double MeanTiming = 0;
    int    OffMeanID[2]={-1,-1};
    
    for( int igamma  = 0; igamma < 5; igamma++){
      for( int jgamma = igamma+1; jgamma < 6; jgamma++){
	double tmpRMS=0;
	double tmpMean=0;

	for( int kgamma = 0; kgamma < 6; kgamma++){
	  if( kgamma == igamma || kgamma == jgamma ){ continue; }
	  tmpRMS += TMath::Power(GammaTime[kgamma]-AllOffset[kgamma],2);
	  tmpMean+= GammaTime[kgamma]-AllOffset[kgamma];
	}
	tmpMean/=4;
	tmpRMS=sqrt(tmpRMS/4-tmpMean*tmpMean);
	if(tmpRMS < MinTimeRMS ){
	  MinTimeRMS = tmpRMS;
	  MeanTiming = tmpMean;
	  OffMeanID[0]=igamma;
	  OffMeanID[1]=jgamma;
	}
      }
    }

    for( int igamma = 1; igamma < 6; igamma++){
      hisTimingDeltaDist[0]->Fill(GammaTime[igamma]-TimeOffset[GammaID[igamma]] -GammaTime[0]+TimeOffset[GammaID[0]]);
      hisTimingDeltaDist[1]->Fill(GammaTime[igamma]-TimeOffset[GammaID[igamma]]-TOFOffset[igamma]-ShowerOffset[igamma]
				  -GammaTime[0]+TimeOffset[GammaID[0]]+TOFOffset[0]+ShowerOffset[0]);
      hisTimingDeltaDist[2]->Fill(GammaTime[igamma]-AllOffset[igamma]-GammaTime[0]+AllOffset[0]);      
      hisEnergyTimingDist[0]->Fill(GammaEnergy[0],GammaTime[igamma]-TimeOffset[GammaID[igamma]]-GammaTime[0]+TimeOffset[GammaID[0]]);
      hisEnergyTimingDist[1]->Fill(GammaEnergy[0],GammaTime[igamma]-TimeOffset[GammaID[igamma]]-TOFOffset[igamma]-ShowerOffset[igamma]
				   -GammaTime[0]+TimeOffset[GammaID[0]]+TOFOffset[0]+ShowerOffset[0]);
      hisEnergyTimingDist[2]->Fill(GammaEnergy[0],GammaTime[igamma]-AllOffset[igamma]-GammaTime[0]+AllOffset[0]);
    }
    for( int igamma = 0; igamma < 6; igamma++){
      bool bLY = LYRegion(GammaX[igamma],GammaY[igamma]);
      if( bLY ){
	hisTimingCutAllDist[0]->Fill(GammaEnergy[igamma],GammaTime[igamma]-AllOffset[igamma]-MeanTiming);
      }else{
	hisTimingCutAllDist[1]->Fill(GammaEnergy[igamma],GammaTime[igamma]-AllOffset[igamma]-MeanTiming);
      }
    }
    if( TMath::Abs(GammaTime[OffMeanID[0]]-AllOffset[OffMeanID[0]]-MeanTiming) >  TMath::Abs(GammaTime[OffMeanID[1]]-AllOffset[OffMeanID[1]]-MeanTiming)){
      hisTimingMaxDelta->Fill(GammaTime[OffMeanID[0]]-AllOffset[OffMeanID[0]]-MeanTiming);
      hisTimingCutDist[0]->Fill(GammaEnergy[OffMeanID[0]],GammaTime[OffMeanID[0]]-AllOffset[OffMeanID[0]]-MeanTiming);
      hisTimingCutDist[1]->Fill(GammaEnergy[OffMeanID[1]],GammaTime[OffMeanID[1]]-AllOffset[OffMeanID[1]]-MeanTiming);
      hisKLMassDelta->Fill(klVec[0].m(),GammaTime[OffMeanID[1]]-AllOffset[OffMeanID[1]]-MeanTiming);
    }else{
      hisTimingMaxDelta->Fill(GammaTime[OffMeanID[1]]-AllOffset[OffMeanID[1]]-MeanTiming);
      hisTimingCutDist[1]->Fill(GammaEnergy[OffMeanID[0]],GammaTime[OffMeanID[0]]-AllOffset[OffMeanID[0]]-MeanTiming);
      hisTimingCutDist[0]->Fill(GammaEnergy[OffMeanID[1]],GammaTime[OffMeanID[1]]-AllOffset[OffMeanID[1]]-MeanTiming);
    }
    hisTimingRMS->Fill(MinTimeRMS);
    hisKLMassvsRMS->Fill(klVec[0].m(),MinTimeRMS);

    if( MinTimeRMS < 5 ){
      hisKLMassRMS[0]->Fill(klVec[0].m());
    }else{
      hisKLMassRMS[1]->Fill(klVec[0].m());
    }
    if( TMath::Abs(klVec[0].m()-KLMass) < 10){
      hisKLMassPeak[0]->Fill(MinTimeRMS);
    }else{
      hisKLMassPeak[1]->Fill(MinTimeRMS);
    }    
  }
std::cout<< "Write" << std::endl;
 hisTimingRMS->Write();
 hisTimingMaxDelta->Write();
 hisKLMassvsRMS->Write();
 hisKLMassDelta->Write();
  for( int i = 0; i< 3; i++){
    hisTimingDeltaDist[i]->Write();
    hisEnergyTimingDist[i]->Write();
  }
  for( int i = 0; i< 2; i++){
    hisTimingCutAllDist[i]->Write();
    hisTimingCutDist[i]->Write();
    hisKLMassRMS[i]->Write();
    hisKLMassPeak[i]->Write();
  }

  tfOut->Close();
}
