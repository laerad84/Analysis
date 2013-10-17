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
  char* name = "WAVNOCV";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"

  tf = new TFile(Form("Kl_Total_%s.root",name));
  tr = (TTree*)tf->Get(Form("trKL"));
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];

  TFile* tfOut = new TFile(Form("TimeResolution_GammaTimeDistribution_%s.root",name),"recreate");
  TTree* trOut    = new TTree("GammaTimeTree","GammaTime");
  double GammaCE[6];
  double GammaTime[6];
  double GammaSignal[6];
  double GammaTimeSigma;
  double GammaTimeSigmaExcept[6];
  double GammaE[6];
  double GammaX[6];
  double GammaY[6];
  double KLMass;
  double KLChisq;
  double KLChisqSec;
  double Pi0Pt[3];
  double MaxHeight;
  double TLowGammaMean;
  double TLowGammaSigma;
  int    nLowGamma;
  trOut->Branch("MaxHeight",&MaxHeight,"MaxHeight/D");
  trOut->Branch("nLowGamma",&nLowGamma,"nLowGamma/I");
  trOut->Branch("TLowGammaMean",&TLowGammaMean,"TLowGammaMean/D");
  trOut->Branch("TLowGammaSigma",&TLowGammaSigma,"TLowGammaSigma/D");
  trOut->Branch("GammaCE",GammaCE,"GammaCE[6]/D");
  trOut->Branch("GammaTime",GammaTime,"GammaTime[6]/D");
  trOut->Branch("GammaTimeSigma",&GammaTimeSigma,"GammaTimeSigma/D");
  trOut->Branch("GammaTimeSigmaExcept",GammaTimeSigmaExcept,"GammaTimeSigmaExcept[6]/D");
  trOut->Branch("GammaE",GammaE,"GammaE[6]/D");
  trOut->Branch("GammaSignal",GammaSignal,"GammaSignal[6]/D");
  trOut->Branch("GammaX",GammaX,"GammaX[6]/D");
  trOut->Branch("GammaY",GammaY,"GammaY[6]/D");
  trOut->Branch("KLMass",&KLMass,"KLMass/D");
  trOut->Branch("KLChisq",&KLChisq,"KLChisq/D");
  trOut->Branch("KLChisqSec",&KLChisqSec,"KLChisqSec/D");
  trOut->Branch("Pi0Pt",Pi0Pt,"Pi0Pt[3]/D");

  TH2D* hisResolution = new TH2D("hisResolution","hisResolution",100,0,400,100,-20,20);
  TH2D* hisResolutionAdj = new TH2D("hisResolutionAdj","hisResolutionAdj",100,0,400,100,-20,20);
  TH2D* hisResolutionLY[4];//0: good/good 1:good/bad 2:bad/good 3:bad/bad
  TH2D* hisX = new TH2D("hisX","hisX",2716,0,2716,80,-1000,1000);
  TH2D* hisY = new TH2D("hisY","hisY",2716,0,2716,80,-1000,1000);
  for( int i = 0; i< 4; i++){
    hisResolutionLY[i] = new TH2D(Form("hisResolutionLY_%d",i),Form("hisResolutionLY_%d",i),100,0,400,100,-20,20);
  }
  TH2D* hisTimingNonlinearity[3];
  char* Name[3] = {"Good","Bad","Large"};
  for( int i = 0; i< 3; i++){
    hisTimingNonlinearity[i] = new TH2D(Form("hisTimingNonlinearity_%d",i),Form("hisTimingNonLinearity_%s",Name[i]),160,0,16000,100,-20,20);
  }

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
  TH1D* hisGammaTimeDistribution[nAngleCut];
  TH1D* hisGammaTime0[nAngleCut];

  for( int i = 0; i< nAngleCut; i++){
    hisGammaTimeDistribution[i] = new TH1D(Form("hisGTimeDist_%d",i),Form("hisGTimeDist_%d",i),-160,-80,80);
    hisGammaTime0[i] = new TH1D(Form("hisGammaTime0_%d",i),Form("hisGammaTime0_%d",i),200,0,400);
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
  int CsiNumber;
  double CsiSignal[3000];
  int CsiModID[3000];
  //double KLMass;
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiSignal",CsiSignal);
  tr->SetBranchAddress("CsiModID",CsiModID);
  //tr->SetBranchAddress("KLMass",&KLMass);
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
  std::ifstream ifs("TimeOffset_ShowerHeight_10.dat");
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
    //std::cout<< ievent << std::endl;
    //if( ievent  >= 100000 ){ break ; } 
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    data.getData(clist);
    data.getData(glist);
    data.getData(klVec);
    

    GammaTimeSigma=0;
    KLChisq = 0;
    KLChisqSec = -1;
    KLMass  = 0; 
    for( int i = 0; i< 3; i++){
      Pi0Pt[i] = 0;
    }
    for( int i = 0; i< 6; i++){
      GammaTime[i] = 0;
      GammaY[i] = 0;
      GammaX[i] = 0;
      GammaE[i] = 0;
      GammaCE[i] = 0;
      GammaTimeSigmaExcept[i] = 0;
    }
    
    MaxHeight = 0;
    nLowGamma = 0;
    TLowGammaMean = 0;
    TLowGammaSigma = 0;





    double klpos[3]={0};
    klpos[0] = klVec[0].vx();
    klpos[1] = klVec[1].vy();
    klpos[2] = klVec[2].vz();
    //if( klVec[0].chisqZ() > 15 ){ continue; }
    //if( klpos[2] > 5500 ){continue; }
    Double_t x0(0),y0(0);
    Double_t x1(0),y1(0);
    Double_t t0(0),t1(0); 
    Double_t L(0);
    Double_t E0(0),E1(0);
    Double_t R(0),D(0);
    Double_t InjectionAngle={0};
    std::list<Gamma>::iterator git = glist.begin();

    Double_t GammaTime0=-100;
    bool bAbortIDEvent = false;
    int gIndex =0; 
    GammaTime0 = 0; 
    for( git = glist.begin();
	 git != glist.end();
	 git++){
      if( gIndex >= 6 ){ break;}
      
      GammaTime0 += (*git).clusterTimeVec()[0]-showerTimeDelayAdj(klVec[0],(*git))-TimeOffset[(*git).clusterIdVec()[0]]-gammaLOF(klVec[0],(*git))/sol;
      gIndex++;
    }
    GammaTime0 /= 6;
    gIndex = 0; 

    for( git = glist.begin();
	 git != glist.end();
	 git++){

      if( gIndex >= 6 ){ break; }
      //if( abs((*git).x()) > 600 || abs((*git).y() > 600 )){ continue; } 
      /*
	if( git == glist.begin() ){
	  GammaTime0 = (*git).t()-showerTimeDelayAdj(klVec[0],(*git))-TimeOffset[(*git).clusterIdVec()[0]];
	}
      */
	GammaTime[gIndex] = (*git).clusterTimeVec()[0]-showerTimeDelayAdj(klVec[0],(*git))-TimeOffset[(*git).clusterIdVec()[0]]-gammaLOF(klVec[0],(*git))/sol;;
	GammaCE[gIndex]   = (*git).clusterEVec()[0];
	GammaE[gIndex]    = (*git).e();
	GammaX[gIndex]    = (*git).x();
	GammaY[gIndex]    = (*git).y();
	GammaTimeSigma    += TMath::Power(GammaTime[gIndex] -GammaTime0,2);
	int tmpID  = (*git).clusterIdVec()[0];
	for( int i = 0; i< CsiNumber; i++){
	  if( tmpID == CsiModID[i] ){
	    GammaSignal[gIndex]= CsiSignal[i];
	    break;
	  }
	}
	if( GammaSignal[gIndex] > MaxHeight ){ MaxHeight = GammaSignal[gIndex];}
	if( GammaSignal[gIndex] < 1000 &&GammaSignal[gIndex] > 100){ 
	  TLowGammaMean += GammaTime[gIndex];
	  TLowGammaSigma += GammaTime[gIndex]*GammaTime[gIndex];
	  nLowGamma++;
	}
	gIndex++;
    }
    TLowGammaMean /= nLowGamma;
    TLowGammaSigma = TMath::Sqrt(TLowGammaSigma/nLowGamma - TLowGammaMean*TLowGammaMean);

    int PositionIndex = 0; 
    if( nLowGamma >= 3 && GammaTimeSigma < 4){
      for( int i = 0; i< 6; i++){
	if( GammaSignal[i] >=1000 && TLowGammaSigma < 4){
	  if( TMath::Abs(GammaX[i]) < 600 ){
	    if( GammaY[i] > 0  ){
	      if( GammaX[i] < -100 ){
		PositionIndex = 1;
	      }else{
		PositionIndex = 0;
	      }
	    }else{
	      if( GammaX[i] < 75 ){
		PositionIndex = 1;
	      }else{
		PositionIndex = 0;
	      }
	    }
	  }else{
	    PositionIndex = 2;
	  }
	  hisTimingNonlinearity[PositionIndex]->Fill(GammaSignal[i], GammaTime[i]-TLowGammaMean );
	}
      }
    }
    for( int i = 0; i< 6; i++){
      for( int j = 0; j< 6; j++){
	if( i == j){continue;}
	GammaTimeSigmaExcept[i] +=TMath::Power(GammaTime[j] -GammaTime[i],2);
      }
      GammaTimeSigmaExcept[i] = TMath::Sqrt(GammaTimeSigmaExcept[i]/5);
    }

    bool gammaMinE = false;
    bool gammapos  = false;
    for( int i = 0; i< 6;i++){
      if( GammaE[i] < 50 ){gammaMinE = true;}
      if( abs(GammaX[i]) < 150 && abs(GammaY[i]) < 150 ){ gammapos = true;}
      if( abs( GammaY[i] ) > 550 ){ gammapos = true; }
      if( sqrt(GammaX[i]*GammaX[i]+GammaY[i]*GammaY[i]) > 850 ){ gammapos = true; }
    }
    if( gammaMinE ){ continue; }
    if( gammapos  ){ continue; }
    GammaTimeSigma =TMath::Sqrt(GammaTimeSigma/5);
    KLChisq = klVec[0].chisqZ();
    KLMass  = klVec[0].m();
    if( klVec.size() == 2 ){
      KLChisqSec = klVec[1].chisqZ(); 
    }
    for( int i = 0; i< 3; i++){
      Pi0Pt[i] = klVec[0].pi0()[i].p3().perp();
    }
    trOut->Fill();




    git = glist.begin();
    int    g0crystalID = (*git).clusterIdVec()[0];
    if( g0crystalID >= 2240 ){continue; }
    double g0time   = (*git).clusterTimeVec()[0];
    double g0length = TMath::Sqrt( TMath::Power(((*git).x() - klVec[0].vx()),2) 
				 + TMath::Power(((*git).y() - klVec[0].vy()),2) 
				 + TMath::Power(((*git).z() - klVec[0].vz()),2));
    double g0Offset = TimeOffset[g0crystalID];//man->GetT0Offset(g0crystalID);
    double g0Shower = showerTimeDelayAdj(klVec[0],(*git));
    double g0Delta  = g0Offset-g0length/sol-g0Shower;//man->GetT0Offset(g0crystalID);

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

    git = glist.begin();    
    double g0Ene = (*git).clusterEVec()[0];
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
  

  hisTimingNonlinearity[0]->Write();
  hisTimingNonlinearity[1]->Write();
  hisTimingNonlinearity[2]->Write();
  trOut->Write();
  tfOut->Close();
}
