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
#include "gamma/GammaFinder.h"
#include "cluster/Cluster.h"
#include "cluster/ClusterFinder.h"
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

#include "CLHEP/Vector/ThreeVector.h"
#include "csimap/CsiMap.h"
#include "CsIPoly.h"

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
  CsiMap* map = CsiMap::getCsiMap();

  const int nFile = 1;
  //char* name = "SumUp";
  char* name = argv[1];

  Float_t OEVVetoEne;
  Float_t CC03VetoEne;
  Float_t OEVTotalVetoEne;
  Float_t CC03TotalVetoEne;
  Float_t CBARVetoEne;
  Float_t NCCVetoEne;
  Float_t CVVetoEne;

  TFile* tfin = new TFile(Form("%s",name));
  TTree* trin = (TTree*)tfin->Get("KLongRecTree");

  Double_t MaxR;
  Double_t MinX;
  Double_t MinY;
  Double_t coex;
  Double_t coey;
  Double_t vtxx;
  Double_t vtxy;
  Double_t vtxz;
  Double_t klp[3];
  Double_t klv[3];
  Double_t klE;
  Double_t MaxGChisq;
  Int_t    EventNumber;
  Double_t GPos[6][3];
  Double_t GE[6];
  Double_t GChisq[6];
  Double_t GMaxChisq;
  Double_t GTimeDelta[6];
  Double_t GTimeMaxDelta;
  Double_t GTimeMaxSigma;
  Double_t MinGDist;
  Double_t klchisqZ;
  Double_t klMass;
  Double_t GMinE;
  Double_t coeR;
  Int_t    VetoCut;
  Int_t    StdCut;
  Int_t    AdvCut;
  trin->SetBranchAddress("EventNumber"     ,&EventNumber);//EventNumber/I");
  trin->SetBranchAddress("OEVVetoEne"      ,&OEVVetoEne);//OEVVetoEne/F");
  trin->SetBranchAddress("OEVTotalVetoEne" ,&OEVTotalVetoEne);//OEVTotalVetoEne/F");
  trin->SetBranchAddress("CC03VetoEne"     ,&CC03VetoEne);//CC03VetoEne/F");
  trin->SetBranchAddress("CC03TotalVetoEne",&CC03TotalVetoEne);//CC03TotalVetoEne/F");
  trin->SetBranchAddress("NCCVetoEne"      ,&NCCVetoEne);//NCCVetoEne/F");
  trin->SetBranchAddress("CVVetoEne"       ,&CVVetoEne);//CVVetoEne/F");
  trin->SetBranchAddress("CBARVetoEne"     ,&CBARVetoEne);//CBARVetoEne/F");
  trin->SetBranchAddress("MaxR"            ,&MaxR);//MaxR/D");
  trin->SetBranchAddress("MinX"            ,&MinX);//MinX/D");
  trin->SetBranchAddress("MinY"            ,&MinY);//MinY/D");
  trin->SetBranchAddress("MinGDist"        ,&MinGDist);//MinDist/D");
  trin->SetBranchAddress("coex"            ,&coex);//coex/D");
  trin->SetBranchAddress("coey"            ,&coey);//coey/D");
  trin->SetBranchAddress("coeR"            ,&coeR);//coeR/D");
  trin->SetBranchAddress("klMass"          ,&klMass);//klMass/D");
  trin->SetBranchAddress("klchisqZ"        ,&klchisqZ);//klchisqZ/D");
  trin->SetBranchAddress("klp"             ,klp);//klp[3]/D");
  trin->SetBranchAddress("klv"             ,klv);//klv[3]/D");
  trin->SetBranchAddress("klE"             ,&klE);//klE/D");
  trin->SetBranchAddress("MaxGChisq"       ,&MaxGChisq);//MaxGChisq/D");
  trin->SetBranchAddress("GPos"            ,GPos);//GPos[6][3]/D");
  trin->SetBranchAddress("GE"              ,GE);//GE[6]/D");
  trin->SetBranchAddress("GChisq"          ,GChisq);//GChisq[6]/D");
  trin->SetBranchAddress("GTimeDelta"      ,GTimeDelta);//GTimeDelta[6]/D");
  trin->SetBranchAddress("GTimeMaxDelta"   ,&GTimeMaxDelta);//GTimeMaxDelta/D");
  trin->SetBranchAddress("GTimeMaxSigma"   ,&GTimeMaxSigma);//GTimeMaxSigma/D");
  trin->SetBranchAddress("GMinE"           ,&GMinE);//GMinE/D");
  trin->SetBranchAddress("VetoCut"         ,&VetoCut);//VetoCut/I");
  trin->SetBranchAddress("StdCut"          ,&StdCut);//StdCut/I");
  trin->SetBranchAddress("AdvCut"          ,&AdvCut);//AdvCut/I");


  TFile* tfOut = new TFile(Form("Dist_%s",name),"recreate");
  // Nocut , StdCut , StdCut+VetoCut;
  
  char* hisname[5] ={ "NoCut","StdCut","StdCutVetoCut","TimingCore","TimingHalo"};

  TH1D* hisklMass[5];
  TH1D* hisklpz[5];
  TH1D* hisklE[5];
  TH2D* hisCoe[5];
  TH1D* hisCoeY[5];
  TH1D* hisCoeX[5];
  TH1D* hisCoeR[5];
  TH1D* hisGTimeDeltaMax[5]; 
  TH1D* hisGChisqMax[5];
  TH1D* hisGE[5];
  TH1D* hisklchisqZ[5];
  TH1D* hisR[5];
  TH1D* CV[5];
  TH1D* NCC[5];
  TH1D* CBAR[5];
  TH1D* CC03[5];
  TH1D* OEV[5];

  for( int i = 0; i< 5; i++){
    hisklMass[i] = new TH1D(Form("hisklMass_%d",i),Form("hisklMass_%s",hisname[i]),200,400,600);
    hisklpz[i]   = new TH1D(Form("hisklpz_%d",i),Form("hisklPz_%s",hisname[i]),200,0,10000);
    hisklE[i]    = new TH1D(Form("hisklE_%d",i),Form("hisklE_%s",hisname[i]),200,0,10000);
    hisCoe[i]    = new TH2D(Form("hisCoe_%d",i),Form("hisCoe_%s",hisname[i]),800,-400,400,800,-400,400);
    hisCoeX[i]    = new TH1D(Form("hisCoeX_%d",i),Form("hisCoeX_%s",hisname[i]),800,-400,400);
    hisCoeY[i]    = new TH1D(Form("hisCoeY_%d",i),Form("hisCoeY_%s",hisname[i]),800,-400,400);
    hisklchisqZ[i] = new TH1D(Form("hisklchisqZ_%d",i),Form("hisklchisqZ_%s",hisname[i]),100,0,50); 
    hisGTimeDeltaMax[i] = new TH1D(Form("hisGTimeDeltaMax_%d",i),Form("hisGTimeDeltaMax_%s",hisname[i]),300,-30,30);
    hisGChisqMax[i]     = new TH1D(Form("hisGChisqMax_%d",i),Form("hisGChisqMax_%s",hisname[i]),150,0,30);
    hisGE[i]     = new TH1D(Form("hisGE_%d",i),Form("hisGE_%s",hisname[i]),300,0,3000);
    hisR[i]             = new TH1D(Form("hisR_%d",i),Form("hisR_%s",hisname[i]),1000,0,1000);
    CV[i]  = new TH1D(Form("hisCV_%d",i),Form("hisCV_%s",hisname[i]),100,0,10);
    NCC[i] = new TH1D(Form("hisNCC_%d",i),Form("hisNCC_%s",hisname[i]),50,0,100);
    CBAR[i]= new TH1D(Form("hisCBAR_%d",i),Form("hisCBAR_%s",hisname[i]),50,0,100);
    CC03[i]= new TH1D(Form("hisCC03_%d",i),Form("hisCC03_%s",hisname[i]),50,0,100);
    OEV[i] = new TH1D(Form("hisOEV_%d",i),Form("hisOEV_%s",hisname[i]),50,0,100);
    
  }
    




  for( int ievent  =0; ievent < trin->GetEntries(); ievent++){
    trin->GetEntry( ievent );
    if( (ievent % 1000) == 0){
      std::cout<< ievent <<"/" << trin->GetEntries()<< std::endl;
    }
    Double_t MaxGE = 0; 
    for( int i = 0; i< 6; i++){
      if( GE[i] > MaxGE ){
	MaxGE = GE[i];
      }
    }



    hisklMass[0]->Fill(klMass);
    hisklpz[0]->Fill(klp[2]);
    hisklchisqZ[0]->Fill(klchisqZ);
    hisklE[0]->Fill(klE);
    hisCoe[0]->Fill(coex,coey);
    hisCoeX[0]->Fill(coex);
    hisCoeY[0]->Fill(coey);
    hisGTimeDeltaMax[0]->Fill(GTimeMaxDelta);
    hisGChisqMax[0]->Fill(MaxGChisq);
    hisR[0]->Fill(coeR);
    CV[0]->Fill(CVVetoEne);
    NCC[0]->Fill(NCCVetoEne);
    CBAR[0]->Fill(CBARVetoEne);
    CC03[0]->Fill(CC03VetoEne);
    OEV[0]->Fill(OEVVetoEne);
    for( int i = 0; i< 6; i++){
      hisGE[0]->Fill(GE[i]);
    }
    if( StdCut == 0 ){
      hisklMass[1]->Fill(klMass);
      hisklpz[1]->Fill(klp[2]);
      hisklchisqZ[1]->Fill(klchisqZ);
      hisklE[1]->Fill(klE);
      hisCoe[1]->Fill(coex,coey);
      hisCoeX[1]->Fill(coex);
      hisCoeY[1]->Fill(coey);
      hisGTimeDeltaMax[1]->Fill(GTimeMaxDelta);
      hisGChisqMax[1]->Fill(MaxGChisq);
      hisR[1]->Fill(coeR);
      CV[1]->Fill(CVVetoEne);
      NCC[1]->Fill(NCCVetoEne);
      CBAR[1]->Fill(CBARVetoEne);
      CC03[1]->Fill(CC03VetoEne);
      OEV[1]->Fill(OEVVetoEne);
      for( int i = 0; i< 6; i++){
	hisGE[1]->Fill(GE[i]);
      }
   
	if( TMath::Abs(GTimeMaxDelta) < 3 ){
	  hisklMass[3]->Fill(klMass);
	  hisklpz[3]->Fill(klp[2]);
	  hisklchisqZ[3]->Fill(klchisqZ);
	  hisklE[3]->Fill(klE);
	  hisCoe[3]->Fill(coex,coey);
	  hisCoeX[3]->Fill(coex);
	  hisCoeY[3]->Fill(coey);
	  hisGTimeDeltaMax[3]->Fill(GTimeMaxDelta);
	  hisGChisqMax[3]->Fill(MaxGChisq);
	  hisR[3]->Fill(coeR);
	  CV[3]->Fill(CVVetoEne);
	  NCC[3]->Fill(NCCVetoEne);
	  CBAR[3]->Fill(CBARVetoEne);
	  CC03[3]->Fill(CC03VetoEne);
	  OEV[3]->Fill(OEVVetoEne);		  
	  for( int i = 0; i< 6; i++){
	    hisGE[3]->Fill(GE[i]);
	  }
	}else if( TMath::Abs( GTimeMaxDelta ) > 17 && TMath::Abs( GTimeMaxDelta ) <= 20 ){
	  hisklMass[4]->Fill(klMass);
	  hisklpz[4]->Fill(klp[2]);
	  hisklchisqZ[4]->Fill(klchisqZ);
	  hisklE[4]->Fill(klE);
	  hisCoe[4]->Fill(coex,coey);
	  hisCoeX[4]->Fill(coex);
	  hisCoeY[4]->Fill(coey);
	  hisGTimeDeltaMax[4]->Fill(GTimeMaxDelta);
	  hisGChisqMax[4]->Fill(MaxGChisq);
	  hisR[4]->Fill(coeR);
	  CV[4]->Fill(CVVetoEne);
	  NCC[4]->Fill(NCCVetoEne);
	  CBAR[4]->Fill(CBARVetoEne);
	  CC03[4]->Fill(CC03VetoEne);
	  OEV[4]->Fill(OEVVetoEne);	
	  for( int i = 0; i< 6; i++){
	    hisGE[4]->Fill(GE[i]);
	  }
	}



      //if( VetoCut ==0){
      if( CVVetoEne < 0.4 && OEVVetoEne < 20 && CBARVetoEne < 30 && CC03VetoEne < 30 && NCCVetoEne<20 &&MaxGE < 1500 && klchisqZ < 10){
	hisklMass[2]->Fill(klMass);
	hisklpz[2]->Fill(klp[2]);
	hisklchisqZ[2]->Fill(klchisqZ);
	hisklE[2]->Fill(klE);
	hisCoe[2]->Fill(coex,coey);
	hisCoeX[2]->Fill(coex);
	hisCoeY[2]->Fill(coey);
	hisGTimeDeltaMax[2]->Fill(GTimeMaxDelta);
	hisGChisqMax[2]->Fill(MaxGChisq);
	hisR[2]->Fill(coeR);
	CV[2]->Fill(CVVetoEne);
	NCC[2]->Fill(NCCVetoEne);
	CBAR[2]->Fill(CBARVetoEne);
	CC03[2]->Fill(CC03VetoEne);
	OEV[2]->Fill(OEVVetoEne);	
	for( int i = 0; i< 6; i++){
	  hisGE[2]->Fill(GE[i]);
	}

      }
    }
  }
  for( int i = 0; i< 5; i++){
    hisklpz[i]->Write();
    hisklchisqZ[i]->Write();
    hisklE[i]->Write();
    hisCoe[i]->Write();
    hisCoeX[i]->Write();
    hisCoeY[i]->Write();
    hisGTimeDeltaMax[i]->Write();
    hisGChisqMax[i]->Write();
    hisGE[i]->Write();
    hisR[i]->Write();
    CV[i]->Write();
    NCC[i]->Write();
    CBAR[i]->Write();
    CC03[i]->Write();
    OEV[i]->Write();
  }
  tfOut->Close();
}
