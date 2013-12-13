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
  TFile* tf;
  TTree* tr;
  //char* name = "SumUp";
  char* name = argv[1];
  tf = new TFile(Form("%s",name));
  tr = (TTree*)tf->Get(Form("RecTree"));
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  Float_t OEVVetoEne;
  Float_t CC03VetoEne;
  Float_t OEVTotalVetoEne;
  Float_t CC03TotalVetoEne;
  tr->SetBranchAddress("OEVVetoEne",&OEVVetoEne);
  tr->SetBranchAddress("CC03VetoEne",&CC03VetoEne);
  tr->SetBranchAddress("OEVTotalVetoEne",&OEVTotalVetoEne);
  tr->SetBranchAddress("CC03TotalVetoEne",&CC03TotalVetoEne);
  

  TFile* tfOut = new TFile(Form("KL_Halo_%s",name),"recreate");
  TTree* trOut = new TTree("KLongRecTree","KLongRecTree");

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
  Double_t GTimeDelta[6];
  Double_t GTimeMaxDelta;
  Double_t GTimeMaxSigma;
  Double_t MinGDist;
  Double_t klchisqZ;
  Double_t klMass;
  Double_t GMinE;
  trOut->Branch("EventNumber"     ,&EventNumber,"EventNumber/I");
  trOut->Branch("OEVVetoEne"      ,&OEVVetoEne,"OEVVetoEne/F");
  trOut->Branch("OEVTotalVetoEne" ,&OEVTotalVetoEne,"OEVTotalVetoEne/F");
  trOut->Branch("CC03VetoEne"     ,&CC03VetoEne,"CC03VetoEne/F");
  trOut->Branch("CC03TotalVetoEne",&CC03TotalVetoEne,"CC03TotalVetoEne/F");
  trOut->Branch("MaxR"            ,&MaxR,"MaxR/D");
  trOut->Branch("MinX"            ,&MinX,"MinX/D");
  trOut->Branch("MinY"            ,&MinY,"MinY/D");
  trOut->Branch("MinGDist"        ,&MinGDist,"MinDist/D");
  trOut->Branch("coex"            ,&coex,"coex/D");
  trOut->Branch("coey"            ,&coey,"coey/D");
  trOut->Branch("klMass"          ,&klMass,"klMass/D");
  trOut->Branch("klchisqZ"        ,&klchisqZ,"klchisqZ/D");
  trOut->Branch("klp"             ,klp,"klp[3]/D");
  trOut->Branch("klv"             ,klv,"klv[3]/D");
  trOut->Branch("klE"             ,&klE,"klE/D");
  trOut->Branch("MaxGChisq"       ,&MaxGChisq,"MaxGChisq/D");
  trOut->Branch("GPos"            ,GPos,"GPos[6][3]/D");
  trOut->Branch("GE"              ,GE,"GE[6]/D");
  trOut->Branch("GChisq"          ,GChisq,"GChisq[6]/D");
  trOut->Branch("GTimeDelta"      ,GTimeDelta,"GTimeDelta[6]/D");
  trOut->Branch("GTimeMaxDelta"   ,&GTimeMaxDelta,"GTimeMaxDelta/D");
  trOut->Branch("GTimeMaxSigma"   ,&GTimeMaxSigma,"GTimeMaxSigma/D");
  trOut->Branch("GMinE"           ,&GMinE,"GMinE/D");

  TH2D* hisECenter = new TH2D("hisECenter","hisECenter",800,-400,400,800,-400,400);
  TH2D* hisKLVtx   = new TH2D("hisKLVtx","hisKLVtx",800,-400,400,800,-400,400);
  TH2D* hisPt      = new TH2D("hisPt","hisPt",800,-400,400,800,-400,400);
  TH1D* hisR       = new TH1D("hisR","hisR",800,0,400);
  TH2D* hisXZ      = new TH2D("hisXZ","hisXZ",600,0,6000,800,-400,400);
  TH2D* hisYZ      = new TH2D("hisYZ","hisYZ",600,0,6000,800,-400,400);
  TH1D* hisKLP     = new TH1D("hisKLP","hisKLP",500,0,5000);
  TH1D* hisKLPHalo = new TH1D("hisKLPHalo","hisKLPHalo",500,0,5000);
  TH2D* hisGammaPos= new TH2D("hisGammaPos","hisGammaPos",80,-1000,1000,80,-1000,1000);
  TH2D* hisGammaPosHalo= new TH2D("hisGammaPosHalo","hisGammaPosHalo",80,-1000,1000,80,-1000,1000);

  TH1D* hisGammaR= new TH1D("hisGammaR","hisGammaR",80,0,1000);
  TH1D* hisGammaRHalo= new TH1D("hisGammaRHalo","hisGammaRHalo",80,0,1000);
  TH1D* hisGammaChisq= new TH1D("hisGammaChisq","hisGammaChisq",200,0,50);
  TH1D* hisGammaChisqHalo= new TH1D("hisGammaChisqHalo","hisGammaChisqHalo",200,0,50);
  TH1D* hisGammaMaxChisq= new TH1D("hisGammaMaxChisq","hisGammaChisq",200,0,50);
  TH1D* hisGammaMaxChisqHalo= new TH1D("hisGammaMaxChisqHalo","hisGammaChisqHalo",200,0,50);
  TH1D* hisGammaTDelta= new TH1D("hisGammaTDelta","hisGammaTDelta",200,-20,20);
  TH1D* hisGammaTDeltaHalo = new TH1D("hisGammaTDeltaHalo","hisGammaTDeltaHalo",200,-20,20);

  TH1D* hisCC03ETHalo = new TH1D("hisCC03ETHalo","hisCC03ETHalo",200,0,25);
  TH1D* hisCC03ET     = new TH1D("hisCC03ET","hisCC03ET",200,0,25);
  TH1D* hisOEVETHalo  = new TH1D("hisOEVETHalo","hisOEVETHalo",200,0,25);
  TH1D* hisOEVET      = new TH1D("hisOEVET","hisOEVET",200,0,25);

  for( int ievent  =0; ievent < tr->GetEntries(); ievent++){
    tr->GetEntry( ievent );
    if( (ievent % 1000) == 0){
      std::cout<< ievent <<"/" << tr->GetEntries()<< std::endl;
    }
    EventNumber = ievent;
    std::list<Gamma>   glist;
    //std::vector<Klong> klVec;
    std::list<Pi0> plist;
    std::vector<Klong> klVec;

    data.getData(glist);
    data.getData(plist);    
    data.getData(klVec);
    std::list<Gamma>::iterator git = glist.begin();
    std::list<Pi0>::iterator   pit = plist.begin();


    Double_t FlightTime[6];
    Double_t baseTime = 0;
    Double_t MaxTimeDelta = 0;
    Double_t TimeDelta=0;
    MaxR = 0;
    MinX = 1000;
    MinY = 1000;
    GMinE = 100000;
    MaxGChisq     = 0;
    GTimeMaxDelta = 0;
    GTimeMaxSigma = 0;
    MinGDist = 0;
    coex = 0;
    coey = 0;
    Double_t SumE = 0;
    int gIndex = 0;        
    Double_t EX=0;
    Double_t EY=0;
    Double_t PX=0;
    Double_t PY=0;    
    for( int i = 0; i< 6; i++,git++){
      double l = TMath::Sqrt( TMath::Power((*git).x()-klVec[0].vx(),2)
			      +TMath::Power((*git).y()-klVec[0].vy(),2)
			      +TMath::Power((*git).z()-klVec[0].vz(),2));
      FlightTime[i] = l/sol;
      baseTime+=(*git).t()-FlightTime[i];
      coex      +=(*git).e()*(*git).x();
      coey      +=(*git).e()*(*git).y();
      SumE      +=(*git).e();
      if( (*git).e() < GMinE ){
	GMinE = (*git).e();
      }
    }

    
    baseTime = baseTime/6;    
    coex = coex/SumE;
    coey = coey/SumE;

    git = glist.begin();
    for( int i = 0; i< 6; i++,git++){
      GPos[i][0]=(*git).x();
      GPos[i][1]=(*git).y();
      GPos[i][2]=(*git).z();
      GE[i]     =(*git).e();
      GChisq[i] =(*git).chisq();
      if( GChisq[i] > MaxGChisq ){
	MaxGChisq =GChisq[i];
      }
      GTimeDelta[i] = (*git).t() -baseTime-FlightTime[i];
      if( TMath::Abs(GTimeMaxDelta) < TMath::Abs((*git).t() -baseTime-FlightTime[i]) ){
	GTimeMaxDelta= (*git).t() -baseTime-FlightTime[i];
      }

      Double_t R = TMath::Sqrt(TMath::Power((*git).x(),2)+TMath::Power((*git).y(),2));
      if( R > MaxR){
	MaxR = R;
      }
      if( TMath::Abs((*git).x())< MinX ){
	MinX  = TMath::Abs((*git).x());
      }
      if( TMath::Abs((*git).y())< MinY ){
	MinY  = TMath::Abs((*git).y());	
      }
    }
    for( int i = 0; i< 6; i++){
      Double_t tmpSigma = 0;
      for( int j = 0; j<6; j++){
	if( i==j){ continue; }
	tmpSigma+= TMath::Power(GTimeDelta[i]-GTimeDelta[j],2);
      }
      tmpSigma = TMath::Sqrt( tmpSigma /5 );
      if( tmpSigma > GTimeMaxSigma ){
	GTimeMaxSigma = tmpSigma ;
      }
    }




    klv[0] = klVec[0].vx();
    klv[1] = klVec[0].vy();
    klv[2] = klVec[0].vz();
    klp[0] = klVec[0].p3()[0];
    klp[1] = klVec[0].p3()[1];
    klp[2] = klVec[0].p3()[2];
    klchisqZ=klVec[0].chisqZ();
    klE    = klVec[0].e();
    klMass = klVec[0].m();
    trOut->Fill();

    if( TMath::Abs(klVec[0].m() - KLMass ) > 10 ){ continue; }
    if( klVec[0].chisqZ() > 5 ){ continue; }
    
    git = glist.begin();
    bool bGamma = false;
    for( int i = 0; i< 6; i++,git++){
      if( TMath::Abs((*git).x()) < 125 && TMath::Abs((*git).y()) < 125){
	bGamma = true;
      }
      Double_t gr = TMath::Sqrt(pow( (*git).x(),2) +pow((*git).y(),2 ));
      if( gr > 850 ){
	bGamma = true;
      }
    }

    if( bGamma ){ continue; }
    //if( CC03TotalVetoEne > 1.5 ){ continue; }
    //if( OEVTotalVetoEne > 1.5 ){ continue; }
    if( GMinE < 200 ){ continue; }
    Double_t R = TMath::Sqrt(pow(coex-5.874,2)+pow(coey-1.501,2));    
    bool bHalo = false;    
    if( TMath::Abs( coex - 5.984 ) > 100 || 
	TMath::Abs( coey - 1.501 ) > 100 ){
      bHalo = true;
    }  

    if( bHalo ){
      hisKLPHalo->Fill( klVec[0].p3().mag());
    }else{
      hisKLP->Fill( klVec[0].p3().mag());
    }

    git = glist.begin();
    double maxGammaChisq=0;
    for( int i = 0; i< glist.size();i++,git++){
      if( i>= 6 ){ break; }
      Double_t gr = TMath::Sqrt(pow( (*git).x(),2) +pow((*git).y(),2 ));
      if( maxGammaChisq < (*git).chisq()){
	maxGammaChisq = (*git).chisq();
      }
      if( bHalo ){
	hisGammaPosHalo->Fill((*git).x(),(*git).y());
	hisGammaRHalo->Fill(gr);
	hisGammaChisqHalo->Fill((*git).chisq());
      }else{
	hisGammaPos->Fill((*git).x(),(*git).y());
	hisGammaR->Fill(gr);
	hisGammaChisq->Fill((*git).chisq());
      }
    }

    if( TMath::Abs(GTimeMaxDelta) > 3 ){ continue; }
    if( maxGammaChisq > 2.5 ){ continue; }
    if( bHalo ){
      hisGammaTDeltaHalo->Fill(TimeDelta);
    }else{
      hisGammaTDelta->Fill(TimeDelta);
    }
    if( bHalo ){
      hisGammaMaxChisqHalo->Fill(maxGammaChisq);
      hisOEVETHalo->Fill( OEVVetoEne );
      hisCC03ETHalo->Fill( CC03VetoEne );
    }else{
      hisGammaMaxChisq->Fill(maxGammaChisq);
      hisOEVET->Fill( OEVVetoEne );
      hisCC03ET->Fill( CC03VetoEne );
    }

    //std::cout<< coex << "\t" << coey << std::endl;
    if( klVec[0].vz() < 5500 ){
      hisR->Fill(R);
      hisKLVtx->Fill(klVec[0].vx(),klVec[0].vy());
      hisECenter->Fill( coex, coey);
    }
    hisXZ->Fill(klVec[0].vz(),coex);
    hisYZ->Fill(klVec[0].vz(),coey);
    if(bHalo){
      CsIPoly* csi = new CsIPoly(Form("CsI_%d",ievent),Form("CsI_%d",ievent));
      git = glist.begin();
      for( int i = 0; i< 6; i++,git++){
	for( int j = 0; j< (*git).clusterIdVec().size();j++){
	  csi->Fill((*git).clusterIdVec()[j],(*git).clusterEVec()[j]);
	}
      }
      csi->Write();
    }
  }

  hisOEVET->Write();
  hisOEVETHalo->Write();
  hisCC03ET->Write();
  hisCC03ETHalo->Write();
  hisGammaTDelta->Write();
  hisGammaTDeltaHalo->Write();
  hisGammaChisq->Write();
  hisGammaChisqHalo->Write();
  hisGammaMaxChisq->Write();
  hisGammaMaxChisqHalo->Write();
  hisGammaPos->Write();
  hisGammaPosHalo->Write();
  hisGammaR->Write();
  hisGammaRHalo->Write();
  hisKLP->Write();
  hisKLPHalo->Write();
  hisXZ->Write();
  hisYZ->Write();
  hisR->Write();
  hisKLVtx->Write();
  hisECenter->Write();
  trOut->Write();
  tfOut->Close();
}
