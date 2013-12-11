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
  
  TFile* tfOut = new TFile(Form("KL_Halo_%s",name),"recreate");

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



  for( int ievent  =0; ievent < tr->GetEntries(); ievent++){
    tr->GetEntry( ievent );
    if( (ievent % 1000) == 0){
      std::cout<< ievent <<"/" << tr->GetEntries()<< std::endl;
    }
    std::list<Gamma>   glist;
    //std::vector<Klong> klVec;
    std::list<Pi0> plist;
    std::vector<Klong> klVec;

    data.getData(glist);
    data.getData(plist);    
    data.getData(klVec);
    if( TMath::Abs(klVec[0].m() - KLMass ) > 10 ){ continue; }
    if( klVec[0].chisqZ() > 5 ){ continue; }

    int gIndex = 0;
    std::list<Gamma>::iterator git = glist.begin();
    std::list<Pi0>::iterator   pit = plist.begin();
    
    
    Double_t EX=0;
    Double_t EY=0;
    Double_t SumE=0;
    Double_t PX=0;
    Double_t PY=0;
    /*
    for( int i = 0; i< klVec[0].pi0().size(); i++){
      EX += klVec[0].pi0()[i].g1().e()*klVec[0].pi0()[i].g1().x();
      EY += klVec[0].pi0()[i].g1().e()*klVec[0].pi0()[i].g1().y();
      EX += klVec[0].pi0()[i].g2().e()*klVec[0].pi0()[i].g2().x();
      EY += klVec[0].pi0()[i].g2().e()*klVec[0].pi0()[i].g2().y();
      SumE+=klVec[0].pi0()[i].g1().e();
      SumE+=klVec[0].pi0()[i].g2().e();
    }
    */
    /*
    for( int i = 0; i< klVec[0].pi0().size(); i++,pit++){
      EX += (*pit).g1().e()*(*pit).g1().x();
      EY += (*pit).g1().e()*(*pit).g1().y();
      EX += (*pit).g2().e()*(*pit).g2().x();
      EY += (*pit).g2().e()*(*pit).g2().y();
      SumE+=(*pit).g1().e();
      SumE+=(*pit).g2().e();
    }
    */


    bool bGamma = false;
    for( int i = 0; i< glist.size(); i++,git++){
      if( TMath::Abs((*git).x()) < 150 && TMath::Abs((*git).y()) < 150){
	bGamma = true;
      }
      Double_t gr = TMath::Sqrt(pow( (*git).x(),2) +pow((*git).y(),2 ));
      if( gr > 800 ){
	bGamma = true;
      }
      if( i >= 6 ){ continue; }
      EX += (*git).e()*(*git).x();
      EY += (*git).e()*(*git).y();
      SumE += (*git).e();
    }
    if( bGamma ){ continue; }
    Double_t ECenterX = EX/SumE;
    Double_t ECenterY = EY/SumE;
    Double_t R = TMath::Sqrt(pow(ECenterX-5.874,2)+pow(ECenterY-1.501,2));

    bool bHalo = false;
    


    if( TMath::Abs( ECenterX - 5.984 ) > 100 || 
	TMath::Abs( ECenterY - 1.501 ) > 100 ){
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
      if( i> 6 ){ continue; }
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
    
    if( bHalo ){
      hisGammaMaxChisqHalo->Fill(maxGammaChisq);
    }else{
      hisGammaMaxChisq->Fill(maxGammaChisq);
    }

    git = glist.begin();
    Double_t baseTime = (*git).t();
    Double_t MaxTimeDelta = 0;
    Double_t TimeDelta=0;
    git++;
    for( int i = 0; i< 6; i++,git++){      
      if( MaxTimeDelta < TMath::Abs((*git).t() - baseTime )){
	MaxTimeDelta = TMath::Abs((*git).t() - baseTime);
	TimeDelta = (*git).t() - baseTime;
      }
    }
    if( bHalo ){
      hisGammaTDeltaHalo->Fill(TimeDelta);
    }else{
      hisGammaTDelta->Fill(TimeDelta);
    }



    if( maxGammaChisq > 2.5 ){ continue; }
    //std::cout<< ECenterX << "\t" << ECenterY << std::endl;
    if( klVec[0].vz() < 5000 ){
      hisR->Fill(R);
      hisKLVtx->Fill(klVec[0].vx(),klVec[0].vy());
      hisECenter->Fill( ECenterX, ECenterY);
    }
    hisXZ->Fill(klVec[0].vz(),ECenterX);
    hisYZ->Fill(klVec[0].vz(),ECenterY);
    if( TMath::Abs( ECenterX ) > 100 || TMath::Abs( ECenterY ) > 100 ){
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
  tfOut->Close();
}
