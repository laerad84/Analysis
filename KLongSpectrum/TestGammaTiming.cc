#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector2.h"
#include "TVector.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TF1.h"

#include "gnana/E14GNAnaDataContainer.h"
#include "gnana/E14GNAnaFunction.h"
#include "gamma/Gamma.h"
#include "gamma/GammaFinder.h"
#include "cluster/Cluster.h"
#include "cluster/ClusterFinder.h"
#include "pi0/Pi0.h"
#include "klong/Klong.h" 
#include "csimap/CsiMap.h"

#include <vector>
#include <list>
#include <cstdlib>
#include <cstdio>
#include <string>

#include "CLHEP/Vector/ThreeVector.h"

#include "GeneralFunctions.h"
#include "LocalFunction.h"
#include "User_Function.h"
#include "User_Functions.h"
#include "CrateIDHandler.h"
#include "CsIPoly.h"

const double KLMass = 497.648;//MeV
double const sol = 299.792458;//[mm/nsec]
double const solc= 80;//[mm/nsec]
double const Pcor[2]={6.49003,0.99254};
double const CsIX0=18.5;//mm

double gammaLOF( Klong kl, Gamma g){
  double length = 0;
  length = sqrt( pow(g.x()-kl.vx(),2)+ pow(g.y()-kl.vy(),2)+pow(g.z()-kl.vz(),2));
  return length;
}
double THFunction( double *x, double *p){
  double value=0;
  value = p[0]+p[1]*TMath::Log( 1 + p[2]*TMath::Exp(x[0]/2000));
  return value;
}

int main( int argc, char** argv){
  TFile* tf;
  TTree* tr;
  char* name = "WAVNOCV";
  ClusterFinder cFinder;
  GammaFinder   gFinder;
  TF1* THCorrFunc = new TF1("THCorrFunc",THCorrectionFunction,0,25000,3);
  THCorrFunc->SetParameters(0,1.672,0.0319651);
  CrateIDHandler * CIDHandler = new CrateIDHandler();

  tf = new TFile(Form("Kl_Total_%s.root",name));
  tr = (TTree*)tf->Get("trKL");
  
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  Int_t CsiNumber;
  Int_t CsiID[3000]    = {-1};
  Double_t CsiSignal[3000]= {0};
  Double_t CsiTime[3000]  = {0};
  Double_t CsiEne[3000]   = {0}; 
  Int_t GamClusNumbers;
  Int_t GamClusSizes[120];//GamClusNumbers
  Double_t GamClusCsiSignal[120][120];//GamClusNumbers
  Double_t GamClusCsiCrate[120][120];//GamClusNumbers
  E14GNAnaDataContainer data;
  //data.setBranchAddress( tr );
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiModID",CsiID);//CsiID
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  tr->SetBranchAddress("CsiTime",CsiTime);//CsiNumber

  //tr->SetBranchAddress("GamClusNumbers",&GamClusNumbers);
  //tr->SetBranchAddress("GamClusCsiSignal",GamClusCsiSignal);//GamClusNumbers
  //tr->SetBranchAddress("GamClusCsiCrate",GamClusCsiCrate);//GamClusNumbers


  TFile* tfOut = new TFile(Form("kl_Total_%s_GammaTime.root",name),"recreate");
  TTree* trOut = new TTree("trKL","GammaTime Adj");
  Double_t GammaCenterTime[6];
  Double_t GammaWeight[6];
  Double_t GammaTOF[6];
  Double_t GammaTotE[6];
  Double_t GammaTSigma[6];
  Double_t GammaTDelta[6];
  Double_t GammaTDeltaMax;
  Double_t EventPeakTime;
  Int_t GammaTDeltaID;
  TH1D*    TimeHist = new TH1D("TimeHist","TimeHist",150,0,300);
  trOut->Branch("CsiL1nTrig",&CsiL1nTrig,"CsiL1nTrig/I");
  trOut->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsiL1TrigCount[20]/D");
  trOut->Branch("GammaCenterTime",GammaCenterTime,"GammaCenterTime[6]/D");
  trOut->Branch("GammaWeight",GammaWeight,"GammaWeight[6]/D");
  trOut->Branch("GammaTOF",GammaTOF,"GammaTOF[6]/D");
  trOut->Branch("GammaTotE",GammaTotE,"GammaTotE[6]/D");
  trOut->Branch("GammaTSigma",GammaTSigma,"GammaTSigma[6]/D");
  trOut->Branch("GammaTDelta",&GammaTDelta,"GammaTDelta[6]/D");
  trOut->Branch("GammaTDeltaMax",&GammaTDeltaMax,"GammaTDeltaMax/D");
  trOut->Branch("GammaTDeltaID",&GammaTDeltaID,"GammaTDeltaID/I");
  trOut->Branch("EventPeakTime",&EventPeakTime,"EventPeakTime/D");
  trOut->Branch("GamClusNumbers",&GamClusNumbers,"GamClusNumbers/I");
  trOut->Branch("GamClusSizes",GamClusSizes,"GamClusSizes[GamClusNumbers]/I");
  trOut->Branch("GmaClusCsiSignal",GamClusCsiSignal,"GamClusCsiSignal[GamClusNumbers][120]/D");

  const int s_arrSize = 120;
  trOut->Branch("GamClusNumbers",&GamClusNumbers,"GamClusNumbers/I");
  trOut->Branch("GamClusSizes",GamClusSizes,"GamClusSizes[GamClusNumbers]/I");  
  trOut->Branch("GamClusCsiSignal",GamClusCsiSignal,Form("GamClusCsiSignal[GamClusNumbers][%d]/D",s_arrSize));
  trOut->Branch("GamClusCsiCrate",GamClusCsiCrate,Form("GamClusCsiCrate[GamClusNumbers][%d]/I",s_arrSize));
  data.branchOfKlong( trOut );
      
  Int_t    CsiTNumber;
  Int_t CsiTCutID[2716];
  Double_t CsiTCutE[2716];
  Double_t CsiTCutT[2716];
  Double_t CsiTCutCID[2716];
  Double_t CsiTCutSig[2716];

  TH1D*  TimeDeltaHist = new TH1D("TimeDeltaHist","TimeDeltaHist",300,-150,150);
  TH2D*  TimeDeltaHistEne = new TH2D("TimeDeltaHistEne","TimeDeltaHistEne",500,0,2000,300,-150,150);
  TH2D*  TimeDeltaHistEneID[2];
  for( int i =0; i< 2; i++){
    TimeDeltaHistEneID[i] = new TH2D(Form("TimeDeltaHistEneID_%d",i),Form("TimeDeltaHistEneID_%d",i),500,0,2000,300,-150,150);
  }
  TH2D*  TimeDeltaHistEneCrate[20];
  TH2D*  TimeDeltaHistSigCrate[20];
  for( int i =0; i< 20; i++){
    TimeDeltaHistEneCrate[i] = new TH2D(Form("TimeDeltaHistEneCrate_%d",i),Form("TimeDeltaHistEneCrate_%d",i),500,0,2000,300,-150,150);
    TimeDeltaHistSigCrate[i] = new TH2D(Form("TimeDeltaHistSigCrate_%d",i),Form("TimeDeltaHistSigCrate_%d",i),500,0,20000,300,-150,150);
  }
  TH2D*  hisIDDelay = new TH2D("hisIDDelay","hisIDDelay",2716,0,2716,20,-10,10);
  CsIPoly* poly = new CsIPoly("csi","csi");

  for( int ievt = 0; ievt< tr->GetEntries(); ievt++){
  //for( int ievt = 0; ievt< 100; ievt++){
    tr->GetEntry(ievt);
    std::vector<Klong> klVec;
    std::list<Cluster> clist;
    std::list<Gamma>   glist;    
    //data.getData(klVec);
    TimeHist->Reset();
    GammaTDeltaMax = 0;
    for( int ig = 0; ig< 6; ig++){
      GammaCenterTime[ig] = 0;
      GammaWeight[ig] = 0;
      GammaTOF[ig] = 0;
      GammaTotE[ig] = 0;      
      GammaTSigma[ig] = 0;
      GammaTDelta[ig] = 0;
    }
    for( int i = 0; i< CsiNumber; i++){
      CsiTime[i] = CsiTime[i] - THCorrFunc->Eval(CsiSignal[i]);
      if( CsiSignal[i] < 8000 ){ 
	TimeHist->Fill( CsiTime[i] ,CsiEne[i] );
      }
    }
    Double_t MaximumBinCenter = TimeHist->GetBinCenter(TimeHist->GetMaximumBin());
    EventPeakTime = MaximumBinCenter;
    CsiTNumber = 0;
    for( int i =0; i< 2716; i++){
      CsiTCutID[i] = 0;
      CsiTCutE[i]  = 0;
      CsiTCutT[i]  = 0;
      CsiTCutCID[i]= 0;
      CsiTCutSig[i]= 0;
    }
    for( int i = 0; i< CsiNumber; i++){
      if( CsiTime[i] < 60 || CsiTime[i] > 200 ){
	continue;
      }
      if( TMath::Abs(CsiTime[i] - EventPeakTime ) >  10){ continue; }
      CsiTCutT[CsiTNumber] = CsiTime[i];
      CsiTCutE[CsiTNumber] = CsiEne[i];
      CsiTCutID[CsiTNumber]= CsiID[i];
      CsiTCutCID[CsiTNumber]= CIDHandler->GetCrate(CsiID[i]);
      CsiTCutSig[CsiTNumber]= CsiSignal[i];
      CsiTNumber++;

      TimeDeltaHist->Fill( CsiTime[i] - MaximumBinCenter );
      if( CsiID[i] < 2240 ){
	TimeDeltaHistEneID[0]->Fill( CsiEne[i], CsiTime[i]-MaximumBinCenter );
      }else{
	TimeDeltaHistEneID[1]->Fill( CsiEne[i], CsiTime[i]-MaximumBinCenter );
      }      
      TimeDeltaHistEne->Fill( CsiEne[i], CsiTime[i]-MaximumBinCenter );
      TimeDeltaHistEneCrate[CIDHandler->GetCrate(CsiID[i])]->Fill( CsiEne[i], CsiTime[i]-MaximumBinCenter );
      TimeDeltaHistSigCrate[CIDHandler->GetCrate(CsiID[i])]->Fill( CsiSignal[i], CsiTime[i]-MaximumBinCenter );
      if( CsiSignal[i] > 8000 ){
	hisIDDelay->Fill( CsiID[i], CsiTime[i]-MaximumBinCenter);
	poly->Fill( CsiID[i]);
      }
    }
    
    //clist = cFinder.findCluster( CsiNumber, CsiID, CsiEne, CsiTime);
    clist = cFinder.findCluster( CsiTNumber, CsiTCutID, CsiTCutE, CsiTCutT);
    gFinder.findGamma( clist, glist );
    if( clist.size() < 6 ){ continue; }
    if( glist.size() !=6 ){ continue; }
    if( !(User_RecG6(glist,klVec))){ continue; }
    user_cut( data,klVec);
    int igamma = 0;
    for( int ipi = 0; ipi < klVec[0].pi0().size(); ipi++){
      double tmpTime[2];
      tmpTime[0] = klVec[0].pi0()[ipi].g1().t();
      tmpTime[1] = klVec[0].pi0()[ipi].g2().t();
      SetGammaTime( klVec[0].pi0()[ipi].g1());
      SetGammaTime( klVec[0].pi0()[ipi].g2());
      /*
      std::cout<< tmpTime[0] << "\t" <<   klVec[0].pi0()[ipi].g1().t() << std::endl;
      std::cout<< tmpTime[1] << "\t" <<   klVec[0].pi0()[ipi].g2().t() << std::endl;
      */
      GammaCenterTime[igamma] = GetTiming( klVec[0].pi0()[ipi].g1());
      GammaWeight[igamma]     = GetWeight( klVec[0].pi0()[ipi].g1());
      GammaTOF[igamma]        = CalGammaTOF( klVec[0],klVec[0].pi0()[ipi].g1());
      GammaTotE[igamma]       = klVec[0].pi0()[ipi].g1().e();
      GammaTSigma[igamma]     = GetClusterTSigma( klVec[0].pi0()[ipi].g1());
      igamma++;
      GammaCenterTime[igamma] = klVec[0].pi0()[ipi].g2().t();
      GammaWeight[igamma]     = GetWeight( klVec[0].pi0()[ipi].g2());
      GammaTOF[igamma]        = CalGammaTOF( klVec[0],klVec[0].pi0()[ipi].g2());
      GammaTotE[igamma]       = klVec[0].pi0()[ipi].g2().e();
      GammaTSigma[igamma]     = GetClusterTSigma( klVec[0].pi0()[ipi].g2());
      igamma++;
    }
    /*
    for( int ig = 0; ig < igamma-1; ig++){
      for( int jg  =ig+1; jg < igamma; jg++){
	if( TMath::Abs(GammaCenterTime[ig]-GammaTOF[ig] - (GammaCenterTime[jg] - GammaTOF[jg])) > GammaTDeltaMax){
	  GammaTDeltaMax = TMath::Abs(GammaCenterTime[ig]-GammaTOF[ig] - (GammaCenterTime[jg] - GammaTOF[jg]));
	}
      }
    }
    */
    GammaTDeltaMax = 0;
    for( int ig = 0; ig < igamma; ig++){
      for( int jg = 0; jg < igamma; jg++){
	if( ig != jg ){
	  GammaTDelta[ig] += GammaCenterTime[ig]-GammaTOF[ig] - (GammaCenterTime[jg] - GammaTOF[jg]);
	}
	GammaTDelta[ig] = GammaTDelta[ig]/5.;
	if( TMath::Abs(GammaTDelta[ig]) > GammaTDeltaMax){
	  GammaTDeltaMax = TMath::Abs(GammaTDelta[ig]);
	}
      }
    }


    data.setData(klVec);
    int clNumber = 0;
    GamClusNumbers = glist.size();
    for( int pt = 0; pt < klVec[0].pi0().size(); pt++){
      GamClusSizes[clNumber] = klVec[0].pi0()[pt].g1().clusterIdVec().size();
      for( int i = 0; i < klVec[0].pi0()[pt].g1().clusterIdVec().size(); i++){
	for( int in = 0; in < CsiNumber; in++){
	  if( CsiID[in] == klVec[0].pi0()[pt].g1().clusterIdVec()[i]){
	    GamClusCsiSignal[clNumber][i] = CsiTCutSig[in];
	    GamClusCsiCrate[clNumber][i]  = CsiTCutCID[in];
	    break;
	  }
	}
      }
      clNumber++;
      GamClusSizes[clNumber] = klVec[0].pi0()[pt].g2().clusterIdVec().size();
      for( int i = 0; i < klVec[0].pi0()[pt].g2().clusterIdVec().size(); i++){
	for( int in = 0; in < CsiNumber; in++){
	  if( CsiID[in] == klVec[0].pi0()[pt].g2().clusterIdVec()[i]){
	    GamClusCsiSignal[clNumber][i] = CsiTCutSig[in];
	    GamClusCsiCrate[clNumber][i]  = CsiTCutCID[in];
	    break;
	  }
	}
      }
      clNumber++;
    }
    

    trOut->Fill();
  }
  trOut->Write();
  TimeDeltaHist->Write();
  TimeDeltaHistEne->Write();
  TimeDeltaHistEneID[0]->Write();
  TimeDeltaHistEneID[1]->Write();
  for( int i = 0; i< 20; i++){
    TimeDeltaHistEneCrate[i]->Write();
    TimeDeltaHistSigCrate[i]->Write();
  }
  hisIDDelay->Write();
  poly->Write();
  tfOut->Close();
}
  
