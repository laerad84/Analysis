#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TH3.h"

#include "TChain.h"
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
#include "T0Manager.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "IDHandler.h"
#include "User_Function.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "csimap/CsiMap.h"


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
  const int cnTrig        = 0;
  const int cgammaPosIn   = 1;
  const int cgammaPosOut  = 2;
  const int cgammaEne     = 3;
  const int cPi0Pt        = 4;
  const int cKLChisqZ     = 5;
  const int cklmass       = 6;
  const int cklVtx        = 7;

  TF1* soltFunc  = new TF1("soltfunc",KLSpectrum,0,12,5);
  TF1* sugarFunc = new TF1("sugarfunc",KLSpectrum,0,12,5);
  const double soltPar[5] = { 1, 1.37,7.48e-1,-2.9e-1,1.68e-2};
  const double sugarPar[5] = {1,1.41991,0.810237,-0.301413,0.017092};
  soltFunc->SetParameters(soltPar);
  sugarFunc->SetParameters(sugarPar);
  IDHandler* handler = new IDHandler();

  CsiMap* map = CsiMap::getCsiMap();

  int mode = std::atoi( argv[1]);
  const int nFile = 1;
  //TFile* tf;
  //TTree* tr;
  TChain* tr;
  std::string name;
  switch( mode ){
  case 0:
    name = "DATA_NONTIMECAL1";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"
    tr = new TChain("trKL");
    tr->Add(Form("kl_KL_%s.root",name.c_str()));
    break;
  case 1:
    name = "DATA_NONTIMECALNOCV";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"
    tr = new TChain("trKL");
    tr->Add(Form("kl_KL_%s.root",name.c_str()));
    break;
  case 2:
    name = "SIM_1e9";
    tr = new TChain("T");
    for( int i = 0; i< 20; i++){
      tr->Add(Form("Sim_e14_KL3pi0_KL_RES_LY_pe_1E8_NON10_%d.root",i));
    }
    break;
  case 3:
    name = "SIM";
    tr = new TChain("trKL");
    tr->Add("Kl_Total_SIM.root");//"Sim_e14_KL3pi0_KL_RES_LY_pe_1E8_NON10_%d.root",i));
    break;

  case 4:
    name = "SIM_SATO";
    tr = new TChain("trKL");
    tr->Add("Kl_Total_SIM_SATO.root");//"Sim_e14_KL3pi0_KL_RES_LY_pe_1E8_NON10_%d.root",i));
    break;

  default:
    return -1; 
  }
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  Int_t CsiNumber;
  Int_t CsiModID[3000];
  Double_t CsiSignal[3000];
  Double_t CsiTime[3000];
  Double_t CsiEne[3000];
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  /*
  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiModID",CsiModID);//CsiNumber
  tr->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  */
  //tr->SetCacheSize(-1); 
  std::cout<< __LINE__ << std::endl;


  std::cout<< __LINE__ << std::endl;
  std::cout<< __LINE__ << std::endl;
  TFile* tfOut = new TFile(Form("kl_%s.root",name.c_str()),"recreate");
  TTree* trOut = new TTree("KLDistribution","Time correction");
  
  int CutConditionKL;
  int EventID;
  E14GNAnaDataContainer dataCp;
  std::cout<< __LINE__ << std::endl;

  std::cout<< __LINE__ << std::endl;
  dataCp.branchOfKlong(trOut);

  trOut->Branch("CsiL1TrigCount",CsiL1TrigCount,"CsiL1TrigCount[20]/D");
  trOut->Branch("CsiL1nTrig",&CsiL1nTrig,"CsiL1nTrig/I");
  trOut->Branch("CutConditionKL",&CutConditionKL,"CutConditionKL/I");
  trOut->Branch("EventID",&EventID,"EventID/I");
  std::cout<< __LINE__ << std::endl;
  
  const int nHist =3; 
  TH1D* hisL1nTrig[nHist];
  TH1D* hisL1TrigCount[nHist][20];
  TH1D* hisKLMass[nHist];
  TH1D* hisGammaE[nHist];
  TH1D* hisGammaX[nHist];
  TH1D* hisGammaY[nHist];
  
  TH1D* hisKLX[nHist];
  TH1D* hisKLY[nHist];
  TH1D* hisKLZ[nHist];
  TH1D* hisKLP[nHist];
  TH1D* hisKLChisqZ[nHist];
  std::cout<< __LINE__ << std::endl;

  for( int i = 0; i< nHist; i++){
    hisL1nTrig[i] = new TH1D(Form("hisL1nTrig_%d",i),Form("hisL1nTrig_%d",i),20,0,20);
    for( int j = 0; j < 20; j++){
      hisL1TrigCount[i][j] = new TH1D(Form("hisL1TrigCount_%d_%d",i,j),Form("hisL1TrigCount_%d_%d",i,j),
				      100,0,10000);
    }
    hisKLMass[i] = new TH1D(Form("hisKlMass_%d",i),Form("hisKlMass_%d",i),400,400,600);
    hisGammaE[i] = new TH1D(Form("hisGammaE_%d",i),Form("hisGammaE_%d",i),400,0,4000);
    hisGammaX[i] = new TH1D(Form("hisGammaX_%d",i),Form("hisGammaX_%d",i),80,-1000,1000);
    hisGammaY[i] = new TH1D(Form("hisGammaY_%d",i),Form("hisGammaY_%d",i),80,-1000,1000);
    hisKLX[i] = new TH1D(Form("hisKLX_%d",i),Form("hisKLX_%d",i),120,-300,300);
    hisKLY[i] = new TH1D(Form("hisKLY_%d",i),Form("hisKLY_%d",i),120,-300,300);
    hisKLZ[i] = new TH1D(Form("hisKLZ_%d",i),Form("hisKLZ_%d",i),60,0,6000);
    hisKLP[i] = new TH1D(Form("hisKLP_%d",i),Form("hisKLP_%d",i),60,0,6000);
    hisKLMass[i] = new TH1D(Form("hisKLMass_%d",i),Form("hisKLMass_%d",i),60,0,6000);
    hisKLChisqZ[i] = new TH1D(Form("hisKLChisqZ_%d",i),Form("hisKLChisqZ_%d",i),60,0,6000);
  }
  std::cout<< "Start" << std::endl;
  for( int ievent  =0; ievent < tr->GetEntries(); ievent++){
    if( (ievent %  1000 ) == 0 ){ 
      std::cout<< (double)(ievent)/tr->GetEntries() << std::endl;
    }
    tr->GetEntry( ievent );
    dataCp.reset();
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;

    data.getData(glist);
    data.getData(klVec);
    if( klVec.size() < 1 ){ continue; }
    CutConditionKL = 0;
    EventID = ievent;
    
    if( CsiL1nTrig < 5 ){ 
      CutConditionKL |= 1 << cnTrig;
    }
    hisL1nTrig[0]->Fill( CsiL1nTrig );
    for( int i = 0; i < 20; i++){
      hisL1TrigCount[0][i]->Fill( CsiL1TrigCount[i] );
    }
    if( (CutConditionKL & 1 ) != 1 ){
      hisL1nTrig[1]->Fill( CsiL1nTrig );
      for( int i = 0; i < 20; i++){
	hisL1TrigCount[1][i]->Fill( CsiL1TrigCount[i] );
      }
    }

    std::list<Gamma>::iterator git = glist.begin();
    for( int i  =0; i< 6; i++, git++){
      double x = (*git).x();
      double y = (*git).y();
      if( (*git).e() < 200 ){ 
	CutConditionKL |= 1 << cgammaEne;
      }
      if( TMath::Abs(x) < 150 && TMath::Abs(y) < 150 ){
	CutConditionKL |= 1 << cgammaPosIn;
      }
      if( TMath::Abs(y) > 550  || sqrt(x*x+y*y) > 850 ){
	CutConditionKL |= 1 << cgammaPosOut;
      }
    }

    /*
    for( int i = 0; i< klVec[0].pi0().size(); i++){
      if( klVec[0].pi0()[i].p3().perp() > 150 ){
	CutConditionKL |= 1 << cPi0Pt;
      }
    }
    */

    if( klVec[0].chisqZ() > 15 ){ 
      CutConditionKL |= 1 << cKLChisqZ;
    }
    if( TMath::Abs(klVec[0].m() -KLMass) > 10 ){
      CutConditionKL |= 1 << cklmass;
    }
    
    if( klVec[0].vz() < 2000 || klVec[0].vz() > 5000 ){
      CutConditionKL |= 1 << cklVtx;
    }

    hisKLChisqZ[0]->Fill(klVec[0].chisqZ());
    if( (CutConditionKL & 1 ) != 1 ){
      hisKLMass[0]->Fill(klVec[0].m());
      hisKLX[0]->Fill(klVec[0].vx() );
      hisKLY[0]->Fill(klVec[0].vy() );
      hisKLZ[0]->Fill(klVec[0].vz() );      
      hisKLP[0]->Fill(klVec[0].p3().mag());
      git = glist.begin();
      for( int i = 0; i< 6; i++){
	hisGammaE[0]->Fill((*git).e());
	hisGammaX[0]->Fill((*git).x());
	hisGammaY[0]->Fill((*git).y());
      }

      //hisKLChisqZ[1]->Fill(klVec[1].chisqZ());
      if( (CutConditionKL & 15)== 0){
	hisKLChisqZ[1]->Fill(klVec[1].chisqZ());
	hisKLMass[1]->Fill(klVec[0].m());
	hisKLX[1]->Fill(klVec[0].vx() );
	hisKLY[1]->Fill(klVec[0].vy() );
	hisKLZ[1]->Fill(klVec[0].vz() );   
	hisKLP[1]->Fill(klVec[0].p3().mag());
	git = glist.begin();
	for( int i = 0; i< 6; i++,git++){
	  hisGammaE[1]->Fill((*git).e());
	  hisGammaX[1]->Fill((*git).x());
	  hisGammaY[1]->Fill((*git).y());
	}
      }

      if( CutConditionKL == 0 ){
	hisKLChisqZ[2]->Fill(klVec[1].chisqZ());
	hisKLMass[2]->Fill(klVec[0].m());
	hisKLX[2]->Fill(klVec[0].vx() );
	hisKLY[2]->Fill(klVec[0].vy() );
	hisKLZ[2]->Fill(klVec[0].vz() );   
	hisKLP[2]->Fill(klVec[0].p3().mag());
	git = glist.begin();
	for( int i = 0; i< 6; i++){
	  hisGammaE[2]->Fill((*git).e());
	  hisGammaX[2]->Fill((*git).x());
	  hisGammaY[2]->Fill((*git).y());
	}
      }
    }
    dataCp.setData(klVec);
    trOut->Fill();
  }
  for( int i = 0; i<2; i++){
    hisL1nTrig[i]->Write();
    for( int j = 0; j<20; j++){
    hisL1TrigCount[i][j]->Write();    
    }
    hisGammaE[i]->Write();
    hisGammaX[i]->Write();
    hisGammaY[i]->Write();
    hisKLX[i]->Write();
    hisKLY[i]->Write();
    hisKLZ[i]->Write();
    hisKLP[i]->Write();
    hisKLMass[i]->Write();
    hisKLChisqZ[i]->Write();
  }
  trOut->Write();
  tfOut->Close();
}
