/// Convert simulation data with trigger simulation.


#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <string>
#include <list>
#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "cluster/ClusterFinder.h"
//#include "ClusterFinder_EDIT.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "User_Function.h"

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include "IDHandler.h"

//#include "E14WavReader.h"
#include "E14WavReader_V1.h"
//#include "E14WaveReader_V2.h"
#include "L1TrigCounter.h"
#include "EnergyConverter.h"
#include "TRandom.h"
#include "TF1.h"

double funcResolutionInvSq( double* x, double* p){
  double value = 0;
  if( x[0] >  0  && p[0] > 0){
    //value = 10000./(1.26*1.26*1e3/(x[0]*p[0])+16900/(x[0]*x[0]) + 0.76*0.76);
    value = 12.7*x[0]*p[0];
  }
  return value;
}

int
main( int argc ,char ** argv ){
  
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string HOME         = std::getenv("HOME");
  
  std::string ROOTFILE_SIMCONV  = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/ConvFile";
  std::string ROOTFILE_SIMPI0   = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/SIMPI0";
  //std::string ROOTFILE_SIMPI0   = "/Volume0/Simulation/Pi0Run/NewPi0Data_2013";
  std::string iFileForm          = "%s/SimPi0_1E6_LYRES_%d.root";        // ROOTFILE_SIMCONV
  std::string oFileForm          = "%s/SimPi0_1E6_LYRES_Merged.root"; // ROOTFILE_SIM3PI0

  TChain* trin = new TChain("T");
  for( int i = 0; i < 5; i++){
    trin->Add(Form(iFileForm.c_str(),ROOTFILE_SIMPI0.c_str(),i));
  }
  int    RunNo;
  int    EventNumber;
  int    CsiNumber;
  int    CsiModID[2716];
  double CsiEne[2716];
  double CsiTime[2716];
  double CsiHHTime[2716];
  double CsiSignal[2716];
  int    CsiL1nTrig;
  double CsiL1TrigCount[20];
  
  trin->SetBranchAddress("RunNumber",&RunNo);
  trin->SetBranchAddress("EventNumber",&EventNumber);
  trin->SetBranchAddress("CsiNumber",&CsiNumber);
  trin->SetBranchAddress("CsiModID",CsiModID);
  trin->SetBranchAddress("CsiEne",CsiEne);
  trin->SetBranchAddress("CsiTime",CsiTime);
  trin->SetBranchAddress("CsiHHTime",CsiHHTime);
  trin->SetBranchAddress("CsiSignal",CsiSignal);
  trin->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  trin->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);


  int    nTrack;
  UShort_t  track[200];
  int    pid[200];
  float mass[200];
  float ek[200];
  float end_ek[200];
  double p[200][3];
  double end_p[200][3];
  double v[200][3];
  double end_v[200][3];

  trin->SetBranchAddress("nTrack",&nTrack);
  trin->SetBranchAddress("track",track);
  trin->SetBranchAddress("pid",pid);
  trin->SetBranchAddress("mass",mass);
  trin->SetBranchAddress("ek",ek);
  trin->SetBranchAddress("end_ek",end_ek);
  trin->SetBranchAddress("p",p);
  trin->SetBranchAddress("end_p",end_p);
  trin->SetBranchAddress("v",v);
  trin->SetBranchAddress("end_v",end_v);

  //double CsiL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
  //					1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};
  double CsiL1TrigCountThreshold[20] = {1000,3000,3000,3400,3400,3400,2200,2200,2400,2400,
					2400,1000,1000,1000,1000,1000,1000,1000,1000,1000};
  double CsiL1TrigHighThreshold = 50000;
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_SIMPI0.c_str()),"recreate");

  E14GNAnaDataContainer data; 
  data.setBranchAddress(trin);
  const int nHist  =5;
  char* Name[nHist] = {"Neutron","Gamma","KL","KLBG","ETC"};

  TH1D* hisPi0[nHist];
  TH1D* hisPi0Trigged[nHist];
  TH1D* hisGammaE[nHist];
  TH1D* hisGammaChi2[nHist];
  TH1D* hisPi0RecZ[nHist];
  TH1D* hisPi0RecZSig2[nHist];
  TH1D* hisCosTheta[nHist];
  TH1D* hisPi0CutMass[nHist];
  for( int i = 0; i< nHist; i++){

    hisPi0[i] = new TH1D(Form("hisPi0_%d",i),Form("hisPi0_%s",Name[i]),150,0,300 );
    hisPi0Trigged[i] = new TH1D(Form("hisPi0Trigged_%d",i),Form("hisPi0Trigged_%s",Name[i]),150,0,300 );
    hisPi0RecZ[i]    = new TH1D(Form("hisPi0RecZ_%d",i),Form("hisPi0RecZ_%s",Name[i]),60,-300,300);
    hisPi0RecZSig2[i] = new TH1D(Form("hisPi0RecZSig2_%d",i),Form("hisPi0RecZSig2_%s",Name[i]),100,0,10000);
    hisGammaE[i]     = new TH1D(Form("hisGammaE_%d",i),
				Form("hisGammaE_%s;GammaEnergy[MeV]",Name[i]),
				100,0,1000);
    hisGammaChi2[i]  = new TH1D(Form("hisGammaChi2_%d",i),
				Form("hisGammaChi2_%s;GammaChi2[MeV]",Name[i]),
				100,0,100);
    hisCosTheta[i]   = new TH1D(Form("hisCosTheta_%d",i),
				Form("hisCosTheta_%s;CosTheta",Name[i]),
				100,0,1);
    hisPi0CutMass[i] = new TH1D(Form("hisPi0CutMass_%d",i),
				Form("hisPi0CutMass_%s;Pi0RecMass[MeV]",Name[i]),150,0,300);
  }


  TH1D* hisL1TrigCount[11];
  for( int i = 0; i < 11 ; i++){
    if( i == 0){
      hisL1TrigCount[i] = new TH1D(Form("hisL1nTrig"),"hisL1nTrig",15,0,15);
    }else{
      hisL1TrigCount[i] = new TH1D(Form("hisL1TrigCount_%d",i),Form("hisL1TrigCount_%d",i),100,0,50000);    
    }
  }

  TH1D* hisL1TrigCountTrigged[11];  
  for( int i = 0; i < 11 ; i++){
    if( i == 0){
      hisL1TrigCountTrigged[i] = new TH1D(Form("hisL1nTrigTrigged"),"hisL1nTrigTrigged",15,0,15);
    }else{
      hisL1TrigCountTrigged[i] = new TH1D(Form("hisL1TrigCountTrigged_%d",i),Form("hisL1TrigCountTrigged_%d",i),100,0,50000);    
    }
  }




  std::cout<< "LOOP" << std::endl;
  long entries = trin->GetEntries();
  std::cout<< entries << std::endl;
  for( int ievent  = 0; ievent < entries ; ievent++){
    //for( int ievent  = 0; ievent < 100 ; ievent++){
    trin->GetEntry(ievent);
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::list<Pi0>     plist;

    data.getData(plist);

    bool bnEvent  = false;
    bool bkEvent  = false;
    bool bgEvent  = false;
    bool bnreact  = false;
    bool bgreact  = false;
    bool bkreact  = false;
    bool bkdecay  = false;
    bool bksdecay = false;
    bool bETC     = false;

    if( pid[0] == 2112 ){
      bnEvent = true;
    }else if( pid[0] == 130 ){
      bkEvent = true;
    }else if( pid[0] == 22 ){
      bgEvent = true;
    }
    int nPi0 = 0;
    int nPipm= 0;
    int nEta = 0;
    int other= 0;
    for( int ip = 1 ; ip < nTrack; ip++){
      if( pid[1]  == 310 ){ bksdecay = true;}
      if( pid[ip] == 211 ){ bksdecay = true;}
      if( pid[ip] == 111 ){ nPi0++; }
      if( pid[ip] == 211 || pid[ip] == -211 ){ nPipm++;}
    }
    
    if( bnEvent ){
      //if( nPi0 ==1 ){ bnreact = true; }
      if( abs(end_v[0][2]-3536) < 10 ){
	bnreact = true;
      }else{
	bETC = true;
      }
    }else if( bgEvent ){
      if( nPi0 == 1 ){ bgreact = true; }
      else{bETC = true;}
    }else if( bkEvent ){
      if( abs(end_v[0][2]-3536) < 10 ){
	bkreact = true;
      }else{
	bkdecay = true; 
      }
    }else{
      bETC = true;
    }
    std::list<Pi0>::iterator pit = plist.begin();
    Int_t hisID = -1;
    if( bnreact ){
      hisID = 0;
    }else if( bgreact ){
      hisID = 1;
    }else if( bkreact ){
      hisID = 2;
    }else if( bkdecay ){       
      hisID = 3;
    }else if( bETC ){
      hisID = 4;
    }
      
    if( hisID >= 0 ){
      hisPi0[hisID]->Fill((*pit).m());    
    }

    int nTrig = 0; 
    for( int i = 1; i< 11; i++){
      if( CsiL1TrigCount[i] > CsiL1TrigCountThreshold[i]  &&
	  CsiL1TrigCount[i] < CsiL1TrigHighThreshold ){
	nTrig++;
      }
    }

    hisL1TrigCount[0]->Fill(nTrig);
    for( int i = 1; i< 11; i++){
      hisL1TrigCount[i]->Fill(CsiL1TrigCount[i]);
    }

    if( nTrig >=2 ){
      if( hisID >= 0){      
	hisPi0Trigged[hisID]->Fill((*pit).m());          
	hisPi0RecZ[hisID]->Fill((*pit).recZ()-(*pit).vz());
	hisPi0RecZSig2[hisID]->Fill((*pit).recZsig2());
	hisGammaE[hisID]->Fill((*pit).g1().e());
	hisGammaE[hisID]->Fill((*pit).g2().e());
	hisGammaChi2[hisID]->Fill((*pit).g1().chisq());
	hisGammaChi2[hisID]->Fill((*pit).g2().chisq());

	double x[2]; 
	x[0] = (*pit).g1().x();
	x[1] = (*pit).g2().x();
	double y[2]; 
	y[0] = (*pit).g1().y();
	y[1] = (*pit).g2().y();
	double cosTheta = TMath::Abs( x[0]*x[1]+y[0]*y[1] )/TMath::Sqrt((x[0]*x[0]+y[0]*y[1])*(x[1]*x[1]+y[1]*y[1]));
	hisCosTheta[hisID]->Fill(cosTheta);
	double Eg1 = (*pit).g1().e();
	double Eg2 = (*pit).g2().e();
	double gchisq_1 = (*pit).g1().chisq();
	double gchisq_2 = (*pit).g2().chisq();
	double pi0pt    = TMath::Sqrt((*pit).p3()[0]*(*pit).p3()[0]+ (*pit).p3()[1]*(*pit).p3()[1]);
	if( Eg1 > 350 &&
	    Eg2 > 150 &&
	    gchisq_1 < 5 && 
	    gchisq_2 < 5 &&
	    pi0pt  > 50  &&
	    cosTheta < 0.9 ){
	  hisPi0CutMass[hisID]->Fill((*pit).m());
	}


      }      
      hisL1TrigCountTrigged[0]->Fill(nTrig);
      for( int i = 1; i< 11; i++){
	hisL1TrigCountTrigged[i]->Fill(CsiL1TrigCount[i]);
      }
    }    
  }

  std::cout<< "Write" << std::endl;
  for( int i = 0; i < nHist; i++){
    hisPi0[i]->Write();
  }
  for( int i = 0; i < nHist; i++){
    hisPi0CutMass[i]->Write();
  }

  for( int i = 0; i< nHist; i++){
    hisPi0Trigged[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisPi0RecZ[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisPi0RecZSig2[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisGammaE[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisGammaChi2[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisCosTheta[i]->Write();
  }
  for( int i = 0; i < 11; i++){
    hisL1TrigCount[i]->Write();
  }
  for( int i = 0; i < 11; i++){
    hisL1TrigCountTrigged[i]->Write();
  }
    
  tfout->Close();
  return 0;
}
