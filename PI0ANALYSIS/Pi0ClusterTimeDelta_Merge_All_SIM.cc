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
  std::string iFileForm          = "%s/SimPi0_1E6_LYRES_KLBEAM_%d.root";        // ROOTFILE_SIMCONV
  std::string oFileForm          = "%s/SimPi0_1E6_LYRES_Merged_NEW_KLBEAM.root"; // ROOTFILE_SIM3PI0

  TChain* trin = new TChain("T");
  for( int i = 0; i < 2; i++){
    trin->Add(Form(iFileForm.c_str(),ROOTFILE_SIMPI0.c_str(),i));
  }
  E14GNAnaDataContainer data; 
  data.setBranchAddress(trin);

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
  
  int CVNumber;
  int CVModID[256];
  Double_t CVEne[256];
  int SciNumber;
  Double_t SciEne[1];
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
  trin->SetBranchAddress("CVNumber",&CVNumber);
  trin->SetBranchAddress("CVModID",CVModID);
  trin->SetBranchAddress("CVEne",CVEne);
  trin->SetBranchAddress("SciEne",SciEne);

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

  const int nHist  =5;
  char* Name[nHist] = {"Neutron","Gamma","KL","KLBG","ETC"};

  TH1D* hisPi0[nHist];
  TH1D* hisPi0Trigged[nHist];
  TH1D* hisGammaE[nHist];
  TH1D* hisGammaECutHigh[nHist];
  TH1D* hisGammaECutLow[nHist];
  TH1D* hisGammaChi2[nHist];
  TH1D* hisPi0RecZ[nHist];
  TH1D* hisPi0RecZSig2[nHist];
  TH1D* hisCosTheta[nHist];
  TH1D* hisPi0CutMass[nHist];
  TH1D* hisPi0E[nHist];
  TH1D* hisPi0ECut[nHist];
  TH2D* hisPi0MassGammaEH[nHist];
  TH2D* hisPi0MassGammaEL[nHist];
  TH2D* hisPi0MassCenterE[nHist];
  TH2D* hisPi0MassHeight[nHist];
  TH2D* hisGammaECompare[nHist];

  for( int i = 0; i< nHist; i++){
    

    hisPi0ECut[i] = new TH1D(Form("hisPi0ECut_%d",i),Form("hisPi0ECut_%s",Name[i]),200,0,5000);
    hisPi0E[i] = new TH1D(Form("hisPi0E_%d",i),Form("hisPi0E_%s",Name[i]),200,0,5000);
    hisPi0[i] = new TH1D(Form("hisPi0_%d",i),Form("hisPi0_%s",Name[i]),150,0,300 );
    hisPi0Trigged[i] = new TH1D(Form("hisPi0Trigged_%d",i),Form("hisPi0Trigged_%s",Name[i]),150,0,300 );
    hisPi0RecZ[i]    = new TH1D(Form("hisPi0RecZ_%d",i),Form("hisPi0RecZ_%s",Name[i]),60,-300,300);
    hisPi0RecZSig2[i] = new TH1D(Form("hisPi0RecZSig2_%d",i),Form("hisPi0RecZSig2_%s",Name[i]),100,0,10000);
    hisGammaE[i]     = new TH1D(Form("hisGammaE_%d",i),
				Form("hisGammaE_%s;GammaEnergy[MeV]",Name[i]),
				150,0,3000);
    hisGammaECutHigh[i]     = new TH1D(Form("hisGammaECutHigh_%d",i),
				       Form("hisGammaECutHigh_%s;GammaEnergy[MeV]",Name[i]),
				       150,0,3000);
    hisGammaECutLow[i]     = new TH1D(Form("hisGammaECutLow_%d",i),
				      Form("hisGammaECutLow_%s;GammaEnergy[MeV]",Name[i]),
				      150,0,3000);
    hisGammaChi2[i]  = new TH1D(Form("hisGammaChi2_%d",i),
				Form("hisGammaChi2_%s;GammaChi2[MeV]",Name[i]),
				100,0,100);
    hisCosTheta[i]   = new TH1D(Form("hisCosTheta_%d",i),
				Form("hisCosTheta_%s;CosTheta",Name[i]),
				100,0,1);
    hisPi0CutMass[i] = new TH1D(Form("hisPi0CutMass_%d",i),
				Form("hisPi0CutMass_%s;Pi0RecMass[MeV]",Name[i]),150,0,300);
    hisPi0MassGammaEH[i] = new TH2D(Form("hisPi0MassGammaEH_%d",i),
				    Form("hisPi0MassGammaEH_%s;GammaE[MeV]",Name[i]),100,0,4000,150,0,300);
    hisPi0MassGammaEL[i] = new TH2D(Form("hisPi0MassGammaEL_%d",i),
				    Form("hisPi0MassGammaEL_%s;GammaE[MeV]",Name[i]),100,0,4000,150,0,300);
    hisPi0MassHeight[i] = new TH2D(Form("hisPi0MassHeight_%d",i),
				Form("hisPi0MassHeight_%s;Height[cnt]",Name[i]),160,0,16000,150,0,300);
    hisPi0MassCenterE[i] = new TH2D(Form("hisPi0MassCenterE_%d",i),
				    Form("hisPi0MassCenterE_%s",Name[i]),100,0,2000,150,0,300);    

    hisGammaECompare[i] = new TH2D(Form("hisGammaECompare_%d",i),
				   Form("hisGammaECompare_%s",Name[i]),100,0,2000,200,0.5,1.5);
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
    
    if(SciEne[0] < 3.86 ){ continue; }
    double CVMax = 0;
    for( int icv  =0; icv < CVNumber; icv++){
      if(CVEne[icv] > CVMax ){
	CVMax = CVEne[icv];
      }
    }
    if( CVMax > 0.4 ){ continue; }

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
      if( TMath::Abs(end_v[0][2]-3536) < 10 ){
	bnreact = true;
      }else{
	bETC = true;
      }
    }else if( bgEvent ){
      if( nPi0 == 1 ){ bgreact = true; }
      else{bETC = true;}
    }else if( bkEvent ){
      if( TMath::Abs(end_v[0][2]-3536) < 10 ){
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

	double R[2];
	
	bool bPosition = true;
	if( TMath::Abs(y[0])> 550  ||
	    TMath::Abs(y[1])> 550  ){ 
	  bPosition = false;
	}
	for( int ig = 0; ig < 2; ig++){
	  R[ig] = TMath::Sqrt( x[ig]*x[ig] + y[ig]*y[ig]);
	  if( R[ig] > 850){ bPosition = false; }
	  if( TMath::Abs(y[ig]) < 150 && TMath::Abs(x[ig])< 150 ){
	    bPosition = false;
	  }	  
	}
	if( !bPosition ){ continue; }


	Double_t GammaEnergy[100];
	Double_t GammaX[100];
	Double_t GammaY[100];
	Double_t GammaZ[100];
	Double_t GammaProductZ[100];
	int nTotalGamma = 0;
	bool bAlProductedGamma = false;
	for( int ip = 0; ip < nTrack; ip++){
	  if( pid[ip] == 22 ){
	    GammaEnergy[nTotalGamma] = ek[ip];
	    GammaX[nTotalGamma]      = end_v[ip][0];
	    GammaY[nTotalGamma]      = end_v[ip][1];
	    GammaZ[nTotalGamma]      = end_v[ip][2];
	    GammaProductZ[nTotalGamma]= v[ip][2];	    
	    nTotalGamma++;
	  }
	}

	int GammaIndex[2] = {-1,-1};
	Double_t Dist[2][2] = {{0}};
	if( nTotalGamma == 2){
	  if( GammaZ[0] > 6148 && GammaZ[1] < 6648 &&
	      GammaZ[1] > 6148 && GammaZ[2] < 6648 &&
	      TMath::Abs(GammaY[0]) < 600 && 
	      TMath::Abs(GammaY[1]) < 600 &&
	      GammaX[0]*GammaX[0]+GammaY[0]*GammaY[0] < 850*850 &&
	      GammaX[1]*GammaX[1]+GammaY[1]*GammaY[1] < 850*850 ){
	      
		       
	    if( TMath::Abs(GammaProductZ[0]-3536)<10  &&
		TMath::Abs(GammaProductZ[1]-3536)<10 ){
	      for( int isg = 0; isg< 2; isg++){
		for( int ig = 0; ig < 2; ig++){
		  Dist[isg][ig] = sqrt((GammaX[isg]-x[ig])*(GammaX[isg]-x[ig])+(GammaY[isg]-y[ig])*(GammaY[isg]-y[ig]));
	      }
	      }	  
	      if( Dist[0][0]*Dist[1][1] < Dist[0][1]*Dist[1][0] ){
	      hisGammaECompare[hisID]->Fill(GammaEnergy[0],(*pit).g1().e()/GammaEnergy[0]);
	      hisGammaECompare[hisID]->Fill(GammaEnergy[1],(*pit).g2().e()/GammaEnergy[1]);
	      }else{
		hisGammaECompare[hisID]->Fill(GammaEnergy[0],(*pit).g1().e()/GammaEnergy[1]);
		hisGammaECompare[hisID]->Fill(GammaEnergy[1],(*pit).g2().e()/GammaEnergy[0]);
	      }
	    }
	  }
	}
      
	


	int    ClusterID[2] ={0};
	double ClusterHeight[2] ={0};
	double MaximumHeight=0;
	ClusterID[0] = (*pit).g1().clusterIdVec()[0];
	ClusterID[1] = (*pit).g2().clusterIdVec()[0];
	Int_t nMatched  = 0;
	for( int iCsi  =0; iCsi < CsiNumber; iCsi++){	  
	  //std::cout << CsiModID[iCsi] << "\t" << CsiSignal[iCsi] << "\t" << ClusterID[0] << "\t" << ClusterID[1] << std::endl;
	  if( CsiModID[iCsi] == ClusterID[0] ){
	    ClusterHeight[0] = CsiSignal[iCsi];
	    nMatched++;
	  }
	  if(CsiModID[iCsi] == ClusterID[1] ){
	    ClusterHeight[1] = CsiSignal[iCsi];
	    nMatched++;
	  }
	  if( nMatched == 2 ){
	    break;
	  }
	}
	if( ClusterHeight[0] > ClusterHeight[1] ){
	  MaximumHeight = ClusterHeight[0];
	}else{
	  MaximumHeight = ClusterHeight[1];
	}
	double cosTheta = TMath::Abs( x[0]*x[1]+y[0]*y[1] )/TMath::Sqrt((x[0]*x[0]+y[0]*y[1])*(x[1]*x[1]+y[1]*y[1]));
	hisCosTheta[hisID]->Fill(cosTheta);
	double Eg1 = (*pit).g1().e();
	double Eg2 = (*pit).g2().e();
	double gchisq_1 = (*pit).g1().chisq();
	double gchisq_2 = (*pit).g2().chisq();
	double pi0pt    = TMath::Sqrt((*pit).p3()[0]*(*pit).p3()[0]+ (*pit).p3()[1]*(*pit).p3()[1]);
	double pi0Mass  = (*pit).m();
	if( Eg1 > 350 &&
	    Eg2 > 200 &&
	    gchisq_1 < 5 && 
	    gchisq_2 < 5 &&
	    pi0pt  > 50  &&
	    cosTheta < 0.9 ){
	  hisPi0CutMass[hisID]->Fill((*pit).m());
	  hisPi0E[hisID]->Fill((*pit).e());
	  hisGammaECutHigh[hisID]->Fill((*pit).g1().e());
	  hisGammaECutLow[hisID]->Fill((*pit).g2().e());

	  hisPi0MassGammaEH[hisID]->Fill( Eg1, pi0Mass );
	  hisPi0MassGammaEL[hisID]->Fill( Eg2, pi0Mass );
	  hisPi0MassHeight[hisID]->Fill( MaximumHeight, pi0Mass);
	  hisPi0MassCenterE[hisID]->Fill( (*pit).g1().clusterEVec()[0], pi0Mass);
	  hisPi0MassCenterE[hisID]->Fill( (*pit).g2().clusterEVec()[0], pi0Mass);

	  if( TMath::Abs((*pit).m()-135) < 10 ){
	    hisPi0ECut[hisID]->Fill((*pit).e());
	  }
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
  for( int i = 0; i < nHist; i++){
    hisPi0E[i]->Write();
  }
  for( int i = 0; i < nHist; i++){
    hisPi0ECut[i]->Write();
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
    hisGammaECutHigh[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisGammaECutLow[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisGammaChi2[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisCosTheta[i]->Write();
  }
  for( int i = 0; i< nHist; i++){
    hisGammaECompare[i]->Write();
  }
  for( int i = 0; i < 11; i++){
    hisL1TrigCount[i]->Write();
  }
  for( int i = 0; i < 11; i++){
    hisL1TrigCountTrigged[i]->Write();
  }
    
  for( int i = 0; i< nHist; i++){
    hisPi0MassGammaEH[i]->Write();
    hisPi0MassGammaEL[i]->Write();
    hisPi0MassHeight[i]->Write();
    hisPi0MassCenterE[i]->Write();
  }
  tfout->Close();
  return 0;
}
