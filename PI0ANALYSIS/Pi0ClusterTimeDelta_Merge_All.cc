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


const double KLMass = 497.648;//MeV
double const sol = 299.792458;//[mm/nsec]
double const solc= 80;//[mm/nsec]
double const Pcor[2]={6.49003,0.99254};
double const CsIX0=18.5;//mm
double showerDepth(double x){
  double L = CsIX0*(Pcor[0]+Pcor[1]*log(x/1000.));//mm  
  return L;
}


double showerTimeDelay(Pi0 kl, Gamma g){
  double depth     = showerDepth(g.e());
  double cosTheta  = abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  double delayTime = 500/solc - depth*(1-cosTheta)/sol;
  return delayTime;
}


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
  std::string ROOTFILE_SIMPI0   = "/gpfs/fs03/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/Pi0Run/SIMPI0";//kekcc
  //std::string ROOTFILE_SIMPI0   = "/Volume0/Simulation/Pi0Run/NewPi0Data_2013";//local
  //std::string iFileForm          = "%s/SimPi0_1E6_LYRES_%d.root";     // ROOTFILE_SIMCONV
  //std::string oFileForm          = "%s/SimPi0_1E6_LYRES_Merged.root"; // ROOTFILE_SIM3PI0

  //std::string iFileForm          = "%s/run_wav_%d_Cal_FNL_COS_newCompensate_pi0.root";
  //std::string oFileForm          = "%s/Pi0_wav_Merged_Data_LaserComp.root";

  std::string iFileForm          = "%s/run_wav_%d_Cal_FNL_COS_newTimeOffset_pi0.root";
  std::string oFileForm          = "%s/Pi0_wav_Merged_Data_IwaiSato_NEW.root";

  Int_t RunN[24]={4502,4503,4504,4505,4506,4507,4508,4509,4510,4511,4512,4513,
		  4514,4515,4516,4517,4518,4519,4520,4521,4522,4523,4524,4525};
  TChain* trin = new TChain("T");
  for( int i = 0; i < 24; i++){
    trin->Add(Form(iFileForm.c_str(),ROOTFILE_WAV.c_str(),RunN[i]));
  }
  E14GNAnaDataContainer data; 

  data.setBranchAddress(trin);

  int    RunNo;
  int    EventNumber;
  int    CsiNumber;

  double CsiEne[2716];
  int    CsiModID[2716];
  double CsiTime[2716];
  double CsiHHTime[2716];
  double CsiSignal[2716];
  int    CsiL1nTrig;
  double CsiL1TrigCount[20];

  int CVNumber;
  Short_t CVModID[256];
  Double_t CVEne[256];
  int SciNumber;
  Double_t SciEne;

  trin->SetBranchAddress("RunNumber",&RunNo);
  trin->SetBranchAddress("EventNumber",&EventNumber);
  trin->SetBranchAddress("CsiNumber",&CsiNumber);

  trin->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  trin->SetBranchAddress("CsiModID",CsiModID);//CsiNumber
  trin->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  trin->SetBranchAddress("CsiHHTime",CsiHHTime);//CsiNumber
  trin->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  trin->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  trin->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  trin->SetBranchAddress("CVNumber",&CVNumber);
  trin->SetBranchAddress("CVModID",CVModID);//CVNumber
  trin->SetBranchAddress("CVEne",CVEne);//CVNumber
  trin->SetBranchAddress("SciNumber",&SciNumber);
  trin->SetBranchAddress("SciEne",&SciEne);//SciNumber
  /*
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
  */

  //double CsiL1TrigCountThreshold[20] = {1000,1800,1800,1800,1800,1800,1200,1200,1200,1200,
  //					1300,1000,1000,1000,1000,1000,1000,1000,1000,1000};
  Double_t TimeOffset[2716]={0};
  std::ifstream ifs("~/local/Analysis/KLongSpectrum/Data/TimeOffset_Shower_10.dat");
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



  double CsiL1TrigCountThreshold[20] = {1000,3000,3000,3400,3400,3400,2200,2200,2400,2400,
					2400,1000,1000,1000,1000,1000,1000,1000,1000,1000};
  //double CsiL1TrigHighThreshold = 50000;
  double CVThreshold[10] = {700,600,600,600,500,600,600,550,600,600};
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  TFile* tfout= new TFile("Pi0_All_Out.root","recreate");
  TTree* trOut= new TTree("GammaTimeTree","GammaTime");
  E14GNAnaDataContainer dataCopy;
  dataCopy.branchOfPi0List(trOut);
  double GammaCE[2];
  double GammaTime[2];
  double GammaTimeSigma;
  double GammaTimeSigmaExcept[2];
  double GammaE[2];
  double GammaX[2];
  double GammaY[2];
  trOut->Branch("GammaCE",GammaCE,"GammaCE[2]/D");
  trOut->Branch("GammaTime",GammaTime,"GammaTime[2]/D");
  trOut->Branch("GammaTimeSigma",&GammaTimeSigma,"GammaTimeSigma/D");
  trOut->Branch("GammaTimeSigmaExcept",GammaTimeSigmaExcept,"GammaTimeSigmaExcept[2]/D");
  trOut->Branch("GammaE",GammaE,"GammaE[2]/D");
  trOut->Branch("GammaX",GammaX,"GammaX[2]/D");
  trOut->Branch("GammaY",GammaY,"GammaY[2]/D");


  //TFile* tfout = new TFile(Form(oFileForm.c_str(),ROOTFILE_SIMPI0.c_str()),"recreate");

  const int nHist  =1;
  char* Name[nHist] = {"Data"};

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



  for( int i = 0; i< nHist; i++){

    hisPi0ECut[i]    = new TH1D(Form("hisPi0ECutData_%d",i),Form("hisPi0CutE_%s",Name[i]),200,0,5000 );
    hisPi0E[i]       = new TH1D(Form("hisPi0EData_%d",i),Form("hisPi0E_%s",Name[i]),200,0,5000 );
    hisPi0[i]        = new TH1D(Form("hisPi0Data_%d",i),Form("hisPi0_%s",Name[i]),150,0,300 );
    hisPi0Trigged[i] = new TH1D(Form("hisPi0TriggedData_%d",i),Form("hisPi0Trigged_%s",Name[i]),150,0,300 );
    hisPi0RecZ[i]    = new TH1D(Form("hisPi0RecZData_%d",i),Form("hisPi0RecZ_%s",Name[i]),60,-300,300);
    hisPi0RecZSig2[i]= new TH1D(Form("hisPi0RecZSig2Data_%d",i),Form("hisPi0RecZSig2_%s",Name[i]),100,0,10000);
    hisGammaE[i]     = new TH1D(Form("hisGammaEData_%d",i),
				Form("hisGammaE_%s;GammaEnergy[MeV]",Name[i]),
				150,0,3000);
    hisGammaECutHigh[i] = new TH1D(Form("hisGammaECutHighData_%d",i),
				   Form("hisGammaECutHigh_%s;GammaEnergy[MeV]",Name[i]),
				   150,0,3000);
    hisGammaECutLow[i]  = new TH1D(Form("hisGammaECutLowData_%d",i),
				   Form("hisGammaECutLow_%s;GammaEnergy[MeV]",Name[i]),
				   150,0,3000);
    hisGammaChi2[i]  = new TH1D(Form("hisGammaChi2Data_%d",i),
				Form("hisGammaChi2_%s;GammaChi2[MeV]",Name[i]),
				100,0,100);
    hisCosTheta[i]   = new TH1D(Form("hisCosThetaData_%d",i),
				Form("hisCosTheta_%s;CosTheta",Name[i]),
				100,0,1);
    hisPi0CutMass[i] = new TH1D(Form("hisPi0CutMassData_%d",i),
				Form("hisPi0CutMass_%s;Pi0RecMass[MeV]",Name[i]),150,0,300);
    hisPi0MassGammaEH[i] = new TH2D(Form("hisPi0MassGammaEHData_%d",i),
				    Form("hisPi0MassGammaEH_%s;GammaE[MeV]",Name[i]),100,0,4000,150,0,300);
    hisPi0MassGammaEL[i] = new TH2D(Form("hisPi0MassGammaELData_%d",i),
				    Form("hisPi0MassGammaEL_%s;GammaE[MeV]",Name[i]),100,0,4000,150,0,300);
    hisPi0MassHeight[i] = new TH2D(Form("hisPi0MassHeightData_%d",i),
				Form("hisPi0MassHeight_%s;Height[cnt]",Name[i]),160,0,16000,150,0,300);
    hisPi0MassCenterE[i] = new TH2D(Form("hisPi0MassCenterEData_%d",i),
				    Form("hisPi0MassCenterE_%s",Name[i]),100,0,2000,150,0,300);    
  }

  TH1D* hisL1TrigCount[11];
  for( int i = 0; i < 11 ; i++){
    if( i == 0){
      hisL1TrigCount[i] = new TH1D(Form("hisL1nTrigData"),"hisL1nTrigData",15,0,15);
    }else{
      hisL1TrigCount[i] = new TH1D(Form("hisL1TrigCountData_%d",i),Form("hisL1TrigCountData_%d",i),100,0,50000);    
    }
  }

  TH1D* hisL1TrigCountTrigged[11];  
  for( int i = 0; i < 11 ; i++){
    if( i == 0){
      hisL1TrigCountTrigged[i] = new TH1D(Form("hisL1nTrigTrigged"),"hisL1nTrigTrigged",15,0,15);
    }else{
      hisL1TrigCountTrigged[i] = new TH1D(Form("hisL1TrigCountTriggedData_%d",i),Form("hisL1TrigCountTriggedData_%d",i),100,0,50000);    
    }
  }

  TH1D* hisCVEneDistrib[20];
  TH1D* hisCVEneDistribTrigged[20];
  for( int i = 0; i< 20; i++){
    hisCVEneDistrib[i] = new TH1D(Form("hisCVEneDistrib_%d",i),Form("hisCVEneDistrib_%d",i),400,0,4000);
    hisCVEneDistribTrigged[i] = new TH1D(Form("hisCVEneDistribTrigged_%d",i),Form("hisCVEneDistribTrigged_%d",i),400,0,4000);
  }
  TH1D* hisSciEneDistrib = new TH1D("hisSciEneDistrib","hisSciEneDistrib",400,0,4000);
  TH1D* hisSciEneDistribTrigged = new TH1D("hisSciEneDistribTrigged","hisSciEneDistribTrigged",400,0,4000);

  std::cout<< "LOOP" << std::endl;
  long entries = trin->GetEntries();
  std::cout<< entries << std::endl;
  for( int ievent  = 0; ievent < entries ; ievent++){
    //for( int ievent  = 0; ievent < 100 ; ievent++){
    trin->GetEntry(ievent);
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::list<Pi0>     plist;
    //if((ievent % 1000)== 0){ std::cout<< ievent << "/" << entries << std::endl; }
    dataCopy.reset();
    data.getData(plist);
    dataCopy.setData(plist);

    Int_t hisID = -1;    
    std::list<Pi0>::iterator pit = plist.begin();
    if( (*pit).vz() < 1000 ){ continue; }

    if( hisID >= 0 ){
      hisPi0[hisID]->Fill((*pit).m());    
    }
    hisSciEneDistrib->Fill(SciEne);
    if(SciEne < 425 ){ continue; }
    bool CVTrig = false;
    for( int icv = 0; icv <CVNumber; icv++){
      if( CVEne[icv] > CVThreshold[CVModID[icv]]){
	CVTrig = true;
      }
      hisCVEneDistrib[CVModID[icv]]->Fill(CVEne[icv]);
    }
    if(CVTrig){continue;}

    hisID = 0;
    int nTrig = 0; 
    for( int i = 1; i< 11; i++){
      if( CsiL1TrigCount[i] > CsiL1TrigCountThreshold[i] ){
	//&&CsiL1TrigCount[i] < CsiL1TrigHighThreshold ){
	nTrig++;
      }
    }

    hisL1TrigCount[0]->Fill(nTrig);
    for( int i = 1; i< 11; i++){
      hisL1TrigCount[i]->Fill(CsiL1TrigCount[i]);
    }
    if( nTrig >=2 ){
      if( hisID >= 0){      
	hisSciEneDistribTrigged->Fill(SciEne);
	for( int icv = 0; icv < CVNumber; icv++){
	  hisCVEneDistribTrigged[CVModID[icv]]->Fill(CVEne[icv]);
	}


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
	// Gamma Position Cut 
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

	//std::cout<< "ID Height" << std::endl;
	//std::cout<< ClusterID[0] << "\t" << ClusterHeight[0] << std::endl;
	//std::cout<< ClusterID[1] << "\t" << ClusterHeight[1] << std::endl;
	//std::cout<< MaximumHeight << std::endl;


	double cosTheta = TMath::Abs( x[0]*x[1]+y[0]*y[1] )/TMath::Sqrt((x[0]*x[0]+y[0]*y[1])*(x[1]*x[1]+y[1]*y[1]));
	hisCosTheta[hisID]->Fill(cosTheta);
	double Eg1 = (*pit).g1().e();
	double Eg2 = (*pit).g2().e();
	double gchisq_1 = (*pit).g1().chisq();
	double gchisq_2 = (*pit).g2().chisq();
	double pi0pt    = TMath::Sqrt((*pit).p3()[0]*(*pit).p3()[0]+ (*pit).p3()[1]*(*pit).p3()[1]);
	double pi0Mass  = (*pit).m();

	GammaTimeSigma = 0;
	for( int i  =0; i< 2; i++){
	  GammaTime[i] = 0;
	  GammaY[i] = 0;
	  GammaX[i] = 0;
	  GammaE[i] = 0;
	  GammaCE[i] = 0;
	  GammaTimeSigmaExcept[i] = 0;
	}
	
	GammaTime[0] = (*pit).g1().clusterTimeVec()[0]-showerTimeDelay((*pit),(*pit).g1())-TimeOffset[(*pit).g1().clusterIdVec()[0]];
	GammaCE[0]   = (*pit).g1().clusterEVec()[0];
	GammaE[0]    = (*pit).g1().e();
	GammaX[0]    = (*pit).g1().x();
	GammaY[0]    = (*pit).g1().y();

	GammaTime[1] = (*pit).g2().clusterTimeVec()[0]-showerTimeDelay((*pit),(*pit).g2())-TimeOffset[(*pit).g2().clusterIdVec()[0]];
	GammaCE[1]   = (*pit).g2().clusterEVec()[0];
	GammaE[1]    = (*pit).g2().e();
	GammaX[1]    = (*pit).g2().x();
	GammaY[1]    = (*pit).g2().y();
	
	GammaTimeSigma = TMath::Abs( GammaTime[0]-GammaTime[1]);
	GammaTimeSigmaExcept[0] = GammaTimeSigma;
	GammaTimeSigmaExcept[1] = GammaTimeSigma;
	
	dataCopy.setData(plist);
	trOut->Fill();

	  


	//std::cout<< pi0Mass << std::endl;
	if( Eg1 > 350 &&
	    Eg2 > 200 &&
	    //gchisq_1 < 5 && 
	    //gchisq_2 < 5 &&
	    pi0pt  > 50  &&
	    cosTheta < 0.9 ){
	  hisPi0CutMass[hisID]->Fill((*pit).m());
	  hisPi0E[hisID]->Fill((*pit).e());
	  hisGammaECutHigh[hisID]->Fill((*pit).g1().e());
	  hisGammaECutLow[hisID]->Fill((*pit).g2().e());

	  hisPi0MassGammaEH[hisID]->Fill(Eg1,pi0Mass);
	  hisPi0MassGammaEL[hisID]->Fill(Eg2,pi0Mass);
	  hisPi0MassHeight[hisID]->Fill(MaximumHeight, pi0Mass);
	  hisPi0MassCenterE[hisID]->Fill((*pit).g1().clusterEVec()[0], pi0Mass);
	  hisPi0MassCenterE[hisID]->Fill((*pit).g2().clusterEVec()[0], pi0Mass);

	  if( TMath::Abs((*pit).m()-135)< 10 ){
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
    hisPi0MassGammaEH[i]->Write();
    hisPi0MassGammaEL[i]->Write();
    hisPi0MassHeight[i]->Write();
    hisPi0MassCenterE[i]->Write();
  }

  for( int i = 0; i < 11; i++){
    hisL1TrigCount[i]->Write();
  }
  for( int i = 0; i < 11; i++){
    hisL1TrigCountTrigged[i]->Write();
  }
  for( int i =0; i< 20; i++){
    hisCVEneDistrib[i]->Write();
    hisCVEneDistribTrigged[i]->Write();
  }
  hisSciEneDistrib->Write();
  hisSciEneDistribTrigged->Write();
  trOut->Write();
  tfout->Close();
  return 0;
}
