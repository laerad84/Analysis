#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"

#include "klong/Klong.h"

#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "TGraphErrors.h"
#include "E14ReadSumFile.h"
#include "ReadWavAna.h"
#include "User_Function.h"
//#include "E14_CALIBRATION_CLUSTER_BUILDER_V1/KL_calibration.h"
#include "Calibration.h"
#include "CalibrationTree.h"
#include "GeneralFunctions.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>w
#include <vector>
#include <list>
#include "TMath.h"

//bool user_rec(std::list<Gamma> const &glist, std::vector<Klong> &klVec);
//void user_cut(E14GNAnaDataContainer &data  , std::vector<Klong> &klVec);

Int_t 
main(int argc,char** argv)
{
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Read Arguement
  /////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<< __LINE__ << std::endl;
  int nloop = -1;
  int runNumber=std::atoi( argv[1]);
  int iterationNumber;
  std::string inputFilename;
  std::string outputFilename;
  std::string calibrationFilename;
  std::string TempCalibrationFilename;
  std::cout<< __LINE__ << std::endl;
  //std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  //std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  //std::string ROOTFILE_3PI0CALIBRATION = std::getenv( "ROOTFILE_3PI0CALIBRATION");
  //std::string ROOTFILE_SIMCONV = "/group/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/3pi0Run/SIM3PI0";
  std::cout<< __LINE__ << std::endl;
  std::string ROOTFILE_SIMCONV = "~/KL3pi0";
  std::string path;
  std::string InitalCalibrationFilename;
  std::cout<< __LINE__ << std::endl;
  inputFilename = Form("%s/Sim_e14_KL3pi0_KL_RES_LY_pe_1E8_%d.root",ROOTFILE_SIMCONV.c_str(), runNumber);
  //inputFilename       = Form("%s/Sim3pi0_wav_fast_KL_RES_LY_pe_5E6_%d_Calibration.root",ROOTFILE_SIMCONV.c_str(),runNumber);
  //inputFilename       = Form("%s/Sim3pi0_wav_fast_KL_RES_LY_pe_5E6_%d_Calibration_mis_1.root",ROOTFILE_SIMCONV.c_str(),runNumber);// Test MisCalibration

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Input RootFile
  /////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout<< __LINE__ << std::endl;
  //E14ReadSumFile* ReadSum = new  E14ReadSumFile();
  TFile* tfinput = new TFile(inputFilename.c_str());
  if(!tfinput->IsOpen()){ std::cout<< Form("%s is not opened",tfinput->GetName()) << std::endl; return -1;}

  TChain* ch = new TChain("T");
  for( int irun = 0; irun< 10; irun++){
    ch->Add(Form("%s/Sim_e14_KL3pi0_KL_RES_LY_pe_1E8_%d.root",ROOTFILE_SIMCONV.c_str(), irun));
  }
  E14GNAnaDataContainer data;

  Int_t CsiNumber;
  Double_t CsiEne[2716];//CsiNumber
  Double_t CsiTime[2716];//CsiNumber
  Double_t CsiSignal[2716];//CsiNumber
  Int_t CsiModID[2716];//CsiNumber

  Int_t nTrack;
  Int_t pid[20];//nTrack
  Double_t end_v[20][3];//nTrack
  Float_t ek[20];//nTrack

  data.setBranchAddress(ch);
  ch->SetBranchAddress("CsiNumber",&CsiNumber );
  ch->SetBranchAddress("CsiTime"  ,CsiTime);//CsiNumber
  ch->SetBranchAddress("CsiModID" ,CsiModID);//CsiNumber
  ch->SetBranchAddress("CsiEne"   ,CsiEne);//CsiNumber
  ch->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber

  ch->SetBranchAddress("nTrack",&nTrack);
  ch->SetBranchAddress("pid",pid);//nTrack
  ch->SetBranchAddress("end_v",end_v);//nTrack
  ch->SetBranchAddress("ek",ek);//nTrack


  TFile* tfcorr = new TFile("GammaECorrection.root");
  TGraphErrors* fcorr = (TGraphErrors*)tfcorr->Get("ECorrection");
  


  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Output RootFile
  /////////////////////////////////////////////////////////////////////////////////////////////////

  TFile* tfout = new TFile("GammaR.root","recreate");
  TH1D* hisGammaR = new TH1D("hisGammaR","hisgammaR",400,0,2000);
  TH2D* hisGammaERatio = new TH2D("hisGammaERatio","hisGammaERatio",50,0,2000,100,0.5,1.5);
  TH2D* hisKlongDecayPosition = new TH2D("hisKlongDecayPosition","hisKlongDecayPosition",120,0,6000,100,-200,200);
  TH2D* hisGammaERatioLog = new TH2D("hisGammaERatioLog","hisGammaERatioLog",25,GenLogArray(25,1,2000),100,GenLinArray(100,0.5,1.5));
  //for(int ievt=0;ievt< 10;ievt++){
  for(int ievt=0;ievt< ch->GetEntries();ievt++){
    ch->GetEntry(ievt);
    if( nTrack != 10 ){ continue; }
    
    //std::cout<< "Build" << std::endl;
    std::list<Cluster> clist;
    std::list<Gamma>   glist; 
    std::vector<Klong> klVec;
    data.getData(klVec);
    Double_t ge[6]={0};
    Double_t gx[6]={0};
    Double_t gy[6]={0};
    Double_t gR[6]={0};
    if ( klVec.size() == 0 || klVec[0].pi0().size() != 3 ){ continue; }
    if( klVec[0].chisqZ() > 5 ){ continue; }
    int ig=0;
    for( int i=0; i< nTrack; i++){
      if( pid[i] !=22 ){ continue; }
      
      //std::cout<< pid[i] << "\t" << ek[i] << "\t" << end_v[i][0] << "\t" << end_v[i][1]<< std::endl;
      ge[ig]=ek[i];
      gx[ig]=end_v[i][0];
      gy[ig]=end_v[i][1];      
      ig++;
    }
    
    hisKlongDecayPosition->Fill(end_v[0][2],klVec[0].vz() - end_v[0][2]);
    
    bool GammaPosCut = false;
    for( int i = 0; i< klVec[0].pi0().size(); i++){
      if( TMath::Abs(klVec[0].pi0()[i].g1().pos()[0]) < 150 &&
	  TMath::Abs(klVec[0].pi0()[i].g1().pos()[1]) < 150 ){
	GammaPosCut = true;
      }
      if( TMath::Abs(klVec[0].pi0()[i].g2().pos()[0]) < 150 &&
	  TMath::Abs(klVec[0].pi0()[i].g2().pos()[1]) < 150 ){
	GammaPosCut = true;
      }
      if( klVec[0].pi0()[i].g1().pos().perp() > 850 ) { GammaPosCut = true; }
      if( klVec[0].pi0()[i].g2().pos().perp() > 850 ) { GammaPosCut = true; }
    }
    if( GammaPosCut ){ continue; }
    for( int i =0;i< klVec[0].pi0().size(); i++){
      //std::cout<< klVec[0].pi0()[i].g1().pos()[0] << "\t" << klVec[0].pi0()[i].g1().pos()[1] << std::endl;
      //std::cout<< klVec[0].pi0()[i].g2().pos()[0] << "\t" << klVec[0].pi0()[i].g2().pos()[1] << std::endl;
      for( int j =0 ; j< 6; j++){
	Double_t R1 = TMath::Sqrt( pow(klVec[0].pi0()[i].g1().pos()[0] - gx[j],2)+pow(klVec[0].pi0()[i].g1().pos()[1] - gy[j],2));
	Double_t R2 = TMath::Sqrt( pow(klVec[0].pi0()[i].g2().pos()[0] - gx[j],2)+pow(klVec[0].pi0()[i].g2().pos()[1] - gy[j],2));
	if( R1 < 80 ){
	  hisGammaERatio->Fill(ge[j],ge[j]/klVec[0].pi0()[i].g1().e()); 
	  hisGammaERatioLog->Fill(ge[j],fcorr->Eval(ge[j],0,"S")*ge[j]/klVec[0].pi0()[i].g1().e()); 
	}
	if( R2 < 80 ){
	  hisGammaERatio->Fill(ge[j],ge[j]/klVec[0].pi0()[i].g2().e()); 
	  hisGammaERatioLog->Fill(ge[j],fcorr->Eval(ge[j],0,"S")*ge[j]/klVec[0].pi0()[i].g2().e()); 
	}
	hisGammaR->Fill(R1);
	hisGammaR->Fill(R2);
      }
    }
  }
  hisKlongDecayPosition->Write();
  hisGammaERatioLog->Write();
  hisGammaR->Write();
  hisGammaERatio->Write();
  tfout->Close();
  std::cout<<"finish!"<<std::endl;  
}
