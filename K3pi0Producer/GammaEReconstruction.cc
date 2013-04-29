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
#include "klong/Klong.h"
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

#include "E14WavReader_V1.h"
#include "L1TrigCounter.h"
#include "EnergyConverter.h"

int main( int argc, char** argv ){
  
  std::string InputFileStr = "/Volume0/Simulation/3pi0/Conv_KL3pi0_FAST_REDUCED_5E6_%d.root";
  TChain* trin = new TChain("T");
  for( int i = 0; i< 10; i++){
    trin->Add(Form(InputFileStr.c_str(),i));
  }
  int EventNumber;
  int CsiNumber;
  int CsiModID[3000];
  double CsiEne[3000];
  double CsiTime[3000];

  int nTrack;
  int pid[20];
  double v[20][3];
  float  ek[20];

  trin->SetBranchAddress("EventNum",&EventNumber);
  trin->SetBranchAddress("CsiNumber",&CsiNumber);
  trin->SetBranchAddress("CsiModID",CsiModID);//CsiNumber
  trin->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  trin->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  trin->SetBranchAddress("nTrack",&nTrack);
  trin->SetBranchAddress("v",v);//nTrack
  trin->SetBranchAddress("pid",pid);//nTrack
  trin->SetBranchAddress("ek",ek);//nTrack

  GammaFinder gFinder;
  ClusterFinder cFinder;
  TFile* tfOut          = new TFile("TestOut_1.root","recreate");
  TH2D*  hisGammaERatio = new TH2D("hisGammaERatio","hisGammaERatio",50,0,2000,200,0,2);
  TH2D*  hisGammaEDist  = new TH2D("hisGammaEDist","hisGammaEDist",50,0,2000,200,0,500);
  TH2D*  hisGammaERatioL = new TH2D("hisGammaERatioL","hisGammaERatioL",50,0,2000,200,0,2);
  TH2D*  hisGammaEDistL  = new TH2D("hisGammaEDistL","hisGammaEDistL",50,0,2000,200,0,500);
  TH2D*  hisGammaERatioS = new TH2D("hisGammaERatioS","hisGammaERatioS",50,0,2000,200,0,2);
  TH2D*  hisGammaEDistS  = new TH2D("hisGammaEDistS","hisGammaEDistS",50,0,2000,200,0,500);

  int nCsI;
  int CsIID[3000];
  double CsIEnergy[3000];
  double CsITime[3000];
  std::list<Cluster>  clist;
  std::list<Gamma>    glist;
  std::vector<Klong>  klvec;
  for( int ievent  = 0; ievent < trin->GetEntries(); ievent++){  
  //for( int ievent  = 267508; ievent < trin->GetEntries(); ievent++){
    if( ievent % 100 == 0 ){
      std::cout<< ievent<<  "/" << trin->GetEntries() << "\n";
    }
    trin->GetEntry(ievent);
    nCsI = 0;
    for( int i = 0; i< 2716; i++){
      CsIID[i] = -1;
      CsIEnergy[i] = 0;
      CsITime[i]   = 0;
    }
    for( int i = 0; i < CsiNumber; i++){
      if( CsiEne[i] < 3 ){ continue; }
      CsIID[nCsI]  = CsiModID[i];
      CsIEnergy[nCsI] = CsiEne[i];
      CsITime[nCsI]   = CsiTime[i];
      nCsI++;
    }
    if( nCsI < 6 ){ continue; }
    if( nCsI > 200 ){ continue; }
    if( nTrack != 6 ){ continue; }
    clist.clear();
    glist.clear();
    klvec.clear();

    /*
    for( int idigi = 0 ; idigi < nCsI; idigi++){
      std::cout << ievent << "\t"<<idigi << ":" <<  nCsI <<"\t"  << CsIID[idigi] << "\t" <<  CsIEnergy[idigi] << "\t" << CsITime[idigi] << std::endl;
    }
    */

    clist = cFinder.findCluster(nCsI,CsIID,CsIEnergy,CsITime);
    //std::cout << "ClFinder" << std::endl;
    gFinder.findGamma(clist,glist);
  
    if( clist.size() < 6 ){ continue; }
    if( glist.size() < 6 ){ continue; }
    if( glist.size() == 6 ){
      if(user_rec(glist,klvec)){
	double tmpGammaE[6]={0};
	double tmpGammaX[6]={0};
	double tmpGammaY[6]={0};
	double tmpGammaR[6][6]={{0}};
	int    minGammaID[6] = {-1};
	double minGammaR[6];
	for( int ip =0; ip <3; ip++){
	  for( int ig = 0; ig < 2; ig++){
	    if( ig ==0 ){
	      tmpGammaE[ip*2+ig] =klvec[0].pi0()[ip].g1().e();
	      tmpGammaX[ip*2+ig] =klvec[0].pi0()[ip].g1().x();
	      tmpGammaY[ip*2+ig] =klvec[0].pi0()[ip].g1().y();
	    }else{
	      tmpGammaE[ip*2+ig] =klvec[0].pi0()[ip].g2().e();
	      tmpGammaX[ip*2+ig] =klvec[0].pi0()[ip].g2().x();
	      tmpGammaY[ip*2+ig] =klvec[0].pi0()[ip].g2().y();
	    }
	    minGammaR[ip*2+ig] = 2000;
	    for( int isg = 0; isg < 6; isg++){
	      tmpGammaR[ip*2+ig][isg] =TMath::Sqrt( TMath::Power(v[isg][0]-tmpGammaX[ip*2+ig],2)+TMath::Power(v[isg][1]-tmpGammaY[ip*2+ig],2));
	      if( minGammaR[ip*2+ig] >  tmpGammaR[ip*2+ig][isg] ){
		minGammaR[ip*2+ig] = tmpGammaR[ip*2+ig][isg];
		minGammaID[ip*2+ig] = isg;
	      }	 
	    }	    
	  }
	}

	/*
	for( int ig = 0; ig< 6; ig++){
	  std::cout<< ig << "\t" << tmpGammaE[ig] << "\t" << ek[ig] << std::endl;
	}
	for( int ig = 0; ig < 6; ig++){
	  std::cout << ig <<"\t";
	  for( int jg = 0; jg < 6; jg++){
	    std::cout <<  tmpGammaR[ig][jg] << "\t";
	  }
	  std::cout << "\n" <<  std::endl;
	}
	*/
	bool tmpFlag = false;
	for( int ig = 0; ig < 6; ig++){
	  if( minGammaR[ig] > 75 ){ tmpFlag  =true; }
	  if( abs(tmpGammaY[ig]) > 550 ){ tmpFlag = true; }
	  double R = sqrt( tmpGammaX[ig]*tmpGammaX[ig] + tmpGammaY[ig]*tmpGammaY[ig] );
	  if( R > 900 ){ tmpFlag =true; }
	  if( abs( tmpGammaX[ig] ) < 150 && abs( tmpGammaY[ig] ) < 150 ){ tmpFlag = true; }
	}
	if( tmpFlag ){ continue; }
	for( int ig = 0; ig < 6; ig++){
	  //std::cout<< ig << "\t" << tmpGammaE[ig] << "\t" <<ek[minGammaID[ig]] << "\t" << minGammaR[ig] << std::endl;
	  hisGammaEDist->Fill(tmpGammaE[ig], minGammaR[ig] );
	  hisGammaERatio->Fill(tmpGammaE[ig], tmpGammaE[ig]/ek[minGammaID[ig]]);
	  if( abs( tmpGammaX[ig] ) < 600 ){
	    hisGammaEDistS->Fill(tmpGammaE[ig], minGammaR[ig] );
	    hisGammaERatioS->Fill(tmpGammaE[ig], tmpGammaE[ig]/ek[minGammaID[ig]]);
	  }else{
	    hisGammaEDistL->Fill(tmpGammaE[ig], minGammaR[ig] );
	    hisGammaERatioL->Fill(tmpGammaE[ig], tmpGammaE[ig]/ek[minGammaID[ig]]);
	  }
	}
	//std::cout<< "/////////////////////////" << std::endl;
	//getchar();	
      }
    }    
  }
  hisGammaERatio->Write();
  hisGammaEDist->Write();
  hisGammaERatioS->Write();
  hisGammaEDistS->Write();
  hisGammaERatioL->Write();
  hisGammaEDistL->Write();
  tfOut->Close();
  return 0;
}


