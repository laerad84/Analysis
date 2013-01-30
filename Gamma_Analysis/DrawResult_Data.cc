#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TApplication.h"
#include "ClusterTimeReader.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TProfile.h"
double AdjFunc( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double value = p0 + p1*exp(p2*x0) + p3*exp(p4*x0);
  return value;
}
double AdjFunc1( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double value = p0 + p1*exp(p2*x0);
  return value;
}
double AdjFuncTest( double* x, double* par){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double value = p0 + p1*exp(p2*x0) + p3*exp(p4*x0);
  return value;

}

int main( int argc, char** argv ){
  
  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string iFileForm = "%s/Data_All_Cal_FNL_COS.root";//ROOTFILE_GAMMACLUS
  std::string oFileForm = "ClusterTimeStructure_Data_FNL_COS_2.root";

  TF1* TimeAdjFunc = new TF1("TimeAdjFunc",AdjFunc, 0, 2000,5);
  //Double_t Par[5] = {-0.0905327,1.54915,-0.114423,0.0758477,0.00487457};
  //Double_t ParErrors[5] = {0.00245834,0.0263153,0.00188467,0.007594,4.81501e-05};
  //Double_t Par[5] = {-0.105097,1.52645,-0.10655,0.0620572,0.00910542};
  //Double_t ParErrors[5]={0.00186334,0.0135129,0.000805812,0.000167933,1.18097e-05};
  //Double_t Par[6] = {-0.0644067,1.1759,-0.165316,0.7*0.0570758,0.0049958,0.001};
  //Double_t ParErrors[5] = {0.00203663,0.0418876,0.00542069,0.00221644,4.08696e-05}; 

  //Double_t Par[6] = {-0.0644067,3*1.1759,-0.165316,2.92809e-1,2.57735e-03};
  //Double_t ParErrors[5] = {0.00203663,0.0418876,0.00542069,0.00221644,4.08696e-05};
  double Par[5] = {-0.261894,1.35789,-0.106939,0.094532,0.00395812};
  TimeAdjFunc->SetParameters(Par);
  //TimeAdjFunc->SetParErrors(ParErrors);
 
  TF1* TimeAdjFunc1 = new TF1("TimeAdjFunc",AdjFunc1,0,2000,3);
  Double_t Par1[3] = {5.70121e-01,-2.92809e-01,2.57735e-03};
  TimeAdjFunc1->SetParameters(Par1);

  TFile* tf = new TFile(Form(iFileForm.c_str(),ROOTFILE_GAMMACLUS.c_str()));
  TTree* tr = (TTree*)tf->Get("trCluster");
  ClusterTimeReader* reader = new ClusterTimeReader(tr);
  Int_t nEntries = reader->fChain->GetEntries();
  
  TFile* tfout = new TFile(oFileForm.c_str(),"recreate");  

  const int nE = 6;
  const int nTheta = 8;
  TH2D* hisPhiPhi[nE-1][nTheta -1 ];
  TH2D* hisRPhi[nE-1][nTheta -1 ];
  TH2D* hisRPhiChisqN[nE-1][nTheta -1];
  TH2D* hisRPhiE[nE-1][nTheta -1];
  TProfile* profET_D_All;
  TProfile* profET_D_All_min;
  TProfile* profET_D_All_max;


  TH2D* hisRT[nE-1][nTheta-1];
  TH2D* hisRE[nE-1][nTheta-1];
  TH2D* hisDT[nE-1][nTheta-1];
  TH2D* hisDE[nE-1][nTheta-1];
  TProfile2D* profRDT[nE-1][nTheta-1];
  TProfile2D* profRDE[nE-1][nTheta-1];
  TProfile2D* profRDT_MIP[nE-1][nTheta-1];
  TProfile2D* profRDE_MIP[nE-1][nTheta-1];
  TProfile2D* profRDT_Low[nE-1][nTheta-1];
  TProfile2D* profRDE_Low[nE-1][nTheta-1];
  TProfile2D* profRDT_High[nE-1][nTheta-1];
  TProfile2D* profRDE_High[nE-1][nTheta-1];
  TProfile2D* profRDE_INV[nE-1][nTheta-1];
  TProfile2D* profRDT_INV[nE-1][nTheta-1];
  TProfile*   profRT_INV[nE-1][nTheta-1];
  TProfile*   profDT_INV[nE-1][nTheta-1];

  TProfile2D* profRDT_NR[nE-1][nTheta-1];
  TProfile* profDT_NR[nE-1][nTheta-1];
  TProfile* profRT_NR[nE-1][nTheta-1];
  TProfile2D* profRDE_NR[nE-1][nTheta-1];

  TProfile2D* profTimeRD[nE-1][nTheta-1];
  TProfile2D* profEET[nE-1][nTheta-1];
  TProfile2D* profEET_D[nE-1][nTheta-1];

  TProfile*   profRTAll[nE-1][nTheta-1];
  TProfile*   profRT_high[nE-1][nTheta-1];
  TProfile*   profRT_low[nE-1][nTheta-1];
  TProfile*   profDTAll[nE-1][nTheta-1];
  TProfile*   profDT_high[nE-1][nTheta-1];
  TProfile*   profDT_low[nE-1][nTheta-1];
  TProfile*   profDEAll[nE-1][nTheta-1];
  TProfile*   profDE_high[nE-1][nTheta-1];
  TProfile*   profDE_low[nE-1][nTheta-1];
  TProfile*   profREAll[nE-1][nTheta-1];
  TProfile*   profRE_high[nE-1][nTheta-1];
  TProfile*   profRE_low[nE-1][nTheta-1];
  
  TProfile*   profET_D[nE-1][nTheta-1];
  TProfile*   profET_D0[nE-1][nTheta-1];
  TProfile*   profET_D1[nE-1][nTheta-1];
  TProfile*   profET_D2[nE-1][nTheta-1];
  TProfile*   profET_D3[nE-1][nTheta-1];

  TProfile*   profRT_MIP[nE-1][nTheta-1];
  TProfile*   profDT_MIP[nE-1][nTheta-1];
  
  TProfile2D* profRTE[nE-1][nTheta-1];
  TProfile2D* profDTE[nE-1][nTheta-1];
  TProfile2D* profRDT_NRCL[nTheta-1];
  TProfile*   profRT_NRCL[nTheta-1];
  TProfile*   profDT_NRCL[nTheta-1];
  

  const int nBinsR   = 41;
  const int nBinsPhi = 41;
  
  TH1D* hisRPhiTime[nE-1][nTheta-1][41][41];
  TH1D* hisRPhiEnergy[nE-1][nTheta-1][41][41];

  Int_t EArr[nE] = {100,300,500,700,900,1100};
  Int_t ThetaArr[nTheta] = {10,15,20,25,30,35,40,45};
  
  Int_t RArr[ nBinsR ];
  Int_t PhiArr[ nBinsPhi ];
  for( int kIndex  = 0; kIndex < nBinsR; kIndex++){
    RArr[ kIndex ] = -102.5+ 5*kIndex;
  }
  for( int lIndex = 0; lIndex < nBinsPhi; lIndex++){
    PhiArr[ lIndex ] = -102.5+5*lIndex;
  }
  
  Double_t width   = 12.5;
  Int_t    nBinsRD = 19;
  Double_t RDMax   = width*nBinsRD/2.;
  Double_t RDMin   = -1*width*nBinsRD/2.;


  for( int iIndex = 0; iIndex < nTheta-1; iIndex++){
    profRDT_NRCL[iIndex] = new TProfile2D(Form("profRDT_NRCL_Theta_%d",iIndex),
					  Form("profRDT_NRCL_Theta_%d",ThetaArr[iIndex]),
					  nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
    profRT_NRCL[iIndex] = new TProfile(Form("profRT_NRCL_Theta_%d",iIndex),
				       Form("profRT_NRCL_Theta_%d",iIndex),
				       nBinsRD,RDMin,RDMax);
    profDT_NRCL[iIndex] = new TProfile(Form("profDT_NRCL_Theta_%d",iIndex),
				       Form("profDT_NRCL_Theta_%d",iIndex),
				       nBinsRD,RDMin,RDMax);
  }

  profET_D_All = new TProfile("profET_D_All","profET_D_All",50,0,800);
  profET_D_All_min = new TProfile("profET_D_All_min","profET_D_All_min",160,0,800);
  profET_D_All_max = new TProfile("profET_D_All_max","profET_D_All_max",32,0,800);
  for( int iIndex = 0; iIndex < nE-1; iIndex++){
    for( int jIndex = 0; jIndex < nTheta-1; jIndex++){
      for( int kIndex  = 0; kIndex  < nBinsR; kIndex++){
	for( int lIndex = 0; lIndex < nBinsPhi; lIndex++){
	  hisRPhiTime[ iIndex][jIndex][kIndex][lIndex] 
	    = new TH1D(Form("hisRPhiTime_%d_%d_%d_%d",iIndex,jIndex,kIndex,lIndex),
		       Form("hisRPhiTime_E_%d_%d_Theta_%d_%d_xy_%d_%d",
			    EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1],
			    RArr[kIndex],PhiArr[lIndex]),
		       200,-10,10);
	  hisRPhiEnergy[iIndex][jIndex][kIndex][lIndex]
	    = new TH1D(Form("hisRPhiEnergy_%d_%d_%d_%d",iIndex,jIndex,kIndex,lIndex),
		       Form("hisRPhiEnergy_E_%d_%d_Theta_%d_%d_xy_%d_%d",
			    EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1],
			    RArr[kIndex],PhiArr[lIndex]),
		       100,0,1000);					       
	}
      }
      profRTE[iIndex][jIndex] = new TProfile2D(Form("profRTE_E_%d_Theta_%d",iIndex,jIndex),
					       Form("profRTE_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
					       nBinsRD,RDMin,RDMax,50,0,400);

      profDTE[iIndex][jIndex] = new TProfile2D(Form("profDTE_E_%d_Theta_%d",iIndex,jIndex),
					       Form("profDTE_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
					       nBinsRD,RDMin,RDMax,50,0,400);


      profRT_MIP[iIndex][jIndex] = new TProfile(Form("profRT_MIP_E_%d_Theta_%d",iIndex,jIndex),
						Form("profRT_MIP_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      
      profDT_MIP[iIndex][jIndex] = new TProfile(Form("profDT_MIP_E_%d_Theta_%d",iIndex,jIndex),
						Form("profDT_MIP_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      
      hisRT[iIndex][jIndex] = new TH2D(Form("hisRT_E_%d_Theta_%d",iIndex,jIndex),
				       Form("hisRT_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
				       nBinsRD,RDMin,RDMax,200,-10,10);
      hisRE[iIndex][jIndex] = new TH2D(Form("hisRE_E_%d_Theta_%d",iIndex,jIndex),
				       Form("hisRE_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
				       nBinsRD,RDMin,RDMax,200,0,800);
      hisDT[iIndex][jIndex] = new TH2D(Form("hisDT_E_%d_Theta_%d",iIndex,jIndex),
				       Form("hisDT_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
				       nBinsRD,RDMin,RDMax,200,-10,10);
      hisDE[iIndex][jIndex] = new TH2D(Form("hisDE_E_%d_Theta_%d",iIndex,jIndex),
				       Form("hisDE_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
				       nBinsRD,RDMin,RDMax,200,0,800);
      profRDT[iIndex][jIndex] = new TProfile2D(Form("profRDT_E_%d_Theta_%d",iIndex,jIndex),
					       Form("profRDT_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
					       nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDE[iIndex][jIndex] = new TProfile2D(Form("profRDE_E_%d_Theta_%d",iIndex,jIndex),
					       Form("profRDE_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
					       nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDT_MIP[iIndex][jIndex] = new TProfile2D(Form("profRDT_MIP_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDT_MIP_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDE_MIP[iIndex][jIndex] = new TProfile2D(Form("profRDE_MIP_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDE_MIP_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDT_High[iIndex][jIndex] = new TProfile2D(Form("profRDT_High_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDT_High_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDE_High[iIndex][jIndex] = new TProfile2D(Form("profRDE_High_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDE_High_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDT_Low[iIndex][jIndex] = new TProfile2D(Form("profRDT_Low_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDT_Low_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDE_Low[iIndex][jIndex] = new TProfile2D(Form("profRDE_Low_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDE_Low_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDT_INV[iIndex][jIndex] = new TProfile2D(Form("profRDT_INV_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDT_INV_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDE_INV[iIndex][jIndex] = new TProfile2D(Form("profRDE_INV_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDE_INV_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRDT_NR[iIndex][jIndex] = new TProfile2D(Form("profRDT_NR_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDT_NR_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profRT_NR[iIndex][jIndex] = new TProfile(Form("profRT_NR_E_%d_Theta_%d",iIndex,jIndex),
					       Form("profRT_NR_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
					       nBinsRD,RDMin,RDMax);
      profDT_NR[iIndex][jIndex] = new TProfile(Form("profDT_NR_E_%d_Theta_%d",iIndex,jIndex),
					       Form("profDT_NR_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
					       nBinsRD,RDMin,RDMax);
      profRDE_NR[iIndex][jIndex] = new TProfile2D(Form("profRDE_NR_E_%d_Theta_%d",iIndex,jIndex),
						   Form("profRDE_NR_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      profDT_INV[iIndex][jIndex] = new TProfile(Form("profDT_INV_E_%d_Theta_%d",iIndex,jIndex),
						Form("profDT_INV_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      profRT_INV[iIndex][jIndex] = new TProfile(Form("profRT_INV_E_%d_Theta_%d",iIndex,jIndex),
						Form("profRT_INV_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      
      profREAll[iIndex][jIndex]  = new TProfile(Form("profREAll_E_%d_Theta_%d",iIndex,jIndex),
						Form("profREAll_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      profRE_low[iIndex][jIndex] = new TProfile(Form("profRE_low_E_%d_Theta_%d",iIndex,jIndex),
						Form("profRE_low_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      profRE_high[iIndex][jIndex]= new TProfile(Form("profRE_high_E_%d_Theta_%d",iIndex,jIndex),
						Form("profRE_high_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      profDEAll[iIndex][jIndex]  = new TProfile(Form("profDEAll_E_%d_Theta_%d",iIndex,jIndex),
						Form("profDEAll_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      profDE_low[iIndex][jIndex] = new TProfile(Form("profDE_low_E_%d_Theta_%d",iIndex,jIndex),
						Form("profDE_low_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      profDE_high[iIndex][jIndex]= new TProfile(Form("profDE_high_E_%d_Theta_%d",iIndex,jIndex),
						Form("profDE_high_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      profRTAll[iIndex][jIndex]  = new TProfile(Form("profRTAll_E_%d_Theta_%d",iIndex,jIndex),
						Form("profRTAll_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);

      profRT_high[iIndex][jIndex]= new TProfile(Form("profRT_high_E_%d_Theta_%d",iIndex,jIndex),
						Form("profRThigh_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      
      profRT_low[iIndex][jIndex] = new TProfile(Form("profRT_low_E_%d_Theta_%d",iIndex,jIndex),
						Form("profRT_low_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      
      profDTAll[iIndex][jIndex]  = new TProfile(Form("profDTAll_E_%d_Theta_%d",iIndex,jIndex),
						Form("profDTAll_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
					       nBinsRD,RDMin,RDMax);
      
      profDT_high[iIndex][jIndex]= new TProfile(Form("profDT_high_E_%d_Theta_%d",iIndex,jIndex),
						Form("profDThigh_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      
      profDT_low[iIndex][jIndex] = new TProfile(Form("profDT_low_E_%d_Theta_%d",iIndex,jIndex),
						Form("profDT_low_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						nBinsRD,RDMin,RDMax);
      
      profEET[iIndex][jIndex]    = new TProfile2D(Form("profEET_E_%d_Theta_%d",iIndex,jIndex),
						  Form("profEET_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						  20,0,800,20,0,800);
      profEET_D[iIndex][jIndex]    = new TProfile2D(Form("profEET_D_E_%d_Theta_%d",iIndex,jIndex),
						  Form("profEET_D_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						  20,0,800,20,0,800);
      profET_D[iIndex][jIndex]    = new TProfile(Form("profET_D_E_%d_Theta_%d",iIndex,jIndex),
						 Form("profET_D_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						 50,0,800);
      profET_D0[iIndex][jIndex]    = new TProfile(Form("profET_D0_E_%d_Theta_%d",iIndex,jIndex),
						 Form("profET_D0_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						 50,0,800);
      profET_D1[iIndex][jIndex]    = new TProfile(Form("profET_D1_E_%d_Theta_%d",iIndex,jIndex),
						 Form("profET_D1_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						 50,0,800);
      profET_D2[iIndex][jIndex]    = new TProfile(Form("profET_D2_E_%d_Theta_%d",iIndex,jIndex),
						 Form("profET_D2_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						 50,0,800);
      profET_D3[iIndex][jIndex]    = new TProfile(Form("profET_D3_E_%d_Theta_%d",iIndex,jIndex),
						 Form("profET_D3_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						 50,0,800);
      profTimeRD[iIndex][jIndex] = new TProfile2D( Form("hisTimeRD_E_%d_Theta_%d",iIndex,jIndex),
						   Form("hisTimeRD_E_%d_%d_Theta_%d_%d",EArr[iIndex],EArr[iIndex+1],ThetaArr[jIndex],ThetaArr[jIndex+1]),
						   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      hisRPhi[iIndex][jIndex]   = new TH2D(Form("hisRPhi_E_%d_Theta_%d",iIndex,jIndex),
					   Form("hisRPhi_E_%d_%d_Theta_%d_%d",EArr[iIndex], EArr[iIndex+1],ThetaArr[jIndex], ThetaArr[jIndex+1]),
					   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      hisRPhiE[iIndex][jIndex]  = new TH2D(Form("hisRPhiE_E_%d_Theta_%d",
						iIndex,jIndex),
					   Form("hisRPhiE_E_%d_%d_Theta_%d_%d",
						EArr[iIndex], EArr[iIndex+1],
						ThetaArr[jIndex], ThetaArr[jIndex+1]),
					   nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      
      hisRPhiChisqN[iIndex][jIndex] = new TH2D(Form("hisRPhiChisqN_E_%d_Theta_%d",
						    iIndex,jIndex),
					       Form("hisRThetaChisqN_E_%d_%d_Theta_%d_%d",
						    EArr[iIndex], EArr[iIndex+1],
						    ThetaArr[jIndex], ThetaArr[jIndex+1]),
					       nBinsRD,RDMin,RDMax,nBinsRD,RDMin,RDMax);
      hisPhiPhi[iIndex][jIndex] = new TH2D(Form("hisPhiPhi_E_%d_Theta_%d",
						iIndex,jIndex),
					   Form("hisPhiPhi_E_%d_%d_Theta_%d_%d",
						EArr[iIndex],EArr[iIndex+1],
						ThetaArr[jIndex],ThetaArr[jIndex+1]),
					   60,-1*TMath::Pi(),TMath::Pi(),60,-1*TMath::Pi(),TMath::Pi());
    }
  }

				
  ///////////////////////////////////////////////////////////////////////
  /// Cut Condition 
  ///////////////////////////////////////////////////////////////////////
  Double_t OuterRadCut = 500;
  Double_t InnerRadCut = 250;
  Double_t LowEnergyCut= 3;
  Int_t    ClusterSizeCut = 6;



  for( int evtIndex = 0; evtIndex < nEntries; evtIndex++){
    reader->GetEntry( evtIndex );
    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      if( reader->ClusterR[clusterIndex] > OuterRadCut ){ continue; }
      if( reader->ClusterR[clusterIndex] < InnerRadCut ){ continue; }
      if( reader->nCrystal[clusterIndex] < ClusterSizeCut ){ continue; }

      bool bAbort = true; 
      Int_t EnergyIndex = 0;
      Int_t ThetaIndex  = 0;
      for( int EIndex = 0; EIndex < nE-1; EIndex++){
	if( reader->ClusterEnergy[clusterIndex] >  EArr[EIndex] &&
	    reader->ClusterEnergy[clusterIndex] <= EArr[EIndex +1] ){
	  EnergyIndex  = EIndex;
	  bAbort = false;
	  break;
	}
      }
      if( bAbort ) { continue;}
      bAbort = true; 
      for( int TIndex = 0; TIndex < nTheta -1 ; TIndex++){
	if( reader->ClusterTheta[clusterIndex]*180./TMath::Pi() >  ThetaArr[TIndex] &&
	    reader->ClusterTheta[clusterIndex]*180./TMath::Pi() <= ThetaArr[TIndex+1]){
	  ThetaIndex = TIndex;
	  bAbort = false;
	  break;
	}
      }
      if( bAbort ){ continue; }

      bool blr = true;// view from downstream , l is false( -1,0 ), r is true( 1 )//
      if( TMath::Cos( reader->ClusterPhi[ clusterIndex ] ) < 0 ){ 
	blr = false;
      }else{ blr = true; }
            
      Double_t RCenterCrystal = 1000;
      Double_t ECenterCrystal = 0;
      Double_t TCenterCrystal = 0;
      Double_t RCenter=0;
      Double_t DCenter=0;
      for( int crystalIndex = 0; crystalIndex < reader->nCrystal[ clusterIndex]; crystalIndex++){
	Double_t RadinCluster = reader->CrystalR[clusterIndex][crystalIndex];
	Double_t RinCluster = reader->CrystalR[clusterIndex][crystalIndex]*TMath::Cos(reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t DinCluster = reader->CrystalR[clusterIndex][crystalIndex]*TMath::Sin(reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t EinCluster = reader->CrystalEnergy[clusterIndex][crystalIndex];
	// Time-Energy relation fixed ? // 
	Double_t TinCluster = reader->CrystalT[clusterIndex][crystalIndex]-TimeAdjFunc->Eval(EinCluster);
	if( RadinCluster < TMath::Abs(RCenterCrystal)){ 
	  RCenterCrystal = RadinCluster;
	  ECenterCrystal = EinCluster;
	  TCenterCrystal = TinCluster-TimeAdjFunc->Eval(ECenterCrystal);
	  RCenter        = RinCluster;
	  DCenter        = DinCluster;
	}
      }
      /*
      if(TMath::Abs(RCenter) > 12.5/2. || TMath::Abs(DCenter) > 12.5/2.){
	continue;
      }
      */
      Double_t CutValue = 25;
      for( int crystalIndex  = 0; crystalIndex < reader->nCrystal[ clusterIndex ]; crystalIndex++){       
	if( reader->CrystalEnergy[clusterIndex][crystalIndex] <= LowEnergyCut ) { continue; }
	//if( reader->CrystalEnergy[clusterIndex][crystalIndex] > 400 ){ continue; }
	Double_t RinCluster = reader->CrystalR[clusterIndex][crystalIndex]*TMath::Cos(reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t DinCluster = reader->CrystalR[clusterIndex][crystalIndex]*TMath::Sin(reader->CrystalPhi[clusterIndex][crystalIndex]);
	Double_t EinCluster = reader->CrystalEnergy[clusterIndex][crystalIndex];
	// Time-Energy relation fixed ? // 
	Double_t TinCluster = reader->CrystalT[clusterIndex][crystalIndex]-TCenterCrystal-TimeAdjFunc->Eval(EinCluster);
	if( EinCluster == 0){ continue; }

	hisPhiPhi[EnergyIndex][ThetaIndex]->Fill( reader->ClusterPhi[clusterIndex], reader->CrystalPhi[clusterIndex][crystalIndex]);
	if(reader->CrystalR[clusterIndex][crystalIndex]  > RCenterCrystal ){
	  profRDT[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);
	  profRDE[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,EinCluster);

	  if( EinCluster > ECenterCrystal*0.8 && EinCluster < ECenterCrystal*1.2){
	    profRDT_NR[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);
	    profRDE_NR[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,EinCluster);
	    profRDT_NRCL[ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);
	    
	    if( TMath::Abs(DinCluster)< CutValue ){
	      profRT_NRCL[ThetaIndex]->Fill(RinCluster,TinCluster);
	      profRT_NR[EnergyIndex][ThetaIndex]->Fill(RinCluster,TinCluster);
	    }
	    if( TMath::Abs(RinCluster)< CutValue ){
	      profDT_NRCL[ThetaIndex]->Fill(DinCluster,TinCluster);
	      profDT_NR[EnergyIndex][ThetaIndex]->Fill(DinCluster,TinCluster);
	    }	    
	  }

	  if( EinCluster > ECenterCrystal ){
	    profRDT_INV[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);
	    profRDE_INV[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,EinCluster);
	    if( TMath::Abs(DinCluster)< CutValue ){
	      profRT_INV[EnergyIndex][ThetaIndex]->Fill(RinCluster,TinCluster);
	    }
	    if( TMath::Abs(RinCluster)< CutValue ){
	      profDT_INV[EnergyIndex][ThetaIndex]->Fill(DinCluster,TinCluster);
	    }
	  }

	  if( EinCluster < 24 ){
	    profRDE_MIP[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,EinCluster);
	    profRDT_MIP[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);	    
	  }
	  if( ECenterCrystal > reader->ClusterEnergy[crystalIndex] *0.3){
	    profRDT_High[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);
	    profRDE_High[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,EinCluster);	    
	  }else{
	    profRDT_Low[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);
	    profRDE_Low[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,EinCluster);	    
	  }

	  ////////////////// Fill RT ///////////////////	  
	  if( TMath::Abs(DinCluster) < CutValue ){
	    hisRE[EnergyIndex][ThetaIndex]->Fill(RinCluster,EinCluster);
	    hisRT[EnergyIndex][ThetaIndex]->Fill(RinCluster,TinCluster);
	    profRTAll[EnergyIndex][ThetaIndex]->Fill(RinCluster,TinCluster);
	    profREAll[EnergyIndex][ThetaIndex]->Fill(RinCluster,EinCluster);
	    profRTE[EnergyIndex][ThetaIndex]->Fill(RinCluster,EinCluster,TinCluster);
	    if( EinCluster < 24 && EinCluster > 6){
	      profRT_MIP[EnergyIndex][ThetaIndex]->Fill(RinCluster,TinCluster); 
	    }
	    if( ECenterCrystal > 0.3*reader->ClusterEnergy[clusterIndex]){
	      profRT_high[EnergyIndex][ThetaIndex]->Fill(RinCluster,TinCluster);
	      profRE_high[EnergyIndex][ThetaIndex]->Fill(RinCluster,EinCluster);
	    }else{
	      profRT_low[EnergyIndex][ThetaIndex]->Fill(RinCluster,TinCluster);
	      profRE_low[EnergyIndex][ThetaIndex]->Fill(RinCluster,EinCluster);
	    }
	  }	  
	  ////////////////// Fill DT ///////////////////
	  if( TMath::Abs(RinCluster) < CutValue ){
	    hisDE[EnergyIndex][ThetaIndex]->Fill(DinCluster,EinCluster);
	    hisDT[EnergyIndex][ThetaIndex]->Fill(DinCluster,TinCluster);
	    profDTAll[EnergyIndex][ThetaIndex]->Fill(DinCluster,TinCluster);
	    profDEAll[EnergyIndex][ThetaIndex]->Fill(DinCluster,EinCluster); 
	    profDTE[EnergyIndex][ThetaIndex]->Fill(DinCluster,EinCluster,TinCluster);
	    if( EinCluster < 40 && EinCluster > 10){
	      profDT_MIP[EnergyIndex][ThetaIndex]->Fill(DinCluster,TinCluster);
	    }
	    if( ECenterCrystal > 0.3*reader->ClusterEnergy[clusterIndex]){
	      profDT_high[EnergyIndex][ThetaIndex]->Fill(DinCluster,TinCluster);
	      profDE_high[EnergyIndex][ThetaIndex]->Fill(DinCluster,EinCluster);
	    }else{
	      profDT_low[EnergyIndex][ThetaIndex]->Fill(DinCluster,TinCluster);
	      profDE_low[EnergyIndex][ThetaIndex]->Fill(DinCluster,EinCluster);
	    }
	    hisDT[EnergyIndex][ThetaIndex]->Fill( DinCluster,TinCluster);
	    if(TMath::Abs( DinCluster) < 25 ){
	      profEET_D[EnergyIndex][ThetaIndex]->Fill(ECenterCrystal,EinCluster,TinCluster);
	      if( ThetaIndex <3 ){
		profET_D_All->Fill(ECenterCrystal,TinCluster);
		profET_D_All_min->Fill(ECenterCrystal,TinCluster);
		profET_D_All_max->Fill(ECenterCrystal,TinCluster);

	      }
	      profET_D[EnergyIndex][ThetaIndex]->Fill(ECenterCrystal,TinCluster);
	      if( EinCluster <100 ){
		profET_D0[EnergyIndex][ThetaIndex]->Fill(ECenterCrystal,TinCluster);
	      }else if( EinCluster <200 ){
		profET_D1[EnergyIndex][ThetaIndex]->Fill(ECenterCrystal,TinCluster);
	      }else if( EinCluster <300 ){
		profET_D2[EnergyIndex][ThetaIndex]->Fill(ECenterCrystal,TinCluster);
	      }else if( EinCluster <400 ){
		profET_D3[EnergyIndex][ThetaIndex]->Fill(ECenterCrystal,TinCluster);
	      }

	    }
	  }	  
	  profTimeRD[EnergyIndex][ThetaIndex]->Fill(RinCluster,DinCluster,TinCluster);
	  if(TMath::Abs(RinCluster) < 25 && TMath::Abs(DinCluster) < 25 ){
	    profEET[EnergyIndex][ThetaIndex]->Fill(ECenterCrystal,EinCluster,TinCluster);
	  }
	}
      }
    }
  }

  
  for( int jIndex = 0; jIndex < nTheta-1; jIndex++){
    profRDT_NRCL[jIndex]->Write();
    profRT_NRCL[jIndex]->Write();
    profDT_NRCL[jIndex]->Write();
  }

  for( int iIndex = 0; iIndex < nE-1; iIndex++){
    for( int jIndex = 0; jIndex < nTheta-1; jIndex++){
      profRT_MIP[iIndex][jIndex]->Write();
      profDT_MIP[iIndex][jIndex]->Write();
      profDTE[iIndex][jIndex]->Write();
      profRTE[iIndex][jIndex]->Write();
      profRDT[iIndex][jIndex]->Write();
      profRDE[iIndex][jIndex]->Write();

      profRDE_INV[iIndex][jIndex]->Write();
      profRDT_INV[iIndex][jIndex]->Write();
      profRT_INV[iIndex][jIndex]->Write();
      profDT_INV[iIndex][jIndex]->Write();

      profRDE_NR[iIndex][jIndex]->Write();
      profRDT_NR[iIndex][jIndex]->Write();
      profRT_NR[iIndex][jIndex]->Write();
      profDT_NR[iIndex][jIndex]->Write();
      
      profRDT_MIP[iIndex][jIndex]->Write();
      profRDE_MIP[iIndex][jIndex]->Write();
      profRDT_High[iIndex][jIndex]->Write();
      profRDE_High[iIndex][jIndex]->Write();
      profRDT_Low[iIndex][jIndex]->Write();
      profRDE_Low[iIndex][jIndex]->Write();
      profREAll[iIndex][jIndex]->Write();
      profRE_high[iIndex][jIndex]->Write();
      profRE_low[iIndex][jIndex]->Write();
      profDEAll[iIndex][jIndex]->Write();
      profDE_high[iIndex][jIndex]->Write();
      profDE_low[iIndex][jIndex]->Write();
      hisRPhi[iIndex][jIndex]->Write();
      hisRPhiChisqN[iIndex][jIndex]->Write();      
      hisRPhiE[iIndex][jIndex]->Write();
      hisPhiPhi[iIndex][jIndex]->Write();
      hisRT[iIndex][jIndex]->Write();
      hisDT[iIndex][jIndex]->Write();
      hisRE[iIndex][jIndex]->Write();
      hisDE[iIndex][jIndex]->Write();
      profTimeRD[iIndex][jIndex]->Write();
      profEET[iIndex][jIndex]->Write();
      profEET_D[iIndex][jIndex]->Write();
      profET_D[iIndex][jIndex]->Write();
      profET_D0[iIndex][jIndex]->Write();
      profET_D1[iIndex][jIndex]->Write();
      profET_D2[iIndex][jIndex]->Write();
      profET_D3[iIndex][jIndex]->Write();
      profRTAll[iIndex][jIndex]->Write();
      profRT_low[iIndex][jIndex]->Write();
      profRT_high[iIndex][jIndex]->Write();
      profDTAll[iIndex][jIndex]->Write();
      profDT_low[iIndex][jIndex]->Write();
      profDT_high[iIndex][jIndex]->Write();
      
    }
  }
  profET_D_All->Write();
  profET_D_All_min->Write();
  profET_D_All_max->Write();
  tfout->Close();
  return 0; 
}
