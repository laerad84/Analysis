#include "TTree.h"
#include "TFile.h"

#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include <iostream>
#include "GeneralFunctions.h"
#include "ChisqCalc.h"
int main( int argc, char** argv){



  TFile* tfShape = new TFile("GammaTimeShapeRD.root");
  const int nAngle = 9;
  Double_t AngleBaseChisq[nAngle] = {0.373699,0.456174,0.509379,0.553928,0.604852,0.647827,0.688128,0.714016,0.747594};
  //Double_t AngleBaseChisq[nAngle]={0.403432,0.4644,0.513773,0.555798,0.606122,0.649179,0.688979,0.71543,0.742123};
//Double_t AngleBaseChisq[nAngle]={0.390709,0.45339,0.507781,0.551798,0.604505,0.649143,0.69142,0.724839,0.753404};


  TH2D* hisGammaTimeShape[nAngle];
  for( int i  =0; i< nAngle; i++){
    hisGammaTimeShape[i] = (TH2D*)tfShape->Get(Form("hisGammaTimeShape_%d",i));
  }
  TH2D* hisGammaTimeShapeBinID = new TH2D("hisGammaTimeShapeBinID","hisGammaTimeShapeBinID",
					  9,-112.5,112.5,
					  9,-112.5,112.5);

  TFile* tf = new TFile("kl_GammaTimeShape_DATA_NONTIMECAL.root");
  TTree* tr = (TTree*)tf->Get("GammaTimeShape");

  Double_t Radius;//Distance from Csi Center 
  Double_t ZVtx;//Distance from Csi surface
  Double_t theta;//InjectionAngle
  //Double_t cot;//cotangent of Injection Angle Z/R
  Double_t phi;//radial angle
  Double_t X;
  Double_t Y;
  Int_t    ClusterSize;
  Double_t E[120];//Energy;
  Double_t T[120];//Timing
  Double_t R[120];//Radial distance
  Double_t D[120];//
  Double_t FractionAngle[120];
  Int_t    CutCondition;
  Int_t    EventID;
  Double_t BaseTime;
  Double_t EGamma;
  Double_t ECenter;
  Double_t GammaChisq;
  tr->SetBranchAddress("EGamma",&EGamma);
  tr->SetBranchAddress("ECenter",&ECenter);
  tr->SetBranchAddress("EventID",&EventID);
  tr->SetBranchAddress("Radius",&Radius);
  tr->SetBranchAddress("X",&X);
  tr->SetBranchAddress("Y",&Y);
  tr->SetBranchAddress("ZVtx",&ZVtx);
  tr->SetBranchAddress("theta",&theta);
  tr->SetBranchAddress("phi",&phi);
  tr->SetBranchAddress("ClusterSize",&ClusterSize);
  tr->SetBranchAddress("BaseTime",&BaseTime);
  tr->SetBranchAddress("E",E);//ClusterSize
  tr->SetBranchAddress("T",T);//ClusterSize
  tr->SetBranchAddress("R",R);//ClusterSize
  tr->SetBranchAddress("D",D);//ClusterSize
  tr->SetBranchAddress("FractionAngle",FractionAngle);//ClusterSize
  tr->SetBranchAddress("CutCondition",&CutCondition);
  tr->SetBranchAddress("GammaChisq",&GammaChisq);

  TFile* tfOut = new TFile("GammaTimeChisq.root","recreate");
  TTree* trOut = new TTree("GammaTimeChisq","");
  Int_t initialAngleID;
  Int_t MinimumAngleID;
  Double_t ChisqTime[nAngle];
  Double_t NDF[nAngle];
  Double_t RegionChisqTime[nAngle];
  Double_t RegionNDF[nAngle];
  Double_t NormalizedChisq[nAngle];
  Double_t MinTimeOffset[nAngle];
  Int_t TestAngle[nAngle]={0,1,2,3,4,5,6,7,8};
  Int_t RegionMinimumAngleID;
  Int_t NormalizedMinimumID;
  trOut->Branch("initialAngleID",&initialAngleID,"initialAngleID/I");
  trOut->Branch("ChisqTime",ChisqTime,"ChisqTime[9]/D");
  trOut->Branch("NDF",NDF,"NDF[9]/D");
  trOut->Branch("RegionChisqTime",RegionChisqTime,"RegionChisqTime[9]/D");
  trOut->Branch("RegionNDF",RegionNDF,"RegionNDF[9]/D");
  trOut->Branch("NormalizedChisq",NormalizedChisq,"NormalizedChisq[9]/D");
  trOut->Branch("MinTimeOffset",MinTimeOffset,"MinTimeOffset[9]/D");
  trOut->Branch("TestAngle",TestAngle,"TestAngle[9]/I");
  trOut->Branch("GammaChisq",&GammaChisq,"GammaChisq/D");
  trOut->Branch("ClusterSize",&ClusterSize,"ClusterSize/I");
  trOut->Branch("EGamma",&EGamma,"EGamma/D");
  trOut->Branch("MinimumAngleID",&MinimumAngleID,"MinimumAngleID/I");
  trOut->Branch("RegionMinimumAngleID",&RegionMinimumAngleID,"RegionMinimumAngleID/I");
  trOut->Branch("NormalizedMinimumID",&NormalizedMinimumID,"NormalizedMinimumID/I");
  Double_t *AngleDist = GenLogArray(nAngle+1,1,40);
  const int      nDiv      = 4;
  const int      nDivTime  = 200;

  TH1D* hisAngle = new TH1D("hisAngle","hisAngle",9,AngleDist);
  for( int ievent = 0; ievent< tr->GetEntries(); ievent++){
    tr->GetEntry(ievent);    
    if( (ievent %10000) == 0 ){ std::cout<< ievent << "/" << tr->GetEntries() << std::endl;}
    for( int i =0; i< nAngle; i++){
      ChisqTime[i] = 0;
      TestAngle[i] = 0;
      NDF[i] = 0;
      RegionChisqTime[i] = 0;
      RegionNDF[i] = 0;
    }

    Int_t AngleIndex = hisAngle->Fill(ZVtx/Radius)-1;
    if( AngleIndex >= nAngle ){ continue; }
    if( AngleIndex <0 ){ continue; }
    if( CutCondition != 0 ){ continue; }
    if( EGamma < 100 ){ continue; }
    initialAngleID =AngleIndex;
    
    MinimumAngleID = -1;
    double MinimumAngleChisq = 100000000;
    RegionMinimumAngleID = -1;
    double RegionMinimumAngleChisq = 100000000;
    NormalizedMinimumID = -1;
    double NormMinChisq = 100000000;
    for( int iangle = 0; iangle <nAngle; iangle++){
      for( int i = 0; i< ClusterSize; i++){
	//if( E[i] < 10 ){ continue;}
	int nBinCount = hisGammaTimeShapeBinID->Fill(R[i],D[i]);
	Double_t BinError = hisGammaTimeShape[iangle]->GetBinError(nBinCount);
	if( BinError == 0 ){ continue; }
	Double_t BinCenter = hisGammaTimeShape[iangle]->GetBinContent(nBinCount);
	//double fracChisq = TMath::Power(BinCenter - (T[i]-BaseTime),2)/GetTimeResSq(E[i]);
	//double sigma = BinError*BinError+GetTimeResSq(E[i]);
	double sigma =GetTimeResSq(E[i]);
	double fracChisq = TMath::Power(BinCenter - (T[i]-BaseTime),2)/sigma;
	ChisqTime[iangle]+=fracChisq;
	NDF[iangle]++;	
	/*
	if( TMath::Abs( D[i] ) < 12.5*5 ){
	  if( R[i] > -12.5* 5 && R[i] < 12.5 *7 ){	   
	    RegionChisqTime[iangle]+=fracChisq;
	    RegionNDF[iangle]++;
	  }
	}
	*/
      }
      


      double testRegionChisq=0;
      double testRegionChisqOld=1000000;
      double testRegionNDF=0;      
      for( int itime = 0; itime< 21; itime++){
	double offsetTime = BaseTime-1+itime*0.1;	
	RegionTimeChisq( hisGammaTimeShape[iangle], hisGammaTimeShapeBinID, ClusterSize, R,D,E,T,offsetTime,testRegionChisq,testRegionNDF);
	if( testRegionChisqOld > testRegionChisq ){
	  testRegionChisqOld = testRegionChisq;
	  RegionChisqTime[iangle] = testRegionChisq;
	  RegionNDF[iangle]       = testRegionNDF;
	  MinTimeOffset[iangle]   = offsetTime;
	  //std::cout<< RegionChisqTime[iangle] << std::endl;
	}
      }

      TestAngle[iangle] = iangle;
      ChisqTime[iangle] = ChisqTime[iangle]/(NDF[iangle]-1);
      
      //RegionChisqTime[iangle] = RegionChisqTime[iangle]/(RegionNDF[iangle]-1);
      NormalizedChisq[iangle] = ChisqTime[iangle]/AngleBaseChisq[iangle];
      if( ChisqTime[iangle] < MinimumAngleChisq ){
	MinimumAngleChisq = ChisqTime[iangle];
	MinimumAngleID = iangle;
      }
      if( RegionChisqTime[iangle] < RegionMinimumAngleChisq && iangle > 0){
	RegionMinimumAngleChisq = RegionChisqTime[iangle];
	RegionMinimumAngleID = iangle;
      }
      if( NormMinChisq > NormalizedChisq[iangle] ){
	NormMinChisq = NormalizedChisq[iangle];
	NormalizedMinimumID = iangle;
      }
		      


    }
    trOut->Fill();
  }
  trOut->Write();
  tfOut->Close();
  return 0; 
}
