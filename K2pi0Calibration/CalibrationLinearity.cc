#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"

#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include "gamma/Gamma.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"

#include <cstdlib>
#include <cstdio>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <list>

#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"


int main( int argc ,char** argv){
  
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("neRMI");
  const int nCSI = 2716;  
  
  int    nCal[nCSI];
  double CalFactor[nCSI];
  double CalRMS[nCSI];
  double preCalFactor[nCSI];  
  
  for( int i = 0; i< nCSI; i++){
    preCalFactor[i] = 1;
  }


  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  TChain* ch = new TChain("trCalibration");
  ch->Add("CalibrationADV_*_10.root");
  //E14GNAnaDataContainer data;
  //data.setBranchAddress(ch);
  Int_t    GamClusNumber;
  Double_t GamClusDepE[12];
  Int_t    CutCondition;
  Double_t GammaEnergy[6];
  Double_t Ratio[6];
  Double_t SecondRatio[6];
  Double_t Corr[6];
  Int_t    FlagCalibrated[6];
  Int_t    CorrID[6];
  Double_t GammaSigma[6];
  Int_t    nCalibrated;
  Int_t    FlagKL_prefit;
  Double_t chisq[6];
  
  Int_t    LeadingChID[6];
  Double_t LeadingHeight[6];
  Double_t LeadingEnergy[6];
  ch->SetBranchAddress("CutCondition",&CutCondition);
  ch->SetBranchAddress("GamClusNumber",&GamClusNumber);
  ch->SetBranchAddress("GamClusDepE",GamClusDepE);//GamClusNumber;
  ch->SetBranchAddress("FlagKL_prefit",&FlagKL_prefit);
  ch->SetBranchAddress("GammaEnergy",GammaEnergy);
  ch->SetBranchAddress("Ratio",Ratio);
  ch->SetBranchAddress("chisq",chisq);
  ch->SetBranchAddress("SecondRatio",SecondRatio);
  ch->SetBranchAddress("Corr",Corr);
  ch->SetBranchAddress("FlagCalibrated", FlagCalibrated);
  ch->SetBranchAddress("CorrID",CorrID);
  ch->SetBranchAddress("GammaSigma",GammaSigma);
  ch->SetBranchAddress("nCalibrated",&nCalibrated);
  ch->SetBranchAddress("LeadingChID",LeadingChID);
  ch->SetBranchAddress("LeadingHeight",LeadingHeight);
  ch->SetBranchAddress("LeadingEnergy",LeadingEnergy);

  long Entries = ch->GetEntries();
  std::cout<< "Nentries:" << Entries  << std::endl;


  TFile* tfOut = new TFile("Output.root","recreate");
  TH1D* hisCalibrationFactor[nCSI];
  TH2D* hisCalibrationFactorEne[nCSI];
  TH2D* hisCalibrationFactorRatio[nCSI];
  TH2D* hisCalibrationFactorSecondRatio[nCSI];
  TH2D* hisCalibrationFactorSigma[nCSI];

  TH2D* hisCalibrationFactorHeight[nCSI];
  TH2D* hisCalibrationFactorHeight_Weighted[nCSI];
  TH2D* hisCalibrationFactorHeight_Weighted_INV[nCSI];

  for( int i =0; i < nCSI; i++){
    
    hisCalibrationFactor[i]
      = new TH1D(Form("hisCalibrationFactor_%d",i),
		 Form("CalibrationFactor:%d",i),
		 500,0,2);
    hisCalibrationFactorEne[i]
      = new TH2D(Form("hisCalibrationFactorEne_%d",i),
		 Form("hisCCalirationFactorEne_%d;Energy[MeV];CalibratioFactor",i),
		 20, 0 ,2000, 100, 0, 2);
    hisCalibrationFactorRatio[i]
      = new TH2D(Form("hisCalibrationRatio_%d",i),
		 Form("hisCalibrationRatio_%d;Ratio;CalibrationFactor", i),
		 20,0,1, 100,0,2);    
    hisCalibrationFactorSecondRatio[i]
      = new TH2D(Form("hisCalibrationSecondRatio_%d",i),
		 Form("hisCalibrationSecondRatio_%d;SecondRatio;CalibrationFactor",i),
		 20,0,1,100,0,2);
    hisCalibrationFactorSigma[i]
      = new TH2D(Form("hisCalibrationSigma_%d",i),
		 Form("hisCalibrationSigma_%d;Sigma of Gamma;CalibrationFactor",i),
		 50,0,50,100,0,2);
    hisCalibrationFactorHeight[i]
      = new TH2D(Form("hisCalibrationFactorHeight_%d",i),
		 Form("hisClaibrationFactorHeight_%d;Height_of_Channel;CalibrationFactor",i),
		 40,0,16000,100,0,2);	       
    hisCalibrationFactorHeight_Weighted[i]
      = new TH2D(Form("hisCalibrationFactorHeight_Weighted_%d",i),
		 Form("hisClaibrationFactorHeight_Weighted_%d;Height_of_Channel;Weighted_CalibrationFactor",i),
		 40,0,16000,100,0,2);
    hisCalibrationFactorHeight_Weighted_INV[i]
      = new TH2D(Form("hisCalibrationFactorHeight_Weighted_INV_%d",i),
		 Form("hisClaibrationFactorHeight_Weighted_INV_%d;Height_of_Channel;Weighted_&_Inversed_CalibrationFactor",i),
		 40,0,16000,100,0,2);	    
  }  

  ////////////////////////////////////////////////////////////////////////////
  // Loop
  ////////////////////////////////////////////////////////////////////////////
  //- Add Histogram for all Run;
  for( int ievet  = 0 ;ievet  < Entries ; ++ievet){
    if( ievet %1000  == 0 ){
      std::cout << ievet << "/" << Entries << std::endl;
    }
    
    ch->GetEntry(ievet);

    /*
      std::cout<< nCalibrated << std::endl;
      for( int i = 0; i< 6; i++ ) {
      std::cout << Corr[i] << " : " << GammaEnergy[i]  << " : " << CorrID[i] << std::endl;
      }
    */
    //- Fill Calibration Sample
    
    //std::vector<Klong> klVec;
    //data.setData(klVec);
    
    if((CutCondition & (1+8)) != 0){continue;}
    if( FlagKL_prefit != 0){continue;}
    if( nCalibrated  == 0 ){continue;} 
    /*
    double depE[6]={0};
    Int_t gIndex = 0;
    for( int i = 0; i< klVec[0].pi0().size(); i++){
      depE[gIndex] = klVec[0].pi0()[i].g1().edep();
      gIndex++;
      depE[gIndex] = klVec[0].pi0()[i].g2().edep();
      gIndex++;
    }
    */
    for( int i = 0; i< 6; i++){
      if( FlagCalibrated[i] == 0 ){
	double ratio = 1+((Corr[i]-1)*GamClusDepE[i]/LeadingEnergy[i]);
	hisCalibrationFactor[CorrID[i]]->Fill(Corr[i]);
	hisCalibrationFactorEne[CorrID[i]]->Fill(GammaEnergy[i],Corr[i]);
	hisCalibrationFactorRatio[CorrID[i]]->Fill(Ratio[i], Corr[i]);
	hisCalibrationFactorSecondRatio[CorrID[i]]->Fill(SecondRatio[i],Corr[i]);
	hisCalibrationFactorSigma[CorrID[i]]->Fill(GammaSigma[i],Corr[i]);;	
	hisCalibrationFactorHeight[LeadingChID[i]]->Fill(LeadingHeight[i],Corr[i]);
	//hisCalibrationFactorHeight_Weighted[LeadingChID[i]]->Fill(LeadingHeight[i],(1+((Corr[i]-1)*GammaEnergy[i]/LeadingEnergy[i])));
	//hisCalibrationFactorHeight_Weighted_INV[LeadingChID[i]]->Fill(LeadingHeight[i], 1./LeadingEnergy[i]*1./(1+((Corr[i]-1)*GammaEnergy[i]/LeadingEnergy[i])));
	hisCalibrationFactorHeight_Weighted[LeadingChID[i]]->Fill(LeadingHeight[i],ratio);
	hisCalibrationFactorHeight_Weighted_INV[LeadingChID[i]]->Fill(LeadingHeight[i],1/ratio);
      }
    }
  }

  TH1D* hisReNormCal = new TH1D("hisReNormCal",
				"Histogram for Renormalize Calibration Factor", 
				60,0.7,1.3);
  for( int i = 0; i < 2716; ++i){
    if(hisCalibrationFactor[i]->GetEntries() == 0){continue;}

    hisCalibrationFactor[i]->Write();
    hisCalibrationFactorEne[i]->Write();
    hisCalibrationFactorRatio[i]->Write();
    hisCalibrationFactorSecondRatio[i]->Write();
    hisCalibrationFactorSigma[i]->Write();
    hisCalibrationFactorHeight[i]->Write();
    hisCalibrationFactorHeight_Weighted[i]->Write();
    hisCalibrationFactorHeight_Weighted_INV[i]->Write();
    if( hisCalibrationFactor[i]->Integral() < 144 ){
      continue;
    }           
    hisCalibrationFactor[i]->Fit("gaus","Q","");
    TF1* calFunction= hisCalibrationFactor[i]->GetFunction("gaus");
    CalFactor[i]    = hisCalibrationFactor[i]->GetMean();
    CalRMS[i]       = hisCalibrationFactor[i]->GetRMS();
    //CalFactor[i] = calFunction->GetParameter(1);
    hisReNormCal->Fill(CalFactor[i]*preCalFactor[i]);
  }

  // Renomalize Calibration Factor, Mean -> 1
  
  tfOut->Close();
  
}
