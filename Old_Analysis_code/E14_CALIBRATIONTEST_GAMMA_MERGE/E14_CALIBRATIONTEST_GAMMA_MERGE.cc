#include <iostream>
#include <string>
#include <cstring>

//#include "gnana/E14GNAnaFunction.h"
//#include "gnana/E14GNAnaDataContainer.h"
//#include "rec2g/Rec2g.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TPostScript.h"

//#include "klong/Klong.h"
//#include "pi0/Pi0.h"

#include "CLHEP/Vector/ThreeVector.h"

//#include "CalibrationTree.h"

//void Merge(int iterationNumber){

int 
main( int argc, char** argv ){
  int iterationNumber  = 0;

  std::string Filename = "SimCalibration_Data/Calibration_%04d_%d.root";
  TChain* ch = new TChain("trCalibration");
  for( int i = 0; i < 200; ++i){
    ch->Add(Form(Filename.c_str(),i*10,iterationNumber));
  }

  Int_t KlongFit;
  Double_t Corr[6];
  Int_t CorrID[6];
  Double_t CorrE[6];
  Double_t chisq[6];
  Int_t FlagCalibrated[6];
  Int_t nCalibrated;
  Double_t SimGammaEne[6];
  Double_t RecGammaEne[12];
  Double_t CalGammaEne[6];
  Double_t CstGammaEne[6];

  Double_t KlongPos[2][3];
  Double_t KlongMom[2][3];
  Double_t KlongChisqZ[2];
  Double_t KlongMass[2];
  Double_t SimKLMom[3];
  Double_t Pi0Mass[6];
  Double_t Pi0RecZ[6];

  Int_t CutCondition;

  ch->SetBranchAddress("Pi0Mass",Pi0Mass);
  ch->SetBranchAddress("Pi0RecZ",Pi0RecZ);

  ch->SetBranchAddress("KlongMass",KlongMass);
  ch->SetBranchAddress("KlongMom",KlongMom);
  ch->SetBranchAddress("KlongPos",KlongPos);
  ch->SetBranchAddress("KlongChisqZ",KlongChisqZ);
  ch->SetBranchAddress("SimKLMom",SimKLMom);
  ch->SetBranchAddress("SimGammaEneMatched",SimGammaEne);
  //ch->SetBranchAddress("GammaEnergy",CalGammaEne);
  ch->SetBranchAddress("GammaE"     ,RecGammaEne);
  ch->SetBranchAddress("GamClusDepE",CstGammaEne);

  ch->SetBranchAddress("CutCondition",&CutCondition);
  ch->SetBranchAddress("KlongFit", &KlongFit);
  ch->SetBranchAddress("Corr",Corr);
  ch->SetBranchAddress("CorrID",CorrID);  
  ch->SetBranchAddress("CorrE",CalGammaEne);
  ch->SetBranchAddress("chisq",chisq);
  ch->SetBranchAddress("FlagCalibrated",FlagCalibrated);
  ch->SetBranchAddress("nCalibrated",&nCalibrated);

  TFile* tf = new TFile("SimCalibration_Data/Merge_All.root","recreate");
  // TO Draw //
  // 
  // Compare Calibration Constant vs KLMass
  // Compare Calibration Constant vs Pi0Mass 
  // Compare Calibration Constant vs Gamma Energy
  // Compare Calibration Constant vs Lower GammaEnergy
  // Compare Calibration Constant vs Delta Of Reconstructed Pi0 Position.
  
  // Compare Pi0Mass vs GammaEnergy
  // Compare Pi0Mass vs LowerGammaE
  // Compare Pi0Mass vs Delta of Reconstructed Pi0 Position - Reconstructed KL Position .

  TH2D* his_Cal_KLMass = new TH2D("his_Cal_KLMass",
				  "KLMass:Calibration Sample;KLMass;Calibration Sample",
				  40,480,520,60,0.7,1.3);
  TH2D* his_Cal_Pi0Mass = new TH2D("his_Cal_Pi0Mass",
				   "Pi0Mass:Calibration Sample;Pi0Mass;Calibration Sample",
				   40,125,145,60,0.7,1.3);
  TH2D* his_Cal_GammaE = new TH2D("his_Cal_GammaE",
				  "GammaE:Calibration Sample;GammaE;Calibration Sample",
				  40,0,2000,60,0.7,1.3);
  TH2D* his_Cal_LGammaE      = new TH2D("his_Cal_LGammaE",
					"LGammaE:Calibration Cample;GammaE;Calibration Sample",
					40,0,2000,60,0.7,1.3);
  TH2D* his_Cal_DeltaRec     = new TH2D("his_Cal_DeltaRec",
					"RecPosition:Calibration Sample;DeltaRecPosition;Calibration Sample",
					40,-200,200,60,0.7,1.3);
  TH2D* his_Pi0Mass_GammaE   = new TH2D("his_Pi0Mass_GammaE",
					"Gamma Energy:Pi0Mass;Gamma Energy;Pi0Mass",
					40,0,2000,40,125,145);
  TH2D* his_Pi0Mass_LGammaE  = new TH2D("his_Pi0Mass_LGammaE",
					"LGammaE:Pi0Mass;Low Gamma Energy;Pi0Mass",
					40,0,2000,40,125,145);
  TH2D* his_Pi0Mass_DeltaRec = new TH2D("his_Pi0Mass_DeltaRec",
					"RecPosition:Pi0Mass;DeltaRecPosition;Pi0Mass",
					40,-200,200,50,120,150);
  TH2D* his_GammaE_DeltaE    = new TH2D("his_GammaE_DeltaE","",
					40,0,2000,200,0,2);
  TH2D* his_HGammaE_LGammaE_KLChisqZ = new TH2D("his_HGammaE_LGammaE_KLChisqZ","",
						40,0,2000,40,0,2000);
  TH1D* his_KLChisqZ_HGLG[40][40];
  for( int i = 0; i< 40; i++){
    for( int j  =0; j< 40; j++){
      his_KLChisqZ_HGLG[i][j] = new TH1D(Form("his_ChisqHGLG_%d_%d",i,j),"",100,-1,1);
    }
  }



  const int nCH = 2716;
  TH1D* his_Corr_All[nCH];
  TH1D* his_Corr_Mis[nCH];
  TH1D* his_Corr_Fit[nCH];
  TH1D* his_Gamma_Ratio_Rec[nCH];
  TH1D* his_Gamma_Ratio_Cal[nCH];
  TH2D* his_Mass_Corr[nCH];
  TH2D* his_Energy_Corr[nCH];
  TH2D* his_DeltaCos_Pi0_All    = new TH2D("his_DeltaCos_Pi0_All","",
					   60,-6,6,80,-0.01,0.01);
  TH2D* his_GammaE_Pi0_All      = new TH2D("his_GammaE_Pi0_All","",
					   60,-6,6,80,-0.2,0.2);
  TH2D* his_GammaE_Min_All      = new TH2D("his_GammaE_Min_All","",
					   60,-6,6,80,-0.2,0.2);
  TH2D* his_GammaE_Max_All      = new TH2D("his_GammaE_Min_Max","",
					   60,-6,6,80,-0.2,0.2);

  TH2D* his_deltaPz_recPz       = new TH2D("his_deltaPz_recPz","",
					   75,0,7500,100,-500,500);
  TH2D* his_KlMass_Corr_All     = new TH2D("his_KlMass_Corr_All","",
					   50,0.5,1.5,100,475,520);
  TH2D* his_Pi0Mass_Corr_All    = new TH2D("his_Pi0Mass_Corr_All","",
					   50,0.5,1.5,100,125,145);
  TH2D* his_ChisqZ_Corr_All     = new TH2D("his_ChisqZ_Corr_All","",
					   50,0.5,1.5,100,0,10);
  TH2D* his_GammaE_Corr_All     = new TH2D("his_GammaE_Corr_All","",
					   50,0.5,1.5,100,0,4000);
  TH2D* his_CalChisq_Corr_All   = new TH2D("his_CalChisq_Corr_All","",
					   50,0.5,1.5,100,0,20);
  TH2D* his_KlMass_Corr_All_CUT = new TH2D("his_KlMass_Corr_All_CUT","",
					   50,0.5,1.5,100,475,520);
  TH2D* his_Pi0Mass_Corr_All_CUT   = new TH2D("his_Pi0Mass_Corr_All_CUT","",
					      50,0.5,1.5,40,125,145);
  TH2D* his_Pi0Mass_GammaE_All_CUT = new TH2D("his_Pi0Mass_GammaE_All_CUT","",
					      50,0,5000,40,125,145);

  for( int i = 0; i< 2716; ++i ){    
    his_Corr_All[i] = new TH1D(Form("hisCalib_All_%04d",i),"", 200,0,2);
    his_Corr_Fit[i] = new TH1D(Form("hisCalib_Fit_%04d",i),"", 200,0,2);
    his_Corr_Mis[i] = new TH1D(Form("hisCalib_Mis_%04d",i),"", 200,0,2);
    
    his_Gamma_Ratio_Rec[i] = new TH1D(Form("his_Gamma_Ratio_Rec_%04d",i),"",200,0,2);
    his_Gamma_Ratio_Cal[i] = new TH1D(Form("his_Gamma_Ratio_Cal_%04d",i),"",200,0,2);
    his_Mass_Corr[i]       = new TH2D(Form("his_Mass_Corr_%04d",i),"",
				      50,0.5,1.5,40,450,550);
    his_Energy_Corr[i]     = new TH2D(Form("his_Energy_Corr_%04d",i),"",
				      50,0,5000,50,0.5,1.5);
  }
  long nevent = ch->GetEntries();
  
  for( int ievt = 0 ; ievt < nevent ; ++ievt){
    if( ievt%1000==0){
      std::cout<< ievt << " / " << nevent << std::endl;
    }
    ch->GetEntry(ievt);    

    if((CutCondition&(17)) == 0 ){
      for( int ipi0 = 0; ipi0 < 3; ++ipi0){
	
	Double_t maxGammaDelta=(RecGammaEne[ipi0*2]-SimGammaEne[ipi0*2])/RecGammaEne[ipi0*2];
	Double_t minGammaDelta=(RecGammaEne[ipi0*2+1]-SimGammaEne[ipi0*2+1])/RecGammaEne[ipi0*2+1];
	
	his_GammaE_Pi0_All->Fill( (Pi0Mass[ipi0]-134.9766),minGammaDelta);
	his_GammaE_Pi0_All->Fill( (Pi0Mass[ipi0]-134.9766),maxGammaDelta);
	his_GammaE_Max_All->Fill( (Pi0Mass[ipi0]-134.9766),maxGammaDelta);
	his_GammaE_Min_All->Fill( (Pi0Mass[ipi0]-134.9766),minGammaDelta);
	
	
	double costhetaSim = 1-134.9766*134.9766/(SimGammaEne[ipi0*2]*SimGammaEne[ipi0*2+1]*2);
	double costhetaRec = 1-134.9766*134.9766/(RecGammaEne[ipi0*2]*RecGammaEne[ipi0*2+1]*2);
	his_DeltaCos_Pi0_All->Fill( (Pi0Mass[ipi0]-134.9766),costhetaRec-costhetaSim);
      }
    }
    
    if( (CutCondition &( 17 ) ) == 0 ){
      if( nCalibrated != 0 ){
	for( int ipi0 = 0; ipi0 <3; ipi0++){
	  int iIndex = (int)(RecGammaEne[ipi0*2]/50);
	  int jIndex = (int)(RecGammaEne[ipi0*2+1]/50);
	  his_KLChisqZ_HGLG[iIndex][jIndex]->Fill((KlongPos[0][2]-Pi0RecZ[ipi0])/TMath::Abs(KlongPos[0][2]-6148));
	}


	for( int igamma = 0; igamma < 6; igamma++){
	  if(FlagCalibrated[igamma] != 0 ){continue;}
	  his_GammaE_DeltaE->Fill(RecGammaEne[igamma],SimGammaEne[igamma]/RecGammaEne[igamma]);
	  his_Cal_KLMass->Fill(KlongMass[0],Corr[igamma]);
	  his_Cal_Pi0Mass->Fill(Pi0Mass[igamma/2],Corr[igamma]);
	  his_Cal_GammaE->Fill(RecGammaEne[igamma],Corr[igamma]);
	  if( igamma%2 !=0 ){// Low energy gamma
	    his_Cal_LGammaE->Fill(RecGammaEne[igamma],Corr[igamma]);
	    his_Pi0Mass_LGammaE->Fill(RecGammaEne[igamma],Pi0Mass[igamma/2]);
	  }
	  his_Cal_DeltaRec->Fill(KlongPos[0][2]-Pi0RecZ[igamma/2],Corr[igamma]);
	  his_Pi0Mass_GammaE->Fill(RecGammaEne[igamma],Pi0Mass[igamma/2]);
	}
      }
    }
    
    
    if( nCalibrated != 0 ){
      if( FlagCalibrated[0] == 0){
	his_Pi0Mass_Corr_All_CUT->Fill(Corr[0], Pi0Mass[0]);
      }
      
      for( int i = 0; i< 6; ++i ){
	if( FlagCalibrated[i] == 0 ){
	  his_KlMass_Corr_All->Fill(Corr[i],KlongMass[0]);
	  int FlagPi0Mass=0;
	  for( int j = 0 ; j< 3 ; ++j ){
	    if( TMath::Abs(Pi0Mass[0]- 134.9766)> 3 ){
	    FlagPi0Mass|=1;
	    }
	  }
	  if( FlagPi0Mass!=1 && TMath::Abs(KlongMass[0]-497.648) < 5 && CorrID[i]==1600){
	    his_KlMass_Corr_All_CUT->Fill(Corr[i],KlongMass[0]);	    
	  }
	  if( i%2 == 1){
	    his_Pi0Mass_GammaE_All_CUT->Fill( RecGammaEne[i],Pi0Mass[i/2]);
	  }
	  //if( (CutCondition & (1+2+4+8))==0){
	  his_Corr_All[CorrID[i]]->Fill(Corr[i]);
	  if(KlongFit==1){
	    his_Corr_Fit[CorrID[i]]->Fill(Corr[i]);
	  }else{
	    his_Corr_Mis[CorrID[i]]->Fill(Corr[i]);
	  }	  
	  
	  his_Gamma_Ratio_Rec[CorrID[i]]->Fill((RecGammaEne[i])/CstGammaEne[i]);
	  his_Gamma_Ratio_Cal[CorrID[i]]->Fill((CalGammaEne[i])/CstGammaEne[i]);
	  his_Mass_Corr[CorrID[i]]->Fill(Corr[i],KlongMass[0]);
	  his_Energy_Corr[CorrID[i]]->Fill(RecGammaEne[i],Corr[i]);
	  his_deltaPz_recPz->Fill(KlongPos[0][2],SimKLMom[2]-KlongMom[0][2]);
	  his_ChisqZ_Corr_All->Fill(Corr[i],KlongChisqZ[0]);
	  his_GammaE_Corr_All->Fill(Corr[i],RecGammaEne[i]);	    
	  his_CalChisq_Corr_All->Fill(Corr[i],chisq[i]);
	  for( int j = 0; j< 3; ++j){
	    his_Pi0Mass_Corr_All->Fill(Corr[i],Pi0Mass[0]);
	  }
	  //}
	}
	
      }     
    }
  }
  for( int i = 0; i<40; i++){
    for( int j = 0; j< 40; j++){
      his_HGammaE_LGammaE_KLChisqZ->SetBinContent(i+1,j+1,his_KLChisqZ_HGLG[i][j]->GetRMS());
    }
  }

  his_HGammaE_LGammaE_KLChisqZ->Write();


  for( int i = 0; i< 2716; ++i){
    std::cout<<i<< ":" << his_Corr_All[i]->GetEntries() << std::endl;
    his_Corr_All[i]->Write();
    his_Corr_Fit[i]->Write();
    his_Corr_Mis[i]->Write();
    his_Gamma_Ratio_Rec[i]->Write();
    his_Gamma_Ratio_Cal[i]->Write();
    his_Mass_Corr[i]->Write();
    his_Energy_Corr[i]->Write();    
  }
  his_DeltaCos_Pi0_All->Write();
  his_GammaE_Pi0_All->Write();
  his_GammaE_Max_All->Write();
  his_GammaE_Min_All->Write();
  his_deltaPz_recPz->Write();
  his_KlMass_Corr_All->Write();
  his_KlMass_Corr_All_CUT->Write();
  his_Pi0Mass_Corr_All->Write();
  his_ChisqZ_Corr_All->Write();
  his_GammaE_Corr_All->Write();
  his_CalChisq_Corr_All->Write();
  his_Pi0Mass_Corr_All_CUT->Write();
  his_Pi0Mass_GammaE_All_CUT->Write();

  his_Cal_KLMass->Write();
  his_Cal_Pi0Mass->Write();
  his_Cal_GammaE->Write();
  his_Cal_LGammaE->Write();
  his_Cal_DeltaRec->Write();
  his_Pi0Mass_GammaE->Write();
  his_Pi0Mass_LGammaE->Write();
  his_GammaE_DeltaE->Write();

  std::cout <<ch->GetEntries() << std::endl;
  std::cout<< __LINE__<< std::endl;
  tf->Close();
  
  return 0;
}

  
