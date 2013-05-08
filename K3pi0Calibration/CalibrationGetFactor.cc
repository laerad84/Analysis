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

#include <cstdlib>
#include <cstdio>
#include <list>

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
  
  Int_t CalibrationNumber  = -1; //iterationNumber 
  std::string runListFile;
  std::string path;
  if( argc == 3){
    CalibrationNumber = atoi(argv[1]);
    runListFile       = argv[2];
  }else if( argc == 4 ){
    CalibrationNumber = atoi(argv[1]);
    runListFile       = argv[2];
    path = argv[3];
  }else{
    std::cerr << "<<<>>> Argument Error <<<>>> " << "\n"
	      << "Usage:" << argv[0] 
	      << "[CalibrationNumber] [Filename of List] "
	      << std::endl;
    return -1;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // Set initial CalibrationFactor and Files 
  //////////////////////////////////////////////////////////////////////////////
  
  std::vector<int> VecRunNum;
  std::ifstream ifsList(runListFile.c_str());
  std::cout<< runListFile << std::endl;
  if( !ifsList.is_open() ){
    return -1; 
  }
  while( !ifsList.eof() ){
    Int_t runNumber;
    if( ifsList >>  runNumber){
      VecRunNum.push_back(runNumber );
      //if(runNumber ==4200){break;}
    }
  }
  
  int    nCal[nCSI];
  double CalFactor[nCSI];
  double CalRMS[nCSI];
  double preCalFactor[nCSI];  
  
  for( int i = 0; i< nCSI; i++){
    preCalFactor[i] = 1;
  }
  
  std::string CalibrationDataPath=std::getenv("ROOTFILE_3PI0CALIBRATION");
  if( argc ==4){
    CalibrationDataPath += "/";
    CalibrationDataPath += path.c_str();
  }

  
  std::string InputCalData      = CalibrationDataPath.c_str();
  std::string InputRootFile     = CalibrationDataPath.c_str();
  std::string OutputCalData     = CalibrationDataPath.c_str();
  std::string OutputCalStatData = CalibrationDataPath.c_str();
  std::string OutputCalHist     = CalibrationDataPath.c_str();

  InputCalData      += "/CalibrationFactorADV_%d.dat";
  InputRootFile     += "/CalibrationADV_%d_%d.root";
  OutputCalData     += "/CalibrationFactorADV_%d.dat";
  OutputCalStatData += "/CalibrationStaticsADV_%d.dat";
  OutputCalHist     += "/CalhistListADV_%s_%d.root";

  std::cout<< InputCalData      << std::endl;
  std::cout<< OutputCalData     << std::endl;
  std::cout<< OutputCalStatData << std::endl;
  std::cout<< OutputCalHist     << std::endl;
  
  //CalibrationDataPath+="/CalibrationFactorADV_%d.dat";
  if( CalibrationNumber!=0 ){
    std::ifstream ifs(Form(InputCalData.c_str(),
			   CalibrationNumber));
    if(!ifs.is_open()){ return -1;}
    Int_t id;
    Double_t precal;
    while( !ifs.eof() ){
      if( ifs >> id >> precal){
	preCalFactor[id] = precal;
      }
    }
  }
  
  int initialRunNumber;
  int finalRunNumber;
  Int_t nextCalNum = CalibrationNumber +1;
  std::ofstream ofs(Form(OutputCalData.c_str(),
			 nextCalNum));
  std::ofstream ofs1(Form(OutputCalStatData.c_str(),
			  nextCalNum));
  std::string listFilename=runListFile.substr(runListFile.find_last_of('/')+1);
  std::cout<< listFilename << std::endl;
  TFile* tfOut 
    = new TFile(Form(OutputCalHist.c_str(),
		     listFilename.substr(0,-4).c_str(),
		     CalibrationNumber),
		"RECREATE");
  
  ////////////////////////////////////////////////////////////////////////////  
  // Calibration Factor Histogram
  ////////////////////////////////////////////////////////////////////////////  
  TH1D* hisKLMassRaw    = new TH1D("his_KL_mass_Raw","kl_mass_raw",400,400,600);
  TH1D* hisKLMassBefore = new TH1D("his_KL_mass_Use","kl_mass_Use",400,400,600);
  TH1D* hisKLMassAfter  = new TH1D("his_KL_mass_Calibrated" ,"kl_mass_Calibrated" ,400,400,600);
  TH1D* hisPi0MassBefore= new TH1D("his_Pi0_mass_Before","pi0_mass_before",400,100,200);
  TH1D* hisPi0MassAfter = new TH1D("his_Pi0_mass_After" ,"pi0_mass_after" ,400,100,200);
  
  TH1D* hisChisq        = new TH1D("his_Chisq","chisq",100,-1,10);
  TH1D* hisChisqDof     = new TH1D("his_Chisq_dof","chist_dof",100,-1,100);
  TH1D* hisChisqZ       = new TH1D("his_ChisqZ","KL_ChisqZ",80,0,40);
  TH1D* hisChisqZ_2nd   = new TH1D("his_ChisqZ_2nd","KL_ChisqZ_2nd",80,0,40);
  TH1D* hisChisqZ_Cal   = new TH1D("his_ChisqZ_Cal","KL_ChisqZ_Cal",80,0,40);
  TH1D* hisChisqZ_Cal_2nd = new TH1D("his_ChisqZ_Cal_2nd","KL_ChisqZ_Cal_2nd",80,0,40);

  TH1D* hisKlongPt      = new TH1D("his_KlongPt","KlongPt",100,0,200);
  TH1D* hisKlongPz      = new TH1D("his_KlongPz","KlongPz",100,0,5000);
  TH1D* hisKlongPt_Cal  = new TH1D("his_KlongPt_Cal","KlongPt_Cal",100,0,200);
  TH1D* hisKlongPz_Cal  = new TH1D("his_KlongPz_Cal","KlongPz_Cal",100,0,5000);  

  TH1D* hisKlongVz      = new TH1D("his_KlongVz","KlongVz",50,1000,6000);
  TH1D* hisKlongVz_Cal  = new TH1D("his_KlongVz_Cal","KlongVz_Cal",50,1000,6000);

  TH1D* hisCalSigmaTotal     = new TH1D("his_SigmaCalSample","Sigma_Of_Calibration_Sample",
					100,0,0.1);  
  TH2D* hisCalFactorEneTotal = new TH2D("hisCalFactorEneTotal",
					"hisCalFactorEneTotal",
					20,0,2000,100,0,2);
  TH2D* hisCalFactorRatioTotal = new TH2D("hisCalFactorRatioTotal",
					  "hisCalFactorRatioTotal",
					  20,0,1,100,0,2);
  TH2D* hisCalFactorSecondRatioTotal = new TH2D("hisCalFactorSecondRatioTotal",
						"hisCalFactorSecondRatioTotal",
						20,0,1,100,0,2);
  TH2D* hisCalFactorSigmaTotal = new TH2D("hisCalFactorSigmaTotal",
					  "hisCalFactorSigmaTotal",
					  50,0,50,100,0,2);
  
  
  TH1D* hisCalibrationFactor[nCSI];
  TH2D* hisCalibrationFactorEne[nCSI];
  TH2D* hisCalibrationFactorRatio[nCSI];
  TH2D* hisCalibrationFactorSecondRatio[nCSI];
  TH2D* hisCalibrationFactorSigma[nCSI];

  TH2D* hisCalibrationFactorHeight[nCSI];
  TH2D* hisCalibrationFactorHeight_Weighted[nCSI];

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
	       
  }  
  



  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////
  TChain* ch = new TChain("trCalibration");
  std::vector<int>::iterator iRun;

  for( iRun  = VecRunNum.begin();
       iRun != VecRunNum.end();
       ++iRun){
    std::string Filename = Form(InputRootFile.c_str(), *iRun, CalibrationNumber);
    int rst = ch->Add(Filename.c_str());
    if(rst ==0 ) {std::cout << Filename << std::endl; }
  }
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

  Int_t    CutCondition;  
  Int_t    KlongNumber;
  Double_t KlongId[2];
  Double_t KlongMass[2];
  Double_t KlongE[2];
  Double_t KlongPos[2][3];
  Double_t KlongMom[2][3];
  Double_t KlongPt[2];
  Double_t KlongDeltaZ[2];
  Double_t KlongChisqZ[2];

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


  ch->SetBranchAddress("CutCondition",&CutCondition);
  ch->SetBranchAddress("KlongNumber",&KlongNumber);
  ch->SetBranchAddress("KlongId",KlongId);//KlongNumber
  ch->SetBranchAddress("KlongMass",KlongMass);//KlongNumber
  ch->SetBranchAddress("KlongPos",KlongPos);//KlongNumber
  ch->SetBranchAddress("KlongPt",KlongPt);//KlongNumber
  ch->SetBranchAddress("KlongMom",KlongMom);//KlongNumber
  ch->SetBranchAddress("KlongDeltaZ",KlongDeltaZ);//KlongNumber
  ch->SetBranchAddress("KlongChisqZ",KlongChisqZ);//KlongNumber
  ch->SetBranchAddress("KlongE",KlongE);//KlongNumber  


  long Entries = ch->GetEntries();
  std::cout<< "Nentries:" << Entries  << std::endl;

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
    hisKLMassRaw->Fill(KlongMass[0]);    
    hisChisqZ->Fill(KlongChisqZ[0]);
    if( KlongNumber >1){
      hisChisqZ_2nd->Fill(KlongChisqZ[1]);
    }
    
    if((CutCondition & (1+8)) != 0){continue;}
    hisKlongPt->Fill(KlongPt[0]);
    hisKlongPz->Fill(KlongMom[0][2]);
    hisKlongVz->Fill(KlongPos[0][2]);
    
    if( FlagKL_prefit != 0){continue;}
    hisKLMassBefore->Fill(KlongMass[0]);

    if( nCalibrated  == 0 ){continue;} 
    hisKLMassAfter->Fill(KlongMass[0]);
    hisChisqZ_Cal->Fill(KlongChisqZ[0]);
    if( KlongNumber >1 ){
      hisChisqZ_Cal_2nd->Fill(KlongChisqZ[1]);
    }
    hisKlongPt_Cal->Fill(KlongPt[0]);
    hisKlongPz_Cal->Fill(KlongMom[0][2]);
    hisKlongVz_Cal->Fill(KlongPos[0][2]);

    for( int i = 0; i< 6; i++){
      if( FlagCalibrated[i] == 0 ){
	hisCalibrationFactor[CorrID[i]]->Fill(Corr[i]);
	hisCalibrationFactorEne[CorrID[i]]->Fill(GammaEnergy[i],Corr[i]);
	hisCalibrationFactorRatio[CorrID[i]]->Fill(Ratio[i], Corr[i]);
	hisCalibrationFactorSecondRatio[CorrID[i]]->Fill(SecondRatio[i],Corr[i]);
	hisCalibrationFactorSigma[CorrID[i]]->Fill(GammaSigma[i],Corr[i]);;	
	hisCalibrationFactorHeight[LeadingChID[i]]->Fill(LeadingHeight[i],Corr[i]);
	hisCalibrationFactorHeight_Weighted[LeadingChID[i]]->Fill(LeadingHeight[i], 1./(1+((Corr[i]-1)*GammaEnergy[i]/LeadingEnergy[i])));
	hisChisq->Fill(chisq[i]);
	hisChisqDof->Fill(chisq[i]/4);
      }
    }
  }
  
  hisKLMassRaw->Write();
  hisKLMassBefore->Write();
  hisKLMassAfter->Write();
  hisPi0MassBefore->Write();
  hisPi0MassAfter->Write();
  hisChisq->Write();
  hisChisqDof->Write();
  hisChisqZ->Write();
  hisChisqZ_2nd->Write();
  hisChisqZ_Cal->Write();
  hisChisqZ_Cal_2nd->Write();
  hisKlongPt->Write();
  hisKlongPt_Cal->Write();
  hisKlongPz->Write();
  hisKlongPz_Cal->Write();
  hisKlongVz->Write();
  hisKlongVz_Cal->Write();

  TH1D* hisReNormCal = new TH1D("hisReNormCal",
				"Histogram for Renormalize Calibration Factor", 
				60,0.7,1.3);
  for( int i = 0; i < 2716; ++i){
    if(hisCalibrationFactor[i]->GetEntries() == 0){continue;}
    hisCalSigmaTotal->Fill(hisCalibrationFactor[i]->GetRMS());
    hisCalFactorRatioTotal->Add(hisCalibrationFactorRatio[i]);
    hisCalFactorSecondRatioTotal->Add(hisCalibrationFactorSecondRatio[i]);
    hisCalFactorEneTotal->Add(hisCalibrationFactorEne[i]);
    hisCalFactorSigmaTotal->Add(hisCalibrationFactorSigma[i]);

    hisCalibrationFactor[i]->Write();
    hisCalibrationFactorEne[i]->Write();
    hisCalibrationFactorRatio[i]->Write();
    hisCalibrationFactorSecondRatio[i]->Write();
    hisCalibrationFactorSigma[i]->Write();
    hisCalibrationFactorHeight[i]->Write();
    hisCalibrationFactorHeight_Weighted[i]->Write();
    ofs1 << i << "\t" << hisCalibrationFactor[i]->Integral() << std::endl;
    if( hisCalibrationFactor[i]->Integral() < 144 ){
      continue;
    }           
    hisCalibrationFactor[i]->Fit("gaus","Q","");
    TF1* calFunction= hisCalibrationFactor[i]->GetFunction("gaus");
    CalFactor[i]    = hisCalibrationFactor[i]->GetMean();
    CalRMS[i]       = hisCalibrationFactor[i]->GetRMS();
    //CalFactor[i] = calFunction->GetParameter(1);
    //CalRMS[i] = calFunction->GetParameter(2);    
    //ofs << i << "\t" << CalFactor[i]*preCalFactor[i] << std::endl;    
    hisReNormCal->Fill(CalFactor[i]*preCalFactor[i]);
  }

  // Renomalize Calibration Factor, Mean -> 1
  for( int i = 0; i< 2716; ++i ){
    if( hisCalibrationFactor[i]->Integral() < 144 ){continue;}
    ofs << i << "\t" 
	<< CalFactor[i]*preCalFactor[i]/hisReNormCal->GetMean()
	<< std::endl;
  }
  
  hisCalSigmaTotal->Write();
  hisCalFactorRatioTotal->Write();
  hisCalFactorSecondRatioTotal->Write();
  hisCalFactorEneTotal->Write();
  hisCalFactorSigmaTotal->Write();  
  
  tfOut->Close();
  ofs.close();      
  ofs1.close();
  
}
