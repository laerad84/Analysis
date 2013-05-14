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
  OutputCalHist     += "/TestDistribution_%s_%d.root";

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
  std::string listFilename=runListFile.substr(runListFile.find_last_of('/')+1);
  std::cout<< listFilename << std::endl;
  TFile* tfOut 
    = new TFile(Form(OutputCalHist.c_str(),
		     listFilename.substr(0,-4).c_str(),
		     CalibrationNumber),
		"RECREATE");
  
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


  E14GNAnaDataContainer data;
  data.setBranchAddress(ch);
  

  ////////////////////////////////////////////////////////////////////////////
  // Loop
  ////////////////////////////////////////////////////////////////////////////
  //- Add Histogram for all Run;


  TH2D* hisPi0DeltaZ[2];
  char* Name[2]={"LowHeight","HighHeight"};
  for( int i = 0 ; i < 5; i++){
    hisPi0DeltaZ[i] = new TH2D(Form("hisPi0DeltaZ_%d",i),Form("hisPi0DeltaZ_%s",Name[i]),40,0,2000,240,-300,300);
  }


  for( int ievet  = 0 ;ievet  < Entries ; ++ievet){
    if( ievet %1000  == 0 ){
      std::cout << ievet << "/" << Entries << std::endl;
    }
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::list<Pi0>     plist;
    std::vector<Klong> klVec;

    ch->GetEntry(ievet);

    data.getData( clist );
    data.getData( glist );
    data.getData( plist );
    data.getData( klVec );
    
    if( LeadingHeight[i*2]<5000){ 
      hisPi0DeltaZ[0]->Fill(klVec[0].pi0()[0].g1().e(),klVec[0].pi0()[0].recZ()-klVec[0].pi0()[0].vz());
    }else{
      hisPi0DeltaZ[1]->Fill(klVec[0].pi0()[0].g1().e(),klVec[0].pi0()[0].recZ()-klVec[0].pi0()[0].vz());
    }
  }

  hisPi0DeltaZ[0]->Write();
  hisPi0DeltaZ[1]->Write();
  tfOut->Close();
}
