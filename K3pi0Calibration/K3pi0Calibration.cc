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

#include "E14ReadSumFile.h"
#include "ReadWavAna.h"
#include "User_Function.h"
//#include "E14_CALIBRATION_CLUSTER_BUILDER_V1/KL_calibration.h"
#include "Calibration.h"
#include "CalibrationTree.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <list>

//bool user_rec(std::list<Gamma> const &glist, std::vector<Klong> &klVec);
//void user_cut(E14GNAnaDataContainer &data  , std::vector<Klong> &klVec);

Int_t 
main(int argc,char** argv)
{
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Read Arguement
  /////////////////////////////////////////////////////////////////////////////////////////////////

  int nloop = -1;
  int runNumber;
  int iterationNumber;
  std::string inputFilename;
  std::string outputFilename;
  std::string calibrationFilename;
  std::string TempCalibrationFilename;
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB  = std::getenv("ANALYSISLIB");
  std::string ROOTFILE_3PI0CALIBRATION = std::getenv( "ROOTFILE_3PI0CALIBRATION");

  std::string path;
  if( argc == 3 ){
    runNumber = atoi(argv[1]);
    iterationNumber = atoi(argv[2]);
  }else if(argc == 4 ){
    runNumber = atoi(argv[1]);
    iterationNumber = atoi(argv[2]);
    path = argv[3];
  }else{
    std::cerr << "<<<>>>Arguement Error<<<>>>" <<"\n"
	      << "Usage:: " << argv[0] 
	      << " : [runNumber] [iterationNumber]" << std::endl;
    return -1;
  }
  
  inputFilename       = Form("%s/run_wav_%04d_cl.root",ROOTFILE_WAV.c_str(),runNumber);

  if( argc  == 3 ){
    outputFilename      = Form("%s/CalibrationADV_%04d_%d.root",ROOTFILE_3PI0CALIBRATION.c_str(),runNumber,iterationNumber);
    calibrationFilename = Form("%s/CalibrationFactorADV_%d.dat",ROOTFILE_3PI0CALIBRATION.c_str(),iterationNumber);
  }else if( argc == 4 ){
    outputFilename      = Form("%s/%s/CalibrationADV_%04d_%d.root",ROOTFILE_3PI0CALIBRATION.c_str(),path.c_str(),runNumber,iterationNumber);
    calibrationFilename = Form("%s/%s/CalibrationFactorADV_%d.dat",ROOTFILE_3PI0CALIBRATION.c_str(),path.c_str(),iterationNumber);
  }

  std::cout<<"Input file        : "<< inputFilename        << std::endl;
  std::cout<<"Output file       : "<< outputFilename       << std::endl;
  std::cout<<"Calibration Number: "<< calibrationFilename  << std::endl;

  TempCalibrationFilename = Form("%s/Data/Temperature_Factor/TemperatureCorrectionFactor.root",ANALYSISLIB.c_str());
  TFile* tfTempCorr  =new TFile(TempCalibrationFilename.c_str());
  TTree* trTempCorr  =(TTree*)tfTempCorr->Get("TemperatureCorrectionCsI");
  Double_t TempCorFactor=0;
  trTempCorr->SetBranchAddress("CorrectionFactor",&TempCorFactor);
  trTempCorr->GetEntry(runNumber);
  if( TempCorFactor == 0){
    TempCorFactor = 1;
  }
  std::cout<< TempCorFactor << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Input RootFile
  /////////////////////////////////////////////////////////////////////////////////////////////////
    
  //E14ReadSumFile* ReadSum = new  E14ReadSumFile();
  TFile* tfinput = new TFile(inputFilename.c_str());
  TTree* ch = (TTree*)tfinput->Get("T");

  Int_t CsiNumber;
  Double_t CsiEne[2716];//CsiNumber
  Double_t CsiTime[2716];//CsiNumber
  Int_t CsiModID[2716];//CsiNumber
  ch->SetBranchAddress("CsiNumber", &CsiNumber );
  ch->SetBranchAddress("CsiTime",CsiTime);//CsiNumber
  ch->SetBranchAddress("CsiModID",CsiModID);//CsiNumber
  ch->SetBranchAddress("CsiEne",CsiEne);//CsiNumber
  

  /*

  TChain* ch = new TChain("T");
  ch->Add( inputFilename.c_str());

  */


  //ReadWavAna *ReadSum = new ReadWavAna(ch);





  //ReadSum->Add(inputFilename.c_str());

  /*
  for( int i = 0; i< 20 ; i++){
    inputFilename       = Form("sumdata/Sum%04d.root",runNumber+i);
  }
  */
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Output RootFile
  /////////////////////////////////////////////////////////////////////////////////////////////////

  TFile *fout = new TFile(outputFilename.c_str(),"RECREATE");  
  TTree *trout = new TTree("trCalibration","output from E14_CALIBRTION_CULSTER_BUILDER"); 
  Calibration* calibrator  = new Calibration();  
  E14GNAnaDataContainer data;
  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  CalibrationTree calData;
  calData.Branch(trout);

  //data.branchOfGammaList( trout );
  data.branchOfKlong( trout );
  GammaFinder gFinder;
  
  // declare  ClusterFinder and variables
  int nCSIDigi=0;
  int CSIDigiID[3000]={0};
  double CSIDigiE[3000]={0},CSIDigiTime[3000]={0};
  double CSICalFactor[3000]={0};
  ClusterFinder clusterFinder;
  double CC03IntegratedADC[32];

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Read Calibration File
  /////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout<< "Read Calibration File " << std::endl;
  
  for(int  ich = 0; ich < 3000; ich++){
    CSICalFactor[ich] = 1;
  }
  if( iterationNumber  != 0){
    std::ifstream ifs(calibrationFilename.c_str());
    if(!ifs.is_open()){
      std::cerr <<"File is not Exist" << std::endl;
      return -1;
    }
    int    id;
    double CalibrationFactor;;
    while( !ifs.eof() ){
      if(ifs >> id >> CalibrationFactor){
	if( !(CalibrationFactor==0 ) ){
	  CSICalFactor[id] = CalibrationFactor;
	}else{
	  CSICalFactor[id] = 1;
	}
      }
    }    
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Prepare Calibration
  /////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout << "Prepare Calibration" << std::endl;
  /*
  TH1D* his_Fit[100];
  TH1D* his_CSI[N_CSI];
  TH2D* his_CSIEne[N_CSI];
  TH2D* his_CSISigma[N_CSI];
  TH2D* his_CSIRatio[N_CSI];
  TH2D* his_CSISecondRatio[N_CSI];
  std::cout<< "Set Histogram" << std::endl;
  for( int icsi = 0; icsi < N_CSI; ++icsi){
    std::cout<< icsi << std::endl; 
    his_CSI[icsi]
      = new TH1D(Form("CalibrationFactor_%d",icsi),
		 Form("CalibrationFacotr_%d",icsi),
		 200, 0, 2 );
    his_CSIEne[icsi]
      = new TH2D(Form("CalibrationFactorEne_%d",icsi),
		 Form("CalibrationFactorEne_%d;Energy[Mev];CalibrationFactor",icsi),
		 200, 0, 4000,200, 0, 2); 
    his_CSISigma[icsi]
      = new TH2D(Form("CalibrationFactorSigma_%d",icsi),
		 Form("CalibrationFactorSigma_%d;Sigma Of Gamma;CalibrationFactor",icsi),
		 200,0,50,200,0,2);
    his_CSIRatio[icsi] 
      = new TH2D(Form("CalibrationFactorRatio_%d",icsi),
		 Form("CalibrationFactorRatio_%d;Ratio OF Crystal;CalibrationFactor",icsi),
		 100,0,1,200,0,2);
    his_CSISecondRatio[icsi]
      = new TH2D(Form("CalibrationFactorSecondRatio_%d",icsi),
		 Form("CalibrationFacotrSecondRatio_%d;Ratio Of SecondCrystal;CalibrationFactor",icsi),
		 100,0,1,200,0,2);
  }
  */

  std::cout << "End Prepare Calibration" << std::endl;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // loop
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //loop analysis
  
  std::cout<< "Loop!!"<< std::endl;
  int nentry = ch->GetEntries();    
  int nKL    = 0; 
  std::cout<<"# of entry in input tree =="<<nentry<<std::endl;
  
  if( nloop<0 || nloop>nentry ) nloop = nentry;  
  std::cout<<"\n start loop analysis for "<<nloop<<" events..."<<std::endl; 
  //nloop = 5000;
  
  for(int ievt=0;ievt<nloop;ievt++){
    //std::cout << ievt << std::endl;
    if(nloop>100&&ievt%(nloop/100)==0)
      std::cout<<ievt/(nloop/100)<<"%"<<std::endl;
    
    //std::cout<< "ievt:"<< ievt << std::endl;
    ch->GetEntry(ievt);
    nCSIDigi = 0; 
    calibrator->InitValue();
    calData.InitValue();
    //std::cout<< "loop" << std::endl;
    for( int i = 0; i< CsiNumber; i++){
      Double_t Energy = CsiEne[i]*CSICalFactor[ (CsiModID[i]) ]/TempCorFactor;      
      if( Energy > 3 ){
	CSIDigiID[nCSIDigi]  = CsiModID[i];
	CSIDigiE[nCSIDigi]   = Energy;
	CSIDigiTime[nCSIDigi]= CsiTime[i];
	nCSIDigi++;
      }
      /*
      std::cout << CsiNumber << " : " << i << " : " 
		<< CsiModID[i] << " : "
		<< CsiEne[i] << " : "
		<< CsiTime[i] << std::endl;
      */
    }

    /*   
    int hitCC03=0;
    for( int i = 0; i< ReadSum->CC03Number; i++){
      CC03IntegratedADC[i] = ReadSum->CC03Ene[i];      
      //assume CC03has Uniform Gain;
      if( CC03IntegratedADC[i] > 1000){
	hitCC03++;
      }      
    }
    if( hitCC03> 0) {
      data.eventID++;
      continue;
    }
    */    
    //std::cout<< "Build" << std::endl;
    std::list<Cluster> clist;
    std::list<Gamma>   glist; 
    std::vector<Klong> klVec;
    
    clist = clusterFinder.findCluster(nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);
    gFinder.findGamma(clist, glist);
    if( glist.size() == 6 ){
      if(user_rec(glist,klVec)){
	data.setData( clist );
	data.setData( glist );    
	user_cut(data,klVec);
	data.setData(klVec);    
	///change here JWLEE
	//Fiducial cut
	bool GammaFlag = false;
	for( int i = 0; i< klVec[0].pi0().size(); i++){
	  if( TMath::Abs(klVec[0].pi0()[i].g1().y()) > 550){
	    GammaFlag = true;
	  }
	  if( TMath::Abs(klVec[0].pi0()[i].g2().y()) > 550){
	    GammaFlag = true; 
	  }
	}
	//if( (data.CutCondition & 2) != 0 ){ continue; }
	//if( GammaFlag ){ continue;} 

	int result = calibrator->CalEnergy_idv(klVec);
	if(result >= 1){
	  nKL++;
	  calibrator->GetResult(calData);
	}
	
	trout->Fill();	    
      }
    }else{
      data.CutCondition = -1;
    }
    data.eventID++;    
  }

  // end of analysis
  /*
  for( int icsi = 0; icsi < N_CSI; ++icsi){
    his_CSI[icsi]->Write();
    his_CSIEne[icsi]->Write();
    his_CSISigma[icsi]->Write();
    his_CSIRatio[icsi]->Write();
    his_CSISecondRatio[icsi]->Write();
  }
  */

  trout->Write();
  fout->Close();
  std::cout<<"finish!"<<std::endl;  
  return 0;

}
