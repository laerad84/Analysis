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
  std::string ROOTFILE_SIMCONV = "/group/had/koto/ps/klea/work/jwlee/RootFiles/Data/Simulation/3pi0Run/SIM3PI0";
  std::string path;
  std::string InitialCalibrationFilename;
  
  


  Double_t ScaleFactor = 0;
  if( argc == 4 ){
    runNumber                  = atoi(argv[1]);
    iterationNumber            = atoi(argv[2]);
    InitialCalibrationFilename = argv[3];
    outputFilename      = Form("%s/CalibrationADV_%d_%d.root",ROOTFILE_3PI0CALIBRATION.c_str(),runNumber,iterationNumber);
    calibrationFilename = Form("%s/CalibrationFactorADV_%d.dat",ROOTFILE_3PI0CALIBRATION.c_str(),iterationNumber);
  }else if(argc == 5 ){
    runNumber           = atoi(argv[1]);
    iterationNumber     = atoi(argv[2]);
    path                = argv[3];
    InitialCalibrationFilename = argv[4];
    outputFilename      = Form("%s/%s/CalibrationADV_%d_%d.root",ROOTFILE_3PI0CALIBRATION.c_str(),path.c_str(),runNumber,iterationNumber);
    calibrationFilename = Form("%s/%s/CalibrationFactorADV_%d.dat",ROOTFILE_3PI0CALIBRATION.c_str(),path.c_str(),iterationNumber);
  }else if( argc  == 6 ){
    runNumber           = atoi(argv[1]);
    iterationNumber     = atoi(argv[2]);
    path                = argv[3];
    InitialCalibrationFilename = argv[4];
    ScaleFactor         = atof(argv[5]);    

   outputFilename      = Form("%s/%s/CalibrationADV_%d_%d.root",ROOTFILE_3PI0CALIBRATION.c_str(),path.c_str(),runNumber,iterationNumber);
    calibrationFilename = Form("%s/%s/CalibrationFactorADV_%d.dat",ROOTFILE_3PI0CALIBRATION.c_str(),path.c_str(),iterationNumber);
  }else{
    std::cerr << "<<<>>>Arguement Error<<<>>>" <<"\n"
	      << "Usage:: " << argv[0] 
	      << " : [runNumber] [iterationNumber]" << std::endl;
    return -1;
  }
  //inputFilename       = Form("%s/run_wav_%04d_cl.root",ROOTFILE_WAV.c_str(),runNumber);
  //inputFilename         = Form("%s/Sim3pi0_wav_fast_KL_RES_LY_pe_5E6_%d_Calibration.root",ROOTFILE_SIMCONV.c_str(),runNumber);
  //inputFilename       = Form("%s/Sim3pi0_wav_fast_KL_RES_LY_pe_5E6_%d_Calibration_mis_1.root",ROOTFILE_SIMCONV.c_str(),runNumber);// Test MisCalibration
  inputFilename       = Form("%s/Sim3pi0_wav_fast_5E6_%d.root",ROOTFILE_SIMCONV.c_str(),runNumber);
  std::cout<<"********************************************************\n";
  std::cout<<"Input file        : "<< inputFilename        <<"\n";
  std::cout<<"Output file       : "<< outputFilename       <<"\n";
  std::cout<<"RunNumber         : "<< runNumber            <<"\n";
  std::cout<<"IterationNumber   : "<< iterationNumber      <<"\n";
  std::cout<<"Calibration Number: "<< calibrationFilename  <<"\n";
  std::cout<<"Path              : "<< path                 <<"\n";
  std::cout<<"ScaleFactor       : "<< ScaleFactor          <<"%\n";
  std::cout<<"********************************************************\n" << std::endl;

  TempCalibrationFilename = Form("%s/Data/Temperature_Factor/TemperatureCorrectionFactor.root",ANALYSISLIB.c_str());
  TFile* tfTempCorr  =new TFile(TempCalibrationFilename.c_str()); 
  if( !tfTempCorr->IsOpen() ){ std::cout<< Form("%s is not opened",tfTempCorr->GetName()) << std::endl; return -1;}
  TTree* trTempCorr  =(TTree*)tfTempCorr->Get("TemperatureCorrectionCsI");
  
  Double_t TempCorFactor=1;
  /*
  trTempCorr->SetBranchAddress("CorrectionFactor",&TempCorFactor);
  trTempCorr->GetEntry(runNumber);
  if( TempCorFactor == 0){
    TempCorFactor = 1;
   }
  */


  std::cout<< "TempCorFactor :" << TempCorFactor << std::endl;
  std::cout<< "ScaleFactor   :" << ScaleFactor << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Read Initial Calibration Factor File
  /////////////////////////////////////////////////////////////////////////////////////////////////
  double CsICalFactor[2716];
  double CsIInitCalFactor[2716];

  std::cout<< "Set Initial Calibration Factor" << std::endl;
  int tmpID;
  double tmpCalFactor;
  for( int ich = 0; ich <2716; ich++){
    CsIInitCalFactor[ich] = 1;
  }
  if( argc >=4 ){
    std::ifstream ifsinitCal(InitialCalibrationFilename.c_str());
    if( !(ifsinitCal.is_open()) ){ std::cout<< InitialCalibrationFilename << " can't open" << std::endl; return -1; }
    while(ifsinitCal >> tmpID >> tmpCalFactor ){
      CsIInitCalFactor[tmpID] = tmpCalFactor;
    }
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Read Calibration File
  /////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout<< "Read Calibration File " << std::endl;
  
  for(int  ich = 0; ich < 2716; ich++){
    CsICalFactor[ich] = 1;
  }
  if( iterationNumber  != 0){
    std::ifstream ifs(calibrationFilename.c_str());
    if(!ifs.is_open()){
      std::cerr << Form("%s is ont opened", calibrationFilename.c_str()) << std::endl;
      std::cerr <<"File is not Exist" << std::endl;
      return -1;
    }

    int    id;
    double CalibrationFactor;;
    while(ifs >> id >> CalibrationFactor){
      if( !(CalibrationFactor==0 ) ){
	CsICalFactor[id] = CalibrationFactor;
      }else{
	CsICalFactor[id] = 1;
      }
    }
  }    

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Input RootFile
  /////////////////////////////////////////////////////////////////////////////////////////////////
    
  //E14ReadSumFile* ReadSum = new  E14ReadSumFile();

  //TFile* tfinput = new TFile(inputFilename.c_str());
  //if(!tfinput->IsOpen()){ std::cout<< Form("%s is not opened",tfinput->GetName()) << std::endl; return -1;}
  //TTree* ch = (TTree*)tfinput->Get("Tree");
  TChain* ch = new TChain("Tree");
  for( int i = 0; i< 25;i++){
    inputFilename       = Form("%s/Sim3pi0_wav_fast_5E6_%d.root",ROOTFILE_SIMCONV.c_str(),runNumber*25+i);
    ch->Add(inputFilename.c_str());
  }

  Int_t    CsiNumber;
  Double_t CsiEne[2716];//CsiNumber
  Double_t CsiTime[2716];//CsiNumber
  Double_t CsiSignal[2716];//CsiNumber
  Int_t    CsiModID[2716];//CsiNumber
  ch->SetBranchAddress("CsiNumber",&CsiNumber );
  ch->SetBranchAddress("CsiTime"  ,CsiTime);//CsiNumber
  ch->SetBranchAddress("CsiModID" ,CsiModID);//CsiNumber
  ch->SetBranchAddress("CsiEne"   ,CsiEne);//CsiNumber
  ch->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber

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
  if( !fout->IsOpen() ){ std::cout<< Form("%s is not opened",fout->GetName()) << std::endl; return -1; }
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
  int nCSIDigi              =0;
  int CSIDigiID[3000]       ={0};
  double CSIDigiE[3000]     ={0};
  double CSIDigiTime[3000]  ={0};
  double CSIHeight[3000]    ={0};
  ClusterFinder clusterFinder;
  double CC03IntegratedADC[32];

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
  std::cout << "ScaleFactor: "<<  ScaleFactor << std::endl;
  double Scale = 1+0.01*ScaleFactor;
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
      Double_t Energy = CsiEne[i]*CsICalFactor[(CsiModID[i])]/TempCorFactor/CsIInitCalFactor[(CsiModID[i])];   
      if( Energy > 3 ){
	CSIDigiID[nCSIDigi]  = CsiModID[i];
	CSIDigiE[nCSIDigi]   = Energy;
	CSIDigiTime[nCSIDigi]= CsiTime[i];
	CSIHeight[nCSIDigi]  = CsiSignal[i];
	//std::cout << CSIHeight[nCSIDigi] <<"\t" << CsiSignal[nCSIDigi] << std::endl;
	nCSIDigi++;
      }
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
	data.setData(klVec);    
	user_cut(data,klVec);
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
	  int tmpCnt = 0; 
	  for( int iDigi = 0; iDigi < nCSIDigi; iDigi++){
	    for( int iGID = 0; iGID < 6; iGID++){
	      if( CSIDigiID[iDigi] == calData.LeadingChID[iGID]){
		calData.LeadingHeight[iGID]=CSIHeight[iDigi];
		//std::cout << iGID << "\t" << CSIHeight[iDigi] << "\t" << calData.LeadingHeight[iGID] << std::endl;
		tmpCnt++;
	      }
	    }
	    if(tmpCnt==6){break;}
	  }
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
