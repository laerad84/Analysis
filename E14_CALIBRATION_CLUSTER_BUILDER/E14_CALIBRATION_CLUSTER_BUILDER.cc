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

#include "E14_CALIBRATION_CLUSTER_BUILDER/E14ReadSumFile.h"
#include "E14_CALIBRATION_CLUSTER_BUILDER/User_Function.h"
#include "E14_CALIBRATION_CLUSTER_BUILDER/KL_calibration.h"

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

  if( argc == 3 ){
    runNumber = atoi(argv[1]);
    iterationNumber = atoi(argv[2]);
  }else{
    std::cerr << "<<<>>>Arguement Error<<<>>>" <<"\n"
	      << "Usage:: " << argv[0] 
	      << " : [runNumber] [iterationNumber]" << std::endl;
    return -1;
  }
  
  inputFilename       = Form("sumup_data/Sum%04d.root",runNumber);
  outputFilename      = Form("Calibration_Data/Calibration_%04d_%d.root",
			     runNumber,iterationNumber);
  calibrationFilename = Form("Calibration_Data/CalibrationFactor_%d.dat",
			     iterationNumber);
  
  std::cout<<"input file: "        << inputFilename        << std::endl;
  std::cout<<"output file: "       << outputFilename       << std::endl;
  std::cout<<"Calibration Number: "<< calibrationFilename  << std::endl;

  TempCalibrationFilename = "TemperatureCorrectionFactor.root";
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
  
  E14ReadSumFile* ReadSum = new  E14ReadSumFile();
  ReadSum->Add(inputFilename.c_str());

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
  E14GNAnaDataContainer data;
  data.branchOfClusterList( trout );
  data.branchOfDigi( trout );
  //data.branchOfGammaList( trout );
  data.branchOfKlong( trout );
  GammaFinder gFinder;
  
  // declare  ClusterFinder and variables
  int nCSIDigi=0;
  int CSIDigiID[3000]={0};
  double CSIDigiE[3000]={0},CSIDigiTime[3000]={0};
  double CSICalFactor[3000]={0};
  ClusterFinder clusterFinder;
  
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
  TH1D* his_Fit[100];
  TH1D* his_CSI[N_CSI];
  TH2D* his_CSIEne[N_CSI];

  user_ana_recg6_init(his_Fit,his_CSI);
  std::cout << his_CSI[0] << std::endl;
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // loop
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //loop analysis
  std::cout<< "Loop!!"<< std::endl;
  int nentry = ReadSum->GetEntries();    
  int nKL    = 0; 
  std::cout<<"# of entry in input tree =="<<nentry<<std::endl;
  
  if( nloop<0 || nloop>nentry ) nloop = nentry;  
  std::cout<<"\n start loop analysis for "<<nloop<<" events..."<<std::endl; 
  
  for(int ievt=0;ievt<nloop;ievt++){
    if(nloop>100&&ievt%(nloop/100)==0)
      std::cout<<ievt/(nloop/100)<<"%"<<std::endl;
    
    ReadSum->GetEntry(ievt);
    nCSIDigi = 0; 
    
    for( int i = 0; i< ReadSum->CsiNumber; i++){
      Double_t Energy = ReadSum->CsiEne[i]*CSICalFactor[ReadSum->CsiModID[i]]*TempCorFactor;      
      if( Energy > 3 ){
	CSIDigiID[nCSIDigi]  = ReadSum->CsiModID[i];
	CSIDigiE[nCSIDigi]   = Energy;
	CSIDigiTime[nCSIDigi]= ReadSum->CsiTime[i];
	nCSIDigi++;
      }
    }

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
	trout->Fill();
	if((data.CutCondition & (1+8)) ==0 ){
	  if( klVec[0].chisqZ() < 10 ){
	    if( klVec.size() >= 2 ){
	      if( (klVec[1].chisqZ() - klVec[0].chisqZ()) > 5 ){		
		int result = CalEnergy_idv(klVec[0],his_Fit,his_CSI);
		if(result == 1) nKL++;
	      }
	    }else{
		int result = CalEnergy_idv(klVec[0],his_Fit,his_CSI);
		if(result == 1) nKL++;
	    }
	  }
	}
      }
    }else{
      data.CutCondition = -1;
    }
    data.eventID++;    
  }
  // end of analysis
  trout->Write();
  std::cout << "nKL:"<<  nKL <<std::endl;
  for( int i = 0; i < N_CSI; i++){
    if( his_CSI[i]->Integral() > 0)
      std::cout <<"Write()" << i <<  " : "<< his_CSI[i]->Integral()<<  std::endl;
    his_CSI[i]->Write();
  }
  for( int i = 0; i< 5; i++){
    his_Fit[i]->Write();
  }
  
  fout->Close();
  std::cout<<"finish!"<<std::endl;  
  return 0;
}
