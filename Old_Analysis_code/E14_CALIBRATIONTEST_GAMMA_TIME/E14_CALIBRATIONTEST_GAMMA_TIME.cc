// Apr. 2nd 2012
// Add Function which set Calibration Factor when Calculating Simulation Data

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"

#include "klong/Klong.h"

#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "E14ReadSumFile.h"
#include "E14_CALIBRATIONTEST_GAMMA/User_Function.h"

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

#include "E14_CALIBRATIONTEST_GAMMA/Calibration.h"
#include "E14_CALIBRATIONTEST_GAMMA/CalibrationTree.h"

//bool user_rec(std::list<Gamma> const &glist, std::vector<Klong> &klVec);
//void user_cut(E14GNAnaDataContainer &data  , std::vector<Klong> &klVec);


int
main(int argc,char** argv)
{
  //////////////////////////////////////////////////////////////////////////////////////
  // Read Arguement
  //////////////////////////////////////////////////////////////////////////////////////

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

  
  std::string BeamCalibrationFilename = "SimRoot/BeamCalibrationFactor.dat";
  Double_t BeamFactor[2716];
  for( int i = 0; i< 2716; ++i ){
    if( i >= 2240 && i< 2240+120){
      BeamFactor[i] = 0.;
    }else if ( i >= 2716-120 && i< 2716 ){
      BeamFactor[i]  =0.;
    }else{
      BeamFactor[i] = 1.;
    }
  }

  /*
  std::ifstream ifsBeamCal(BeamCalibrationFilename.c_str());
  if( ifsBeamCal.is_open()){
    Int_t chID;
    Double_t CalFactor;
    while( ifsBeamCal >> chID >> CalFactor ){
      if( CalFactor > 1.3 || CalFactor < 0.7 || CalFactor ==0 ){
	BeamFactor[chID] = 1.;	
      }else{
	BeamFactor[chID] = 1./CalFactor;
      }
    }
  }else{
    std::cerr <<" Error: BeamCalibration File is not exist" << std::endl;
    return -1;
  }
  */
  
  /*
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
  */
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Input RootFile
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  E14ReadSumFile* ReadSum = new  E14ReadSumFile(1);

  for( int subindex = 0; subindex < 10 ;++subindex){    
    inputFilename       = Form("SimRoot/Conv_KL3pi0.mac_1000000_%d.root",runNumber+subindex);
    ReadSum->Add(inputFilename.c_str());
  }

  outputFilename      = Form("SimCalibration_Data/TimeCalibration_%04d_%d.root",
			     runNumber,iterationNumber);
  calibrationFilename = Form("SimCalibration_Data/CalibrationFactor_%d.dat",
			     iterationNumber);

  std::cout<< "Total Event Number : " << ReadSum->GetEntries() << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Output RootFile
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  TFile *fout = new TFile(outputFilename.c_str(),"RECREATE");
  TTree *trout = new TTree("trCalibration","output from E14_CALIBRTION_CULSTER_BUILDER");
  Calibration* calibrator = new Calibration();

  Double_t SimGammaEne[6];
  Double_t SimGammaEneMatched[6];
  Double_t SimGammaPosx[6];
  Double_t SimGammaPosy[6];
  Double_t SimGammaPosz[6];
  Double_t SimGammaR[6];
  Double_t SimKLVtx[3];
  Double_t SimKLMom[3];
  Double_t SimKLE;

  Int_t    SimGammaCombination[6];  
  Int_t    KlongFit;
  
  trout->Branch("SimGammaEne" ,SimGammaEne ,"SimGammaEne[6]/D");
  trout->Branch("SimGammaEneMatched" ,SimGammaEneMatched,"SimGammaEneMatched[6]/D");
  trout->Branch("SimGammaPosx",SimGammaPosx,"SimGammaPosx[6]/D");
  trout->Branch("SimGammaPosy",SimGammaPosy,"SimGammaPosy[6]/D");
  trout->Branch("SimGammaPosz",SimGammaPosz,"SimGammaPosz[6]/D");
  trout->Branch("SimGammaR"   ,SimGammaR   ,"SimGammaR[6]/D");
  trout->Branch("SimKLVtx"    ,SimKLVtx    ,"SimKLVtx[3]/D");
  trout->Branch("SimKLMom"    ,SimKLMom    ,"SimKLMom[3]/D");
  trout->Branch("SimKLE"      ,&SimKLE     ,"SimKLE/D");
  trout->Branch("SimGammaCombination",SimGammaCombination,"SimGammaCombination[6]/I");
  trout->Branch("KlongFit"    ,&KlongFit,"KlongFit/I");
  
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
  
  TH1D* hisStage    = new TH1D("hisStage","Distribution of breakPoint",100,0,100);
  TH1D* hisStageSum = new TH1D("hisStageSum","Residual number of event untion break Point", 10, 0, 10);
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Read Calibration File
  /////////////////////////////////////////////////////////////////////////////////////////////////
  
  std::cout<< "Read Calibration File " << std::endl;
  
  for(int  ich = 0; ich < 3000; ich++){
    CSICalFactor[ich] = 1;
  }

  //  Test about to Calibration, Don't use calibration.. 
  /*
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
  */
  
  /////////////////////////////////////////////////////////////////////////////////////////////////
  // loop
  /////////////////////////////////////////////////////////////////////////////////////////////////
  //loop analysis  

  std::cout<< "Loop!!"<< std::endl;
  int nentry = ReadSum->GetEntries();    
  int nKL    = 0;
  int nCalib = 0;
  std::cout<<"# of entry in input tree =="<<nentry<<std::endl;  
  if( nloop<0 || nloop>nentry ) nloop = nentry;  
  std::cout<<"\n start loop analysis for "<<nloop<<" events..."<<std::endl; 
  
  //for(int ievt=0;ievt< 10;ievt++){
  int Stage;
  int SubStage;
  for(int ievt=0;ievt<nloop;ievt++){
    if(nloop>100&&ievt%(nloop/100)==0){
      std::cout<<ievt/(nloop/100)<<"%"<<std::endl;    
    }
    ReadSum->GetEntry(ievt);
    nCSIDigi = 0; 
    Stage    = 0;
    SubStage = 0; 
    data.eventID = ievt;
    calibrator->InitValue();
    calData.InitValue();
    
    // Do Something 
    hisStageSum->Fill(Stage);
    ++Stage;SubStage = 0;

    ///////////////////////////////////////////////////////////////////////////
    // Read Csi Data::Stage1
    ///////////////////////////////////////////////////////////////////////////
    for( int i = 0; i< ReadSum->CsiNumber; i++){
      Double_t Energy = ReadSum->CsiEne[i]*CSICalFactor[ReadSum->CsiModID[i]]*BeamFactor[ReadSum->CsiModID[i]];      
      if( Energy > 3 ){
	CSIDigiID[nCSIDigi]  = ReadSum->CsiModID[i];
	CSIDigiE[nCSIDigi]   = Energy;
	CSIDigiTime[nCSIDigi]= ReadSum->CsiTime[i];
	nCSIDigi++;
      }
    }
    
    hisStageSum->Fill(Stage);
    ++Stage;SubStage=0;
    ///////////////////////////////////////////////////////////////////////////
    // Clustering ::Stage2
    ///////////////////////////////////////////////////////////////////////////
    std::list<Cluster> clist;
    std::list<Gamma>   glist; 
    std::vector<Klong> klVec;
    clist = clusterFinder.findCluster(nCSIDigi,CSIDigiID,CSIDigiE,CSIDigiTime);
    gFinder.findGamma(clist, glist);    
    if( glist.size() != 6 ){ 
      hisStage->Fill(Stage*10+SubStage);
      continue; 
    }else{
      SubStage++;
    }
    
    if(!user_rec(glist,klVec)){
      hisStage->Fill(Stage*10+SubStage);
      continue;
    }else{
      SubStage++;
      nKL++;
    }
        
    data.setData( clist );
    data.setData( glist );    
    user_cut(data,klVec);
    data.setData(klVec);
    
    ////////////////////////////////////////////////////////////////////////
    /////// Set SimGammaData ::Stage3
    /////// Init SimGammaData .......
    ////////////////////////////////////////////////////////////////////////
    /*
      if((data.CutCondition & (1+4)) !=0 ){      
      hisStage->Fill(Stage*10+SubStage);
      continue;
    }else{
    SubStage = 1;
    }
    */

    
    for( int i = 0; i< 6; i++){
      SimGammaEne[i] = 0.;
      SimGammaEneMatched[i] = 0.;
      SimGammaPosx[i] = 0.;
      SimGammaPosy[i] = 0.;
      SimGammaPosz[i] = 0.;
      SimGammaCombination[i] = 0;  
      SimGammaR[i]    = 0.;
    }   
    SimKLVtx[0] = 0.;
    SimKLVtx[1] = 0.;
    SimKLVtx[2] = 0.;
    SimKLE      = 0.; 
    SimKLMom[0] = 0.; 
    SimKLMom[1] = 0.; 
    SimKLMom[2] = 0.; 
    KlongFit    = 0;
    
    if( ReadSum->nTrack != 10 ){ 
      hisStage->Fill(Stage*10+SubStage);
      continue; 
    }else{
      SubStage=2;
    }
    
    std::vector<int> pi0TrackNum;
    std::vector<CLHEP::Hep3Vector> gammaPos;
    std::vector<int> gammaTrackNum;
    
    int gammaCombinationTrackNum[3][2]={0};
    int gammaIndex[3]      = {0};
    int gammaFindNum[3][2] = {0};
    int gammaAgree[3]      = {0};
    int klongAgree[2]      = {0};     
    int nGamma=0;
    
    for( int i = 0 ;i < ReadSum->nTrack; ++i){
      
      // Set vertex of KL

      if( ReadSum->pid[i] == 130 ){
	SimKLVtx[0] = ReadSum->end_v[i][0];
	SimKLVtx[1] = ReadSum->end_v[i][1];
	SimKLVtx[2] = ReadSum->end_v[i][2];
	SimKLMom[0] = ReadSum->end_p[i][0];
	SimKLMom[1] = ReadSum->p[i][1];
	SimKLMom[2] = ReadSum->p[i][2];
	SimKLE      = ReadSum->ek[i];
      }
      
      if( ReadSum->pid[i] == 111 ){
	pi0TrackNum.push_back(i);
      }      
      if( ReadSum->pid[i] == 22 ){
	gammaTrackNum.push_back(i);
	if( nGamma > 6 ){ 
	  std::cerr << "Number of Gamma is more than 6" <<std::endl;	  
	  hisStage->Fill(Stage*10+SubStage);
	  continue;
	}else{
	  SubStage=3;
	}
	
	SimGammaEne[nGamma] = ReadSum->ek[i];
	SimGammaPosx[nGamma]= ReadSum->end_v[i][0];
	SimGammaPosy[nGamma]= ReadSum->end_v[i][1];
	SimGammaPosz[nGamma]= ReadSum->end_v[i][2];
	
	nGamma++;
	
      }
    }    
    
    if( pi0TrackNum.size() != 3 || gammaTrackNum.size() != 6 ){
      hisStage->Fill(Stage*10+SubStage);
      continue;
    }else{
      SubStage=4;
    }
    
    hisStageSum->Fill(Stage);
    ++Stage;SubStage = 0;
    ////////////////////////////////////////////////////////////////////////////
    // Set Combination of Gamma ::Stage4
    ////////////////////////////////////////////////////////////////////////////
    // find mother pi0 
    for( int i = 0; i< gammaTrackNum.size(); ++i){
      for( int j = 0; j < pi0TrackNum.size() ; ++j){
	if( ReadSum->mother[gammaTrackNum[i]] == ReadSum->track[pi0TrackNum[j]]){
	  gammaCombinationTrackNum[j][gammaIndex[j]] = gammaTrackNum[i];
	  gammaIndex[j]++;	      
	}
      }	  
    }
    
    for( int i = 0; i< pi0TrackNum.size() ; ++i){
      Double_t g1 = ReadSum->ek[gammaCombinationTrackNum[i][0]];
      Double_t g2 = ReadSum->ek[gammaCombinationTrackNum[i][1]];
      // Sort by GammaEnergy;
      if( g1 < g2){
	Int_t temp = gammaCombinationTrackNum[i][0];
	gammaCombinationTrackNum[i][0] = gammaCombinationTrackNum[i][1];
	gammaCombinationTrackNum[i][1] = temp;
      }
    }	
    
    hisStageSum->Fill(Stage);
    ++Stage;SubStage = 0;
    ////////////////////////////////////////////////////////////////////////////
    // Match with KL Data. ::Stage5
    ////////////////////////////////////////////////////////////////////////////
    
    int nKlongLoop= klVec.size();
    if( nKlongLoop >2){
      nKlongLoop  = 2; 
    }      
    for( int iklong = 0; iklong < nKlongLoop; iklong++){	  

      gammaAgree[0] = 0;
      gammaAgree[1] = 0;
      gammaAgree[2] = 0;
      
      for( int i = 0; i< klVec[iklong].pi0().size(); ++i){
	Double_t deltaR = 0.;
	Double_t minR   = 2000.; 	  
	Int_t    minimumID = 0;//invalid ID because this is track ID of Klong
	Double_t x;
	Double_t y;
	
	for( int j = 0; j < gammaTrackNum.size(); ++j){	    
	  
	  x = ReadSum->end_v[gammaTrackNum[j]][0];
	  y = ReadSum->end_v[gammaTrackNum[j]][1];
	  
	  deltaR 
	    = TMath::Sqrt(TMath::Power(TMath::Abs(klVec[iklong].pi0()[i].g1().x() - x),2)+
			  TMath::Power(TMath::Abs(klVec[iklong].pi0()[i].g1().y() - y),2));
	  if( deltaR < minR ) {
	    minR  = deltaR; 
	    minimumID = gammaTrackNum[j];
	    if( iklong ==0 ){
	      SimGammaEneMatched[i*2] = ReadSum->ek[gammaTrackNum[j]];
	    }	      
	  }
	}
	gammaFindNum[i][0] = minimumID;
	
	deltaR = 0.;
	minR   = 2000.;
	minimumID = 0;
	
	for( int j =0; j< gammaTrackNum.size(); ++j){
	  x = ReadSum->end_v[gammaTrackNum[j]][0];
	  y = ReadSum->end_v[gammaTrackNum[j]][1];
	  deltaR 
	    = TMath::Sqrt(TMath::Power(TMath::Abs(klVec[iklong].pi0()[i].g2().x() - x),2)+
			  TMath::Power(TMath::Abs(klVec[iklong].pi0()[i].g2().y() - y),2));
	  if( deltaR < minR ) { 
	    minR  = deltaR;
	    minimumID = gammaTrackNum[j];
	    if( iklong ==0 ){
	      SimGammaEneMatched[i*2+1] = ReadSum->ek[gammaTrackNum[j]];
	    }	      
	  }
	}	  
	gammaFindNum[i][1] = minimumID;
      }     
      
      for( int i = 0; i< 3; ++i){	  
	for( int j = 0; j< 3; ++j){
	  if((gammaFindNum[i][0] == gammaCombinationTrackNum[j][0] && 
	      gammaFindNum[i][1] == gammaCombinationTrackNum[j][1] ) ||
	     (gammaFindNum[i][0] == gammaCombinationTrackNum[j][1] &&
	      gammaFindNum[i][1] == gammaCombinationTrackNum[j][0] )){
	    gammaAgree[i] = 1;
	  }
	}
      }
      if(gammaAgree[0]==1 &&
	 gammaAgree[1]==1 &&
	 gammaAgree[2]==1 ){	
	klongAgree[iklong] = 1; 
	KlongFit |= 1<< iklong;
      }
    }//end match with kldata::Stage 5
    
    hisStageSum->Fill(Stage);
    ++Stage;SubStage = 0;    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Stage 6 : Test Calculation Timing of Gamma
    ////////////////////////////////////////////////////////////////////////////////////////////



    
    
    Stage=9;
    hisStage->Fill(Stage*10+9);     
    trout->Fill();    
  }
  
  // end of analysis
  trout->Write();
  hisStage->Write();
  hisStageSum->Write();
  std::cout << "nKL:"<<  nKL <<std::endl;
  std::cout << "nCalib:"<<nCalib << std::endl;
  fout->Close();
  std::cout<<"finish!"<<std::endl;    
  return 0;
  
}
