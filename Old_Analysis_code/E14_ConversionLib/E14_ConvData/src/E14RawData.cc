#include "E14RawData.h"
#include <iostream>
#include <iomanip>

#if !defined(__CINT__)
ClassImp(E14RawData)
#endif

E14RawData::E14RawData()
{
  nFADC = 0;
  nSamples = 0;
  EventNo = 0;
  //  BufferLoopNo = 0;

  initializeDataValues();
}


E14RawData::~E14RawData()
{
  ;
}

void E14RawData::Clear(Option_t* )
{
  initializeDataValues();
}

void E14RawData::SetRunNum( int val )
{
  RunNum = val;
}

void E14RawData::InitTree( char* s )
{
  o_rootFile = s;
  //  hfile = new TFile("/dev/shm/a.root", "RECREATE","",1);
  hfile = new TFile( o_rootFile.c_str(), "RECREATE");
  tree=new TTree("EventTree","Conv data for 2012 run");

  tree->Branch("nFADC",&nFADC,"nFADC/S");// set 20 as nFADC
  tree->Branch("nSamples",&nSamples,"nSamples/S");
  tree->Branch("EventNo",&EventNo,"EventNo/I");
  tree->Branch("Data",Data,"Data[20][16][48]/S");
  tree->Branch("Pedestal",Pedestal,"Pedestal[20][16]/F");
  tree->Branch("PeakHeight",PeakHeight,"PeakHeight[20][16]/S");
  tree->Branch("PeakTime",PeakTime,"PeakTime[20][16]/S");
  tree->Branch("IntegratedADC",IntegratedADC,"IntegratedADC[20][16]/F");
  tree->Branch("EtSum_FADC",EtSum_FADC,"EtSum_FADC[20][48]/I");

  tree->Branch("BufferLoopNo",&BufferLoopNo,"BufferLoopNo/L");
  tree->Branch("Error",Error,"Error[20]/S");
  tree->Branch("TimeStamp",TimeStamp,"TimeStamp[20]/I");
  tree->Branch("TrigNo",TrigNo,"TrigNo[20]/S");
  tree->Branch("SpillNo",SpillNo,"SpillNo[20]/S");
  tree->Branch("SlotNo",SlotNo,"SlotNo[20]/S");
  tree->Branch("Compression_flag",Compression_flag,"Compression_flag[20][16]/S");
    
}

void E14RawData::Fill()
{
  tree -> Fill();
}

void E14RawData::Write()
{
  //  hfile = new TFile( o_rootFile.c_str(), "RECREATE");
  hfile -> cd();
  tree->Write();
  hfile -> Close();
}
  
void E14RawData::initializeDataValues()
{

  for(int i=0;i<20;i++){
    for(int j=0;j<16;j++){
      for(int k=0;k<48;k++) Data[i][j][k] = 0;
      Pedestal[i][j] = 0;
      PeakHeight[i][j] = 0;
      PeakTime[i][j] = 0;
      IntegratedADC[i][j] = 0;
      //      pe[i][j] = 0;
      Compression_flag[i][j] = 0;
    }
    for(int k=0;k<48;k++) EtSum_FADC[i][k] = 0;
    Error[i] = 0;
    TimeStamp[i] = 0;
    TrigNo[i] = 0;
    SpillNo[i] = 0;
    SlotNo[i] = 0;
  }

}


void E14RawData::dump(int i=0,int j=0)
{
  std::cout << std::setw(4) << Data[i][j][0];
  std::cout << std::setw(6) << Pedestal[i][j];
  std::cout << std::setw(6) << IntegratedADC[i][j];
  std::cout << std::setw(6) << PeakHeight[i][j];
  std::cout << std::setw(6) << PeakTime[i][j];
  std::cout << std::endl;
}

int E14RawData::CalcAll(){

  if( nFADC == 0 ){
    std::cout << "nFAC is 0. Need to be set first." << std::endl;
    return 1;
  }

  CalcPedestal();
  CalcPeakHeight();
  CalcIntegratedADC();
  CalcEtSumFADC();

  return 0;
}

int E14RawData::CalcPedestal(){

  for(int i=0;i<nFADC;i++){
    for(int j=0;j<16;j++){
      Pedestal[i][j] = (float)( Data[i][j][0] + Data[i][j][1] + Data[i][j][2] ) / 3;
    }
  }
    
  return 0;
}

int E14RawData::CalcPeakHeight(){

  int MaxValue;
  int MaxValTime;
  
  for(int i=0;i<nFADC;i++){
    for(int j=0;j<16;j++){

      MaxValue   = 0;
      MaxValTime = 0;
      for(int k=0;k<nSamples;k++){	
	if( MaxValue < Data[i][j][k] ){
	  MaxValue = Data[i][j][k];
	  MaxValTime = k;
	}
      }
      
      PeakHeight[i][j] = MaxValue - (int)Pedestal[i][j];
      PeakTime[i][j] = MaxValTime;
      
    }
  }    
  
  return 0;
}

int E14RawData::CalcIntegratedADC(){

  for(int i=0;i<nFADC;i++){
    for(int j=0;j<16;j++){
      for(int k=0;k<nSamples;k++){
	IntegratedADC[i][j] += Data[i][j][k] - Pedestal[i][j];
      }
    }
  }
  
  return 0;
}

int E14RawData::CalcEtSumFADC(){

  for(int i=0;i<nFADC;i++){
    for(int j=0;j<16;j++){
      for(int k=0;k<nSamples;k++){
	EtSum_FADC[i][k] += Data[i][j][k];
      }
    }
  }

  return 0;
}
