#include "E14RawData.h"
#include <iostream>
#include <iomanip>

#if !defined(__CINT__)
ClassImp(E14RawData)
#endif

E14RawData::E14RawData() : TObject()
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
  
void E14RawData::initializeDataValues()
{

  for(int i=0;i<20;i++){
    for(int j=0;j<16;j++){
      for(int k=0;k<48;k++) Data[i][j][k] = 0;
      Pedestal[i][j] = 0;
      PeakHeight[i][j] = 0;
      PeakTime[i][j] = 0;
      ADCSum[i][j] = 0;
      pe[i][j] = 0;
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
  std::cout << std::setw(6) << ADCSum[i][j];
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
  CalcADCsum();
  CalcEtSumFADC();

  return 0;
}

int E14RawData::CalcPedestal(){

  for(int i=0;i<nFADC;i++){
    for(int j=0;j<16;j++){
      Pedestal[i][j] = (float)( Data[i][j][0] + Data[i][j][1] ) / 2;
    }
  }
    
  return 0;
}

int E14RawData::CalcPeakHeight(){

  float MaxValue;
  float MaxValTime;
  
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
      
      PeakHeight[i][j] = MaxValue - Pedestal[i][j];
      PeakTime[i][j] = MaxValTime;
      
    }
  }    
  
  return 0;
}

int E14RawData::CalcADCsum(){

  for(int i=0;i<nFADC;i++){
    for(int j=0;j<16;j++){
      for(int k=0;k<nSamples;k++){
	ADCSum[i][j] += Data[i][j][k] - Pedestal[i][j];
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
