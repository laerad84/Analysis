#include "E14Calibrator.h"
#include <iostream>
#include <iomanip>

E14Calibrator::E14Calibrator() : TObject()
{

  AllNum = 0;
  DataStoreTag = false;

  ThisID = 9999;
}


E14Calibrator::~E14Calibrator()
{
  ;
}

void E14Calibrator::initializeDataValues(){

  for(int i=0;i<4096;i++){
    GainMean[i]  = 9999;
    GainSigma[i] = 9999;
  }

}

int E14Calibrator::dump(){

  printf("\nDump calibration data\n");
  for(int i=0;i<AllNum;i++){
    printf("e14 ID = %d : ",i);
    printf("Gain Mean = %f : ",GainMean[i]);
    printf("Gain Sigma = %f\n",GainSigma[i]);
  }

  return AllNum;
}

int E14Calibrator::GetAllNum(){

  if( DataStoreTag == false ){
    printf("Data is not stored yet... aborting...\n");
    return -1;
  }

  return AllNum;
}

void E14Calibrator::GetDataFromText( char *InFileName ){

  int detID = 9999;
  
  initializeDataValues();

  FILE *fp;
  fp = fopen( InFileName ,"rt");
  if( fp == NULL ){
    printf("Calibrtaion file %s is not found.\n", InFileName);
    return;
  }
  printf("Calibration file is open as %s.\n",InFileName);
  
  fscanf(fp,"%d\n",&AllNum);  
  for(int i=0;i<AllNum;i++){
    fscanf(fp,"%d",&detID);
    fscanf(fp,"%f %f\n", &GainMean[detID], &GainSigma[detID] );    
  }
  fclose( fp );
  
  DataStoreTag = true;
}

// Not implement yet
void E14Calibrator::GetDataFromDB( ){

  if( ThisID == 9999 ){
    printf("Set ID for calibration file. e.g. = 0 for CsI\n");
    return;
  }
  
  DataStoreTag = true;
}

void E14Calibrator::ShowGainInfo( int e14id ){
  
  if( e14id < 0 || AllNum < e14id ){
    printf("e14id should be less than %d\n",e14id);
  }

  printf("e14id = %d : GainMean = %f : GainSigma = %f\n",
	 e14id, GainMean[e14id], GainSigma[e14id]);
}

float E14Calibrator::CalcEnergy( int n, float ADCSum ){

  if( GainMean[n] == 9999 ){
    printf("Gain is not read for e14id=%d. Return w/o calibraion.\n",n);
    return (float)( ADCSum );
  }

  // ??? 14.07, 2334.63 from 2010 Autumn run file???
  return (float)( ADCSum * GainMean[n] * 14.07 / 2334.63 );

}
