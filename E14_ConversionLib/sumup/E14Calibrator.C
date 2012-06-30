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
    //    printf("%d %f %f\n", detID, GainMean[detID], GainSigma[detID] );    
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

double E14Calibrator::CalcEnergy( int ich, double IntegratedADC, short PeakHeight ){

  if( GainMean[ich] == 9999 ){
    printf("Gain is not read for e14id=%d. Return w/o calibraion.\n",ich);
    return (double)( IntegratedADC );
  }
  
  double energy = 0;

  energy = IntegratedADC / GetLinearityFactor( PeakHeight );

  /*
  // ??? 14.07, 2334.63 from 2010 Autumn run file???
  //  return (double)( ADCSum * GainMean[n] * 14.07 / 2334.63 );
  if( ich < 2240 ){ // Small CsI
    energy = (double)( IntegratedADC / GainMean[ich] * 14.07 );
  }else{  // Large CsI
    energy = (double)( IntegratedADC / GainMean[ich] * 14.07 * 2 );    
  }
  */

  // 2012 run, 14.07 is used for large too.
  energy = (double)( energy / GainMean[ich] * 14.07 );

  return energy;

}

double E14Calibrator::GetLinearityFactor(double peak){
  int mode = 0;
  static TSpline3 *spl[3]={0};
  static double xmin[3]={0};
  static double xmax[3]={0};
  static double ymax[3]={0};
  if(spl[0]==0){
    std::ifstream ifs("/share/apps/production/2012VME/sumup/calib/getLinearityFunction.dat");
    if(!ifs){
      std::cout<<"error: can't find linearity data"<<std::endl;
    }
    for(int imode=0;imode<3;imode++){
      int thismode,num;
      ifs>>thismode>>num;
      //      std::cout<<"linearity function::  mode="<<thismode<<" npoint="<<num<std::endl;    std::cout<<"linearity function::  mode="<<thismode<<" npoint="<<num<      
      double xspl[1000]={0},yspl[1000]={0};  
      for(int i=0;i<num;i++){ 
        ifs>>xspl[i]>>yspl[i]; 
      }
      xmin[imode] = xspl[0];
      xmax[imode] = xspl[num-1];  
      ymax[imode] = yspl[num-1];   
      //      std::cout<<"range:"<<xmin[imode]<<" to  "<<xmax[imode]<<std::endl;
      spl[imode] = new TSpline3(Form("spl%d",imode),xspl,yspl,num); 
    }
    ifs.close();
  }
  if(peak<xmin[mode]) return 1;
  else if(peak>=xmax[mode]) return ymax[mode];
  else return spl[mode]->Eval(peak);
}
