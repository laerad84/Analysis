#include "E14ConvWriterModule.h"

E14ConvWriterModule::E14ConvWriterModule( TTree* tr, char* name ){
  strcpy(this->m_DetectorName, name );
  std::cout<< m_DetectorName << std::endl; 
  this->m_Tree = tr;
  //this->Branch();  
}

E14ConvWriterModule::~E14ConvWriterModule(){
  ;  
}

bool E14ConvWriterModule::InitData(){
  for( int i = 0; i< 4096; i++){
    m_ID[i]            = -1;
    m_Pedestal[i]      = 0.;
    m_Signal[i]        = 0.;
    m_Time[i]          = -9999;
    m_FitTime[i]       = -9999;    
    m_HHTime[i]        = -9999;
    
    m_ParA[i]          = -9999;
    m_ParB[i]          = -9999;
    
    m_FitTime[i]       = -9999;
    m_FitShape[i]      = 0.;
    m_DeltaDiff[i]     = 0.;    
    m_ADC[i]           = 0.;
    m_FitADC[i]        = 0.;
  }
  this->m_nDigi = 0;

  return true;
}

bool E14ConvWriterModule::SetBranchAddress(){

  m_Tree->SetBranchAddress(Form("%sNumber"   ,m_DetectorName),&m_nDigi);
  m_Tree->SetBranchAddress(Form("%sID"       ,m_DetectorName),m_ID);

  m_Tree->SetBranchAddress(Form("%sPedestal" ,m_DetectorName),m_Pedestal);
  m_Tree->SetBranchAddress(Form("%sSignal"   ,m_DetectorName),m_Signal);
  m_Tree->SetBranchAddress(Form("%sTime"     ,m_DetectorName),m_Time);
  m_Tree->SetBranchAddress(Form("%sHHTime"   ,m_DetectorName),m_HHTime);
  m_Tree->SetBranchAddress(Form("%sParA"     ,m_DetectorName),m_ParA);
  m_Tree->SetBranchAddress(Form("%sParB"     ,m_DetectorName),m_ParB);
  m_Tree->SetBranchAddress(Form("%sADC"      ,m_DetectorName),m_ADC);
  m_Tree->SetBranchAddress(Form("%sFitADC"   ,m_DetectorName),m_FitADC);
  m_Tree->SetBranchAddress(Form("%sFitHeight",m_DetectorName),m_FitHeight);
  m_Tree->SetBranchAddress(Form("%sFitTime"  ,m_DetectorName),m_FitTime);
  m_Tree->SetBranchAddress(Form("%sFitShape" ,m_DetectorName),m_FitShape);
  m_Tree->SetBranchAddress(Form("%sDeltaDiff",m_DetectorName),m_DeltaDiff);
  
  m_Tree->SetBranchAddress(Form("%sChisq"    ,m_DetectorName),m_Chisq);
  m_Tree->SetBranchAddress(Form("%sNDF"      ,m_DetectorName),m_NDF);

  return kTRUE;  
}

bool E14ConvWriterModule::Branch(){
  m_Tree->Branch(Form("%sNumber"  ,m_DetectorName),&m_nDigi  ,
		 Form("%sNumber/I"             ,m_DetectorName));
  m_Tree->Branch(Form("%sID"      ,m_DetectorName) ,m_ID      ,
		 Form("%sID[%sNumber]/S"       ,m_DetectorName,m_DetectorName));//m_nDigi

  ///// Branches of Fit Quality /////
  m_Tree->Branch(Form("%sFitHeight",m_DetectorName),m_FitHeight,
		 Form("%sFitHeight[%sNumber]/S",m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFitTime",m_DetectorName),m_FitTime   ,
		 Form("%sFitTime[%sNumber]/S"  ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sChisq",m_DetectorName),m_Chisq       ,
		 Form("%sChisq[%sNumber]/D"    ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sNDF",m_DetectorName),m_NDF,
		 Form("%sNDF[%sNumber]/S"    ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFitShape",m_DetectorName),m_FitShape,
		 Form("%sFitShape[%sNumber]/D",m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sDeltaDiff",m_DetectorName),m_DeltaDiff,
		 Form("%sDeltaDiff[%sNumber]/D",m_DetectorName,m_DetectorName));//m_nDigi
  ///// Raw Fit Result /////
  m_Tree->Branch(Form("%sPedestal",m_DetectorName) ,m_Pedestal,
		 Form("%sPed[%sNumber]/D"      ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sSignal"  ,m_DetectorName) ,m_Signal  ,
		 Form("%sSignal[%sNumber]/D"   ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sTime"  ,m_DetectorName) ,m_Time  ,
		 Form("%sTime[%sNumber]/D"   ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sHHTime",m_DetectorName) ,m_HHTime,
		 Form("%sHHTime[%sNumber]/D" ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sParA"    ,m_DetectorName) ,m_ParA    ,
		 Form("%sParA[%sNumber]/D"     ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sParB"    ,m_DetectorName) ,m_ParB    ,
		 Form("%sParB[%sNumber]/D"     ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sADC"     ,m_DetectorName) ,m_ADC     ,
		 Form("%sADC[%sNumber]/D"      ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFitADC"  ,m_DetectorName) ,m_FitADC  ,
		 Form("%sFitADC[%sNumber]/D"   ,m_DetectorName,m_DetectorName));//m_nDigi
  return kTRUE;
}

bool E14ConvWriterModule::Close(){
  m_Tree->Write();
  return true;
}
