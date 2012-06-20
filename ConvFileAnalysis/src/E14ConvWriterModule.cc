#include "E14ConvWriterModule.h"

E14ConvWriterModule::E14ConvWriterModule( TTree* otr, char* name ){
  strcpy(this->m_DetectorName, name );
  std::cout<< m_DetectorName << std::endl; 
  this->m_OutputTree = otr;
  this->Branch();  
}

E14ConvWriterModule::E14ConvWriterModule(){
  ;  
}
Bool_t E14ConvWriterModule::InitData(){
  m_nDigi = 0;
  for( int i = 0; i< 4096; i++){
    m_ID[i]        = -1;
    m_Pedestal[i]  = 0;
    m_Signal[i]    = 0;
    m_Timing[i]    = -9999;
    m_HHTimging[i] = -9999;
    m_ParA[i]      = -9999;
    m_ParB[i]      = -9999;
    m_Fit[i]       = -9999;
  }
  return true;
}

Bool_t E14ConvWriterModule::SetBranchAddress(){
  m_Tree->SetBranchAddress(Form("%sNumber"  ,m_DetectorName),&m_nDigi);
  m_Tree->SetBranchAddress(Form("%sID"      ,m_DetectorName),m_ID);
  m_Tree->SetBranchAddress(Form("%sPedestal",m_DetectorName),m_Pedestal);
  m_Tree->SetBranchAddress(Form("%sSignal"  ,m_DetectorName),m_Signal);
  m_Tree->SetBranchAddress(Form("%sTiming"  ,m_DetectorName),m_Timing);
  m_Tree->SetBranchAddress(Form("%sHHTiming",m_DetectorName),m_HHTiming);
  m_Tree->SetBranchAddress(Form("%sParA"    ,m_DetectorName),m_ParA);
  m_Tree->SetBranchAddress(Form("%sParB"    ,m_DetectorName),m_ParB);
  m_Tree->SetBranchAddress(Form("%sFit"     ,m_DetectorName),m_Fit);
  return kTRUE;  
}

Bool_t E14ConvWriterModule::Branch(){
  m_Tree->Branch(Form("%sNumber"  ,m_DetectorName),&m_nDigi  ,
		 Form("%sNumber/I"            ,m_DetectorName));
  m_Tree->Branch(Form("%sID"      ,m_DetectorName),m_ID      ,
		 Form("%sID[%sNumber]/I"      ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sPedestal",m_DetectorName),m_Pedestal,
		 Form("%sPed[%sNumber]/D"     ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sSignal"  ,m_DetectorName),m_Signal  ,
		 Form("%sSignal[%sNumber]/D"  ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sTiming"  ,m_DetectorName),m_Timing  ,
		 Form("%sTiming[%sNumber]/D"  ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sHHTiming",m_DetectorName),m_HHTiming,
		 Form("%sHHTiming[%sNumber]/D",m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sParA"    ,m_DetectorName),m_ParA    ,
		 Form("%sParA[%sNumber]/D"    ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sParB"    ,m_DetectorName),m_ParB    ,
		 Form("%sParB[%sNumber]/D"    ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit"     ,m_DetectorName),m_Fit     ,
		 Form("%sFit[%sNumber]/I"     ,m_DetectorName,m_DetectorName));//m_nDigi
  return kTRUE;
}


