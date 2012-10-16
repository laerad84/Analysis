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
  m_nTimeCluster = 0;
  for( int i  =0; i< 32 ; i++){
    m_TimeClusterHead[i] = 0.;
    m_TimeClusterTail[i] = 0.;
    m_TotalEnergyInTimeCluster[i] = 0.;
  }

  m_TotalEnergy = 0.;
  m_nOverFlow   = 0;
  m_nUnderFlow  = 0;  
  for( int i = 0; i< 2716; i++){
    m_IDOverFlow[i]  = 0;
    m_IDUnderFlow[i] = 0;

    m_ID[i]            = -1;
    m_TimeClusterID[i] = 0;
    m_Pedestal[i]      = 0.;
    m_Signal[i]        = 0.;
    m_Time[i]          = 0.;
    m_FitTime[i]       = 0.;    
    m_HHTime[i]        = 0.;
    m_Energy[i]        = 0.;
    
    m_ParA[i]          = 0.;
    m_ParB[i]          = 0.;
    
    m_FitTime[i]       = 0.;
    m_FitShape[i]      = 0.;
    m_DeltaDiff[i]     = 0.;    
    m_ADC[i]           = 0.;
    m_FitADC[i]        = 0.;

    m_wav_SlopeDelta[i] = 0.;
    m_wav_Height[i]     = 0.;
    m_wav_PeakTime[i]   = 0.;
    m_wav_Pedestal[i]   = 0.;
    m_wav_Width[i]      = 0.;
    m_wav_FrontHalfTime[i] = 0.;
    m_wav_RearHalfTime[i]  = 0.;

    m_Fit_Pedestal[i]  = 0.;
    m_Fit_Time[i]      = 0.;
    m_Fit_Height[i]    = 0.;
    m_Fit_HHTime[i]    = 0.;
    m_Fit_ChisqNDF[i]  = 0.;
    m_Fit_ChisqPed[i]  = 0.;
    m_Fit_ChisqFront[i]= 0.;
    m_Fit_ChisqRear[i] = 0.;
    m_Fit_ChisqTail[i] = 0.;
    
    m_Conv_Energy[i]   = 0.;
    m_Conv_Time[i]     = 0.;
  }
  this->m_nDigi = 0;

  return true;
}

bool E14ConvWriterModule::SetBranchAddress(){

  m_Tree->SetBranchAddress(Form("%snTimeCluster",m_DetectorName),&m_nTimeCluster);
  m_Tree->SetBranchAddress(Form("%sTimeClusterHead",m_DetectorName),m_TimeClusterHead);
  m_Tree->SetBranchAddress(Form("%sTimeClusterTail",m_DetectorName),m_TimeClusterTail);
  m_Tree->SetBranchAddress(Form("%sTotalEnergyInTimeCluster",m_DetectorName),m_TotalEnergyInTimeCluster);

  m_Tree->SetBranchAddress(Form("%sTotalEnergy",m_DetectorName),&m_TotalEnergy);
  // Under/OverFlow 
  m_Tree->SetBranchAddress(Form("%snOverFlow",m_DetectorName),&m_nOverFlow);
  m_Tree->SetBranchAddress(Form("%snUnderFlow",m_DetectorName),&m_nUnderFlow);
  m_Tree->SetBranchAddress(Form("%sIDOverFlow",m_DetectorName),m_IDOverFlow);
  m_Tree->SetBranchAddress(Form("%sIDUnderFlow",m_DetectorName),m_IDUnderFlow);			   

  // Values of waveform // 
  m_Tree->SetBranchAddress(Form("%sNumber"   ,m_DetectorName),&m_nDigi);
  m_Tree->SetBranchAddress(Form("%sID"       ,m_DetectorName),m_ID);
  m_Tree->SetBranchAddress(Form("%sTimeClusterID",m_DetectorName),m_TimeClusterID);

  m_Tree->SetBranchAddress(Form("%sPedestal" ,m_DetectorName),m_Pedestal);
  m_Tree->SetBranchAddress(Form("%sSignal"   ,m_DetectorName),m_Signal);
  m_Tree->SetBranchAddress(Form("%sTime"     ,m_DetectorName),m_Time);
  m_Tree->SetBranchAddress(Form("%sHHTime"   ,m_DetectorName),m_HHTime);
  m_Tree->SetBranchAddress(Form("%sEne"      ,m_DetectorName),m_Energy);
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
  
  if(strcmp( m_DetectorName, "Csi" ) == 0){
    m_Tree->SetBranchAddress(Form("%swav_SlopeDelta"   ,m_DetectorName),m_wav_SlopeDelta);
    m_Tree->SetBranchAddress(Form("%swav_Height"       ,m_DetectorName),m_wav_Height);
    m_Tree->SetBranchAddress(Form("%swav_PeakTime"     ,m_DetectorName),m_wav_PeakTime);
    m_Tree->SetBranchAddress(Form("%swav_Pedestal"     ,m_DetectorName),m_wav_Pedestal);
    m_Tree->SetBranchAddress(Form("%swav_Width"        ,m_DetectorName),m_wav_Width);
    m_Tree->SetBranchAddress(Form("%swav_FrontHalfTime",m_DetectorName),m_wav_FrontHalfTime);
    m_Tree->SetBranchAddress(Form("%swav_RearHalfTime" ,m_DetectorName),m_wav_RearHalfTime);
    
    m_Tree->SetBranchAddress(Form("%sFit_Pedestal"  ,m_DetectorName),m_Fit_Pedestal);
    m_Tree->SetBranchAddress(Form("%sFit_Time"      ,m_DetectorName),m_Fit_Time);
    m_Tree->SetBranchAddress(Form("%sFit_Height"    ,m_DetectorName),m_Fit_Height);
    m_Tree->SetBranchAddress(Form("%sFit_HHTime"    ,m_DetectorName),m_Fit_HHTime);
    m_Tree->SetBranchAddress(Form("%sFit_ChisqNDF"  ,m_DetectorName),m_Fit_ChisqNDF);
    m_Tree->SetBranchAddress(Form("%sFit_ChisqPed"  ,m_DetectorName),m_Fit_ChisqPed);
    m_Tree->SetBranchAddress(Form("%sFit_ChisqFront",m_DetectorName),m_Fit_ChisqFront);
    m_Tree->SetBranchAddress(Form("%sFit_ChisqRear" ,m_DetectorName),m_Fit_ChisqRear);
    m_Tree->SetBranchAddress(Form("%sFit_ChisqTail" ,m_DetectorName),m_Fit_ChisqTail);
    
    m_Tree->SetBranchAddress(Form("%sConv_Energy",m_DetectorName),m_Conv_Energy);
    m_Tree->SetBranchAddress(Form("%sConv_Time"  ,m_DetectorName),m_Conv_Time);
  }

  return kTRUE;  
}

bool E14ConvWriterModule::Branch(){
  m_Tree->Branch(Form("%snTimeCluster"  ,m_DetectorName),&m_nTimeCluster,
		 Form("%snTimeCluster/I",m_DetectorName));
  m_Tree->Branch(Form("%sTimeClusterHead"                           ,m_DetectorName),m_TimeClusterHead,
		 Form("%sTimeClusterHead[%snTimeCluster]/D"         ,m_DetectorName,m_DetectorName));
  m_Tree->Branch(Form("%sTimeClusterTail"                           ,m_DetectorName),m_TimeClusterTail,
		 Form("%sTimeClusterTail[%snTimeCluster]/D"         ,m_DetectorName,m_DetectorName));
  m_Tree->Branch(Form("%sTotalEnergyInTimeCluster"                  ,m_DetectorName),m_TotalEnergyInTimeCluster,
		 Form("%sTotalEnergyInTimeCluster[%snTimeCluster]/I",m_DetectorName,m_DetectorName));

  m_Tree->Branch(Form("%sTotalEnergy",m_DetectorName),&m_TotalEnergy,
		 Form("%sTotalEnergy/D",m_DetectorName));

  m_Tree->Branch(Form("%snOverFlow"                    ,m_DetectorName),&m_nOverFlow,
		 Form("%snOverFlow/I"                  ,m_DetectorName));
  m_Tree->Branch(Form("%sIDOverFlow"                   ,m_DetectorName),m_IDOverFlow,
		 Form("%sIDOverFlow[%snOverFlow]/I"    ,m_DetectorName,m_DetectorName));//m_nOverFlow
  m_Tree->Branch(Form("%snUnderFlow"                   ,m_DetectorName),&m_nUnderFlow,
		 Form("%snUnderFlow/I"                 ,m_DetectorName));
  m_Tree->Branch(Form("%sIDUnderFlow"                  ,m_DetectorName),m_IDUnderFlow,
		 Form("%sIDUnderFlow[%snUnderFlow]/I"  ,m_DetectorName,m_DetectorName));//m_nUnderFlow

  m_Tree->Branch(Form("%sNumber"                       ,m_DetectorName),&m_nDigi,
		 Form("%sNumber/I"                     ,m_DetectorName));
  m_Tree->Branch(Form("%sID"                           ,m_DetectorName),m_ID,
		 Form("%sID[%sNumber]/S"               ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sTimeClusterID"                ,m_DetectorName),m_TimeClusterID,
		 Form("%sTimeClusterID[%sNumber]/I"    ,m_DetectorName,m_DetectorName));//m_nDigi
  ///// Branches of Fit Quality /////
  m_Tree->Branch(Form("%sFitHeight"                    ,m_DetectorName),m_FitHeight,
		 Form("%sFitHeight[%sNumber]/S"        ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFitTime"                      ,m_DetectorName),m_FitTime,
		 Form("%sFitTime[%sNumber]/S"          ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sChisq"                        ,m_DetectorName),m_Chisq,
		 Form("%sChisq[%sNumber]/D"            ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sNDF"                          ,m_DetectorName),m_NDF,
		 Form("%sNDF[%sNumber]/S"              ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFitShape"                     ,m_DetectorName),m_FitShape,
		 Form("%sFitShape[%sNumber]/D"         ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sDeltaDiff"                    ,m_DetectorName),m_DeltaDiff,
		 Form("%sDeltaDiff[%sNumber]/D"        ,m_DetectorName,m_DetectorName));//m_nDigi
  ///// Raw Fit Result /////
  m_Tree->Branch(Form("%sPedestal"                     ,m_DetectorName) ,m_Pedestal,
		 Form("%sPedestal[%sNumber]/D"         ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sSignal"                       ,m_DetectorName) ,m_Signal,
		 Form("%sSignal[%sNumber]/D"           ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sEne"                          ,m_DetectorName)    ,m_Energy,
		 Form("%sEne[%sNumber]/D"              ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sTime"                         ,m_DetectorName) ,m_Time,
		 Form("%sTime[%sNumber]/D"             ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sHHTime"                       ,m_DetectorName) ,m_HHTime,
		 Form("%sHHTime[%sNumber]/D"           ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sParA"                         ,m_DetectorName) ,m_ParA,
		 Form("%sParA[%sNumber]/D"             ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sParB"                         ,m_DetectorName) ,m_ParB,
		 Form("%sParB[%sNumber]/D"             ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sADC"                          ,m_DetectorName) ,m_ADC,
		 Form("%sADC[%sNumber]/D"              ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFitADC"                       ,m_DetectorName)     ,m_FitADC,
		 Form("%sFitADC[%sNumber]/D"           ,m_DetectorName,m_DetectorName));//m_nDigi

  if( strcmp( m_DetectorName, "Csi" ) == 0 ){
  m_Tree->Branch(Form("%swav_SlopeDelta"               ,m_DetectorName),m_wav_SlopeDelta,
		 Form("%swav_SlopeDelta[%sNumber]/D"   ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%swav_Height"                   ,m_DetectorName),m_wav_Height,
		 Form("%swav_Height[%sNumber]/D"       ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%swav_PeakTime"                 ,m_DetectorName),m_wav_PeakTime,
		 Form("%swav_PeakTime[%sNumber]/D"     ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%swav_Pedestal"                 ,m_DetectorName),m_wav_Pedestal,
		 Form("%swav_Pedestal[%sNumber]/D"     ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%swav_Width"                    ,m_DetectorName),m_wav_Width,
		 Form("%swav_Width[%sNumber]/D"        ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%swav_FrontHalfTime"            ,m_DetectorName),m_wav_FrontHalfTime,
		 Form("%swav_FrontHalfTime[%sNumber]/D",m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%swav_RearHalfTime"             ,m_DetectorName),m_wav_RearHalfTime,
		 Form("%swav_RearHalfTime[%sNumber]/D" ,m_DetectorName,m_DetectorName));//m_nDigi

  m_Tree->Branch(Form("%sFit_Pedestal"                 ,m_DetectorName),m_Fit_Pedestal,
		 Form("%sFit_Pedestal[%sNumber]/D"     ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit_Time"                     ,m_DetectorName),m_Fit_Time,
		 Form("%sFit_Time[%sNumber]/D"         ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit_Height"                   ,m_DetectorName),m_Fit_Height,
		 Form("%sFit_Height[%sNumber]/D"       ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit_HHTime"                   ,m_DetectorName),m_Fit_HHTime,
		 Form("%sFit_HHTime[%sNumber]/D"       ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit_ChisqNDF"                 ,m_DetectorName),m_Fit_ChisqNDF,
		 Form("%sFit_ChisqNDF[%sNumber]/D"     ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit_ChisqPed"                 ,m_DetectorName),m_Fit_ChisqPed,
		 Form("%sFit_ChisqPed[%sNumber]/D"     ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit_ChisqFront"               ,m_DetectorName),m_Fit_ChisqFront,
		 Form("%sFit_ChisqFront[%sNumber]/D"   ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit_ChisqRear"                ,m_DetectorName),m_Fit_ChisqRear,
		 Form("%sFit_ChisqRear[%sNumber]/D"    ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sFit_ChisqTail"                ,m_DetectorName),m_Fit_ChisqTail,
		 Form("%sFit_ChisqTail[%sNumber]/D"    ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sConv_Energy"                  ,m_DetectorName),m_Conv_Energy,
		 Form("%sConv_Energy[%sNumber]/D"      ,m_DetectorName,m_DetectorName));//m_nDigi
  m_Tree->Branch(Form("%sConv_Time"                    ,m_DetectorName),m_Conv_Time,
		 Form("%sConv_Time[%sNumber]/D"        ,m_DetectorName,m_DetectorName));//m_nDigi
  }
  return kTRUE;
}

bool E14ConvWriterModule::Close(){
  m_Tree->Write();
  return true;
}
