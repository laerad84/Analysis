/*
 * ConvertDigiData.cc
 *
 *  Created on: Jun 16, 2011
 *      Author: jwlee
 */

/*
 * ANAconvsimCsI.cc
 *  mArr[detectorIndex]                           = this->m_eventData->digi;

 *  Created on: Apr 1, 2011
 *      Author: jwlee
 */
#ifndef CONVERTDIGIDATA_H_
#include "ConvertDigiData.h"
#endif
#include <iostream>
ConvertDigiData::ConvertDigiData(TTree* itr, TTree* otr, char* name)
{
	strcpy(this->m_DetectorName,name);
	std::cout << m_DetectorName << std::endl;
	this->m_InputTree    = itr;
	this->m_OutputTree   = otr;
	this->m_DetectorData = new GsimDetectorData();
	this->m_EventData    = new GsimDetectorEventData();
	this->SetBranchAddress();
	this->Branch();

}

ConvertDigiData::~ConvertDigiData(){
}


Bool_t ConvertDigiData::Branch(){
  //mOutputTree->Branch(Form("%seventNum",detectorName),&eventNum,Form("%seventNum/I",detectorName));
  //Before Change
  /*
  m_OutputTree->Branch(Form("%s.nDigi"   ,m_DetectorName),&m_nDigi ,
		       Form("%s.nDigi/I" ,m_DetectorName));
  m_OutputTree->Branch(Form("%s.totalE"  ,m_DetectorName),&m_TotalE,
		       Form("%s.totalE/F",m_DetectorName));
  m_OutputTree->Branch(Form("%s.ID"    ,m_DetectorName),m_ID     ,
		       Form("%s.ID[%s.nDigi]/s",m_DetectorName,m_DetectorName));//m_nDigi
  m_OutputTree->Branch(Form("%s.energy",m_DetectorName),m_energy ,
		       Form("%s.energy[%s.nDigi]/F",m_DetectorName,m_DetectorName));//m_nDigi
  m_OutputTree->Branch(Form("%s.time"  ,m_DetectorName),m_time   ,
		       Form("%s.time[%s.nDigi]/F"  ,m_DetectorName,m_DetectorName));//m_nDigi
  */
  //After Change
  if( strcmp( "CSI", m_DetectorName )==0){
    m_OutputTree->Branch(Form("%sNumber"   ,"Csi"),&m_nDigi ,
			 Form("%sNumber/I" ,"Csi"));
    m_OutputTree->Branch(Form("%sTotalE"  ,"Csi"),&m_TotalE,
			 Form("%sTotalE/D","Csi"));
    m_OutputTree->Branch(Form("%sModID"    ,"Csi"),m_ID     ,
			 Form("%sModID[%sNumber]/I","Csi","Csi"));//m_nDigi
    m_OutputTree->Branch(Form("%sEne","Csi"),m_energy ,
			 Form("%sEne[%sNumber]/D","Csi","Csi"));//m_nDigi
    m_OutputTree->Branch(Form("%sTime"  ,"Csi"),m_time   ,
			 Form("%sTime[%sNumber]/D"  ,"Csi","Csi"));//m_nDigi
  }else{
    m_OutputTree->Branch(Form("%sNumber"   ,m_DetectorName),&m_nDigi ,
			 Form("%sNumber/I" ,m_DetectorName));
    m_OutputTree->Branch(Form("%sTotalE"  ,m_DetectorName),&m_TotalE,
			 Form("%sTotalE/D",m_DetectorName));
    m_OutputTree->Branch(Form("%sModID"    ,m_DetectorName),m_ID     ,
			 Form("%sModID[%sNumber]/I",m_DetectorName,m_DetectorName));//m_nDigi
    m_OutputTree->Branch(Form("%sEne",m_DetectorName),m_energy ,
			 Form("%sEne[%sNumber]/D",m_DetectorName,m_DetectorName));//m_nDigi
    m_OutputTree->Branch(Form("%sTime"  ,m_DetectorName),m_time   ,
			 Form("%sTime[%sNumber]/D"  ,m_DetectorName,m_DetectorName));//m_nDigi
  }
  return true;
}

Bool_t ConvertDigiData::SetBranchAddress(){
  this->m_InputTree->SetBranchStatus(Form("%s.*"          ,m_DetectorName),0);
  this->m_InputTree->SetBranchStatus(Form("%s.nDigi"      ,m_DetectorName));
  this->m_InputTree->SetBranchStatus(Form("%s.digi.modID" ,m_DetectorName));
  this->m_InputTree->SetBranchStatus(Form("%s.digi.energy",m_DetectorName));
  this->m_InputTree->SetBranchStatus(Form("%s.digi.time"  ,m_DetectorName));
  std::cout << this->m_InputTree->SetBranchAddress(Form("%s.",m_DetectorName),&(this->m_EventData)) << std::endl;
  return true;
}

Bool_t ConvertDigiData::Convert(){
  Int_t m_nDigi  = this->m_EventData->nDigi;
  this->m_nDigi  = this->m_EventData->nDigi;
  m_Arr          = this->m_EventData->digi;
  this->m_TotalE = 0;

  //std::cout<< m_nDigi << std::endl;
  for( Int_t digiIndex = 0; digiIndex < m_nDigi; digiIndex++){
    GsimDigiData* digi        = (GsimDigiData*)m_Arr->UncheckedAt(digiIndex);	    
    this->m_ID[digiIndex]     = (int)digi->modID;
    this->m_time[digiIndex]   = (double)digi->time;
    this->m_energy[digiIndex] = (double)digi->energy;
    this->m_TotalE           += digi->energy;
  }
  return true;
}
