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
#ifndef CONVERTDIGIDATA_V1_H_
#include "ConvertDigiData_v1.h"
#endif
#include <iostream>
ConvertDigiData_v1::ConvertDigiData_v1(char* name)
{
	strcpy(this->m_DetectorName,name);
	std::cout << m_DetectorName << std::endl;
	this->m_DetectorData = new GsimDetectorData();
	this->m_EventData    = new GsimDetectorEventData();

}

ConvertDigiData_v1::~ConvertDigiData_v1(){
}


Bool_t ConvertDigiData_v1::Branch(TTree* tree){

  tree->Branch(Form("%sModuleNumber"   ,m_DetectorName),&m_nDigi ,
	       Form("%sModuleNumber/I" ,m_DetectorName));
  tree->Branch(Form("%sTotalVetoEne"  ,m_DetectorName),&m_TotalE,
	       Form("%sTotalVetoEne/F",m_DetectorName));
  tree->Branch(Form("%sModuleModID"    ,m_DetectorName),m_ID     ,
	       Form("%sModuleModID[%sNumber]/I",m_DetectorName,m_DetectorName));//m_nDigi
  tree->Branch(Form("%sModuleEne",m_DetectorName),m_energy ,
	       Form("%sModuleEne[%sNumber]/F",m_DetectorName,m_DetectorName));//m_nDigi
  tree->Branch(Form("%sModuleTime"  ,m_DetectorName),m_time   ,
	       Form("%sModuleTime[%sNumber]/F"  ,m_DetectorName,m_DetectorName));//m_nDigi
  tree->Branch(Form("%sVetoEne"       ,m_DetectorName),&m_VetoEne,Form("%sVetoEne/F",m_DetectorName));
  tree->Branch(Form("%sVetoTime"      ,m_DetectorName),&m_VetoTime,Form("%sVetoTime/F",m_DetectorName));

  return true;
}
Bool_t ConvertDigiData_v1::setBranchOfDigi(TTree* tree){

  tree->SetBranchAddress(Form("%sModuleNumber"  ,m_DetectorName),&m_nDigi );
  tree->SetBranchAddress(Form("%sTotalVetoEne"  ,m_DetectorName),&m_TotalE);
  tree->SetBranchAddress(Form("%sModuleModID"   ,m_DetectorName),m_ID     );//m_nDigi
  tree->SetBranchAddress(Form("%sModuleEne"     ,m_DetectorName),m_energy );//m_nDigi
  tree->SetBranchAddress(Form("%sModuleTime"    ,m_DetectorName),m_time   );//m_nDigi
  tree->SetBranchAddress(Form("%sVetoEne"       ,m_DetectorName),&m_VetoEne);
  tree->SetBranchAddress(Form("%sVetoTime"      ,m_DetectorName),&m_VetoTime);
  return true;
}

Bool_t ConvertDigiData_v1::SetBranchAddress(TTree* tree){
  tree->SetBranchStatus(Form("%s.*"          ,m_DetectorName),0);
  tree->SetBranchStatus(Form("%s.nDigi"      ,m_DetectorName));
  tree->SetBranchStatus(Form("%s.digi.modID" ,m_DetectorName));
  tree->SetBranchStatus(Form("%s.digi.energy",m_DetectorName));
  tree->SetBranchStatus(Form("%s.digi.time"  ,m_DetectorName));
  std::cout << tree->SetBranchAddress(Form("%s.",m_DetectorName),&(this->m_EventData)) << std::endl;
  return true;
}

Bool_t ConvertDigiData_v1::Convert(){
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
  m_VetoTime = 0;
  m_VetoEne  = 0;
  for( int digiIndex  = 0; digiIndex < m_nDigi; digiIndex++ ){
    if( m_energy[digiIndex] > m_VetoEne ){
      m_VetoEne = m_energy[digiIndex];
      m_VetoTime = m_time[digiIndex];
    }
  }
  


  return true;
}
