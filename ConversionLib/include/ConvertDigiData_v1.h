/*
 * ConvertDigiData.h
 *
 *  Created on: Jun 16, 2011
 *      Author: jwlee
 */

#ifndef CONVERTDIGIDATA_V1_H_
#define CONVERTDIGIDATA_V1_H_

#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <GsimData/GsimDetectorEventData.h>
#include <GsimData/GsimDetectorData.h>
#include <GsimData/GsimDigiData.h>

class ConvertDigiData_v1
{
private:
	static const Int_t      s_MaxCH = 3000; 
	GsimDetectorData*      	m_DetectorData;
	char                    m_DetectorName[128];

public:
	TClonesArray*           m_Arr;
	GsimDetectorEventData*  m_EventData;
	Int_t                   m_nDigi;
	Float_t                m_TotalE;
	Int_t                   m_ID[s_MaxCH];//m_nDigi
	Float_t                m_energy[s_MaxCH];//m_nDigi
	Float_t                m_time[s_MaxCH];//m_nDigi
	Float_t                m_VetoEne;
	Float_t                m_VetoTime;
	//Int_t    eventNum;

	ConvertDigiData_v1(char* name);
	~ConvertDigiData_v1();

	Bool_t SetBranchAddress(TTree* tr);
	Bool_t Branch(TTree* tr);
	Bool_t setBranchOfDigi(TTree* tr);
	Bool_t Convert();
};


#endif /* CONVERTDIGIDATA_H_ */
