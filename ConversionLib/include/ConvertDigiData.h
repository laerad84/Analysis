/*
 * ConvertDigiData.h
 *
 *  Created on: Jun 16, 2011
 *      Author: jwlee
 */

#ifndef CONVERTDIGIDATA_H_
#define CONVERTDIGIDATA_H_

#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <GsimData/GsimDetectorEventData.h>
#include <GsimData/GsimDetectorData.h>
#include <GsimData/GsimDigiData.h>

class ConvertDigiData
{
private:
	static const Int_t      s_MaxCH = 3000;
	TTree*                  m_InputTree;
	TTree*                  m_OutputTree;
 
	GsimDetectorData*      	m_DetectorData;
	char                    m_DetectorName[128];

public:
	TClonesArray*           m_Arr;
	GsimDetectorEventData*  m_EventData;
	Int_t                   m_nDigi;
	Double_t                m_TotalE;
	Int_t                   m_ID[s_MaxCH];//m_nDigi
	Double_t                m_energy[s_MaxCH];//m_nDigi
	Double_t                m_time[s_MaxCH];//m_nDigi

	//Int_t    eventNum;

	ConvertDigiData( TTree* , TTree* ,char*);
	~ConvertDigiData();

	Bool_t SetBranchAddress();
	Bool_t Branch();
	Bool_t Convert();
};


#endif /* CONVERTDIGIDATA_H_ */
