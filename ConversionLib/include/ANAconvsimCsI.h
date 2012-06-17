/*
 * ANAconvsimCsI.h
 *
 *  Created on: Apr 1, 2011
 *      Author: jwlee
 */

#ifndef ANACONVSIMCSI_H_
#define ANACONVSIMCSI_H_

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <GsimData/GsimDetectorEventData.h>
#include <GsimData/GsimDetectorData.h>
#include <GsimData/GsimDigiData.h>

class ANAconvSim
{
private:

	static const Int_t nTotalCH = 3000;
	TFile* 					m_itfile;
	TFile* 					m_otfile;
	TTree*   			   	m_eventTree;
	TTree*                 	m_outTree;
	GsimDetectorData*      	m_detectorData;
	GsimDetectorEventData* 	m_eventData;

	Int_t    eventNum;
	Int_t    nDigi;
	Double_t totalE;
	Int_t    ID[nTotalCH];//nDigi
	Double_t depE[nTotalCH];//nDigi

public:
	ANAconvSim(char*, char* ,char*);
	~ANAconvSim();
	Bool_t SetBranchAddress();
	Bool_t Branch();
	Long_t Convert();
	Bool_t SetInputFile(char*);
	Bool_t SetOutputFile(char*,char*);
	Bool_t Clear();
};

#endif /* ANACONVSIMCSI_H_ */
