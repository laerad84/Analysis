/*
 * ConvertData.h
 *
 *  Created on: Jun 17, 2011
 *      Author: jwlee
 */

#ifndef CONVERTDATA_H_
#define CONVERTDATA_H_

#include <iostream>
#include <list>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <GsimData/GsimDetectorEventData.h>
#include <GsimData/GsimDetectorData.h>
#include <GsimData/GsimDigiData.h>

#include "ConvertDigiData.h"
#include "ConvertTrackData.h"
class ConvertData{
private:
	static const Int_t mNmaxDetector = 128;
	ConvertDigiData*   mDetectorDigiData[mNmaxDetector];
	ConvertTrackData*  mTrackData;
	TClonesArray*      mArr[mNmaxDetector];
	TClonesArray       mDetectorArray;
	TClonesArray*      mTrackArr;

	Int_t  mNDetector;
	TFile* mInputFile;
	TFile* mOutputFile;
	TTree* mInputTree;
	TTree* mOutputTree;
	Long_t mEntries;
	Long_t mEventNum;
	Bool_t mfConvertTrack;

public:

	ConvertData(char* ,char* );
	~ConvertData();
	Bool_t Convert(Long_t);
	Long_t ConvertAll();
	Int_t  GetNDetector();
	Int_t  Add(char*);
	Bool_t AddTrack();

};


#endif /* CONVERTDATA_H_ */
