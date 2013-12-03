/*
 * ConvertData.cc
 *
 *  Created on: Jun 17, 2011
 *      Author: jwlee
 */

#ifndef CONVERTDATA_H_
#include "ConvertData.h"
#define CONVERTDATA_H_

ConvertData::ConvertData(char* iFilename, char* oFilename):mfConvertTrack(false)
{
	this->mInputFile  = new TFile(iFilename);
	this->mInputTree  = (TTree*)this->mInputFile->Get("eventTree00");
	this->mOutputFile = new TFile(oFilename,"RECREATE");
	this->mOutputTree = new TTree("T","");
	mNDetector        = 0;
	this->mEntries    = this->mInputTree->GetEntries();
	this->mOutputTree->Branch("EventNum",&(this->mEventNum),"EventNum/I");
	this->mOutputTree->SetAutoSave();
}

ConvertData::~ConvertData()
{
	this->mOutputTree->Write();
	this->mOutputFile->Close();
	//this->mDetectorArray.Delete();
}

Int_t ConvertData::Add(char* detectorName){
	//this->mDetectorDigiData[mNDetector] = new(this->mDetectorArray[mNDetector]) ConvertDigiData(this->mInputTree,this->mOutputTree,detectorName);
	this->mDetectorDigiData[mNDetector] = new ConvertDigiData(this->mInputTree,this->mOutputTree,detectorName);
	mNDetector++;
	return mNDetector;
}

Bool_t ConvertData::AddTrack(){
  if( !mfConvertTrack ){
    this->mTrackData   = new ConvertTrackData(this->mInputTree,this->mOutputTree);
    mfConvertTrack = true;
    return mfConvertTrack;
  }else{
    return false;
  }
}

Int_t ConvertData::GetNDetector(){
	return mNDetector;
}

Bool_t ConvertData::Convert(Long_t eventNum){
  if(eventNum > this->mInputTree->GetEntries()){return false;}
  this->mInputTree->GetEntry(eventNum);
  //std::cout << eventNum << std::endl;
  this->mEventNum = eventNum;  
  /*
    for( Int_t detectorIndex  = 0 ; detectorIndex < this->mNDetector; detectorIndex++){
    Int_t nDigi                                   = this->mDetectorDigiData[detectorIndex]->m_eventData->nDigi;
    this->mDetectorDigiData[detectorIndex]->nDigi = this->mDetectorDigiData[detectorIndex]->m_eventData->nDigi;
    mArr[detectorIndex]                           = this->mDetectorDigiData[detectorIndex]->m_eventData->digi;
    mDetectorDigiData[detectorIndex]->totalE     = 0;
    
    for( Int_t digiIndex = 0; digiIndex < nDigi; digiIndex++){
    GsimDigiData* digi                                = (GsimDigiData*)mArr[detectorIndex]->UncheckedAt(digiIndex);	    
    mDetectorDigiData[detectorIndex]->ID[digiIndex]   = digi->modID;
    mDetectorDigiData[detectorIndex]->depE[digiIndex] = digi->energy;
    mDetectorDigiData[detectorIndex]->totalE         += digi->energy;
    }
    }
    
    if( mConvertTrack ){
    Int_t nTracks            = this->mTrackData->m_trackData->GetSize();
    this->mTrackData->nTrack = this->mTrackData->m_trackData->GetSize();
    for( Int_t trackIndex = 0; trackIndex < nTracks; trackIndex++){
    
    }
    }
  */

  if( this->mNDetector > 0){
    for( Int_t detectorIndex = 0; detectorIndex < this->mNDetector; detectorIndex++){
      this->mDetectorDigiData[detectorIndex]->Convert();
    }
  }
  if( this->mfConvertTrack ){
    this->mTrackData->Convert();
  }
  //std::cout << "Convert End"  << std::endl;
  return true;
}




Long_t ConvertData::ConvertAll(){

	Long_t nEvent = this->mInputTree->GetEntries();
	//nEvent=10000;
	for( Long_t ievent = 0; ievent < nEvent; ievent++){
		if(ievent%1000==0){std::cout << ievent << "/" << nEvent << std::endl;}
		//this->mInputTree->GetEntry(ievent);
		this->Convert(ievent);

		this->mOutputTree->Fill();
		//std::cout << "FillEnd" << std::endl;
	}


	return nEvent;
}
#endif //CONVERTDATA_H_
