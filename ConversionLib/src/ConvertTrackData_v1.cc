/*
 * ConvertTrackData_v1.cc
 *
 *  Created on: Jun 18, 2011
 *      Author: jwlee
 */

#ifndef CONVERTDIGIDATA_V1_H_
#include "ConvertTrackData_v1.h"
#endif
#include <iostream>
ConvertTrackData_v1::ConvertTrackData_v1()
{
	this->m_particleData = new GsimGenParticleData();
	this->m_trackData    = new GsimTrackData();
}

ConvertTrackData_v1::~ConvertTrackData_v1(){
  ;
}

Bool_t ConvertTrackData_v1::Branch(TTree* OutputTree){
  OutputTree->Branch("nTrack",&nTrack,"nTrack/I");
  OutputTree->Branch("track" ,track  ,"track[nTrack]/s");//nTrack
  OutputTree->Branch("mother",mother ,"mother[nTrack]/S");//nTrack
  OutputTree->Branch("pid"   ,pid    ,"pid[nTrack]/I");//nTrack
  OutputTree->Branch("mass"  ,mass   ,"mass[nTrack]/F");//nTrack
  OutputTree->Branch("ek"    ,ek     ,"ek[nTrack]/F");//nTrack
  OutputTree->Branch("end_ek",end_ek ,"end_ek[nTrack]/F");//nTrack
  OutputTree->Branch("p"     ,p      ,"p[nTrack][3]/D");//nTrack
  OutputTree->Branch("end_p" ,end_p  ,"end_p[nTrack][3]/D");//nTrack
  OutputTree->Branch("v"     ,v      ,"v[nTrack][3]/D");//nTrack
  OutputTree->Branch("end_v" ,end_v  ,"end_v[nTrack][3]/D");//nTrack
}

Bool_t ConvertTrackData_v1::setBranchOfTrack(TTree* InputTree){
  InputTree->SetBranchAddress("nTrack",&nTrack);
  InputTree->SetBranchAddress("track" ,track  );//nTrack
  InputTree->SetBranchAddress("mother",mother );//nTrack
  InputTree->SetBranchAddress("pid"   ,pid    );//nTrack
  InputTree->SetBranchAddress("mass"  ,mass   );//nTrack
  InputTree->SetBranchAddress("ek"    ,ek     );//nTrack
  InputTree->SetBranchAddress("end_ek",end_ek );//nTrack
  InputTree->SetBranchAddress("p"     ,p      );//nTrack
  InputTree->SetBranchAddress("end_p" ,end_p  );//nTrack
  InputTree->SetBranchAddress("v"     ,v      );//nTrack
  InputTree->SetBranchAddress("end_v" ,end_v  );//nTrack
}


Bool_t ConvertTrackData_v1::SetBranchAddress(TTree* InputGsimTree){
  std::cout << InputGsimTree->SetBranchAddress("GenParticle.",&(this->m_particleData)) << std::endl;
  return true;
}

Bool_t ConvertTrackData_v1::Convert(){

  this->m_trackArr = (TClonesArray*)(this->m_particleData->briefTracks);
  this->nTrack     = this->m_trackArr->GetEntries();
  //std::cout << "nTrack : " << this->nTrack << std::endl;
  for( Int_t trackIndex = 0; trackIndex < this->nTrack; trackIndex++){
    GsimTrackData* trackData = (GsimTrackData*)m_trackArr->UncheckedAt(trackIndex);
    this->track[trackIndex]  = trackData->track;
    this->mother[trackIndex] = trackData->mother;
    this->pid[trackIndex]    = trackData->pid;
    this->mass[trackIndex]   = trackData->mass;
    this->ek[trackIndex]     = trackData->ek;
    this->end_ek[trackIndex] = trackData->end_ek;
    trackData->p.GetXYZ    (this->p[trackIndex]);
    trackData->v.GetXYZ    (this->v[trackIndex]);
    trackData->end_p.GetXYZ(this->end_p[trackIndex]);
    trackData->end_v.GetXYZ(this->end_v[trackIndex]);
  }
  return true;
}
