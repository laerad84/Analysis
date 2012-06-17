/*
 * ConvertTrackData.cc
 *
 *  Created on: Jun 18, 2011
 *      Author: jwlee
 */

#ifndef CONVERTDIGIDATA_H_
#include "ConvertTrackData.h"
#endif
#include <iostream>
ConvertTrackData::ConvertTrackData(TTree* itr, TTree* otr)
{

	this->mInputTree = itr;
	this->mOutputTree = otr;
	this->m_particleData = new GsimGenParticleData();
	this->m_trackData    = new GsimTrackData();
	this->SetBranchAddress();
	this->Branch();

}

ConvertTrackData::~ConvertTrackData(){
}


Bool_t ConvertTrackData::Branch(){

	mOutputTree->Branch("nTrack",&nTrack,"nTrack/I");
	mOutputTree->Branch("track" ,track  ,"track[nTrack]/s");//nTrack
	mOutputTree->Branch("mother",mother ,"mother[nTrack]/S");//nTrack
	mOutputTree->Branch("pid"   ,pid    ,"pid[nTrack]/I");//nTrack
	mOutputTree->Branch("mass"  ,mass   ,"mass[nTrack]/F");//nTrack
	mOutputTree->Branch("ek"    ,ek     ,"ek[nTrack]/F");//nTrack
	mOutputTree->Branch("end_ek",end_ek ,"end_ek[nTrack]/F");//nTrack
	mOutputTree->Branch("p"     ,p      ,"p[nTrack][3]/D");//nTrack
	mOutputTree->Branch("end_p" ,end_p  ,"end_p[nTrack][3]/D");//nTrack
	mOutputTree->Branch("v"     ,v      ,"v[nTrack][3]/D");//nTrack
	mOutputTree->Branch("end_v" ,end_v  ,"end_v[nTrack][3]/D");//nTrack

	/*
	mOutputTree->Branch("px"    ,px     ,"px[nTrack]/D");//nTrack
	mOutputTree->Branch("py"    ,py     ,"py[nTrack]/D");//nTrack
	mOutputTree->Branch("pz"    ,pz     ,"pz[nTrack]/D");//nTrack
	mOutputTree->Branch("vx"    ,vx     ,"vx[nTrack]/D");//nTrack
	mOutputTree->Branch("vy"    ,vy     ,"vy[nTrack]/D");//nTrack
	mOutputTree->Branch("vz"    ,vz     ,"vz[nTrack]/D");//nTrack
	mOutputTree->Branch("end_px",end_px ,"end_px[nTrack]/D");//nTrack
	mOutputTree->Branch("end_py",end_py ,"end_py[nTrack]/D");//nTrack
	mOutputTree->Branch("end_pz",end_pz ,"end_pz[nTrack]/D");//nTrack
	mOutputTree->Branch("end_vx",end_vx ,"end_vx[nTrack]/D");//nTrack
	mOutputTree->Branch("end_vy",end_vy ,"end_vy[nTrack]/D");//nTrack
	mOutputTree->Branch("end_vz",end_vz ,"end_vz[nTrack]/D");//nTrack
	*/
	return true;
}


Bool_t ConvertTrackData::SetBranchAddress(){
  std::cout << this->mInputTree->SetBranchAddress("GenParticle.",&(this->m_particleData)) << std::endl;
  return true;
}

Bool_t ConvertTrackData::Convert(){

  this->m_trackArr = (TClonesArray*)(this->m_particleData->briefTracks);
  this->nTrack     = this->m_trackArr->GetEntries();

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
    
    /*
    this->px[trackIndex]     = trackData->p.fX;
    this->py[trackIndex]     = trackData->p.fY;
    this->pz[trackIndex]     = trackData->p.fZ;
    this->vx[trackIndex]     = trackData->v.fX;
    this->vy[trackIndex]     = trackData->v.fY;
    this->vz[trackIndex]     = trackData->v.fZ;

    this->end_px[trackIndex] = trackData->end_p.fX;
    this->end_py[trackIndex] = trackData->end_p.fY;
    this->end_pz[trackIndex] = trackData->end_p.fZ;
    this->end_vx[trackIndex] = trackData->end_v.fX;
    this->end_vy[trackIndex] = trackData->end_v.fY;
    this->end_vz[trackIndex] = trackData->end_v.fZ;
    */
  }
  return true;
}
