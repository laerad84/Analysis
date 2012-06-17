/*
 * ConvertTrackData.h
 *
 *  Created on: Jun 18, 2011
 *      Author: jwlee
 */

#ifndef CONVERTTRACKDATA_H_
#define CONVERTTRACKDATA_H_

#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>

#include <GsimData/GsimGenParticleData.h>
#include <GsimData/GsimTrackData.h>

class ConvertTrackData
{
private:
	static const Int_t nMaxTrack =400;
	TTree* mInputTree;
	TTree* mOutputTree;
public:
	GsimTrackData*          m_trackData;
	GsimGenParticleData* 	m_particleData;
	TClonesArray*           m_trackArr;
	Int_t 		       	nTrack;
	UShort_t                track[nMaxTrack];//nTrack
	Short_t                 mother[nMaxTrack];//nTrack
	Int_t                   pid[nMaxTrack];//nTrack
	Float_t 		mass[nMaxTrack];//nTrack
	Float_t			ek[nMaxTrack];//nTrack
	Float_t			end_ek[nMaxTrack];//nTrack

	Double_t                p[nMaxTrack][3];//nTrack
	Double_t                v[nMaxTrack][3];//nTrack
	Double_t                end_p[nMaxTrack][3];//nTrack
	Double_t                end_v[nMaxTrack][3];//nTrack

	/*
	Double_t		px[nMaxTrack][3];//nTrack
	Double_t		py[nMaxTrack][3];//nTrack
	Double_t		pz[nMaxTrack][3];//nTrack
	Double_t 		vx[nMaxTrack][3];//nTrack
	Double_t 		vy[nMaxTrack][3];//nTrack
	Double_t 		vz[nMaxTrack][3];//nTrack

	Double_t		end_px[nMaxTrack][3];//nTrack
	Double_t		end_py[nMaxTrack][3];//nTrack
	Double_t		end_pz[nMaxTrack][3];//nTrack
	Double_t 		end_vx[nMaxTrack][3];//nTrack
	Double_t 		end_vy[nMaxTrack][3];//nTrack
	Double_t 		end_vz[nMaxTrack][3];//nTrack
	*/
	ConvertTrackData(TTree* ,TTree* );
	~ConvertTrackData();

	Bool_t SetBranchAddress();
	Bool_t Branch();
	Bool_t Convert();
};



#endif /* CONVERTTRACKDATA_H_ */
