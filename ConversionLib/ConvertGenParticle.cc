#include "ConvertTrackData.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "GsimData/GsimGenParticleData.h"
#include "GsimData/GsimTrackData.h"
#include "TMath.h"

int main( int argc, char** argv){
  //argv[1] = InputFilename
  //argv[2] = outputFilename

  if( argc != 3 ){ return -1 ;}
  TFile* tf = new TFile(argv[1]);
  TTree* tr = (TTree*)tf->Get("eventTree00");
  GsimTrackData* trackData          = new GsimTrackData();
  GsimGenParticleData* particleData = new GsimGenParticleData();
  tr->SetBranchAddress("GenParticle.",&(particleData));
  UInt_t NCCnhit;
  //tr->SetBranchAddress("NCC.nHit",&NCCnhit);
  
  TFile* tfout = new TFile(argv[2],"recreate");
  TTree* trout = new TTree("HEPEvt","HEPEvt");
  int nMax = 300;
  int NHEP;
  int ISTHEP[nMax];
  int IDHEP[nMax];
  int JDAHEP1[nMax];
  int JDAHEP2[nMax];
  double PHEP1[nMax];
  double PHEP2[nMax];
  double PHEP3[nMax];
  double PHEP5[nMax];
  double VHEP1[nMax];
  double VHEP2[nMax];
  double VHEP3[nMax];
  double VHEP4[nMax];
  double PHEP0[nMax];
  double PHEP4[nMax];

  trout->Branch("NHEP",&NHEP,"NHEP/I");
  trout->Branch("IDHEP",IDHEP,"IDHEP[NHEP]/I");//NHEP
  trout->Branch("JDAHEP1",JDAHEP1,"JDAHEP1[NHEP]/I");//NHEP
  trout->Branch("JDAHEP2",JDAHEP2,"JDAHEP2[NHEP]/I");//NHEP
  trout->Branch("ISTHEP",ISTHEP,"ISTHEP[NHEP]/I");//NHEP
  trout->Branch("PHEP1",PHEP1,"PHEP1[NHEP]/D");//NHEP
  trout->Branch("PHEP2",PHEP2,"PHEP2[NHEP]/D");//NHEP
  trout->Branch("PHEP3",PHEP3,"PHEP3[NHEP]/D");//NHEP
  trout->Branch("PHEP5",PHEP5,"PHEP5[NHEP]/D");//NHEP
  trout->Branch("VHEP1",VHEP1,"VHEP1[NHEP]/D");//NHEP
  trout->Branch("VHEP2",VHEP2,"VHEP2[NHEP]/D");//NHEP
  trout->Branch("VHEP3",VHEP3,"VHEP3[NHEP]/D");//NHEP
  trout->Branch("VHEP4",VHEP4,"VHEP4[NHEP]/D");//NHEP
  trout->Branch("PHEP0",PHEP0,"PHEP0[NHEP]/D");//NHEP
  trout->Branch("PHEP4",PHEP4,"PHEP4[NHEP]/D");//NHEP

  TClonesArray* trackArr;
  for( int i = 0; i< tr->GetEntries(); i++){
    tr->GetEntry(i);    
    //std::cout<< NCCnhit << std::endl;
    //if(NCCnhit == 0){ continue;}
    for( int j = 0; j< 100; j++){//100 times duplicate
    NHEP=0;
    for( int it = 0; it < nMax; it++){
      ISTHEP[it] = 0;
      IDHEP[it] = 0;
      JDAHEP1[it] = 0;
      JDAHEP2[it] = 0;
      PHEP1[it] = 0;
      PHEP2[it] = 0;
      PHEP3[it] = 0;
      PHEP5[it] = 0;
      VHEP1[it] = 0;
      VHEP2[it] = 0;
      VHEP3[it] = 0;
      PHEP0[it] = 0;
      PHEP4[it] = 0;
    }

    trackArr = (TClonesArray*)particleData->briefTracks;
    int nTrack = trackArr->GetEntries();
    bool bEtaTrack = true;
    for( int itrack= 0; itrack < nTrack; itrack++){
      GsimTrackData* trackData = (GsimTrackData*)trackArr->UncheckedAt(itrack);
      /*
      if( trackData->mother == 1 ){
	if( trackData->pid != 221 ){
	  bEtaTrack = false;
	}
      }
      */
      //if(trackData->mother == 1 && trackData->pid == 221 ){
      if( trackData->mother != 1 ){ continue; }
      if( trackData->v[2] > 6080 ){ continue; }
      ISTHEP[NHEP] = 1;
	JDAHEP1[NHEP] = 0;
	JDAHEP2[NHEP] = 0;
	IDHEP[NHEP] = trackData->pid;
	VHEP1[NHEP] = trackData->v[0];
	VHEP2[NHEP] = trackData->v[1];
	VHEP3[NHEP] = trackData->v[2];
	VHEP4[NHEP] = 0;
	PHEP0[NHEP] = sqrt(pow(trackData->p[0]/1000,2)
			   +pow(trackData->p[1]/1000,2)
			   +pow(trackData->p[2]/1000,2));
	PHEP1[NHEP] = trackData->p[0]/1000;
	PHEP2[NHEP] = trackData->p[1]/1000;
	PHEP3[NHEP] = trackData->p[2]/1000;
	PHEP4[NHEP] = trackData->ek/1000;
	  /*
	  sqrt(pow(trackData->p[0]/1000,2)
			+pow(trackData->p[1]/1000,2)
			+pow(trackData->p[2]/1000,2)
			+pow(trackData->mass/1000,2))-trackData->mass/1000;
	  */
	PHEP5[NHEP] = trackData->mass/1000;
	NHEP++;
	//}
    }
    //if( !bEtaTrack ){ continue; }
    
    trout->Fill();

  }
  }

  trout->Write();
  tfout->Close();
  return 0;
}

