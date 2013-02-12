//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov  5 15:48:44 2012 by ROOT version 5.30/04
// from TTree eventTree00/eventTree00
// found on file: /Volume0/gamma/template_gamma_210MeV_10deg-1E5-0.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

const Int_t kMaxCSI = 1;
const Int_t kMaxCSI_hits = 6000;
const Int_t kMaxCSI_digi = 1;
const Int_t kMaxCSI_mtime = 1;
const Int_t kMaxCSI_trig = 1;
const Int_t kMaxEvent = 1;
const Int_t kMaxGenParticle = 1;
const Int_t kMaxGenParticle_briefTracks = 1000;

class EventTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
 //GsimDetectorEventData *CSI_;
   UInt_t          CSI_TObject_fUniqueID;
   UInt_t          CSI_TObject_fBits;
   UShort_t        CSI_status;
   UShort_t        CSI_nHit;
   UShort_t        CSI_nDigi;
   UShort_t        CSI_nTrig;
   Float_t         CSI_totalEnergy;
   Float_t         CSI_firstHitTime;
   Float_t         CSI_lastHitTime;
   Int_t           CSI_hits_;
   UInt_t          CSI_hits_fUniqueID[kMaxCSI_hits];   //[CSI.hits_]
   UInt_t          CSI_hits_fBits[kMaxCSI_hits];   //[CSI.hits_]
   UShort_t        CSI_hits_thisID[kMaxCSI_hits];   //[CSI.hits_]
   UShort_t        CSI_hits_track[kMaxCSI_hits];   //[CSI.hits_]
   UShort_t        CSI_hits_stop[kMaxCSI_hits];   //[CSI.hits_]
   UShort_t        CSI_hits_hitChannel[kMaxCSI_hits];   //[CSI.hits_]
   Float_t         CSI_hits_time[kMaxCSI_hits];   //[CSI.hits_]
   Float_t         CSI_hits_edep[kMaxCSI_hits];   //[CSI.hits_]
   Int_t           CSI_hits_pid[kMaxCSI_hits];   //[CSI.hits_]
   UInt_t          CSI_hits_r_fUniqueID[kMaxCSI_hits];   //[CSI.hits_]
   UInt_t          CSI_hits_r_fBits[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_r_fX[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_r_fY[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_r_fZ[kMaxCSI_hits];   //[CSI.hits_]
   Float_t         CSI_hits_ek[kMaxCSI_hits];   //[CSI.hits_]
   UInt_t          CSI_hits_p_fUniqueID[kMaxCSI_hits];   //[CSI.hits_]
   UInt_t          CSI_hits_p_fBits[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_p_fX[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_p_fY[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_p_fZ[kMaxCSI_hits];   //[CSI.hits_]
   Int_t           CSI_digi_;
   UInt_t          CSI_digi_fUniqueID[kMaxCSI_digi];   //[CSI.digi_]
   UInt_t          CSI_digi_fBits[kMaxCSI_digi];   //[CSI.digi_]
   UShort_t        CSI_digi_detID[kMaxCSI_digi];   //[CSI.digi_]
   UShort_t        CSI_digi_modID[kMaxCSI_digi];   //[CSI.digi_]
   Float_t         CSI_digi_energy[kMaxCSI_digi];   //[CSI.digi_]
   Float_t         CSI_digi_time[kMaxCSI_digi];   //[CSI.digi_]
   UShort_t        CSI_digi_thisID[kMaxCSI_digi];   //[CSI.digi_]
   UShort_t        CSI_digi_status[kMaxCSI_digi];   //[CSI.digi_]
   UShort_t        CSI_digi_track[kMaxCSI_digi];   //[CSI.digi_]
   UShort_t        CSI_digi_mtimeEntry[kMaxCSI_digi];   //[CSI.digi_]
   UShort_t        CSI_digi_mtimeSize[kMaxCSI_digi];   //[CSI.digi_]
   Int_t           CSI_mtime_;
   UInt_t          CSI_mtime_fUniqueID[kMaxCSI_mtime];   //[CSI.mtime_]
   UInt_t          CSI_mtime_fBits[kMaxCSI_mtime];   //[CSI.mtime_]
   UShort_t        CSI_mtime_modID[kMaxCSI_mtime];   //[CSI.mtime_]
   Float_t         CSI_mtime_energy[kMaxCSI_mtime];   //[CSI.mtime_]
   Float_t         CSI_mtime_time[kMaxCSI_mtime];   //[CSI.mtime_]
   Int_t           CSI_trig_;
   UInt_t          CSI_trig_fUniqueID[kMaxCSI_trig];   //[CSI.trig_]
   UInt_t          CSI_trig_fBits[kMaxCSI_trig];   //[CSI.trig_]
   UShort_t        CSI_trig_modID[kMaxCSI_trig];   //[CSI.trig_]
   Float_t         CSI_trig_energy[kMaxCSI_trig];   //[CSI.trig_]
   Float_t         CSI_trig_time[kMaxCSI_trig];   //[CSI.trig_]
 //GsimEventData   *Event_;
   UInt_t          Event_TObject_fUniqueID;
   UInt_t          Event_TObject_fBits;
   UShort_t        Event_expMC;
   UShort_t        Event_run_number;
   UShort_t        Event_spill_number;
   UInt_t          Event_event_number;
   UInt_t          Event_trigger;
   UInt_t          Event_time_stamp;
   UInt_t          Event_status;
   UShort_t        Event_version;
 //GsimGenParticleData *GenParticle_;
   UInt_t          GenParticle_TObject_fUniqueID;
   UInt_t          GenParticle_TObject_fBits;
   Int_t           GenParticle_briefTracks_;
   UInt_t          GenParticle_briefTracks_fUniqueID[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_fBits[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UShort_t        GenParticle_briefTracks_track[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Short_t         GenParticle_briefTracks_mother[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Int_t           GenParticle_briefTracks_pid[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_p_fUniqueID[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_p_fBits[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_p_fX[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_p_fY[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_p_fZ[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Float_t         GenParticle_briefTracks_ek[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Float_t         GenParticle_briefTracks_mass[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Float_t         GenParticle_briefTracks_time[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_v_fUniqueID[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_v_fBits[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_v_fX[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_v_fY[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_v_fZ[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_end_p_fUniqueID[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_end_p_fBits[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_end_p_fX[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_end_p_fY[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_end_p_fZ[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Float_t         GenParticle_briefTracks_end_ek[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Float_t         GenParticle_briefTracks_end_time[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_end_v_fUniqueID[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UInt_t          GenParticle_briefTracks_end_v_fBits[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_end_v_fX[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_end_v_fY[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Double_t        GenParticle_briefTracks_end_v_fZ[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   Short_t         GenParticle_briefTracks_mech[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UShort_t        GenParticle_briefTracks_status[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UShort_t        GenParticle_briefTracks_thisID[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   ULong64_t       GenParticle_briefTracks_history[kMaxGenParticle_briefTracks];   //[GenParticle.briefTracks_]
   UShort_t        GenParticle_version;

   // List of branches
   TBranch        *b_CSI_TObject_fUniqueID;   //!
   TBranch        *b_CSI_TObject_fBits;   //!
   TBranch        *b_CSI_status;   //!
   TBranch        *b_CSI_nHit;   //!
   TBranch        *b_CSI_nDigi;   //!
   TBranch        *b_CSI_nTrig;   //!
   TBranch        *b_CSI_totalEnergy;   //!
   TBranch        *b_CSI_firstHitTime;   //!
   TBranch        *b_CSI_lastHitTime;   //!
   TBranch        *b_CSI_hits_;   //!
   TBranch        *b_CSI_hits_fUniqueID;   //!
   TBranch        *b_CSI_hits_fBits;   //!
   TBranch        *b_CSI_hits_thisID;   //!
   TBranch        *b_CSI_hits_track;   //!
   TBranch        *b_CSI_hits_stop;   //!
   TBranch        *b_CSI_hits_hitChannel;   //!
   TBranch        *b_CSI_hits_time;   //!
   TBranch        *b_CSI_hits_edep;   //!
   TBranch        *b_CSI_hits_pid;   //!
   TBranch        *b_CSI_hits_r_fUniqueID;   //!
   TBranch        *b_CSI_hits_r_fBits;   //!
   TBranch        *b_CSI_hits_r_fX;   //!
   TBranch        *b_CSI_hits_r_fY;   //!
   TBranch        *b_CSI_hits_r_fZ;   //!
   TBranch        *b_CSI_hits_ek;   //!
   TBranch        *b_CSI_hits_p_fUniqueID;   //!
   TBranch        *b_CSI_hits_p_fBits;   //!
   TBranch        *b_CSI_hits_p_fX;   //!
   TBranch        *b_CSI_hits_p_fY;   //!
   TBranch        *b_CSI_hits_p_fZ;   //!
   TBranch        *b_CSI_digi_;   //!
   TBranch        *b_CSI_digi_fUniqueID;   //!
   TBranch        *b_CSI_digi_fBits;   //!
   TBranch        *b_CSI_digi_detID;   //!
   TBranch        *b_CSI_digi_modID;   //!
   TBranch        *b_CSI_digi_energy;   //!
   TBranch        *b_CSI_digi_time;   //!
   TBranch        *b_CSI_digi_thisID;   //!
   TBranch        *b_CSI_digi_status;   //!
   TBranch        *b_CSI_digi_track;   //!
   TBranch        *b_CSI_digi_mtimeEntry;   //!
   TBranch        *b_CSI_digi_mtimeSize;   //!
   TBranch        *b_CSI_mtime_;   //!
   TBranch        *b_CSI_mtime_fUniqueID;   //!
   TBranch        *b_CSI_mtime_fBits;   //!
   TBranch        *b_CSI_mtime_modID;   //!
   TBranch        *b_CSI_mtime_energy;   //!
   TBranch        *b_CSI_mtime_time;   //!
   TBranch        *b_CSI_trig_;   //!
   TBranch        *b_CSI_trig_fUniqueID;   //!
   TBranch        *b_CSI_trig_fBits;   //!
   TBranch        *b_CSI_trig_modID;   //!
   TBranch        *b_CSI_trig_energy;   //!
   TBranch        *b_CSI_trig_time;   //!
   TBranch        *b_Event_TObject_fUniqueID;   //!
   TBranch        *b_Event_TObject_fBits;   //!
   TBranch        *b_Event_expMC;   //!
   TBranch        *b_Event_run_number;   //!
   TBranch        *b_Event_spill_number;   //!
   TBranch        *b_Event_event_number;   //!
   TBranch        *b_Event_trigger;   //!
   TBranch        *b_Event_time_stamp;   //!
   TBranch        *b_Event_status;   //!
   TBranch        *b_Event_version;   //!
   TBranch        *b_GenParticle_TObject_fUniqueID;   //!
   TBranch        *b_GenParticle_TObject_fBits;   //!
   TBranch        *b_GenParticle_briefTracks_;   //!
   TBranch        *b_GenParticle_briefTracks_fUniqueID;   //!
   TBranch        *b_GenParticle_briefTracks_fBits;   //!
   TBranch        *b_GenParticle_briefTracks_track;   //!
   TBranch        *b_GenParticle_briefTracks_mother;   //!
   TBranch        *b_GenParticle_briefTracks_pid;   //!
   TBranch        *b_GenParticle_briefTracks_p_fUniqueID;   //!
   TBranch        *b_GenParticle_briefTracks_p_fBits;   //!
   TBranch        *b_GenParticle_briefTracks_p_fX;   //!
   TBranch        *b_GenParticle_briefTracks_p_fY;   //!
   TBranch        *b_GenParticle_briefTracks_p_fZ;   //!
   TBranch        *b_GenParticle_briefTracks_ek;   //!
   TBranch        *b_GenParticle_briefTracks_mass;   //!
   TBranch        *b_GenParticle_briefTracks_time;   //!
   TBranch        *b_GenParticle_briefTracks_v_fUniqueID;   //!
   TBranch        *b_GenParticle_briefTracks_v_fBits;   //!
   TBranch        *b_GenParticle_briefTracks_v_fX;   //!
   TBranch        *b_GenParticle_briefTracks_v_fY;   //!
   TBranch        *b_GenParticle_briefTracks_v_fZ;   //!
   TBranch        *b_GenParticle_briefTracks_end_p_fUniqueID;   //!
   TBranch        *b_GenParticle_briefTracks_end_p_fBits;   //!
   TBranch        *b_GenParticle_briefTracks_end_p_fX;   //!
   TBranch        *b_GenParticle_briefTracks_end_p_fY;   //!
   TBranch        *b_GenParticle_briefTracks_end_p_fZ;   //!
   TBranch        *b_GenParticle_briefTracks_end_ek;   //!
   TBranch        *b_GenParticle_briefTracks_end_time;   //!
   TBranch        *b_GenParticle_briefTracks_end_v_fUniqueID;   //!
   TBranch        *b_GenParticle_briefTracks_end_v_fBits;   //!
   TBranch        *b_GenParticle_briefTracks_end_v_fX;   //!
   TBranch        *b_GenParticle_briefTracks_end_v_fY;   //!
   TBranch        *b_GenParticle_briefTracks_end_v_fZ;   //!
   TBranch        *b_GenParticle_briefTracks_mech;   //!
   TBranch        *b_GenParticle_briefTracks_status;   //!
   TBranch        *b_GenParticle_briefTracks_thisID;   //!
   TBranch        *b_GenParticle_briefTracks_history;   //!
   TBranch        *b_GenParticle_version;   //!

   EventTree(TTree *tree=0);
   virtual ~EventTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EventTree_cxx
EventTree::EventTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

EventTree::~EventTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("CSI.TObject.fUniqueID", &CSI_TObject_fUniqueID, &b_CSI_TObject_fUniqueID);
   fChain->SetBranchAddress("CSI.TObject.fBits", &CSI_TObject_fBits, &b_CSI_TObject_fBits);
   fChain->SetBranchAddress("CSI.status", &CSI_status, &b_CSI_status);
   fChain->SetBranchAddress("CSI.nHit", &CSI_nHit, &b_CSI_nHit);
   fChain->SetBranchAddress("CSI.nDigi", &CSI_nDigi, &b_CSI_nDigi);
   fChain->SetBranchAddress("CSI.nTrig", &CSI_nTrig, &b_CSI_nTrig);
   fChain->SetBranchAddress("CSI.totalEnergy", &CSI_totalEnergy, &b_CSI_totalEnergy);
   fChain->SetBranchAddress("CSI.firstHitTime", &CSI_firstHitTime, &b_CSI_firstHitTime);
   fChain->SetBranchAddress("CSI.lastHitTime", &CSI_lastHitTime, &b_CSI_lastHitTime);
   fChain->SetBranchAddress("CSI.hits", &CSI_hits_, &b_CSI_hits_);
   fChain->SetBranchAddress("CSI.hits.fUniqueID", CSI_hits_fUniqueID, &b_CSI_hits_fUniqueID);
   fChain->SetBranchAddress("CSI.hits.fBits", CSI_hits_fBits, &b_CSI_hits_fBits);
   fChain->SetBranchAddress("CSI.hits.thisID", CSI_hits_thisID, &b_CSI_hits_thisID);
   fChain->SetBranchAddress("CSI.hits.track", CSI_hits_track, &b_CSI_hits_track);
   fChain->SetBranchAddress("CSI.hits.stop", CSI_hits_stop, &b_CSI_hits_stop);
   fChain->SetBranchAddress("CSI.hits.hitChannel", CSI_hits_hitChannel, &b_CSI_hits_hitChannel);
   fChain->SetBranchAddress("CSI.hits.time", CSI_hits_time, &b_CSI_hits_time);
   fChain->SetBranchAddress("CSI.hits.edep", CSI_hits_edep, &b_CSI_hits_edep);
   fChain->SetBranchAddress("CSI.hits.pid", CSI_hits_pid, &b_CSI_hits_pid);
   fChain->SetBranchAddress("CSI.hits.r.fUniqueID", CSI_hits_r_fUniqueID, &b_CSI_hits_r_fUniqueID);
   fChain->SetBranchAddress("CSI.hits.r.fBits", CSI_hits_r_fBits, &b_CSI_hits_r_fBits);
   fChain->SetBranchAddress("CSI.hits.r.fX", CSI_hits_r_fX, &b_CSI_hits_r_fX);
   fChain->SetBranchAddress("CSI.hits.r.fY", CSI_hits_r_fY, &b_CSI_hits_r_fY);
   fChain->SetBranchAddress("CSI.hits.r.fZ", CSI_hits_r_fZ, &b_CSI_hits_r_fZ);
   fChain->SetBranchAddress("CSI.hits.ek", CSI_hits_ek, &b_CSI_hits_ek);
   fChain->SetBranchAddress("CSI.hits.p.fUniqueID", CSI_hits_p_fUniqueID, &b_CSI_hits_p_fUniqueID);
   fChain->SetBranchAddress("CSI.hits.p.fBits", CSI_hits_p_fBits, &b_CSI_hits_p_fBits);
   fChain->SetBranchAddress("CSI.hits.p.fX", CSI_hits_p_fX, &b_CSI_hits_p_fX);
   fChain->SetBranchAddress("CSI.hits.p.fY", CSI_hits_p_fY, &b_CSI_hits_p_fY);
   fChain->SetBranchAddress("CSI.hits.p.fZ", CSI_hits_p_fZ, &b_CSI_hits_p_fZ);
   fChain->SetBranchAddress("CSI.digi", &CSI_digi_, &b_CSI_digi_);
   fChain->SetBranchAddress("CSI.digi.fUniqueID", &CSI_digi_fUniqueID, &b_CSI_digi_fUniqueID);
   fChain->SetBranchAddress("CSI.digi.fBits", &CSI_digi_fBits, &b_CSI_digi_fBits);
   fChain->SetBranchAddress("CSI.digi.detID", &CSI_digi_detID, &b_CSI_digi_detID);
   fChain->SetBranchAddress("CSI.digi.modID", &CSI_digi_modID, &b_CSI_digi_modID);
   fChain->SetBranchAddress("CSI.digi.energy", &CSI_digi_energy, &b_CSI_digi_energy);
   fChain->SetBranchAddress("CSI.digi.time", &CSI_digi_time, &b_CSI_digi_time);
   fChain->SetBranchAddress("CSI.digi.thisID", &CSI_digi_thisID, &b_CSI_digi_thisID);
   fChain->SetBranchAddress("CSI.digi.status", &CSI_digi_status, &b_CSI_digi_status);
   fChain->SetBranchAddress("CSI.digi.track", &CSI_digi_track, &b_CSI_digi_track);
   fChain->SetBranchAddress("CSI.digi.mtimeEntry", &CSI_digi_mtimeEntry, &b_CSI_digi_mtimeEntry);
   fChain->SetBranchAddress("CSI.digi.mtimeSize", &CSI_digi_mtimeSize, &b_CSI_digi_mtimeSize);
   fChain->SetBranchAddress("CSI.mtime", &CSI_mtime_, &b_CSI_mtime_);
   fChain->SetBranchAddress("CSI.mtime.fUniqueID", &CSI_mtime_fUniqueID, &b_CSI_mtime_fUniqueID);
   fChain->SetBranchAddress("CSI.mtime.fBits", &CSI_mtime_fBits, &b_CSI_mtime_fBits);
   fChain->SetBranchAddress("CSI.mtime.modID", &CSI_mtime_modID, &b_CSI_mtime_modID);
   fChain->SetBranchAddress("CSI.mtime.energy", &CSI_mtime_energy, &b_CSI_mtime_energy);
   fChain->SetBranchAddress("CSI.mtime.time", &CSI_mtime_time, &b_CSI_mtime_time);
   fChain->SetBranchAddress("CSI.trig", &CSI_trig_, &b_CSI_trig_);
   fChain->SetBranchAddress("CSI.trig.fUniqueID", &CSI_trig_fUniqueID, &b_CSI_trig_fUniqueID);
   fChain->SetBranchAddress("CSI.trig.fBits", &CSI_trig_fBits, &b_CSI_trig_fBits);
   fChain->SetBranchAddress("CSI.trig.modID", &CSI_trig_modID, &b_CSI_trig_modID);
   fChain->SetBranchAddress("CSI.trig.energy", &CSI_trig_energy, &b_CSI_trig_energy);
   fChain->SetBranchAddress("CSI.trig.time", &CSI_trig_time, &b_CSI_trig_time);
   fChain->SetBranchAddress("Event.TObject.fUniqueID", &Event_TObject_fUniqueID, &b_Event_TObject_fUniqueID);
   fChain->SetBranchAddress("Event.TObject.fBits", &Event_TObject_fBits, &b_Event_TObject_fBits);
   fChain->SetBranchAddress("Event.expMC", &Event_expMC, &b_Event_expMC);
   fChain->SetBranchAddress("Event.run_number", &Event_run_number, &b_Event_run_number);
   fChain->SetBranchAddress("Event.spill_number", &Event_spill_number, &b_Event_spill_number);
   fChain->SetBranchAddress("Event.event_number", &Event_event_number, &b_Event_event_number);
   fChain->SetBranchAddress("Event.trigger", &Event_trigger, &b_Event_trigger);
   fChain->SetBranchAddress("Event.time_stamp", &Event_time_stamp, &b_Event_time_stamp);
   fChain->SetBranchAddress("Event.status", &Event_status, &b_Event_status);
   fChain->SetBranchAddress("Event.version", &Event_version, &b_Event_version);
   fChain->SetBranchAddress("GenParticle.TObject.fUniqueID", &GenParticle_TObject_fUniqueID, &b_GenParticle_TObject_fUniqueID);
   fChain->SetBranchAddress("GenParticle.TObject.fBits", &GenParticle_TObject_fBits, &b_GenParticle_TObject_fBits);
   fChain->SetBranchAddress("GenParticle.briefTracks", &GenParticle_briefTracks_, &b_GenParticle_briefTracks_);
   fChain->SetBranchAddress("GenParticle.briefTracks.fUniqueID", GenParticle_briefTracks_fUniqueID, &b_GenParticle_briefTracks_fUniqueID);
   fChain->SetBranchAddress("GenParticle.briefTracks.fBits", GenParticle_briefTracks_fBits, &b_GenParticle_briefTracks_fBits);
   fChain->SetBranchAddress("GenParticle.briefTracks.track", GenParticle_briefTracks_track, &b_GenParticle_briefTracks_track);
   fChain->SetBranchAddress("GenParticle.briefTracks.mother", GenParticle_briefTracks_mother, &b_GenParticle_briefTracks_mother);
   fChain->SetBranchAddress("GenParticle.briefTracks.pid", GenParticle_briefTracks_pid, &b_GenParticle_briefTracks_pid);
   fChain->SetBranchAddress("GenParticle.briefTracks.p.fUniqueID", GenParticle_briefTracks_p_fUniqueID, &b_GenParticle_briefTracks_p_fUniqueID);
   fChain->SetBranchAddress("GenParticle.briefTracks.p.fBits", GenParticle_briefTracks_p_fBits, &b_GenParticle_briefTracks_p_fBits);
   fChain->SetBranchAddress("GenParticle.briefTracks.p.fX", GenParticle_briefTracks_p_fX, &b_GenParticle_briefTracks_p_fX);
   fChain->SetBranchAddress("GenParticle.briefTracks.p.fY", GenParticle_briefTracks_p_fY, &b_GenParticle_briefTracks_p_fY);
   fChain->SetBranchAddress("GenParticle.briefTracks.p.fZ", GenParticle_briefTracks_p_fZ, &b_GenParticle_briefTracks_p_fZ);
   fChain->SetBranchAddress("GenParticle.briefTracks.ek", GenParticle_briefTracks_ek, &b_GenParticle_briefTracks_ek);
   fChain->SetBranchAddress("GenParticle.briefTracks.mass", GenParticle_briefTracks_mass, &b_GenParticle_briefTracks_mass);
   fChain->SetBranchAddress("GenParticle.briefTracks.time", GenParticle_briefTracks_time, &b_GenParticle_briefTracks_time);
   fChain->SetBranchAddress("GenParticle.briefTracks.v.fUniqueID", GenParticle_briefTracks_v_fUniqueID, &b_GenParticle_briefTracks_v_fUniqueID);
   fChain->SetBranchAddress("GenParticle.briefTracks.v.fBits", GenParticle_briefTracks_v_fBits, &b_GenParticle_briefTracks_v_fBits);
   fChain->SetBranchAddress("GenParticle.briefTracks.v.fX", GenParticle_briefTracks_v_fX, &b_GenParticle_briefTracks_v_fX);
   fChain->SetBranchAddress("GenParticle.briefTracks.v.fY", GenParticle_briefTracks_v_fY, &b_GenParticle_briefTracks_v_fY);
   fChain->SetBranchAddress("GenParticle.briefTracks.v.fZ", GenParticle_briefTracks_v_fZ, &b_GenParticle_briefTracks_v_fZ);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_p.fUniqueID", GenParticle_briefTracks_end_p_fUniqueID, &b_GenParticle_briefTracks_end_p_fUniqueID);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_p.fBits", GenParticle_briefTracks_end_p_fBits, &b_GenParticle_briefTracks_end_p_fBits);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_p.fX", GenParticle_briefTracks_end_p_fX, &b_GenParticle_briefTracks_end_p_fX);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_p.fY", GenParticle_briefTracks_end_p_fY, &b_GenParticle_briefTracks_end_p_fY);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_p.fZ", GenParticle_briefTracks_end_p_fZ, &b_GenParticle_briefTracks_end_p_fZ);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_ek", GenParticle_briefTracks_end_ek, &b_GenParticle_briefTracks_end_ek);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_time", GenParticle_briefTracks_end_time, &b_GenParticle_briefTracks_end_time);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_v.fUniqueID", GenParticle_briefTracks_end_v_fUniqueID, &b_GenParticle_briefTracks_end_v_fUniqueID);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_v.fBits", GenParticle_briefTracks_end_v_fBits, &b_GenParticle_briefTracks_end_v_fBits);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_v.fX", GenParticle_briefTracks_end_v_fX, &b_GenParticle_briefTracks_end_v_fX);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_v.fY", GenParticle_briefTracks_end_v_fY, &b_GenParticle_briefTracks_end_v_fY);
   fChain->SetBranchAddress("GenParticle.briefTracks.end_v.fZ", GenParticle_briefTracks_end_v_fZ, &b_GenParticle_briefTracks_end_v_fZ);
   fChain->SetBranchAddress("GenParticle.briefTracks.mech", GenParticle_briefTracks_mech, &b_GenParticle_briefTracks_mech);
   fChain->SetBranchAddress("GenParticle.briefTracks.status", GenParticle_briefTracks_status, &b_GenParticle_briefTracks_status);
   fChain->SetBranchAddress("GenParticle.briefTracks.thisID", GenParticle_briefTracks_thisID, &b_GenParticle_briefTracks_thisID);
   fChain->SetBranchAddress("GenParticle.briefTracks.history", GenParticle_briefTracks_history, &b_GenParticle_briefTracks_history);
   fChain->SetBranchAddress("GenParticle.version", &GenParticle_version, &b_GenParticle_version);
   Notify();
}

Bool_t EventTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EventTree_cxx
