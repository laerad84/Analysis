//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 13 18:37:12 2013 by ROOT version 5.30/04
// from TTree eventTree00/eventTree00
// found on file: /Volume0/gamma/template_gamma_450MeV_10deg-1E5-0.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
   const Int_t kMaxCSI = 1;
   const Int_t kMaxCSI_hits = 4000;
   const Int_t kMaxCSI_digi = 1;
   const Int_t kMaxCSI_mtime = 1;
   const Int_t kMaxCSI_trig = 1;
   const Int_t kMaxEvent = 1;
   const Int_t kMaxGenParticle = 1;
   const Int_t kMaxGenParticle_briefTracks = 10;

class EventTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UShort_t        CSI_nHit;
   Float_t         CSI_totalEnergy;
   Int_t           CSI_hits_;
   Float_t         CSI_hits_time[kMaxCSI_hits];   //[CSI.hits_]
   Float_t         CSI_hits_edep[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_r_fX[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_r_fY[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_r_fZ[kMaxCSI_hits];   //[CSI.hits_]
   Float_t         CSI_hits_ek[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_p_fX[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_p_fY[kMaxCSI_hits];   //[CSI.hits_]
   Double_t        CSI_hits_p_fZ[kMaxCSI_hits];   //[CSI.hits_]


   // List of branches

   TBranch        *b_CSI_nHit;   //!
   TBranch        *b_CSI_totalEnergy;   //!
   TBranch        *b_CSI_hits_;   //!
   TBranch        *b_CSI_hits_time;   //!
   TBranch        *b_CSI_hits_edep;   //!
   TBranch        *b_CSI_hits_r_fX;   //!
   TBranch        *b_CSI_hits_r_fY;   //!
   TBranch        *b_CSI_hits_r_fZ;   //!
   TBranch        *b_CSI_hits_ek;   //!
   TBranch        *b_CSI_hits_p_fX;   //!
   TBranch        *b_CSI_hits_p_fY;   //!
   TBranch        *b_CSI_hits_p_fZ;   //!

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
  /*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Volume0/gamma/template_gamma_450MeV_10deg-1E5-0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Volume0/gamma/template_gamma_450MeV_10deg-1E5-0.root");
      }
      f->GetObject("eventTree00",tree);

   }
  */
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
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchAddress("CSI.nHit", &CSI_nHit, &b_CSI_nHit);
   fChain->SetBranchAddress("CSI.totalEnergy", &CSI_totalEnergy, &b_CSI_totalEnergy);
   fChain->SetBranchAddress("CSI.hits", &CSI_hits_, &b_CSI_hits_);
   fChain->SetBranchAddress("CSI.hits.time", CSI_hits_time, &b_CSI_hits_time);
   fChain->SetBranchAddress("CSI.hits.edep", CSI_hits_edep, &b_CSI_hits_edep);
   fChain->SetBranchAddress("CSI.hits.r.fX", CSI_hits_r_fX, &b_CSI_hits_r_fX);
   fChain->SetBranchAddress("CSI.hits.r.fY", CSI_hits_r_fY, &b_CSI_hits_r_fY);
   fChain->SetBranchAddress("CSI.hits.r.fZ", CSI_hits_r_fZ, &b_CSI_hits_r_fZ);
   fChain->SetBranchAddress("CSI.hits.ek", CSI_hits_ek, &b_CSI_hits_ek);
   fChain->SetBranchAddress("CSI.hits.p.fX", CSI_hits_p_fX, &b_CSI_hits_p_fX);
   fChain->SetBranchAddress("CSI.hits.p.fY", CSI_hits_p_fY, &b_CSI_hits_p_fY);
   fChain->SetBranchAddress("CSI.hits.p.fZ", CSI_hits_p_fZ, &b_CSI_hits_p_fZ);
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
