//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul  2 02:10:05 2012 by ROOT version 5.32/01
// from TTree CosmicOut/
// found on file: testTiming80k.root
//////////////////////////////////////////////////////////

#ifndef E14CosmicOut_h
#define E14CosmicOut_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class E14CosmicOut {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   static const int MaxChannel = 2716;
   // Declaration of leaf types
   Int_t           nDigi;
   Int_t           CsIID[2716];   //[nDigi]
   Int_t           CsIADC[2716];   //[nDigi]
   Double_t        CsIdepE[2716];   //[nDigi]
   Double_t        CsITiming[2716];   //[nDigi]
   Double_t        CsIHHTiming[2716];   //[nDigi]
   Double_t        CsIFitTiming[2716];   //[nDigi]
   Double_t        CsISplTiming[2716];   //[nDigi]
   Int_t           nHitUp;
   Int_t           nHitDn;
   Int_t           HitCoinUp;
   Int_t           HitCoinDn;
   Int_t           HitUp;
   Int_t           HitDn;
   Int_t           UpperID[1];   //[nHitUp]
   Int_t           DownID[1];   //[nHitDn]
   Int_t           Trigger;
   Double_t        CalFactor;
   Double_t        roh;
   Double_t        theta;
   Int_t           CosmicFit;
   Int_t           CosmicBoolUp;
   Int_t           CosmicBoolDn;

   // List of branches
   TBranch        *b_nDigi;   //!
   TBranch        *b_CsIID;   //!
   TBranch        *b_CsIADC;   //!
   TBranch        *b_CsIdepE;   //!
   TBranch        *b_CsITiming;   //!
   TBranch        *b_CsIHHTiming;   //!
   TBranch        *b_CsIFitTiming;   //!
   TBranch        *b_CsISplTiming;   //!
   TBranch        *b_nHitUp;   //!
   TBranch        *b_nHitDn;   //!
   TBranch        *b_HitCoinUp;   //!
   TBranch        *b_HitCoinDn;   //!
   TBranch        *b_HitUp;   //!
   TBranch        *b_HitDn;   //!
   TBranch        *b_UpperID;   //!
   TBranch        *b_DownID;   //!
   TBranch        *b_Trigger;   //!
   TBranch        *b_CalFactor;   //!
   TBranch        *b_roh;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_CosmicFit;   //!
   TBranch        *b_CosmicBoolUp;   //!
   TBranch        *b_CosmicBoolDn;   //!

   E14CosmicOut(TTree *tree=0);
   virtual ~E14CosmicOut();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef E14CosmicOut_cxx
E14CosmicOut::E14CosmicOut(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

E14CosmicOut::~E14CosmicOut()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t E14CosmicOut::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t E14CosmicOut::GetEntries(){
  if (!fChain) return 0;
  return fChain->GetEntries();
}
Long64_t E14CosmicOut::LoadTree(Long64_t entry)
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

void E14CosmicOut::Init(TTree *tree)
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

   fChain->SetBranchAddress("nDigi", &nDigi, &b_nDigi);
   fChain->SetBranchAddress("CsIID", CsIID, &b_CsIID);
   fChain->SetBranchAddress("CsIADC", CsIADC, &b_CsIADC);
   fChain->SetBranchAddress("CsIdepE", CsIdepE, &b_CsIdepE);
   fChain->SetBranchAddress("CsITiming", CsITiming, &b_CsITiming);
   fChain->SetBranchAddress("CsIHHTiming", CsIHHTiming, &b_CsIHHTiming);
   fChain->SetBranchAddress("CsIFitTiming", CsIFitTiming, &b_CsIFitTiming);
   fChain->SetBranchAddress("CsISplTiming", CsISplTiming, &b_CsISplTiming);
   fChain->SetBranchAddress("nHitUp", &nHitUp, &b_nHitUp);
   fChain->SetBranchAddress("nHitDn", &nHitDn, &b_nHitDn);
   fChain->SetBranchAddress("HitCoinUp", &HitCoinUp, &b_HitCoinUp);
   fChain->SetBranchAddress("HitCoinDn", &HitCoinDn, &b_HitCoinDn);
   fChain->SetBranchAddress("HitUp", &HitUp, &b_HitUp);
   fChain->SetBranchAddress("HitDn", &HitDn, &b_HitDn);
   fChain->SetBranchAddress("UpperID", &UpperID, &b_UpperID);
   fChain->SetBranchAddress("DownID", &DownID, &b_DownID);
   fChain->SetBranchAddress("Trigger", &Trigger, &b_Trigger);
   fChain->SetBranchAddress("CalFactor", &CalFactor, &b_CalFactor);
   fChain->SetBranchAddress("roh", &roh, &b_roh);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("CosmicFit", &CosmicFit, &b_CosmicFit);
   fChain->SetBranchAddress("CosmicBoolUp", &CosmicBoolUp, &b_CosmicBoolUp);
   fChain->SetBranchAddress("CosmicBoolDn", &CosmicBoolDn, &b_CosmicBoolDn);
   Notify();
}

Bool_t E14CosmicOut::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void E14CosmicOut::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t E14CosmicOut::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef E14CosmicOut_cxx
