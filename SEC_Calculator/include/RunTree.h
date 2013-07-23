//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 20 05:40:05 2013 by ROOT version 5.30/04
// from TTree runTree/
// found on file: run4200.root
//////////////////////////////////////////////////////////

#ifndef RunTree_h
#define RunTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class RunTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runID;
   UInt_t          StartDate;
   UInt_t          EndDate;
   Int_t           bufSize;
   Int_t           TrigType;
   Int_t           nTrigMAX;
   Int_t           trigdelay;
   Int_t           Threshold;
   Int_t           Threshold_offset;
   Int_t           SpillEtRecordCycle;
   Int_t           SECTMonRecordCycle;
   UInt_t          TotalnSpill;
   UInt_t          TotalnSEC;
   UInt_t          TotalnTMon;
   UInt_t          TotalnTrigRequested;
   UInt_t          TotalnTrigAccepted;
   UInt_t          TotalnEt;
   UInt_t          TotalnEt_withVeto;
   UInt_t          TotalnTrigOR;
   UInt_t          TotalnCalibTrig;
   UInt_t          TotalnChamberTrig;
   UInt_t          TotalnVetoCount;

   // List of branches
   TBranch        *b_runID;   //!
   TBranch        *b_StartDate;   //!
   TBranch        *b_EndDate;   //!
   TBranch        *b_bufSize;   //!
   TBranch        *b_TrigType;   //!
   TBranch        *b_nTrigMAX;   //!
   TBranch        *b_trigdelay;   //!
   TBranch        *b_Threshold;   //!
   TBranch        *b_Threshold_offset;   //!
   TBranch        *b_SpillEtRecordCycle;   //!
   TBranch        *b_SECTMonRecordCycle;   //!
   TBranch        *b_TotalnSpill;   //!
   TBranch        *b_TotalnSEC;   //!
   TBranch        *b_TotalnTMon;   //!
   TBranch        *b_TotalnTrigRequested;   //!
   TBranch        *b_TotalnTrigAccepted;   //!
   TBranch        *b_TotalnEt;   //!
   TBranch        *b_TotalnEt_withVeto;   //!
   TBranch        *b_TotalnTrigOR;   //!
   TBranch        *b_TotalnCalibTrig;   //!
   TBranch        *b_TotalnChamberTrig;   //!
   TBranch        *b_TotalnVetoCount;   //!

   RunTree(TTree *tree=0);
   virtual ~RunTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RunTree_cxx
RunTree::RunTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

RunTree::~RunTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RunTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RunTree::LoadTree(Long64_t entry)
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

void RunTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("StartDate", &StartDate, &b_StartDate);
   fChain->SetBranchAddress("EndDate", &EndDate, &b_EndDate);
   fChain->SetBranchAddress("bufSize", &bufSize, &b_bufSize);
   fChain->SetBranchAddress("TrigType", &TrigType, &b_TrigType);
   fChain->SetBranchAddress("nTrigMAX", &nTrigMAX, &b_nTrigMAX);
   fChain->SetBranchAddress("trigdelay", &trigdelay, &b_trigdelay);
   fChain->SetBranchAddress("Threshold", &Threshold, &b_Threshold);
   fChain->SetBranchAddress("Threshold_offset", &Threshold_offset, &b_Threshold_offset);
   fChain->SetBranchAddress("SpillEtRecordCycle", &SpillEtRecordCycle, &b_SpillEtRecordCycle);
   fChain->SetBranchAddress("SECTMonRecordCycle", &SECTMonRecordCycle, &b_SECTMonRecordCycle);
   fChain->SetBranchAddress("TotalnSpill", &TotalnSpill, &b_TotalnSpill);
   fChain->SetBranchAddress("TotalnSEC", &TotalnSEC, &b_TotalnSEC);
   fChain->SetBranchAddress("TotalnTMon", &TotalnTMon, &b_TotalnTMon);
   fChain->SetBranchAddress("TotalnTrigRequested", &TotalnTrigRequested, &b_TotalnTrigRequested);
   fChain->SetBranchAddress("TotalnTrigAccepted", &TotalnTrigAccepted, &b_TotalnTrigAccepted);
   fChain->SetBranchAddress("TotalnEt", &TotalnEt, &b_TotalnEt);
   fChain->SetBranchAddress("TotalnEt_withVeto", &TotalnEt_withVeto, &b_TotalnEt_withVeto);
   fChain->SetBranchAddress("TotalnTrigOR", &TotalnTrigOR, &b_TotalnTrigOR);
   fChain->SetBranchAddress("TotalnCalibTrig", &TotalnCalibTrig, &b_TotalnCalibTrig);
   fChain->SetBranchAddress("TotalnChamberTrig", &TotalnChamberTrig, &b_TotalnChamberTrig);
   fChain->SetBranchAddress("TotalnVetoCount", &TotalnVetoCount, &b_TotalnVetoCount);
   Notify();
}

Bool_t RunTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RunTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RunTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RunTree_cxx
