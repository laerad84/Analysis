//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 23 14:22:26 2012 by ROOT version 5.26/00
// from TTree EventTree/Conv data for 2012 run
// found on file: /disk/compute-1-11/partition2/conv_data/crate1/run1331_conv.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class EventTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Short_t         nFADC;
   Short_t         nSamples;
   Int_t           EventNo;
   Short_t         Data[20][16][48];
   Float_t         Pedestal[20][16];
   Short_t         PeakHeight[20][16];
   Short_t         PeakTime[20][16];
   Float_t         IntegratedADC[20][16];
   Int_t           EtSum_FADC[20][48];
   Long64_t        BufferLoopNo;
   Short_t         Error[20];
   Int_t           TimeStamp[20];
   Short_t         TrigNo[20];
   Short_t         SpillNo[20];
   Short_t         SlotNo[20];
   Short_t         Compression_flag[20][16];

   // List of branches
   TBranch        *b_nFADC;   //!
   TBranch        *b_nSamples;   //!
   TBranch        *b_EventNo;   //!
   TBranch        *b_Data;   //!
   TBranch        *b_Pedestal;   //!
   TBranch        *b_PeakHeight;   //!
   TBranch        *b_PeakTime;   //!
   TBranch        *b_IntegratedADC;   //!
   TBranch        *b_EtSum_FADC;   //!
   TBranch        *b_BufferLoopNo;   //!
   TBranch        *b_Error;   //!
   TBranch        *b_TimeStamp;   //!
   TBranch        *b_TrigNo;   //!
   TBranch        *b_SpillNo;   //!
   TBranch        *b_SlotNo;   //!
   TBranch        *b_Compression_flag;   //!

   EventTree(TTree *tree=0);
   virtual ~EventTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual Int_t    FileOpen( char *InFileName );

   TTree *tree;
};

#endif

#ifdef EventTree_cxx
EventTree::EventTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
/*
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/disk/compute-1-11/partition2/conv_data/crate1/run1331_conv.root");
      if (!f) {
         f = new TFile("/disk/compute-1-11/partition2/conv_data/crate1/run1331_conv.root");
      }
      tree = (TTree*)gDirectory->Get("EventTree");

   }
   Init(tree);
*/
}

EventTree::~EventTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventTree::FileOpen( char *InFileName )
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject( InFileName );
      if (!f) {
         f = new TFile( InFileName );
      }
      tree = (TTree*)gDirectory->Get("EventTree");

   }
   Init(tree);
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
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
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

   fChain->SetBranchAddress("nFADC", &nFADC, &b_nFADC);
   fChain->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
   fChain->SetBranchAddress("EventNo", &EventNo, &b_EventNo);
   fChain->SetBranchAddress("Data", Data, &b_Data);
   fChain->SetBranchAddress("Pedestal", Pedestal, &b_Pedestal);
   fChain->SetBranchAddress("PeakHeight", PeakHeight, &b_PeakHeight);
   fChain->SetBranchAddress("PeakTime", PeakTime, &b_PeakTime);
   fChain->SetBranchAddress("IntegratedADC", IntegratedADC, &b_IntegratedADC);
   fChain->SetBranchAddress("EtSum_FADC", EtSum_FADC, &b_EtSum_FADC);
   fChain->SetBranchAddress("BufferLoopNo", &BufferLoopNo, &b_BufferLoopNo);
   fChain->SetBranchAddress("Error", Error, &b_Error);
   fChain->SetBranchAddress("TimeStamp", TimeStamp, &b_TimeStamp);
   fChain->SetBranchAddress("TrigNo", TrigNo, &b_TrigNo);
   fChain->SetBranchAddress("SpillNo", SpillNo, &b_SpillNo);
   fChain->SetBranchAddress("SlotNo", SlotNo, &b_SlotNo);
   fChain->SetBranchAddress("Compression_flag", Compression_flag, &b_Compression_flag);
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
