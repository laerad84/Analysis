//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 29 17:29:06 2013 by ROOT version 5.30/04
// from TTree EventTree/
// found on file: run_conv_mode5.root
//////////////////////////////////////////////////////////

#ifndef Waveform_Reader_h
#define Waveform_Reader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class Waveform_Reader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           xaxis[48];
   Double_t        volt;
   Int_t           runID;
   Int_t           ientry;
   Double_t        data[48];
   Double_t        maxdata;
   Double_t        chi2;
   Double_t        pedestal;
   Double_t        param[5];
   Double_t        EventSum;

   // List of branches
   TBranch        *b_xaxis;   //!
   TBranch        *b_volt;   //!
   TBranch        *b_runid;   //!
   TBranch        *b_ientry;   //!
   TBranch        *b_data;   //!
   TBranch        *b_maxdata;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_pedestal;   //!
   TBranch        *b_param;   //!
   TBranch        *b_EventSum;   //!

   Waveform_Reader(TTree *tree=0);
   virtual ~Waveform_Reader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Waveform_Reader_cxx
Waveform_Reader::Waveform_Reader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

Waveform_Reader::~Waveform_Reader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Waveform_Reader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Waveform_Reader::LoadTree(Long64_t entry)
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

void Waveform_Reader::Init(TTree *tree)
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

   fChain->SetBranchAddress("xaxis", xaxis, &b_xaxis);
   fChain->SetBranchAddress("volt", &volt, &b_volt);
   fChain->SetBranchAddress("runID", &runID, &b_runid);
   fChain->SetBranchAddress("ientry", &ientry, &b_ientry);
   fChain->SetBranchAddress("data", data, &b_data);
   fChain->SetBranchAddress("maxdata", &maxdata, &b_maxdata);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("pedestal", &pedestal, &b_pedestal);
   fChain->SetBranchAddress("param", param, &b_param);
   fChain->SetBranchAddress("EventSum", &EventSum, &b_EventSum);
   Notify();
}

Bool_t Waveform_Reader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Waveform_Reader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Waveform_Reader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Waveform_Reader_cxx
