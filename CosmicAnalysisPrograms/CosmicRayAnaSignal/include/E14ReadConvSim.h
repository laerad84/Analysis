//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  4 23:37:36 2012 by ROOT version 5.32/01
// from TTree T/
// found on file: CosmicSim/Conv_Cosmic/cosmic_conv_0.root
//////////////////////////////////////////////////////////

#ifndef E14ReadConvSim_h
#define E14ReadConvSim_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class E14ReadConvSim {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EventNum;
   Int_t           CsiNumber;
   Double_t        CsiTotalE;
   Int_t           CsiModID[2716];   //CsiNumber
   Double_t        CsiEne[2716];   //CsiNumber
   Double_t        CsiTime[2716];   //CsiNumber
   Int_t           CC03Number;
   Double_t        CC03TotalE;
   Int_t           CC03ModID[32];   //CC03Number
   Double_t        CC03Ene[32];   //CC03Number
   Double_t        CC03Time[32];   //CC03Number
   Int_t           CVNumber;
   Double_t        CVTotalE;
   Int_t           CVModID[44];   //CVNumber
   Double_t        CVEne[44];   //CVNumber
   Double_t        CVTime[44];   //CVNumber
   Int_t           OEVNumber;
   Double_t        OEVTotalE;
   Int_t           OEVModID[44];   //OEVNumber
   Double_t        OEVEne[44];   //OEVNumber
   Double_t        OEVTime[44];   //OEVNumber
   Int_t           LaserNumber;
   Double_t        LaserTotalE;
   Int_t           LaserModID[5];   //LaserNumber
   Double_t        LaserEne[5];   //LaserNumber
   Double_t        LaserTime[5];   //LaserNumber
   Int_t           CosmicNumber;
   Double_t        CosmicTotalE;
   Int_t           CosmicModID[20];   //CosmicNumber
   Double_t        CosmicEne[20];   //CosmicNumber
   Double_t        CosmicTime[20];   //CosmicNumber
   Int_t           sciNumber;
   Double_t        sciTotalE;
   Int_t           sciModID[10];   //sciNumber
   Double_t        sciEne[10];   //sciNumber
   Double_t        sciTime[10];   //sciNumber
   Int_t           nTrack;
   UShort_t        track[256];   //[nTrack]
   Short_t         mother[256];   //[nTrack]
   Int_t           pid[256];   //[nTrack]
   Float_t         mass[256];   //[nTrack]
   Float_t         ek[256];   //[nTrack]
   Float_t         end_ek[256];   //[nTrack]
   Double_t        p[256][3];   //[nTrack]
   Double_t        end_p[256][3];   //[nTrack]
   Double_t        v[256][3];   //[nTrack]
   Double_t        end_v[256][3];   //[nTrack]

   // List of branches
   TBranch        *b_EventNum;   //!
   TBranch        *b_CsiNumber;   //!
   TBranch        *b_CsiTotalE;   //!
   TBranch        *b_CsiModID;   //!
   TBranch        *b_CsiEne;   //!
   TBranch        *b_CsiTime;   //!
   TBranch        *b_CC03Number;   //!
   TBranch        *b_CC03TotalE;   //!
   TBranch        *b_CC03ModID;   //!
   TBranch        *b_CC03Ene;   //!
   TBranch        *b_CC03Time;   //!
   TBranch        *b_CVNumber;   //!
   TBranch        *b_CVTotalE;   //!
   TBranch        *b_CVModID;   //!
   TBranch        *b_CVEne;   //!
   TBranch        *b_CVTime;   //!
   TBranch        *b_OEVNumber;   //!
   TBranch        *b_OEVTotalE;   //!
   TBranch        *b_OEVModID;   //!
   TBranch        *b_OEVEne;   //!
   TBranch        *b_OEVTime;   //!
   TBranch        *b_LaserNumber;   //!
   TBranch        *b_LaserTotalE;   //!
   TBranch        *b_LaserModID;   //!
   TBranch        *b_LaserEne;   //!
   TBranch        *b_LaserTime;   //!
   TBranch        *b_CosmicNumber;   //!
   TBranch        *b_CosmicTotalE;   //!
   TBranch        *b_CosmicModID;   //!
   TBranch        *b_CosmicEne;   //!
   TBranch        *b_CosmicTime;   //!
   TBranch        *b_sciNumber;   //!
   TBranch        *b_sciTotalE;   //!
   TBranch        *b_sciModID;   //!
   TBranch        *b_sciEne;   //!
   TBranch        *b_sciTime;   //!
   TBranch        *b_nTrack;   //!
   TBranch        *b_track;   //!
   TBranch        *b_mother;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_ek;   //!
   TBranch        *b_end_ek;   //!
   TBranch        *b_p;   //!
   TBranch        *b_end_p;   //!
   TBranch        *b_v;   //!
   TBranch        *b_end_v;   //!

   E14ReadConvSim(TTree *tree=0);
   virtual ~E14ReadConvSim();
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

#ifdef E14ReadConvSim_cxx
E14ReadConvSim::E14ReadConvSim(TTree *tree) : fChain(0) 
{
   Init(tree);
}

E14ReadConvSim::~E14ReadConvSim()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t E14ReadConvSim::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Int_t E14ReadConvSim::GetEntries(){
  if( !fChain ){ return 0; }
  Int_t Entries = fChain->GetEntries();
  return Entries;
}
Long64_t E14ReadConvSim::LoadTree(Long64_t entry)
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

void E14ReadConvSim::Init(TTree *tree)
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

   fChain->SetBranchAddress("EventNum", &EventNum, &b_EventNum);
   fChain->SetBranchAddress("CsiNumber", &CsiNumber, &b_CsiNumber);
   fChain->SetBranchAddress("CsiTotalE", &CsiTotalE, &b_CsiTotalE);
   fChain->SetBranchAddress("CsiModID", CsiModID, &b_CsiModID);
   fChain->SetBranchAddress("CsiEne", CsiEne, &b_CsiEne);
   fChain->SetBranchAddress("CsiTime", CsiTime, &b_CsiTime);
   fChain->SetBranchAddress("CC03Number", &CC03Number, &b_CC03Number);
   fChain->SetBranchAddress("CC03TotalE", &CC03TotalE, &b_CC03TotalE);
   fChain->SetBranchAddress("CC03ModID", CC03ModID, &b_CC03ModID);
   fChain->SetBranchAddress("CC03Ene", CC03Ene, &b_CC03Ene);
   fChain->SetBranchAddress("CC03Time", CC03Time, &b_CC03Time);
   fChain->SetBranchAddress("CVNumber", &CVNumber, &b_CVNumber);
   fChain->SetBranchAddress("CVTotalE", &CVTotalE, &b_CVTotalE);
   fChain->SetBranchAddress("CVModID", &CVModID, &b_CVModID);
   fChain->SetBranchAddress("CVEne", &CVEne, &b_CVEne);
   fChain->SetBranchAddress("CVTime", &CVTime, &b_CVTime);
   fChain->SetBranchAddress("OEVNumber", &OEVNumber, &b_OEVNumber);
   fChain->SetBranchAddress("OEVTotalE", &OEVTotalE, &b_OEVTotalE);
   fChain->SetBranchAddress("OEVModID", &OEVModID, &b_OEVModID);
   fChain->SetBranchAddress("OEVEne", &OEVEne, &b_OEVEne);
   fChain->SetBranchAddress("OEVTime", &OEVTime, &b_OEVTime);
   fChain->SetBranchAddress("LaserNumber", &LaserNumber, &b_LaserNumber);
   fChain->SetBranchAddress("LaserTotalE", &LaserTotalE, &b_LaserTotalE);
   fChain->SetBranchAddress("LaserModID", &LaserModID, &b_LaserModID);
   fChain->SetBranchAddress("LaserEne", &LaserEne, &b_LaserEne);
   fChain->SetBranchAddress("LaserTime", &LaserTime, &b_LaserTime);
   fChain->SetBranchAddress("CosmicNumber", &CosmicNumber, &b_CosmicNumber);
   fChain->SetBranchAddress("CosmicTotalE", &CosmicTotalE, &b_CosmicTotalE);
   fChain->SetBranchAddress("CosmicModID", &CosmicModID, &b_CosmicModID);
   fChain->SetBranchAddress("CosmicEne", &CosmicEne, &b_CosmicEne);
   fChain->SetBranchAddress("CosmicTime", &CosmicTime, &b_CosmicTime);
   fChain->SetBranchAddress("sciNumber", &sciNumber, &b_sciNumber);
   fChain->SetBranchAddress("sciTotalE", &sciTotalE, &b_sciTotalE);
   fChain->SetBranchAddress("sciModID", sciModID, &b_sciModID);
   fChain->SetBranchAddress("sciEne", sciEne, &b_sciEne);
   fChain->SetBranchAddress("sciTime", sciTime, &b_sciTime);
   fChain->SetBranchAddress("nTrack", &nTrack, &b_nTrack);
   fChain->SetBranchAddress("track", track, &b_track);
   fChain->SetBranchAddress("mother", mother, &b_mother);
   fChain->SetBranchAddress("pid", pid, &b_pid);
   fChain->SetBranchAddress("mass", mass, &b_mass);
   fChain->SetBranchAddress("ek", ek, &b_ek);
   fChain->SetBranchAddress("end_ek", end_ek, &b_end_ek);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("end_p", end_p, &b_end_p);
   fChain->SetBranchAddress("v", v, &b_v);
   fChain->SetBranchAddress("end_v", end_v, &b_end_v);
   Notify();
}

Bool_t E14ReadConvSim::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void E14ReadConvSim::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t E14ReadConvSim::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef E14ReadConvSim_cxx
