//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul  2 14:35:39 2013 by ROOT version 5.30/04
// from TTree T/
// found on file: Conv_KL3pi0_FAST_REDUCED_5E6_0.root
//////////////////////////////////////////////////////////

#ifndef ConvReader_h
#define ConvReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ConvReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EventNum;
   Int_t           CsiNumber;
   Double_t        CsiTotalE;
   Int_t           CsiModID[502];   //[CsiNumber]
   Double_t        CsiEne[502];   //[CsiNumber]
   Double_t        CsiTime[502];   //[CsiNumber]
   Int_t           CC03Number;
   Double_t        CC03TotalE;
   Int_t           CC03ModID[16];   //[CC03Number]
   Double_t        CC03Ene[16];   //[CC03Number]
   Double_t        CC03Time[16];   //[CC03Number]
   Int_t           CVNumber;
   Double_t        CVTotalE;
   Int_t           CVModID[1];   //[CVNumber]
   Double_t        CVEne[1];   //[CVNumber]
   Double_t        CVTime[1];   //[CVNumber]
   Int_t           OEVNumber;
   Double_t        OEVTotalE;
   Int_t           OEVModID[8];   //[OEVNumber]
   Double_t        OEVEne[8];   //[OEVNumber]
   Double_t        OEVTime[8];   //[OEVNumber]
   Int_t           LaserNumber;
   Double_t        LaserTotalE;
   Int_t           LaserModID[1];   //[LaserNumber]
   Double_t        LaserEne[1];   //[LaserNumber]
   Double_t        LaserTime[1];   //[LaserNumber]
   Int_t           CosmicNumber;
   Double_t        CosmicTotalE;
   Int_t           CosmicModID[1];   //[CosmicNumber]
   Double_t        CosmicEne[1];   //[CosmicNumber]
   Double_t        CosmicTime[1];   //[CosmicNumber]
   Int_t           SciNumber;
   Double_t        SciTotalE;
   Int_t           SciModID[1];   //[SciNumber]
   Double_t        SciEne[1];   //[SciNumber]
   Double_t        SciTime[1];   //[SciNumber]
   Int_t           nTrack;
   UShort_t        track[15];   //[nTrack]
   Short_t         mother[15];   //[nTrack]
   Int_t           pid[15];   //[nTrack]
   Float_t         mass[15];   //[nTrack]
   Float_t         ek[15];   //[nTrack]
   Float_t         end_ek[15];   //[nTrack]
   Double_t        p[15][3];   //[nTrack]
   Double_t        end_p[15][3];   //[nTrack]
   Double_t        v[15][3];   //[nTrack]
   Double_t        end_v[15][3];   //[nTrack]

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
   TBranch        *b_SciNumber;   //!
   TBranch        *b_SciTotalE;   //!
   TBranch        *b_SciModID;   //!
   TBranch        *b_SciEne;   //!
   TBranch        *b_SciTime;   //!
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

   ConvReader(TTree *tree=0);
   virtual ~ConvReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ConvReader_cxx
ConvReader::ConvReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Conv_KL3pi0_FAST_REDUCED_5E6_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Conv_KL3pi0_FAST_REDUCED_5E6_0.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

ConvReader::~ConvReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ConvReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ConvReader::LoadTree(Long64_t entry)
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

void ConvReader::Init(TTree *tree)
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
   fChain->SetBranchAddress("OEVModID", OEVModID, &b_OEVModID);
   fChain->SetBranchAddress("OEVEne", OEVEne, &b_OEVEne);
   fChain->SetBranchAddress("OEVTime", OEVTime, &b_OEVTime);
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
   fChain->SetBranchAddress("SciNumber", &SciNumber, &b_SciNumber);
   fChain->SetBranchAddress("SciTotalE", &SciTotalE, &b_SciTotalE);
   fChain->SetBranchAddress("SciModID", &SciModID, &b_SciModID);
   fChain->SetBranchAddress("SciEne", &SciEne, &b_SciEne);
   fChain->SetBranchAddress("SciTime", &SciTime, &b_SciTime);
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

Bool_t ConvReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ConvReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ConvReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ConvReader_cxx
