//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 27 15:18:01 2012 by ROOT version 5.30/04
// from TTree trCluster/
// found on file: /Volume0/gamma/Data_All.root
//////////////////////////////////////////////////////////

#ifndef ClusterTimeReader_h
#define ClusterTimeReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ClusterTimeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           EventNumber;
   Int_t           nCluster;
   Int_t           nCrystal[120];   //[nCluster]
   Int_t           ClusterID[120];   //[nCluster]
   Double_t        ClusterEnergy[120];   //[nCluster]
   Double_t        ClusterR[120];   //[nCluster]
   Double_t        ClusterT[120];   //[nCluster]
   Double_t        ClusterTheta[120];   //[nCluster]
   Double_t        ClusterPhi[120];   //[nCluster]
   Double_t        CrystalT[120][120];   //[nCluster]
   Double_t        CrystalEnergy[120][120];   //[nCluster]
   Double_t        CrystalR[120][120];   //[nCluster]
   Double_t        CrystalPhi[120][120];   //[nCluster]
   Int_t           CrystalID[120][120];   //[nCluster]
   Double_t        CrystalSignal[120][120];//[nCluster]
   // List of branches
   TBranch        *b_EventNumber;   //!
   TBranch        *b_nCluster;   //!
   TBranch        *b_nCrystal;   //!
   TBranch        *b_ClusterID;   //!
   TBranch        *b_ClusterEnergy;   //!
   TBranch        *b_ClusterR;   //!
   TBranch        *b_ClusterT;   //!
   TBranch        *b_ClusterTheta;   //!
   TBranch        *b_ClusterPhi;   //!
   TBranch        *b_CrystalT;   //!
   TBranch        *b_CrystalEnergy;   //!
   TBranch        *b_CrystalR;   //!
   TBranch        *b_CrystalPhi;   //!
   TBranch        *b_CrystalID;  //!
   TBranch        *b_CrystalSignal; //!

   ClusterTimeReader(TTree *tree=0);
   virtual ~ClusterTimeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ClusterTimeReader_cxx
ClusterTimeReader::ClusterTimeReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

ClusterTimeReader::~ClusterTimeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ClusterTimeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ClusterTimeReader::LoadTree(Long64_t entry)
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

void ClusterTimeReader::Init(TTree *tree)
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

   fChain->SetBranchAddress("EventNumber"  ,&EventNumber , &b_EventNumber);
   fChain->SetBranchAddress("nCluster"     ,&nCluster    , &b_nCluster);
   fChain->SetBranchAddress("nCrystal"     ,nCrystal     , &b_nCrystal);
   fChain->SetBranchAddress("ClusterID"    ,ClusterID    , &b_ClusterID);
   fChain->SetBranchAddress("ClusterEnergy",ClusterEnergy, &b_ClusterEnergy);
   fChain->SetBranchAddress("ClusterR"     ,ClusterR     , &b_ClusterR);
   fChain->SetBranchAddress("ClusterT"     ,ClusterT     , &b_ClusterT);
   fChain->SetBranchAddress("ClusterTheta" ,ClusterTheta , &b_ClusterTheta);
   fChain->SetBranchAddress("ClusterPhi"   ,ClusterPhi   , &b_ClusterPhi);
   fChain->SetBranchAddress("CrystalT"     ,CrystalT     , &b_CrystalT);
   fChain->SetBranchAddress("CrystalEnergy",CrystalEnergy, &b_CrystalEnergy);
   fChain->SetBranchAddress("CrystalR"     ,CrystalR     , &b_CrystalR);
   fChain->SetBranchAddress("CrystalPhi"   ,CrystalPhi   , &b_CrystalPhi);
   fChain->SetBranchAddress("CrystalID"    ,CrystalID    , &b_CrystalID);
   fChain->SetBranchAddress("CrystalSignal",CrystalSignal, &b_CrystalSignal);
   Notify();
}

Bool_t ClusterTimeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ClusterTimeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ClusterTimeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ClusterTimeReader_cxx
