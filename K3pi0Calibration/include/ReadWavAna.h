//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov  6 13:09:40 2012 by ROOT version 5.30/04
// from TTree T/Output from Time zero
// found on file: /Volume0/ExpData/2012_Feb_Beam/RootFile_wav/run_wav_4200_cl.root
//////////////////////////////////////////////////////////

#ifndef ReadWavAna_h
#define ReadWavAna_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ReadWavAna {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           RunNumber;
   Int_t           EventNumber;
   Int_t           CsiNumber;
   Int_t           CsiModID[280];   //[CsiNumber]
   Double_t        CsiEne[280];   //[CsiNumber]
   Double_t        CsiTime[280];   //[CsiNumber]
   Double_t        CsiHHTime[280];   //[CsiNumber]
   Int_t           eventID;
   Int_t           OrigEventID;
   Int_t           CutCondition;
   Int_t           VetoCondition;
   Int_t           ClusterNumber;
   Int_t           ClusterId[10];   //[ClusterNumber]
   Int_t           ClusterStatus[10];   //[ClusterNumber]
   Double_t        ClusterThreshold[10];   //[ClusterNumber]
   Double_t        ClusterDepE[10];   //[ClusterNumber]
   Double_t        ClusterCoePos[10][3];   //[ClusterNumber]
   Double_t        ClusterTime[10];   //[ClusterNumber]
   Double_t        ClusterRMS[10];   //[ClusterNumber]
   Int_t           ClusterSize[10];   //[ClusterNumber]
   Int_t           ClusterCsiId[10][120];   //[ClusterNumber]
   Double_t        ClusterCsiE[10][120];   //[ClusterNumber]
   Double_t        ClusterCsiTime[10][120];   //[ClusterNumber]
   Int_t           GamClusNumber;
   Int_t           GamClusId[7];   //[GamClusNumber]
   Int_t           GamClusStatus[7];   //[GamClusNumber]
   Double_t        GamClusThreshold[7];   //[GamClusNumber]
   Double_t        GamClusDepE[7];   //[GamClusNumber]
   Double_t        GamClusCoePos[7][3];   //[GamClusNumber]
   Double_t        GamClusTime[7];   //[GamClusNumber]
   Double_t        GamClusRMS[7];   //[GamClusNumber]
   Int_t           GamClusSize[7];   //[GamClusNumber]
   Int_t           GamClusCsiId[7][120];   //[GamClusNumber]
   Double_t        GamClusCsiE[7][120];   //[GamClusNumber]
   Double_t        GamClusCsiTime[7][120];   //[GamClusNumber]
   Int_t           GammaNumber;
   Int_t           GammaId[12];   //[GammaNumber]
   Int_t           GammaStatus[12];   //[GammaNumber]
   Double_t        GammaE[12];   //[GammaNumber]
   Double_t        GammaPos[12][3];   //[GammaNumber]
   Double_t        GammaTime[12];   //[GammaNumber]
   Double_t        GammaMom[12][3];   //[GammaNumber]
   Double_t        GammaSigmaE[12];   //[GammaNumber]
   Double_t        GammaSigmaPos[12][3];   //[GammaNumber]
   Double_t        GammaChi2[12];   //[GammaNumber]
   Double_t        GammaAnn[12];   //[GammaNumber]
   Int_t           Gamma_clusIndex[12];   //[GammaNumber]
   Int_t           Pi0Number;
   Int_t           Pi0Id[6];   //[Pi0Number]
   Int_t           Pi0Status[6];   //[Pi0Number]
   Double_t        Pi0E[6];   //[Pi0Number]
   Double_t        Pi0Pos[6][3];   //[Pi0Number]
   Double_t        Pi0Mom[6][3];   //[Pi0Number]
   Double_t        Pi0Pt[6];   //[Pi0Number]
   Double_t        Pi0Mass[6];   //[Pi0Number]
   Double_t        Pi0RecZ[6];   //[Pi0Number]
   Double_t        Pi0RecZsig2[6];   //[Pi0Number]
   Int_t           Pi0_gamIndex[6][2];   //[Pi0Number]
   Int_t           KlongNumber;
   Int_t           KlongId[2];   //[KlongNumber]
   Int_t           KlongStatus[2];   //[KlongNumber]
   Double_t        KlongE[2];   //[KlongNumber]
   Double_t        KlongPos[2][3];   //[KlongNumber]
   Double_t        KlongMom[2][3];   //[KlongNumber]
   Double_t        KlongPt[2];   //[KlongNumber]
   Double_t        KlongMass[2];   //[KlongNumber]
   Double_t        KlongDeltaZ[2];   //[KlongNumber]
   Double_t        KlongChisqZ[2];   //[KlongNumber]
   Int_t           KlongVertFlag[2];   //[KlongNumber]
   Int_t           KlongSortFlag[2];   //[KlongNumber]
   Int_t           Klong_piIndex[2][3];   //[KlongNumber]

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_CsiNumber;   //!
   TBranch        *b_CsiModID;   //!
   TBranch        *b_CsiEne;   //!
   TBranch        *b_CsiTime;   //!
   TBranch        *b_CsiHHTime;   //!
   TBranch        *b_eventID;   //!
   TBranch        *b_OrigEventID;   //!
   TBranch        *b_CutCondition;   //!
   TBranch        *b_VetoCondition;   //!
   TBranch        *b_ClusterNumber;   //!
   TBranch        *b_ClusterId;   //!
   TBranch        *b_ClusterStatus;   //!
   TBranch        *b_ClusterThreshold;   //!
   TBranch        *b_ClusterDepE;   //!
   TBranch        *b_ClusterCoePos;   //!
   TBranch        *b_ClusterTime;   //!
   TBranch        *b_ClusterRMS;   //!
   TBranch        *b_ClusterSize;   //!
   TBranch        *b_ClusterCsiId;   //!
   TBranch        *b_ClusterCsiE;   //!
   TBranch        *b_ClusterCsiTime;   //!
   TBranch        *b_GamClusNumber;   //!
   TBranch        *b_GamClusId;   //!
   TBranch        *b_GamClusStatus;   //!
   TBranch        *b_GamClusThreshold;   //!
   TBranch        *b_GamClusDepE;   //!
   TBranch        *b_GamClusCoePos;   //!
   TBranch        *b_GamClusTime;   //!
   TBranch        *b_GamClusRMS;   //!
   TBranch        *b_GamClusSize;   //!
   TBranch        *b_GamClusCsiId;   //!
   TBranch        *b_GamClusCsiE;   //!
   TBranch        *b_GamClusCsiTime;   //!
   TBranch        *b_GammaNumber;   //!
   TBranch        *b_GammaId;   //!
   TBranch        *b_GammaStatus;   //!
   TBranch        *b_GammaE;   //!
   TBranch        *b_GammaPos;   //!
   TBranch        *b_GammaTime;   //!
   TBranch        *b_GammaMom;   //!
   TBranch        *b_GammaSigmaE;   //!
   TBranch        *b_GammaSigmaPos;   //!
   TBranch        *b_GammaChi2;   //!
   TBranch        *b_GammaAnn;   //!
   TBranch        *b_Gamma_clusIndex;   //!
   TBranch        *b_Pi0Number;   //!
   TBranch        *b_Pi0Id;   //!
   TBranch        *b_Pi0Status;   //!
   TBranch        *b_Pi0E;   //!
   TBranch        *b_Pi0Pos;   //!
   TBranch        *b_Pi0Mom;   //!
   TBranch        *b_Pi0Pt;   //!
   TBranch        *b_Pi0Mass;   //!
   TBranch        *b_Pi0RecZ;   //!
   TBranch        *b_Pi0RecZsig2;   //!
   TBranch        *b_Pi0_gamIndex;   //!
   TBranch        *b_KlongNumber;   //!
   TBranch        *b_KlongId;   //!
   TBranch        *b_KlongStatus;   //!
   TBranch        *b_KlongE;   //!
   TBranch        *b_KlongPos;   //!
   TBranch        *b_KlongMom;   //!
   TBranch        *b_KlongPt;   //!
   TBranch        *b_KlongMass;   //!
   TBranch        *b_KlongDeltaZ;   //!
   TBranch        *b_KlongChisqZ;   //!
   TBranch        *b_KlongVertFlag;   //!
   TBranch        *b_KlongSortFlag;   //!
   TBranch        *b_Klong_piIndex;   //!

   ReadWavAna(TTree *tree=0);
   virtual ~ReadWavAna();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ReadWavAna_cxx
ReadWavAna::ReadWavAna(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  Init(tree);
}

ReadWavAna::~ReadWavAna()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ReadWavAna::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ReadWavAna::LoadTree(Long64_t entry)
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

void ReadWavAna::Init(TTree *tree)
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

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("CsiNumber", &CsiNumber, &b_CsiNumber);
   fChain->SetBranchAddress("CsiModID", CsiModID, &b_CsiModID);
   fChain->SetBranchAddress("CsiEne", CsiEne, &b_CsiEne);
   fChain->SetBranchAddress("CsiTime", CsiTime, &b_CsiTime);
   fChain->SetBranchAddress("CsiHHTime", CsiHHTime, &b_CsiHHTime);
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("OrigEventID", &OrigEventID, &b_OrigEventID);
   fChain->SetBranchAddress("CutCondition", &CutCondition, &b_CutCondition);
   fChain->SetBranchAddress("VetoCondition", &VetoCondition, &b_VetoCondition);
   fChain->SetBranchAddress("ClusterNumber", &ClusterNumber, &b_ClusterNumber);
   fChain->SetBranchAddress("ClusterId", ClusterId, &b_ClusterId);
   fChain->SetBranchAddress("ClusterStatus", ClusterStatus, &b_ClusterStatus);
   fChain->SetBranchAddress("ClusterThreshold", ClusterThreshold, &b_ClusterThreshold);
   fChain->SetBranchAddress("ClusterDepE", ClusterDepE, &b_ClusterDepE);
   fChain->SetBranchAddress("ClusterCoePos", ClusterCoePos, &b_ClusterCoePos);
   fChain->SetBranchAddress("ClusterTime", ClusterTime, &b_ClusterTime);
   fChain->SetBranchAddress("ClusterRMS", ClusterRMS, &b_ClusterRMS);
   fChain->SetBranchAddress("ClusterSize", ClusterSize, &b_ClusterSize);
   fChain->SetBranchAddress("ClusterCsiId", ClusterCsiId, &b_ClusterCsiId);
   fChain->SetBranchAddress("ClusterCsiE", ClusterCsiE, &b_ClusterCsiE);
   fChain->SetBranchAddress("ClusterCsiTime", ClusterCsiTime, &b_ClusterCsiTime);
   fChain->SetBranchAddress("GamClusNumber", &GamClusNumber, &b_GamClusNumber);
   fChain->SetBranchAddress("GamClusId", GamClusId, &b_GamClusId);
   fChain->SetBranchAddress("GamClusStatus", GamClusStatus, &b_GamClusStatus);
   fChain->SetBranchAddress("GamClusThreshold", GamClusThreshold, &b_GamClusThreshold);
   fChain->SetBranchAddress("GamClusDepE", GamClusDepE, &b_GamClusDepE);
   fChain->SetBranchAddress("GamClusCoePos", GamClusCoePos, &b_GamClusCoePos);
   fChain->SetBranchAddress("GamClusTime", GamClusTime, &b_GamClusTime);
   fChain->SetBranchAddress("GamClusRMS", GamClusRMS, &b_GamClusRMS);
   fChain->SetBranchAddress("GamClusSize", GamClusSize, &b_GamClusSize);
   fChain->SetBranchAddress("GamClusCsiId", GamClusCsiId, &b_GamClusCsiId);
   fChain->SetBranchAddress("GamClusCsiE", GamClusCsiE, &b_GamClusCsiE);
   fChain->SetBranchAddress("GamClusCsiTime", GamClusCsiTime, &b_GamClusCsiTime);
   fChain->SetBranchAddress("GammaNumber", &GammaNumber, &b_GammaNumber);
   fChain->SetBranchAddress("GammaId", GammaId, &b_GammaId);
   fChain->SetBranchAddress("GammaStatus", GammaStatus, &b_GammaStatus);
   fChain->SetBranchAddress("GammaE", GammaE, &b_GammaE);
   fChain->SetBranchAddress("GammaPos", GammaPos, &b_GammaPos);
   fChain->SetBranchAddress("GammaTime", GammaTime, &b_GammaTime);
   fChain->SetBranchAddress("GammaMom", GammaMom, &b_GammaMom);
   fChain->SetBranchAddress("GammaSigmaE", GammaSigmaE, &b_GammaSigmaE);
   fChain->SetBranchAddress("GammaSigmaPos", GammaSigmaPos, &b_GammaSigmaPos);
   fChain->SetBranchAddress("GammaChi2", GammaChi2, &b_GammaChi2);
   fChain->SetBranchAddress("GammaAnn", GammaAnn, &b_GammaAnn);
   fChain->SetBranchAddress("Gamma_clusIndex", Gamma_clusIndex, &b_Gamma_clusIndex);
   fChain->SetBranchAddress("Pi0Number", &Pi0Number, &b_Pi0Number);
   fChain->SetBranchAddress("Pi0Id", Pi0Id, &b_Pi0Id);
   fChain->SetBranchAddress("Pi0Status", Pi0Status, &b_Pi0Status);
   fChain->SetBranchAddress("Pi0E", Pi0E, &b_Pi0E);
   fChain->SetBranchAddress("Pi0Pos", Pi0Pos, &b_Pi0Pos);
   fChain->SetBranchAddress("Pi0Mom", Pi0Mom, &b_Pi0Mom);
   fChain->SetBranchAddress("Pi0Pt", Pi0Pt, &b_Pi0Pt);
   fChain->SetBranchAddress("Pi0Mass", Pi0Mass, &b_Pi0Mass);
   fChain->SetBranchAddress("Pi0RecZ", Pi0RecZ, &b_Pi0RecZ);
   fChain->SetBranchAddress("Pi0RecZsig2", Pi0RecZsig2, &b_Pi0RecZsig2);
   fChain->SetBranchAddress("Pi0_gamIndex", Pi0_gamIndex, &b_Pi0_gamIndex);
   fChain->SetBranchAddress("KlongNumber", &KlongNumber, &b_KlongNumber);
   fChain->SetBranchAddress("KlongId", KlongId, &b_KlongId);
   fChain->SetBranchAddress("KlongStatus", KlongStatus, &b_KlongStatus);
   fChain->SetBranchAddress("KlongE", KlongE, &b_KlongE);
   fChain->SetBranchAddress("KlongPos", KlongPos, &b_KlongPos);
   fChain->SetBranchAddress("KlongMom", KlongMom, &b_KlongMom);
   fChain->SetBranchAddress("KlongPt", KlongPt, &b_KlongPt);
   fChain->SetBranchAddress("KlongMass", KlongMass, &b_KlongMass);
   fChain->SetBranchAddress("KlongDeltaZ", KlongDeltaZ, &b_KlongDeltaZ);
   fChain->SetBranchAddress("KlongChisqZ", KlongChisqZ, &b_KlongChisqZ);
   fChain->SetBranchAddress("KlongVertFlag", KlongVertFlag, &b_KlongVertFlag);
   fChain->SetBranchAddress("KlongSortFlag", KlongSortFlag, &b_KlongSortFlag);
   fChain->SetBranchAddress("Klong_piIndex", Klong_piIndex, &b_Klong_piIndex);
   Notify();
}

Bool_t ReadWavAna::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ReadWavAna::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ReadWavAna::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ReadWavAna_cxx
