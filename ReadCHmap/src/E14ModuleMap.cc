#define E14ModuleMap_cxx
#include "E14ModuleMap.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

E14ModuleMap::E14ModuleMap(const char* name, TTree *tree)
{
  fName = name;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     std::cout<< "Please Set Tree" << std::endl;
   }
   Init(tree);
}

E14ModuleMap::~E14ModuleMap()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t E14ModuleMap::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t E14ModuleMap::LoadTree(Long64_t entry)
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

void E14ModuleMap::Init(TTree *tree)
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

   fChain->SetBranchAddress(Form("%s.TObject.fUniqueID",fName.c_str()), &TObject_fUniqueID, &b_TObject_fUniqueID);
   fChain->SetBranchAddress(Form("%s.TObject.fBits"    ,fName.c_str()), &TObject_fBits, &b_TObject_fBits);
   fChain->SetBranchAddress(Form("%s.ThisID"           ,fName.c_str()), &ThisID, &b_ThisID);
   fChain->SetBranchAddress(Form("%s.AllNum"           ,fName.c_str()), &AllNum, &b_AllNum);
   fChain->SetBranchAddress(Form("%s.DataStoreTag"     ,fName.c_str()), &DataStoreTag, &b_DataStoreTag);
   fChain->SetBranchAddress(Form("%s.CrateID[4096]"    ,fName.c_str()), CrateID, &b_CrateID);
   fChain->SetBranchAddress(Form("%s.FADCID[4096]"     ,fName.c_str()), FADCID, &b_FADCID);
   fChain->SetBranchAddress(Form("%s.CHID[4096]"       ,fName.c_str()), CHID, &b_CHID);
   Notify();
}

Bool_t E14ModuleMap::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   return kTRUE;
}

void E14ModuleMap::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t E14ModuleMap::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void E14ModuleMap::PrintAll(){
  for( int i = 0; i< AllNum; i++){
    std::cout << CrateID[i] << " : " << FADCID[i] << " : " << CHID[i] << std::endl;
  }
}

void E14ModuleMap::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;
   }
}

bool E14ModuleMap::GetCFC(int modCHID, int& crateID, int& fadcID, int& chID){
  if( modCHID >= AllNum || modCHID < 0 ){ return false; }
  crateID = CrateID[modCHID];
  fadcID  = FADCID[modCHID];
  chID    = CHID[modCHID];
  return true;
}
