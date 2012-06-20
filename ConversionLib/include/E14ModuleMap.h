//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun  5 15:46:34 2012 by ROOT version 5.28/00b
// from TTree map/Mapping data
// found on file: Sum4557.root
//////////////////////////////////////////////////////////

#ifndef E14ModuleMap_h
#define E14ModuleMap_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <cstring>
#include <iostream>

const Int_t kMaxModule = 1;
class E14ModuleMap {
 public :
  std::string     fName;    // Name of Module
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Declaration of leaf types
  //E14Mapper       *Csi_;
  UInt_t          TObject_fUniqueID;
  UInt_t          TObject_fBits;
  Int_t           ThisID;
  Int_t           AllNum;
  Bool_t          DataStoreTag;
  Int_t           CrateID[4096];
  Int_t           FADCID[4096];
  Int_t           CHID[4096];
  // List of branches
  TBranch        *b_TObject_fUniqueID;   //!
  TBranch        *b_TObject_fBits;   //!
  TBranch        *b_ThisID;   //!
  TBranch        *b_AllNum;   //!
  TBranch        *b_DataStoreTag;   //!
  TBranch        *b_CrateID;   //!
  TBranch        *b_FADCID;   //!
  TBranch        *b_CHID;   //!
  
  E14ModuleMap(const char*, TTree *tree=0);
  virtual ~E14ModuleMap();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
  virtual Bool_t   GetCFC(int,int&,int&,int&);
  virtual void     PrintAll();

  virtual std::string GetName(){ return fName; }
  virtual Int_t       GetAllNum(){return AllNum;}
};

#endif

