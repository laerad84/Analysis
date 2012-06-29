#ifndef E14CONVREADER__H__
#define E14CONVREADER__H__

#include "TROOT.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

class E14ConvReader 
{
 private:
 public:

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


  TTree* m_ch;

 public:
  /*
  E14ConvReader(TFile*);
  ~E14ConvReader();
  
  virtual bool SetBranchAddress();
  virtual long GetEntries();
  virtual int  GetEntry(int ientry);
  virtual int  AddFile( const char *filename);
  */

   E14ConvReader(TTree *tree=0);
   virtual ~E14ConvReader();
   virtual Int_t    GetEntries();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
