#ifndef EventTree_h
#define EventTree_h

#include "TObject.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

class EventTree 
{
 public :
  
  EventTree();
  virtual ~EventTree();
  
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Int_t    GetEntriesFast();
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Int_t    FileOpen( char *InFileName );

 public:

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
   
 private:
   TTree *tree;

   ClassDef(EventTree,1)

};

#endif // EventTree_h
