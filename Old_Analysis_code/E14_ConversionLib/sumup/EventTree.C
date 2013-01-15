#include "EventTree.h"

#if !defined(__CINT__)
ClassImp(EventTree)
#endif

EventTree::EventTree()
{
  ;
}

EventTree::~EventTree()
{
  ;
}

void EventTree::Loop()
{

   int jentry = 0;
   tree->GetEntry(jentry);
   
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

  return 1;
}

Int_t EventTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.                                                      
  if (!tree) return 0;
  return tree->GetEntry(entry);
}

Int_t EventTree::GetEntriesFast()
{
  if (!tree) return 0;
  return tree->GetEntriesFast();
}

void EventTree::Init(TTree *tree)
{
  if (!tree) return;

  tree->SetBranchAddress("nFADC", &nFADC);
  tree->SetBranchAddress("nSamples", &nSamples);
  tree->SetBranchAddress("EventNo", &EventNo);
  tree->SetBranchAddress("Data", Data);
  tree->SetBranchAddress("Pedestal", Pedestal);
  tree->SetBranchAddress("PeakHeight", PeakHeight);
  tree->SetBranchAddress("PeakTime", PeakTime);
  tree->SetBranchAddress("IntegratedADC", IntegratedADC);
  tree->SetBranchAddress("EtSum_FADC", EtSum_FADC);
  tree->SetBranchAddress("BufferLoopNo", &BufferLoopNo);
  tree->SetBranchAddress("Error", Error);
  tree->SetBranchAddress("TimeStamp", TimeStamp);
  tree->SetBranchAddress("TrigNo", TrigNo);
  tree->SetBranchAddress("SpillNo", SpillNo);
  tree->SetBranchAddress("SlotNo", SlotNo);
  tree->SetBranchAddress("Compression_flag", Compression_flag);

  tree->SetBranchStatus("*",0);  // disable all branches                  
  tree->SetBranchStatus("nFADC",1);
  //  tree->SetBranchStatus("nSamples",1);
  tree->SetBranchStatus("EventNo",1);
  //  tree->SetBranchStatus("Data",1);
  tree->SetBranchStatus("PeakHeight",1);
  tree->SetBranchStatus("PeakTime",1);
  tree->SetBranchStatus("IntegratedADC",1);
  tree->SetBranchStatus("TimeStamp",1);
  tree->SetBranchStatus("SpillNo", SpillNo);
  tree->SetBranchStatus("Error",1);
  tree->SetBranchStatus("EtSum_FADC",1);

}
