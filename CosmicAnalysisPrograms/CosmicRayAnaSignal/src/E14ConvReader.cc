#ifndef E14CONVREADER__H__
#include "E14ConvReader.h"
#include <iostream>
#endif

/*
E14ConvReader::E14ConvReader(TFile* tf){
  //m_ch = NULL;
  //m_ch = new TChain("EventTree");
  m_ch= (TTree*)tf->Get("EventTree");
}

E14ConvReader::~E14ConvReader(){
  std::cout  << "Delete E14ConvReader" << std::endl;
  m_ch->Clear();
  m_ch = NULL;
}

bool E14ConvReader::SetBranchAddress(){

  m_ch->SetBranchAddress("EventNo"         ,&EventNo);
  m_ch->SetBranchAddress("nFADC"           ,&nFADC);
  m_ch->SetBranchAddress("nSamples"        ,&nSamples);
  m_ch->SetBranchAddress("Data"            ,Data);//nFADC,nSamples
  m_ch->SetBranchAddress("PeakHeight"      ,PeakHeight);//nFADC
  m_ch->SetBranchAddress("PeakTime"        ,PeakTime);
  m_ch->SetBranchAddress("Pedestal"        ,Pedestal);
  m_ch->SetBranchAddress("IntegratedADC"   ,IntegratedADC);
  m_ch->SetBranchAddress("EtSum_FADC"      ,EtSum_FADC);
  //m_ch->SetBranchAddress("BufferLoopNo"    ,&BufferLoopNo);
  m_ch->SetBranchAddress("Error"           ,Error);
  m_ch->SetBranchAddress("TimeStamp"       ,TimeStamp);
  m_ch->SetBranchAddress("TrigNo"          ,&TrigNo);
  m_ch->SetBranchAddress("SlotNo"          ,&SlotNo);
  m_ch->SetBranchAddress("Compression_flag",&Compression_flag);

  return true;
}

long E14ConvReader::GetEntries(){
  return m_ch->GetEntries();
}

int E14ConvReader::GetEntry(int ientry){
  return m_ch->GetEntry(ientry);
}

int E14ConvReader::AddFile( const char* filename){
  //return m_ch->Add(filename);
  return 0;
}

*/

E14ConvReader::E14ConvReader(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   Init(tree);
}

E14ConvReader::~E14ConvReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t E14ConvReader::GetEntries(){
  if( !fChain ) return 0;
  return fChain->GetEntries();
}

Int_t E14ConvReader::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void E14ConvReader::Init(TTree *tree)
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
   fChain->SetBranchAddress("nFADC", &nFADC, &b_nFADC);
   fChain->SetBranchAddress("nSamples", &nSamples, &b_nSamples);
   fChain->SetBranchAddress("EventNo", &EventNo, &b_EventNo);
   fChain->SetBranchAddress("Data", Data, &b_Data);
   fChain->SetBranchAddress("Pedestal", Pedestal, &b_Pedestal);
   fChain->SetBranchAddress("PeakHeight", PeakHeight, &b_PeakHeight);
   fChain->SetBranchAddress("PeakTime", PeakTime, &b_PeakTime);
   fChain->SetBranchAddress("IntegratedADC", IntegratedADC, &b_IntegratedADC);
   fChain->SetBranchAddress("EtSum_FADC", EtSum_FADC, &b_EtSum_FADC);
   fChain->SetBranchAddress("BufferLoopNo", &BufferLoopNo, &b_BufferLoopNo);
   fChain->SetBranchAddress("Error", Error, &b_Error);
   fChain->SetBranchAddress("TimeStamp", TimeStamp, &b_TimeStamp);
   fChain->SetBranchAddress("TrigNo", TrigNo, &b_TrigNo);
   fChain->SetBranchAddress("SpillNo", SpillNo, &b_SpillNo);
   fChain->SetBranchAddress("SlotNo", SlotNo, &b_SlotNo);
   fChain->SetBranchAddress("Compression_flag", Compression_flag, &b_Compression_flag);
   Notify();
}

Bool_t E14ConvReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void E14ConvReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
