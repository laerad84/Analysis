#ifndef E14CONVREADER__H__
#include "ConvFileAnalysis/E14ConvReader.h"
#include <iostream>
#endif
E14ConvReader::E14ConvReader(){
  m_ch = new TChain("EventTree");
}

E14ConvReader::~E14ConvReader(){
  std::cout  << "Delete E14ConvReader" << std::endl;
  //delete m_ch;
}

bool E14ConvReader::SetBranchAddress(){

  m_ch->SetBranchAddress("EventNo"         ,&EventNo);
  m_ch->SetBranchAddress("nFADC"           ,&nFADC);
  m_ch->SetBranchAddress("nSamples"        ,&nSamples);
  m_ch->SetBranchAddress("Data"            ,Data);//nFADC,nSamples
  m_ch->SetBranchAddress("PeakHeight"      ,PeakHeight);
  m_ch->SetBranchAddress("PeakTime"        ,PeakTime);
  m_ch->SetBranchAddress("Pedestal"        ,Pedestal);
  m_ch->SetBranchAddress("IntegratedADC"   ,IntegratedADC);
  m_ch->SetBranchAddress("EtSum_FADC"      ,EtSum_FADC);
  m_ch->SetBranchAddress("BufferLoopNo"    ,&BufferLoopNo);
  m_ch->SetBranchAddress("Error"           ,Error);
  m_ch->SetBranchAddress("TimeStamp"       ,TimeStamp);
  m_ch->SetBranchAddress("TrigNo"          ,&TrigNo);
  m_ch->SetBranchAddress("SlotNo"          ,&SlotNo);
  m_ch->SetBranchAddress("Compression_flag",&Compression_flag);

  return true;
}

long E14ConvReader::GetChainEntries(){
  return m_ch->GetEntries();
}

int E14ConvReader::GetChainEntry(int ientry){
  return m_ch->GetEntry(ientry);
}

int E14ConvReader::AddFile( const char* filename){
  return m_ch->Add(filename);
}



