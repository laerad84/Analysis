#ifndef __E14ReadConvFile__h__
#define __E14ReadConvFile__h__

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"

static const int  NSLOT_PER_CRATE = 20;
static const int  NCH_PER_SLOT    = 16;
static const int  NPOINT_PER_CH   = 48;

class E14ReadConv{
  
 public:
  TChain* ch;  
  short EventNo;
  short nFADC;
  short nSamples;
  short Data[NSLOT_PER_CRATE][NCH_PER_SLOT][NPOINT_PER_CH];
  short PeakHeight[NSLOT_PER_CRATE][NCH_PER_SLOT];
  short PeakTime[NSLOT_PER_CRATE][NCH_PER_SLOT];  
  float Pedestal[NSLOT_PER_CRATE][NCH_PER_SLOT];
  float IntegratedADC[NSLOT_PER_CRATE][NCH_PER_SLOT];
  int   EtSum_FADC[NSLOT_PER_CRATE][NPOINT_PER_CH];
  long  BufferLoopNo;
  short Error[NSLOT_PER_CRATE];
  int   TimeStamp[NSLOT_PER_CRATE];
  short TrigNo[NSLOT_PER_CRATE];
  short SpillNo[NSLOT_PER_CRATE];
  short Compression_flag[NSLOT_PER_CRATE][NCH_PER_SLOT];
  



  
};



    
