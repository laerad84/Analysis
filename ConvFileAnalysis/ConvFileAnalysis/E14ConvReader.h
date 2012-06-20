#ifndef E14CONVREADER__H__
#define E14CONVREADER__H__

#include "TChain.h"

class E14ConvReader 
{
 private:
  TChain* m_ch;
 public:
  int   EventNo;
  short nFADC;
  short nSamples;
  short Data[20][16][48];
  short PeakHeight[20][16];
  short PeakTime[20][16];
  float Pedestal[20][16];
  float IntegratedADC[20][16];
  int   EtSum_FADC[20][48];
  long  BufferLoopNo;
  short Error[20];
  int   TimeStamp[20];
  short TrigNo[20];
  short SpillNo[20];
  short SlotNo[20];
  short Compression_flag;
  
 public:
  E14ConvReader();
  ~E14ConvReader();
  
  virtual bool SetBranchAddress();
  virtual long GetChainEntries();
  virtual int  GetChainEntry(int ientry);
  virtual int  AddFile( const char *filename);
  
};
#endif
