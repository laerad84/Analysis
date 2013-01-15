#ifndef E14DataReadVME_h
#define E14DataReadVME_h

#include "E14RawData.h"

#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>

#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TString.h"
#include <cstdlib>
#include <cstdio>

#include <iostream>
#include <sstream>
#include <iomanip>

class E14DataReadVME
{
 public:

  E14DataReadVME();
  virtual ~E14DataReadVME();
  
  virtual void  SetDebugMode() { DebugMode = true; }
  virtual void  SetInputFile(  std::string s) { i_rootFile = s; }
  virtual void  SetInputDirectory(  std::string s) { id_rootFile = s; }
  virtual void  SetOutputFile(  std::string s) { o_rootFile = s; }
  virtual void  SetRunNumber(  int num ) { RunNumber = num; }

  /* Raw data get from binary file */
  virtual void  FileOpen( char *InFileName );
  virtual int   HeaderRead();
  virtual int   GetBuffer();
  virtual int   DataRead( int sn );

  /* Raw data get from conv file */
  virtual int  FileOpenConv( char *InFileName );
  virtual void  DataReadConv( int ievt );

  /* Data accessor for storage */
  virtual short GetCrateID(){ return (short)CrateID; }
  virtual int   GetEventNo(){ return EventNo; }
  virtual short GetnFADC(){ return nFADC; }
  virtual short GetnSamples(){ return nSamples; }

  virtual short GetData( int slot_n, int ch_n, int sample_n){ return Data[slot_n][ch_n][sample_n]; }
  virtual short GetError( int slot_n ){ return Error[slot_n]; }
  virtual int   GetTimeStamp( int slot_n ){ return TimeStamp[slot_n]; }
  virtual short GetTrigNo( int slot_n ){ return TrigNo[slot_n]; }
  virtual int   GetSpillNo( int slot_n ){ return SpillNo[slot_n]; }
  virtual short GetSlotNo( int slot_n ){ return SlotNo[slot_n]; }
  virtual int   GetCompressionFlag( int slot_n, int ch_n ){ return Compression_flag[slot_n][ch_n]; }

 private:  
  virtual void initializeDataValues();
  virtual void TimeStampJudge();

 private: /* for data storage */
  int CrateID;

  int nFADC;
  int nSamples;
  int EventNo;

  int Data[20][16][48];
  int Error[20];
  int TimeStamp[20];
  int TrigNo[20];
  int SpillNo[20];
  int SlotNo[20];
  int Compression_flag[20][16];

 private:

  E14RawData rawdata[6];

  char InFileName[256];
  std::string i_rootFile;
  std::string id_rootFile;
  std::string o_rootFile;

  int RunNumber;

  FILE *fout;
  unsigned int dp;
  int nCrate;

  TFile *m_file;
  TTree *m_tree;

  int bufSize;
  unsigned short bufData[20][35000];

  int ich,count;
  int eventhead,amthead,hodohead,adchead;
  int nread,len;
  int readend;
  int stoptime,hitnumber,hitch,hittime;

  bool DebugMode;

  ClassDef(E14DataReadVME,1)
    };

#endif // E14DataReadVME_h
