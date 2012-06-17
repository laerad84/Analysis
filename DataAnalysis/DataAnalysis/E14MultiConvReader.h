#ifndef E14MULTICONVREADER__H__
#define E14MULTICONVREADER__H__
static const int N_MAXIMUM_CRATE_NUMBER = 20;

#include "DataAnalysis/E14ConvReader.h"
#include "DataAnalysis/E14IDHandler.h"

class E14MultiConvReader{
 private:

  short crateNo;
  short slotNo;
  short channelNo;

  int    NumberOfCrate;
  short* Data[N_MAXIMUM_CRATE_NUMBER];
  short* PeakHeight[N_MAXIMUM_CRATE_NUMBER];
  short* PeakTime[N_MAXIMUM_CRATE_NUMBER];
  float* Pedestal[N_MAXIMUM_CRATE_NUMBER];
  float* IntegratedADC[N_MAXIMUM_CRATE_NUMBER];
  short* Error[N_MAXIMUM_CRATE_NUMBER];  


 public:
  short* CsIData[E14_NCsI];
  E14ConvReader* conv[N_MAXIMUM_CRATE_NUMBER];
  E14IDHandler * IDhandler;

  E14MultiConvReader(int nCrate = 11);
  ~E14MultiConvReader();

  virtual bool  Init();
  virtual bool  SetAddress();
  virtual int   GetChainEntry(int ientry);
  virtual long  GetChainEntries();
  virtual int   AddFile( const char *ConvFileDir,int RunNumber);
  
};

#endif //E14MULTICONVREADER__H__
