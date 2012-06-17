#ifndef E14RawData_h
#define E14RawData_h

#include "TObject.h"

#include "TFile.h"
#include "TTree.h"

class E14RawData
{
 public:
  
  E14RawData();
  virtual ~E14RawData();
  
  virtual void   Clear(Option_t* opt="");
  
  void SetRunNum( int val );
  void initializeDataValues();
  
  void dump(int i, int j);
  
 public:
  //storage
  short nFADC;
  short nSamples; 
  int   EventNo;
  
  short Data[20][16][48];
  float Pedestal[20][16];
  short PeakHeight[20][16];
  short PeakTime[20][16];
  float IntegratedADC[20][16];
  int   EtSum_FADC[20][48];
  
  short Error[20];
  int   TimeStamp[20];
  short TrigNo[20];
  short SpillNo[20];
  short SlotNo[20];
  short Compression_flag[20][16];
  long  BufferLoopNo;

  TFile *hfile;
  TTree *tree;
  
 public: 
  virtual void Fill();
  virtual void Write();
  virtual int CalcAll();
  virtual void InitTree( char* s );
  
 private:
    
  virtual int CalcPedestal();
  virtual int CalcPeakHeight();
  virtual int CalcIntegratedADC();
  virtual int CalcEtSumFADC();

  int RunNum;

  std::string o_rootFile;
  
  ClassDef(E14RawData,1)

    };

#endif // E14RawData_h
