#ifndef __E14ReadSumFile__h__
#define __E14ReadSumFile__h__
#include "TChain.h"
#include "TFile.h"

static const int N_TOTAL_CC03   = 32;
static const int N_TOTAL_OEV    = 44;
static const int N_TOTAL_CV     = 10;
static const int N_TOTAL_Laser  = 5;
static const int N_TOTAL_Cosmic = 20;
static const int N_TOTAL_CSI    = 2716;
static const int N_TOTAL_CRATE  = 20;

class E14ReadSumFile{
 public:

  TChain* ch;
  TTree*  trOut;
  TFile*  tfOut;

  int EventNumber;

  int CsiNumber;
  int CsiModID[N_TOTAL_CSI];//CsiNumber
  double CsiTotalE;
  double CsiADC[N_TOTAL_CSI];//CsiNumber
  double CsiEne[N_TOTAL_CSI];//CsiNumber
  double CsiTime[N_TOTAL_CSI];//CsiNumber

  int CC03Number;
  int CC03HitNumber;
  int CC03ModID[N_TOTAL_CC03];//CC03Number
  double CC03TotalE;
  double CC03ADC[N_TOTAL_CC03];//CC03Number
  double CC03Ene[N_TOTAL_CC03];//CC03Number
  double CC03Time[N_TOTAL_CC03];//CC03Number
  
  int OEVNumber;
  int OEVHitNumber;
  int OEVModID[N_TOTAL_OEV];//OEVNumber
  double OEVTotalE;
  double OEVADC[N_TOTAL_OEV];//OEVNumber
  double OEVEne[N_TOTAL_OEV];//OEVNumber
  double OEVTime[N_TOTAL_OEV];//OEVNumber

  int CVNumber;
  int CVHitNumber;
  int CVModID[N_TOTAL_CV];//CVNumber
  double CVTotalE;
  double CVADC[N_TOTAL_CV];//CVNumber
  double CVEne[N_TOTAL_CV];//CVNumber
  double CVTime[N_TOTAL_CV];//CVNumber

  int CrateNumber;
  int CrateHitNumber;
  int CrateModID[N_TOTAL_CRATE];//CrateNumber
  double CrateTotalE;
  double CrateADC[N_TOTAL_CRATE];//CrateNumber
  double CrateEne[N_TOTAL_CRATE];//CrateNumber
  double CrateTime[N_TOTAL_CRATE];//CrateNumber

  int CosmicNumber;
  int CosmicModID[N_TOTAL_Cosmic];//CosmicNumber
  double CosmicTotalE;  
  double CosmicEne[N_TOTAL_Cosmic];//CosmicNumber
  double CosmicTime[N_TOTAL_Cosmic];//CosmicNumber

  int LaserNumber;
  int LaserModID[N_TOTAL_Laser];//LaserNumber
  double LaserTotalE;
  double LaserEne[N_TOTAL_Laser];//LaserNumber
  double LaserTime[N_TOTAL_Laser];//LaserNumber

  short LaserBit;
  short CosmicBitUp;
  short CosmicBitDn;
  short CVBit;
  short CC03Bit;

  double LaserThreshold[N_TOTAL_Laser];
  double CosmicThreshold[N_TOTAL_Cosmic];
  double CVThreshold[N_TOTAL_CV];
  double CC03Threshold[N_TOTAL_CC03];
  
  E14ReadSumFile();
  ~E14ReadSumFile();
  
  bool SetBranchAddress();
  int  GetEntry(long eventNumber);
  long GetEntries();
  int  Add(const char* filename);
  bool SetOutputFile(const  char* filename);
  bool BranchtoTree(TTree* tr);
  bool SumData();
  bool SetTriggerBit();
  bool SumCrate();
  bool Write();
  bool Fill();
};

#endif

