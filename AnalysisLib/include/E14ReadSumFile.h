#ifndef __E14ReadSumFile__h__
#define __E14ReadSumFile__h__
#include "TChain.h"

static const int N_TOTAL_CC03   = 32;
static const int N_TOTAL_OEV    = 44;
static const int N_TOTAL_CV     = 10;
static const int N_TOTAL_Laser  = 5;
static const int N_TOTAL_Cosmic = 20;
static const int N_TOTAL_CSI    = 2716;
static const int N_TOTAL_Etc    = 16;
static const int N_TOTAL_Sci    = 1; 
static const int nMaxTrack      = 400;

class E14ReadSumFile{
 public:

  TChain* ch;

  int EventNumber;

  int CsiNumber;
  int CsiModID[N_TOTAL_CSI];//CsiNumber
  double CsiEne[N_TOTAL_CSI];//CsiNumber
  double CsiTime[N_TOTAL_CSI];//CsiNumber
  double CsiIntegratedADC[N_TOTAL_CSI];//CsiNumber

  int CC03Number;
  int CC03ModID[N_TOTAL_CC03];//CC03Number
  double CC03Ene[N_TOTAL_CC03];//CC03Number
  double CC03Time[N_TOTAL_CC03];//CC03Number
  
  int OEVNumber;
  int OEVModID[N_TOTAL_OEV];//OEVNumber
  double OEVEne[N_TOTAL_OEV];//OEVNumber
  double OEVTime[N_TOTAL_OEV];//OEVNumber

  int CVNumber;
  int CVModID[N_TOTAL_CV];//CVNumber
  double CVEne[N_TOTAL_CV];//CVNumber
  double CVTime[N_TOTAL_CV];//CVNumber

  int CosmicNumber;
  int CosmicModID[N_TOTAL_Cosmic];//CosmicNumber
  double CosmicEne[N_TOTAL_Cosmic];//CosmicNumber
  double CosmicTime[N_TOTAL_Cosmic];//CosmicNumber

  int LaserNumber;
  int LaserModID[N_TOTAL_Laser];//LaserNumber
  double LaserEne[N_TOTAL_Laser];//LaserNumber
  double LaserTime[N_TOTAL_Laser];//LaserNumber

  int EtcNumber;
  int EtcModID[N_TOTAL_Etc];//EtcNumber
  double EtcEne[N_TOTAL_Etc];//EtcNumber
  double EtcTime[N_TOTAL_Etc];//EtcNumber

  int SciNumber;
  int SciModID[N_TOTAL_Sci];//SciNumber
  double SciEne[N_TOTAL_Sci];//SciNumber
  double SciTime[N_TOTAL_Sci];//SciNumber


  Int_t    nTrack;
  UShort_t track[nMaxTrack];//nTrack
  Short_t  mother[nMaxTrack];//nTrack
  Int_t    pid[nMaxTrack];//nTrack
  Float_t  mass[nMaxTrack];//nTrack
  Float_t  ek[nMaxTrack];//nTrack
  Float_t  end_ek[nMaxTrack];//nTrack
  Double_t p[nMaxTrack][3];//nTrack
  Double_t v[nMaxTrack][3];//nTrack
  Double_t end_p[nMaxTrack][3];//nTrack
  Double_t end_v[nMaxTrack][3];//nTrack
  
  E14ReadSumFile(Int_t Flag = 0);
  ~E14ReadSumFile();
  
  bool SetBranchAddress();
  bool ReadTrack();
  int  GetEntry(long eventNumber);
  long GetEntries();
  int  Add(const char* filename);
};

#endif

