#ifndef SemiOfflineHist_h
#define SemiOfflineHist_h

#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"

#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TBox.h"

#include "E14Mapper.h"

class SemiOfflineHist
{
 public:
  
  SemiOfflineHist();
  virtual ~SemiOfflineHist();

  virtual void Fill();
  virtual void SetCrateID( int num ){ CrateID = num; } 
  virtual void SetRunNumber( int num ){ RunNum = num; } 
  virtual void InitHist();
  virtual void WriteAll();

  virtual void SpikeCheck( );

 public:
  //storage
  short nFADC;
  short nCH;
  short nSamples;
  int   EventNo;
  float Pedestal[20][16];
  float IntegratedADC[20][16];
  short PeakHeight[20][16];
  float PeakTime[20][16];
  short Data[20][16][48];

 private:

  int CrateID;
  int RunNum;

  int NDet;
  std::string DetName[6];
  E14Mapper     *map[6];

  TBox *box;

  TFile *hfile;
  TH1F *hLaserADC[20][16];  
  TH1F *hLaserTime[20][16];  
  TH1F *hPed[20][16];  
  TH2F *hRaw[20][2];  
  TCanvas *cLaserRaw;
  TCanvas *cPedestal;
  char psFileName1[256];
  char psFileName2[256];
  TGraphErrors *gLaserADC;  
  TGraphErrors *gLaserTime[20];  
  TH1F *hSpikeHist;  

  int SpikeEvent;
  bool SPIKEID[20][16];

  ClassDef(SemiOfflineHist,1)

    };

#endif // SemiOfflineHist_h
