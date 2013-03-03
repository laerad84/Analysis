#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"



int main( int argc, char** argv ){
  Int_t  iRun  = atoi(argv[1]);
  TFile* tfOut = new TFile(Form("LaserTimeStability_%d.root",iRun),"RECREATE");
  
  TH1D* hisLaserHeight = new TH1D(Form("hisLaserHeight_%d",iRun),Form("hisLaserHeight_%d",iRun),
			    160,0,4800);
  TH1D* hisLaserTime = new TH1D(Form("hisLaserTime_%d",iRun),Form("hisLaserTime_%d",iRun),
			  50,0,500);
  TH1D* hisTimeDeltaLaser[2716];
  TH1D* hisTimeDelta[2716];
  TH1D* hisHeight[2716];
  for( int ich = 0; ich < 2716; ich++)
      hisHeight[tmpCsIID]->Fill(CsiSignal[ich]);
    }
  }
  
  tfOut->cd();

  for( int ich = 0; ich < 2716; ich++){
      hisTimeDeltaLaser[ich]->Write();
      hisTimeDelta[ich]->Write();
      hisHeight[ich]->Write();
  }
}


