#include <ifstream>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TH2D.h"
#include "IDHandler.h"
#include "CsIImage.h"

int main( int argc, char** argv){
  
  IDHandler* handler = new IDHandler("Data/crystal.txt");
  CsIImage*  image   = new CsIImage(handler);
  
  std::vector<int> runNumberList;
  std::vector<int>::iterator itRun;
  
  std::ifstream ifs(argv[1]);
  while(1){
    int runNumber;
    ifs >> runNumber;
    if( ifs.eof){break;}
    runNumberList.push_back(runNumber);
  }
  
  TGraphErrors* grChange[2716];
  for( int i = 0; i< 2716; i++){
    grChange[i] = new TGraphErrors();
    grChange[i]->SetNameTitle(Form("GainChange%04d",i),
			      Form("GainChange%04d",i));
  }
  TTree* tr = new TTree("LaserOut","");
  Double_t laserOut[2716];
  Double_t PINout[4];
  
  for( int itRun = runNumber.begin(); itRun != runNumber.end(); itRun++){
    TFile* tfCam = new TFile(Form("laserRootFiles/laserdata_%04d.root",*itRun));
    TFile* tfVME = new TFile(Form("laserRootFiles/laser_%04d.root",*itRun));
    
    TTree* trCam = (TTree*)tfCam->Get("laserEventData");
    TH1D* hisCHLaser[2716];
    for( int i = 0; i< 2716; i++){
      hisCHLaser[i] = (TH1D*)tfVME->Get(Form("laser%04d",i));
    }
    TH1D* hisOut[4];
    TH1D* hisPed[4];
    Double_t PINOut[4];
    Double_t PINErr[4];
    for( int ipin = 0; ipin < 4; ipin++){
      hisOut[ipin] = new TH1D(Form("hisOut%04d",i),
			      Form("hisOut%04d",i),
			      500,0,4000);
      hisPed[ipin] = new TH1D(Form("hisOut%04d",i),
			      Form("hisOut%04d",i),
			      500,0,4000);
      trCam->Project(hisOut[ipin],
		     Form("PINData[%d]",ipin),
		     Form("PINData[%d]>1500",ipin));
      trCam->Project(hisPed[ipin],
		     Form("PINData[%d]",ipin),
		     Form("PINData[%d]<1500",ipin));
      PINOut[ipin] = hisOut[ipin]->GetMean() - hisPed[ipin]->GetMean();
      PINErr[ipin] = hisOut[ipin]->GetRMS();
    }
  }  
}  

