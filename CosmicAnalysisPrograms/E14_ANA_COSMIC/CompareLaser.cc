#include <fstream>
#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TH2D.h"
#include "IDHandler.h"
#include "CsIImage.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TFile.h"

int main( int argc, char** argv){
  
  IDHandler* handler = new IDHandler();
  CsIImage*  image   = new CsIImage(handler);
  
  std::vector<int> runNumberList;
  std::vector<int>::iterator itRun;
  
  std::ifstream ifs(argv[1]);
  int runNumber;
  while(1){
    ifs >> runNumber;
    if( ifs.eof()){break;}
    runNumberList.push_back(runNumber);
    std::cout<< runNumber << std::endl;
  }
  
  TGraphErrors* grChange[2716];
  for( int i = 0; i< 2716; i++){
    grChange[i] = new TGraphErrors();
    grChange[i]->SetNameTitle(Form("GainChange%04d",i),
			      Form("GainChange%04d",i));
  }
  //TTree* tr = new TTree("LaserOut","");
  //Double_t laserOut[2716];
  //Double_t PINout[4];
  
  for( itRun = runNumberList.begin();
       itRun!= runNumberList.end(); 
       ++itRun){
    TFile* tfCam = new TFile(Form("laserROOTFiles/laserdata_%04d.root",*itRun));
    TFile* tfVME = new TFile(Form("laserROOTFiles/laser_%04d.root",*itRun));

    TH1D* hisOut[4];
    TH1D* hisPed[4];
    Double_t PINOut[4];
    Double_t PINErr[4];
    Double_t x,y;
    TTree* trCam = (TTree*)tfCam->Get("laserEventData");
    for( int ipin = 0; ipin < 4; ipin++){
      
      hisOut[ipin] = new TH1D(Form("hisOut%04d",ipin),
			      Form("hisOut%04d",ipin),
			      500,0,4000);
      hisPed[ipin] = new TH1D(Form("hisPed%04d",ipin),
			      Form("hisPed%04d",ipin),
			      500,0,4000);
      trCam->Project(hisOut[ipin]->GetName(),
		     Form("PINData[%d]",ipin),
		     Form("PINData[%d]>1500",ipin));
      trCam->Project(hisPed[ipin]->GetName(),
		     Form("PINData[%d]",ipin),
		     Form("PINData[%d]<1500",ipin));
      PINOut[ipin] = hisOut[ipin]->GetMean() - hisPed[ipin]->GetMean();
      PINErr[ipin] = hisOut[ipin]->GetRMS();            
      hisOut[ipin]->Clear();
      hisPed[ipin]->Clear();
      std::cout << PINOut[ipin] << std::endl;
    }
    std::cout<<"PIN" << std::endl;
    
    TH1D* hisCHLaser[2716];        
    for( int ich = 0; ich< 2716; ich++){      
      hisCHLaser[ich] = NULL;
      hisCHLaser[ich] = (TH1D*)tfVME->Get(Form("laser%d",ich));
      if(hisCHLaser[ich] == NULL ){
	grChange[ich]->SetPoint(grChange[ich]->GetN(),*itRun,0);
      }
      double out = hisCHLaser[ich]->GetMean();

      handler->GetMetricPosition(ich,x,y);      
      int pinid;
      if( x < 0 ){
	if( y>0 ){ 
	  pinid = 0;
	}else{
	  pinid  =1 ;
	}
      }else{
	if( y>0){
	  pinid = 2;
	}else{
	  pinid = 3;
	}
      }      
      grChange[ich]->SetPoint(grChange[ich]->GetN(),*itRun,out/PINOut[pinid]);
      hisCHLaser[ich]->Clear();
    }       
  }
  std::cout<< "OUT" << std::endl;
  TFile* tfOUTPUT = new TFile("test.root","recreate");
  for( int ich = 0; ich < 2716; ich++){
    grChange[ich]->Write();
  }
  tfOUTPUT->Close();
    
    
}  

