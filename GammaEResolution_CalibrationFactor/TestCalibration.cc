#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"

#include "TCanvas.h"
#include "TMath.h"

#include "cluster/ClusterFinder.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/DigiReader.h"
#include "gamma/GammaFinder.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <list>
int main( int argc, char** argv){

  Int_t InjectionEnergy    = atoi(argv[1]);
  Int_t InjectionAngle     = atoi(argv[2]);
  Int_t CalPrecision = atoi( argv[3]);

  //Int_t RunNumber = atoi(argv[3]);

  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string iFileForm = Form("%s/Cluster_%dMeV_%ddeg-1E5-1.root",
			       ROOTFILE_GAMMACLUS.c_str(),InjectionEnergy,InjectionAngle);
  
  TFile* tf = new TFile(iFileForm.c_str());
  TTree* tr = (TTree*)tf->Get("T");
  Int_t    nDigi;
  Int_t    CsiID[2716];
  Double_t CsiEne[2716];
  Double_t CsiTime[2716];
  tr->SetBranchAddress("nDigi",&nDigi);
  tr->SetBranchAddress("CsiID",CsiID);
  tr->SetBranchAddress("CsiEne",CsiEne);
  tr->SetBranchAddress("CsiTime",CsiTime);



  std::string ofname = Form("Output_%dMeV_%ddeg_%d.root",InjectionEnergy,InjectionAngle,CalPrecision);
  TFile* tfout = new TFile(ofname.c_str(),"RECREATE");
  TTree* trout = new TTree("Tree","Tree");

  Double_t CalFactor[2716];  
  std::ifstream ifs(Form("CalFactor_%d.txt",CalPrecision));
  Int_t tmpID;
  Double_t tmpCal; 
  TH1D* hisCalFactor = new TH1D("hisCalFactor","hisCalFactor",100,0.5,1.5);
  TH2D* hisGammaE    = new TH2D("hisGammaE","hisGammaE",100,0,1000,100,0,1000);
  TH1D* hisGammaERatio = new TH1D("hisGammaERatio","hisGammaERatio",400,0,2);
  TH1D* hisGammaEResolution = new TH1D(Form("hisGammaEResolution_%d_%d_%d",InjectionEnergy,InjectionAngle,CalPrecision),"hisGammaEResolution",400,0,2);

  while( ifs >> tmpID >> tmpCal ){
    CalFactor[tmpID] = tmpCal;
  }


  E14GNAnaDataContainer data;
  data.branchOfClusterList( trout );
  data.branchOfGammaList( trout );
  data.branchOfDigi( trout);

  ClusterFinder clusterFinder;
  GammaFinder   gFinder;

  Int_t nEntries = tr->GetEntries();
  for( int ievt = 0; ievt < nEntries; ievt++){
    tr->GetEntry(ievt);
    Int_t CsiDigiID[2716];
    Double_t CsiDigiE[2716];
    Double_t CsiDigiTime[2716];
    Int_t nCsiDigi=0;
    for( int idigi = 0; idigi < nDigi; idigi++){
      if( CsiEne[idigi]*CalFactor[CsiID[idigi]] > 3 ){
	CsiDigiID[nCsiDigi] = CsiID[idigi];
	CsiDigiE[nCsiDigi] = CsiEne[idigi]*CalFactor[CsiID[idigi]];
	CsiDigiTime[nCsiDigi] = CsiTime[idigi]; 
	nCsiDigi++;
      }
    }    
    std::list<Cluster> clist = clusterFinder.findCluster(nCsiDigi,CsiDigiID,CsiDigiE,CsiDigiTime);
    std::list<Gamma> glist;
    if( clist.size() != 1 ){ continue;}
    gFinder.findGamma( clist , glist );
    if( glist.size() != 1 ){ continue;}
    double gEne = glist.front().e();
    double Theta = InjectionAngle*TMath::Pi()/180;
    double x   = glist.front().x();
    double y   = glist.front().y();
    double Phi = TMath::ATan2(x,y);
    double px = gEne*sin(Theta)*cos(Phi);
    double py = gEne*sin(Theta)*sin(Phi);
    double pz = gEne*cos(Theta);
    glist.front().setP3(px,py,pz);
    E14GNAnaFunction::getFunction()->correctEnergyWithAngle(glist.front());
    double gEneFix = glist.front().e();
    hisGammaE->Fill(gEneFix,gEne);
    hisGammaERatio->Fill(gEneFix/gEne);
    hisGammaEResolution->Fill(gEneFix/InjectionEnergy);

    data.setData( clist );
    data.setData( glist );
    trout->Fill();
    data.eventID++;
  }
  hisGammaE->Write();
  hisGammaERatio->Write();
  hisGammaEResolution->Write();
  trout->Write();
  tfout->Close();
  std::cout<< "Finished" << std::endl;
  return 0;

}
