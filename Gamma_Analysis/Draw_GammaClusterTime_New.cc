#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TApplication.h"
#include "ClusterTimeReader.h"
#include "TMath.h"
#include "TF1.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "CsIPoly.h"

double AdjFunc( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double value = p0 + 1.2*p1*exp(p2*x0) + 0*p3*exp(p4*x0);
  return value;
}
double AdjFunc1( double* x, double* par ){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double value = p0 + p1*exp(p2*x0);
  return value;
}
double AdjFuncTest( double* x, double* par){
  double x0 = x[0];
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double value = p0 + p1*exp(p2*x0) + p3*exp(p4*x0);
  return value;

}
double AdjfuncTH( double*x ,double*par ){  
  double value = par[0]*(-1*TMath::Log(1 + par[1])+TMath::Log( 1 + par[1]*exp( par[2]*x[0] )));
  return value;
}

int main( int argc, char** argv) {

  std::string ROOTFILE_GAMMACLUS = std::getenv("ROOTFILE_GAMMACLUS");
  std::string ROOTFILE_WAV       = std::getenv("ROOTFILE_WAV");
  std::string ANALYSISLIB        = std::getenv("ANALYSISLIB");
  
  TF1* TimeAdjFuncEnergy = new TF1("TimeAdjFuncEnergy",AdjFunc,0,2000,5);
  //TimeAdjFuncEnergy->SetParameters(-7.51860e-01,9.57348e-01,-9.55972e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-9.69218e-01,1.12202e+00,-7.52470e-02,0,0);
  //TimeAdjFuncEnergy->SetParameters(-6.32987e-01,7.64508e-01,-9.74328e-02,0,0);//E:200-300;
  //TimeAdjFuncEnergy->SetParameters(-6.38508e-01,7.55330e-01,-9.05449e-02,0,0);//E:200-300,T:10-14
  //TimeAdjFuncEnergy->SetParameters(0,0,0,0,0);

  TF1* TimeAdjFuncHeight = new TF1("TimeAdjFuncHeight",AdjfuncTH,0,16000,3);
  TimeAdjFuncHeight->SetParameters(3.308,0.005644,0.0004385);


  TChain* ch = new TChain("T");
  std::ifstream ifs(Form("%s/Data/RunList/KLRunList_2.txt",ANALYSISLIB.c_str()));
  int tmpRunNumber;
  while( ifs >> tmpRunNumber ){
    std::cout<<  tmpRunNumber  << std::endl;
    ch->Add(Form("%s/run_wav_%d_Cal_FNL_COS.root", ROOTFILE_WAV.c_str(),tmpRunNumber));
  }

  GammaFinder gFinder;


  Int_t EventNumber;
  Int_t RunNumber;
  Int_t nCrystal;
  Int_t ClusterID;
  
  E14GNAnaDataContainer data;
  std::cout<< ch->GetEntries() << std::endl;
  for( Int_t eventIndex = 0; eventIndex  < ch->GetEntries(); eventIndex++){
    ch->GetEntry( eventIndex );
    EventNumber = eventIndex;

    std::list<Cluster> clist;
    std::list<Gammma>  glist;
    std::vector<Klong> klVec;
    
    data.getData( clist );
    data.getData( klVec );
    gFinder.findGamma( clist, glist );
    if( clist.size() != 6 ){ continue; }
    if( glist.size() != 6 ){ continue; }

    std::list<Cluster>::iterator itcl;

  }

  


  /*
  double ene = gam.e();  
  double sinTheta = sin(gam.p3().theta());
  double newx = gam.coex()-L*sinTheta*cos(gam.p3().phi());
  double newy = gam.coey()-L*sinTheta*sin(gam.p3().phi());
  */

  static double const Pcor[2]= {6.49003,0.99254};
  static double const CsIX0  = 18.5;//mm

  for( int ievent = 0; ievent < nEntries; ievent++){
    reader->GetEntry(ievent);
    hisCut->Fill(0);
    Int_t    MEGIndex    = -1;
    Double_t MEGEnergy   = 0;
    Double_t MEGTime     = 0;
    Double_t MEGTimeNon  = 0;
    Double_t MEGTimeNew  = 0;
    Double_t MEGR        = 1000;
    Double_t GTime[6]    = {0.};
    Double_t GRad[6]     = {0.};    
    Double_t GTimeNon[6] = {0.};
    Double_t GTimeNew[6] = {0.};
    Double_t GPhi[6]    = {0.};
    Double_t GE[6]      = {0.};
    Double_t GCE[6]     = {0.};
    Double_t GCT[6]     = {0.};
    Double_t GL[6]      = {0.};
    if( reader->nCluster !=6 ){ continue; }
    bool ClusterRegion = true;
    bool ClusterMinRBool = false;//  < 300 mm at least one crystal 

    for( int clusterIndex = 0; clusterIndex < reader->nCluster; clusterIndex++){
      Double_t clusterR      = reader->ClusterR[clusterIndex];
      Double_t clusterTheta  = reader->ClusterTheta[clusterIndex];
      Double_t clusterPhi    = reader->ClusterPhi[clusterIndex];
      Double_t clusterEnergy = reader->ClusterEnergy[clusterIndex];
      Double_t clusterT      = reader->ClusterT[clusterIndex];
      Double_t clusterZ      = clusterR/TMath::Tan(clusterTheta);
      Double_t clusterL      = clusterR/TMath::Sin(clusterTheta);
      Double_t MECHeight     = 0;
      Double_t MECEnergy     = 0;      
      
      Double_t clusterx      = reader->ClusterPos[clusterIndex][0];
      Double_t clustery      = reader->ClusterPos[clusterIndex][1];
      Double_t clusterz      = reader->ClusterPos[clusterIndex][2];
      Double_t klongz        = reader->KlongPos[2];

      Double_t Length        = sqrt( (klongz-clusterz)*(klongz-clusterz)+ clusterx*clusterx+clustery*clustery);
      GRad[clusterIndex] = clusterR;
      GL[clusterIndex]   = Length;
      GTime[clusterIndex]   = clusterT;
      GE[clusterIndex]   = clusterEnergy;
    }    
    for( int i = 0; i< 6; i++){
      for( int j = i+1; j< 6; j++){
	hisRT->Fill(GRad[i],GTime[i]-GL[i]/299.8-(GTime[j]-GL[i]/299.8));
      }
    }
  }
  hisWT->Write();
  hisHT->Write();
  profTime->Write();
  hisCut->Write();
  hisRTNon->Write();
  hisRTNew->Write();
  hisRT->Write();
  tfout->Close();
}
