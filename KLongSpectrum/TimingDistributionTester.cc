#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"

#include "TCanvas.h"
#include <iostream>
#include "TMath.h"
#include "TVector2.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "gamma/Gamma.h"
#include "cluster/Cluster.h"
#include "klong/Klong.h"
#include "pi0/Pi0.h"
#include <vector>
#include <list>
#include "TH2.h"
#include <string>
#include <cstdlib>
#include <cstdio>
#include "T0Manager.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TSpline.h"
#include "TGraphErrors.h"
const double KLMass = 497.648;//MeV

double KLSpectrum(double* x,double*p){
  //double klp   = TMath::Sqrt( x[0]*x[0]-KLMass*KLMass);  
  double klp = x[0];
  double sigma = p[2]*(1-(p[3]+p[4]*klp)*(klp-p[1]));
  double value = p[0]*TMath::Exp(-0.5*(klp-p[1])*(klp-p[1])/(sigma*sigma));
  return value;
}

double const sol = 299.792458;//[mm/nsec]
double const solc= 80;//[mm/nsec]
double const Pcor[2]={6.49003,0.99254};
double const CsIX0=18.5;//mm

double showerDepth(double x){
  double L = CsIX0*(Pcor[0]+Pcor[1]*log(x/1000.));//mm  
  return L;
}
double showerTimeDelay(Klong kl, Gamma g){
  double depth     = showerDepth(g.e());
  double cosTheta  = TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  double delayTime = 500/solc - depth*(1-cosTheta)/sol;
  return delayTime;
}
double showerTimeDelayAdj(Klong kl, Gamma g){
  double clusterTimingConstant = 0.541812;
  double depth     = showerDepth( g.e() );
  double cosTheta  = TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  //double delayTime = 500/solc - depth*(1-cosTheta)/sol + 0.541812*TMath::Log10(TMath::E())*log(g.e()/1000);
  double delayTime = depth*(1/sol-cosTheta/solc);
  return delayTime;
}
double ThetaConstant(Klong kl, Gamma g){
  double cosTheta  =1-sol/solc*TMath::Abs(TMath::Sqrt(TMath::Power(g.x() - kl.vx(),2)+TMath::Power(g.y() - kl.vy(),2))/(g.z() - kl.vz()));
  return cosTheta;
}

double gammaLOF( Klong kl, Gamma g){
  double length = 0;
  length = sqrt( pow(g.x()-kl.vx(),2)+ pow(g.y()-kl.vy(),2)+pow(g.z()-kl.vz(),2));
  return length;
}

double HeightDelay( double height ){
  double value = TMath::Log( 1 + 0.03566*TMath::Exp( height/1621));
  return value;
}

/*
TFile* tflin = new TFile("TimingLinearityFuncLaser.root");
TGraphErrors* grTimingLinearity = (TGraphErrors*)tflin->Get("grTimingLinearityLaser");
TSpline3* spl = new TSpline3("spl",grTimingLinearity);
double HeightAdjFunc(double* x, double* par){
  return spl->Eval(x[0]);
}
*/
//void DistributionTester(){
int main( int argc, char** argv){
  int nTimeIteration =  atoi( argv[1] );
  TApplication* app = new TApplication("app",&argc,argv);
  TF1* soltFunc  = new TF1("soltfunc",KLSpectrum,0,12,5);
  TF1* sugarFunc = new TF1("sugarfunc",KLSpectrum,0,12,5);
  const double soltPar[5] = { 1, 1.37,7.48e-1,-2.9e-1,1.68e-2};
  const double sugarPar[5] = {1,1.41991,0.810237,-0.301413,0.017092};
  soltFunc->SetParameters(soltPar);
  sugarFunc->SetParameters(sugarPar);

  const int nFile = 1;
  TFile* tf;
  TTree* tr;
  char* name = "DATA_NONTIMECAL";//"SIM","3pi0_OldComp","WAVNOCV","3pi0_OldComp_wopi0","3pi0_noCompNoCal","3pi0_LaserComp_NOCV"

  tf = new TFile(Form("Kl_Total_%s.root",name));
  tr = (TTree*)tf->Get(Form("trKL"));
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];
  int      CsiNumber;
  double   CsiSignal[2716];
  int      CsiModID[2716];
  
  TFile* tfOut = new TFile(Form("DistributionTimeData_ShowerHeight_%s_%d.root",name,nTimeIteration),"recreate");
  TH2D* hisTimeID;
  TH2D* hisAdjTimeID;
  hisTimeID    = new TH2D(Form("hisTimeID_%d",nTimeIteration),Form("hisTimeID_%d",nTimeIteration),
			  2716,0, 2716,500,-40,60);
  hisAdjTimeID = new TH2D(Form("hisAdjTimeID_%d",nTimeIteration),Form("hisAdjTimeID_%d",nTimeIteration),
			  2716,0, 2716,500,-40,60);
  //double sol = 299.792458;//[mm/nsec]
  E14GNAnaDataContainer data;
  data.setBranchAddress( tr );
  tr->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
  tr->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
  tr->SetBranchAddress("CsiNumber",&CsiNumber);
  tr->SetBranchAddress("CsiSignal",CsiSignal);//CsiNumber
  tr->SetBranchAddress("CsiModID",CsiModID);//CsiNumber
  /*T0Manager* man = new T0Manager();
  if( nTimeIteration > 0 ){
    if( !(man->ReadFile(Form("TimeOffset_%d.dat",nTimeIteration)))){
      std::cout<< "No file exist" << std::endl;
      return -1;
    }
  }
  */
  //man->PrintOffset();

  Double_t TimeOffset[2716]={0};
  if( nTimeIteration > 1 ){ 
    std::ifstream ifs(Form("TimeOffset_ShowerHeight_%d.dat",nTimeIteration));
    if( !ifs.is_open() ){
      std::cout<< "No CalibrationFile" << std::endl;
      return -1;
    }
    int tmpID;
    double tmpOffset;
    while( ifs >> tmpID >> tmpOffset ){
      TimeOffset[tmpID] = tmpOffset;
      std::cout<<tmpID << "\t" <<  TimeOffset[tmpID] << std::endl;
    }
    ifs.close();
  }
  for( int ievent = 0; ievent < tr->GetEntries(); ievent++){      
    tr->GetEntry(ievent);
    //if( ievent  >= 100000 ){ break ; } 
    if(CsiNumber > 500 ){ continue; }
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    data.getData(clist);
    data.getData(glist);
    data.getData(klVec);
    
    
    double klpos[3]={0};
    klpos[0] = klVec[0].vx();
    klpos[1] = klVec[1].vy();
    klpos[2] = klVec[2].vz();
    if( klVec[0].chisqZ() > 10 ){ continue; }
    std::list<Gamma>::iterator git = glist.begin();
    int    g0crystalID = (*git).clusterIdVec()[0];
    double g0time   = (*git).clusterTimeVec()[0];
    double g0length = TMath::Sqrt( TMath::Power(((*git).x() - klVec[0].vx()),2) 
				 + TMath::Power(((*git).y() - klVec[0].vy()),2) 
				 + TMath::Power(((*git).z() - klVec[0].vz()),2));
    double g0Offset = TimeOffset[g0crystalID];//man->GetT0Offset(g0crystalID);
    double g0Shower = showerTimeDelayAdj(klVec[0],(*git));
    double g0Delta  = g0Offset+g0length/sol+g0Shower;//man->GetT0Offset(g0crystalID);

    if( g0crystalID >= 2240 ){continue; }
    
    double g0Ene = (*git).clusterEVec()[0];
    if(g0Ene > 400 ){ continue; }
    if(g0Ene < 100 ){ continue; }

    double GammaTimeSigma=0;
    double GammaTimeMean = 0; 
    int    CrystalID[6]={-1,-1,-1,-1,-1,-1};
    double CrystalHeight[6]={0,0,0,0,0,0};
    double HeightOffset[6]={0,0,0,0,0,0};

    git = glist.begin();
    for( int igamma = 0;igamma < 6; igamma++,git++){
      GammaTimeSigma += (*git).t()*(*git).t();
      GammaTimeMean  += (*git).t();
      CrystalID[igamma] =  (*git).clusterIdVec()[0];
      for( int icsi = 0; icsi < CsiNumber; icsi++){
	if( CrystalID[igamma] == CsiModID[icsi]){
	  CrystalHeight[igamma] = CsiSignal[icsi];
	  HeightOffset[igamma] = HeightDelay(CrystalHeight[igamma]);
	  break;
	}
      }
    }
    g0Delta+=HeightOffset[0];

    GammaTimeMean /= 6;
    GammaTimeSigma = sqrt((GammaTimeSigma/6) - (GammaTimeMean*GammaTimeMean));
    if( nTimeIteration > 10 ){ 
      if( GammaTimeSigma > 10 ){ continue; }
    }

    git = glist.begin();
    git++;
    if(CrystalHeight[0] > 4000 ){ continue; }
    if(CrystalHeight[0] < 200 ){continue; }
    for( int igamma = 1; igamma < 6; igamma++,git++){
      if( CrystalHeight[igamma] > 4000 ){ continue; }
      if( CrystalHeight[igamma] < 200  ){ continue; }
      int crystalID = (*git).clusterIdVec()[0];
      double Ene = (*git).clusterEVec()[0];
      if(Ene > g0Ene){ continue; }
      if(Ene < 100 ){ continue; }
      double Offset = TimeOffset[crystalID];//man->GetT0Offset(crystalID);
      double length = TMath::Sqrt( TMath::Power(((*git).x() - klVec[0].vx()),2) 
				   + TMath::Power(((*git).y() - klVec[0].vy()),2) 
				   + TMath::Power(((*git).z() - klVec[0].vz()),2));
      double Shower = showerTimeDelayAdj(klVec[0],(*git));
      double Delta = Offset+length/sol+Shower+HeightOffset[igamma];//man->GetT0Offset(crystalID);
      hisTimeID->Fill(crystalID,((*git).clusterTimeVec()[0]-Offset)-(g0time-g0Offset));
      hisAdjTimeID->Fill(crystalID,((*git).clusterTimeVec()[0]-Delta)-(g0time-g0Delta));
    }
  }
  /// Evaluation T0 /// 
  double gausFit[2] ={0};
  if( nTimeIteration < 5 ){ 
    gausFit[0] = -10*(5-nTimeIteration);
    gausFit[1] = 10*(5-nTimeIteration);
  }else{
    gausFit[0] = -5;
    gausFit[1] = 5;
  }
  TF1* func = new TF1("func","gaus",gausFit[0],gausFit[1]);

  hisTimeID->FitSlicesY(func);
  hisAdjTimeID->FitSlicesY(func);
  TProfile* profTimeID = hisTimeID->ProfileX();
  TProfile* profAdjTimeID = hisAdjTimeID->ProfileX();

  TH1D *hisTimeID_0 = (TH1D*)gDirectory->Get(Form("%s_0",hisTimeID->GetName()));
  TH1D *hisTimeID_1 = (TH1D*)gDirectory->Get(Form("%s_1",hisTimeID->GetName()));
  TH1D *hisTimeID_2 = (TH1D*)gDirectory->Get(Form("%s_2",hisTimeID->GetName()));
  TH1D *hisAdjTimeID_0 = (TH1D*)gDirectory->Get(Form("%s_0",hisAdjTimeID->GetName()));
  TH1D *hisAdjTimeID_1 = (TH1D*)gDirectory->Get(Form("%s_1",hisAdjTimeID->GetName()));
  TH1D *hisAdjTimeID_2 = (TH1D*)gDirectory->Get(Form("%s_2",hisAdjTimeID->GetName()));

  
  std::cout << "Slice" << std::endl;
  int calibrated[2716] = {-1};
  int nCalibrated = 0;
  double mean=0;
  std::ofstream ofs(Form("TimeOffset_ShowerHeight_%d.dat",nTimeIteration+1));
  double timeCalFactor[2716] = {0};
  double tmptimeCalFactor[2716] = {0};
  TH1D* tmpTimeOffsetDistribution = new TH1D("tmpTimeOffset","tmpTimeOffset",200,-20,20);
  /*
  for( int i = 1; i<= hisAdjTimeID_1->GetNbinsX(); i++){
    tmptimeCalFactor[i-1] = hisAdjTimeID_1->GetBinContent(i);
    if( profAdjTimeID->GetBinEntries(i) < 1 ){continue;}
    tmpTimeOffsetDistribution->Fill( tmptimeCalFactor[i-1] - mean);
    std::cout<< i-1 << "\t" << tmptimeCalFactor[i-1] << std::endl;
  }
*/

  TH1D* TimeOffsetDistribution = new TH1D(Form("TimeOffsetDistribution_%d",nTimeIteration),
					  Form("TimeOffsetDistribution_%d",nTimeIteration),
					  200,-20,20);
  TGraphErrors* grTimeOffset = new TGraphErrors();
  grTimeOffset->SetNameTitle("grTimeOffset","grTimeOffset");
  TCanvas* can  =new TCanvas("can","can",800,800);
  if( nTimeIteration <= 10 ){
    for( int i = 0; i< 2716; i++){
      if( profAdjTimeID->GetBinEntries(i+1) < 1 ){ continue;}
      calibrated[i] = 1; 
      tmptimeCalFactor[i] = profAdjTimeID->GetBinContent(i+1);
      tmpTimeOffsetDistribution->Fill( tmptimeCalFactor[i] );
    }
    
    mean = tmpTimeOffsetDistribution->GetMean();
    
    for( int i = 0; i< 2716; i++){
      if( calibrated[i] >0){      
	timeCalFactor[i] = profAdjTimeID->GetBinContent(i+1) - mean;
	TimeOffsetDistribution->Fill( profAdjTimeID->GetBinContent(i+1) - mean);
      }else{
	timeCalFactor[i] = 0;
      }
      ofs << i  << "\t" << timeCalFactor[i]+TimeOffset[i] << "\n";
    }
    ofs.close();
  }else{
    for( int i = 0; i< 2716; i++){
      TH1D* hist = hisAdjTimeID->ProjectionY(Form("hist_%d",i),i+1,i+1);
      if( hist->GetEntries() < 1 ){
	grTimeOffset->SetPoint(grTimeOffset->GetN(),i,0);
	grTimeOffset->SetPointError(grTimeOffset->GetN()-1,0,0);
	continue;
      }
      /*
      hist->Draw();
      can->Update();
      can->Modified();
      getchar();
      */

      double FitWidth = 20-(nTimeIteration-10)*2;
      if(FitWidth < 2 ){ FitWidth = 2;}
      hist->Fit("gaus","","",-hist->GetRMS()*5,hist->GetRMS()*5);
      TF1* func = hist->GetFunction("gaus");
      calibrated[i]    = 1;
      timeCalFactor[i] = func->GetParameter(1);
      double error     = func->GetParError(1);
      tmpTimeOffsetDistribution->Fill( tmptimeCalFactor[i] );
      grTimeOffset->SetPoint(grTimeOffset->GetN(),i,timeCalFactor[i] );
      grTimeOffset->SetPointError( grTimeOffset->GetN()-1,0,error );      
      std::cout<< i << "\t" << tmptimeCalFactor[i] << std::endl;
    }
    mean = tmpTimeOffsetDistribution->GetMean();
    for( int i = 0 ; i< 2716; i++){
      if( calibrated[i] > 0 ){
	timeCalFactor[i] = tmptimeCalFactor[i] - mean;
      }else{
	timeCalFactor[i] = 0;
      }

      ofs << i << "\t" << timeCalFactor[i]+TimeOffset[i] << "\n";
    }
  }
  
  
  /*
  for( int ibin = 1; ibin < hisAdjTimeID_1->GetNbinsX();ibin++){
    if( profAdjTimeID->GetBinEntries(ibin) > 0 ){
      timeCalFactor[ibin-1] = hisAdjTimeID_1->GetBinContent(ibin)-1;
      TimeOffsetDistribution->Fill( timeCalFactor[ibin-1] - mean);
    }else{
      timeCalFactor[ibin-1] = 0;
    }
    ofs << ibin-1  << "\t" << timeCalFactor[ibin-1]+TimeOffset[ibin-1] << "\n";
  }
  */
  ofs.close();
  

  hisTimeID->Write();
  hisAdjTimeID->Write();
  TimeOffsetDistribution->Write();
  /*
  hisTimeID_0->Write();
  hisTimeID_1->Write();
  hisTimeID_2->Write();
  hisAdjTimeID_0->Write();
  hisAdjTimeID_1->Write();
  hisAdjTimeID_2->Write();
  */
  
  grTimeOffset->Write();
  tfOut->Close();
  //app->Run();
}
