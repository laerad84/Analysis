#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include <fstream>
#include <iostream>
#include "TProfile.h"

double fcn( double *x ,double *par ){
  double x0     = x[0];
  double L      = 2622;
  double l      = 100;
  //double deltaT = -1*(l*TMath::Sqrt(L*L+x0*x0) - 3*l*L)/300./TMath::Sqrt(L*L+x0*x0);
  //double deltaT0=( ( L * L) + l*TMath::Sqrt(L*L) - 3*l*L)/299.7/TMath::Sqrt(L*L);
  //double deltaT = deltaT0-1*( ( L * L + x0 * x0) + l*TMath::Sqrt(L*L+x0*x0) - 3*l*L)/299.7/TMath::Sqrt(L*L+x0*x0);
  double deltaT0 = L/299.7;
  double deltaT = deltaT0-1*( ( L * L + x0 * x0) )/299.7/TMath::Sqrt(L*L+x0*x0);
  return deltaT;
}

void DrawDelta(){
  std::string ANALISYSLIB = std::getenv("ANALYSISLIB"); 
  gSystem->Load(Form("%s/lib/libAnalysisLib.so",ANALISYSLIB.c_str()));
  IDHandler* handler = new IDHandler();
  CsIImage*  image   = new CsIImage(handler);
  CsIImage*  image0   = new CsIImage(handler);
  std::ifstream ifs("CosmicOut_V4_TimeDeltaResolution.dat");
  int tmpID;
  double tmpDelta;
  double tmpResolution;

  double Delta[2716]={0};
  double Resolution[2716]={0xFFFF};

  while( ifs >> tmpID >> tmpDelta >> tmpResolution ){
    Delta[ tmpID ] = tmpDelta;
    Resolution[ tmpID ] = tmpResolution;
  }
  std::ifstream ifs1("Pi0Peak.dat");
  double Delta1[2716];
  double Resolution1[2716]={-1};
  while( ifs1 >> tmpID >> tmpDelta >> tmpResolution ){
    Delta1[ tmpID ] = tmpDelta;
    Resolution1[ tmpID ] = tmpResolution; 
  }
  
  double x,y;
  
  TGraph* gr        = new TGraph();
  TH2D*   hisDelta  = new TH2D("hisDelta" ,"hisDelta" , 20,0,1000,200,-10,10);
  TH2D*   hisDeltaL = new TH2D("hisDeltaL","hisDeltaL", 20,0,1000,200,-10,10);
  TH1D*   hisDelta1 = new TH1D("hisDelta1","hisDelta1",160,-20,20);
  TH1D*   hisDelta2 = new TH1D("hisDelta2","hisDelta2",160,-20,20);
  TH1D*   hisDelta3 = new TH1D("hisDelta3","hisDelta3",160,-20,20);

  for( int i = 0; i< 2716; i++){
    if( Resolution[i] > 0xFFFE || Resolution1[i] < 0){ continue; }
    handler->GetMetricPosition(i, x, y);
    double Radius = TMath::Sqrt( x*x +y*y  );
    Double_t L = (TMath::Sqrt( 2624*2624 + x*x + y*y ) - 2624)/299.7;
    if(( i >=2240 && i < 2240+120 )|| ( i>=2716-120 )){ continue; }
    image->Fill( i, Delta[i]-L);
    //image->Fill( i, Delta[i]);

    if( i < 2240){
      image0->Fill(i, Delta1[i]);
      gr->SetPoint( gr->GetN(), Radius, Delta[i] );
      hisDelta->Fill( Radius, Delta[i]);
      hisDelta1->Fill(Delta[i]);
      hisDelta2->Fill(Delta1[i]);
      hisDelta3->Fill(Delta[i]-Delta1[i]);
    }else{
      image0->Fill(i,Delta1[i]-16);
      hisDeltaL->Fill( Radius, Delta[i]);
    }
  }
  TF1* fcnt   = new TF1("fcn",fcn, 0,1000,0);

  TCanvas* can = new TCanvas("can","",1000 ,1000);
  can->Divide(2,2);
  can->cd(1);
  hisDelta->Draw("colz");
  TProfile* prof = hisDelta->ProfileX();
  prof->Draw("same");
  fcnt->Draw("same");
  can->cd(2);
  image->DrawWithRange("colz",-2,2);
  can->cd(3);
  /*
  hisDelta1->Draw();
  hisDelta2->SetLineColor(2);
  hisDelta2->Draw("same");
  */
  hisDeltaL->Draw("colz");
  TProfile* profL = hisDeltaL->ProfileX();
  profL->Draw("same");
  
  fcnt->Draw("same");
  can->cd(4);
  image0->DrawWithRange("colz",-5,5);
  //hisDelta3->Draw();
}

    
