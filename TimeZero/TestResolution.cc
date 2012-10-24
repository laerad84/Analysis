#include <iostream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include "EnergyConverter.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <fstream>
#include "TApplication.h"
#include "IDHandler.h"
#include "CsIImage.h"

#include "TProfile.h"

int 
main( int argc ,char** argv ){
  std::string ANALIB = std::getenv("ANALYSISLIB");
  TApplication* app = new TApplication("app",&argc, argv );
  EnergyConverter* EConverter = new EnergyConverter();
  EConverter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",ANALIB.c_str()));

  std::ifstream ifs("TimeResolutionCosmic.dat");
  
  double Resolution[ 2716 ]={0xFFFF}; 
  double Delta[2716]       ={0xFFFF};

  IDHandler* handler = new IDHandler();
  CsIImage* imageRes = new CsIImage( handler );
  CsIImage* imageOut = new CsIImage( handler );
  imageRes->SetTitle("Resolution[ns]");
  imageOut->SetTitle("Calibration Constant[ MeV / cnt ]");
  int tmpID;
  double tmpRes;
  double tmpDelta;
  while ( ifs >> tmpID >> tmpDelta >> tmpRes ){
    Delta[ tmpID ] = tmpDelta; 
    Resolution[ tmpID ] = tmpRes;
  }

  double LY[2716]  = {0};
  double Gain[2716]= {0};

  TFile* tfCsIData = new TFile("CsI_Data.root");
  TTree* trCsIData = (TTree*)tfCsIData->Get("pmt_crystal");

  int    ID;
  int    PMTID;
  int    CrystalID;
  double PMTGain;
  double CrystalLY;
  trCsIData->SetBranchAddress("ID"       ,&ID);
  trCsIData->SetBranchAddress("PMTID"    ,&PMTID);
  trCsIData->SetBranchAddress("CrystalID",&CrystalID);
  trCsIData->SetBranchAddress("PMTGain"  ,&PMTGain);
  trCsIData->SetBranchAddress("CrystalLY",&CrystalLY);
  for( int i = 0; i< trCsIData->GetEntries(); i++){
    trCsIData->GetEntry( i );
    LY[ID] = CrystalLY;
    Gain[ID] = PMTGain;
  }
      
  TGraph* gr            = new TGraph();
  TGraph* grSmall       = new TGraph();
  TGraph* grLarge       = new TGraph();
  TGraph* grLYSmall     = new TGraph();
  TGraph* grLYLarge     = new TGraph();
  TGraph* grGainSmall   = new TGraph();
  TGraph* grGainLarge   = new TGraph();
  TGraph* grLYSmallGood = new TGraph();
  TGraph* grLYSmallBad  = new TGraph();
  TH2D*   hisResYDep    = new TH2D("hisResYDep","hisResYDep",48,-600,600,40,0,2);

  grSmall->SetMarkerColor( 2 );
  grLarge->SetMarkerColor( 3 );
  grSmall->SetMarkerStyle( 6 );
  grLarge->SetMarkerStyle( 6 );

  
  grLYSmall->SetMarkerStyle( 6 );
  grLYLarge->SetMarkerStyle( 6 );
  grLYSmall->SetMarkerColor( 2 );
  grLYLarge->SetMarkerColor( 3 );
  grLYSmallGood->SetMarkerStyle( 6 );
  grLYSmallBad->SetMarkerStyle( 6 );
  grLYSmallGood->SetMarkerColor( 2 );
  grLYSmallBad->SetMarkerColor( 3 );
  

  grGainSmall->SetMarkerStyle( 6 );
  grGainLarge->SetMarkerStyle( 6 );
  grGainSmall->SetMarkerColor( 2 );
  grGainLarge->SetMarkerColor( 3 );

  grLYSmall    ->SetNameTitle("grLYSmall"    ,"LYSmall;Resolution;LY");
  grLYLarge    ->SetNameTitle("grLYLarge"    ,"LYLarge;Resolution;LY");
  grLYSmallGood->SetNameTitle("grLYSmallGood","LYSmallGood;Resolution;LY");
  grLYSmallBad ->SetNameTitle("grLYSmallBad" ,"LYSmallBad;Resolution;LY");
  grGainSmall  ->SetNameTitle("grGainSmall"  ,"GainSmall;Resolution;Gain");
  grGainLarge  ->SetNameTitle("grGainLarge"  ,"GainLarge;Resolution;Gain");

  double x,y;
  for( int i = 0; i< 2716; i++){
    if( Resolution[i] >= 0xFFFF ){ continue; }
    gr->SetPoint(gr->GetN() , Resolution[ i ], EConverter->GetCalibrationConstant( i ));    
    imageRes->Fill( i, Resolution[ i ] ); 
    imageOut->Fill( i, EConverter->GetCalibrationConstant( i ));
    
    if( i < 2240 ){
      grSmall    ->SetPoint( grSmall->GetN()    , Resolution[ i ], EConverter->GetCalibrationConstant( i ));    
      grLYSmall  ->SetPoint( grLYSmall->GetN()  , Resolution[ i ], LY[i] );
      grGainSmall->SetPoint( grGainSmall->GetN(), Resolution[ i ], Gain[i] );
      handler->GetMetricPosition( i , x, y ); 
      if( x < -10 ){
	grLYSmallGood->SetPoint( grLYSmallGood->GetN(), LY[i],Resolution[i]);
	hisResYDep->Fill( y, Resolution[i] );
      }else if( x > 10 ){
	grLYSmallBad->SetPoint( grLYSmallBad->GetN(), LY[i], Resolution[i] );
      }else{
	;
      }
    }else{
      grLarge    ->SetPoint( grLarge->GetN()    , Resolution[ i ], EConverter->GetCalibrationConstant( i ));    
      grLYLarge  ->SetPoint( grLYLarge->GetN()  , Resolution[ i ], LY[i]);
      grGainLarge->SetPoint( grGainLarge->GetN(), Resolution[ i ], Gain[i]);
    }
  }

  TCanvas* can = new TCanvas("can","",1600,800);  
  can->Divide( 4,2 );
  can->cd(1);
  imageRes->Draw("colz");
  can->cd(2);
  imageOut->Draw("colz");
  can->cd(3);
  gr->SetMarkerStyle( 6);
  gr->Draw("AP");
  can->cd(4);
  grSmall->Draw("AP");
  grLarge->Draw("P");
  can->cd(5);
  //grLYSmall->Draw("AP");
  grLYSmallGood->Draw("AP");
  grLYSmallGood->Fit("pol1","","",0,2);
  can->cd(6);
  //grLYLarge->Draw("AP");
  grLYSmallBad->Draw("AP");
  grLYSmallBad->Fit("pol1","","",0,2);
  can->cd(7);
  hisResYDep->Draw("colz");  
  TProfile* prof = hisResYDep->ProfileX();
  prof->SetLineColor( 1 );
  prof->SetLineWidth( 2 ); 
  prof->Draw("same");
  prof->Fit("pol1","","",-600,600);
  //grGainSmall->Draw("AP");
  can->cd(8);
  //grGainLarge->Draw("AP");
  app->Run();
  return 0; 

}
