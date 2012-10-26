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
#include "TMultiGraph.h"

int 
main( int argc ,char** argv ){
  std::string ANALIB = std::getenv("ANALYSISLIB");
  TApplication* app = new TApplication("app",&argc, argv );
  EnergyConverter* EConverter = new EnergyConverter();
  EConverter->ReadCalibrationRootFile(Form("%s/Data/Cosmic_Calibration_File/CosmicResult_20120209.root",ANALIB.c_str()));

  

  IDHandler* handler = new IDHandler();
  CsIImage* imageRes = new CsIImage( handler );
  CsIImage* imageOut = new CsIImage( handler );
  imageRes->SetTitle("Resolution[ns]");
  imageOut->SetTitle("Calibration Constant[ MeV / cnt ]");

  std::cout<< __PRETTY_FUNCTION__ << std::endl;

  std::ifstream ifs("Delta_Resolution_by_Energy.dat");
  double Resolution[ 2716 ][4]={{0xFFFF}}; 
  double Delta[2716][4]       ={{0.}};
  int tmpID;
  double tmpRes[4];
  double tmpDelta[4];
  
  while ( ifs >> tmpID >> tmpDelta[0] >> tmpDelta[1] >> tmpDelta[2] >> tmpDelta[3] >> tmpRes[0] >> tmpRes[1] >> tmpRes[2] >> tmpRes[3] ){
    for( int i = 0 ; i < 4; i++){
      Delta[ tmpID ][i]      = tmpDelta[i]; 
      Resolution[ tmpID ][i] = tmpRes[i];
    }
  }

  std::cout<< __PRETTY_FUNCTION__ << std::endl;

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

  std::cout<< __PRETTY_FUNCTION__ << std::endl;

  TGraph* gr            = new TGraph();
  TGraph* grSmall       = new TGraph();
  TGraph* grLarge       = new TGraph();
  TGraph* grLYSmall     = new TGraph();
  TGraph* grLYLarge     = new TGraph();
  TGraph* grGainSmall   = new TGraph();
  TGraph* grGainLarge   = new TGraph();
  TGraph* grLYSmallGood = new TGraph();
  TGraph* grLYSmallBad  = new TGraph();
  TGraph* grLYLargeGood = new TGraph();
  TGraph* grLYLargeBad  = new TGraph();

  //TH2D*   hisResYDep    = new TH2D("hisResYDep","hisResYDep1",48,-600,600,40,0,2);


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
  grLYLargeGood->SetMarkerStyle( 7 );
  grLYLargeBad->SetMarkerStyle( 7 );
  grLYLargeGood->SetMarkerColor( 4 );
  grLYLargeBad->SetMarkerColor( 6 );
  

  grGainSmall->SetMarkerStyle( 6 );
  grGainLarge->SetMarkerStyle( 6 );
  grGainSmall->SetMarkerColor( 2 );
  grGainLarge->SetMarkerColor( 3 );

  grLYSmall    ->SetNameTitle("grLYSmall"    ,"LYSmall;LY;Resolution");
  grLYLarge    ->SetNameTitle("grLYLarge"    ,"LYLarge;LY;Resolution");
  grLYSmallGood->SetNameTitle("grLYSmallGood","LYSmallGood;LY;Resolution");
  grLYSmallBad ->SetNameTitle("grLYSmallBad" ,"LYSmallBad;LY;Resolution");
  grLYLargeGood->SetNameTitle("grLYLargeGood","LYLargeGood;LY;Resolution");
  grLYLargeBad ->SetNameTitle("grLYLargeBad" ,"LYLargeBad;LY;Resolution");
  grGainSmall  ->SetNameTitle("grGainSmall"  ,"GainSmall;Resolution;Gain");
  grGainLarge  ->SetNameTitle("grGainLarge"  ,"GainLarge;Resolution;Gain");

  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  std::cout<< "LOOP" << std::endl;

  int IndexNumber = 2;
  double x,y;
  for( int i = 0; i< 2716; i++){
    std::cout << i << std::endl;
    if( Resolution[i][IndexNumber] >= 0xFFFF ){ continue; }
    gr->SetPoint(gr->GetN() , Resolution[ i ][IndexNumber], EConverter->GetCalibrationConstant( i ));    
    imageRes->Fill( i, Resolution[ i ][IndexNumber] ); 
    imageOut->Fill( i, EConverter->GetCalibrationConstant( i ));
    handler->GetMetricPosition( i , x, y );     
    if( i < 2240 ){
      grSmall    ->SetPoint( grSmall->GetN()    , Resolution[ i ][IndexNumber], EConverter->GetCalibrationConstant( i ));    
      grLYSmall  ->SetPoint( grLYSmall->GetN()  , Resolution[ i ][IndexNumber], LY[i] );
      grGainSmall->SetPoint( grGainSmall->GetN(), Resolution[ i ][IndexNumber], Gain[i] );

      if( x < -10 ){
	grLYSmallGood->SetPoint( grLYSmallGood->GetN(), LY[i],Resolution[i][IndexNumber]);
	//hisResYDep->Fill( y, Resolution[i][IndexNumber] );
      }else if( x > 10 ){
	grLYSmallBad->SetPoint( grLYSmallBad->GetN(), LY[i]/2, Resolution[i][IndexNumber] );
      }else{
	;
      }
    }else{
      if( y <-600 || y > 600 ){ 
	;
      }else{
	if( x < -10 ){
	  grLYLargeGood->SetPoint( grLYLargeGood->GetN(), 2*LY[i],Resolution[i][IndexNumber]);
	}else if( x > 10 ){
	  grLYLargeBad->SetPoint( grLYLargeBad->GetN(), 2*LY[i], Resolution[i][IndexNumber] );
	}else{
	  ;
	} 
      }
      grLarge    ->SetPoint( grLarge->GetN()    , Resolution[ i ][IndexNumber], EConverter->GetCalibrationConstant( i ));    
      grLYLarge  ->SetPoint( grLYLarge->GetN()  , Resolution[ i ][IndexNumber], LY[i]);
      grGainLarge->SetPoint( grGainLarge->GetN(), Resolution[ i ][IndexNumber], Gain[i]);
    }
  }

  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  std::cout<< grLarge->GetN() << std::endl;

  TCanvas* can = new TCanvas("can","",1200,800);  
  can->Divide( 3,2 );
  can->cd(1);
  imageRes->Draw("colz");
  can->cd(2);
  imageOut->Draw("colz");
  can->cd(3);
  gr->SetMarkerStyle( 6);
  //gr->Draw("AP");
  grSmall->Draw("AP");
  grLarge->Draw("P");
  can->cd(4);
  //grLYSmall->Draw("AP");
  grLYSmallGood->Draw("AP");
  grLYSmallBad->Draw("P");
  grLYSmallGood->GetXaxis()->SetRangeUser(0,4);
  can->cd(5);
  grLYLargeGood->Draw("AP");
  grLYLargeBad->Draw("P");
  grLYLargeGood->GetXaxis()->SetRangeUser(0,4);

  //grLYLarge->Draw("AP");
  //grLYSmallBad->Fit("pol1","","",0,2);
  can->cd(6);
  gPad->SetGridx();
  gPad->SetGridy();
  TMultiGraph* mgr  = new TMultiGraph();
  mgr->Add(grLYSmallGood);
  mgr->Add(grLYSmallBad);
  mgr->Add(grLYLargeGood);
  mgr->Add(grLYLargeBad);
  mgr->Draw("AP");
  can->cd(8);

  /*
  can->cd(7);
  hisResYDep->Draw("colz");  

  TProfile* prof = hisResYDep->ProfileX();
  prof->SetLineColor( 1 );
  prof->SetLineWidth( 2 ); 
  prof->Draw("same");
  //prof->Fit("pol1","","",-600,600);
  //grGainSmall->Draw("AP");
  can->cd(8);
  //grGainLarge->Draw("AP");
  */
  app->Run();
  return 0; 
}
