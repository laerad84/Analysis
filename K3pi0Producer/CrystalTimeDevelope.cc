#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TFile.h"
#include "TF1.h"
#include "TProfile.h"

//void CrystalTimeDevelope(){
int main( int argc , char** argv ){
  
  int ZHead = atoi( argv[1] );
  int ZTail = atoi( argv[2] );
  
  TFile* tf = new TFile(Form("hist_time_output_%d_%d.root",ZHead,ZTail));
  TH2D* hist_D_Theta_Time[10];
  TH2D* hist_R_Theta_Time[10];
  TH2D* hist_L_Theta_Time[10];
  
  for( int i = 0; i < 10; i++){
    hist_D_Theta_Time[i] = (TH2D*)tf->Get(Form("hist_D_Theta_Time_%d",i));
    hist_R_Theta_Time[i] = (TH2D*)tf->Get(Form("hist_R_Theta_Time_%d",i));
    hist_L_Theta_Time[i] = (TH2D*)tf->Get(Form("hist_L_Theta_Time_%d",i));
  }
  TH1D* his_Slice_R_Plus[10][10];
  TH1D* his_Slice_R_Minus[10][10];
  TH1D* his_Slice_D[10][10];
  TH1D* his_Slice_L[10][10];
  for( int i = 0; i< 10; i++){
    for( int j = 0; j< 10; j++){
      his_Slice_R_Plus[i][j] = new TH1D(Form("his_Slice_R_PLUS_%d_%d",i,j),
					Form("his_Slice_R_PLUS_%d_%d",i,j),
					100,-5,5);
      his_Slice_R_Minus[i][j] = new TH1D(Form("his_Slice_R_Minus_%d_%d",i,j),
					 Form("his_Slice_R_Minus_%d_%d",i,j),
					 100,-5,5);
    }
    
    for( int j = 0; j< 10; j++){
      his_Slice_D[i][j] = new TH1D(Form("his_Slice_D_%d_%d",i,j),
				   Form("his_Slice_D_%d_%d",i,j),
				   100,-5,5);
      his_Slice_L[i][j] = new TH1D(Form("his_Slice_L_%d_%d",i,j),
				   Form("his_Slice_L_%d_%d",i,j),
				   100,-5,5);
    }
  }    
  TPostScript* ps = new TPostScript(Form("output_%d_%d.ps",ZHead,ZTail),111);
  ps->NewPage();
  TCanvas* can = new TCanvas("can","can",600,800);
  TGraphErrors* gr_Time_R_Plus[10];
  TGraphErrors* gr_Time_R_Minus[10]; 
  TGraphErrors* gr_Time_D[10];
  TGraphErrors* gr_Time_L[10];
  for( int index = 0; index < 10 ; index++){
    gr_Time_R_Plus[index] = new TGraphErrors();
    gr_Time_R_Minus[index] = new TGraphErrors();
    gr_Time_D[index] = new TGraphErrors();
    gr_Time_L[index] = new TGraphErrors();
  }

  for (int index = 0; index <10 ; index++){

    for( int ibin  = 0; ibin < 10; ibin++){
      his_Slice_R_Plus[index][ibin] = hist_R_Theta_Time[index]->ProjectionY(his_Slice_R_Plus[index][ibin]->GetName(),20+ibin+1,20+ibin+2);
      his_Slice_R_Minus[index][ibin] = hist_R_Theta_Time[index]->ProjectionY(his_Slice_R_Minus[index][ibin]->GetName(),20-ibin-1,20-ibin);
      his_Slice_D[index][ibin] = hist_R_Theta_Time[index]->ProjectionY(his_Slice_D[index][ibin]->GetName(),ibin+1,ibin+2);
      his_Slice_L[index][ibin] = hist_L_Theta_Time[index]->ProjectionY(his_Slice_L[index][ibin]->GetName(),ibin+1,ibin+2);
    }
    
    for( int ibin = 0; ibin < 10 ; ibin++){
      double mean_RP = his_Slice_R_Plus[index][ibin]->GetMean();
      double RMS_RP  = his_Slice_R_Plus[index][ibin]->GetRMS();
      if( his_Slice_R_Plus[index][ibin]->GetEntries() <10 ){ continue; }      
      his_Slice_R_Plus[index][ibin]->Fit("gaus","Q","",
					 mean_RP-RMS_RP, mean_RP+RMS_RP);
      TF1* func  = his_Slice_R_Plus[index][ibin]->GetFunction("gaus");
      double mean_RP_Fit = func->GetParameter(1);
      double RMS_RP_Fit  = func->GetParameter(2);
      his_Slice_R_Plus[index][ibin]->Fit("gaus","Q","",mean_RP_Fit-RMS_RP_Fit,
					 mean_RP_Fit+RMS_RP_Fit);
      func = his_Slice_R_Plus[index][ibin]->GetFunction("gaus");
      gr_Time_R_Plus[index]->SetPoint( gr_Time_R_Plus[index]->GetN(),ibin*10+5 , func->GetParameter(1));
      gr_Time_R_Plus[index]->SetPointError( gr_Time_R_Plus[index]->GetN()-1,
					    0, func->GetParameter( 2 ));
      can->Clear();
      his_Slice_R_Plus[index][ibin]->Draw();
      can->Modified();      
      can->Update();
      //ps->NewPage();
      
    }

    for( int ibin = 0; ibin < 10 ; ibin++){
      double mean_RP = his_Slice_R_Minus[index][ibin]->GetMean();
      double RMS_RP  = his_Slice_R_Minus[index][ibin]->GetRMS();
      if( his_Slice_R_Minus[index][ibin]->GetEntries() <10 ){ continue; }

      his_Slice_R_Minus[index][ibin]->Fit("gaus","Q","",
					 mean_RP-RMS_RP, mean_RP+RMS_RP);
      TF1* func  = his_Slice_R_Minus[index][ibin]->GetFunction("gaus");
      double mean_RP_Fit = func->GetParameter(1);
      double RMS_RP_Fit  = func->GetParameter(2);
      his_Slice_R_Minus[index][ibin]->Fit("gaus","Q","",mean_RP_Fit-RMS_RP_Fit,
					 mean_RP_Fit+RMS_RP_Fit);
      func = his_Slice_R_Minus[index][ibin]->GetFunction("gaus");
      gr_Time_R_Minus[index]->SetPoint( gr_Time_R_Minus[index]->GetN(),ibin*10+5 , func->GetParameter(1));
      gr_Time_R_Minus[index]->SetPointError( gr_Time_R_Minus[index]->GetN()-1,
					  0, func->GetParameter( 2 ));
      can->Clear();
      his_Slice_R_Minus[index][ibin]->Draw();
      can->Modified();
      can->Update();
      //ps->NewPage();
    }

    for( int ibin = 0; ibin < 10 ; ibin++){
      double mean_DP = his_Slice_D[index][ibin]->GetMean();
      double RMS_DP  = his_Slice_D[index][ibin]->GetRMS();
      if( his_Slice_D[index][ibin]->GetEntries() <10 ){ continue; }

      his_Slice_D[index][ibin]->Fit("gaus","Q","",
					 mean_DP-RMS_DP, mean_DP+RMS_DP);
      TF1* func  = his_Slice_D[index][ibin]->GetFunction("gaus");
      double mean_DP_Fit = func->GetParameter(1);
      double RMS_DP_Fit  = func->GetParameter(2);
      his_Slice_D[index][ibin]->Fit("gaus","Q","",mean_DP_Fit-RMS_DP_Fit,
					 mean_DP_Fit+RMS_DP_Fit);
      func = his_Slice_D[index][ibin]->GetFunction("gaus");
      gr_Time_D[index]->SetPoint( gr_Time_D[index]->GetN(),ibin*10+5 , func->GetParameter(1));
      gr_Time_D[index]->SetPointError( gr_Time_D[index]->GetN()-1,
					  0, func->GetParameter( 2 ));
      can->Clear();
      his_Slice_D[index][ibin]->Draw();
      can->Modified();
      can->Update();
      //ps->NewPage();
    }

    for( int ibin = 0; ibin < 10 ; ibin++){
      double mean_LP = his_Slice_L[index][ibin]->GetMean();
      double RMS_LP  = his_Slice_L[index][ibin]->GetRMS();
      if( his_Slice_L[index][ibin]->GetEntries() <10 ){ continue; }

      his_Slice_L[index][ibin]->Fit("gaus","Q","",
					 mean_LP-RMS_LP, mean_LP+RMS_LP);
      TF1* func  = his_Slice_L[index][ibin]->GetFunction("gaus");
      double mean_LP_Fit = func->GetParameter(1);
      double RMS_LP_Fit  = func->GetParameter(2);
      his_Slice_L[index][ibin]->Fit("gaus","Q","",mean_LP_Fit-RMS_LP_Fit,
					 mean_LP_Fit+RMS_LP_Fit);
      func = his_Slice_L[index][ibin]->GetFunction("gaus");
      gr_Time_L[index]->SetPoint( gr_Time_L[index]->GetN(),ibin*10+5 , func->GetParameter(1));
      gr_Time_L[index]->SetPointError( gr_Time_L[index]->GetN()-1,
					  0, func->GetParameter( 2 ));
      can->Clear();
      his_Slice_L[index][ibin]->Draw();
      can->Modified();
      can->Update();

      //ps->NewPage();
    }

  }


  for( int i = 0; i< 10; i++){
    gPad->SetLogz(1);
    hist_R_Theta_Time[i]->Draw("colz");
    hist_R_Theta_Time[i]->ProfileX()->Draw("same");
    can->Update();
    can->Modified();
  }
  for( int i = 0; i< 10; i++){
    gPad->SetLogz(1);
    hist_D_Theta_Time[i]->Draw("colz");
    hist_D_Theta_Time[i]->ProfileX()->Draw("same");
    can->Update();
    can->Modified();
  }
  for( int i = 0; i< 10; i++){
    gPad->SetLogz(1);
    hist_L_Theta_Time[i]->Draw("colz");
    hist_L_Theta_Time[i]->ProfileX()->Draw("same");
    can->Update();
    can->Modified();
  }
  gPad->SetLogz(0);


  can->Clear();
  ps->NewPage();
  can->Divide( 4,3);
  for( int index = 0; index< 10; index++){
    can->cd(index+1);
    gr_Time_R_Plus[index]->SetMarkerStyle(8);
    gr_Time_R_Plus[index]->SetMarkerColor( 1 );
    gr_Time_R_Plus[index]->Draw("AP");
    gr_Time_R_Minus[index]->SetMarkerStyle(7);
    gr_Time_R_Minus[index]->SetMarkerColor( 2 );
    gr_Time_R_Minus[index]->Draw("P");
  }
  can->Update();
  can->Modified();
  ps->NewPage();
  for( int index = 0; index< 10; index++){
    can->cd(index+1);
    gr_Time_D[index]->SetMarkerStyle(8);
    gr_Time_D[index]->SetMarkerColor( 1 );
    gr_Time_D[index]->Draw("AP");
    gr_Time_L[index]->SetMarkerStyle(7);
    gr_Time_L[index]->SetMarkerColor( 2 );
    gr_Time_L[index]->Draw("P");
  }
  can->Update();
  can->Modified();
  ps->NewPage();
  ps->Close();
  
}


				   

  
