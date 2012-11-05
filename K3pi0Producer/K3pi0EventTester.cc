#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <list>
#include <vector>

#include "cluster/Cluster.h"
#include "gamma/Gamma.h"
#include "pi0/Pi0.h"
#include "klong/Klong.h"

#include "gnana/E14GNAnaDataContainer.h"
#include "gnana/E14GNAnaFunction.h"
#include "rec2g/Rec2g.h"

#include "User_Function.h"

#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TApplication.h"

#include "CsIImage.h"
#include "IDHandler.h"

int main( int argc , char** argv ){

  TChain* ch  = new TChain("T");
  std::string ROOTFILE_WAV = std::getenv("ROOTFILE_WAV");

  Double_t ZHead = atof( argv[1] );
  Double_t ZTail = atof( argv[2] );

  //Int_t RunNumber =  atoi( argv[1] );

  TApplication* app = new TApplication("app",&argc, argv );

  std::ifstream ifs( "KLRunList_2.txt");
  int tmpRunNumber; 
  while(ifs >> tmpRunNumber ){
    ch->Add(Form("%s/run_wav_%d_cl.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
  }	    
	
  /*
  for( int i = 0; i< 10; i++){
    ch->Add(Form("%s/run_wav_%d_cl.root",ROOTFILE_WAV.c_str(),RunNumber +i));
  }
    
  */
  std::ifstream ifsPi0Peak("Data/Pi0Peak.dat");
  int tmpID;
  double tmpDelta;
  double tmpResolution;
  double Pi0Delta[2716]={0};
  double Pi0Resolution[2716]={0xFFFF};
  while( ifsPi0Peak >> tmpID >> tmpDelta >> tmpResolution ){
    if( tmpResolution > 0xFFFE || tmpResolution < 0 ){
      continue;
    }
    Pi0Delta[ tmpID ] = tmpDelta;
    Pi0Resolution[ tmpID ] = tmpResolution; 
  }


  E14GNAnaDataContainer data;
  data.setBranchAddress( ch );

  IDHandler* handler = new IDHandler();
  CsIImage* imageTime = new CsIImage(handler);
  CsIImage* imageEnergy = new CsIImage( handler); 
  

  double x,y;
  TGraph* grTimeCluster = new TGraph(); 
  grTimeCluster->SetMarkerStyle(6);
  TCanvas* can  = new TCanvas("can","",800,800);
  can->Divide(2,2);

  TFile* tfout  = new TFile(Form("hist_time_output_%d_%d.root",(int)ZHead,(int)ZTail), "recreate");
  TH2D* hist_R_Theta_Time[10];
  TH2D* hist_D_Theta_Time[10];
  TH2D* hist_L_Theta_Time[10];
  TH2D* hist_R_Theta_MeanTime[10];

  for( int i = 0; i< 10; i++){
    hist_R_Theta_Time[i] = new TH2D(Form("hist_R_Theta_Time_%d",i),Form("hist_R_Theta_Time_%d",i),
				    40,-200,200,100,-5,5);
    hist_D_Theta_Time[i] = new TH2D(Form("hist_D_Theta_Time_%d",i),Form("hist_D_Theta_Time_%d",i),
				    20,0,200,100,-5,5);
    hist_L_Theta_Time[i] = new TH2D(Form("hist_L_Theta_Time_%d",i),Form("hist_L_Theta_Time_%d",i),
				    20,0,200,100,-5,5);
    hist_R_Theta_MeanTime[i] = new TH2D(Form("hist_R_Theta_MeanTime_%d",i),Form("hist_R_Theta_MeanTime_%d",i),
					3,0,3,100,-5,5);
  }
 
  for( int ievent  = 0; ievent < ch->GetEntries(); ievent++){
    ch->GetEntry( ievent );
    std::list<Cluster> clist;
    std::list<Gamma>   glist;
    std::vector<Klong> klVec;
    data.getData( clist );
    data.getData( glist );
    data.getData( klVec );
    //std::cout<< clist.size() << " : " << glist.size() << " : " << klVec.size() << std::endl;
    if( klVec.size() ==  0){ continue; }
    if( clist.size() != 6 ){ continue; }    
    if( klVec[0].m() > 505 || klVec[0].m() < 485 ){ continue; }
    if( klVec[0].p3()[2] < ZHead || klVec[0].p3()[2] > ZTail) {continue; }

    for( std::list<Cluster>::iterator itCluster = clist.begin();
	 itCluster !=clist.end();
	 itCluster++){
      //std::cout<< klVec[0].p3()[2] << " : "<<  (*itCluster).id() << std::endl;      
      

      Double_t clusterPosX = -1*(*itCluster).pos()[0];
      Double_t clusterPosY = (*itCluster).pos()[1];
      Double_t clusterPosR = TMath::Sqrt( clusterPosX*clusterPosX + clusterPosY*clusterPosY );
      Double_t clex = clusterPosX/clusterPosR;
      Double_t cley = clusterPosY/clusterPosR;
      Int_t    Theta_Time_hist_ID = (int)(clusterPosR/100);
      Double_t ClusterCenterEnergy=0;
      //std::cout << clusterPosX  << " : " << clusterPosY << std::endl;
      if( clusterPosY > 0 ){
	if( clusterPosX > 100 ){ continue; }
      }else{
	if( clusterPosX >-100 ){ continue; }
      }
      Double_t ClusterMeanTime=0;
      Int_t    nCrystalTime  = 0;
      Int_t    MaxCHID;
      Double_t MaxCHTime;
      Double_t MaxEnergy = 0; 
      Double_t RadiusMax;
      Double_t Radius;
      grTimeCluster->Set(0);
      imageTime->Reset();
      imageEnergy->Reset();


      Double_t weightX=0.;
      Double_t weightY=0.;
      Double_t TotalEInCluster=0.;

      for( int iCrystal = 0; iCrystal < (*itCluster).clusterIdVec().size(); iCrystal++){
	Int_t    IDInCluster    = (*itCluster).clusterIdVec()[iCrystal]; 
	Double_t TimeInCluster  = (*itCluster).clusterTimeVec()[iCrystal]; 
	Double_t EnergyInCluster= (*itCluster).clusterEVec()[iCrystal];
	handler->GetMetricPosition( (*itCluster).clusterIdVec()[iCrystal] , x, y);
	//std::cout << (*itCluster).clusterIdVec()[iCrystal] << std::endl;
	if( (*itCluster).clusterEVec()[iCrystal] > MaxEnergy ){
	  MaxEnergy = (*itCluster).clusterEVec()[iCrystal];
	  MaxCHID   = (*itCluster).clusterIdVec()[iCrystal];
	  MaxCHTime = (*itCluster).clusterTimeVec()[iCrystal]-Pi0Delta[(*itCluster).clusterIdVec()[iCrystal]];
	  RadiusMax = TMath::Sqrt( x*x +y*y);
	}
	weightX  += x*EnergyInCluster;
	weightY  += y*EnergyInCluster;
	TotalEInCluster += EnergyInCluster;
      }

      weightX = weightX / TotalEInCluster;
      weightY = weightY / TotalEInCluster;
      if( weightX > -100 ){continue;}
      for( int iCrystal = 0; iCrystal < (*itCluster).clusterIdVec().size(); iCrystal++){
	Int_t    IDInCluster    = (*itCluster).clusterIdVec()[iCrystal]; 
	Double_t TimeInCluster  = (*itCluster).clusterTimeVec()[iCrystal]; 
	Double_t EnergyInCluster= (*itCluster).clusterEVec()[iCrystal];
	handler->GetMetricPosition( IDInCluster, x, y);
	if( TMath::Sqrt( ( x-weightX )*( x-weightX ) + ( y - weightY )*( y-weightY ) ) < 25 ){
	  ClusterCenterEnergy += EnergyInCluster;
	  if( EnergyInCluster > 30 && EnergyInCluster < 500 ){
	    ClusterMeanTime += TimeInCluster - Pi0Delta[ IDInCluster ]; 
	    nCrystalTime++;
	  }	  
	}
      }      
      if( ClusterCenterEnergy < 100 ){ continue; }
      if( nCrystalTime == 0){ continue; }
      if( TotalEInCluster < 100 ){ continue; }
      if( MaxEnergy  > 600 ){ continue ; }
      ClusterMeanTime /= nCrystalTime;
      Int_t RegionNumber=-1;// 0-2
      Double_t MeanTime[3]={0};
      Double_t MeanWeight[3]= {0};
      for( int iCrystal = 0; iCrystal < (*itCluster).clusterIdVec().size(); iCrystal++){
	Int_t    IDInCluster    = (*itCluster).clusterIdVec()[iCrystal]; 
	Double_t TimeInCluster  = (*itCluster).clusterTimeVec()[iCrystal]; 
	Double_t EnergyInCluster= (*itCluster).clusterEVec()[iCrystal];

	handler->GetMetricPosition( IDInCluster, x, y);
	Radius = TMath::Sqrt( x*x +y*y );
	double CenterRadius  = TMath::Sqrt(weightX*weightX + weightY*weightY);
	double CenterDistance= TMath::Sqrt( (weightX - x)*(weightX - x) + (weightY - y)*(weightY - y ));
	if( CenterDistance < 30 ){
	  RegionNumber = 1;
	}else if( Radius > CenterRadius ){
	  RegionNumber = 2;
	}else if( Radius < CenterRadius ){
	  RegionNumber = 0;
	} 

	MeanTime[RegionNumber]   += TimeInCluster*EnergyInCluster;
	MeanWeight[RegionNumber] += EnergyInCluster;

	Double_t RInCluster = clex*(x-clusterPosX) + cley*(y-clusterPosY);
	Double_t DInCluster = TMath::Abs(cley*(x-clusterPosX) - clex*(y-clusterPosY));
	Double_t LInCluster = TMath::Abs(TMath::Sqrt((x-clusterPosX)*(x-clusterPosX) + (y-clusterPosY)*(y-clusterPosY)));
	if( EnergyInCluster < 12 ) continue ;	
	grTimeCluster->SetPoint( grTimeCluster->GetN(), 
				 Radius-RadiusMax,
				 TimeInCluster - Pi0Delta[IDInCluster]);
	imageTime->Fill( IDInCluster , TimeInCluster- Pi0Delta[ IDInCluster ]);
	imageEnergy->Fill( IDInCluster , EnergyInCluster );
	
	/*
	if( TMath::Abs(RInCluster) < 50){
	  hist_D_Theta_Time[Theta_Time_hist_ID]->Fill( DInCluster, TimeInCluster - MaxCHTime);
	}

	if( DInCluster < 50 ){
	  hist_R_Theta_Time[Theta_Time_hist_ID]->Fill( RInCluster, TimeInCluster - MaxCHTime);
	}
	hist_L_Theta_Time[ Theta_Time_hist_ID]->Fill( LInCluster, TimeInCluster - MaxCHTime);
	*/

	if( TMath::Abs(RInCluster) < 50){
	  hist_D_Theta_Time[Theta_Time_hist_ID]->Fill( DInCluster, TimeInCluster - Pi0Delta[ IDInCluster] - ClusterMeanTime);
	}

	if( DInCluster < 50 ){
	  hist_R_Theta_Time[Theta_Time_hist_ID]->Fill( RInCluster, TimeInCluster - Pi0Delta[ IDInCluster] - ClusterMeanTime);
	}
	hist_L_Theta_Time[ Theta_Time_hist_ID]->Fill( LInCluster, TimeInCluster - Pi0Delta[ IDInCluster]- ClusterMeanTime );

      }

      for( int iRegion = 0; iRegion < 3; iRegion++){
	if( MeanWeight[iRegion] <= 0){ continue; }
	MeanTime[iRegion] /= MeanWeight[iRegion];
	//std::cout<< iRegion <<  " : " << MeanTime[iRegion] -ClusterMeanTime << std::endl ;
	hist_R_Theta_MeanTime[ Theta_Time_hist_ID ]->Fill( iRegion, MeanTime[iRegion] -ClusterMeanTime);
      }


      /*
      can->cd(1);
      gPad->SetLogz();
      imageEnergy->Draw("colz");
      can->cd(2);
      imageTime->Draw("colz");
      can->cd(3);
      grTimeCluster->Draw("AP");
      can->Update();
      can->Modified();
      getchar();
      */
    }
  }

  for( int i = 0; i< 10; i++){
    hist_D_Theta_Time[i]->Write();
    hist_R_Theta_Time[i]->Write();
    hist_L_Theta_Time[i]->Write();
    hist_R_Theta_MeanTime[i]->Write();
  }
  tfout->Close();
  
  //app->Run();
  return 0;
}
