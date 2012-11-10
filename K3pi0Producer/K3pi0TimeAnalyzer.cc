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

  TApplication* app = new TApplication("app",&argc, argv );

  std::ifstream ifs( "KLRunList_2.txt");
  int tmpRunNumber; 
  int nFiles = 0;
  while(ifs >> tmpRunNumber ){
    ch->Add(Form("%s/run_wav_%d_Cal.root",ROOTFILE_WAV.c_str(),tmpRunNumber));
    nFiles++; 
    //    if( nFiles > 10 ){ break; }
  }	    
	
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

  Double_t sol = 299.792458;// [mm/ns]

  double x,y;

  TCanvas* can  = new TCanvas("can","",800,800);
  can->Divide(2,2);

  TFile* tfout  = new TFile(Form("hist_KlongZ_TimeDelta.root"),"recreate");
  TH2D* hisKlongZTimeDelta    = new TH2D("hisKlongZTimeDelta","hisKlongZTimeDelta", 300,3000,6000,200,-5,5);
  TH2D* hisKlongZTimeDeltaCal = new TH2D("hisKlongZTimeDeltaCal","hisKlongZTimeDeltaCal", 300,3000,6000,200,-5,5);
  // Fill by 10Degrees // 
  TH2D* hisKlongZGammaInjectDirection[5] ;
  TH2D* hisKlongZGammaEnergy[5];
  for( int i = 0; i< 5; i++){
    hisKlongZGammaInjectDirection[i] =  new TH2D(Form("hisKlongGammaInjectionDirection%d",i),
						 Form("hisKlongGammaInjectionDirection%d",i),
						 300, 3000,6000,200,-5,5);
    hisKlongZGammaEnergy[i] = new TH2D(Form("hisKlongZGammaEnergy%d",i),
				       Form("hisKlongZGammaEnergy%d",i),
				       300,3000,6000,200,-5,5);
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
 

    Double_t InnerR= 1000.;//mm
    Double_t OuterR= 0.;//mm
    Double_t InnerTime = 0.;
    Double_t OuterTime = 0.;
    Double_t InnerRx;
    Double_t InnerRy;
    Double_t OuterRx;
    Double_t OuterRy;

    Int_t GammaID = 0;
    Int_t MinRGammaID;
    for( std::list<Gamma>::iterator itGamma = glist.begin();
	 itGamma != glist.end() ;
	 itGamma++,GammaID++){
      if( GammaID >=6 ){ 
	GammaID = 0; 
	break; 
      }
      Double_t gx,gy,tmpR;
      gx = (*itGamma).x();
      gy = (*itGamma).y();
      tmpR = TMath::Sqrt( gx*gx + gy*gy );
      if( tmpR > OuterR ){
	OuterR    = tmpR;
	OuterTime = (*itGamma).t();
	OuterRx = gx;
	OuterRy = gy;
      }
      if( tmpR < InnerR ){
	MinRGammaID = GammaID;
	InnerR = tmpR;
	InnerTime = (*itGamma).t();
	InnerRx = gx;
	InnerRy = gy;
      }      
    }
    
    if( InnerR > 300 ){ continue; }
    if( InnerR < 150 ){ continue; }
    if( OuterR < 300 ){ continue; }

    Double_t KlongZ = klVec[0].vz();
    Double_t CsISurface = 6148;
    Double_t DeltaZ = CsISurface - KlongZ;
    Double_t CalInnerR = TMath::Sqrt( ( InnerRx - klVec[0].vx() )*( InnerRx - klVec[0].vx() ) + ( InnerRy - klVec[0].vy() )*( InnerRy - klVec[0].vy() ));
    Double_t CalOuterR = TMath::Sqrt( ( OuterRx - klVec[0].vx() )*( OuterRx - klVec[0].vx() ) + ( OuterRy - klVec[0].vy() )*( OuterRy - klVec[0].vy() ));

    Double_t DeltaLengthInner = TMath::Sqrt( (KlongZ -CsISurface)*(KlongZ - CsISurface ) + InnerR*InnerR) - TMath::Abs( CsISurface - KlongZ );
    Double_t DeltaLengthOuter = TMath::Sqrt( (KlongZ -CsISurface)*(KlongZ - CsISurface ) + OuterR*OuterR) - TMath::Abs( CsISurface - KlongZ );
    Double_t DeltaLengthInnerCal = TMath::Sqrt( (KlongZ -CsISurface)*(KlongZ - CsISurface ) + CalInnerR*CalInnerR) - TMath::Abs( CsISurface - KlongZ );
    Double_t DeltaLengthOuterCal = TMath::Sqrt( (KlongZ -CsISurface)*(KlongZ - CsISurface ) + CalOuterR*CalOuterR) - TMath::Abs( CsISurface - KlongZ );

    //std::cout << KlongZ << std::endl;
    Double_t TimeZero = InnerTime - DeltaLengthInner/sol;
    Double_t DeltaTime = OuterTime - DeltaLengthOuter/sol - TimeZero; 
    Double_t TimeZeroCal  = InnerTime - DeltaLengthInnerCal/sol;
    Double_t DeltaTimeCal = OuterTime - DeltaLengthOuterCal/sol - TimeZeroCal; 

    hisKlongZTimeDelta->Fill( KlongZ, DeltaTime );
    hisKlongZTimeDeltaCal->Fill( KlongZ, DeltaTimeCal );


    // Search minimum/Maximum Radius Gamma    
    for( std::list<Gamma>::iterator itGamma = glist.begin();
	 itGamma != glist.end(); 
	 itGamma++, GammaID++){
      if( GammaID >=6  ) { 
	GammaID = 0; 
	break; 
      }
      if( GammaID == MinRGammaID ){ continue; }
      Double_t gx, gy, tmpR, calR, Theta, Energy, GammaTime;
      Double_t DeltaLength;
      gx      = (*itGamma).x();
      gy      = (*itGamma).y();
      Energy  = (*itGamma).e();
      tmpR    = TMath::Sqrt( gx*gx + gy*gy );
      calR    = TMath::Sqrt( (gx - klVec[0].vx())*(gx - klVec[0].vx()) + (gy - klVec[0].vy())*(gy - klVec[0].vy()));
      DeltaLength = TMath::Abs(TMath::Sqrt( (KlongZ - CsISurface)*(KlongZ - CsISurface) + calR*calR ) - TMath::Abs(CsISurface - KlongZ));
      GammaTime   = (*itGamma).t() - DeltaLength/sol - TimeZeroCal;
      Theta       = TMath::ATan(calR/TMath::Abs(KlongZ - CsISurface))/TMath::Pi()*180;
      
      Int_t ThetaID = (int)(Theta/10);
      Int_t EnergyID = (int)(Energy/100);
      if( EnergyID < 5 ){
	hisKlongZGammaEnergy[EnergyID]->Fill(KlongZ,GammaTime);
      }
      if( ThetaID < 5 ){
	hisKlongZGammaInjectDirection[ ThetaID ]->Fill(KlongZ,GammaTime);
      }
    }

    
    /*
    // Cluster  Analysis // 
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
      imageTime->Reset();
      imageEnergy->Reset();
      Double_t weightX=0.;
      Double_t weightY=0.;
      Double_t TotalEInCluster=0.;      
    }
    */

  }// Event Process //

  hisKlongZTimeDelta->Write();
  hisKlongZTimeDeltaCal->Write();
  for( int i = 0; i< 5; i++){
    hisKlongZGammaEnergy[i]->Write();
    hisKlongZGammaInjectDirection[i]->Write();
  }
  tfout->Close();
  
  //app->Run();
  return 0;
}
