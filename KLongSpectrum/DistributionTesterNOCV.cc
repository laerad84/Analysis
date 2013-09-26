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

const double KLMass = 497.648;//MeV

double KLSpectrum(double* x,double*p){
  //double klp   = TMath::Sqrt( x[0]*x[0]-KLMass*KLMass);  
  double klp = x[0];
  double sigma = p[2]*(1-(p[3]+p[4]*klp)*(klp-p[1]));
  double value = p[0]*TMath::Exp(-0.5*(klp-p[1])*(klp-p[1])/(sigma*sigma));
  return value;
}


//void DistributionTester(){
int main( int argc, char** argv){

  TF1* soltFunc  = new TF1("soltfunc",KLSpectrum,0,12,5);
  TF1* sugarFunc = new TF1("sugarfunc",KLSpectrum,0,12,5);
  const double soltPar[5] = { 1, 1.37,7.48e-1,-2.9e-1,1.68e-2};
  const double sugarPar[5] = {1,1.41991,0.810237,-0.301413,0.017092};
  soltFunc->SetParameters(soltPar);
  sugarFunc->SetParameters(sugarPar);

  const int nFile = 2;
  TFile* tf[nFile]; 
  TTree* tr[nFile];  
  char* name[nFile] = {"SIMFULL","WAVNOCV"};
  
  for( int i = 0; i < nFile; i++){
    tf[i] = new TFile(Form("Kl_Total_%s.root",name[i]));
    tr[i] = (TTree*)tf[i]->Get(Form("trKL"));
  }
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];

  TH1D* hisKLP[nFile];
  TH1D* hisKLZ[nFile];
  TH1D* hisKLE[nFile];
  TH1D* hisKLMass[nFile];
  TH1D* hisKLChisq[nFile];
  TH1D* hisKLSecChisq[nFile];
  TH1D* hisKLX[nFile];
  TH1D* hisKLY[nFile];

  TH1D* hisGammaE[nFile];
  TH1D* hisGammaTime[nFile];
  TH1D* hisGammaPosX[nFile];
  TH1D* hisGammaPosY[nFile];
  TH1D* hisKLZAcceptance[nFile][10];
  TH1D* hisKLZPosition[nFile];
  TH1D* hisVETOHistory[nFile];
  TH1D* hisPi0E[nFile];
  TH1D* hisPi0Mass[nFile];
  TH1D* hisPi0Z[nFile];
  TH1D* hisPi0ZSig[nFile];
  TH1D* hisPi0P[nFile];
  TH1D* hisPi0Pt[nFile];
  TH1D* hisL1TrigCount[nFile][20];

  TH1D* hisKLEN[nFile];
  TH1D* hisKLMassN[nFile];
  TH1D* hisKLXN[nFile];
  TH1D* hisKLYN[nFile];
  TH1D* hisPi0MassN[nFile];
  TH1D* hisKLChisqN[nFile];
  TH1D* hisKLSecChisqN[nFile];
  
  for( int i = 0; i< nFile; i++){
    for( int j = 0; j< 10; j++){
      hisKLZAcceptance[i][j] = new TH1D(Form("hisKLZAcceptance_%d_%d",i,j),Form("hisKLZAcceptance_%d_%d_%d",i,500*j,500*(j+1)),70,0,7000);
    }
    hisKLZPosition[i]= new TH1D(Form("hisKLZPosition_%d",i),Form("hisKLZPosition_%s",name[i]),70,0,7000);
    hisKLZ[i]        = new TH1D(Form("hisKLZ_%d",i),Form("hisKLZ_%s",name[i]),70,0,7000);
    hisKLP[i]        = new TH1D(Form("hisKLP_%d",i),Form("hisKLP_%s",name[i]),100,0,10000);
    hisKLE[i]        = new TH1D(Form("hisKLE_%d",i),Form("hisKLE_%s",name[i]),100,0,10000);
    hisKLX[i]        = new TH1D(Form("hisKLX_%d",i),Form("hisKLX_%s",name[i]),160,-400,400);
    hisKLY[i]        = new TH1D(Form("hisKLY_%d",i),Form("hisKLY_%s",name[i]),160,-400,400);
    hisKLMass[i]     = new TH1D(Form("hisKLMass_%d",i),Form("hisKLMass_%s",name[i]),100,450,550);
    hisKLChisq[i]    = new TH1D(Form("hisKLChisq_%d",i),Form("hisKLChisq_%s",name[i]),100,0,100);
    hisKLSecChisq[i] = new TH1D(Form("hisKLSecChisq_%d",i),Form("hisKLSecChisq_%s",name[i]),100,0,100);
    hisGammaE[i]     = new TH1D(Form("hisGammaE_%d",i),Form("hisGammaE_%s",name[i]),100,0,2000);
    hisGammaTime[i]  = new TH1D(Form("hisGammaTime_%d",i),Form("hisGammaTime_%s",name[i]),100,0,400);
    hisGammaPosX[i]  = new TH1D(Form("hisGammaPosX_%d",i),Form("hisGammaPosX_%s",name[i]),80,-1000,1000);
    hisGammaPosY[i]  = new TH1D(Form("hisGammaPosY_%d",i),Form("hisGammaPosY_%s",name[i]),80,-1000,1000);
    hisVETOHistory[i]= new TH1D(Form("hisVETOHistory_%d",i),Form("hisVETOHistory_%s",name[i]),20,0,20);
    hisPi0E[i]       = new TH1D(Form("hisPi0E_%d",i),Form("hisPi0E_%s",name[i]),100,0,3000);
    hisPi0P[i]       = new TH1D(Form("hisPi0P_%d",i),Form("hisPi0P_%s",name[i]),100,0,3000);
    hisPi0Pt[i]      = new TH1D(Form("hisPi0Pt_%d",i),Form("hisPi0Pt_%s",name[i]),100,0,400);
    hisPi0Mass[i]    = new TH1D(Form("hisPi0Mass_%d",i),Form("hisPi0Mass_%s",name[i]),150,50,200);
    hisPi0Z[i]       = new TH1D(Form("hisPi0Z_%d",i),Form("hisPi0Z_%s",name[i]),70,0,7000);
    hisPi0ZSig[i]    = new TH1D(Form("hisPi0ZSig_%d",i),Form("hisPi0ZSig_%s",name[i]),150,0,15000);
    for( int k = 0; k<20; k++){
	hisL1TrigCount[i][k]= new TH1D(Form("hisL1TrigCount_%d_%d",i,k),Form("hisL1TrigCount_%s_%d",name[i],k),500,0,30000);
    }
    hisKLEN[i]        = new TH1D(Form("hisKLEN_%d",i),Form("hisKLEN_%s",name[i]),100,0,10000);
    hisKLMassN[i]     = new TH1D(Form("hisKLMassN_%d",i),Form("hisKLMassN_%s",name[i]),100,450,550);
    hisKLXN[i]        = new TH1D(Form("hisKLXN_%d",i),Form("hisKLXN_%s",name[i]),160,-400,400);
    hisKLYN[i]        = new TH1D(Form("hisKLYN_%d",i),Form("hisKLYN_%s",name[i]),160,-400,400);
    hisPi0MassN[i]    = new TH1D(Form("hisPi0MassN_%d",i),Form("hisPi0MassN_%s",name[i]),150,50,200);
    hisKLChisqN[i]    = new TH1D(Form("hisKLChisqN_%d",i),Form("hisKLChisqN_%s",name[i]),100,0,100);
    hisKLSecChisqN[i] = new TH1D(Form("hisKLSecChsiqN_%d",i),Form("hisKLSecChisqN_%s",name[i]),100,0,100);    
  }

  TFile* tfOut = new TFile("DistributionTesterNOCV.root","recreate");
  for( int iFile = 0; iFile < nFile; iFile++){
    E14GNAnaDataContainer data;
    data.setBranchAddress( tr[iFile] );
    tr[iFile]->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
    tr[iFile]->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);

    for( int ievent = 0; ievent < tr[iFile]->GetEntries(); ievent++){      
      tr[iFile]->GetEntry(ievent);
      hisVETOHistory[iFile]->Fill(0);
      //if( ievent  >= 100000 ){ break ; } 
      std::list<Cluster> clist;
      std::list<Gamma>   glist;
      std::vector<Klong> klVec;
      data.getData(clist);
      data.getData(glist);
      data.getData(klVec);

      Double_t klptSq  = klVec[0].p3()[0]*klVec[0].p3()[0]+klVec[0].p3()[1]*klVec[0].p3()[2];
      Double_t klMom   = TMath::Sqrt(klVec[0].e()*klVec[0].e()-klVec[0].m()*klVec[0].m());
      //Double_t Ratio   = 1; 
      Double_t Ratio   = sugarFunc->Eval(klMom)/soltFunc->Eval(klMom);

      if( iFile != 2 ){
	if( CsiL1nTrig< 5 ){ continue; }
      }
      hisVETOHistory[iFile]->Fill(1);
      hisVETOHistory[iFile]->Fill(2);
      bool  bInnerGamma = false;
      bool  bOuterGamma = false;
      bool  bClusterMaxE= false;
      // Cut on Gamma // 
      Double_t MinGammaE = 200;
      Int_t nEGamma = 0;      
      Double_t ClusterMaxE = 0;
      std::list<Gamma>::iterator git = glist.begin();
      for( int igamma = 0; igamma < 6; igamma++,git++){

	if( (*git).e() > MinGammaE ){ nEGamma++;}
	if( (*git).clusterEVec()[0] > ClusterMaxE){
	  ClusterMaxE = (*git).clusterEVec()[0];
	  if(ClusterMaxE > 400 ){ bClusterMaxE = true; }
	}
	if( TMath::Abs((*git).x()) < 150 &&
	    TMath::Abs((*git).y()) < 150 ){
	  bInnerGamma = true; 	  
	}
	if( TMath::Sqrt((*git).x()*(*git).x()+(*git).y()*(*git).y()) > 850 ){
	  bOuterGamma = true;
	}
	if( TMath::Abs((*git).y())>550){ bOuterGamma = true; }
      }


      /////////////////////////////
      /// Basic Cut condition   ///
      //if( bClusterMaxE){ continue; }
      hisVETOHistory[iFile]->Fill(4);
      if( nEGamma < 6 ){ continue; }
      hisVETOHistory[iFile]->Fill(5);
      if( bInnerGamma ){ continue; }
      hisVETOHistory[iFile]->Fill(6);
      if( bOuterGamma ){ continue; }
      hisVETOHistory[iFile]->Fill(7);

      hisKLChisq[iFile]->Fill(klVec[0].chisqZ());
      if( klVec.size() >1 ){
	hisKLSecChisq[iFile]->Fill(klVec[1].chisqZ());
      }

      if( klVec[0].chisqZ() > 15 ){continue;} 
      //////////////////////////////
      bool bEGammaN = false;
      bool bKLZ     = false;
      bool bKLP     = false;
      bool bPi0Pt   = false; 
      for( git = glist.begin();
	   git != glist.end();
	   git++){
	if( (*git).e() > 1300 ){ bEGammaN = true; }
      }
      std::vector<Pi0>::iterator pit = klVec[0].pi0().begin();      
      if( klVec[0].p3().mag() > 4500 || klVec[0].p3().mag() < 1500 ){
	bKLP = true;
      }
      if( klVec[0].vz() > 4500 || klVec[0].vz() < 2000 ){
	bKLZ = true; 
      }
      for(pit= klVec[0].pi0().begin();
	  pit!=klVec[0].pi0().end();
	  pit++){
	if( (*pit).p3().perp() > 140 ){ 
	  bPi0Pt = true;
	}
      }
      
      for( int il1 = 0; il1 < 20; il1++){
	hisL1TrigCount[iFile][il1]->Fill(CsiL1TrigCount[il1]);
      }

      for(git= glist.begin();
	  git!= glist.end();
	  git++){
	hisGammaE[iFile]->Fill((*git).e());
	hisGammaTime[iFile]->Fill((*git).t());
	hisGammaPosX[iFile]->Fill((*git).x());
	hisGammaPosY[iFile]->Fill((*git).y());
      }
      for(pit= klVec[0].pi0().begin();
	  pit!=klVec[0].pi0().end();
	  pit++){
	hisPi0E[iFile]->Fill((*pit).e());
	hisPi0Z[iFile]->Fill((*pit).vz());
	hisPi0P[iFile]->Fill((*pit).p3().mag());
	hisPi0Pt[iFile]->Fill((*pit).p3().perp());
	hisPi0Mass[iFile]->Fill((*pit).m());
	hisPi0ZSig[iFile]->Fill((*pit).recZsig2());	
      }
      
      /*
      if( klVec[0].vz() < 5000 && klVec[0].vz() > 3000){
	if(TMath::Abs(klVec[0].m()-KLMass )< 10 ){
	  hisKLP[iFile]->Fill(klMom);
	}
      }
      */
      hisKLP[iFile]->Fill(klVec[0].p3().mag());
      hisKLE[iFile]->Fill(klVec[0].e());
      hisKLZ[iFile]->Fill(klVec[0].vz());
      hisKLX[iFile]->Fill(klVec[0].vx());
      hisKLY[iFile]->Fill(klVec[0].vy());
      hisKLMass[iFile]->Fill(klVec[0].m());

      if( TMath::IsNaN(Ratio) == 0 ){ 
	hisKLZPosition[iFile]->Fill(klVec[0].vz(),Ratio);
      }



      if( !bEGammaN && !bKLZ && !bKLP && !bPi0Pt ){

	hisKLEN[iFile]->Fill(klVec[0].e());
	hisKLMassN[iFile]->Fill(klVec[0].m());
	hisKLXN[iFile]->Fill(klVec[0].vx());
	hisKLYN[iFile]->Fill(klVec[0].vy());
	hisKLChisqN[iFile]->Fill(klVec[0].chisqZ());
	if( klVec.size() >1 ){
	  hisKLSecChisqN[iFile]->Fill(klVec[1].chisqZ());
	}
	for(pit= klVec[0].pi0().begin();
	    pit!=klVec[0].pi0().end();
	    pit++){
	  hisPi0MassN[iFile]->Fill((*pit).m());       
	}
      }

      Int_t KLEIndex = (int)(klVec[0].e()/500);
      if( KLEIndex <0|| KLEIndex >= 10){ continue;}
      hisKLZAcceptance[iFile][KLEIndex]->Fill(klVec[0].vz());            
    }
  }
  for( int iFile = 0; iFile < nFile; iFile++){
    hisKLE[iFile]->SetLineColor(iFile+1);
    hisKLP[iFile]->SetLineColor(iFile+1);
    hisKLZ[iFile]->SetLineColor(iFile+1);
    hisKLMass[iFile]->SetLineColor(iFile+1);
    hisKLChisq[iFile]->SetLineColor(iFile+1);
    hisKLSecChisq[iFile]->SetLineColor(iFile+1);
    hisGammaE[iFile]->SetLineColor(iFile+1);
  }
 
  for( int i = 0; i < nFile; i++){
    hisKLChisq[i]->Write();
    hisKLSecChisq[i]->Write();
    hisGammaE[i]->Write();
    hisGammaTime[i]->Write();
    hisGammaPosX[i]->Write();
    hisGammaPosY[i]->Write();

    hisPi0E[i]->Write();
    hisPi0P[i]->Write();
    hisPi0Pt[i]->Write();
    hisPi0Z[i]->Write();
    hisPi0ZSig[i]->Write();
    hisPi0Mass[i]->Write();
    
    hisKLE[i]->Write();
    hisKLP[i]->Write();
    hisKLZ[i]->Write();
    hisKLX[i]->Write();
    hisKLY[i]->Write();
    hisKLZPosition[i]->Write();
    hisKLMass[i]->Write();
    hisVETOHistory[i]->Write();

    hisKLEN[i]->Write();
    hisKLMassN[i]->Write();
    hisKLXN[i]->Write();
    hisKLYN[i]->Write();
    hisPi0MassN[i]->Write();
    hisKLChisqN[i]->Write();
    hisKLSecChisqN[i]->Write();
   
  }
  for( int i = 0; i< nFile; i++){
    for( int j = 0; j< 10; j++){
      hisKLZAcceptance[i][j]->Write();
    }
    for( int j = 0; j< 20; j++){
      hisL1TrigCount[i][j]->Write();
    }
  }
  tfOut->Write();
  tfOut->Close();
}
