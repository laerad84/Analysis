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

  const int nFile = 6;
  TFile* tf[nFile]; 
  TTree* tr[nFile];  
  char* name[nFile] = {"SIMFAST","3pi0_LaserComp","WAV","SIM","3pi0_OldComp","WAVNOCV"};
  
  for( int i = 0; i < nFile; i++){
    tf[i] = new TFile(Form("Kl_Total_%s.root",name[i]));
    tr[i] = (TTree*)tf[i]->Get(Form("trKL"));
  }
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];

  /*
  Double_t KLMass;
  Double_t KLChisq;
  Double_t KLSecChisq;
  Double_t KLE;
  Double_t KLPos[nFile];
  Double_t KLMom[nFile];
  Double_t GammaE[6];
  Double_t GammaPos[6][nFile]; 
  Double_t GammaTime[6];

  for( int i = 0; i< nFile; i++){
    tr[i]->SetBranchAddress("KLMass",&KLMass);
    tr[i]->SetBranchAddress("KLChisq",&KLChisq);
    tr[i]->SetBranchAddress("KLSecChisq",&KLSecChisq);
    tr[i]->SetBranchAddress("KLE",&KLE);
    tr[i]->SetBranchAddress("KLPos",KLPos);
    tr[i]->SetBranchAddress("KLMom",KLMom);
    tr[i]->SetBranchAddress("GammaE",GammaE);
    tr[i]->SetBranchAddress("GammaPos",GammaPos);
    tr[i]->SetBranchAddress("GammaTime",GammaTime);
  }
*/
  TH1D* hisKLP[nFile];
  TH1D* hisKLZ[nFile];
  TH1D* hisKLE[nFile];
  TH1D* hisKLMass[nFile];
  TH1D* hisKLChisq[nFile];
  TH1D* hisKLSecChisq[nFile];
  TH1D* hisGammaE[nFile];
  TH1D* hisKLZAcceptance[nFile][10];
  TH1D* hisKLZPosition[nFile];

  for( int i = 0; i< nFile; i++){
    for( int j = 0; j< 10; j++){
      hisKLZAcceptance[i][j] = new TH1D(Form("hisKLZAcceptance_%d_%d",i,j),Form("hisKLZAcceptance_%d_%d_%d",i,500*j,500*(j+1)),70,0,7000);
    }
    hisKLZPosition[i]= new TH1D(Form("hisKLZPosition_%d",i),Form("hisKLZPosition_%s",name[i]),70,0,7000);
    hisKLZ[i]        = new TH1D(Form("hisKLZ_%d",i),Form("hisKLZ_%s",name[i]),70,0,7000);
    hisKLP[i]        = new TH1D(Form("hisKLP_%d",i),Form("hisKLP_%s",name[i]),100,0,10000);
    hisKLE[i]        = new TH1D(Form("hisKLE_%d",i),Form("hisKLE_%s",name[i]),100,0,10000);
    hisKLMass[i]     = new TH1D(Form("hisKLMass_%d",i),Form("hisKLMass_%s",name[i]),50,400,800);
    hisKLChisq[i]    = new TH1D(Form("hisKLChisq_%d",i),Form("hisKLChisq_%s",name[i]),100,0,100);
    hisKLSecChisq[i] = new TH1D(Form("hisKLSecChisq_%d",i),Form("hisKLSecChisq_%s",name[i]),100,0,100);
    hisGammaE[i]     = new TH1D(Form("hisGammaE_%d",i),Form("hisGammaE_%s",name[i]),100,0,2000);
  }

  TFile* tfOut = new TFile("DistributionSIMDATA.root","recreate");
  for( int iFile = 0; iFile < nFile; iFile++){
    E14GNAnaDataContainer data;
    data.setBranchAddress( tr[iFile] );
    tr[iFile]->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
    tr[iFile]->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);

    for( int ievent = 0; ievent < tr[iFile]->GetEntries(); ievent++){      
      tr[iFile]->GetEntry(ievent);
      //if( ievent  >= 100000 ){ break ; } 
      std::list<Cluster> clist;
      std::list<Gamma>   glist;
      std::vector<Klong>   klVec;
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
      if( klVec[0].chisqZ() > 6 ){continue;} 
      bool  bInnerGamma = false;
      bool  bOuterGamma = false;
      bool  bClusterMaxE= false;
      // Cut on Gamma // 
      Double_t MinGammaE = 200;
      Int_t nEGamma = 0;      
      if( klVec.size() > 1 ){
	if( klVec[1].chisqZ() - klVec[0].chisqZ() < 6 ){ continue; }
      }


      Double_t ClusterMaxE = 0;
      std::list<Gamma>::iterator git = glist.begin();
      for( int igamma = 0; igamma < 6; igamma++,git++){

	hisGammaE[iFile]->Fill((*git).e());
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
	if( TMath::Abs((*git).y())>575){ bOuterGamma = true; }
      }


      if( bClusterMaxE){ continue; }
      if( nEGamma < 6 ){ continue; }
      if( bInnerGamma ){ continue; }
      if( bOuterGamma ){ continue; }
      //if( bInnerGamma || bOuterGamma ) { continue; }

      /*
      // CutValues // 
      // true->Cut // 
      bool  bNGamma     = false;
      bool  bRWGamma    = false;
      bool  bTGamma     = false;
      bool  bGammaX     = false;
      bool  bGammaGood  = false;

      Int_t nGamma = 0;
      Int_t nGammaGood = 0;
      Int_t GEIndex = -1;
      Double_t MaximumR = 0; 
      Double_t GammaMinDist = 1000;
      Double_t GammaTMean   = 0; 
      Double_t MaximumTDelta= 0;
      Double_t MinGammaE = 1000000;
      
      // Cut on Gamma // 
      for( int igamma = 0; igamma < 6; igamma++){
	if( GammaE[igamma] < MinGammaE ){MinGammaE = GammaE[igamma];}
	if( TMath::Abs(GammaPos[igamma][0]) < 150 &&
	    TMath::Abs(GammaPos[igamma][1]) < 150 ){
	  bInnerGamma = true; 	  
	}
	if( TMath::Sqrt(GammaPos[igamma][0]*GammaPos[igamma][0]+GammaPos[igamma][1]*GammaPos[igamma][1]) > 800){
	  bOuterGamma = true;
	}	
	if( TMath::Abs(GammaPos[igamma][1])>550){ bOuterGamma = true; }
	if( GammaE[igamma] > 200 ){
	  nGamma++;
	}
	if( GammaPos[igamma][1]>0&&GammaPos[igamma][0] > -200 ){ nGammaGood++;}
	else if( GammaPos[igamma][1]<=0&&GammaPos[igamma][0] > 175 ){ nGammaGood++;}
	
	if( TMath::Abs( GammaPos[igamma][0] ) > 600 ){ bGammaX = true ;}
	TVector2 vec( GammaPos[igamma][0], GammaPos[igamma][1]);
	if( vec.Mod() > MaximumR ){ MaximumR = vec.Mod();}
	
	for( int jgamma = igamma+1; jgamma < 6; jgamma++){
	  Double_t R = TMath::Sqrt((GammaPos[igamma][0]-GammaPos[jgamma][0])
				   *(GammaPos[igamma][0]-GammaPos[jgamma][0])
				   +(GammaPos[igamma][1]-GammaPos[jgamma][1])
				   *(GammaPos[igamma][1]-GammaPos[jgamma][1]));
	  if( GammaMinDist > R ) { GammaMinDist = R;}
	}
      }
      //if( nGamma !=6 ){ continue; }

      GEIndex =  (int)(MinGammaE/50);
      if( GammaMinDist < 200 ){bRWGamma = true;}
      //if( nGamma != 6 ){ bNGamma = true; }
      //if( nGammaGood <4 ){bNGamma = true; }
      //if( bGammaX ){ continue; }
      // Cut on KL // 
      bool bklChisq = false;
      bool bkle     = false;
      bool bklz     = false;
      bool bklmass  = false;
      bool bklpt    = false;
      Double_t Klmass = 497.648;
      if( KLE > 5000 ) bkle = true; 
      if( KLMass   < Klmass-10  || KLMass > Klmass+10  ){ bklmass = true; }
      if( KLPos[2] < 1000 || KLPos[2]> 5500){ bklz    = true; } 
      if( klptSq   > 100  ){ bklpt = true; }
      if( KLChisq  > 5 ){ bklChisq = true; }
      
      //if( KLSecChisq - KLChisq < 5 ){ continue; }
      // Initial Cut // 
            
      //if( bInnerGamma || bOuterGamma || bNGamma ){ continue; }      
      //if( bInnerGamma || bOuterGamma || bNGamma || bRWGamma ){ continue; }      
      //if( bklmass || bklpt || bklChisq ){ continue; }
      //if( bTGamma ){ continue; }
      if( bInnerGamma || bOuterGamma ){ continue; }
      std::cout<< GEIndex << std::endl;
      Double_t KLECutMin = 2700;
      Double_t KLECutMax = 3300;
      */
      if( klVec[0].vz() < 5000 && klVec[0].vz() > 3000){
	if(TMath::Abs(klVec[0].m()-KLMass )< 10 ){
	  hisKLP[iFile]->Fill(klMom);
	}
      }
      hisKLE[iFile]->Fill(klVec[0].e());
      hisKLZ[iFile]->Fill(klVec[0].vz());
      /*
      std::cout <<klMom  << "\tSolt:\t:" 
		<< soltFunc->Eval(klMom) << "\tSugar:\t" 
		<< sugarFunc->Eval(klMom)<< "\t"
		<< Ratio << "\t" 
		<< TMath::IsNaN(Ratio) << std::endl;
      */
      if( TMath::IsNaN(Ratio) == 0 ){ 
	hisKLZPosition[iFile]->Fill(klVec[0].vz(),Ratio);
      }
      hisKLMass[iFile]->Fill(klVec[0].m());
      hisKLChisq[iFile]->Fill(klVec[0].chisqZ());
      if( klVec.size() >1 ){
	hisKLSecChisq[iFile]->Fill(klVec[1].chisqZ());
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

  std::cout << "Draw" << std::endl;
  TCanvas* can = new TCanvas("can","can",1200,800);
  can->Divide(3,2);
  Double_t ScaleFactor = hisKLZ[0]->Integral()/hisKLZ[1]->Integral();
  can->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  hisKLE[0]->Draw();
  hisKLE[1]->Scale(ScaleFactor);
  hisKLE[1]->Draw("same");
  can->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  hisKLZ[0]->Draw();
  hisKLZ[1]->Scale(ScaleFactor);
  hisKLZ[1]->Draw("same");
  /*
  hisKLZ[2]->SetLineColor(3);
  hisKLZ[2]->Scale(ScaleFactor);
  hisKLZ[2]->Draw("same");
  */
  can->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  hisKLMass[0]->Draw();
  hisKLMass[1]->Draw("same");
  can->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  hisKLChisq[0]->Draw();
  hisKLChisq[1]->Scale(ScaleFactor);
  hisKLChisq[2]->Scale(ScaleFactor);
  hisKLChisq[1]->Draw("same");
  hisKLSecChisq[0]->Draw("same");
  hisKLSecChisq[1]->Draw("same");
  can->cd(5);
  gPad->SetGridx();
  gPad->SetGridy();
  gPad->SetLogy();
  hisGammaE[0]->Draw();
  hisGammaE[1]->Scale(ScaleFactor);
  hisGammaE[1]->Draw("same");
  /*
  hisKLE[0]->Write();
  hisKLE[1]->Write();
  hisKLP[0]->Write();
  hisKLP[1]->Write();
  hisKLZ[0]->Write();
  hisKLZ[1]->Write();
  hisKLZPosition[0]->Write();
  hisKLZPosition[1]->Write();
  hisKLMass[0]->Write();
  hisKLMass[1]->Write();
  */
  for( int i = 0; i < nFile; i++){
    hisKLE[i]->Write();
    hisKLP[i]->Write();
    hisKLZ[i]->Write();
    hisKLZPosition[i]->Write();
    hisKLMass[i]->Write();
  }
  for( int i = 0; i< nFile-1; i++){
    for( int j = 0; j< 10; j++){
      hisKLZAcceptance[i][j]->Write();
    }
  }
  tfOut->Write();
  tfOut->Close();
}
