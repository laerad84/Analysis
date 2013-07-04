#include "Calibration.h"
#include "CalibrationTree.h"
#include "ConvReader.h"
#include <string>
#include <iostream>
#include <fstream>

#include "E14Fsim/E14FsimFunction.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "cluster/Cluster.h"
#include "cluster/ClusterFinder.h"
#include "gamma/Gamma.h"
#include "gamma/GammaFinder.h"
#include "rec2g/Rec2g.h"

#include "User_Function.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "TApplication.h"
#include "TCanvas.h"
int main(int argc, char** argv){

  TApplication* app = new TApplication("app",&argc, argv);
  std::string filename = "/Volume0/Simulation/Conv_KL3pi0_FAST_REDUCED_5E6_0.root";
  TFile* tf = new TFile(filename.c_str());
  TTree* tr = (TTree*)tf->Get("T");
  ConvReader* reader = new ConvReader(tr);  
  Int_t RearrangeID[6] = {0,1,2,3,4,5};
  int nCSIDigi;
  Double_t CsiEne[2716];
  Int_t CsiID[2716];
  Double_t CsiTime[2716];

  Calibration* calibrator = new Calibration();
  for( int i = 0; i< 6; i++){
    //calibrator->grChisq[i]->SetMarkerStyle(21);
    calibrator->grChisq[i]->SetMarkerSize(1);
  }
  CalibrationTree calData; 
  ClusterFinder clusterFinder;
  GammaFinder   gFinder;

  TGraph* grChisquareTest[6];
  TGraph* grCorrTest[6];
  for( int i = 0; i< 6; i++){
    grChisquareTest[i] = new TGraph();
    grChisquareTest[i]->SetNameTitle(Form("ChisquareTest_%d",i),Form("ChisquareTest_%d;Energy[Mev];#Chi^{2}",i));
    grCorrTest[i] = new TGraph();
    grCorrTest[i]->SetMarkerColor(2);
  }

  TCanvas* can = new TCanvas("can","",1200,800);
  can->Divide(3,2);
  for( int ievt = 0; ievt < 2000; ievt++){
    tr->GetEntry(ievt);
    nCSIDigi = 0;
    for( int i = 0; i< 6; i++){
	grChisquareTest[i]->Set(0);
	grCorrTest[i]->Set(0);
      }
    for( int idigi = 0; idigi < reader->CsiNumber; idigi++){
      if( reader->CsiEne[idigi] > 3 ){
	CsiEne[nCSIDigi]  = reader->CsiEne[idigi];
	CsiTime[nCSIDigi] = reader->CsiTime[idigi];
	CsiID[nCSIDigi]   = reader->CsiModID[idigi];
	nCSIDigi++;
      }
    }

    if( reader->nTrack != 6 ){ continue; }
    Int_t nGamma = 0;
    for( int itrack = 0; itrack < reader->nTrack; itrack++){
      if( reader->pid[itrack] == 22  ){
	nGamma++;
      }
    }
    if( nGamma != 6 ){ 
      continue; 
    }
    
    std::list<Cluster> clist; 
    std::list<Gamma>   glist;    
    std::vector<Klong> klVec;
    
    clist = clusterFinder.findCluster(nCSIDigi, CsiID, CsiEne, CsiTime );
    gFinder.findGamma( clist, glist );
    if( glist.size() != 6 ){ continue; }
    Int_t grIndex  = 0;
    std::list<Gamma>::iterator it = glist.begin();
    for( ; it != glist.end(); it++){
      double InitialEnergy = (*it).e();
      double InitialSigma  = (*it).sigmaE();
      for( int iCalIndex  =0; iCalIndex < 21; iCalIndex++){
	double TestEnergy = InitialEnergy*(1-0.1*(1-0.1*iCalIndex));
	(*it).setEnergy(TestEnergy);
	double TestSigma  = E14FsimFunction::getFunction()->csiEnergyRes(TestEnergy/1000.)*1000;
	(*it).setSigmaE(TestSigma);

	if( user_rec( glist,klVec)){
	  bool GammaFlag = false; 
	  for( unsigned int i = 0; i< klVec[0].pi0().size(); i++){
	    Double_t x,y;
	    x = klVec[0].pi0()[i].g1().x();
	    y = klVec[0].pi0()[i].g1().y();
	    if( x*x +y*y > 850*850 ){ GammaFlag = true; }
	    if( TMath::Abs(x) < 150 && TMath::Abs(y) < 150 ){ GammaFlag = true;}
	    x = klVec[0].pi0()[i].g2().x();
	    y = klVec[0].pi0()[i].g2().y();
	    if( x*x +y*y > 850*850 ){ GammaFlag = true; }
	    if( TMath::Abs(x) < 150 && TMath::Abs(y) < 150 ){ GammaFlag = true;}
	  }
	  if( GammaFlag ){continue;}
	  //std::cout<<"EventNumber" <<  ievt << std::endl;
	  int result = calibrator->CalEnergy_idv(klVec);
	  std::cout<< calibrator->GetChisq(grIndex)  << "\t" <<  calibrator->GetCorr(grIndex) << std::endl;
	  grChisquareTest[grIndex]->SetPoint( grChisquareTest[grIndex]->GetN(),
					      (1-0.1*(1-0.1*iCalIndex)),
					      calibrator->GetFirstChisq(grIndex));
	  grCorrTest[grIndex]->SetPoint( grCorrTest[grIndex]->GetN(),
					 calibrator->GetCorr(grIndex),
					 calibrator->GetChisq(grIndex));
					 
	  /*
	    if( result >= 1 ){
	    for( int iGamma = 0; iGamma < 6; iGamma++){
	    can->cd(iGamma+1);
	    calibrator->grChisq[iGamma]->Draw("APL");
	    std::cout<< iGamma << ":" 
	    << calibrator->GetCorr(RearrangeID[iGamma]) << "\t"
	    << calibrator->GetCorrE(RearrangeID[iGamma])<< "\t" 
	    << calibrator->GetChisq(RearrangeID[iGamma])<< std::endl;
	  }
	  can->Update();
	  can->Modified();
	  can->Update();
	  getchar();
	  }
	  */
	}
      }
      grIndex++;
      (*it).setEnergy(InitialEnergy);
      (*it).setSigmaE(InitialSigma);
    }
    std::cout << calibrator->GetNCalibrated() << std::endl;
    for( int iGamma = 0; iGamma < 6; iGamma++){
    grChisquareTest[iGamma]->GetYaxis()->SetRangeUser(0,4);
    can->cd(iGamma+1);
    grChisquareTest[iGamma]->Draw("AP");
    grCorrTest[iGamma]->Draw("P");
    }
    can->Update();
    can->Modified();
    can->Update();
    getchar();

  }
   


  app->Run();
  return 0; 
}
