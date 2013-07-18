#include <iostream>
#include <fstream>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <list>

#include "gnana/DigiReader.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "klong/Klong.h"
#include "cluster/ClusterFinder.h"
#include "rec2g/Rec2g.h"
#include "gamma/GammaFinder.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "TFile.h"
#include "TTree.h"

#include "TApplication.h"
#include "TSystem.h"
#include "CsIPoly.h"

int main( int argc, char** argv ){
  
  std::string ROOTFILE_SUMUP = std::getenv("ROOTFILE_SUMUP");
  std::string ROOTFILE_WAV   = std::getenv("ROOTFILE_WAV");
  const int nFiles = 3;
  char*  filename[nFiles] = {"SIM","SIMADJ","WAV"};
  std::cout<< __PRETTY_FUNCTION__ << std::endl;

  TFile* tfout = new TFile("ResultOutput.root","recreate");
  TH2D*  hisClusterSizeEne[nFiles][2716];
  TH2D*  hisGChisqEne[nFiles][2716];
  TH2D*  hisClusterSizeZ[nFiles][2716];

  for( int i = 0; i< nFiles; i++){
    for( int j = 0; j< 2716; j++){
      hisClusterSizeEne[i][j] = new TH2D(Form("his_%s_%d",filename[i],j),
					 Form("his_%s_%d",filename[i],j),
					 40,0,2000,40,0,40);
      hisGChisqEne[i][j] = new TH2D(Form("his_chisq_%s_%d",filename[i],j),
				    Form("his_chisq_%s_%d",filename[i],j),
				    40,0,2000,40,0,40);
      hisClusterSizeZ[i][j] = new TH2D(Form("his_z_%s_%d",filename[i],j),
				       Form("his_z_%s_%d",filename[i],j),
				       40,0,6000,40,0,40);
    }
  }

  
  TFile* tfRoot[nFiles]; 
  TTree* tr[nFiles];
  E14GNAnaDataContainer data[nFiles];
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  for( int i = 0; i<nFiles; i++){
    tfRoot[i] = new TFile(Form("Kl_Total_%s.root",filename[i]));
    tr[i] = (TTree*)tfRoot[i]->Get("trKL");
    data[i].setBranchAddress(tr[i]);    
  }
  Int_t CsiL1nTrig;
  Double_t CsiL1TrigCount[20];

  CsIPoly* csi = new CsIPoly("csi","csi");
  std::cout<< __PRETTY_FUNCTION__ << std::endl;
  for( int ifile  =0;  ifile < nFiles; ifile++){
    tr[ifile]->SetBranchAddress("CsiL1TrigCount",CsiL1TrigCount);
    tr[ifile]->SetBranchAddress("CsiL1nTrig",&CsiL1nTrig);
    for( int ievent  = 0; ievent < tr[ifile]->GetEntries(); ievent++){
      //data[ifile].reset();
      tr[ifile]->GetEntry(ievent);
      std::list<Cluster> clist;
      std::list<Gamma>   glist;
      std::vector<Klong> klVec;
      //data.getData( clist );
      //data[ifile].getData( glist );
      data[ifile].getData( klVec );

      
      if( klVec.size() == 0 ){ continue; }
      if( CsiL1nTrig< 5 ){ continue; }

      /*
      if( klVec[0].chisqZ()>6 ){ continue; }
      if( klVec.size() >1 ){
	if( klVec[1].chisqZ() - klVec[0].chisqZ() > 6 ){ continue; }
      }
      */
      Double_t KLZ = klVec[0].vz();
      for( unsigned int ipi = 0; ipi < klVec[0].pi0().size(); ipi++){
	
	Double_t GammaCOEX[nFiles];
	Double_t GammaCOEY[nFiles];
	Double_t GammaDepE[nFiles];
	Double_t GammaChisq[nFiles];
	Int_t    ClusterSize[nFiles];
	Int_t    CenterID[nFiles]; 
	
	GammaCOEX[0]=klVec[0].pi0()[ipi].g1().coex();
	GammaCOEY[0]=klVec[0].pi0()[ipi].g1().coey();
	GammaDepE[0]=klVec[0].pi0()[ipi].g1().edep();
	GammaChisq[0]=klVec[0].pi0()[ipi].g1().chisq();
	ClusterSize[0] = klVec[0].pi0()[ipi].g1().clusterIdVec().size();
	CenterID[0] = csi->Fill( -GammaCOEX[0], GammaCOEY[0], GammaDepE[0])-1;
	if( CenterID[0] < 0 || CenterID[0] >= 2716){ continue; }
	hisClusterSizeEne[ifile][CenterID[0]]->Fill( GammaDepE[0], ClusterSize[0] );
	hisGChisqEne[ifile][CenterID[0]]->Fill( GammaDepE[0], GammaChisq[0] );
	hisClusterSizeZ[ifile][CenterID[0]]->Fill( KLZ,ClusterSize[0] );

	GammaCOEX[1]=klVec[0].pi0()[ipi].g2().coex();
	GammaCOEY[1]=klVec[0].pi0()[ipi].g2().coey();
	GammaDepE[1]=klVec[0].pi0()[ipi].g2().edep();
	GammaChisq[1]=klVec[0].pi0()[ipi].g2().chisq();
	ClusterSize[1] = klVec[0].pi0()[ipi].g2().clusterIdVec().size();
	CenterID[1] = csi->Fill( -GammaCOEX[1], GammaCOEY[1], GammaDepE[1])-1;
			    
	if( CenterID[1] < 0 || CenterID[1] >= 2716){ continue; }
	hisClusterSizeEne[ifile][CenterID[1]]->Fill( GammaDepE[1], ClusterSize[1] );
	hisGChisqEne[ifile][CenterID[1]]->Fill( GammaDepE[1], GammaChisq[1] );
	hisClusterSizeZ[ifile][CenterID[1]]->Fill( KLZ, ClusterSize[1]);

	//std::cout << CenterID[1] << "\t" << GammaDepE[1] << "\t" << ClusterSize[1] << std::endl;
      } 
    }
  }
  tfout->cd();
  for( int j = 0; j< 2716; j++){
    for( int i = 0; i < nFiles; i++){
      hisClusterSizeEne[i][j]->Write();
      hisGChisqEne[i][j]->Write();
      hisClusterSizeZ[i][j]->Write();
    }
  }
  tfout->Close();
}
