#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <list>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TROOT.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "rec2g/Rec2g.h"
#include "cluster/Cluster.h"

#include "gamma/Gamma.h"
#include "pi0/Pi0.h"
#include "klong/Klong.h"
#include "klong/RecKlong.h"
#include "e14_Calibration/KL_calibration.h"

#include "TApplication.h"
#include "TCanvas.h"

class Gamma;
class Klong;

/*
static int        csi_cal_factor_count[N_CSI];
static double     csi_cal_factor[N_CSI];
static double     csi_cal_factor_sum[N_CSI];
static double     csi_cal_factor_weight[N_CSI];
static double     csi_energy[N_CSI];
static double     csi_time[N_CSI];
static double     csi_dxy[N_CSI];
static double     csi_energy[N_CSI];
static double     csi_time[N_CSI];
static CLHEP::Hep3Vector csi_pos[N_CSI];
*/

//double Refit_KL_full(Klong &KL);
//double Refit_KL_freelast(Klong &KL);
//int    CalEnergy_idv(Klong &KL_prefit);
//void   user_ana_recg6_init();

int
main(int argc,char** argv){
  if( argc != 3 ){
    std::cout << "::Usage::\n"
	      << Form("bin/%s inputFilename outputFilename\n",argv[0])
	      << std::endl;
    return 0; 
  }
  TH1D* his_Fit[100];
  TH1D* his_CSI[N_CSI];
  user_ana_recg6_init(his_Fit,his_CSI);
  std::cout << his_CSI[0] << std::endl;
  
  
  std::string ifilename = argv[1];
  std::string ofilename = argv[2];
  E14GNAnaDataContainer data;
  // prepare for Read out
  TChain* InputTree = new TChain("Tree");
  //InputTree->Add(ifilename.c_str());
  std::string dir = "/home/had/jwlee/workdir/RootFiles/out_KL3pi0.mac_500000_%d_FEB_CL_KL.root";
  
  for( int seedIndex = 0; seedIndex < 99; seedIndex++){
    InputTree->Add(Form(dir.c_str(),seedIndex));
  }

  std::cout << InputTree->GetEntries() << std::endl;
  std::cout << "Input File name:" << ifilename << std::endl;
  data.setBranchAddress( InputTree );
  
  // prepare for Write out
  TFile* OutputFile = new TFile(ofilename.c_str(),"recreate");
  TTree* OutputTree = new TTree("tr_Calibration","Tree for Calibration");
  std::cout << " Output File name:" << ofilename << std::endl;
  data.branchOfKlong( OutputTree );
  std::vector<Klong> KLvec;
  
  int nKL = 0;
  int nloop = InputTree->GetEntries();
  std::cout << nloop << std::endl;
  for( int ievent = 0; ievent < nloop ; ievent++){
    if( ievent %( nloop/10 ) == 0 && nloop >100 ){
      std::cout << ievent/ (nloop/100) << "%" << std::endl;
    }
    InputTree->GetEntry( ievent );
    data.getData( KLvec );
    
    /////////////////////
    // iteration loop //
    int result = CalEnergy_idv(KLvec[0],his_Fit,his_CSI );
    if( result ==1 ) nKL++;
    
    OutputTree->Fill();    
    //std::cout << "Fill " << std::endl;
  }
  std::cout << "nKL:"<<  nKL <<std::endl;
  for( int i = 0; i < N_CSI; i++){
    if( his_CSI[i]->Integral() > 0)
      std::cout <<"Write()" << i <<  " : "<< his_CSI[i]->Integral()<<  std::endl;
    his_CSI[i]->Write();
  }
  for( int i = 0; i< 5; i++){
    his_Fit[i]->Write();
  }
  
  
  OutputTree->Write();
  OutputFile->Close();
  
}

