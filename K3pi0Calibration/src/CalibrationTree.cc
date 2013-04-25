#include "CalibrationTree.h"

CalibrationTree::CalibrationTree(){
  m_Tree = NULL;
  InitValue();
}

CalibrationTree::~CalibrationTree(){
}

void CalibrationTree::Branch(TTree* tr){
  m_Tree = tr;
  m_Tree->Branch("FlagKL_prefit" ,&FlagKL_prefit,"FlagKL_prefit/I");
  m_Tree->Branch("FlagKL"        ,FlagKL        ,"FlagKL[6]/I");
  m_Tree->Branch("FlagCluster"   ,FlagCluster   ,"FlagCluster[6]/I");
  m_Tree->Branch("FlagCalibrated",FlagCalibrated,"FlagCalibrated[6]/I");
  m_Tree->Branch("CorrID"        ,CorrID        ,"CorrID[6]/I");
  m_Tree->Branch("Corr"          ,Corr          ,"Corr[6]/D");
  m_Tree->Branch("CorrE"         ,CorrE         ,"CorrE[6]/D");
  m_Tree->Branch("Ratio"         ,Ratio         ,"Ratio[6]/D");
  m_Tree->Branch("SecondRatio"   ,SecondRatio   ,"SecondRatio[6]/D");
  m_Tree->Branch("GammaEnergy"   ,GammaEnergy   ,"GammaEnergy[6]/D");
  m_Tree->Branch("GammaSigma"    ,GammaSigma    ,"GammaSigma[6]/D");
  m_Tree->Branch("chisq"         ,chisq         ,"chisq[6]/D");
  m_Tree->Branch("LeadingChID"   ,LeadingChID   ,"LeadingChID[6]/I");
  m_Tree->Branch("LeadingEnergy" ,LeadingEnergy ,"LeadingEnergy[6]/D");
  m_Tree->Branch("LeadingHeight" ,LeadingHeight ,"LeadingHeight[6]/D");
  m_Tree->Branch("nCalibrated"   ,&nCalibrated  ,"nCalibrated/I");
}
void CalibrationTree::SetBranchAddress(TTree* tr){
  m_Tree = tr;
  m_Tree->SetBranchAddress("FlagKL_prefit" ,&FlagKL_prefit);
  m_Tree->SetBranchAddress("FlagKL"        ,FlagKL);
  m_Tree->SetBranchAddress("FlagCluster"   ,FlagCluster);
  m_Tree->SetBranchAddress("FlagCalibrated",FlagCalibrated);
  m_Tree->SetBranchAddress("CorrID"        ,CorrID);
  m_Tree->SetBranchAddress("Corr"          ,Corr);
  m_Tree->SetBranchAddress("CorrE"         ,CorrE);
  m_Tree->SetBranchAddress("Ratio"         ,Ratio);
  m_Tree->SetBranchAddress("SecondRatio"   ,SecondRatio);
  m_Tree->SetBranchAddress("GammaEnergy"   ,GammaEnergy);
  m_Tree->SetBranchAddress("GammaSigma"    ,GammaSigma);
  m_Tree->SetBranchAddress("chisq"         ,chisq);
  m_Tree->SetBranchAddress("LeadingChID"   ,LeadingChID);
  m_Tree->SetBranchAddress("LeadingEnergy" ,LeadingEnergy);
  m_Tree->SetBranchAddress("LeadingHeight" ,LeadingHeight);
  m_Tree->SetBranchAddress("nCalibrated"   ,&nCalibrated);
}

int CalibrationTree::GetEntries(){
  int rst = m_Tree->GetEntries();
  return rst;
}
int CalibrationTree::GetEntry( int ientry ){
  int rst = m_Tree->GetEntry( ientry );
  return rst;
}

int CalibrationTree::InitValue(){
  FlagKL_prefit = 0;
  for (int i = 0; i< 6; ++i){
    FlagKL[i]         = 0;
    FlagCluster[i]    = 0;
    FlagCalibrated[i] = -1;
    CorrID[i] = -1;
    Corr[i]   = 0.;
    CorrE[i]  = 0.;
    Ratio[i]  = 0.;
    SecondRatio[i] = 0.; 
    GammaEnergy[i] = 0.;
    GammaSigma[i]  = 0.;
    chisq[i]       = 0.;
    LeadingChID[i] = -1;
    LeadingEnergy[i] = 0;
    LeadingHeight[i] = 0;
  }
  nCalibrated= 0;
  return 0;
}


