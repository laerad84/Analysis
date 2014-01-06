#include "L1TrigCounter.h"
L1TrigCounter::L1TrigCounter(){
  hisTrigCounter = new TH1D("hisTrigCounter","hisTrigCounter",nMaxCrate,0,nMaxCrate);
  CsiL1Map    = new CsIPoly("CsiL1Map","CsiL1Map");
  m_Threshold = 10000;
  for( int i = 0; i< 2716; i++){
    m_L1Map[i] = -1;
  }
}

L1TrigCounter::~L1TrigCounter(){
  CsiL1Map->Delete();
  hisTrigCounter->Delete();
}

bool L1TrigCounter::ReadMapFile( char* filename ){ 
  std::cout<< filename << std::endl;
  std::string ANALYSISLIB = std::getenv("ANALYSISLIB");
  std::ifstream ifs(Form("%s/Data/%s",ANALYSISLIB.c_str(),filename));
  int tmpid;
  int tmpcrate;
  int tmpslot;
  int tmpchannel;
  int tmpL1;
  if(!ifs.is_open()){ return false;}
  while( ifs >> tmpid >> tmpcrate >> tmpslot >> tmpchannel >> tmpL1 ){
    if( tmpL1 < 0 || tmpL1 > nMaxCrate ){
      m_L1Map[tmpid] =-1;
    }else{
      m_L1Map[tmpid] = tmpL1;
      CsiL1Map->Fill(tmpid,(double)tmpL1);
    }
  }
  return true;
}

int L1TrigCounter::Fill( int id, double value){
  if( id <0 || id >= 2716 ){ return -1; }
  if( m_L1Map[id] < 0 || m_L1Map[id] >= nMaxCrate ){ return -1; }
  hisTrigCounter->Fill( m_L1Map[id], value );
  m_Et += value;

  return m_L1Map[id];
}

std::vector<int> L1TrigCounter::GetTriggedCrate(){
  std::vector<int> vec;
  for( int i = 0; i< nMaxCrate; i++){
    if( hisTrigCounter->GetBinContent(i+1) > m_Threshold ){
      vec.push_back(i);
    }
  }
  return vec;
}
std::vector<double> L1TrigCounter::GetCount(){
  std::vector<double> vec;
  for( int i = 0; i< nMaxCrate; i++){
    vec.push_back( hisTrigCounter->GetBinContent(i+1));
  }
  return vec;
}

double L1TrigCounter::GetEtCount(){
  return m_Et;
}

void L1TrigCounter::DrawTrig(){
  hisTrigCounter->Draw();
}
void L1TrigCounter::DrawTriggerMap(){
  CsiL1Map->Draw("colz");
}
