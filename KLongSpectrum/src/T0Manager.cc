#include "T0Manager.h"

T0Manager::T0Manager(){
  Reset();
}

T0Manager::~T0Manager(){
  Reset();
  std::cout<< "T0 Deleted" << std::endl;
}

bool T0Manager::ReadFile(char* filename){
  std::ifstream tmpifs( filename );
  if( !tmpifs.is_open() ){ std::cout<< Form("%s is not exist", filename ) << std::endl;return false;}
  int tmpID;
  double tmpOffset;
  while(tmpifs >> tmpID >> tmpOffset){
    std::cout<< tmpID <<"\t" <<  tmpOffset << std::endl;
    m_T0[tmpID] = tmpOffset;
  }
  Normalizer();
  tmpifs.close();
  return true;
}
double T0Manager::GetT0Offset( int id ){
  return m_T0[id];
}

void T0Manager::Normalizer(){
  int nT0;
  double sumT0;
  for( int iloop = 0; iloop < 10; iloop++){
    for( int i = 0; i< 2716; i++){
      sumT0 += m_T0[i];
      nT0++;
    }
    sumT0/=nT0;
    for( int i = 0; i< 2716; i++){
      m_T0[i]= m_T0[i] - sumT0;
    }
  }
}
void T0Manager::Reset(){
  for( int i = 0; i< 2716; i++){
    m_T0[i] = 0;
  }
}
    
