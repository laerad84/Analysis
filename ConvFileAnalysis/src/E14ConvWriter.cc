#include "ConvFileAnalysis/E14ConvWriter.h"

E14ConvWriter::E14ConvWriter(char* SumupFile,TTree* tr){  
  m_mapFilename = SumupFile;
  m_tr = tr;
  map  = new E14MapReader( m_mapFilename.c_str() );
  m_nModule   = 0; 
  bInitialize = false;
}

E14ConvWriter::~E14ConvWriter(){
  //delete map;
}

bool E14ConvWriter::AddModule( char* ModuleName ){

  if( bInitialize ){ return false; }
  if( m_nModule > nMaxModule ){return false; }

  mod[m_nModule] = new E14ConvWriterModule(m_tr, ModuleName);
  map->Add(ModuleName);
  m_nModule++;
  return true;

}

bool E14ConvWriter::Set(){
  if( bInitialize ){
    return false;
  }else{
    bInitialize = true;
    return true;
  }
}

bool E14ConvWriter::SetMap(){
  map->SetMap();
  for( int iMod =0 ; iMod< map->GetNmodule(); iMod++){
    map->CopyMap( iMod, ModMap[iMod]);
  }
  std::cout << "Setup Map is Done" << std::endl;
  std::cout << map->GetNmodule()   << std::endl;			 
  return true;
}

bool E14ConvWriter::Branch(){
  for( int i = 0; i< m_nModule; i++){
    mod[i]->Branch();
  }
  return true;
}

bool E14ConvWriter::SetBranchAddress(){
  for( int i = 0; i< m_nModule; i++){
    mod[i]->SetBranchAddress();
  }
  return true;
}

bool E14ConvWriter::InitData(){
  for( int i = 0; i < m_nModule; i++){
    mod[i]->InitData();
  }
  return true;
}

int E14ConvWriter::GetNmodule() {
  return m_nModule;
}
