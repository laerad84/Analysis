#include "E14ConvWriter.h"

E14ConvWriter::E14ConvWriter(char* SumupFile,TTree* tr){  
  m_mapFilename = SumupFile;
  m_tr = tr;
  map  = new E14MapReader( m_mapFilename.c_str() );
  m_nModule   = 0; 
  bInitialize = false;
  m_ModID_Cosmic = GetModuleNumber( "Cosmic" );
  m_ModID_CV     = GetModuleNumber( "CV" );
  m_ModID_Laser  = GetModuleNumber( "Laser" );
  m_gr           = new TGraph();
}
E14ConvWriter::~E14ConvWriter(){
  //delete map;
}
bool E14ConvWriter::AddModule( char* ModuleName ){
  std::string ModuleNameStr = ModuleName;
  m_modList.push_back( Modulename );
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
bool E14ConvWriter::ScanMod(char* modName){
  std::list<std::string>iterator it;
  std::string modNameStr  = modName;
  for( it  = m_modList.begin();
       it != m_modList.end();
       it++){
    int rst = modNameStr.compare( *it );
    if( rst == 0 ){
      return true;
    }
  }
  return false;
}
bool E14ConvWriter::TrigJudgement(){
}
bool E14ConvWriter::InitTriggerFlag(){
  m_TrigFlag   = 0;
  m_CosmicTrig = 0;
  m_LaserTrig  = 0; 
  m_CVTrig     = 0;

  m_CosmicTrigFlagUp = 0;
  m_CosmicTrigFlagDn = 0;
  
  m_CVTrigFlag = 0;
  m_LaserFlag  = 0;

  return true; 
}
int  E14ConvWriter::GetNmodule(){
  return m_nModule;
}
int  E14ConvWriter::GetModuleID( char* ){
  std::list<std::string>iterator it;
  std::string modNameStr  = modName;
  int ModID = 0; 
  for( it  = m_modList.begin();
       it != m_modList.end();
       it++){
    int rst = modNameStr.compare( *it );
    if( rst == 0 ){
      return ModID;
    }
    ModID++;
  }
  return -1;
}
int  E14ConvWriter::Fit( int ModID){
  for( iSubMod = 0; iSubMod < (this->ModMap[iMod]).nMod; iSubMod++){
    int iCrate = 9999;
    int iSlot  = 9999;
    int iCh    = 9999;
    iCrate = (this->ModMap[iMod]).Map[iSubMod][0];
    iSlot  = (this->ModMap[iMod]).Map[iSubMod][1];
    iCh    = (this->ModMap[iMod]).Map[iSubMod][2];
    if( iCrate == 9999 || iSlot == 9999 || iCh == 9999){ continue; }
    m_gr->Set(0);
  }
  return 0;
}
int  E14ConvWriter::FitAll(){
}
