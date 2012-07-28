#include "E14ConvWriter.h"


E14ConvWriter::E14ConvWriter( char* mapFilename, TTree* tr){  
  std::string SUMUPFILEDIR = std::getenv("SUMUPFILEDIR");
  std::string CONVFILEDIR  = std::getenv("ROOTFILE_CONV");
  //m_mapFilename = Form("%s/Sum%d.root",SUMUPFILEDIR,RunNumber); 
  m_mapFilename = mapFilename;
  m_tr        = tr;
  map         = new E14MapReader( m_mapFilename.c_str() );
  m_nModule   = 0; 
  bInitialize = false;
  m_gr        = new TGraph();
 
}

E14ConvWriter::E14ConvWriter( int RunNumber, TTree* tr){  
  std::string SUMUPFILEDIR = std::getenv("SUMUPFILEDIR");
  std::string CONVFILEDIR  = std::getenv("ROOTFILE_CONV");
  m_mapFilename = Form("%s/Sum%d.root",SUMUPFILEDIR.c_str(),RunNumber); 
  //m_mapFilename = SumupFile;
  m_tr        = tr;
  map         = new E14MapReader( m_mapFilename.c_str() );
  m_nModule   = 0; 
  bInitialize = false;
  m_gr        = new TGraph();
 
}

E14ConvWriter::~E14ConvWriter(){
  //delete map;
}

bool E14ConvWriter::AddModule( char* ModuleName ){
  std::string ModuleNameStr = ModuleName;
  m_modList.push_back( ModuleNameStr );
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

  std::string modNameStr  = modName;
  for( std::list<std::string>::iterator it = m_modList.begin();
       it != m_modList.end();
       it++){
    int rst = modNameStr.compare( *it );
    if( rst == 0 ){
      return true;
    }
  }
  return false;
}


int  E14ConvWriter::GetNmodule(){
  return m_nModule;
}
int  E14ConvWriter::GetModuleID( char* modName){
  std::list<std::string>::iterator it;
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

