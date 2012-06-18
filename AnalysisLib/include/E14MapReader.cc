#include "E14MapReader.h"

E14MapReader::E14MapReader( const char* filename ){
  mfile = new TFile( filename );
  mtree = (TTree*)mfile->Get("map");
  mNmodule = 0;
  Init();
}

E14MapReader::~E14MapReader(){

}

void E14MapReader::Init(){
  bAdd = true;
  bSet = false;
}

bool E14MapReader::Add( const char* moduleName ){
  std::string newModulename = moduleName; 
  if( !bAdd ){ return false; }
  if( mNmodule >= nMaxModule || mNmodule < 0 ){bAdd = false; return false;}
  if( bSet ){return false; }
  for( int i = 0; i< mNmodule ;i++){
    if( newModulename.compare( module[i]->GetName() ) == 0 ){ return false; }
  }
  std::cout<<"Make module : " <<  moduleName << " : " << mNmodule << std::endl;
  module[mNmodule] = new E14ModuleMap(moduleName, mtree);
  mNmodule++;
  return true;
}

bool E14MapReader::SetMap(){
  if( bSet ) { return false; } 
  for( int i = 0; i< mNmodule ; i++){
    module[i]->GetEntry(0);
    //std::cout<< module[i]->GetAllNum() << std::endl;
  }
  bSet = true; 
  return bSet;
}

bool E14MapReader::VerifyModuleID(int moduleID){
  if( moduleID >= nMaxModule || moduleID < 0 ){bAdd = false; return false;}
  return true;
}

bool E14MapReader::GetCFC( const char* moduleName, int modulechID, int& crateID, int& fadcID, int& CHID){
  int modid = FindModuleNumber(moduleName);
  if( !VerifyModuleCHID( modid, modulechID ) ){ return false;}
  module[modid]->GetCFC(modulechID, crateID, fadcID, CHID );
  return true;
}

bool E14MapReader::VerifyModuleCHID(int moduleID, int modulechID){  
  if( !VerifyModuleID( moduleID ) ){ return false; }
  if( modulechID > module[moduleID]->GetAllNum() || modulechID < 0){ return false; }
  return true;
}

int E14MapReader::FindModuleNumber( const char* moduleName ){
  std::string name = moduleName;
  for( int i = 0; i< mNmodule; i++){
    if( name.compare( module[i]->GetName() ) == 0 ){ return i;}
  }
  return -1;
}

int E14MapReader::CopyCFCtoMap( int modid , int chmap[4096][3] ){
  if( !VerifyModuleID( modid )){ return 0;}
  for( int i = 0; i< module[modid]->GetAllNum() ; i++){
    module[modid]->GetCFC( i, chmap[i][0], chmap[i][1], chmap[i][2]);
  }
  return module[modid]->GetAllNum();
}
int E14MapReader::CopyCFCtoMap(  const char* moduleName, int chmap[4096][3]){  
  int modid = FindModuleNumber( moduleName );
  if( !VerifyModuleID( modid )){ return 0;}
  for( int i = 0; i< module[modid]->GetAllNum() ; i++){
    module[modid]->GetCFC( i, chmap[i][0], chmap[i][1], chmap[i][2]);
  }
  return module[modid]->GetAllNum();
}

const char* E14MapReader::GetModuleName( int modid ){  
  return (module[modid]->GetName()).c_str();
}

  
