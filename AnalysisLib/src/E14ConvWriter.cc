#include "E14ConvWriter.h"

E14ConvWriter::E14ConvWriter( char* mapFilename, TTree* tr){  
  std::string SUMUPFILEDIR = std::getenv("ROOTFILE_SUMUP");
  std::string CONVFILEDIR  = std::getenv("ROOTFILE_CONV");
  //m_mapFilename = Form("%s/Sum%d.root",SUMUPFILEDIR,RunNumber); 
  m_mapFilename = mapFilename;
  m_tr        = tr;
  map         = new E14MapReader( m_mapFilename.c_str() );
  m_nModule   = 0; 
  bInitialize = false;
 
}
E14ConvWriter::E14ConvWriter( int RunNumber, TTree* tr){  
  std::string SUMUPFILEDIR = std::getenv("ROOTFILE_SUMUP");
  std::string CONVFILEDIR  = std::getenv("ROOTFILE_CONV");
  m_mapFilename = Form("%s/Sum%d.root",SUMUPFILEDIR.c_str(),RunNumber); 
  //m_mapFilename = SumupFile;
  m_tr        = tr;
  map         = new E14MapReader( m_mapFilename.c_str() );
  m_nModule   = 0; 
  bInitialize = false;
 
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
  m_tr->Branch("RunNo"           ,&m_RunNo           ,"RunNo/I");
  m_tr->Branch("EventNo"         ,&m_EventNo         ,"EventNo/I");
  m_tr->Branch("TrigFlag"        ,&m_TrigFlag        ,"TrigFlag/I");
  m_tr->Branch("CosmicTrig"      ,&m_CosmicTrig      ,"CosmicTrig/I");
  m_tr->Branch("LaserTrig"       ,&m_LaserTrig       ,"LaserTrig/I");
  m_tr->Branch("CVTrig"          ,&m_CVTrig          ,"CVTrig/I");
  m_tr->Branch("CosmicTrigFlagUp",&m_CosmicTrigFlagUp,"CosmicTrigFlagUp/I");
  m_tr->Branch("CosmicTrigFlagDn",&m_CosmicTrigFlagDn,"CosmicTrigFlagDn/I");
  m_tr->Branch("CVTrigFlag"      ,&m_CVTrigFlag      ,"CVTrigFlag/I");
  m_tr->Branch("LasertrigFlag"   ,&m_LaserTrigFlag   ,"LaserTrigFlag/I");
  m_tr->Branch("TimePeak"        ,&m_TimePeak        ,"TimePeak/D");
  m_tr->Branch("TimeSigma"       ,&m_TimeSigma       ,"TimeSigma/D");

  for( int i = 0; i< m_nModule; i++){
    mod[i]->Branch();
  }
  return true;
}
bool E14ConvWriter::SetBranchAddress(){
  std::cout<< "SetBranchAddress" << std::endl;
  m_tr->SetBranchAddress("RunNo"           ,&(this->m_RunNo));
  m_tr->SetBranchAddress("EventNo"         ,&m_EventNo);
  m_tr->SetBranchAddress("TrigFlag"        ,&m_TrigFlag);
  m_tr->SetBranchAddress("CosmicTrig"      ,&m_CosmicTrig);
  m_tr->SetBranchAddress("LaserTrig"       ,&m_LaserTrig);
  m_tr->SetBranchAddress("CVTrig"          ,&m_CVTrig);
  m_tr->SetBranchAddress("CosmicTrigFlagUp",&m_CosmicTrigFlagUp);
  m_tr->SetBranchAddress("CosmicTrigFlagDn",&m_CosmicTrigFlagDn);
  m_tr->SetBranchAddress("CVTrigFlag"      ,&m_CVTrigFlag);
  m_tr->SetBranchAddress("LaserTrigFlag"   ,&m_LaserTrigFlag);
  m_tr->SetBranchAddress("TimePeak"        ,&m_TimePeak);
  m_tr->SetBranchAddress("TimeSigma"       ,&m_TimeSigma);

  for( int i = 0; i< m_nModule; i++){
    std::cout << i << std::endl;
    mod[i]->SetBranchAddress();
  }
  return true;
}
bool E14ConvWriter::InitData(){

  m_TrigFlag         = 0;
  m_CosmicTrig       = 0;
  m_LaserTrig        = 0;
  m_CVTrig           = 0;
  m_CosmicTrigFlagUp = 0; 
  m_CosmicTrigFlagDn = 0;
  m_CVTrigFlag       = 0; 
  m_LaserTrigFlag    = 0;  
  m_TimePeak         = 0;
  m_TimeSigma        = 0;

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

int  E14ConvWriter::GetNmodule() const {
  return m_nModule;
}

int  E14ConvWriter::GetNsubmodule(int ModID) const {
  if( ModID >= m_nModule ){ return -1;}
  return (this->ModMap[ModID]).nMod;
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

bool E14ConvWriter::GetCFC( int ModID, int SubModID, int& CrateID, int& FADCID, int& ChannelID ){
  if( (ModID >= this->GetNmodule()) || 
      (SubModID >= (this->ModMap[ModID]).nMod) ){
    CrateID   = 9999;
    FADCID    = 9999;
    ChannelID = 9999; 
  }else{
    CrateID   = (this->ModMap[ModID]).Map[SubModID][0]; 
    FADCID    = (this->ModMap[ModID]).Map[SubModID][1]; 
    ChannelID = (this->ModMap[ModID]).Map[SubModID][2]; 
  }
  if( CrateID == 9999 || FADCID == 9999 || ChannelID == 9999 ){
    return false;
  }else{
    return true; 
  }
}
int E14ConvWriter::SetGraph( int ModID, int SubModID, E14ConvReader* conv[], TGraph* gr){
  gr->Set(0);
  int npoint  = 0;
  if( GetCFC( ModID, SubModID, m_tempCrateID, m_tempFADCID, m_tempChannelID) ){
    for( int ipoint = 0; ipoint < 48; ipoint++){
      gr->SetPoint( gr->GetN(), ipoint*8, 
		    conv[ m_tempCrateID ]->Data[ m_tempFADCID ][ m_tempChannelID ][ ipoint ] );
      npoint++;
    }
    return npoint;
  }else{
    return 0;
  }
}

bool E14ConvWriter::SetData( int iModID, int SubModID, E14ConvReader* conv[], Double_t Data[]){
  int npoint  = 0; 
  if( GetCFC( iModID, SubModID, m_tempCrateID, m_tempFADCID, m_tempChannelID ) ){
    for( int ipoint = 0; ipoint < 48; ipoint++){
      Data[ipoint] = conv[ m_tempCrateID ]->Data[ m_tempFADCID ][ m_tempChannelID ][ ipoint ];
	npoint++;
    }
  }
  if( npoint == 48 ) return true; 
  else return false; 
}


