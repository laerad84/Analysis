#include "EnergyConverter.h"

EnergyConverter::EnergyConverter( ){
  std::cerr << "Converter Build" << std::endl;
  m_fInit = false;
  m_fRead = false;
  Init();
  m_nChannel = 2716;
  std::string ANALIBDIR ="";
  ANALIBDIR = std::getenv("ANALYSISLIB");
  m_CalFilename = Form("%s/Data/Calibration_Result_0000.dat",ANALIBDIR.c_str());
  bool rst = ReadCalibrationFile( m_CalFilename.c_str() );
  if( !rst ){
    std::cerr << __FUNCTION__ << ":" << __LINE__ << std::endl; 
    std::cerr << "Fail to read file" << std::endl;
  }
}  

EnergyConverter::EnergyConverter(const char* CalFilename , int nChannel ){
  std::cerr << "Converter Build" << std::endl;
  m_fInit        = false;
  m_fRead        = false; 
  Init(); 
  m_nChannel     = nChannel;
  m_CalFilename  = CalFilename; 
  bool rst = ReadCalibrationFile( CalFilename );
  if( !rst ){
    std::cerr << __FUNCTION__ << ":" << __LINE__ << std::endl; 
    std::cerr << "Fail to read file" << std::endl;
  }
}

EnergyConverter::~EnergyConverter(){
  std::cerr << "Converter Deleted" << std::endl; 
  delete m_Compensater;
  m_Compensater = NULL;
}

void EnergyConverter::Init(){
  if( m_fInit ){ return; }
  m_CalFilename.clear();
  m_nChannel = 0;
  for( int i = 0; i< 4096; i++){
    m_CalibrationConstant[ i ]  = 0;
    m_fCalibrationConstant[ i ] = false;
  }
  m_Compensater = new PeakCompensater();
}

void EnergyConverter::Reset(){
  m_nChannel = 0;
  for( int i = 0; i< 4096; i++){
    m_CalibrationConstant[ i ] = 0;
    m_fCalibrationConstant[ i ] = false;
  }  
}
bool EnergyConverter::ReadCalibrationFile(const  char* ){
  Reset();
  TFile* tfCal = new TFile( m_CalFilename.c_str());
  TTree* trCal = (TTree*)tfCal->Get("GainFitPar");
  int ID;
  double Peak;
  trCal->SetBranchAddress("ID",&ID);
  trCal->SetBranchAddress("Peak",&Peak);
  for( int i = 0; i< trCal->GetEntries(); i++){
    trCal->GetEntry(i);
    if( Peak > 0){
      m_CalibrationConstant[ID] = 14.014/Peak;
      m_fCalibrationConstant[ID] = true;
    }
  }
  trCal->Delete();
  tfCal->Close();
  /*
  std::ifstream ifs( m_CalFilename.c_str());
  if( !ifs.is_open()) { m_fRead = false; return m_fRead;}
  else{ m_fRead = true; }
  int    ID;
  double CalConstant;
  while( ifs >> ID >> CalConstant){
    m_CalibrationConstant[ID]  = CalConstant;
    m_fCalibrationConstant[ID] = true;
    }
  */
  return m_fRead;
}

double EnergyConverter::GetCalibrationConstant( int channelNumber ) const {
  if( channelNumber > m_nChannel || channelNumber < 0 || !m_fRead ){
    return 0;
  }else{
    return m_CalibrationConstant[ channelNumber ];    
  }  
}

bool EnergyConverter::IsGoodChannel( int channelNumber ) const {  
  return m_fCalibrationConstant[ channelNumber ];
}

double EnergyConverter::ConvertToEnergy( int channelNumber, double PeakHeight ) const {
  if( m_Compensater == NULL ){ 
    std::cerr << "Compensater is not exist." << std::endl;
    return 0;
  }
  double CompensatedHeight = m_Compensater->Compensate( channelNumber , PeakHeight );
  double Energy            = this->GetCalibrationConstant( channelNumber )*CompensatedHeight;
  return Energy;
}
