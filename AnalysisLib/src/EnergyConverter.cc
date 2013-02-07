#include "EnergyConverter.h"

EnergyConverter::EnergyConverter(){
  // For CsI Main Detector // 
  std::cerr << "Converter Build" << std::endl;
  m_fInit = false;
  m_fRead = false;
  Init();
  m_nChannel = 2716;
  m_DetectorName = "CsI";
  m_Compensater  = new PeakCompensater();
}  
EnergyConverter::EnergyConverter(const char* DetectorName , int nChannel ){
  // For Other Detector 
  std::cerr << "Converter Build" << std::endl;
  m_fInit        = false;
  m_fRead        = false; 
  Init(); 
  m_nChannel     = nChannel;
  m_DetectorName = DetectorName;
  m_Compensater  = NULL;
}
EnergyConverter::~EnergyConverter(){
  std::cerr << "Converter Deleted" << std::endl; 
  delete m_Compensater;
  m_Compensater = NULL;
}
void   EnergyConverter::Init(){
  if( m_fInit ){ return; }
  m_CalFilename.clear();
  for( int i = 0; i< 4096; i++){
    m_CalibrationConstant[ i ]  = 0;
    m_fCalibrationConstant[ i ] = false;
  }
  m_Compensater = new PeakCompensater();
}
void   EnergyConverter::Reset(){
  for( int i = 0; i< 4096; i++){
    m_CalibrationConstant[ i ]  = 0;
    m_fCalibrationConstant[ i ] = false;
  }  
}
bool   EnergyConverter::ReadCalibrationTextFile( const char* CalibrationFilename){
  Reset();
  m_CalFilename = CalibrationFilename;
  std::ifstream ifs( m_CalFilename.c_str());
  if( !ifs.is_open()) { m_fRead = false; return m_fRead;}
  else{ m_fRead = true; }
  int    ID;
  double CalConstant;
  while( ifs >> ID >> CalConstant){
    m_CalibrationConstant[ID]  = CalConstant;
    m_fCalibrationConstant[ID] = true;
  }
  return m_fRead;
}
bool   EnergyConverter::ReadCalibrationRootFile( const char* CalibrationFilename){
  Reset();
  m_CalFilename = CalibrationFilename;
  TFile* tfCal = new TFile( m_CalFilename.c_str());
  TTree* trCal = NULL;
  trCal        = (TTree*)tfCal->Get("GainFitPar");
  if( tfCal->IsOpen() ){
    m_fRead = true;
  }else{ 
    m_fRead = false;
    return  m_fRead;
  }
  if( trCal != NULL){
    m_fRead = true; 
  }else{
    m_fRead = false;
    return m_fRead;
  }
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
  return m_fRead;
}
double EnergyConverter::GetCalibrationConstant ( int channelNumber ) const {
  if( channelNumber > m_nChannel || channelNumber < 0 || !m_fRead ){    
    return 0;
  }else{
    return m_CalibrationConstant[ channelNumber ];    
  }  
}
bool   EnergyConverter::IsGoodChannel          ( int channelNumber ) const {
  return m_fCalibrationConstant[ channelNumber ];
}
double EnergyConverter::ConvertToEnergy        ( int channelNumber, double PeakHeight ) const {
  double Energy = 0;
  double CompensatedHeight = 0;
  if( m_Compensater == NULL ){ 
    Energy            = this->GetCalibrationConstant( channelNumber )*PeakHeight;
  }else{
    // In Case of CsI // 
    CompensatedHeight = m_Compensater->Compensate( channelNumber , PeakHeight );
    Energy            = this->GetCalibrationConstant( channelNumber ) * CompensatedHeight;
  }
  return Energy;
}
double EnergyConverter::ConvertToHeight       ( int channelNumber, double Energy ) const  {
  double PeakHeight = 0;
  double CompensatedHeight = 0;
  if( Energy < 0 ) { return 0;}
  if( m_Compensater == NULL ){
    PeakHeight  = Energy/(this->GetCalibrationConstant( channelNumber ));
  }else{
    // In Case of CsI // 

    CompensatedHeight = Energy /(this->GetCalibrationConstant( channelNumber ));
    PeakHeight = m_Compensater->InvCompensate( channelNumber, CompensatedHeight);
    //std::cout<< CompensatedHeight << "\t" << PeakHeight << std::endl;
  }
  return PeakHeight;
}
