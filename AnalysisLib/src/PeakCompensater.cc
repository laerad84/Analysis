#ifndef PEAKCOMPENSATER__H__
#include "PeakCompensater.h"
#endif //PEAKCOMPENSATER__H__


PeakCompensater::PeakCompensater(){
  if( !Init() ){
    std::cerr << " File is not exist" << std::endl;
  }
}

PeakCompensater::~PeakCompensater(){
  ;
}
bool PeakCompensater::Init(){
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  std::ifstream ifs(Form("%s/Data/getLinearityFunction.dat",ANALIBDIR.c_str()));
  if( !ifs.is_open()){ return false; }
  double x[3][50]={{0}};
  double y[3][50]={{0}};
  int ID;
  int nPoint;
  for( int i = 0; i < 3; i++){
    m_gr[i]  = new TGraph();
    //spl[i]->SetTitle(Form("Linearity%d",i));
  }
  for( int i = 0; i < 3; i++){
    ifs >> ID >> nPoint;
    for( int j = 0; j < nPoint; j++){
      ifs >> x[i][j] >> y[i][j];
      m_gr[i]->SetPoint(m_gr[i]->GetN(), x[i][j], y[i][j]);      
    }
  }
  for( int i = 0; i< 3; i++){
    m_spl[i]= new TSpline3(Form("Linearity%d",i), m_gr[i]);    
  }
  ifs.close();
  return true;
}

double PeakCompensater::Compensate(int id , double Peak ){
  int splID=-1;
  double CompensateOut=0;
  if( id == -1 ){
    splID = 1;
  }else if( id < 0 || id >2716){
    return -1;
  }else if( id < 2240 ){
    splID = 0;
  }else if( id >= 2240){
    splID = 2;
  }else{
    return -1;
  }
  if( Peak < 0. ){
    CompensateOut = -1;
  }else if( Peak > 15840 ){
    CompensateOut = Peak/m_spl[splID]->Eval(15840);
  }else{
    CompensateOut = Peak/m_spl[splID]->Eval(Peak);
  }
  return CompensateOut;
}
void PeakCompensater::Draw(int id, char* DrawOption){
  if( id >=3 ){ return ;}
  m_gr[id]->Draw(DrawOption);
}
