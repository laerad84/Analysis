#ifndef PEAKCOMPENSATER__H__
#include "PeakCompensater.h"
#endif //PEAKCOMPENSATER__H__


PeakCompensater::PeakCompensater(){
  if( !Init() ){
    std::cerr << " File is not exist" << std::endl;
  }
  m_version = 0;
}
PeakCompensater::PeakCompensater( int version ){
  if( !Init( version ) ){
    std::cerr << " File is not exist" << std::endl;
  }
  m_version = version;
}
PeakCompensater::~PeakCompensater(){
  ;
}
bool PeakCompensater::Init( int version ){
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  SetMap();
  if( version == 0){ 
    std::ifstream ifs(Form("%s/Data/getLinearityFunction.dat",ANALIBDIR.c_str()));
    if( !ifs.is_open()){ return false; }
    double x[3][50]={{0}};
    double y[3][50]={{0}};
    int ID;
    int nPoint;
    for( int i = 0; i < 3; i++){
      m_gr[i]    = new TGraph();
      m_grInv[i] = new TGraph();
      //spl[i]->SetTitle(Form("Linearity%d",i));
    }
    for( int i = 0; i < 3; i++){
      ifs >> ID >> nPoint;
    for( int j = 0; j < nPoint; j++){
      ifs >> x[i][j] >> y[i ][j];
      m_gr[i]->SetPoint(m_gr[i]->GetN(), x[i][j], y[i][j]);      
      m_grInv[i]->SetPoint(m_grInv[i]->GetN(),x[i][j]/y[i][j],x[i][j]);
    }
    }
    for( int i = 0; i< 3; i++){
      m_spl[i]= new TSpline3(Form("Linearity%d",i), m_gr[i]);    
      m_splInv[i] = new TSpline3(Form("LineartyInv%d",i),m_grInv[i]);
    }
    ifs.close();
  }else if ( version == 1 ){
    TFile* tf = new TFile(Form("%s/Data/HeightLinearity_Laser.root",ANALIBDIR.c_str()));
    m_gr[0] = (TGraph*)tf->Get("heightLinearity_0");//for most of Small 
    m_gr[1] = (TGraph*)tf->Get("heightLinearity_1");//for some of Small
    m_gr[2] = (TGraph*)tf->Get("heightLinearity_2");//for All  of Large
    m_grInv[0] = new TGraph();
    m_grInv[2] = new TGraph();
    m_grInv[1] = new TGraph();
    m_grInv[0]->SetNameTitle("heigtLinearityInv_0","heigtLinearityInv_0");
    m_grInv[1]->SetNameTitle("heigtLinearityInv_1","heigtLinearityInv_1");
    m_grInv[2]->SetNameTitle("heigtLinearityInv_2","heigtLinearityInv_2");
    for( int i = 0; i< 3; i++){
      m_spl[i] = new TSpline3(Form("Linearity%d",i),m_gr[i]);      
      for( int j = 0; j< 160; j++){
	double x = j*100;
	double y = m_spl[i]->Eval(x);
	m_grInv[i]->SetPoint( m_grInv[i]->GetN(),x/y, x );	
      }
      m_splInv[i] = new TSpline3(Form("LinearityInv%d",i),m_grInv[i]);
    }
  }else if ( version == 2 ){
    //version 2 : no adjust // 
    for( int i = 0; i< 3; i++){
      m_gr[i] = new TGraph();
      m_grInv[i]= new TGraph();
    }
    for( int ipoint =0 ;ipoint < 160; ipoint++){
      m_gr[0]->SetPoint( ipoint , ipoint*100, 1);
      m_gr[1]->SetPoint( ipoint , ipoint*100, 1);
      m_gr[2]->SetPoint( ipoint , ipoint*100, 1);
    }
    for( int i = 0; i< 3; i++){
      m_spl[i] = new TSpline3(Form("Linearity%d",i),m_gr[i]);
      for( int j = 0; j< 160; j++){
	double x = j*100;
	double y = m_spl[i]->Eval(x);
	m_grInv[i]->SetPoint( m_grInv[i]->GetN(),x/y, x );	
      }
      m_splInv[i] = new TSpline3(Form("LinearityInv%d",i),m_grInv[i]);
    }
  }else if( version == 3){
    //version 3 : using graph from Calibration 
    TFile* tf0 = new TFile(Form("%s/Data/HeightLinearity_Calibration.root",ANALIBDIR.c_str()));
    m_gr[0] = (TGraph*)tf0->Get("heightLinearityCal");
    m_gr[1] = (TGraph*)tf0->Get("heightLinearityCal");

    TFile* tf1 = new TFile(Form("%s/Data/HeightLinearity_Laser.root",ANALIBDIR.c_str()));
    //m_gr[0] = (TGraph*)tf1->Get("heightLinearity_1");//for most of Small
    //m_gr[1] = (TGraph*)tf1->Get("heightLinearity_1");//for some of Small
    m_gr[2] = (TGraph*)tf1->Get("heightLinearity_2");//for All  of Large
    m_grInv[0] = new TGraph();
    m_grInv[2] = new TGraph();
    m_grInv[1] = new TGraph();
    m_grInv[0]->SetNameTitle("heigtLinearityInv_0","heigtLinearityInv_0");
    m_grInv[1]->SetNameTitle("heigtLinearityInv_1","heigtLinearityInv_1");
    m_grInv[2]->SetNameTitle("heigtLinearityInv_2","heigtLinearityInv_2");
    for( int i = 0; i< 3; i++){
      m_spl[i] = new TSpline3(Form("Linearity%d",i),m_gr[i]);      
      for( int j = 0; j< 160; j++){
	double x = j*100;
	double y = m_spl[i]->Eval(x);
	m_grInv[i]->SetPoint( m_grInv[i]->GetN(),x/y, x );	
      }
      m_splInv[i] = new TSpline3(Form("LinearityInv%d",i),m_grInv[i]);
    }
  }else if( version == 4 ){
        //version 3 : using graph from Calibration 
    TFile* tf0 = new TFile(Form("%s/Data/HeightLinearity_Calibration_nonInv.root",ANALIBDIR.c_str()));
    m_gr[0] = (TGraph*)tf0->Get("heightLinearityCal");
    m_gr[1] = (TGraph*)tf0->Get("heightLinearityCal");

    TFile* tf1 = new TFile(Form("%s/Data/HeightLinearity_Laser.root",ANALIBDIR.c_str()));
    //m_gr[0] = (TGraph*)tf1->Get("heightLinearity_1");//for most of Small
    //m_gr[1] = (TGraph*)tf1->Get("heightLinearity_1");//for some of Small
    m_gr[2] = (TGraph*)tf1->Get("heightLinearity_2");//for All  of Large
    m_grInv[0] = new TGraph();
    m_grInv[2] = new TGraph();
    m_grInv[1] = new TGraph();
    m_grInv[0]->SetNameTitle("heigtLinearityInv_0","heigtLinearityInv_0");
    m_grInv[1]->SetNameTitle("heigtLinearityInv_1","heigtLinearityInv_1");
    m_grInv[2]->SetNameTitle("heigtLinearityInv_2","heigtLinearityInv_2");
    for( int i = 0; i< 3; i++){
      m_spl[i] = new TSpline3(Form("Linearity%d",i),m_gr[i]);
      for( int j = 0; j< 160; j++){
	double x = j*100;
	double y = m_spl[i]->Eval(x);
	m_grInv[i]->SetPoint( m_grInv[i]->GetN(),x/y, x );	
      }
      m_splInv[i] = new TSpline3(Form("LinearityInv%d",i),m_grInv[i]);
    }
  }else{ return false; }
 return true;
}
double PeakCompensater::Compensate(int id , double Peak ){
  int splID=-1;
  double CompensateOut=0;  
  switch (m_version){
  case 0 :
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
    break;
  case 1 :
    if( id ==  -1  ){ splID = -1;}
    if( id >= 2716){ splID = -1;}
    if( id < 2240 ){ 
      if( m_map[id][0] == 7 ){
	splID = 1;
      }else{
	splID = 0;
      }
    }else{
      splID = 2;
    }
    //splID = m_map[id][3];    
    if( splID <0 || splID > 2 ){ 
      CompensateOut = 0;
    }
    if( Peak < 0. ){ CompensateOut =  0;}
    else if( Peak > 15840 ){
      CompensateOut = Peak/m_spl[splID]->Eval(15840);
    }else{
      CompensateOut = Peak/m_spl[splID]->Eval(Peak);
    }
    break;
  case 2 :
    if( id ==  -1  ){ splID = -1;}
    if( id >= 2716){ splID = -1;}
    if( id < 2240 ){ 
      if( m_map[id][0] == 7 ){
	splID = 1;
      }else{
	splID = 0;
      }
    }else{
      splID = 2;
    }
    //splID = m_map[id][3];    
    if( splID <0 || splID > 2 ){ 
      CompensateOut = 0;
    }
    if( Peak < 0. ){ CompensateOut =  0;}
    else if( Peak > 15840 ){
      CompensateOut = Peak/m_spl[splID]->Eval(15840);
    }else{
      CompensateOut = Peak/m_spl[splID]->Eval(Peak);
    }
    break;    
  case 3 :
    if( id == -1 ){ splID = -1;}
    if( id > 2716){ splID = -1;}
    if( id < 2240){ splID = 0; }
    if( id >=2240 && id< 2716 ){ splID = 2;}
    if( Peak < 0. ){ CompensateOut = 0; }
    else if( Peak > 15840 ){ 
      CompensateOut = Peak/m_spl[splID]->Eval(15840);
    }else{
      CompensateOut = Peak/m_spl[splID]->Eval(Peak);
    }
    break;
  case 4 :
    if( id == -1 ){ splID = -1;}
    if( id > 2716){ splID = -1;}
    if( id < 2240){ splID = 0; }
    if( id >=2240 && id< 2716 ){ splID = 2;}
    if( Peak < 0. ){ CompensateOut = 0; }
    else if( Peak > 15840 ){ 
      CompensateOut = Peak*m_spl[splID]->Eval(15840);
    }else{
      CompensateOut = Peak*m_spl[splID]->Eval(Peak);
    }
    break;

  default :
    CompensateOut = 0;
    break; 
  }
  return CompensateOut;
}

double PeakCompensater::InvCompensate( int id, double CompensateOut){
  //This function was made for recover peakHeight.
  //result is equal ~0.1% level.
  //No Compensatable for waveforms which height is over 15000 cnt.
  int splID=-1;
  double Peak = 0.;
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
  if( CompensateOut < 0. ){
    CompensateOut = 0;
  }else{
    Peak = m_splInv[splID]->Eval(CompensateOut);
  }
  return Peak;
}

void PeakCompensater::Draw(int id, char* DrawOption){
  if( id >=3 ){ return ;}
  m_gr[id]->Draw(DrawOption);
}
void PeakCompensater::SetMap(){
  std::string ANALIBDIR = std::getenv("ANALYSISLIB");
  std::ifstream ifs(Form("%s/Data/ch_map_CsI_L1.txt",ANALIBDIR.c_str()));
  int tmpID;
  int tmpCrate;
  int tmpSlot;
  int tmpChannel;
  int tmpL1;
  for( int im_ch = 0; im_ch < 2716; im_ch++){
    m_map[im_ch][0] = -1;
    m_map[im_ch][1] = -1;
    m_map[im_ch][2] = -1;
    m_map[im_ch][3] = -1;
  }
  while(ifs >> tmpID >> tmpCrate >> tmpSlot >> tmpChannel >> tmpL1){
    if( tmpID >= 2716 || tmpCrate > 50 ){ continue; }
    m_map[tmpID][0] = tmpCrate;
    m_map[tmpID][1] = tmpSlot;
    m_map[tmpID][2] = tmpChannel;
    m_map[tmpID][3] = tmpL1;
  }
}
