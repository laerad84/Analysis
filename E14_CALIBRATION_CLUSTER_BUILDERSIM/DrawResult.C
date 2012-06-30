#include <fstream>
void DrawResult(){
  char* ANALIBDIR = std::getenv("ANALYSISLIB");
  char* libFile   = Form("%s/lib/libAnalysisLib.so",ANALIBDIR);
  gSystem->Load(libFile);
  IDHandler* handler = new IDHandler();
  CsIImage* image = new CsIImage(handler);

  std::ifstream ifsInit("SimRoot/BeamCalibrationFactor.dat");
  
  Double_t BeamCalibrationFactor[2716]={0};
  Double_t SimCalibrationFactor[2716]={0};
  Double_t cal;
  Int_t id;

  while( ifsInit>> id >> cal){
    BeamCalibrationFactor[id] = cal;
  }
  std::ifstream ifsResult("SimCalibration_Data/CalibrationFactor_15.dat");

  while( ifsResult >> id >> cal){
    SimCalibrationFactor[id] = cal;
  }

  TH1D* his = new TH1D( "his","",100,0.8,1.2);
  for( int i = 0; i< 2716; i++){

    if( SimCalibrationFactor[i] == 0 ) continue;
    if( BeamCalibrationFactor[i] == 0 ) continue;

    double x,y;
    handler->GetMetricPosition(i, x,y);
    double R = TMath::Sqrt(x*x+y*y);
    double ax = TMath::Abs(x);
    double ay = TMath::Abs(y);
    if( ax <150 && ay<150 ) continue;
    if( R >850 ) continue;
    if( ay > 550 )continue;
    if( ax > 550 && ax < 600 ) continue;
    his->Fill(SimCalibrationFactor[i]/BeamCalibrationFactor[i]);
    image->Fill(i,SimCalibrationFactor[i]/BeamCalibrationFactor[i]);
  }


  const int nIteration = 15;


  std::string   FilenameFactor = "CalibrationFactorADV_%d.dat";
  std::string   FilenameStatic = "CalibrationStaticsADV_%d.dat";
  std::ifstream ifsFactor[nIteration];
  std::ifstream ifsStatic[nIteration];

  Double_t Factor[nIteration][2716];
  Double_t Static[nIteration][2716];
  for( int i = 0; i< nIteration; ++i ){
    for( int j = 0; j< 2716; ++j){
      Factor[i][j] = 0;
      Static[i][j] = 0;
    }
  }


  for( int i =0 ;i< nIteration; ++i){
    std::cout<< Form(FilenameFactor.c_str(), i+1) << std::endl;
    ifsFactor[i].open(Form(FilenameFactor.c_str(), i+1));
    if( !ifsFactor[i].is_open() ){
      std::cout<< "File Open Error:"<<  i+1  << std::endl;
      return;
    }



  TCanvas* can = new TCanvas("can","",800,800);
  can->Divide(2,2);
  can->cd(1);
  his->Draw();
  can->cd(2);
  image->DrawWithRange("colz",0.95,1.05);
}
    
  
