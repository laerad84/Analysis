#include <fstream>

void ResultDraw(){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(111111111);
  gStyle->SetOptStat("neRMI");
  char* ANALIBDIR = std::getenv("ANALYSISLIB");
  char* libFile   = Form("%s/lib/libAnalysisLib.so",ANALIBDIR);
  gSystem->Load(libFile);

  IDHandler* handler = new IDHandler();
  CsIImage* image  =new CsIImage(handler);
  
  const int nIteration = 11;
  
  std::ifstream ifsChamberFactor("/group/had/koto/ps/klea/work/jwlee/Analysis/E14_CALIBRATION_CLUSTER_BUILDERSIM/SimRoot/BeamCalibrationFactor.dat" );
  //std::ifstream ifsChamberFactor("ke3CalibConst.dat");
  Double_t FactorChamber[2716]={0};
  if( !ifsChamberFactor.is_open() ){
    std::cout<< "File Open Error" << std::endl;
    return;
  }

  int chId;
  double chamberCalFactor;
  while( ifsChamberFactor >> chId >> chamberCalFactor){
    std::cout<< chId << ":" << chamberCalFactor << std::endl;
    FactorChamber[chId] = chamberCalFactor;
  }
  
  //std::string   FilenameFactor = "SimCalibration_Data/CalibrationFactor_%d.dat";
  //std::string   FilenameStatic = "SimCalibration_Data/CalibrationStatics_%d.dat";  
  std::string   FilenameFactor = "CalibrationFactor_%d.dat";
  std::string   FilenameStatic = "CalibrationStatics_%d.dat";  
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
    ifsStatic[i].open(Form(FilenameStatic.c_str(),i+1));
    if( !ifsStatic[i].is_open() ){
      std::cout<< "File Open Error:" << i+1 << std::endl;
      return;
    }
  }



  int ID;
  Double_t CalFactor;
  
  for( int i = 0; i< nIteration; ++i){
    while( ifsFactor[i] >> ID >> CalFactor ){
      //std::cout << CalFactor << std::endl;      
      Factor[i][ID] = CalFactor; 
    }
  }
  for( int i = 0; i< nIteration; ++i){
    while( ifsStatic[i] >> ID >> CalFactor ){
      Static[i][ID] = CalFactor;
    }
  }
  
  TH2D* hisStaticsCal = new TH2D("his_StaticsCal","",30,0,1500,60,0.9,1.1);
  TH1D* hisStaticsCal_proj[30];
  TH1D* hisCalibChange[nIteration];
  TH1D* hisCalibrationDist[nIteration];
  TH1D* hisCalibrationCompare[nIteration];
  for( int i =0; i< nIteration; ++i){
    hisCalibChange[i] = new TH1D(Form("hisCalibration%d",i),"",100,0.95,1.05);    
    hisCalibrationDist[i] = new TH1D(Form("hisCalibrationDist%d",i),"",100,0.9,1.1);
    hisCalibrationCompare[i] = new TH1D(Form("hisCalibrationCompare%d",i),"",100,0.8,1.2);
  }
  TGraph* grAspect[nIteration];
  TGraph* grCompare[nIteration];
  for( int i = 0; i< nIteration; ++i){
    grCompare[i] = new TGraph();
    grAspect[i]  = new TGraph();
  }
    
  for( int i = 0; i<nIteration; ++i ){
    for( int j =0; j < 2716; ++j){
      double x,y;
      handler->GetMetricPosition( j , x,y );
      double R = TMath::Sqrt( x*x+ y*y);
      double ax = TMath::Abs(x);
      double ay = TMath::Abs(y);
      //Fiducial 
      if ( ax <150 && ay<150 )	continue;
      if( R> 850 ) continue;
      if( ay > 550 ) continue;
      int nStatics= 144;

      if( Factor[i][j] != 0 && i== nIteration-1 && Static[i][j] >= nStatics){
	hisStaticsCal->Fill(Static[i][j],Factor[i][j]/FactorChamber[j]);
	image->Fill(j,Factor[i][j]/FactorChamber[j]);
      }

      if( Factor[i][j] != 0 && Static[i][j]>= nStatics){
	hisCalibrationDist[i]->Fill(Factor[i][j]);
	hisCalibrationCompare[i]->Fill(Factor[i][j]/FactorChamber[j]);
	grCompare[i]->SetPoint(grCompare[i]->GetN(), FactorChamber[j],Factor[i][j]);
	grAspect[i]->SetPoint(grAspect[i]->GetN(),Static[i][j],Factor[i][j]/FactorChamber[j]);

	if( i==0 ){
	  hisCalibChange[i]->Fill(Factor[i][j]/1);
	}else{
	  hisCalibChange[i]->Fill(Factor[i][j]/Factor[i-1][j]);
	}
      }
    }
  }
  

  TGraph* grRMS  = new TGraph();
  TGraphErrors* grMean = new TGraphErrors();
  TGraphErrors* grMeanChange = new TGraphErrors();
  TGraphErrors* grMeanCompare = new TGraphErrors();
  TGraphErrors* grEntries     = new TGraphErrors();
  grRMS->SetMarkerStyle(5);
  grMean->SetMarkerStyle(5);
  grMeanCompare->SetMarkerStyle(4);

  for( int i = 0 ;i< nIteration; ++i){
    
    grMean->SetPoint(grMean->GetN(), i, hisCalibrationDist[i]->GetMean());
    grMean->SetPointError(grMean->GetN()-1,0,hisCalibrationDist[i]->GetMeanError());
    
    grMeanChange->SetPoint(grMeanChange->GetN(), i, hisCalibChange[i]->GetRMS()/hisCalibChange[i]->GetMean());
    //grMeanChange->SetPointError(grMeanChange->GetN()-1,0,hisCalibChange[i]->GetMeanError());
    grMeanCompare
      ->SetPoint(grMeanCompare->GetN(), i, 
		 hisCalibrationCompare[i]->GetRMS()/hisCalibrationCompare[i]->GetMean());
    grEntries->SetPoint(grEntries->GetN(), i, hisCalibrationCompare[i]->GetEntries());
    //grMeanCompare->SetPointError(grMeanCompare->GetN()-1, 0, hisCalibrationCompare[i]->GetMeanError());			    
  }
  
  for( int i = 0; i< 30; i++){
    hisStaticsCal_proj[i] = hisStaticsCal->ProjectionY(Form("his%d",i),i+1, i+2);
    grRMS->SetPoint(grRMS->GetN(),(i+1)*50,hisStaticsCal_proj[i]->GetRMS()); 
  }
  grMeanChange->SetMarkerStyle(7);
    

  TCanvas* can = new TCanvas("can","",1600,1200);  
  can->Divide(4,3);  
  can->cd(1);
  /*
  hisCalibChange[0]->Draw();
  for( int i = 0 ; i < nIteration; ++i ){
    hisCalibChange[i]->SetLineColor(i%5+1);
    hisCalibChange[i]->Draw("same");
    }*/
  hisCalibrationCompare[nIteration -1]->Draw();
  can->cd(2);
  grRMS->Draw("AP");
  can->cd (3);
  //grMean->Draw("AP");
  grAspect[nIteration-1]->Draw("AP");
  can->cd(4);
  /*
    hisCalibrationDist[0]->Draw();
    for( int i = 0; i< nIteration; ++i){
    hisCalibrationDist[i]->Draw("same");
    }
  */
  //grMeanChange->Draw("APL");
  grMeanCompare->Draw("AP");
  //grCompare[nIteration-1]->Draw("AP");
  can->cd(5);
  /*
  */
  hisStaticsCal->Draw("colz");
  TProfile* prof = hisStaticsCal->ProfileX();
  prof->Draw("same");
  can->cd(6);
  grEntries->SetMarkerStyle(4);
  grEntries->Draw("AP");
  can->cd(7);
  grMean->Draw("AP");
  can->cd(8);
  image->DrawWithRange("colz",0.95,1.05);

  can->SaveAs("Result.pdf");

  /*  TCanvas* canvas  =new TCanvas("canvas","",600,600);
  image->Draw();
  canvas->SaveAs("Image.gif");
  */
}
