void DrawWaveform(){
  gStyle->SetOptFit(1111111111);
  TChain* ch = new TChain("Waveform","");
  ch->Add("/Volume0/ExpData/2012_Feb_Beam/RootFile_wav/TEMPLATE_FIT_RESULT_DUMP_4503.root");
  ch->Add("/Volume0/ExpData/2012_Feb_Beam/RootFile_wav/TEMPLATE_FIT_RESULT_DUMP_4504.root");
  
  TH2D* his1 = new TH2D("his1","his1",600,-250,350,150,-0.25,1.25);
  TH2D* his2 = new TH2D("his2","his2",600,-250,350,150,-0.25,1.25);
  ch->Project(his1->GetName(),"(Waveform-Pedestal)/Height:TimeInfo-PeakTime",
	      "PeakTime>100&&PeakTime<325&&Height>100&&Height<250&&ChisqNDF<30&&ModuleNumber==1220");
  ch->Project(his2->GetName(),"(Waveform-Pedestal)/Height:TimeInfo-PeakTime",
	      "PeakTime>100&&PeakTime<325&&Height<1000&&Height>500&&ChisqNDF<30&&ModuleNumber==1220");
  TH1D* hisHeight = new TH1D("hisHeight","hisHeight",1000,0,5000);
  ch->Project(hisHeight->GetName(),"Height","PeakTime>100&&PeakTime<325&&ModuleNumber==1220&&ChisqNDF<30");
  
  TProfile* prof1 = his1->ProfileX();
  TProfile* prof2 = his2->ProfileX();
  TH1D*     hisFitHeight1= his1->ProjectionY("hisp1",250,251);
  TH1D*     hisFitHeight2= his2->ProjectionY("hisp2",250,251);
  hisFitHeight1->Fit("gaus","","",0.95,1.05);
  hisFitHeight2->Fit("gaus","","",0.95,1.05);
  TCanvas* can = new TCanvas("can","",1200,800);
  can->Divide(3,2);
  can->cd(1);
  his1->Draw("col");
  can->cd(2);
  his2->Draw("col");
  can->cd(3);
  prof1->Draw();
  prof2->SetLineColor(2);
  prof2->Draw("same");
  can->cd(4);
  hisFitHeight1->SetLineColor(1);
  hisFitHeight2->SetLineColor(2);
  hisFitHeight1->Draw();
  hisFitHeight2->Draw("same");
  can->cd(5);
  hisHeight->Draw();
}

  
  
