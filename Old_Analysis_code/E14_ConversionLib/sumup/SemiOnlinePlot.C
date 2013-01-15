 #include "SemiOnlinePlot.h"

 #if !defined(__CINT__)
 ClassImp(SemiOnlinePlot)
 #endif

 SemiOnlinePlot::SemiOnlinePlot(){
   threshold = 2300/5;   // 2300count = 1MIP (15MeV)
   EneThreshold = 3;   // 3 MeV
   //  threshold = 0;
   TargetCH = 3;  

   NDet = 7;
   DetNumber[0] = 2716;
   DetNumber[1] = 32;
   DetNumber[2] = 44;
   DetNumber[3] = 10;
   DetNumber[4] = 20;
   DetNumber[5] = 5;
   DetNumber[6] = 1;
   
   ReadLaserRefValue();
 }

 SemiOnlinePlot::~SemiOnlinePlot(){

   delete adcDistribution;
   delete adcDistributionLaser;
   delete energyDistribution;
   delete hitPosition;
   delete eventDisplay;
 }


 void SemiOnlinePlot::Init(int runID){

   RunNum = runID;

   hfile = new TFile( Form("/disk/kotodaq/2012Feb/run%d/RootFile/SemiOnlinePlot_%d.root",runID,runID),"recreate" );

   handler = new IDHandler("/home/koto/togawa/production/2012VME/sumup/crystal.txt");

   TStyle* gStyle = new TStyle();
   gStyle->SetPalette(1);

   adcDistribution         = new CsIImage(handler);
   adcDistributionLaser    = new CsIImage(handler);
   energyDistribution      = new CsIImage(handler);
   hitPosition             = new CsIImage(handler);
   eventDisplay            = new CsIImage(handler);

   // For eventdisplay
   adcDistribution->SetTitle("ADC distribution (except for Laser);;;Integrated ADC");
   adcDistributionLaser->SetTitle("ADC distribution (Laser);;;Integrated ADC");
   energyDistribution->SetTitle("Energy distribution (except for Laser);;;Energy (MeV)");
   hitPosition->SetTitle("Hit position");

   // For Laser Monitor
   for(int ipin=0;ipin<4;ipin++){
     hPinDiode[ipin] = new TH1F( Form("hPinDiode_%d",ipin),Form("Pin-diode out %d",ipin),400,0,4000);
   }
   for(int icrystal=0;icrystal<DetNumber[0];icrystal++){
     hLaserRef[icrystal] = new TH1F( Form("hLaserRef_%d",icrystal),"",5000,0,10000);
     hLaserRefNorm[icrystal] = new TH1F( Form("hLaserRefNorm_%d",icrystal),"",5000,0,5);
     hPedNoise[icrystal] = new TH1F( Form("hPedNoise_%d",icrystal),Form("Energy distribution above 1 MeV for crystal#%d ; Energy (MeV)",icrystal),100,0,10);
   }

   // For CV monitor
   for(int icv=0;icv<10;icv++){
     hCV[icv] = new TH1F( Form("hCV_%d",icv),Form("CV ch:%d; ADC",icv),500,0,50000);
   }

   // For Et Monitor
   hIntegratedADCSum = new TH1F("hIntegratedADCSum","Sum of Integrated ADC; ADC",300,0,300000);
   hEneSum = new TH1F("hEneSum","Sum of Energy; Energy (MeV)",300,0,3000);
   hEtSum = new TH2D("hEtSum","Sum of ADC in each sample (Except for Laser); Sample; ADC",48,0,48,2000,1.1e+6,1.2e+6);
   hEtSumLaser = new TH2D("hEtSumLaser","Sum of ADC in each sample (Laser); Sample; ADC",48,0,48,2000,1.1e+6,2.0e+6);
   //  hEtSumLaser = new TH2D("hEtSumLaser","Sum of ADC in each sample (Laser)",48,0,48,2000,1.1e+6,30e+6);

   c1 = new TCanvas( "c1", "c1", 800, 600 );
   c2 = new TCanvas( "c2", "c2", 800, 600 );

   sprintf( psFileName1, "/disk/kotodaq/2012Feb/run%d/SemiOnlinePlot_%d.ps", runID , runID);
   c1->Print( Form( "%s[", psFileName1 ) );
   c1->SetRightMargin(0.2);

   sprintf( psFileName2, "/disk/kotodaq/2012Feb/run%d/PedNoise_%d.ps", runID , runID);
   c2->Print( Form( "%s[", psFileName2 ) );
   c2->SetRightMargin(0.2);

   tline = new TLine();
   tline -> SetLineStyle(2);
   tline -> SetLineColor(2);

   ttex = new TLatex();

 }

 void SemiOnlinePlot::SetIntegratedADC( int idet, int ich, Double_t Ene ){
   DetIntegratedADC[idet][ich] = Ene;
   NOfIntegratedADC[idet]++;
 }

 void SemiOnlinePlot::SetPeakHeight( int idet, int ich, Short_t Peak ){
   DetPeakHeight[idet][ich] = Peak;
 }

 void SemiOnlinePlot::SetDetEne( int idet, int ich, Double_t Ene ){
   //  printf("%d %d %f\n",NDet,DetNumber[NDet],Ene);
   DetEne[idet][ich] = Ene;
   NOfDetEne[idet]++;
 }

 void SemiOnlinePlot::SetEtSum( int icrate, int ifadc, Int_t *ADC ){
   //  printf("%d %d %f\n",NDet,DetNumber[NDet],Ene);
   for(int isample=0;isample<48;isample++){
     EtSum[icrate][ifadc][isample] = ADC[isample];
   }
 }

 void SemiOnlinePlot::InitValue()
 {
   for(int idet=0;idet<NDet;idet++){
     NOfIntegratedADC[idet] = 0;
     NOfDetEne[idet] = 0;
     for(int icrystal=0; icrystal<DetNumber[idet]; icrystal++){
       DetEne[idet][icrystal] = 0;
     }
   }
 }

 void SemiOnlinePlot::Fill( int event )
 {

   //  printf("%d %d\n",NOfIntegratedADC[0],NOfDetEne[0]);

   bool LASERTAG = DetEne[5][0] > 5000;
   bool BRANKTAG = DetEne[5][1] > 10000;

   if( LASERTAG && !BRANKTAG ){

     for(int ipin=0;ipin<4;ipin++){
       hPinDiode[ipin] -> Fill( DetEne[5][ipin+1]);
     }
     for(int icrystal=0;icrystal<DetNumber[0];icrystal++){
       if( LaserRefMean[icrystal] == 0 ) continue;
       double x = 0;
       double y = 0;
       handler -> GetMetricPosition( icrystal, x, y);
       if( x<0 && y>0 ){
	 int pind = DetPeakHeight[5][4];
	 if( pind == 0 ) continue;
	 hLaserRef[icrystal] -> Fill( DetPeakHeight[0][icrystal] / pind );
	 hLaserRefNorm[icrystal]
	   -> Fill( DetPeakHeight[0][icrystal] / pind / LaserRefMean[icrystal] );
       }
       if( x<0 && y<0 ){
	 int pind = DetPeakHeight[5][3];
	 if( pind == 0 ) continue;
	 hLaserRef[icrystal] -> Fill( DetPeakHeight[0][icrystal] / pind );
	 hLaserRefNorm[icrystal]
	   -> Fill( DetPeakHeight[0][icrystal] / pind / LaserRefMean[icrystal] );
       }
       if( x>0 && y>0 ){
	 int pind = DetPeakHeight[5][2];
	 if( pind == 0 ) continue;
	 hLaserRef[icrystal] -> Fill( DetPeakHeight[0][icrystal] / pind );
	 hLaserRefNorm[icrystal] 
	   -> Fill( DetPeakHeight[0][icrystal] / pind / LaserRefMean[icrystal] );
       }
       if( x>0 && y<0 ){
	 int pind = DetPeakHeight[5][1];
	 if( pind == 0 ) continue;
	 hLaserRef[icrystal] -> Fill( DetPeakHeight[0][icrystal] / pind );
	 hLaserRefNorm[icrystal]
	   -> Fill( DetPeakHeight[0][icrystal] / pind / LaserRefMean[icrystal] );
       }
     }

     for(int isample=0;isample<48;isample++){
       double EtSumSample = 0;
       for(int icrate=0;icrate<10;icrate++){
	 for(int ifadc=0;ifadc<16;ifadc++){	
	   EtSumSample += EtSum[icrate][ifadc][isample];
	 }
       }
       //    printf("%f\n",EtSumSample);
       hEtSumLaser -> Fill(isample, EtSumSample);
     }

   }else{

     // Veto Laser trigger twice
     if( DetEne[5][0] > -1000 ){
       
       double IntegratedADCSum = 0;
       double EneSum = 0;
       for(int icrystal=0; icrystal<DetNumber[0]; icrystal++){
	 IntegratedADCSum += DetIntegratedADC[0][icrystal];
       }
       for(int icrystal=0; icrystal<DetNumber[0]; icrystal++){
	 if( DetEne[0][icrystal] > 3 ) EneSum += DetEne[0][icrystal];
	 if( DetEne[0][icrystal] > 1 ) hPedNoise[icrystal] -> Fill( DetEne[0][icrystal] );
       }
       
       hIntegratedADCSum -> Fill( IntegratedADCSum );
       hEneSum -> Fill( EneSum );
       //  printf("%f\n",IntegratedADCSum);
       
       for(int isample=0;isample<48;isample++){
	 double EtSumSample = 0;
	 for(int icrate=0;icrate<10;icrate++){
	   for(int ifadc=0;ifadc<16;ifadc++){	
	     EtSumSample += EtSum[icrate][ifadc][isample];
	   }
	 }
	 //    printf("%f\n",EtSumSample);
	 hEtSum -> Fill(isample, EtSumSample);

       }
     }
   }

   for(int icv=0;icv<10;icv++){
     hCV[icv] -> Fill( DetIntegratedADC[3][icv]);
   }

 }

 void SemiOnlinePlot::RunDisplay()
 {
   bool LASERTAG = DetEne[5][0] > 5000;

   // CsI
   for( int crystalID=0; crystalID<DetNumber[0]; crystalID++ ){
     if( LASERTAG ){
       if( DetIntegratedADC[0][crystalID] > threshold ){
	 adcDistributionLaser -> Fill( crystalID, DetIntegratedADC[0][crystalID] );
       }
     }else{

       if( DetEne[5][0] > -1000){
	 if( DetIntegratedADC[0][crystalID] > threshold ){
	   adcDistribution -> Fill( crystalID, DetIntegratedADC[0][crystalID] );
	 }
	 if( DetEne[0][crystalID] > 3 ){
	   energyDistribution -> Fill( crystalID, DetEne[0][crystalID] );
	   hitPosition        -> Fill( crystalID, 1 );
	 }
       }

     }
   }
   
   // CC03
   for( int crystalID=0; crystalID<DetNumber[1]; crystalID++ ){
     if( DetIntegratedADC[1][crystalID] > 800 ){
       if( LASERTAG ){
	 adcDistributionLaser -> Fill( (crystalID+28)%32+3000, DetIntegratedADC[1][crystalID]/2 );
       }else{

	 if( DetEne[5][0] > -1000){
	   adcDistribution -> Fill( (crystalID+28)%32+3000, DetIntegratedADC[1][crystalID]/2 );
	 }

       }
       if( DetEne[5][0] > -1000){
	 hitPosition        -> Fill( (crystalID+28)%32+3000, 1 );
       }
     }
   }

 }

 void SemiOnlinePlot::EventDisplay(int event )
 {
   eventDisplay->Reset();
   eventDisplay->SetTitle( Form("Event %d (Z axis = energy in MeV)", event) );

   // CsI
   for( int crystalID=0; crystalID<DetNumber[0]; crystalID++ ){
     if( DetEne[0][crystalID] > EneThreshold ){
       eventDisplay->Fill( crystalID, DetEne[0][crystalID] );      
     }
   }
   // CC03
   for( int crystalID=0; crystalID<DetNumber[1]; crystalID++ ){
     if( DetIntegratedADC[1][crystalID] > 800 ){
       //       eventDisplay->Fill( crystalID+3000, DetIntegratedADC[1][crystalID] );
       eventDisplay->Fill( (crystalID+28)%32+3000, DetIntegratedADC[1][crystalID]/2/100 );
     }
   }

   c1 -> cd();
   eventDisplay -> Draw();
   c1 -> Print( psFileName1 );
   c1 -> Clear();

 }

 void SemiOnlinePlot::DrawAll()
 {

   double x[2716]  = {0};
   double xE[2716] = {0};
   double Mean[2716]  = {0};
   double MeanE[2716] = {0};

   // For Laser gain
   int LaserEventNum = (int)hLaserRefNorm[0] -> GetEntries();
   printf("Laser event num = %d\n",LaserEventNum);
   for(int icrystal=0;icrystal<DetNumber[0];icrystal++){
     x[icrystal] = icrystal+0.5;
     Mean[icrystal]  = hLaserRefNorm[icrystal] -> GetMean();
     if( LaserEventNum == 0 ) continue;
     Mean[icrystal]  = hLaserRefNorm[icrystal] -> GetMean();
     MeanE[icrystal] = hLaserRefNorm[icrystal] -> GetRMS() / sqrt(LaserEventNum);
   }
   for(int i=0;i<35;i++){
     gLaserStabilityS[i] = new TGraphErrors(64,x,&Mean[i*64],xE,&MeanE[i*64]);
     gLaserStabilityS[i] -> SetName( Form("gLaserStabilityS%d",i) );
     gLaserStabilityS[i] -> SetTitle( Form("Laser normalized by pin-diode for Small in run %d crystal#%d-%d ; Sample# ; Mean (Error = RMS/#sqrt{%d} )",RunNum,i*64,(i+1)*64-1, LaserEventNum ));
     gLaserStabilityS[i] -> SetMaximum(1.05);
     gLaserStabilityS[i] -> SetMinimum(0.95);
     gLaserStabilityS[i] -> SetMarkerStyle(20);
     gLaserStabilityS[i] -> SetMarkerSize(1);
   }
   for(int i=0;i<4;i++){
     gLaserStabilityL[i] = new TGraphErrors(64,x,&Mean[i*64+2360],xE,&MeanE[i*64+2360]);
     gLaserStabilityL[i] -> SetName( Form("gLaserStabilityL%d",i) );
     gLaserStabilityL[i] -> SetTitle( Form("Laser normalized by pin-diode for Large in run %d crystal#%d-%d ; Sample# ; Mean (Error = RMS/#sqrt{%d} )",RunNum,i*64+2360,(i+1)*64-1+2360, LaserEventNum ));
     gLaserStabilityL[i] -> SetMaximum(1.05);
     gLaserStabilityL[i] -> SetMinimum(0.95);
     gLaserStabilityL[i] -> SetMarkerStyle(20);
     gLaserStabilityL[i] -> SetMarkerSize(1);
   }

   c1 -> cd();
   c1 -> Divide(2,2);
   c1 -> GetPad(1) -> SetLogy();
   c1 -> GetPad(2) -> SetLogy();
   c1 -> cd(1);
   hIntegratedADCSum -> Draw();
   c1 -> cd(2);
   hEneSum -> Draw();
   c1 -> cd(3);
   hEtSum -> Draw("COLTZ");
   c1 -> cd(4);
   hEtSumLaser -> Draw("COLTZ");
   c1 -> Print( psFileName1 );
   c1 -> Clear();

   c1 -> cd();
   c1 -> Divide(4,3);
   for(int icv=0;icv<10;icv++){
     c1 -> GetPad(icv+1) -> SetLogy();
     c1 -> cd( icv+1 );
     hCV[icv] -> Draw();
   }
   c1 -> Print( psFileName1 );
   c1 -> Clear();

   c1 -> cd();
   c1 -> Divide(2,2);
   for(int ipin=0;ipin<4;ipin++){
     c1 -> cd( ipin+1 );
     hPinDiode[ipin] -> Draw();
   }
   c1 -> Print( psFileName1 );
   c1 -> Clear();

   for(int i=0;i<9;i++){
     c1 -> cd();
     c1 -> Divide(1,4);
     for(int j=0;j<4;j++){
       if( i==8 && j==3 ) continue;
       c1 -> cd( j+1 );
       c1 -> GetPad( j+1 ) -> SetGridy();;
       gLaserStabilityS[i*4+j] -> Draw("AP");
     }
     c1 -> Print( psFileName1 );
     c1 -> Clear();
   }
   c1 -> cd();
   c1 -> Divide(1,4);
   for(int j=0;j<4;j++){
     c1 -> cd( j+1 );
     c1 -> GetPad( j+1 ) -> SetGridy();;
     gLaserStabilityL[j] -> Draw("AP");
   }
   c1 -> Print( psFileName1 );
   c1 -> Clear();

   c1 -> cd();
   adcDistribution -> Draw();
   c1 -> Print( psFileName1 );
   c1 -> Clear();

   c1 -> cd();
   energyDistribution -> Draw();
   c1 -> Print( psFileName1 );
   c1 -> Clear();

   c1 -> cd();
   adcDistributionLaser -> Draw();
   c1 -> Print( psFileName1 );
   c1 -> Clear();

   c1 -> cd();
   hitPosition -> Draw();
   c1 -> Print( psFileName1 );
   c1 -> Clear();

   c1->Print( Form( "%s)", psFileName1 ) );

   // Pedestal noise
   // Small
   for(int i=0;i<140;i++){
     c2 -> cd();
     c2 -> Divide(4,4);
     for(int j=0;j<16;j++){
       c2 -> cd( j+1 );
       hPedNoise[i*16+j] -> Draw();
       tline -> DrawLine(3,0,3,hPedNoise[i*16+j]->GetMaximum() );
       if( i*16+j==553 | i*16+j==506 | i*16+j==959 | i*16+j==419 |
	   i*16+j==471 | i*16+j==1898 | i*16+j==1899 | i*16+j==191 | 
	   i*16+j==1688 | i*16+j==1790 | i*16+j==1799 | i*16+j==2113){
	 ttex -> DrawLatex(5,hPedNoise[i*16+j]->GetMaximum()*0.7,
			   "Known as hot cahnnel");
       }
     }
     c2 -> Print( psFileName2 );
     c2 -> Clear();
   }
   // Large
   for(int i=0;i<15;i++){
     c2 -> cd();
     c2 -> Divide(4,4);
     for(int j=0;j<16;j++){
       if( i==14 && 12<=j ) continue;
       c2 -> cd( j+1 );
       hPedNoise[i*16+j+2360] -> Draw();
       tline -> DrawLine(3,0,3,hPedNoise[i*16+j]->GetMaximum() );
     }
     c2 -> Print( psFileName2 );
     c2 -> Clear();
   }

   c2->Print( Form( "%s)", psFileName2 ) );

   // Save histos
   hfile -> cd();
   for(int icrystal=0;icrystal<DetNumber[0];icrystal++){
     hLaserRef[icrystal] -> Write();
     hPedNoise[icrystal] -> Write();
   }
   for(int ipin=0;ipin<4;ipin++) hPinDiode[ipin] -> Write();
   for(int icv=0;icv<10;icv++) hCV[icv] -> Write();
   hIntegratedADCSum -> Write();
   hEneSum -> Write();
   hEtSum -> Write();
   hEtSumLaser -> Write();
   hfile -> Close();

 }

 void SemiOnlinePlot::ReadLaserRefValue(){

   int tmpi;

   FILE *fp;
   fp = fopen("/share/apps/production/2012VME/sumup/param/LaserRefMean.txt","rt");
   
   if( fp == NULL ){
     printf("File not found for Laser gain reference... set 1.\n");
     for(int icrystal=0;icrystal<DetNumber[0];icrystal++){
       LaserRefMean[icrystal]  = 1;
       LaserRefSigma[icrystal] = 1;
     }     

   }else{
     printf("File for Laser gain reference was opend.\n");
     for(int icrystal=0;icrystal<DetNumber[0];icrystal++){
       fscanf(fp,"%d %f %f\n",
	      &tmpi,&LaserRefMean[icrystal],&LaserRefSigma[icrystal]);
       //       printf("%d %f %f\n",
       //	      icrystal,LaserRefMean[icrystal],LaserRefSigma[icrystal]);
     }
   }

   return;

}
