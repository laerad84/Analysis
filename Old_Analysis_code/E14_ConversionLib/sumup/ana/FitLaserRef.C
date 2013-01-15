void FitLaserRef(int RunNum=4107){

  TCanvas *pct = new TCanvas("pict","pict",800,600);

  TFile *hfile = new TFile( Form("/disk/kotodaq/2012Feb/run%d/RootFile/SemiOnlinePlot_%d.root",RunNum,RunNum) );
  
  TH1F *hLaserRef[2716];
  for(int icrystal=0;icrystal<2716;icrystal++){
    hLaserRef[icrystal] = (TH1F *)hfile -> Get( Form("hLaserRef_%d",icrystal) );
  }
  
  TF1 *fg = new TF1("fg","gaus");
  fg -> SetLineColor(2);

  float Mean[2716]  = {0};
  float Sigma[2716] = {0};
  for(int icrystal=0;icrystal<2716;icrystal++){
    hLaserRef[icrystal] = (TH1F *)hfile -> Get( Form("hLaserRef_%d",icrystal) );
    printf("Fitted for %d %d\n",icrystal,hLaserRef[icrystal]);

    pict -> cd();
    hLaserRef[icrystal] -> Fit("fg","NQ");
    //    hLaserRef[icrystal] -> Draw();
    //    pict -> Update();
    //    getchar();

    //    Mean[icrystal]  = fg -> GetParameter(1);
    //    Sigma[icrystal] = fg -> GetParameter(2);
    Mean[icrystal]  = hLaserRef[icrystal] -> GetMean();
    Sigma[icrystal] = hLaserRef[icrystal] -> GetRMS();

  }  

  FILE *fp;
  fp = fopen( Form("LaserRefMean_%d.txt",RunNum),"wt");
  for(int icrystal=0;icrystal<2716;icrystal++){
    fprintf(fp,"%d %f %f\n",icrystal,Mean[icrystal],Sigma[icrystal]);
  }

}
