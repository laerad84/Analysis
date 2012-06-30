
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <cstdio>

void simpleViewer(){
  std::string ANAPATH = std::getenv("ANALISISLIB");
  ANAPATH += "/lib/libtest.so";
  gSystem->Load(ANAPATH.c_str());

  const int nCH = 2716;
  const int nPoint  = 2; 
  E14ReadSumFile* Reader[nPoint];
  int RunNumber[nPoint] = {4796,4797};
  std::string filename = "/Volume0/ExpData/2012_Feb_Beam/Sum%s.root";

  for(int ipoint  =0 ;ipoint < nPoint;++ipoint){
    Reader[ipoint] = new E14ReadSumFile();
    Reader[ipoint]->Add(filename.c_str(),RunNumber[ipoint]);     
  }

  /*
  E14ReadSumFile* Reader2  = new E14ReadSumFile();
  Reader->Add("/Volume0/ExpData/2012_Feb_Beam/Sum4798.root");
  */

  TH1D* hist[nPoint][nCH];
  for( int ipoint  = 0; ipoint < nPoint ; ++ipoint){
    for( int ich = 0; ich < nCH; ++ich){
      hist[ipoint][ich] = new TH1D(Form("hist_%d_%d",ich,ipoint),
				   Form("hist_%d_%d",ich,ipoint),
				   1600,
				   0,
				   16000);
    }
  }
  
  for( int ipoint  = 0; ipoint < nPoint; ++ipoint){
    for( int ievent = 0; ievent < Reader[ipoint]->GetEntries(); ++ievent){
      Reader[ipoint]->GetEntry(ievent);
      for(int idigi = 0; idigi < Reader[ipoint]->CsiNumber; idigi++)F{
	  if(Reader[ipoint]->CsiEne[idigi] > 10){
	    hist[ipoint][ Reader[ipoint]->CsiModID[idigi] ]->Fill(Reader[ipoint]->CsiEne[idigi]);
	  }
	}
    }
  }

  TGraph* gr[nPoint-1];
  for( int i = 0; i< nPoint; ++i){
    gr[i]= new TGraph();
  }
  
  for( int ipoint =1; ipoint < nPoint; ipoint++){
    for( int ich = 0; ich < nCH; ich++){
      if(hist[ipoint][ich]->GetMean()> 10 &&
	 hist[0][ich]->GetMean()     > 10 ){
	gr[ipoint-1]->SetPoint(gr->GetN(),ich,hist[ipoint][ich]->GetMean()/hist[ipoint][ich]->GetMean());
      }
    }
  }
  
  for( int ipoint = 1; ipoint < nPoint; ipoint++){
    gr[ipoint]->SetMarkerStyle(3);
    gr[ipoint]->SetMarkerColor(ipoint);
  }
  

  gr[0]->Draw("AP");
}



	  

