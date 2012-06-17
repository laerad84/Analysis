//============================================================================
// Name        : CsIImage.cpp
// Author      : Takahiko Masuda
// Version     : 1.0
// Copyright   : taka
// Date        : 2010/10/25
//============================================================================
#include "CsIImage.h"
#if !defined(__CINT__)
ClassImp(CsIImage)
#endif

CsIImage::CsIImage(IDHandler* idHandler) {
  
  smallHist   = new TH2D( "", "", 48, -600., 600., 48, -600., 600.);
  largeHist   = new TH2D( "", "", 36, -900., 900., 36, -900., 900.);
  CC03Hist[0] = new TH2D( "", "",  1,   80., 100.,  8,  -80., 100.);
  CC03Hist[1] = new TH2D( "", "",  8,  -80., 100.,  1, -100., -80.);
  CC03Hist[2] = new TH2D( "", "",  1, -100., -80.,  8, -100.,  80.);
  CC03Hist[3] = new TH2D( "", "",  8, -100.,  80.,  1,   80., 100.);
  
  handler = idHandler;
  
  int crystalID = 0;
  double x, y;
  for( ; crystalID<numberOfSmall; crystalID++ ){
    handler->GetMetricPosition( crystalID, x, y);
    box[crystalID] = new TBox( x-12.5, y-12.5, x+12.5, y+12.5 );
  }
  for( ; crystalID<numberOfCrystals; crystalID++ ){
    handler->GetMetricPosition( crystalID, x, y);
    box[crystalID] = new TBox( x-25.0, y-25.0, x+25.0, y+25.0 );
  }
  for( int index=0; index<numberOfCC03; index++ ){
    handler->GetMetricPosition( 3000+index, x, y );
    switch( index/8 ){
    case 0:
    case 2:
      box[index+numberOfCrystals] = new TBox( x-10.0, y-11.25, x+10.0, y+11.25 );
      break;
    case 1:
    case 3:
      box[index+numberOfCrystals] = new TBox( x-11.25, y-10.0, x+11.25, y+10.0 );
      break;
    }
  }
  
  
  for( crystalID=0; crystalID<numberOfCrystals+numberOfCC03; crystalID++ ){
    box[crystalID]->SetFillStyle(0);
  }
  
  smallHist->SetStats(0);
  largeHist->SetStats(0);
  
}






	
CsIImage::~CsIImage() {
  delete smallHist;
  delete largeHist;
  delete [] box;
}




bool CsIImage::Fill( int crystalID, double w ){
  if( !(handler->isChannel(crystalID)) ) return false;
	
  double x, y;
  handler->GetMetricPosition( crystalID, x, y);
  //  std::cout << x << " " << y << std::endl;
  
  if( handler->isSmall(crystalID) ){
    smallHist->Fill( x, y, w );
  }else if( handler->isLarge(crystalID) ){
    largeHist->Fill( x, y, w );
  }
  else if( handler->isCC03(crystalID) ){
    CC03Hist[(crystalID-3000)/8]->Fill(x,y,w);
  }

  return true;
  
}





bool CsIImage::SetFillColor( int crystalID, Color_t color){
  if( !(handler->isChannel(crystalID)) ) return false;
	
  box[crystalID]->SetFillColor( color );
  box[crystalID]->SetFillStyle(1001);
  return true;
}



bool CsIImage::GetBin( int crystalID, int& xbin, int& ybin ){
  if( !(handler->isChannel(crystalID)) ){
    xbin = 0;
    ybin = 0;
    return false;
  }

  double x, y;
  handler->GetMetricPosition( crystalID, x, y );

  if( handler->isSmall(crystalID) ){
    xbin = (int)((x+587.5)/25+1);
    ybin = (int)((y+587.5)/25+1);
  }else if( handler->isLarge(crystalID) ){
    xbin = (int)((x+875)/50+1);
    ybin = (int)((y+875)/50+1);
  }else if( handler->isCC03(crystalID) ){
    switch( (crystalID-3000)/8 ){
    case 0:
      xbin = 1;
      ybin = (int)((y+68.75)/22.5+1);
      break;
    case 1:
      xbin = (int)((x+68.75)/22.5+1);
      ybin = 1;
      break;
    case 2:
      xbin = 1;
      ybin = (int)((y+88.75)/22.5+1);
      break;
    case 3:
      xbin = (int)((x+88.75)/22.5+1);
      ybin = 1;
    }
  }

  return true;
}



const double CsIImage::GetCrystalContent( int crystalID ){
  if( !(handler->isChannel(crystalID)) ) return -1;

  int xbin, ybin;
  GetBin( crystalID, xbin, ybin );
  if( handler->isSmall(crystalID) ){
    return smallHist->GetCellContent( xbin, ybin );
  }else if( handler->isLarge(crystalID) ){
    return largeHist->GetCellContent( xbin, ybin );
  }
  else if( handler->isCC03(crystalID) ){
    return CC03Hist[(crystalID-3000)/8]->GetCellContent( xbin, ybin );
  }

  return -1;
}
  




void CsIImage::FillRandom( const char* fname, Int_t ntimes) {
	smallHist->FillRandom( fname, ntimes );
	largeHist->FillRandom( fname, ntimes );
	for( int index=0; index<4; index++ ){
	  CC03Hist[index]->FillRandom( fname, ntimes );
	}
}

void CsIImage::Draw(char* m_option){
  double maxValue=0, minValue=100000000, content=0;
  
  for( int crystalID=0; crystalID<3032; crystalID++ ){
    content = GetCrystalContent( crystalID );
    if( content!=-1 ){
      if( maxValue<content ) maxValue=content;
      if( minValue>content ) minValue=content;
    }
  }

  smallHist->SetMaximum( maxValue );
  smallHist->SetMinimum( minValue );
  largeHist->SetMaximum( maxValue );
  largeHist->SetMinimum( minValue );
  for( int index=0; index<4; index++ ){
    CC03Hist[index]->SetMaximum( maxValue );
    CC03Hist[index]->SetMinimum( minValue );
  }

  largeHist->Draw(m_option);
  char* sameOption = Form("%s,same",m_option);
  smallHist->Draw(sameOption);
  for( int index=0; index<4; index++ ){
    CC03Hist[index]->Draw(sameOption);
  }
  for( int crystalID=0; crystalID<numberOfCrystals+numberOfCC03; crystalID++ ){
    box[crystalID]->Draw("1");
  }
}

void CsIImage::DrawWithRange(char* m_option,Double_t minUser, Double_t maxUser){
  double maxValue=0, minValue=100000000, content=0;

  for( int crystalID=0; crystalID<3032; crystalID++ ){
    content = GetCrystalContent( crystalID );
    if( content!=-1 ){
      if( maxValue<content ) maxValue=content;
      if( minValue>content ) minValue=content;
    }
  }

  smallHist->SetMaximum( maxUser );
  smallHist->SetMinimum( minUser );
  largeHist->SetMaximum( maxUser );
  largeHist->SetMinimum( minUser );
  for( int index=0; index<4; index++ ){
    CC03Hist[index]->SetMaximum( maxUser );
    CC03Hist[index]->SetMinimum( minUser );
  }

  largeHist->Draw(m_option);
  char* sameOption = Form("%s,same",m_option);
  smallHist->Draw(sameOption);
  for( int index=0; index<4; index++ ){
    CC03Hist[index]->Draw(sameOption);
  }
  for( int crystalID=0; crystalID<numberOfCrystals+numberOfCC03; crystalID++ ){
    box[crystalID]->Draw("1");
  }
}

void CsIImage::SetTitle( const char* title ){
  largeHist->SetTitle(title);
}



void CsIImage::Reset( void ){
  int crystalID = 0;
  int xbin, ybin;
  for( crystalID=0; crystalID<numberOfSmall; crystalID++ ){
    GetBin( crystalID, xbin, ybin );
    smallHist->SetBinContent( xbin, ybin, 0 );
  }
  for( ; crystalID<numberOfCrystals; crystalID++ ){
    GetBin( crystalID, xbin, ybin );
    largeHist->SetBinContent( xbin, ybin, 0 );
  }
  for( crystalID=3000; crystalID<3000+numberOfCC03; crystalID++){
    GetBin( crystalID, xbin, ybin );
    CC03Hist[(crystalID-3000)/8]->SetBinContent( xbin, ybin, 0 );
  }

}



		

/*


using namespace std;


int testCsIImage( void ){
	CsIImage* csi = new CsIImage();
	cerr << " CsIImage object is created." << endl;
	
	double x, y;
	
	for( int crystalID=0; crystalID<2716; crystalID++ ){
		csi->GetPosition( crystalID, x, y);
		csi->Fill( crystalID, (double)(crystalID%100) );
//		cerr << crystalID << " " << x << " " << y << endl;
	}
	
	
	for( int crystalID=100; crystalID<110; crystalID++ ){
		csi->SetFillColor( crystalID, 1 );
	}
		
	csi->Draw();
	return 0;
}





*/




