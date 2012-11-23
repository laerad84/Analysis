#include "TApplication.h"
#include "TCanvas.h"
#include "CsIPoly.h"
#include "OEVPoly.h"
#include "CC03Poly.h"
#include "MBPoly.h"
#include "FBPoly.h"

int main( int argc, char** argv ){
  TApplication* app = new TApplication("app", &argc, argv );
  TCanvas* can = new TCanvas("can", "" ,800,800);
  can->Divide( 2,2 );
  
  CsIPoly*  Csi  = new CsIPoly("csi","csi");
  CC03Poly* Cc03 = new CC03Poly("cc03","cc03");
  OEVPoly*  Oev  = new OEVPoly("oev","oev");
  MBPoly*   Mb   = new MBPoly("mb","mb");
  FBPoly*   Fb   = new FBPoly("fb","fb");

  CsIPoly*  Csi1  = new CsIPoly("csi1","csi");
  CC03Poly* Cc031 = new CC03Poly("cc031","cc03");
  OEVPoly*  Oev1  = new OEVPoly("oev1","oev");
  MBPoly*   Mb1   = new MBPoly("mb1","mb");
  FBPoly*   Fb1   = new FBPoly("fb1","fb");

  for( int index = 0;index <  CsIPoly::numberOfCsI ; index++ ){
    Csi->Fill( index, (double)( index %10 ));
  }
  for( int index = 0; index < CC03Poly::numberOfCC03; index++){
    Cc03->Fill( index, (double)( index %10 ));
  }
  for( int index = 0; index < MBPoly::numberOfMB; index++){
    Mb->Fill( index , (double)( index % 10 ));
  }
  for( int index = 0; index < FBPoly::numberOfFB; index++){
    Fb->Fill( index, (double)( index % 10 ));
  }
  for( int index = 0; index < OEVPoly::numberOfOEV; index++){
    Oev->Fill( index, (double)(index % 10));
  }

  for( int index = 0;index <  CsIPoly::numberOfCsI ; index++ ){
    Csi1->Fill( index, (double)( index ));
  }
  for( int index = 0; index < CC03Poly::numberOfCC03; index++){
    Cc031->Fill( index, (double)( index ));
  }
  for( int index = 0; index < MBPoly::numberOfMB; index++){
    Mb1->Fill( index , (double)( index ));
  }
  for( int index = 0; index < FBPoly::numberOfFB; index++){
    Fb1->Fill( index, (double)( index ));
  }
  for( int index = 0; index < OEVPoly::numberOfOEV; index++){
    Oev1->Fill( index, (double)(index));
  }

  can->cd(1);
  Oev->Draw("col L");
  Csi->Draw("same col L");
  Cc03->Draw("same col L");
  can->cd(2);
  Mb->Draw("col L");
  Fb->Draw("same col L");

  can->cd(3);
  Oev1->Draw("col L");
  Csi1->Draw("same col L");
  Cc031->Draw("same col L");
  can->cd(4);
  Mb1->Draw("col L");
  Fb1->Draw("same col L");
  
  app->Run();
  return 0; 
}
