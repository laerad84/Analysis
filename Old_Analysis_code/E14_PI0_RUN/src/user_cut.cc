#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "csimap/CsiMap.h"
#include "TMath.h"
#include <iostream>
#include <list>

enum{ ENERGY_CUT=0, CSI_FIDUCIAL_CUT, VERTEX_CUT, PT_CUT, 
      COLLINEAR_CUT, DISTANCE_CUT, E_TOTAL_CUT, E_THETA_CUT,
      E_RATIO_CUT, PI_KINE_CUT, SHAPE_CHI2_CUT };

int standardCut(Pi0 const &pi0);
int shapeCut(Pi0 const &pi);
int vetoCut(E14GNAnaDataContainer const &data);
int csiVetoCut(E14GNAnaDataContainer const &data,Pi0 const &pi);

bool pi0kine_cut(Pi0 const &pi0);
double getCsiThreshold(double distance);
bool cutline(double x1, double y1, double x2, double y2,double var1, double var2);


void user_cut(E14GNAnaDataContainer &data,std::list<Pi0> const &piList){
  Pi0 const &pi = piList.front();

  int &iCut = data.CutCondition;
  iCut = 0;
  iCut += shapeCut(pi);
  iCut += standardCut(pi);

  int &vCut = data.VetoCondition;
  vCut = 0;
  vCut += vetoCut(data);
  vCut += csiVetoCut(data,pi);
}


int standardCut(Pi0 const &pi0){
  
  Gamma const &g1 = pi0.g1();
  Gamma const &g2 = pi0.g2();
  
  int iCut  = 0;

  
  double Egamma_Min = 0.1 * 1000; //GeV to MeV
  double Egamma_Max = 2.0 * 1000; //GeV to MeV
  //  double Egamma_Max = 10.0;
  
  //  double R_Min =175;
  double XY_Min =150;
  double R_Max =850;
  //  double R_Max =1350;  // for step2

  double Zvert_Min =3000;
  double Zvert_Max =5000;

  double Pt_Min =0.13 * 1000; //GeV to MeV
  double Pt_Max =0.25 * 1000; //GeV to MeV

  double Acop_Min =30.; // deg

  double gDist_Min =300.; //mm

  double TotE_Min =0.5 * 1000; // GeV to MeV

  double Etheta_Min =2.5 * 1000  ; //GeV to MeV 
  
  double Eratio_Min = 0.2;

  //Energy Cut
  if( g1.e()  < Egamma_Min || g2.e() < Egamma_Min ||
      g1.e() > Egamma_Max || g2.e() > Egamma_Max)
    iCut += 1<<ENERGY_CUT ;
  
  //Fiducial in CsI
  double R1 = g1.pos().perp();
  double R2 = g2.pos().perp();
  //  if ( R1 < R_Min || R1 > R_Max || R2 < R_Min || R2 > R_Max)  
  if( R1>R_Max ||  (fabs(g1.x())<XY_Min && fabs(g1.y())<XY_Min) ||
      R2>R_Max ||  (fabs(g2.x())<XY_Min && fabs(g2.y())<XY_Min) )
    iCut += 1<<CSI_FIDUCIAL_CUT;
  
  // Vertex Cut
  double RecVertZ = pi0.recZ();
  if ( RecVertZ < Zvert_Min || RecVertZ > Zvert_Max)
    iCut += 1<<VERTEX_CUT ;

  // PT Cut
  double RecPi0Pt = pi0.p3().perp();
  if ( RecPi0Pt < Pt_Min || RecPi0Pt > Pt_Max) 
    iCut += 1<<PT_CUT;

  // Acp_angle  Cut
  double acop;
  double  ct = (g1.x()*g2.x()+g1.y()*g2.y())/R1/R2;
  if(TMath::Abs(ct) >= 1){
    acop = 0.0;
  }
  else{
    acop = TMath::Pi() - TMath::ACos(ct);
    acop = acop * 180. / TMath::Pi();
  }
  
  double Acop_angle = acop;
  if (Acop_angle < Acop_Min) iCut += 1<<COLLINEAR_CUT;
  
  //Two gamma distance
  CLHEP::Hep3Vector gdist = g1.pos()-g2.pos();
  double gamma_dist = gdist.mag();
                              
  if (gamma_dist < gDist_Min) iCut += 1<<DISTANCE_CUT;
    
  //Total energy of two gammas
  double E_tot = g1.e() + g2.e();
  if (E_tot < TotE_Min)  iCut += 1<<E_TOTAL_CUT;
  
  // In_ang cut : Odd-pair cut
  double  In_Ang1 = TMath::ATan(R1/(g1.z()-pi0.recZ()));
  In_Ang1 = In_Ang1 * 180 / TMath::Pi();
  double  Etheta1 = g1.e() * In_Ang1;
  double  In_Ang2 = TMath::ATan(R2/(g2.z()-pi0.recZ()));
  In_Ang2 = In_Ang2 * 180 / TMath::Pi();
  double  Etheta2 = g2.e() * In_Ang2;


  if (Etheta1 < Etheta_Min || Etheta2 < Etheta_Min)
    iCut += 1<<E_THETA_CUT;
  
  // E_ratio_cut
  double Eratio = TMath::Min(g1.e(),g2.e())
    / TMath::Max(g1.e(),g2.e());
  if (Eratio < Eratio_Min)    iCut += 1<<E_RATIO_CUT;

  // Pi0 Kinematic Cut
  if(pi0kine_cut(pi0))  {
    iCut += 1<<PI_KINE_CUT;
  }

  return iCut;  
}


bool cutline(double x1, double y1, double x2, double y2,
		     double var1, double var2)
{
  double slope = ( y1 - y2 ) / ( x1 - x2 );
  double intercept = y1 - slope * x1 ;
  double y = slope * var1 + intercept;
  if ( var2 > y ) return kTRUE; // above the line
  return kFALSE;
}


bool pi0kine_cut(Pi0 const &pi0)
{
  double cpi0pt=pi0.p3().perp()*1e-3;//RecPi0Pt*1e-3;
  double cpi0pz=pi0.p3().z()*1e-3;//RecPi0Pz*1e-3;
  double cpi0recz=pi0.recZ()*1e-1;//RecVertZ*1e-1;
  double cpi0e=pi0.e()*1e-3;//RecPi0E*1e-3;

  bool flag = kFALSE;
  double pratio = cpi0pt/cpi0pz;
  if (( cpi0recz < 400 )&&( pratio < 0.1 )) return kTRUE;
  if (( cpi0recz > 400 )&&
      !( cutline(400., 0.1, 500., 0.15, cpi0recz, pratio) )) return kTRUE;
  if ( !(cutline(300., 0.8, 500., 0.4, cpi0recz, cpi0e)) ) return kTRUE;
  if ( (cutline(300., 0.2, 500., 0.34, cpi0recz, pratio)) ) return kTRUE;
  return kFALSE;

}


int shapeCut(Pi0 const &pi){
  double cutVal = 2.5;
  if(pi.g1().chisq() > cutVal || pi.g2().chisq() > cutVal )
    return 1<<SHAPE_CHI2_CUT;
  
  return 0;
}

int vetoCut(E14GNAnaDataContainer const &data){
  double threshold[100];
  for(int i=0;i<100;i++) threshold[i] = -1;
  threshold[DigiReader::CC01]=3;
  threshold[DigiReader::CC02]=3;
  threshold[DigiReader::CC03]=3;
  threshold[DigiReader::CC04]=3;
  threshold[DigiReader::CC05]=3;
  threshold[DigiReader::CC06]=3;
  threshold[DigiReader::CBAR]=3;
  threshold[DigiReader::FBAR]=3;
  threshold[DigiReader::BHPV]=3;
  threshold[DigiReader::CV]=3;
  threshold[DigiReader::BCV]=3;
  threshold[DigiReader::BHCV]=3;
  
  int iVeto = 0;
  for(int i=0;i<data.VetoNumber;i++){
    int detid = data.VetoDetID[i];
    double thre = threshold[detid];
    if(thre<0){
      static bool isFirstOne = true;
      if(isFirstOne){
	GsimMessage::getInstance()->
	  report("warning",Form("there is the unknown detector (detID==%d)",detid));
	isFirstOne = false;
      }
    }else if(thre<data.VetoEne[i]){
      iVeto = ( iVeto | (1<<detid) );
    }
  }
  return iVeto;
}

double getCsiThreshold(double distance){
  double constant = 3.73275;
  double slope = 0.0074758;
  double threshold = exp(constant-slope*distance);

  if(threshold>10) return 10;
  if(threshold<1.5) return 1.5;
  return threshold;
}


int csiVetoCut(E14GNAnaDataContainer const &data,Pi0 const &pi){
  bool isVeto = false;

  for(int icsi = 0; icsi < data.CsiNumber; icsi++){
    int id = data.CsiModID[icsi];
    double edep  = data.CsiEne[icsi];
    
    bool find = false;
    double cx[2]={0},cy[2]={0};
    for(int igam=0; igam<2 ; igam++){
      Gamma const &gam = (igam==0)? pi.g1() : pi.g2();
      cx[igam] = gam.x();
      cy[igam] = gam.y();
      for(int i=0; i<gam.clusterIdVec().size()&&!find; i++){
	if(gam.clusterIdVec().at(i)==id)  find = true; 
      }
    }
    
    if(find) continue;    
    double x,y,w;
    CsiMap::getCsiMap()->getXYW(id,x,y,w);
    double dist = std::min(sqrt(pow(x-cx[0],2)+pow(y-cy[0],2)),
			   sqrt(pow(x-cx[1],2)+pow(y-cy[1],2)));

    if(edep > getCsiThreshold(dist)){
      isVeto = true;
      //      break;
    }
  }
  if(isVeto) return 1<<DigiReader::CSI;
  return 0;
}


