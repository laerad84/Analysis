#include "klong/Klong.h"
#include "gnana/E14GNAnaFunction.h"
#include "gnana/E14GNAnaDataContainer.h"
#include "csimap/CsiMap.h"
#include <iostream>
#include <list>

enum{ ENERGY_CUT=0, CSI_FIDUCIAL_CUT, DISTANCE_CUT, 
      KL_VERTEX_CUT, KL_PT_CUT,KL_CHI2_CUT,KL_2NDBEST_CHI2_CUT,
      PI0_MASS_CUT, PI0_Z_CUT, SHAPE_CHI2_CUT};

double getCsiThreshold(double distance);
int kinematicCut( std::vector<Klong> const &klVec);
int shapeCut(std::vector<Klong> const &klVec);
int vetoCut(E14GNAnaDataContainer const &data);
int csiVetoCut(E14GNAnaDataContainer const &data,
	       std::vector<Klong> const &klVec);


void user_cut(E14GNAnaDataContainer &data,std::vector<Klong> const &klVec){
  
  int &iCut = data.CutCondition;
  iCut = 0;
  iCut += shapeCut(klVec);
  iCut += kinematicCut(klVec);
  
  int &vCut = data.VetoCondition;
  vCut = 0;
  vCut += vetoCut(data);
  vCut += csiVetoCut(data,klVec);
}



int kinematicCut( std::vector<Klong> const &klVec){
  int iCut  = 0;

  double Egamma_Min = 0.10 * 1000; //GeV to MeV
  double Egamma_Max = 2.0 * 1000; //GeV to MeV
  
  double XY_Min =150;
  double R_Max =850;

  double gDist_Min =150.; //mm

  double Zvert_Min =3000;
  double Zvert_Max =5000;

  double Pt_Max = sqrt(1.25e-4) * 1000; //GeV to MeV
  
  double BestChi2_Max = 3.0;

  double Chi2_2ndChi2_Min = 4.0;

  double pi0Mass_diff_Max = 5.125; //MeV

  double pi0VertZ_diff_Max = 100; // mm

  int const nPi0 = klVec[0].pi0().size();
  double clusterPos[10][2]={{0}};
  for(int ipi=0;ipi<nPi0;ipi++){ 
    Pi0 const & pi0 = klVec[0].pi0().at(ipi);
    //pi0 mass cut
    double const pi0Mass = 134.976; //MeV
    if(fabs(pi0.m()-pi0Mass)>pi0Mass_diff_Max)
      iCut = (iCut | 1<<PI0_MASS_CUT );

    for(int igam = 0;igam<2;igam++){
      Gamma const &g = (igam==0) ? pi0.g1() : pi0.g2();      
      //Energy Cut
      if( g.e()  < Egamma_Min || g.e() > Egamma_Max )
	iCut = (iCut | 1<<ENERGY_CUT );

      //Fiducial in CsI      
      if( g.pos().perp()>R_Max ||  (fabs(g.x())<XY_Min && fabs(g.y())<XY_Min) )
	iCut = ( iCut | 1<<CSI_FIDUCIAL_CUT );

      clusterPos[ipi*2+igam][0] = g.x();
      clusterPos[ipi*2+igam][1] = g.y();
    }
  }

  //Two gamma distance
  for(int i=0; i<nPi0*2; i++){
    for(int j=i+1; j<nPi0*2; j++){
      double dist = sqrt( pow(clusterPos[i][0]-clusterPos[i][1],2)
			  + pow(clusterPos[i][0]-clusterPos[i][1],2) );
      if ( dist < gDist_Min ) iCut = ( iCut | 1<<DISTANCE_CUT );      
    }
  }

  // KL Vertex Cut
  double RecVertZ = klVec[0].vz();
  if ( RecVertZ < Zvert_Min || RecVertZ > Zvert_Max )
    iCut = ( iCut | 1<<KL_VERTEX_CUT );
  
  // KL PT Cut
  double RecPt = klVec[0].p3().perp();
  if ( RecPt > Pt_Max )   iCut = ( iCut | 1<<KL_PT_CUT ); 

  // KL Best Chi2 cut
  double Chi2 = klVec[0].chisqZ();
  if ( Chi2 > BestChi2_Max ) iCut = ( iCut | 1<<KL_CHI2_CUT );
    
  // KL 2nd Best - Best Chi2
  if(klVec.size()>1){
    double Chi2_2nd = klVec[1].chisqZ();
    if( Chi2_2nd - Chi2 < Chi2_2ndChi2_Min ) 
      iCut = ( iCut | 1<<KL_2NDBEST_CHI2_CUT );
  }
  
  if(klVec[0].deltaZ()>pi0VertZ_diff_Max)
    iCut = ( iCut | 1<<PI0_Z_CUT );
  
  return iCut;  
}


int shapeCut(std::vector<Klong> const &klVec){
  double cutVal = 2.5;
  int iCut = 0;
  for(std::vector<Pi0>::const_iterator it = klVec[0].pi0().begin();
      it!=klVec[0].pi0().end(); it++){
    if(it->g1().chisq() > cutVal || it->g2().chisq() > cutVal )
      iCut = ( iCut | 1<<SHAPE_CHI2_CUT );
  }
  
  return iCut;
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


int csiVetoCut(E14GNAnaDataContainer const &data,std::vector<Klong> const &klVec){
  bool isVeto = false;
  int nClusterCsi=0;
  int clusterCsiId[3000];
  int const nPi0 = klVec[0].pi0().size();
  double cx[10]={0},cy[10]={0};

  for(int ipi=0;ipi<nPi0;ipi++){ 
    for(int iclus=0; iclus<2 ; iclus++){
      Gamma const &clus = (iclus==0)
	? klVec[0].pi0().at(ipi).g1() : klVec[0].pi0().at(ipi).g2();
      cx[ipi*2+iclus] = clus.x();
      cy[ipi*2+iclus] = clus.y();
      for(int i=0; i<clus.clusterIdVec().size(); i++){
	clusterCsiId[nClusterCsi++] = clus.clusterIdVec().at(i);
      }
    }
  }
  

  for(int icsi = 0; icsi < data.CsiNumber; icsi++){
    int id = data.CsiModID[icsi];
    double edep  = data.CsiEne[icsi];
    
    bool find = false;
    for(int i=0; i<nClusterCsi; i++){
      if(clusterCsiId[i]==id){
	find = true;
	break;
      }
    }

    if(find) continue;
    
    double x,y,w;
    CsiMap::getCsiMap()->getXYW(id,x,y,w);

    double dist = sqrt(pow(x-cx[0],2)+pow(y-cy[0],2));
    for(int i = 1; i<2*nPi0; i++){
      double tmpdist = sqrt(pow(x-cx[i],2)+pow(y-cy[i],2));
      if(dist>tmpdist) dist = tmpdist;
    }
    
    if(edep > getCsiThreshold(dist)){
      isVeto = true;
      break;
    }
  }
  if(isVeto) return 1<<DigiReader::CSI;
  return 0;
}

double getCsiThreshold(double distance){
  double constant = 3.73275;
  double slope = 0.0074758;
  double threshold = exp(constant-slope*distance);

  if(threshold>10) return 10;
  if(threshold<1.5) return 1.5;
  return threshold;
}

