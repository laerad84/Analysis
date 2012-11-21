#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <list>
#include "cluster/Cluster.h"
#include "gamma/Gamma.h"
#include "klong/Klong.h"


#include "gnana/E14GNAnaDataContainer.h"
#include "gnana/E14GNAnaFunction.h"
#include "rec2g/Rec2g.h"

#include "User_Function.h"

#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"



class ClusterTime{
 public:
  clusterTime();
  ~ClusterTime();
  bool operator  < (const ClusterTime& ) const;
  bool operator == (const ClusterTime& ) const;
  
  int id() const {return m_id;}
  int status() const {return m_status; }
  double RThreshold()  const { return m_Rthreshold; }
  double ClusterTime() const { returm m_MeanTime; }

  const std::vector<int>& clusterIDVec()           const {return m_cidVec; }
  const std::vector<double>& clusterTimeDeltaVec() const {return m_cTimeDeltaVec; }
  const std::vector<double>& clusterPiVec()        const {return m_cPhiVec; }
  const std::vector<double>& clusterRVec()         const {return m_cRVec; }
  
 private:
  int    m_id;
  int    m_status;
  double m_Rthreshold;
  double m_MeanTime;
  double m_RCluster;
  double m_PhiCluster;

  std::vector<double> m_cidVec;
  std::vector<double> m_cTimeDeltaVec;
  std::vector<double> m_cRvec;
  std::vector<double> m_cPhiVec;

};

ClusterTime::ClusterTime() 
: m_id(-1),
  m_status = -1,
  m_Rthreshold = TMath::Sqrt(2)*5./6.*25.,
  m_MeanTime = 0,
  m_nMeanCrystal = 0 
{
  ;
}
  ClusterTime::~ClusterTime()
{
    
    ;
  }
bool ClusterTime::operator <  (const ClusterTime& clus ) const
{
  if( compare( clus ) < 0){
    return true;
  }else{
    return false;
  }
}
bool ClusterTime::operator == (const ClusterTime& clus ) const 
{
  if( compare( clus ) == 0){
    return true; 
  }else{ 
    return false;
  }
}
std::ostream& operator << ( std::ostream& out, const Cluster& clus )
{
  out << "ClusterTime"  << "\n"
      << "id : " << clus.m_id << "\n"
      << "status : " << clus.m_status << "\n"
      << "Rthreshold :" << clus.m_Rthreshold << "\n"
      << "MeanTime: " << clus.m_MeanTime << "\n"
      << std::flush;
  out << "idVec : " ;
  for( int i = 0; i < clus.m_cidVec.size() ;i++){
    out << clus.m_cidVec[i] << " " ;
  }
  out << "\n";
  out << "RVec : " ;
  for( int i = 0; i < clus.m_cRVec.size(); i++){
    out << clus.m_cRVec[i] << " " ;
  }
  out << "\n";
  out << "PhiVec : " ; 
  for( int i = 0; i < clus.m_cPhiVec.size(); i++){
    out << clus.m_cPhiVec[i] << " ";
  }
  out << "\n" << std::flush;
  return out; 
}
int ClusterTime::Compare( const ClusterTime& clust ) const 
{
  if (m_MeanTime < clus.m_MeanTime ){
    return 1; 
  }else if( m_energy > clus.m_energy ){
    return -1;
  }
  return 0; 
}

class ClusterTimeAnalyzer {

  ClusterTimeAnalyzer();
  virtual ~ClusterTimeAnalyzer();
  
  const Double_t RadiusThreshold;

  void Init();
  void Reset();
  int  EventProcess ( int ievent );
  long Loop();

  Double_t clusterPosX;
  Double_t clusterPosY;
  Double_t clusterPosR;
  Double_t clex;
  Double_t cley;
  Int_t    Theta_Index;
  Double_t ClusterCenterEnergy;
  Double_t ChannelTimeOffset[ 2716 ];
  
  Double_t RInCluster;
  Double_t PhiInCluster;
  Double_t RInCsI;
  Double_t PhiInCsI;

  std::list<ClusterTime>& clusterTimeList() { return m_clusterTimeList; }
  std::list<ClusterTime>& AnlyzeClusterTime(){ std::list<ClusterTime> clist }

};

std::list<ClusterTime>& ClusterTimeAnalyzer::Convert(std::list<Cluster> clist)
{
  m_clusterTimeList.clear();
  for( std::list<Cluster>::iterator itCluster;
       itCluster != clist.end();
       itCluster++){

    ClusterTime clusTime;    
    





    m_clusterTimeList.push_back(clusTime);
    
  }
  return m_clusterTimeList;
}

Cluster& ClusterTimeAnalyzer::Process( Cluster clus ){
  

}






ClusterTimeAnalyzer::ClusterTimeAnalyzer: RadiusThreshold( 25.*5./6.*TMath::Sqrt(2)){
  Init();
}
ClusterTimeAnalyzer::~ClusterTimeAnalyzer(){
  ;
}
void ClusterTimeAnalyzer::Reset(){
  clusterPosX = 0;
  clusterPosY = 0;
  clusterPosR = 0;
  clex        = 0;
  cley        = 0;
  ClusterCenterEnergy = 0;
  ClusterMeanTime = 0.;
  nCrystalTime = 0;
}
void ClusterTimeAnalyzer::Init(){
  for( int i = 0; i< 2716; i++){
    ChannelTimeOffset[i] = 0; 
  }
  Reset();
}
int ClusterTimeAnalyzer::EventProcess( int ievent ){
  ch->GetEntry( ievent );
  std::list<Cluster> clist; 
  data->GetData( clist );
  Double_t x,y;  
  for( std::list<Cluster>::iteartor itCluster = clist.begin();
       itCluster != clist.end();
       itCluster++ ){
    Double_t weightX = 0; 
    Double_t weightY = 0;
    Double_t TotalEInCluster = 0;
    
    for( int iCrystal = 0; iCrystal < (*itCluster).clusterIdVec().size(); iCrystal++){
      Int_t IDInCluster = (*ItCluster).clusterIdVec()[iCrystal];
      Double_t TimeInCluster = (*itCluster).clusterTimeVec()[iCrystal];
      Double_t EnergyInCluster = (*itCluster).clusterEVec()[iCrystal];
      handler->GetMetricPosition( IDInCluster, x, y);
      weightX += x*EnergyInCluster;
      weightY += y*EnergyInCluster;
      TotalEInCluster += EnergyInCluster;
    }
    RInCsI = TMath::Sqrt( weightX*weightX + weightY*weightY );
    PhiInCsI = TMath::TMath::ATan2( weightY, weightX);    
    for( int iCrystal = 0; iCrystal < (*itcluster).clusterIDVec().size(); iCrystal++){
      Int_t IDInCluster = (*itCluster).clusterIdVec()[iCrystal];
      Double_t TimeInCluster = (*itCluster).clusterTimeVec()[iCrystal];
      Double_t EnergyInCluster = (*itCluster).clusterEVec()[iCrystal];
      handler->GetMetricPosition( IDInCluster, x, y );
      // Define ClusterMeanTime := Mean Time of near Energy Center Crystal( R <=25*(6./5.)*TMath::Sqrt(2) );
      // Range set to statistically static condition nCrystalTime = 4-6;;
      if( TMath::Sqrt( (x-weightX)*(x-weightX ) + (y-weightY)*(y-weightY)) <= RadiusThreshold ){
	ClusterMeanTime += TimeInCluster - ChannelTimeOffSet[ IDInCluster ];
	nCrystalTime++;
      }
    }     
  }   
}



